/*
 *  _____ _____ ________  _
 * /  __//  __//  __/\  \//
 * | |  _|  \  |  \   \  / 
 * | |_//|  /_ |  /_  /  \ 
 * \____\\____\\____\/__/\\
 *
 * Graphics Environment for EXperimentations.
 *  Copyright (C) 2006 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: 
 *
 *     ALICE Project - INRIA
 *     INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 */

#include "sparse_matrix.h"
#include <algorithm>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <math.h>
// OpenMP is supported
#ifdef _OPENMP
   #include <omp.h>
#endif
#include "timing.h"

namespace Geex {

//_________________________________________________________

    bool SparseRowColumn::coeff_bin_search(index_t index, index_t* res) {
        // indexes are shifted by +1 to avoid going negative
        index_t first = 1;
        index_t last = nb_coeffs_;
        index_t mid;
        while (first <= last) {
            mid = (first + last) / 2;
            
            if (index > coeff_[mid-1].index) { 
                first = mid + 1; 
            } else if (index < coeff_[mid-1].index) {
                last = mid - 1;
            } else {
                *res = mid-1;
                return true;
            }
        }
        *res = first-1;
        return false;
    }
    
    bool SparseRowColumn::coeff_back_search(index_t index, index_t* res) {
        // indexes are shifted by +1 to avoid going negative
        index_t ii;
        index_t ind;
        for(ii = nb_coeffs_; ii>0; ii--){
            ind = coeff_[ii-1].index;
            if (ind <= index){
                if(ind == index){
                    *res = ii-1;
                    return true;
                }
                break;
            }
        }
        *res = ii;
        return false;
    }
    
    bool SparseRowColumn::coeff_search(index_t index, index_t* res) {
        index_t ii;
        index_t ind;
        for(ii=0; ii < nb_coeffs_; ii++){
            ind = coeff_[ii].index;
            if (ind >= index){
                if(ind == index){
                    *res = ii;
                    return true;
                }
                break;
            }
        }
        *res = ii;
        return false;
    }
    
    void SparseRowColumn::add(index_t index, double val, search_method_ptr search_method ) {
        index_t ii;
        if( (this->*search_method)(index, &ii)){
            coeff_[ii].a += val ;
        } else {
            nb_coeffs_++ ;
            if(nb_coeffs_ > capacity_) {
                grow() ;
            }
            // move up by one place all the elements with greater indexes
            size_t n = nb_coeffs_-ii-1;
            if( n != 0 ){
                Coeff*      dest = coeff_+ii+1;
                const Coeff* src = coeff_+ii;
                memmove(dest, src, n * sizeof(Coeff));
            }
             
            coeff_[ii].index = index;
            coeff_[ii].a = val;
        }
    }
    
    
    void SparseRowColumn::zero(index_t index, search_method_ptr search_method ) {
        index_t ii;
        if ((this->*search_method)(index, &ii)){
            size_t n = nb_coeffs_-ii-1;
            if( n != 0 ){
                Coeff*      dest = coeff_+ii;
                const Coeff* src = coeff_+ii+1;
                memmove(dest, src, n * sizeof(Coeff));
            }
            nb_coeffs_--;
        }
    }
    
    
    void SparseRowColumn::zeros(size_t nbIndexes, index_t* indexes) {
        
        // reverse order to avoid needing to sort the column after 
        // each removed element
        index_t coeff_cursor = nb_coeffs_; // element AFTER last examinated one
        bool need_sorting = false;
        
        for(index_t i = nbIndexes; i>0 ; i--){
            // Search for a_{index} in reverse order too. We start where the
            // last iteration left us, since we know that the elements already
            // seen are greater than the one we are looking for.
            // All these -1 are here because index_t is unsigned...
            index_t cur_index = indexes[i-1];
            index_t ii;
            for( ii=coeff_cursor; ii > 0; ii--) {
                index_t ind = coeff_[ii-1].index;
                if(ind <= cur_index) {
                    if (ind == cur_index){
                        Coeff& coeff = coeff_[ii-1] ;
                        
                        // replace by the last element and forget the last element's 
                        // existance...
                        coeff.a = coeff_[nb_coeffs_-1].a;
                        coeff.index = coeff_[nb_coeffs_-1].index;
                        nb_coeffs_--;
                        need_sorting = true;
                        //found : next loop begin with the next coeff
                        coeff_cursor = ii-1;
                    } else {
                        // not found : next loop begin with this one
                        coeff_cursor = ii; 
                    }
                    break;
                }
            }
            
            // all existing coeffs have been tested, we can skip the rest
            if (ii == 0){
                break;
            }
        }
        if ( need_sorting ){
            this->sort();
        }
    }
        
    double SparseRowColumn::value(index_t i, search_method_ptr search_method ){
        double val = 0.0;
        index_t ii;
        if ((this->*search_method)(i, &ii)){
            val = coeff_[ii].a;
        }
        return val;
    }
    
     
    void SparseRowColumn::dump_indexes_values(index_t* indexes, double* values ){
        for( index_t ii=0; ii < nb_coeffs_; ii++) {
            indexes[ii] = coeff_[ii].index;
            values[ii] = coeff_[ii].a;
        }
    }
        
    void SparseRowColumn::grow(size_t new_capacity) {
        gx_assert(new_capacity >= capacity_);
        Allocator<Coeff>::reallocate(
            coeff_, capacity_, new_capacity
        ) ;
        capacity_ = new_capacity;
    }

    void SparseRowColumn::grow() {
        grow(capacity_ *2);
    }

    class CoeffIndexCompare {
    public:
        bool operator()(const Coeff& c1, const Coeff& c2) {
            return c1.index < c2.index ;
        }
    } ;

    void SparseRowColumn::sort() {
        Coeff* begin = coeff_ ;
        Coeff* end   = coeff_ + nb_coeffs_ ;
        std::sort(begin, end, CoeffIndexCompare()) ;
    }

//_________________________________________________________

    SparseMatrix::SparseMatrix(size_t m, size_t n, Storage storage) {
        storage_ = NONE ;
        allocate(m,n,storage,false) ;
    }        

    SparseMatrix::SparseMatrix(const SparseMatrix* m, Storage storage){
        storage_ = NONE ;
        allocate(m->m(),m->n(),storage,false) ;
        if (m->has_symmetric_storage()){
            if ( m->columns_are_stored() ){
                for(index_t j = 0; j<m->n() ; j++ ){
                    const SparseRowColumn& col = m->column(j);
                    index_t k;
                    if (col.nb_coeffs() > 0){
                        for( k = 0; k < col.nb_coeffs()-1; k++ ){
                            this->add( col.coeff(k).index, j, col.coeff(k).a );
                            this->add( j, col.coeff(k).index, col.coeff(k).a );
                        }
                        this->add( col.coeff(k).index, j, col.coeff(k).a );
                        if( col.coeff(k).index != j ){
                            this->add( j, col.coeff(k).index, col.coeff(k).a );
                        }
                    }
                }
                return;
            }
            if ( m->rows_are_stored() ){
                for(index_t i = 0; i<m->m() ; i++ ){
                    const SparseRowColumn& row = m->row(i);
                    index_t k;
                    if(row.nb_coeffs() > 0){
                        for( k = 0; k < row.nb_coeffs()-1; k++){
                            this->add( i, row.coeff(k).index, row.coeff(k).a );
                            this->add( row.coeff(k).index, i, row.coeff(k).a );
                        }
                        this->add( i, row.coeff(k).index, row.coeff(k).a );
                        if( row.coeff(k).index != i){
                            this->add( row.coeff(k).index, i, row.coeff(k).a );
                        }
                    }
                }
                return;
            }
        } else {
            if ( m->columns_are_stored() ){
                for(index_t j = 0; j<m->n() ; j++ ){
                    const SparseRowColumn& col = m->column(j);
                    for( index_t k = 0; k < col.nb_coeffs(); k++ ){
                        this->add( col.coeff(k).index, j, col.coeff(k).a );
                    }
                }
                return;
            }
            if ( m->rows_are_stored() ){
                for(index_t i = 0; i<m->m() ; i++ ){
                    const SparseRowColumn& row = m->row(i);
                    for(index_t k = 0; k < row.nb_coeffs(); k++){
                        this->add( i, row.coeff(k).index, row.coeff(k).a );
                    }
                }
                return;
            }
        }
    }


    SparseMatrix::SparseMatrix(
        size_t n, Storage storage, bool symmetric_storage
    ) {
        storage_ = NONE ;
        allocate(n,n,storage,symmetric_storage) ;
    }        

    SparseMatrix::~SparseMatrix() {
        deallocate() ;
    }

    SparseMatrix::SparseMatrix() {
        m_ = 0 ;
        n_ = 0 ;
        diag_size_ = 0 ;

        row_ = nil ;
        column_ = nil ;
        diag_ = nil ;

        storage_ = NONE ;
        rows_are_stored_ = false ;
        columns_are_stored_ = false ;
        symmetric_storage_ = false ;
        symmetric_tag_ = false ;
    }

    size_t SparseMatrix::nnz() const {
        size_t result = 0 ;
        if(rows_are_stored()) {
            for(index_t i=0; i<m(); i++) {
                result += row(i).nb_coeffs() ;
            }
        } else if(columns_are_stored()) {
            for(index_t j=0; j<n(); j++) {
                result += column(j).nb_coeffs() ;
            }
        } else {
            gx_assert_not_reached ;
        }
        return result ;
    }
        
    void SparseMatrix::deallocate() {
        m_ = 0 ;
        n_ = 0 ;
        diag_size_ = 0 ;

        delete[] row_ ;
        delete[] column_ ;
        delete[] diag_ ;
        row_ = nil ;
        column_ = nil ;
        diag_ = nil ;

        storage_ = NONE ;
        rows_are_stored_ = false ;
        columns_are_stored_ = false ;
        symmetric_storage_ = false ;
    }

    void SparseMatrix::allocate(
        size_t m, size_t n, Storage storage, bool symmetric_storage
    ) {
        if(storage_ != NONE) {
            deallocate() ;
        }

        m_ = m ;
        n_ = n ;
        diag_size_ = gx_min(m,n) ;
        symmetric_storage_ = symmetric_storage ;
        symmetric_tag_ = false ;
        storage_ = storage ;
        switch(storage) {
        case NONE:
            gx_assert(false) ;
            break ;
        case ROWS:
            rows_are_stored_    = true ;
            columns_are_stored_ = false ;
            break ;
        case COLUMNS:
            rows_are_stored_    = false ;
            columns_are_stored_ = true ;
            break ;
        case ROWS_AND_COLUMNS:
            rows_are_stored_    = true ;
            columns_are_stored_ = true ;
            break ;
        }
        diag_ = new double[diag_size_] ;
        for(index_t i=0; i<diag_size_; i++) {
            diag_[i] = 0.0 ;
        }

        if(rows_are_stored_) {
            row_ = new Row[m] ;
        } else {
            row_ = nil ;
        }

        if(columns_are_stored_) {
            column_ = new Column[n] ;
        } else {
            column_ = nil ;
        }
    }

    void SparseMatrix::sort() {
        if(rows_are_stored_) {
            for(index_t i=0; i<m_; i++) {
                row(i).sort() ;
            }
        }
        if(columns_are_stored_) {
            for(index_t j=0; j<n_; j++) {
                column(j).sort() ;
            }
        }
    }

    void SparseMatrix::zero() {
        if(rows_are_stored_) {
            for(index_t i=0; i<m_; i++) {
                row(i).zero() ;
            }
        }
        if(columns_are_stored_) {
            for(index_t j=0; j<n_; j++) {
                column(j).zero() ;
            }
        }
        for(index_t i=0; i<diag_size_; i++) {
            diag_[i] = 0.0 ;
        }
    }
    
    
    // indexes MUST be sorted in ascending order
    void SparseMatrix::zero_rows(size_t nbIndexes, index_t* indexes){
        #ifdef _OPENMP
	omp_set_num_threads(omp_get_num_procs());
	#endif
	if(rows_are_stored_) {
	    #pragma omp parallel for
	    for(index_t i = 0; i<nbIndexes ; i++){
                row(indexes[i]).zero() ;
            }
        }
        if(columns_are_stored_) {
            #pragma omp parallel for
	    for(index_t j=0; j<n_; j++) {
                column(j).zeros(nbIndexes,indexes) ;
            }
        }
        for(index_t i = 0; i<nbIndexes ; i++){
            diag_[indexes[i]] = 0.0 ;
        }
    }
    
    
    // indexes MUST be sorted in ascending order
    void SparseMatrix::zero_columns(size_t nbIndexes, index_t* indexes){
	#ifdef _OPENMP
	omp_set_num_threads(omp_get_num_procs());
        #endif
	if(rows_are_stored_) {
            #pragma omp parallel for
	    for(index_t i=0; i<m_; i++) {
                row(i).zeros(nbIndexes,indexes) ;

	    }
	    
        }
        if(columns_are_stored_) {
            #pragma omp parallel for
	    for(index_t j = 0; j<nbIndexes ; j++){
                column(indexes[j]).zero() ;
            }
        }
        for(index_t i = 0; i<nbIndexes ; i++){
            diag_[indexes[i]] = 0.0 ;
        }
    }
    
    
    void SparseMatrix::zero_row(index_t i){
        if(rows_are_stored_) {
            row(i).zero() ;
        }
        if(columns_are_stored_) {
            for(index_t j=0; j<n_; j++) {
                column(j).zero(i) ;
            }
        }
        diag_[i] = 0.0 ;
    }
    
    void SparseMatrix::zero_column(index_t j){
        if(rows_are_stored_) {
            for(index_t i=0; i<m_; i++) {
                row(i).zero(j) ;
            }
        }
        if(columns_are_stored_) {
            column(j).zero() ;
        }
        diag_[j] = 0.0 ;
    }
    
   
    void SparseMatrix::clear() {
        if(rows_are_stored_) {
            for(index_t i=0; i<m_; i++) {
                row(i).clear() ;
            }
        }
        if(columns_are_stored_) {
            for(index_t j=0; j<n_; j++) {
                column(j).clear() ;
            }
        }
        for(index_t i=0; i<diag_size_; i++) {
            diag_[i] = 0.0 ;
        }
    }



//----------------------------------------------------------------------
//     Helpers for sparse matrix * vector product
//     according to storage modes:
//           rows or columns
//           symmetric
//----------------------------------------------------------------------

    static void mult_rows_symmetric(
        const SparseMatrix& A, const double* x, double* y
    ) {
        size_t m = A.m() ;
        memset(y, 0, m * sizeof(double)) ;
        for(index_t i=0; i<m; i++) {
            const SparseMatrix::Row& Ri = A.row(i) ;
            for(index_t  ij=0; ij<Ri.nb_coeffs(); ij++) {
                const Coeff& c = Ri.coeff(ij) ;
                y[i] += c.a * x[c.index] ;
                if(i != c.index) {
                    y[c.index] += c.a * x[i] ;
                }
            }
        }
    }

    static void mult_rows(
        const SparseMatrix& A, const double* x, double* y
    ) {
        size_t m = A.m() ;
        memset(y, 0, m * sizeof(double)) ;
        for(index_t i=0; i<m; i++) {
            y[i] = 0 ;
            const SparseMatrix::Row& Ri = A.row(i) ;
            for(index_t ij=0; ij<Ri.nb_coeffs(); ij++) {
                const Coeff& c = Ri.coeff(ij) ;
                y[i] += c.a * x[c.index] ;
            }
        }
    }


    static void mult_cols_symmetric(
        const SparseMatrix& A, const double* x, double* y
    ) {
        size_t n = A.n() ;
        memset(y, 0, n * sizeof(double)) ;
        for(index_t j=0; j<n; j++) {
            const SparseMatrix::Column& Cj = A.column(j) ;
            for(index_t ii=0; ii<Cj.nb_coeffs(); ii++) {
                const Coeff& c = Cj.coeff(ii) ;
                y[c.index] += c.a * x[j]       ;
                if(j != c.index) {
                    y[j]   += c.a * x[c.index] ;
                }
            }
        }
    }

    static void mult_cols(
        const SparseMatrix& A, const double* x, double* y
    ) {
        size_t n = A.n() ;
        memset(y, 0, n * sizeof(double)) ;
        for(index_t j=0; j<n; j++) {
            const SparseMatrix::Column& Cj = A.column(j) ;
            for(index_t ii=0; ii<Cj.nb_coeffs(); ii++) {
                const Coeff& c = Cj.coeff(ii) ;
                y[c.index] += c.a * x[j] ;
            }
        }
    }


    // general driver routine for 
    // sparse matrix * vector product
    void mult(const SparseMatrix& A, const double* x, double* y) {
        if(A.rows_are_stored()) {
            if(A.has_symmetric_storage()) {
                mult_rows_symmetric(A, x, y) ;
            } else {
                mult_rows(A, x, y) ;
            }
        } else {
            if(A.has_symmetric_storage()) {
                mult_cols_symmetric(A, x, y) ;
            } else {
                mult_cols(A, x, y) ;
            }
        }
    }

//----------------------------------------------------------------------
//     Helpers for transpose sparse matrix * vector product
//     according to storage modes:
//           rows or columns
//           symmetric
//     Note: in the symmetric case, A = A^t (and we can reuse the
//      non-transpose product)        
//----------------------------------------------------------------------

    static void mult_xpose_rows(
        const SparseMatrix& A, const double* x, double* y
    ) {
        size_t m = A.m() ;
        memset(y, 0, m * sizeof(double)) ;
        for(index_t i=0; i<m; i++) {
            const SparseMatrix::Row& Ri = A.row(i) ;
            for(index_t ij=0; ij<Ri.nb_coeffs(); ij++) {
                const Coeff& c = Ri.coeff(ij) ;
                y[c.index] += c.a * x[i] ;
            }
        }
    }

    static void mult_xpose_cols(
        const SparseMatrix& A, const double* x, double* y
    ) {
        size_t n = A.n() ;
        memset(y, 0, n * sizeof(double)) ;
        for(index_t j=0; j<n; j++) {
            y[j] = 0.0 ;
            const SparseMatrix::Column& Cj = A.column(j) ;
            for(index_t ii=0; ii<Cj.nb_coeffs(); ii++) {
                const Coeff& c = Cj.coeff(ii) ;
                y[j] += c.a * x[c.index] ;
            }
        }
    }


    // general driver routine for 
    // transpose sparse matrix * vector product
    void mult_transpose(
        const SparseMatrix& A, const double* x, double* y
    ) {
        if(A.rows_are_stored()) {
            if(A.has_symmetric_storage()) {
                mult_rows_symmetric(A, x, y) ;
            } else {
                mult_xpose_rows(A, x, y) ;
            }
        } else {
            if(A.has_symmetric_storage()) {
                mult_cols_symmetric(A, x, y) ;
            } else {
                mult_xpose_cols(A, x, y) ;
            }
        }
    }

    //_____________________________________________________________


/*************************************/

//----------------------------------------------------------------------
//     Helpers for sparse matrix * sparse matrix vector product
//     according to storage modes and product type.
//     Note: in some cases, we need to rebuild the matrix to have a suitable
//           storage type (non symmetric, and adapted to the multiplication
//           methode).
//----------------------------------------------------------------------

void mult_R_R(SparseRowColumn& (SparseMatrix::*rowcol_a)(index_t), SparseRowColumn& (SparseMatrix::*rowcol_b)(index_t), SparseMatrix& A, SparseMatrix& B, SparseMatrix& C){
    index_t nb_row = C.m();
    index_t nb_col = C.n();
    
    double val;
    index_t jj;
    
    size_t nb_coeffs_Ca, nb_coeffs_Cb;
    index_t ka, kb;
    
    index_t* look_up = new index_t[nb_col];
    for(index_t li = 0; li<nb_col ; li++ ){
        look_up[li] = 0;
    };
    index_t icoeff;
    
    index_t maxnnz = 8 + nb_col/8 ;
    Coeff* tmp_row = Allocator<Coeff>::allocate(maxnnz) ;
    index_t len = 0;
    
    for(index_t ii = 0; ii < nb_row; ii++ ){
        const SparseRowColumn&  RCa = (A.*rowcol_a)(ii);
        nb_coeffs_Ca = RCa.nb_coeffs();
        const Coeff* ca = &(RCa.coeff(0));
        for( ka = 0; ka < nb_coeffs_Ca ; ka++ ){
            val = ca->a;
            jj = ca->index;
            const SparseRowColumn&  RCb = (B.*rowcol_b)(jj);
            nb_coeffs_Cb = RCb.nb_coeffs();
            const Coeff* cb = &(RCb.coeff(0));
            for( kb = 0; kb < nb_coeffs_Cb ; kb++ ){
                icoeff = look_up[cb->index];
                if (icoeff == 0){
                    if ( len >= maxnnz ){
                        Allocator<Coeff>::reallocate(tmp_row, maxnnz, maxnnz*2 );
                        maxnnz = maxnnz*2;
                    }
                    tmp_row[len].index = cb->index;
                    tmp_row[len].a = val*cb->a;
                    look_up[cb->index] = len+1;
                    len++;
                } else {
                    tmp_row[icoeff-1].a += val*cb->a;
                }
                cb++;
            }
            ca++;
        }
        
        if(C.rows_are_stored() && !C.columns_are_stored() && !C.has_symmetric_storage()){
            std::sort( tmp_row, tmp_row+len, CoeffIndexCompare()) ;
            C.row(ii).grow(len+1);
            for(index_t i = 0; i < len ;i++){
                C.row(ii).add( tmp_row[i].index, tmp_row[i].a, &SparseRowColumn::coeff_back_search);
                look_up[tmp_row[i].index] = 0; 
            }
        }else{
            if (C.rows_are_stored()){
                std::sort( tmp_row, tmp_row+len, CoeffIndexCompare()) ;
            }
            for(index_t i = 0; i < len ;i++){
                C.add( ii, tmp_row[i].index, tmp_row[i].a, SparseMatrix::BACKWARD_LINEAR_SEARCH);
                look_up[tmp_row[i].index] = 0; 
            }
        }
        len = 0;
    }
    
    Allocator<Coeff>::deallocate(tmp_row) ;
    delete [] look_up;
}


void mult_C_C(SparseRowColumn& (SparseMatrix::*rowcol_a)(index_t), SparseRowColumn& (SparseMatrix::*rowcol_b)(index_t), SparseMatrix& A, SparseMatrix& B, SparseMatrix& C){
    index_t nb_row = C.m();
    index_t nb_col = C.n();
    
    double val;
    index_t ii;
    
    size_t nb_coeffs_Ca, nb_coeffs_Cb;
    index_t ka, kb;
    
    index_t* look_up = new index_t[nb_row];
    for(index_t li = 0; li<nb_row ; li++ ){
        look_up[li] = 0;
    };
    index_t icoeff;
    
    index_t maxnnz = 8 + nb_row/8 ;
    Coeff* tmp_col = Allocator<Coeff>::allocate(maxnnz) ;
    index_t len = 0;
    
    for(index_t jj = 0; jj < nb_col; jj++ ){
        const SparseRowColumn&  RCb = (B.*rowcol_b)(jj);
        nb_coeffs_Cb = RCb.nb_coeffs();
        const Coeff* cb = &(RCb.coeff(0));
        for( kb = 0; kb < nb_coeffs_Cb ; kb++ ){
            val = cb->a;
            ii = cb->index;
            const SparseRowColumn&  RCa = (A.*rowcol_a)(ii);
            nb_coeffs_Ca = RCa.nb_coeffs();
            const Coeff* ca = &(RCa.coeff(0));
            for( ka = 0; ka < nb_coeffs_Ca ; ka++ ){
                icoeff = look_up[ca->index];
                if (icoeff == 0){
                    if ( len >= maxnnz ){
                        Allocator<Coeff>::reallocate(tmp_col, maxnnz, maxnnz*2 );
                        maxnnz = maxnnz*2;
                    }
                    tmp_col[len].index = ca->index;
                    tmp_col[len].a = val*ca->a;
                    look_up[ca->index] = len+1;
                    len++;
                } else {
                    tmp_col[icoeff-1].a += val*ca->a;
                }
                ca++;
            }
            cb++;
        }
        
        if(C.columns_are_stored() && !C.rows_are_stored() && !C.has_symmetric_storage()){
            std::sort( tmp_col, tmp_col+len, CoeffIndexCompare()) ;
            C.column(jj).grow(len+1);
            for(index_t i = 0; i < len ;i++){
                C.column(jj).add( tmp_col[i].index, tmp_col[i].a, &SparseRowColumn::coeff_back_search);
                look_up[tmp_col[i].index] = 0; 
            }
        }else{
            if (C.columns_are_stored()){
                std::sort( tmp_col, tmp_col+len, CoeffIndexCompare()) ;
            }
            for(index_t i = 0; i < len ;i++){
                C.add( tmp_col[i].index, jj, tmp_col[i].a, SparseMatrix::BACKWARD_LINEAR_SEARCH);
                look_up[tmp_col[i].index] = 0; 
            }
        }
        len = 0;
    }
    
    Allocator<Coeff>::deallocate(tmp_col) ;
    delete [] look_up;
}

/*
--- AB Multiplication ---
Karnaugh table of the multiplication type, according to the storage modes 
of the matrix :

AS (or BS) = A(or B) is symetric
AR (or BR) = A(or B) has rows stored
AC (or BC) = A(or B) has columns stored

R : multiplication will be performed with the row x row algorithm
C : multiplication will be performed with the column x column algorithm
X : case that should never happend (can thus be treated as R and/or C)
white space = case when any option will do (can be treated as R _OR_ C)


              _ _ _ _  AS
          _ _ _ _      AR
        _ _     _ _    AC
      X X X X X X X X  
    | X C C   C C C X
  | | X C   R       X 
  |   X   R R R R R X 
| |   X C   R       X
| | | X C   R       X   
|   | X C   R       X  
|     X X X X X X X X
BS BC
 BR

So we can use the following expressions to define C and R: 
    C = (!AS && !AR) || (!BS && !BR)
and 
    R = !C


Karnaugh table of the necessity to rebuild A, for the multiplication type 
    defined above, and taking profit of the cases that never happend:

1 : A needs to be rebuilt
0 : A don't need to be rebuilt.

              _ _ _ _  AS
          _ _ _ _      AR
        _ _     _ _    AC
      1 0 0 1 1 1 1 1  
    | 1 0 0 1 1 1 1 1
  | | 0 0 0 0 1 1 1 1 
  |   0 0 0 0 1 1 1 1 
| |   0 0 0 0 1 1 1 1
| | | 0 0 0 0 1 1 1 1   
|   | 0 0 0 0 1 1 1 1  
|     0 0 0 0 1 1 1 1
BS BC
 BR

So rebuild_a = AS || (!AC && !BS && !BC)


Karnaugh table of the necessity to rebuild B, for the multiplication type 
    defined above, and taking profit of the cases that never happend:

1 : B needs to be rebuilt
0 : B don't need to be rebuilt.

              _ _ _ _  AS
          _ _ _ _      AR
        _ _     _ _    AC
      1 1 0 0 0 0 0 0  
    | 0 0 0 0 0 0 0 0
  | | 0 0 0 0 0 0 0 0 
  |   1 1 0 0 0 0 0 0 
| |   1 1 1 1 1 1 1 1
| | | 1 1 1 1 1 1 1 1   
|   | 1 1 1 1 1 1 1 1  
|     1 1 1 1 1 1 1 1
BS BC
 BR

So rebuild_b = BS || (!BR && !AS && AR)

--- Other products ---

when transposed :
    At.S <= A.S
    At.R <= A.C
    At.C <= A.R
    
*/

/** C = AB */
void mult_AB( SparseMatrix* A, SparseMatrix* B, SparseMatrix* C ){
    gx_assert(C->m() == A->m());
    gx_assert(C->n() == B->n());
    gx_assert(C->nnz() == 0);
    gx_assert( A != C );
    gx_assert( A != B );
    
    
    SparseMatrix *a, *b;
    SparseMatrix::Storage st;

    bool del_a = false, del_b = false;
    
    bool cc = ( !A->has_symmetric_storage() && !A->rows_are_stored() )
           || ( !B->has_symmetric_storage() && !B->rows_are_stored() ) ;
    
    if (cc){
        st = SparseMatrix::COLUMNS;
    } else {
        st = SparseMatrix::ROWS;
    }
    
    if ( A->has_symmetric_storage() 
        || ( !A->columns_are_stored()
            && !B->has_symmetric_storage() 
            && !B->rows_are_stored() ) ){
        fprintf(stderr,"Warning : AB : A matrix duplication !\n" );
        a = new SparseMatrix( A, st );
        del_a = true;
    } else {
        a = A;
    }
    
    if ( B->has_symmetric_storage() 
        || ( !B->columns_are_stored() 
            && !A->has_symmetric_storage() 
            && !A->rows_are_stored() ) ){
        fprintf(stderr,"Warning : AB : B matrix duplication !\n" );
        b = new SparseMatrix( B, st );
        del_b = true;
    } else {
        b = B;
    }
    
    if (cc){
        mult_C_C(&SparseMatrix::column, &SparseMatrix::column, *a, *b, *C);
    } else {
        mult_R_R(&SparseMatrix::row, &SparseMatrix::row, *a, *b, *C);
    }
    
    if (del_a)
        delete a;
    if (del_b)
        delete b;
}


/** C = A^t B */
void mult_AtB( SparseMatrix* A, SparseMatrix* B, SparseMatrix* C ){
    gx_assert(C->m() == A->n());
    gx_assert(C->n() == B->n());
    gx_assert(C->nnz() == 0);
    gx_assert( A != C );
    gx_assert( A != B );
    
    SparseMatrix *a, *b;
    bool del_a = false, del_b = false;
    
    bool cc = ( !A->has_symmetric_storage() && !A->columns_are_stored() )
           || ( !B->has_symmetric_storage() && !B->rows_are_stored() ) ;
    
    if ( A->has_symmetric_storage() 
        || ( !A->rows_are_stored()
            && !B->has_symmetric_storage() 
            && !B->rows_are_stored() ) ){
        fprintf(stderr,"Warning : AtB : A matrix duplication !\n" );
        del_a = true;
        if (cc){
            a = new SparseMatrix( A, SparseMatrix::ROWS );
        } else {
            a = new SparseMatrix( A, SparseMatrix::COLUMNS );
        }
    } else {
        a = A;
    }
    
    if ( B->has_symmetric_storage() 
        || ( !B->columns_are_stored() 
            && !A->has_symmetric_storage() 
            && !A->columns_are_stored() ) ){
        fprintf(stderr,"Warning : AtB : B matrix duplication !\n" );
        if (cc){
            b = new SparseMatrix( B, SparseMatrix::COLUMNS );
        } else {
            b = new SparseMatrix( B, SparseMatrix::ROWS );
        }
        del_b = true;
	} else {
        b = B;
	}
    
    if (cc){
        mult_C_C(&SparseMatrix::row, &SparseMatrix::column, *a, *b, *C);
    } else {
        mult_R_R(&SparseMatrix::column, &SparseMatrix::row, *a, *b, *C);
    }
    
    if (del_a)
        delete a;
    if (del_b)
        delete b;
}


/** C = AB^t */
void mult_ABt( SparseMatrix* A, SparseMatrix* B, SparseMatrix* C ){
    gx_assert(C->m() == A->m());
    gx_assert(C->n() == B->m());
    gx_assert(C->nnz() == 0);
    gx_assert( A != C );
    gx_assert( A != B );
    
    SparseMatrix *a, *b;
    bool del_a = false, del_b = false;
    
    bool cc = ( !A->has_symmetric_storage() && !A->rows_are_stored() )
           || ( !B->has_symmetric_storage() && !B->columns_are_stored() ) ;
    
    if ( A->has_symmetric_storage() 
        || ( !A->columns_are_stored()
            && !B->has_symmetric_storage() 
            && !B->columns_are_stored() ) ){
        fprintf(stderr,"Warning : ABt : A matrix duplication !\n" );
        if (cc){
            a = new SparseMatrix( A, SparseMatrix::COLUMNS );
        } else {
            a = new SparseMatrix( A, SparseMatrix::ROWS );
        }
        del_a = true;
    } else {
        a = A;
    }
    
    if ( B->has_symmetric_storage() 
        || ( !B->rows_are_stored() 
            && !A->has_symmetric_storage() 
            && !A->rows_are_stored() ) ){
        fprintf(stderr,"Warning : ABt : B matrix duplication !\n" );
        if (cc){
            b = new SparseMatrix( B, SparseMatrix::ROWS );
        } else {
            b = new SparseMatrix( B, SparseMatrix::COLUMNS );
        }
        del_b = true;
    } else {
        b = B;
    }
    
    if (cc){
        mult_C_C(&SparseMatrix::column, &SparseMatrix::row, *a, *b, *C);
    } else {
        mult_R_R(&SparseMatrix::row, &SparseMatrix::column, *a, *b, *C);
    }
    
    if (del_a)
        delete a;
    if (del_b)
        delete b;
}


/** C = A^tB^t */
void mult_AtBt( SparseMatrix* A, SparseMatrix* B, SparseMatrix* C ){
    gx_assert(C->m() == A->n());
    gx_assert(C->n() == B->m());
    gx_assert(C->nnz() == 0);
    gx_assert( A != C );
    gx_assert( A != B );
    
    SparseMatrix *a, *b;
    SparseMatrix::Storage st;
    bool del_a = false, del_b = false;
    
    bool cc = ( !A->has_symmetric_storage() && !A->columns_are_stored() )
           || ( !B->has_symmetric_storage() && !B->columns_are_stored() ) ;
    
    if (cc){
        st = SparseMatrix::ROWS;
    } else {
        st = SparseMatrix::COLUMNS;
    }
    
    if ( A->has_symmetric_storage() 
        || ( !A->rows_are_stored()
            && !B->has_symmetric_storage() 
            && !B->columns_are_stored() ) ){
        fprintf(stderr,"Warning : AtBt : A matrix duplication !\n" );
        a = new SparseMatrix( A, st );
        del_a = true;
    } else {
        a = A;
    }
    
    if ( B->has_symmetric_storage() 
        || ( !B->rows_are_stored() 
            && !A->has_symmetric_storage() 
            && !A->columns_are_stored() ) ){
        fprintf(stderr,"Warning : AtBt : B matrix duplication !\n" );
        b = new SparseMatrix( B, st );
        del_b = true;
    } else {
        b = B;
    }
    
    if (cc){
        mult_C_C(&SparseMatrix::row, &SparseMatrix::row, *a, *b, *C);
    } else {
        mult_R_R(&SparseMatrix::column, &SparseMatrix::column, *a, *b, *C);
    }
    
    if (del_a)
        delete a;
    if (del_b)
        delete b;
}
    


}

