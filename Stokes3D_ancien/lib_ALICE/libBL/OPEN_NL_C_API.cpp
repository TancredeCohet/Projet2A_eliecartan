#include "OPEN_NL_C_API.h"

#include <stdio.h>
#include <stdlib.h>
#include "sparse_matrix.h"
#include "geex_defs.h"


class MatrixTable {
public:
    MatrixTable(unsigned int initialSize){
        mMatrixTable = new pSparseMatrix[initialSize];
        size=initialSize;
        count = 0;
        for (unsigned int i = 0; i<size; i++){
            mMatrixTable[i] = NULL;
        }
    }
    
    ~MatrixTable(){
        clear();
        delete [] mMatrixTable;
    }
    
    MatrixID_T add(pSparseMatrix m){
        // If table is full, double its size. Could be optimized...
        if (count > size-1){
            pSparseMatrix* tmp = new pSparseMatrix[2*size] ;
            for(unsigned int i=0; i<size; i++) {
                tmp[i] = mMatrixTable[i] ;
            }
            for(unsigned int i=size; i<2*size; i++) {
                tmp[i] = NULL ;
            }
            delete[] mMatrixTable ;
            mMatrixTable = tmp ;
            size = size*2;
        } 
    
        for (unsigned int i = 0; i<size; i++){
            if (mMatrixTable[i] == NULL){
                mMatrixTable[i] = m;
                count++;
                return i;
            }
        }
        return -1;
    }
    
    void del(MatrixID_T i){
        delete mMatrixTable[i];
        mMatrixTable[i] = NULL;
    }
    
    pSparseMatrix operator [] (MatrixID_T i) const{
        if (i>=0 && (unsigned int)i<size){
            return mMatrixTable[i];
        }
        return NULL;
    }
    
    void clear(){
        for (unsigned int ind=0; ind<size;ind++){
            if (mMatrixTable[ind] != NULL){
                delete mMatrixTable[ind];
                mMatrixTable[ind] = NULL;
            }
        }
        count = 0;
    }
    
private:
    pSparseMatrix* mMatrixTable;
    unsigned int size;
    unsigned int count;
};


/********************************************/

static MatrixTable GlobalMatrixTable(50);

MatrixID_T NewSparseMatrix(gx_size_T m,gx_size_T n){
    assert(m>0);
    assert(n>0);

    return GlobalMatrixTable.add(new Geex::SparseMatrix(m,n,Geex::SparseMatrix::COLUMNS));
}

MatrixID_T CopySparseMatrix(MatrixID_T MatrixID){
    pSparseMatrix p = GlobalMatrixTable[MatrixID];
    assert(p!=NULL);
    return GlobalMatrixTable.add(new Geex::SparseMatrix(p,Geex::SparseMatrix::COLUMNS));
}

void DeleteSparseMatrix(MatrixID_T MatrixID){
    assert(GlobalMatrixTable[MatrixID]!=NULL);
    
    GlobalMatrixTable.del(MatrixID);
}


void ClearSparseMatrix(){
    GlobalMatrixTable.clear();
}

void AddSparseMatrix(MatrixID_T MatrixID,gx_index_T i,gx_index_T j,double Valeur){
    assert(GlobalMatrixTable[MatrixID]!=NULL);
    assert(i>=0);
    assert(j>=0);
    
    if (Valeur != 0){
        GlobalMatrixTable[MatrixID]->add(i,j,Valeur);
    }
}

gx_size_T NnzSparseMatrix(MatrixID_T MatrixID){
    assert(GlobalMatrixTable[MatrixID]!=NULL);
    
    return GlobalMatrixTable[MatrixID]->nnz();
}

double coefficient(MatrixID_T MatrixID, gx_index_T i, gx_index_T j){
    assert(GlobalMatrixTable[MatrixID]!=NULL);
    assert(i>=0);
    assert(j>=0);
    
    return GlobalMatrixTable[MatrixID]->value(i,j);
}

gx_size_T SparseMatrixN(MatrixID_T MatrixID){
    assert(GlobalMatrixTable[MatrixID]!=NULL);
    
    return GlobalMatrixTable[MatrixID]->n();
}

gx_size_T SparseMatrixM(MatrixID_T MatrixID){
    assert(GlobalMatrixTable[MatrixID]!=NULL);
    
    return GlobalMatrixTable[MatrixID]->m();
}

void SparseMatrixToMatlab(MatrixID_T MatrixID, size_t* irs, size_t* jcs, double* sr){
    assert(GlobalMatrixTable[MatrixID]!=NULL);
    
    size_t k=0;
    for (size_t j = 0; j < GlobalMatrixTable[MatrixID]->n(); j++) {

        jcs[j] = k;

        GlobalMatrixTable[MatrixID]->column(j).dump_indexes_values(irs+k, sr+k);
        
        k += GlobalMatrixTable[MatrixID]->column(j).nb_coeffs();
    }

    jcs[GlobalMatrixTable[MatrixID]->n()] = k;
}

void RemoveRowsSparseMatrix(MatrixID_T MatrixID, gx_size_T nbIndexes, gx_index_T* indexes){
    assert(GlobalMatrixTable[MatrixID]!=NULL);
    assert(nbIndexes>0);
    
    GlobalMatrixTable[MatrixID]->zero_rows(nbIndexes,indexes);
}

void RemoveColsSparseMatrix(MatrixID_T MatrixID, gx_size_T nbIndexes, gx_index_T* indexes){
    assert(GlobalMatrixTable[MatrixID]!=NULL);
    assert(nbIndexes>0);
    
    GlobalMatrixTable[MatrixID]->zero_columns(nbIndexes,indexes);
}


bool SparseMatrixExists(MatrixID_T MatrixID){
    return GlobalMatrixTable[MatrixID]!=NULL;
}


void VectorMultSparseMatrix(MatrixID_T MatrixID, const double* x, double* y ){
    assert(GlobalMatrixTable[MatrixID]!=NULL);
    
    mult(*GlobalMatrixTable[MatrixID], x, y);
}
void VectorMultTransposeSparseMatrix(MatrixID_T MatrixID, const double* x, double* y ){
    assert(GlobalMatrixTable[MatrixID]!=NULL);
    
    mult_transpose(*GlobalMatrixTable[MatrixID], x, y);
}

// C += AB
void MultABSparseMatrix(MatrixID_T A_ID, MatrixID_T B_ID, MatrixID_T C_ID ){
    assert(GlobalMatrixTable[A_ID]!=NULL);
    assert(GlobalMatrixTable[B_ID]!=NULL);
    assert(GlobalMatrixTable[C_ID]!=NULL);
    
    mult_AB(GlobalMatrixTable[A_ID],GlobalMatrixTable[B_ID],GlobalMatrixTable[C_ID]);
}

// C += At B
void MultAtBSparseMatrix(MatrixID_T A_ID, MatrixID_T B_ID, MatrixID_T C_ID ){
    assert(GlobalMatrixTable[A_ID]!=NULL);
    assert(GlobalMatrixTable[B_ID]!=NULL);
    assert(GlobalMatrixTable[C_ID]!=NULL);
    
    mult_AtB(GlobalMatrixTable[A_ID],GlobalMatrixTable[B_ID],GlobalMatrixTable[C_ID]);
}

// C += A Bt
void MultABtSparseMatrix(MatrixID_T A_ID, MatrixID_T B_ID, MatrixID_T C_ID ){
    assert(GlobalMatrixTable[A_ID]!=NULL);
    assert(GlobalMatrixTable[B_ID]!=NULL);
    assert(GlobalMatrixTable[C_ID]!=NULL);
    
    mult_ABt(GlobalMatrixTable[A_ID],GlobalMatrixTable[B_ID],GlobalMatrixTable[C_ID]);
}

// C += At Bt
void MultAtBtSparseMatrix(MatrixID_T A_ID, MatrixID_T B_ID, MatrixID_T C_ID ){
    assert(GlobalMatrixTable[A_ID]!=NULL);
    assert(GlobalMatrixTable[B_ID]!=NULL);
    assert(GlobalMatrixTable[C_ID]!=NULL);
    
    mult_AtBt(GlobalMatrixTable[A_ID],GlobalMatrixTable[B_ID],GlobalMatrixTable[C_ID]);
}










