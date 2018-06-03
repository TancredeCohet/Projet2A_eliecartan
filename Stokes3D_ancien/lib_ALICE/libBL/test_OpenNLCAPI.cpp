#include "OPEN_NL_C_API.h"

#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#define PRECISION 0.0000001
#define MAXMATRIX 1000


#define TEST_PROFILE_BEGIN(TEST_OBJECT_NAME)    {  \
    struct timeval TEST_OBJECT_NAME##_t0, TEST_OBJECT_NAME##_t1; \
    struct timezone TEST_OBJECT_NAME##_tz;                      \
	gettimeofday(&TEST_OBJECT_NAME##_t0, &TEST_OBJECT_NAME##_tz);
	
	/* truc a timer*/
	
	
#define TEST_PROFILE_END_BIAS(TEST_OBJECT_NAME,BIAS)                                  \
	gettimeofday(&TEST_OBJECT_NAME##_t1, &TEST_OBJECT_NAME##_tz);           \
	size_t TEST_OBJECT_NAME##_diff=(TEST_OBJECT_NAME##_t1.tv_sec-TEST_OBJECT_NAME##_t0.tv_sec) * 1000000L + (TEST_OBJECT_NAME##_t1.tv_usec-TEST_OBJECT_NAME##_t0.tv_usec); \
	printf ("TEST " #TEST_OBJECT_NAME " : %.7f s\n",( (double)TEST_OBJECT_NAME##_diff/1000000L )- BIAS );           \
	}
	
#define TEST_PROFILE_END(TEST_OBJECT_NAME) TEST_PROFILE_END_BIAS(TEST_OBJECT_NAME,0.0)

	
	
int main(){
    
    
    
#define MATRIX_I 4
#define MATRIX_J 5

    
    int m1= NewSparseMatrix(MATRIX_I, MATRIX_J);
    assert(m1==0);
    
    for (int i = 0; i < MATRIX_I; i++ ){
        for(int j = 0; j < MATRIX_J; j++){
            AddSparseMatrix(m1, i, j,  0.34*i+5.1*j);
        }
    }
    assert(coefficient(m1,0,0)<PRECISION);
    assert(coefficient(m1,1,0)-0.34<PRECISION);
    assert(coefficient(m1,2,0)-0.68<PRECISION);
    assert(coefficient(m1,3,0)-1.02<PRECISION);
    assert(coefficient(m1,0,1)-5.1<PRECISION);
    assert(coefficient(m1,1,1)-5.44<PRECISION);
    assert(coefficient(m1,2,1)-5.78<PRECISION);
    assert(coefficient(m1,3,1)-6.12<PRECISION);
    assert(coefficient(m1,0,2)-10.2<PRECISION);
    assert(coefficient(m1,1,2)-10.54<PRECISION);
    assert(coefficient(m1,2,2)-10.88<PRECISION);
    assert(coefficient(m1,3,2)-11.22<PRECISION);
    assert(coefficient(m1,0,3)-15.3<PRECISION);
    assert(coefficient(m1,1,3)-15.64<PRECISION);
    assert(coefficient(m1,2,3)-15.98<PRECISION);
    assert(coefficient(m1,3,3)-16.32<PRECISION);
    assert(coefficient(m1,0,4)-20.4<PRECISION);
    assert(coefficient(m1,1,4)-20.74<PRECISION);
    assert(coefficient(m1,2,4)-21.08<PRECISION);
    assert(coefficient(m1,3,4)-21.42<PRECISION);
    
    
    /*
    printf("Initial Matrix : \n");
    for (int i = 0; i < MATRIX_I; i++ ){
        for(int j = 0; j < MATRIX_J; j++){
            printf("%f\t", coefficient(m1, i, j));
        }
        printf("\n");
    }
    */
    assert(NnzSparseMatrix(m1)==20);
    assert(SparseMatrixM(m1)==4);
    assert(SparseMatrixN(m1)==5);
    
    //printf("NNZ = %d; M = %d, N = %d\n", NnzSparseMatrix(m1), SparseMatrixM(m1), SparseMatrixN(m1));
    
    #define NB_INDEX 2
    gx_index_T indexes[NB_INDEX];
    indexes[0] = 1;
    indexes[1] = 3;
    
    RemoveColsSparseMatrix(m1, NB_INDEX, indexes); 
    RemoveRowsSparseMatrix(m1, NB_INDEX, indexes); 
    for(int i = 0; i< NB_INDEX; i++){
        AddSparseMatrix(m1, indexes[i],indexes[i],1);
    }
    
    /*
    printf("After removing some lines... : \n nnz = %zd\n",NnzSparseMatrix(m1));
    for (int i = 0; i < MATRIX_I; i++ ){
        for(int j = 0; j < MATRIX_J; j++){
            printf("%f\t", coefficient(m1, i, j));
        }
        printf("\n");
    }
    */
    
    assert(coefficient(m1,0,0)<PRECISION);
    assert(coefficient(m1,2,0)-0.68<PRECISION);
    assert(coefficient(m1,0,2)-10.2<PRECISION);
    assert(coefficient(m1,2,2)-10.88<PRECISION);
    assert(coefficient(m1,0,4)-20.4<PRECISION);
    assert(coefficient(m1,2,4)-21.08<PRECISION);
    
    assert(coefficient(m1,1,1)-1.0<PRECISION);
    assert(coefficient(m1,3,3)-1.0<PRECISION);
    assert(coefficient(m1,1,0)<PRECISION);
    assert(coefficient(m1,1,2)<PRECISION);
    assert(coefficient(m1,1,3)<PRECISION);
    assert(coefficient(m1,1,4)<PRECISION);
    assert(coefficient(m1,0,1)<PRECISION);
    assert(coefficient(m1,2,1)<PRECISION);
    assert(coefficient(m1,3,1)<PRECISION);
    assert(coefficient(m1,3,0)<PRECISION);
    assert(coefficient(m1,3,2)<PRECISION);
    assert(coefficient(m1,3,4)<PRECISION);
    assert(coefficient(m1,0,3)<PRECISION);
    assert(coefficient(m1,2,3)<PRECISION);
    
    
    
    
    assert(NnzSparseMatrix(m1)==8);
    assert(SparseMatrixM(m1)==4);
    assert(SparseMatrixN(m1)==5);
    //printf("NNZ = %d; M = %d, N = %d\n", NnzSparseMatrix(m1), SparseMatrixM(m1), SparseMatrixN(m1));
    
    
    
    
    //printf("Suppression de toutes les matrices !\n");
    ClearSparseMatrix();
    assert(SparseMatrixExists(m1)==false);


    int m;
    for (int i = 0; i < MAXMATRIX; i++ ){
        m = NewSparseMatrix(MATRIX_I, MATRIX_J);
        assert(m!= -1);
    }
    assert(m = MAXMATRIX-1);
    
    ClearSparseMatrix();
    
    /*
    
    int m2= NewSparseMatrix(MATRIX_I, MATRIX_J);
    
    double vectorX[MATRIX_J];
    double vectorY[MATRIX_I];
    
    assert(m2==0);
    
    for (int i = 0; i < MATRIX_I; i++ ){
        for(int j = 0; j < MATRIX_J; j++){
            AddSparseMatrix(m1, i, j,  i);
        }
    }
    
    for (int i = 0; i < MATRIX_I; i++ ){
        vectorY[i] = 0.0;
    }
    for (int j = 0; j < MATRIX_J; j++){
        vectorX[j] = j;
    }
    
       
    printf("Initial Matrix : \n");
    for (int i = 0; i < MATRIX_I; i++ ){
        for(int j = 0; j < MATRIX_J; j++){
            printf("%f\t", coefficient(m1, i, j));
        }
        printf("\n");
    }
    
    printf("Vector Y : \n");
    for (int i = 0; i < MATRIX_I; i++ ){
        printf("%f\n",vectorY[i]);
    }
    
    printf("Vector X : \n");
    for (int j = 0; j < MATRIX_J; j++){
        printf("%f\n",vectorX[j]);
    }
    
    VectorMultSparseMatrix(m2, vectorX, vectorY );
    
    VectorMultTransposeSparseMatrix(m2, vectorY, vectorX );
    
    
    printf("Vector Y  AFTER: \n");
    for (int i = 0; i < MATRIX_I; i++ ){
        printf("%f\n",vectorY[i]);
    }
    
    printf("Vector X  AFTER : \n");
    for (int j = 0; j < MATRIX_J; j++){
        printf("%f\n",vectorX[j]);
    }
    */
    
    
    

#undef MATRIX_I
#undef MATRIX_J




    
    
#ifdef USE_TEST_PROFILE

#define MATRIX_I 30
#define MATRIX_J 30
    
    
    TEST_PROFILE_BEGIN(construction)
    for(int wxd = 0; wxd < 10000; wxd++ ){
        int m1= NewSparseMatrix(MATRIX_I, MATRIX_J);
        for (int i = 0; i < MATRIX_I; i++ ){
            for(int j = 0; j < MATRIX_J; j++){
                AddSparseMatrix(m1, i, j,  0.34*i+5.1*j);
            }
        }
    
        gx_index_T indexes[10];
        indexes[0] = 6;
        indexes[1] = 12;   
        indexes[2] = 15;   
        indexes[3] = 18;   
        indexes[4] = 22;   
        indexes[5] = 24;   
        indexes[6] = 25;   
        indexes[7] = 26;   
        indexes[8] = 27;   
        indexes[9] = 28;   
        
       ClearSparseMatrix();
    }
    TEST_PROFILE_END(construction)
    
    TEST_PROFILE_BEGIN(construction_et_remov)
    for(int wxd = 0; wxd < 10000; wxd++ ){
        int m1= NewSparseMatrix(MATRIX_I, MATRIX_J);
        for (int i = 0; i < MATRIX_I; i++ ){
            for(int j = 0; j < MATRIX_J; j++){
                AddSparseMatrix(m1, i, j,  0.34*i+5.1*j);
            }
        }
        
        #define NB_INDEX 10
        gx_index_T indexes[NB_INDEX];
        indexes[0] = 6;
        indexes[1] = 12;   
        indexes[2] = 15;   
        indexes[3] = 18;   
        indexes[4] = 22;   
        indexes[5] = 24;   
        indexes[6] = 25;   
        indexes[7] = 26;   
        indexes[8] = 27;   
        indexes[9] = 28;   
        
        RemoveColsSparseMatrix(m1, NB_INDEX, indexes); 
        RemoveRowsSparseMatrix(m1, NB_INDEX, indexes); 
        for(int i = 0; i< NB_INDEX; i++){
            AddSparseMatrix(m1, indexes[i],indexes[i],1);
        }
        
        ClearSparseMatrix();
    
    }
    TEST_PROFILE_END(construction_et_remov)
    /*
    TEST_PROFILE_BEGIN(construction_et_remov_old)
    for(int wxd = 0; wxd < 10000; wxd++ ){
        int m1= NewSparseMatrix(MATRIX_I, MATRIX_J);
        for (int i = 0; i < MATRIX_I; i++ ){
            for(int j = 0; j < MATRIX_J; j++){
                AddSparseMatrix(m1, i, j,  0.34*i+5.1*j);
            }
        }
    
        gx_index_T indexes[10];
        indexes[0] = 6;
        indexes[1] = 12;   
        indexes[2] = 15;   
        indexes[3] = 18;   
        indexes[4] = 22;   
        indexes[5] = 24;   
        indexes[6] = 25;   
        indexes[7] = 26;   
        indexes[8] = 27;    
        indexes[9] = 28;   
        
        RemoveRowsColsSparseMatrix_old(m1, 10, indexes); 
        
        ClearSparseMatrix();
    
    }
    TEST_PROFILE_END(construction_et_remov_old)   */ 
#endif


}


