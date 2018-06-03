//matrice P2 entierement en MEXfunction
//[DxDx,DyDy,DxDy,DxW,DyW,SS]=MatriceP2_cpp(v1,t1,ne1)

//ATTENTION : il faut passer v1 et t1 en ligne (ils doivent etre transpos�)

#include "mex.h"
#include "libBL.h"

#include "math.h"



/****************** LINUX TEST  *******************/
#ifdef _WIN32
//WINDOWS

#define TEST_PROFILE_BEGIN(TEST_OBJECT_NAME) 
#define TEST_PROFILE_END_BIAS(TEST_OBJECT_NAME,BIAS)
#define TEST_PROFILE_END(TEST_OBJECT_NAME) TEST_PROFILE_END_BIAS(TEST_OBJECT_NAME,0.0)

#else
// LINUX

#include "sys/time.h"

#define TEST_PROFILE_BEGIN(TEST_OBJECT_NAME)    {  \
    struct timeval TEST_OBJECT_NAME##_t0, TEST_OBJECT_NAME##_t1; \
    struct timezone TEST_OBJECT_NAME##_tz;                      \
	gettimeofday(&TEST_OBJECT_NAME##_t0, &TEST_OBJECT_NAME##_tz);
		
	
#define TEST_PROFILE_END_BIAS(TEST_OBJECT_NAME,BIAS)                                  \
	gettimeofday(&TEST_OBJECT_NAME##_t1, &TEST_OBJECT_NAME##_tz);           \
	size_t TEST_OBJECT_NAME##_diff=(TEST_OBJECT_NAME##_t1.tv_sec-TEST_OBJECT_NAME##_t0.tv_sec) * 1000000L + (TEST_OBJECT_NAME##_t1.tv_usec-TEST_OBJECT_NAME##_t0.tv_usec); \
	printf ("TEST " #TEST_OBJECT_NAME " : %.7f s\n",( (double)TEST_OBJECT_NAME##_diff/1000000L )- BIAS );           \
	}
	
#define TEST_PROFILE_END(TEST_OBJECT_NAME) TEST_PROFILE_END_BIAS(TEST_OBJECT_NAME,0.0)
#endif
/*************************************************/



void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    double *v1, *t1;
    
    double SS_elem[6][6], DxDx_elem[6][6], DyDy_elem[6][6], DxDy_elem[6][6], DxW_elem[6][3], DyW_elem[6][3];
    
    size_t ne1, nt1, nv1, nv2, ind1, ind2, ind_triangle;
    
    size_t i1[3], i2[6];
    
    double p[3][2], aire, val, signe;
    
    size_t SS, DxDx, DxDy, DyDy, DxW, DyW;
    
    double *sr1, *sr2, *sr3, *sr4, *sr5, *sr6;
    
    mwIndex *irs1, *jcs1,*irs2, *jcs2,*irs3, *jcs3,*irs4, *jcs4,*irs5, *jcs5,*irs6, *jcs6;
    
    size_t m,n,nzmax;
    
    v1=mxGetPr(prhs[0]);
    
    t1=mxGetPr(prhs[1]);
    
    ne1=(size_t)mxGetScalar(prhs[2]);
    
    nt1=mxGetN(prhs[1]);
    
    nv1=mxGetN(prhs[0]);
    
    nv2=nv1+ne1;
    
    
TEST_PROFILE_BEGIN(CreationSparseMatrix)

    SS=NewSparseMatrix(nv2,nv2);
    
    DxDx=NewSparseMatrix(nv2,nv2);
    
    DyDy=NewSparseMatrix(nv2,nv2);
    
    DxDy=NewSparseMatrix(nv2,nv2);
    
    DxW=NewSparseMatrix(nv2,nv1);
    
    DyW=NewSparseMatrix(nv2,nv1);
    
TEST_PROFILE_END(CreationSparseMatrix)
    
    
    //remplissage de la matrice elementaire SS_elem

    SS_elem[0][0]=6./180;
    SS_elem[0][1]=-1./180;
    SS_elem[0][2]=-1./180;
    SS_elem[0][3]=-4./180;
    SS_elem[0][4]=0./180;
    SS_elem[0][5]=0./180;
    SS_elem[1][0]=-1./180;
    SS_elem[1][1]=6./180;
    SS_elem[1][2]=-1./180;
    SS_elem[1][3]=0./180;
    SS_elem[1][4]=-4./180;
    SS_elem[1][5]=0./180;
    SS_elem[2][0]=-1./180;
    SS_elem[2][1]=-1./180;
    SS_elem[2][2]=6./180;
    SS_elem[2][3]=0./180;
    SS_elem[2][4]=0./180;
    SS_elem[2][5]=-4./180;
    SS_elem[3][0]=-4./180;
    SS_elem[3][1]=0./180;
    SS_elem[3][2]=0./180;
    SS_elem[3][3]=32./180;
    SS_elem[3][4]=16./180;
    SS_elem[3][5]=16./180;
    SS_elem[4][0]=0./180;
    SS_elem[4][1]=-4./180;
    SS_elem[4][2]=0./180;
    SS_elem[4][3]=16./180;
    SS_elem[4][4]=32./180;
    SS_elem[4][5]=16./180;
    SS_elem[5][0]=0./180;
    SS_elem[5][1]=0./180;
    SS_elem[5][2]=-4./180;
    SS_elem[5][3]=16./180;
    SS_elem[5][4]=16./180;
    SS_elem[5][5]=32./180;

    //remplissage de la amtrice DxDx_elem

    DxDx_elem[0][0]=3./6;
    DxDx_elem[0][1]=1./6;
    DxDx_elem[0][2]=0./6;
    DxDx_elem[0][3]=0./6;
    DxDx_elem[0][4]=0./6;
    DxDx_elem[0][5]=-4./6;
    DxDx_elem[1][0]=1./6;
    DxDx_elem[1][1]=3./6;
    DxDx_elem[1][2]=0./6;
    DxDx_elem[1][3]=0./6;
    DxDx_elem[1][4]=0./6;
    DxDx_elem[1][5]=-4./6;
    DxDx_elem[2][0]=0./6;
    DxDx_elem[2][1]=0./6;
    DxDx_elem[2][2]=0./6;
    DxDx_elem[2][3]=0./6;
    DxDx_elem[2][4]=0./6;
    DxDx_elem[2][5]=0./6;
    DxDx_elem[3][0]=0./6;
    DxDx_elem[3][1]=0./6;
    DxDx_elem[3][2]=0./6;
    DxDx_elem[3][3]=8./6;
    DxDx_elem[3][4]=-8./6;
    DxDx_elem[3][5]=0./6;
    DxDx_elem[4][0]=0./6;
    DxDx_elem[4][1]=0./6;
    DxDx_elem[4][2]=0./6;
    DxDx_elem[4][3]=-8./6;
    DxDx_elem[4][4]=8./6;
    DxDx_elem[4][5]=0./6;
    DxDx_elem[5][0]=-4./6;
    DxDx_elem[5][1]=-4./6;
    DxDx_elem[5][2]=0./6;
    DxDx_elem[5][3]=0./6;
    DxDx_elem[5][4]=0./6;
    DxDx_elem[5][5]=8./6;
    
    //remplssage de DyDy_elem
    
    DyDy_elem[0][0]=3./6;
    DyDy_elem[0][1]=0./6;
    DyDy_elem[0][2]=1./6;
    DyDy_elem[0][3]=0./6;
    DyDy_elem[0][4]=-4./6;
    DyDy_elem[0][5]=0./6;
    DyDy_elem[1][0]=0./6;
    DyDy_elem[1][1]=0./6;
    DyDy_elem[1][2]=0./6;
    DyDy_elem[1][3]=0./6;
    DyDy_elem[1][4]=0./6;
    DyDy_elem[1][5]=0./6;
    DyDy_elem[2][0]=1./6;
    DyDy_elem[2][1]=0./6;
    DyDy_elem[2][2]=3./6;
    DyDy_elem[2][3]=0./6;
    DyDy_elem[2][4]=-4./6;
    DyDy_elem[2][5]=0./6;
    DyDy_elem[3][0]=0./6;
    DyDy_elem[3][1]=0./6;
    DyDy_elem[3][2]=0./6;
    DyDy_elem[3][3]=8./6;
    DyDy_elem[3][4]=0./6;
    DyDy_elem[3][5]=-8./6;
    DyDy_elem[4][0]=-4./6;
    DyDy_elem[4][1]=0./6;
    DyDy_elem[4][2]=-4./6;
    DyDy_elem[4][3]=0./6;
    DyDy_elem[4][4]=8./6;
    DyDy_elem[4][5]=0./6;
    DyDy_elem[5][0]=0./6;
    DyDy_elem[5][1]=0./6;
    DyDy_elem[5][2]=0./6;
    DyDy_elem[5][3]=-8./6;
    DyDy_elem[5][4]=0./6;
    DyDy_elem[5][5]=8./6;
    
    //remplissage de DxDy_elem
    
    DxDy_elem[0][0]=3./6;
    DxDy_elem[0][1]=0./6;
    DxDy_elem[0][2]=1./6;
    DxDy_elem[0][3]=0./6;
    DxDy_elem[0][4]=-4./6;
    DxDy_elem[0][5]=0./6;
    DxDy_elem[1][0]=1./6;
    DxDy_elem[1][1]=0./6;
    DxDy_elem[1][2]=-1./6;
    DxDy_elem[1][3]=4./6;
    DxDy_elem[1][4]=0./6;
    DxDy_elem[1][5]=-4./6;
    DxDy_elem[2][0]=0./6;
    DxDy_elem[2][1]=0./6;
    DxDy_elem[2][2]=0./6;
    DxDy_elem[2][3]=0./6;
    DxDy_elem[2][4]=0./6;
    DxDy_elem[2][5]=0./6;
    DxDy_elem[3][0]=0./6;
    DxDy_elem[3][1]=0./6;
    DxDy_elem[3][2]=4./6;
    DxDy_elem[3][3]=4./6;
    DxDy_elem[3][4]=-4./6;
    DxDy_elem[3][5]=-4./6;
    DxDy_elem[4][0]=0./6;
    DxDy_elem[4][1]=0./6;
    DxDy_elem[4][2]=-4./6;
    DxDy_elem[4][3]=-4./6;
    DxDy_elem[4][4]=4./6;
    DxDy_elem[4][5]=4./6;
    DxDy_elem[5][0]=-4./6;
    DxDy_elem[5][1]=0./6;
    DxDy_elem[5][2]=0./6;
    DxDy_elem[5][3]=-4./6;
    DxDy_elem[5][4]=4./6;
    DxDy_elem[5][5]=4./6;

    //remplissage de DxW_elem
    
    DxW_elem[0][0]=-1./6;
    DxW_elem[0][1]=0./6;
    DxW_elem[0][2]=0./6;
    DxW_elem[1][0]=0./6;
    DxW_elem[1][1]=1./6;
    DxW_elem[1][2]=0./6;
    DxW_elem[2][0]=0./6;
    DxW_elem[2][1]=0./6;
    DxW_elem[2][2]=0./6;
    DxW_elem[3][0]=1./6;
    DxW_elem[3][1]=1./6;
    DxW_elem[3][2]=2./6;
    DxW_elem[4][0]=-1./6;
    DxW_elem[4][1]=-1./6;
    DxW_elem[4][2]=-2./6;
    DxW_elem[5][0]=1./6;
    DxW_elem[5][1]=-1./6;
    DxW_elem[5][2]=0./6;
    
    //remplissage de DyW_elem
    
    DyW_elem[0][0]=-1./6;
    DyW_elem[0][1]=0./6;
    DyW_elem[0][2]=0./6;
    DyW_elem[1][0]=0./6;
    DyW_elem[1][1]=0./6;
    DyW_elem[1][2]=0./6;
    DyW_elem[2][0]=0./6;
    DyW_elem[2][1]=0./6;
    DyW_elem[2][2]=1./6;
    DyW_elem[3][0]=1./6;
    DyW_elem[3][1]=2./6;
    DyW_elem[3][2]=1./6;
    DyW_elem[4][0]=1./6;
    DyW_elem[4][1]=0./6;
    DyW_elem[4][2]=-1./6;
    DyW_elem[5][0]=-1./6;
    DyW_elem[5][1]=-2./6;
    DyW_elem[5][2]=-1./6;
    
    //boucle sur les triangles
    

TEST_PROFILE_BEGIN(RemplirLesTriangles)

    for (ind_triangle=0; ind_triangle<nt1; ind_triangle++){
        
        //on commence par remplir les variables du triangle consid�r�
        
        for (ind1=0; ind1<3; ind1++){
            
            i1[ind1]=(size_t) (t1[8*ind_triangle+ind1]-1);
            
            i2[ind1]= i1[ind1];
            
            i2[3+ind1]=(size_t) (t1[8*ind_triangle+ind1+5]+nv1-1);
            
            p[ind1][0]=v1[i1[ind1]*3];
            
            
            p[ind1][1]=v1[i1[ind1]*3+1];
            
        }
        
        aire=t1[8*ind_triangle+4];
  
        
        //boucle sur les cases � remplir
        
        //on calcule le signe des valeurs � inserer dans DxW et DyW
        
        if ( ((p[1][0]-p[0][0])*(p[2][1]-p[0][1])-(p[2][0]-p[0][0])*(p[1][1]-p[0][1]))<0 ){
                        
            signe=-1;

        }else{

            signe=1;

        }
        
        
        for (ind1=0; ind1<6; ind1++){
            
            for (ind2=0; ind2<6;ind2++){
         
                if(ind2<3){
                    
                    //remplissge de DxW et DyW
                    
                    //on calcule la valeur � inserer a l'emplacement [ind1][ind2]
                    
                    val=signe*((p[2][1]-p[0][1])*DxW_elem[ind1][ind2]-(p[1][1]-p[0][1])*DyW_elem[ind1][ind2]);
                    
                    AddSparseMatrix(DxW,i2[ind1],i1[ind2],val);
                    
                    val=signe*(-(p[2][0]-p[0][0])*DxW_elem[ind1][ind2]+(p[1][0]-p[0][0])*DyW_elem[ind1][ind2]);
                    
                    AddSparseMatrix(DyW,i2[ind1],i1[ind2],val);
                    
                }
                
                //remplissage de SS, DxDx, DxDy, DyDy
                
                AddSparseMatrix(SS,i2[ind1],i2[ind2],aire*SS_elem[ind1][ind2]);
                
                val=(pow(p[2][1]-p[0][1],2))*DxDx_elem[ind1][ind2]-(p[2][1]-p[0][1])*(p[1][1]-p[0][1])*(DxDy_elem[ind1][ind2]+DxDy_elem[ind2][ind1])+(pow(p[1][1]-p[0][1],2))*DyDy_elem[ind1][ind2];
                
                AddSparseMatrix(DxDx,i2[ind1],i2[ind2],val/(2*aire));
                
                val=(pow(p[2][0]-p[0][0],2))*DxDx_elem[ind1][ind2]-(p[2][0]-p[0][0])*(p[1][0]-p[0][0])*(DxDy_elem[ind1][ind2]+DxDy_elem[ind2][ind1])+(pow(p[1][0]-p[0][0],2))*DyDy_elem[ind1][ind2];
                
                AddSparseMatrix(DyDy,i2[ind1],i2[ind2],val/(2*aire));
                
                val=-(p[2][0]-p[0][0])*(p[2][1]-p[0][1])*DxDx_elem[ind1][ind2]+(p[1][0]-p[0][0])*(p[2][1]-p[0][1])*DxDy_elem[ind1][ind2]+(p[1][1]-p[0][1])*(p[2][0]-p[0][0])*DxDy_elem[ind2][ind1]-(p[1][0]-p[0][0])*(p[1][1]-p[0][1])*DyDy_elem[ind1][ind2];
                
                AddSparseMatrix(DxDy,i2[ind1],i2[ind2],val/(2*aire));
                
            }
            
        }
        
    }
            
    
TEST_PROFILE_END(RemplirLesTriangles)
    
    
    
    //DxDx
    
TEST_PROFILE_BEGIN(DxDx)
    m=SparseMatrixM(DxDx);

    n=SparseMatrixN(DxDx);
        
    nzmax=NnzSparseMatrix(DxDx);
   
    plhs[0] = mxCreateSparse(m,n,nzmax,mxREAL);

    sr1  = mxGetPr(plhs[0]);

    irs1 = mxGetIr(plhs[0]);

    jcs1 = mxGetJc(plhs[0]);

    SparseMatrixToMatlab(DxDx, irs1, jcs1, sr1);
    
TEST_PROFILE_END(DxDx)


TEST_PROFILE_BEGIN(DyDy)
    //DyDy
    
    m=SparseMatrixM(DyDy);

    n=SparseMatrixN(DyDy);
        
    nzmax=NnzSparseMatrix(DyDy);
   
    plhs[1] = mxCreateSparse(m,n,nzmax,mxREAL);

    sr2  = mxGetPr(plhs[1]);

    irs2 = mxGetIr(plhs[1]);

    jcs2 = mxGetJc(plhs[1]);

    SparseMatrixToMatlab(DyDy, irs2, jcs2, sr2);
    
TEST_PROFILE_END(DyDy)
    
    
    //DxDy
TEST_PROFILE_BEGIN(DxDy)

    m=SparseMatrixM(DxDy);

    n=SparseMatrixN(DxDy);
        
    nzmax=NnzSparseMatrix(DxDy);
   
    plhs[2] = mxCreateSparse(m,n,nzmax,mxREAL);

    sr3  = mxGetPr(plhs[2]);

    irs3 = mxGetIr(plhs[2]);

    jcs3 = mxGetJc(plhs[2]);

    SparseMatrixToMatlab(DxDy, irs3, jcs3, sr3);
TEST_PROFILE_END(DxDy)
   
    //DxW
TEST_PROFILE_BEGIN(DxW)
    
    m=SparseMatrixM(DxW);

    n=SparseMatrixN(DxW);
        
    nzmax=NnzSparseMatrix(DxW);
   
    plhs[3] = mxCreateSparse(m,n,nzmax,mxREAL);

    sr4  = mxGetPr(plhs[3]);

    irs4 = mxGetIr(plhs[3]);

    jcs4 = mxGetJc(plhs[3]);

    SparseMatrixToMatlab(DxW, irs4, jcs4, sr4);
    
TEST_PROFILE_END(DxW)
    
    
    //DyW
    
TEST_PROFILE_BEGIN(DyW)
    m=SparseMatrixM(DyW);

    n=SparseMatrixN(DyW);
        
    nzmax=NnzSparseMatrix(DyW);
   
    plhs[4] = mxCreateSparse(m,n,nzmax,mxREAL);

    sr5  = mxGetPr(plhs[4]);

    irs5 = mxGetIr(plhs[4]);

    jcs5 = mxGetJc(plhs[4]);

    SparseMatrixToMatlab(DyW, irs5, jcs5, sr5);
    
TEST_PROFILE_END(DyW)
    
    
    
    //SS
    
TEST_PROFILE_BEGIN(SS)
    m=SparseMatrixM(SS);

    n=SparseMatrixN(SS);
        
    nzmax=NnzSparseMatrix(SS);
   
    plhs[5] = mxCreateSparse(m,n,nzmax,mxREAL);

    sr6  = mxGetPr(plhs[5]);

    irs6 = mxGetIr(plhs[5]);

    jcs6 = mxGetJc(plhs[5]);

    SparseMatrixToMatlab(SS, irs6, jcs6, sr6);

TEST_PROFILE_END(SS)


    ClearSparseMatrix();

}

