#include "mex.h"
#include <vtkCell.h>
#include <vtkTetra.h>
#include <vtkPointSet.h>
#include <vtkUnstructuredGrid.h>
#include <iostream>
#include <ostream>
#include <fstream>
#include <vtkCellLocator.h>

/*
 *  mexfile FindTetra_VTK 
 *
 *  Mars 2012 
 *
 *  Contributors:   Jeremy Sinoir - INRIA - CORIDA Project Team
 *                  Marc Fuentes - INRIA - SED 
 *
 *  MATLAB Usage : [ind, V_interp, R_interp] = FindTetra_VTK(P, F, T, V, R)
 *
 *       ENTREE : 
 *       --------
 * 	 P : nb_P x 3 les coordonnées des points du maillage volumique de l'objet
 *	 F : nb_F x 4 les indices des points composant les mailles tetrahedriques de l'objet
 *	 T : nb_T x 3 les points du maillage d'origine dont sur lesquels on veut interpoler V et R 
 *	 V : nb_P x 3 le champ des vitesses sur le maillage volumique
 *	 R : nb_P x 1 le champ de densite (rho) sur le maillage volumique
 *	 SORTIE :
 *	 --------
 *	 indices  : nb_F x 1 indices de la maille contenant le point courant (-1 sinon)
 *	 V_interp : nb_T x 3 le champ des vitesses interpolé sur le maillage de test
 *	 R_interp : nb_T x 1 le champ de densité sur le maillage de test
*/


vtkUnstructuredGrid * buildShape_VTK(int nb_P,double * P,int nb_F, double* F) {

	vtkPoints * volume_points = vtkPoints::New();
	// Creation de la liste de points du maillage volumique 
	for (int i = 0 ; i< nb_P; i ++) {
		double temp[3];
		// La lecture des matrices récupérées de Matlab se fait colonne par colonne 
		temp[0] = P[i + nb_P * 0];
		temp[1] = P[i + nb_P * 1];
		temp[2] = P[i + nb_P * 2];
		volume_points->InsertPoint(i, temp);
	}

	//insertion des points dans un unstructured grid représentant notre maillage volumique
	vtkUnstructuredGrid * shape = vtkUnstructuredGrid::New();
	shape->SetPoints(volume_points);
	
	volume_points->Delete();

	for (int i = 0 ; i< nb_F; i ++){
		//creation de la liste de etrahedres reliant les points 	
		vtkIdList * list = vtkIdList::New();
		// La lecture des matrices récupérées de Matlab se fait colonne par colonne 
		list->InsertNextId(static_cast<int>(F[i + nb_F * 0]));
		list->InsertNextId(static_cast<int>(F[i + nb_F * 1]));
		list->InsertNextId(static_cast<int>(F[i + nb_F * 2]));
		list->InsertNextId(static_cast<int>(F[i + nb_F * 3]));
		
		shape->InsertNextCell(VTK_TETRA,list);
	        
		list->Delete();
	}

	
	return shape;


}

void getBarycentricCoordinates_VTK(vtkUnstructuredGrid * forme, int nb_T, double* T, int nb_V, double * v,int nb_R,double * R, double *indices, double* vitesse, double * densite) {
// cette fonction va localiser les tetrahedres de la forme décrite par un ensemble de tetrahedres (unstructured grid) qui contiennent les
// différents points qu'on veut tester (contenus dans la liste T)
// Une fois la cellule trouvée, on calcule les coordonnées barycentriques de chacun de ces points dans sa cellule englobante  

	// le cell locator est l'objet qui contient la méthode FindCell
	vtkCellLocator * cellLocator = vtkCellLocator::New();
	// Il contient l unstructured grid en argument
	cellLocator->SetDataSet(forme);
	cellLocator->Update();

       // boucles sur les points a tester
        for (int i =0;i<nb_T;i++) {
		double temp[3];
		// La lecture des matrices récupérées de Matlab se fait colonne par colonne 
		
                // on recupere le point courant
                temp[0] = T[i + nb_T * 0];
		temp[1] = T[i + nb_T * 1];
		temp[2] = T[i + nb_T * 2];
		
                // on trouve la cellule qui le contient
                vtkIdType ind = cellLocator->FindCell(temp);
		double res[] = {0., 0., 0., 0.}; // initialisation des coordonnées barycentriques
                // vitesse interpolee 
		double u[] = {0., 0., 0.};
		double rho=0;
		if (ind != -1)  {// il y a bien une cellule englobante
		       
			vtkCell * ma_cellule = forme->GetCell(ind);// rrecuperation de la cellule englobante
		        vtkPoints* pointsLimite = ma_cellule->GetPoints();
	 		
			double point0[3] ; pointsLimite->GetPoint(0, point0); 
			double point1[3] ; pointsLimite->GetPoint(1, point1); 
			double point2[3] ; pointsLimite->GetPoint(2, point2); 
			double point3[3] ; pointsLimite->GetPoint(3, point3); 
			
		        //calcul des coordonnées barycentriques
		        vtkTetra::BarycentricCoords(temp,point0,point1,point2,point3,res);
		        
			for (int j = 0; j < 4; j++) {
			    vtkIdType curr = ma_cellule->GetPointId(j);
			    double u_x = v[curr + nb_V * 0];
			    double u_y = v[curr + nb_V * 1];
			    double u_z = v[curr + nb_V * 2];

			    // on interpole la vitesse  u^j = sum_{i=1}^4 \lambda_i u^j_i
			    u[0] += res[j] * u_x;
			    u[1] += res[j] * u_y;
			    u[2] += res[j] * u_z;

                            rho += res[j] *R[curr];
			  
			   
			}
                } else {
			rho=1.;			
		}
		
		
		
		// affectation des coordonnées barycentriques et de l'indice de cellule dans la variable de retour
                indices[i] = static_cast<double>(ind);
		
		vitesse[i+ nb_T * 0] = u[0];
		vitesse[i+ nb_T * 1] = u[1];
		vitesse[i+ nb_T * 2] = u[2];  
        	densite[i] =rho;  
        
		
	}
	//liberation de la memoire
	cellLocator->Delete();

}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
      
        // nbr of arguments
        if (nrhs < 4){
          mexErrMsgIdAndTxt("FindTetra_VTK:Args","attend 5 arguments [ind, V_interp, R_interp] =  FindTetra_VTK(P,F,T,V,R) ");
        }

	// récupération des paramètres de la mexfunction
	size_t nb_P = mxGetM(prhs[0]);
	double * P = mxGetPr(prhs[0]);

	size_t nb_F = mxGetM(prhs[1]);
	double * F = mxGetPr(prhs[1]);

	size_t nb_T = mxGetM(prhs[2]);
	double * T = mxGetPr(prhs[2]);
	
	size_t nb_V = mxGetM(prhs[3]);
	double * V = mxGetPr(prhs[3]);

	size_t nb_R = mxGetM(prhs[4]);
	double * R = mxGetPr(prhs[4]);

	
        if (nb_V != nb_P){
          mexErrMsgIdAndTxt("FindTetra_VTK:Args","P and V doivent avoir le meme nombre de lignes");
        }
	if (nb_P != nb_R){
          mexErrMsgIdAndTxt("FindTetra_VTK:Args","P and R doivent avoir le meme nombre de lignes");
        }
	// contruction de l'objet volumique
	vtkUnstructuredGrid * forme = buildShape_VTK(nb_P, P, nb_F, F);

	plhs[0] = mxCreateDoubleMatrix(nb_T, 1,mxREAL);
	double * indices = mxGetPr(plhs[0]); 
	plhs[1] = mxCreateDoubleMatrix(nb_T, 3,mxREAL);
	double * v_interp = mxGetPr(plhs[1]); 
	plhs[2] = mxCreateDoubleMatrix(nb_T, 1,mxREAL);
	double * rho_interp = mxGetPr(plhs[2]); 
	
	
	// Appel à la fonction de recherche du tetrahedre englobant 
	// et d'extraction des coordonnees barycentriques du point  
	getBarycentricCoordinates_VTK(forme, nb_T, T, nb_V, V,nb_R,R, indices, v_interp, rho_interp);
	// il faut detruire le rigid alloue a l'interieur de la fonction  buildShape_VTK
	forme->Delete();
 
}

