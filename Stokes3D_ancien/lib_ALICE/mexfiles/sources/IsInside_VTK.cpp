#include "mex.h"
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkIntArray.h>
#include <vtkDataArray.h>
#include <vtkCellArray.h>
#include <vtkPointSet.h>
#include <iostream>
#include <ostream>
#include <fstream>
// compile with  :
// for 32 bits architectures :
// mex -v -O -I/usr/include/paraview/ IsInside_VTK.cpp -L/usr/lib/paraview/ -lvtkFiltering -lvtkCommon -lvtkGraphics;

// for 64 bits architectures :
// mex -v -O -largeArrayDims -I/usr/include/paraview/ IsInside_VTK.cpp -L/usr/lib/paraview/ -lvtkFiltering -lvtkCommon -lvtkGraphics;

using namespace std;






// fonction d'affichage d'une matrice (outil de debugging)
void affiche_matrice(size_t nb, double * m) {
	for(int i = 0 ; i < nb; ++i) {
		 for (int j = 0; j < 3 ; ++j)
		    mexPrintf("%f ",m[nb * j + i]);
	     mexPrintf("\n"); 
	}
}







// fonction qui va construire le rigide au format VTK 
// ou :
// P est la matrice des points de la surface
// F est la matrice de connectivité des points de P pour former un maillage surfacique de notre rigide 
vtkPolyData * buildShape_VTK(int nb_P, double *P , int nb_F, double *F) {
	
	
	// Contient les information relatives aux points au format VTK
	vtkPoints * surface_points = vtkPoints::New();
	
	for (int i = 0 ; i< nb_P; i ++) {
		double temp[3];
		// La lecture des matrices récupérées de Matlab se fait colonne par colonne 
		temp[0] = P[i + nb_P * 0];
		temp[1] = P[i + nb_P * 1];
		temp[2] = P[i + nb_P * 2];
		surface_points->InsertPoint(i, temp);
	}
	
	
	// Contient les information de connectivité au format VTK
	vtkCellArray * surface_faces = vtkCellArray::New();

	//insertion des cells tests dans la figure
	for (int i = 0 ; i< nb_F; i ++){
		vtkIdList * list = vtkIdList::New();
		// La lecture des matrices récupérées de Matlab se fait colonne par colonne 
		list->InsertNextId(static_cast<int>(F[i + nb_F * 0]));
		list->InsertNextId(static_cast<int>(F[i + nb_F * 1]));
		list->InsertNextId(static_cast<int>(F[i + nb_F * 2]));
		
		surface_faces->InsertNextCell(list);
	        
		list->Delete();
	}

	// Création du rigide avec les deux infos précédentes
	vtkPolyData * shape = vtkPolyData::New();
	shape->SetPoints(surface_points);
	shape->SetPolys(surface_faces);
        
	// liberation des objets intermediaires
	surface_faces->Delete();
	surface_points->Delete();

	return shape;
}







// Cette fonction va encapsuler la fonction VTK IsInside qui permet de dire si un point est 
// dans une enveloppe 3D (meme si l'enveloppe est non-convexe)
// T est la liste des points à tester
// rigid est l'enveloppe du rigide précédemment construite au format VTK
// Cette fonction retourne un tableau de doubles de la meme taille que T et qui valent : 
// 1 si le point de T correspondant est dans le rigide, 0 sinon
void testPointsInside_VTK(vtkPolyData * rigid, int nb_T, double *T, double * result) {
	
	// Construction de la liste des points à tester
    vtkPoints * testPoints = vtkPoints::New();
    
    double test_point[3];
	for (int i = 0 ; i<nb_T; i++) {
		// La lecture des matrices récupérées de Matlab se fait colonne par colonne 
		test_point[0] = T[i + 0 * nb_T];
		test_point[1] = T[i + 1 * nb_T];
		test_point[2] = T[i + 2 * nb_T];
	    testPoints->InsertNextPoint(test_point);
	}	
	
	vtkPolyData* pointsPolydata = vtkPolyData::New();
	pointsPolydata->SetPoints(testPoints);
	
	
	vtkSelectEnclosedPoints * selectEnclosedPoints = vtkSelectEnclosedPoints::New();
	
	// insertion des points a tester et de la surface dans une structure VTK 
	// qui va nous permettre d'utiliser IsInside
    #if VTK_MAJOR_VERSION <= 5
	selectEnclosedPoints->SetInput(pointsPolydata);
	#else
	selectEnclosedPoints->SetInputData(pointsPolydata);
	#endif
	#if VTK_MAJOR_VERSION <= 5
	selectEnclosedPoints->SetSurface(rigid);
	#else
	selectEnclosedPoints->SetSurfaceData(rigid);
	#endif
	// Update des informations de la structure VTK 
	// (c'est ici qu'est fait le calcul des points intérieurs)
	selectEnclosedPoints->Update();

	
	
	for (int i = 0 ; i<nb_T; i++) {
		// boucle de récupération des points intérieurs 
		result[i] = static_cast<double>(selectEnclosedPoints->IsInside(i));
	}
	
	// libération de la mémoire
	testPoints->Delete();
	pointsPolydata->Delete();
	selectEnclosedPoints->Delete();

}







void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
      
	// P : nb_P x 3 les coordonnées des points de l'enveloppe
	// F : nb_F x 3 les indices des points composant les faces de l'enveloppe 
	// T : nb_T x 3 les points à tester 

	// récupération des paramètres de la mexfunction
	size_t nb_P = mxGetM(prhs[0]);
	double * P = mxGetPr(prhs[0]);

	size_t nb_F = mxGetM(prhs[1]);
	double * F = mxGetPr(prhs[1]);

	size_t nb_T = mxGetM(prhs[2]);
	double * T = mxGetPr(prhs[2]);

	// contruction de l'enveloppe rigide au format VTK
	vtkPolyData * rigid = buildShape_VTK(nb_P, P, nb_F, F);


	plhs[0] = mxCreateDoubleMatrix(nb_T, 1,mxREAL);
	double * sortie = mxGetPr(plhs[0]); 
	// Appel à la fonction de test des points intérieurs 
	testPointsInside_VTK(rigid, nb_T, T, sortie);

	// il faut detruire le rigid alloue a l'interieur de la fonction  buildShape_VTK
	rigid->Delete();
 
}

