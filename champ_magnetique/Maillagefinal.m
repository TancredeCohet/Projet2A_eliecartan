%-------------------------------------------------------------------------
% Solveur de Laplace - Elements finis P1
%
%    - Delta u = f dans Omega=[0,L]x[0,H]
%
% Conditions limites :
%   u=g sur le bord exterieur du domaine {x=0},{x=L},{y=0} et {y=H}. 
%
% Maillage par mesh2d -- GNU GPL - Copyright (C) 2007 Darren Engwirda
%
%%-------------------------------------------------------------------------
%addpath('code prof'i)
addpath('code prof');
addpath('code prof/Mesh2d_v23');
addpath('code prof/Mesh2d_v23\private');
addpath('code prof/Mesh2d_v23\lib');


%==========================================================================
%                     elaboration du maillage
%==========================================================================

%taille de l'aimant carre de taille a
a = 0.75;

%taille de la cuve
L = 10; %longueur
H = 10; %hauteur
node1 = [-a,-a; a,-a; a,a; -a, a];          %liste des noeuds int�rieurs
node2 = [-L,-H; L,-H; L,H; -L, H];          %liste des noeuds ext�rieurs
edge1 = [(1:size(node1,1))',[(2:size(node1,1))'; 1]];       %liste des ar�tes int�rieures
edge2 = [1,2; 2,3; 3,4; 4,1];                               %liste des ar�tes ext�rieures

edge = [edge1; edge2+size(node1,1)];        %liste de toutes les ar�tes mises dans l'ordre
node = [node1; node2];                      %liste de tous les noeuds mises dans l'ordre

pas_h=0.3;                                  %pas maximal du maillage
hdata = [];
hdata.hmax  = pas_h;
face{1} = 1:size(edge1,1);                  %face de maillage interieure
face{2} = 1:size(edge,1);                   %face de maillage exterieure
 
[v,t,fnum] = meshfaces(node,edge,face,hdata); %construction du maillage

%==========================================================================
%                              Calcul des matrices
%
%                       A est la matrice de rigidit�
%                       M la matrice du terme source
%                             nv2 nombre de noeuds
%                             ibint noeuds interieurs
%                             nc2 noeuds du bord
%==========================================================================
%    Aretes du maillage P1
t = [t ones(size(t,1),1) ];
t = area_triangle(v,t);

[M, nn, ibint, ic2] = matrixP1final(v,t,fnum,node1,edge1,node2,edge2);

[A1, M1]= matrixP1vect(v,t);

%==========================================================================
%                        resolution du systeme lineaire
%==========================================================================
A1 = A1(ic2,ic2);

sol = A1\M;


%--------------------------------------------------------------------------
%reconstituion de la solution avec les noeuds int�rieur car la solution
%actuelle ne porte que sur les noeuds int�rieurs du maillage
%--------------------------------------------------------------------------
close all; 
u = zeros(nn,1);
u(ic2) = sol;
u(ibint) = 0;  

%u = zeros(nn,1);


B = gradient(u,v,t);
%on observe que pour u = 1 sur tout le maillage le gradient calcul� sur les
%noeuds est de l'ordre de 1e-14

B_X = B(:,1);
B_Y = B(:,2);

%on recupere les coordonn�es de chaque point du maillage
X = v(:,1);
Y = v(:,2);



%==========================================================================
%                                 affichage de la solution
%==========================================================================

 figure;
 
 M2 = zeros(nn,1);
 M2(ic2) = M;
 
 patch('faces',t(:,1:3),'vertices',v(:,1:2),'FaceVertexCData',sqrt(B_X.^2+B_Y.^2),'facecolor','interp','edgecolor','black');
 axis equal; 
 colorbar;
 figure;
 
 %patch('faces',t(:,1:3),'vertices',v(:,1:2),'facecolor','red');
 quiver(X,Y,B_X,B_Y);
 xlim([-2 2]);
 ylim([-2 2]);
 

