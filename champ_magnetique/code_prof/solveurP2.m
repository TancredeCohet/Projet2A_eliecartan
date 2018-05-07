%-------------------------------------------------------------------------
% Solveur de Laplace - Elements finis P2
%
%    - Delta u = f dans Omega=[0,L]x[0,H]
%
% Conditions limites :
%   u=g sur le bord exterieur du domaine {x=0},{x=L},{y=0} et {y=H}. 
%
% Maillage par mesh2d -- GNU GPL - Copyright (C) 2007 Darren Engwirda
%
%%-------------------------------------------------------------------------
close all;
% 
% Path pour mesh2d (maillage)
addpath('Mesh2d_v23');
addpath 'lib'

fprintf('------------------------------------------\n');
fprintf('        Solveur Laplacien EF-P2           \n');
fprintf('------------------------------------------\n');
%----------------------------
% Définition de la géométrie
%----------------------------
H=1.05;           % hauteur du domaine
L=1;                 % longueur du domaine
pas_h = 4.02;   % pas du maillage "interieur" %0.04 intialement
%------------
% Maillage
%------------
% Bord exterieur du domaine (rectangle) 
node = [0 0; L 0; L H; 0 H];
edge = [1 2; 2 3; 3 4; 4 1];

fprintf('-----------------------------------------------------------------\n');
fprintf('Maillage par ''mesh2d'' v2.3 (Copyright (C) 2007 Darren Engwirda)\n');
fprintf('-----------------------------------------------------------------\n');
tstart=tic;

hdata = [];
hdata.hmax  = pas_h; 

optionsmesh=[];
optionsmesh.output=false;

outmesh=meshP1P2(node,edge,hdata,optionsmesh);
    
v1=outmesh.v1; t1=outmesh.t1; e1=outmesh.e1;
v2=outmesh.v2; t2=outmesh.t2; e2=outmesh.e2;

ibd2=outmesh.ib2; iin2=outmesh.ic2;

% nv1, nt1, ne1 : nb de noeuds, de triangles et d'aretes du maillage P1
% nv2, nt2, ne2 : nb de noeuds, de triangles et d'aretes du maillage P2
nv1=size(v1,1); nt1=size(t1,1); ne1=size(e1,1);
nv2=size(v2,1); nt2=size(t2,1); ne2=size(e2,1);

fprintf('Maillage P1 : nombre de noeuds    (nv1) = %d\n',nv1);
fprintf('              nombre de triangles (nt1) = %d\n',nt1);
fprintf('Maillage P2 : nombre de noeuds    (nv2) = %d\n',nv2);
fprintf('              nombre de triangles (nt2) = %d\n',nt2);

fprintf('Temps de construction du maillage: %f s\n',toc(tstart));

% Affichage du maillage
clf()
patch('faces',t1(:,1:3),'vertices',v1(:,1:2),'facecolor','white','edgecolor','blue');
axis equal; axis on
drawnow 

%-----------------------------------------------
% Calcul des matrices de rigidité et de masse P2
%-----------------------------------------------
tstart=tic;
fprintf('------------------------------------\n');
fprintf('Construction des matrices P2 \n');
fprintf('------------------------------------\n');
%[A,M]=matrixP2(v1,t1,ne1);
%fprintf('Temps d''assemblage des matrices P2: %f s\n',toc(tstart));
%tstart=tic;
[A,M]=matrixP2vect(v1,t1);
fprintf('Temps d''assemblage des matrices P2: %f s\n',toc(tstart));

%-----------------------------------
% Traitement des conditions limites
%-----------------------------------
tstart=tic;
% conditions limites u=g pour le bord
x=v2(ibd2,1);  
y=v2(ibd2,2);
ubd=gcl(x,y);

% Terme source
% x=v2(:,1);  
% y=v2(:,2);
% f=fscmb(x,y); 
% b = M*f;
b = intquadR52(@fscmb,t1,v1,t2,v2);

% Second membre complet
b = b - A(:,ibd2)*ubd;
b(ibd2)=[];

% Mise à jour de la matrice A : on ne retient que les noeuds interieurs 
A=A(iin2,iin2);
fprintf('Taille de la matrice A : %d x %d (nnz = %d)\n',size(A), nnz(A));
fprintf('Temps de traitement des CL: %f s\n',toc(tstart));

%--------------------------------
% Résolution du systeme lineaire
%--------------------------------
tstart=tic;
sol=A\b;

% reconstruction de la solution dans tout le domaine
u=zeros(nv2,1);
u(iin2)=sol;
u(ibd2)=ubd;   

fprintf('Temps de résolution du système linéaire: %f s\n',toc(tstart));

%-----------
% Affichage
%-----------
clf;
patch('faces',t2(:,1:3),'vertices',v2(:,1:2),'FaceVertexCData',u,...
      'facecolor','interp','edgecolor','black');
axis equal; 
%axis off
colorbar;


