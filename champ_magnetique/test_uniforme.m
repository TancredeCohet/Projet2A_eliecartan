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
close all;
% 
% Path pour mesh2d (maillage)
addpath('Mesh2d_v23');
addpath 'lib'

fprintf('------------------------------------------\n');
fprintf('        Solveur Laplacien EF-P1           \n');
fprintf('------------------------------------------\n');
%----------------------------
% Définition de la géometrie
%----------------------------
H = 2;         % hauteur du domaine
L = 3;            % longueur du domaine
pas_h = 0.1;   % pas du maillage "interieur"

%--------------------------------------------------------------------------
% Maillage
%----------------------------------------------------------------------------
% Bord exterieur du domaine (rectangle) 
node = [0 0; L 0; L H; 0 H];
edge = [1 2; 2 3; 3 4; 4 1];

fprintf('-----------------------------------------------------------------\n');
fprintf('Maillage par ''mesh2d'' v2.3 (Copyright (C) 2007 Darren Engwirda)\n');
fprintf('-----------------------------------------------------------------\n');
tstart=tic; %lance chrono temps de calcul

%implementation longueur des arete
hdata = [];
hdata.hmax  = pas_h; 

%definition des options du maillage
optionsmesh = [];
optionsmesh.output = false;


outmesh = meshP1(node,edge,hdata,optionsmesh);

v1=outmesh.v1; 
t1=outmesh.t1; 
e1=outmesh.e1;
ibd1=outmesh.ib1; %noeud du bord
iin1=outmesh.ic1; %noeud intérieur

% nv1, nt1, ne1 : nb de noeuds, de triangles et d'aretes du maillage P1
nv1 = size(v1,1); %nombre de noeuds
nt1 = size(t1,1); %nombre de trizngle
ne1 = size(e1,1); %nombre d'arette

fprintf('Maillage P1 : nombre de noeuds    (nv1) = %d\n',nv1);
fprintf('              nombre de triangles (nt1) = %d\n',nt1);
fprintf('Temps de construction du maillage: %f s\n',toc(tstart));

%affichage du maillage
patch('faces',t1(:,1:3),'vertices',v1(:,1:2),'facecolor','white','edgecolor','blue');
axis equal; axis on
drawnow 


%-----------------------------------------------
% Calcul des matrices de rigiditÃ© et de masse P1
%-----------------------------------------------
tstart=tic;
fprintf('------------------------------------\n');
fprintf('Construction des matrices P1 \n');
fprintf('------------------------------------\n');

[A, M] = matrixP1vect(v1,t1);
fprintf('Temps d''assemblage des matrices P1: %f s\n',toc(tstart));
fprintf('Temps d''assemblage des matrices P1: %f s\n',toc(tstart));
%tstart=tic
%-------------------------------------------------------------------------
% Traitement des conditions limites
%--------------------------------------------------------------------------
tstart = tic;
% conditions limites u=g pour le bord

x = v1(ibd1,1);  %selectionne absisses des noeuds du maillage sur le bord
y = v1(ibd1,2);  %selectionne ordonnés des noeuds du maillage sur le bord
ubd = 0*gcl(x,y);  %créer actuellement matrice  de 1 de taille nb_noeud_bord*1

% Terme source
x = v1(:,1);     %selectionne absisses des noeuds du maillage
y = v1(:,2);     %selectionne ordonnés des noeuds du maillage
f = fscmb(x,y);  %creation terme source actuellement matrice de 1

% Second membre complet
b = M*f - A(:,ibd1)*ubd;
b(ibd1) = [];

% Mise Ã  jour de la matrice A : on ne retient que les noeuds interieurs 
A = A(iin1,iin1);
fprintf('Taille de la matrice A : %d x %d (nnz = %d)\n',size(A), nnz(A));
fprintf('Temps de traitement des CL: %f s\n',toc(tstart));

%--------------------------------------------------------------------------
% Résolution du systeme lineaire
%--------------------------------------------------------------------------
tstart = tic;
sol = A\b;

% reconstruction de la solution dans tout le domaine
u = zeros(nv1,1);
u(iin1) = sol;
u(ibd1) = ubd;   

fprintf('Temps de résolution du système linéaire: %f s\n',toc(tstart));

%--------------------------------------------------------------------------
% Affichage du maillage et de la solution
%--------------------------------------------------------------------------

%calcul gradient
%u = ones(nv1,1);
B = gradient(u,v1,t1);
B(ibint) = 0;
B_X = B(:,1);
B_Y = B(:,2);
X = v1(:,1);
Y = v1(:,2);

figure;
patch('faces',t1(:,1:3),'vertices',v1(:,1:2),'FaceVertexCData',u,...
      'facecolor','interp','edgecolor','black');
axis equal;
colorbar;

figure;
quiver(X,Y,B_X,B_Y);
axis equal;
colorbar;



