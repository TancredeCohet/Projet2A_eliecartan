%-------------------------------------------------------------------------
% Solveur de Laplace - Elements finis P1
%
%    Calcul du champ magnétique sur un maillage 2D donnée
%%-------------------------------------------------------------------------
function [B_X,B_Y] = champ_magnetique_fct(vz,tz,L,P,taille_aimant);

%   v liste des cordonnées des neuds
%   t matrice 3*N des triangle

% L largueur de la cuve
% P prfondeur de la cuve

addpath('../champ magnétique/code prof');
addpath('../champ magnétique/code prof/Mesh2d_v23');
addpath('../champ magnétique/code prof/Mesh2d_v23/private');
addpath('../champ magnétique/code prof/Mesh2d_v23/lib');

tz = [tz ones(size(tz,1),1) ];
tz = area_triangle(vz,tz); %on ajoute l'aire des triangles en 5ieme colonne

%taille de l'aimant carre de taille a
a = taille_aimant;

node1 = [-a,-a; a,-a; a,a; -a, a];          %liste des noeuds interieurs
node2 = [-L,-P; L,-P; L,P; -P, P];          %liste des noeuds exterieurs
edge1 = [(1:size(node1,1))',[(2:size(node1,1))'; 1]];       %liste des aretes interieures
edge2 = [1,2; 2,3; 3,4; 4,1];                               %liste des aretes exterieures

edge = [edge1; edge2+size(node1,1)];        %liste de toutes les aretes mises dans l'ordre
node = [node1; node2];                      %liste de tous les noeuds mises dans l'ordre

nt=size(tz,1);
index_suitable_triangle=ones(nt,1);
%[v,t,fnum] = meshfaces(node,edge,face,hdata); %construction du maillage

fnum = 2*ones(size(vz,1),1);
int = find(vz(:,1) >= taille_aimant & vz(:,1) <= 2*taille_aimant & vz(:,2) >= taille_aimant & vz(:,2) <= 2*taille_aimant);
for k=1:nt;
        
        if ismember(tz(k,1),int) & ismember(tz(k,2),int) & ismember(tz(k,3),int);
           fnum(k) = 1;
        end
end


%==========================================================================
%                              Calcul des matrices
%==========================================================================

%                       A est la matrice de rigidite
%                       M la matrice du terme source
%                             nv2 nombre de noeuds
%                             ibint noeuds interieurs
%                             nc2 noeuds du bord


[M, nn, ibint, ic2] = matrixP1final(vz,tz,fnum,node1,edge1,node2,edge2);

[A1, M1]= matrixP1vect(vz,tz);

%==========================================================================
%                        resolution du systeme lineaire
%==========================================================================
A1 = A1(ic2,ic2);

sol = A1\M;


%--------------------------------------------------------------------------
%reconstituion de la solution avec les noeuds intérieurs car la solution
%actuelle ne porte que sur les noeuds intérieurs du maillage
%--------------------------------------------------------------------------
close all; 
u = zeros(nn,1);
u(ic2) = sol;
u(ibint) = 0;  

%u = zeros(nn,1);


B = gradient(u,vz,tz);
%on observe que pour u = 1 sur tout le maillage le gradient calcule sur les
%noeuds est de l'ordre de 1e-14

B_X = B(:,1);
B_Y = B(:,2);

%on recupere les coordonn�es de chaque point du maillage
X = vz(:,1);
Y = vz(:,2);



%==========================================================================
%                                 affichage de la solution
%==========================================================================

 figure;
 
 M2 = zeros(nn,1);
 M2(ic2) = M;
 
 %affichage norme B
 patch('faces',tz(:,1:3),'vertices',vz(:,1:2),'FaceVertexCData',sqrt(B_X.^2+B_Y.^2),'facecolor','interp','edgecolor','black');
 axis equal; 
 colorbar;
 figure;
 
 %patch('faces',t(:,1:3),'vertices',v(:,1:2),'facecolor','red');
 quiver(X,Y,B_X,B_Y);
 xlim([-2 2]);
 ylim([-2 2]);
end
 

