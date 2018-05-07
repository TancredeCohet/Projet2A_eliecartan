function [B_X,B_Y] = champ_magnetique_fct_Tancrede(v_total_2D,t_total_2D,node_aimant,edge_aimant,node_total,edge_total)
%-------------------------------------------------------------------------
% Solveur de Laplace - Elements finis P1
%
%    Calcul du champ magnetique sur un maillage 2D donnee
%%-------------------------------------------------------------------------
% L largueur de la cuve
% P prfondeur de la cuve

addpath('../champ_magnetique/code_prof');
addpath('../champ_magnetique/code_prof/Mesh2d_v23');
addpath('../champ_magnetique/code_prof/Mesh2d_v23/private');
addpath('../champ_magnetique/code_prof/Mesh2d_v23/lib');
addpath('../Stokes3D')

nt = size(t_total_2D,1); %nombre de triangle



% ajoute aire triangle
t_total_2D = [t_total_2D ones(size(t_total_2D,1),1) ];
t_total_2D = area_triangle(v_total_2D,t_total_2D); 

%taille de l'aimant carre de taille a
% a = taille_aimant;


% definiton de fnum pour situer pour chaque triangle du maillage sur quelle
% face du maillage il est4% face 1 = aimant
% face 2 = cuve
% face 3 = le reste avec le maillage de résolution total

fnum = 3*ones(size(t_total_2D,1),1); %fnum a 3 de base apres on cherche la cuve et l'aimant

int_aimant = find(v_total_2D(:,1) >= 5 & v_total_2D(:,1) <= 7 & v_total_2D(:,2) >=  2 & v_total_2D(:,2) <= 5);
int_cuve = find(v_total_2D(:,1) >= 3 & v_total_2D(:,1) <= 7 & v_total_2D(:,2) >=  4 & v_total_2D(:,2) <= 8);

for k = 1:nt
        
        if ismember(t_total_2D(k,1),int_aimant) & ismember(t_total_2D(k,2),int_aimant) & ismember(t_total_2D(k,3),int_aimant);
           fnum(k) = 1;
        end 
        if ismember(t_total_2D(k,1),int_cuve) & ismember(t_total_2D(k,2),int_cuve) & ismember(t_total_2D(k,3),int_cuve);
           fnum(k) = 2;
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


[M, nn, ibint, ic2] = matrixP1final(v_total_2D,t_total_2D,fnum,node_aimant,edge_aimant,node_total,edge_total);

[A1, M1]= matrixP1vect(v_total_2D,t_total_2D);

%==========================================================================
%                        resolution du systeme lineaire
%==========================================================================
A1 = A1(ic2,ic2);

sol = A1\M;


%--------------------------------------------------------------------------
%reconstituion de la solution avec les noeuds interieurs car la solution
%actuelle ne porte que sur les noeuds interieurs du maillage
%--------------------------------------------------------------------------
close all; 
u = zeros(nn,1);
u(ic2) = sol;
u(ibint) = 0;  


B = gradient(u,v_total_2D,t_total_2D);
%on observe que pour u = 1 sur tout le maillage le gradient calcule sur les
%noeuds est de l'ordre de 1e-14



% recuperation de la solution sur la cuve
% t_cuve = []; %triangle cuve
j = 1;
for k = 1:nt
        
        if fnum(k)==2
            t_cuve(j,:) = t_total_2D(k,1:3);
            j = j + 1;
        end
end
numeros_noeuds_cuve = unique(sort(t_cuve(:)));
v_cuve = v_total_2D(numeros_noeuds_cuve,:);

u_cuve = u(numeros_noeuds_cuve,:);
B_cuve = B(numeros_noeuds_cuve,:);

B_X = B_cuve(:,1);
B_Y = B_cuve(:,2);

%on recupere les coordonnees de chaque point du maillage
% X_2D = v_total_2D(:,1);
% Y_2D = v_total_2D(:,2);
X_2D = v_cuve(:,1);
Y_2D = v_cuve(:,2);
%==========================================================================
%                                 affichage de la solution
%==========================================================================

 figure;
 
 M2 = zeros(nn,1);
 M2(ic2) = M;
 
 %affichage norme B
 patch('faces',t_total_2D(:,1:3),'vertices',v_total_2D(:,1:2),'FaceVertexCData',sqrt(B_X.^2+B_Y.^2),'facecolor','interp','edgecolor','black');
 axis equal; 
 colorbar;
 figure;
 
 %patch('faces',t(:,1:3),'vertices',v(:,1:2),'facecolor','red');
 quiver(X_2D,Y_2D,B_X,B_Y);
% xlim([-2 2]);
% ylim([-2 2]);
end
 

