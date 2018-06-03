function [B_cuve_X, B_cuve_Y] = champ_magnetique_fct_Tancrede_1(v_total_2D,t_total_2D,node_aimant,edge_aimant,node_ensemble,edge_ensemble,aimant_centre)
%============================================================================================================
% Solveur de Laplace - Elements finis P1
%
% Calcul du champ magnetique sur un maillage 2D donnee
%============================================================================================================
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

%------------------------------------------------------------------------------------------------------------
%                   definiton de fnum pour situer pour chaque triangle du maillage sur quelle
%                   face du maillage il est
%------------------------------------------------------------------------------------------------------------
% face 1 = aimant
% face 2 = cuve
% face 3 = le reste avec le maillage de resolution total

fnum = 3 * ones(size(t_total_2D,1),1); %fnum a 3 de base apres on cherche la cuve et l'aimant

if aimant_centre == true %si l'aimant est en position centree
    int_aimant = find(v_total_2D(:,1) >= 4.75 & v_total_2D(:,1) <= 5.25 & v_total_2D(:,2) >=  3.5 & v_total_2D(:,2) <= 4.5);
end
if aimant_centre == false %si l'aimant est sur le coté
    int_aimant = find(v_total_2D(:,1) >= 5.1 & v_total_2D(:,1) <= 5.6 & v_total_2D(:,2) >=  3.5 & v_total_2D(:,2) <= 4.5);
end
int_cuve = find(v_total_2D(:,1) >= 4.5 & v_total_2D(:,1) <= 5.5 & v_total_2D(:,2) >=  4.5 & v_total_2D(:,2) <= 5.5);

for k = 1:nt
        
        if ismember(t_total_2D(k,1),int_aimant) && ismember(t_total_2D(k,2),int_aimant) && ismember(t_total_2D(k,3),int_aimant);
           fnum(k) = 1; %aimant
        end
        if ismember(t_total_2D(k,1),int_cuve) && ismember(t_total_2D(k,2),int_cuve) && ismember(t_total_2D(k,3),int_cuve);
           fnum(k) = 2; %cuve
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


[M, nn, ibint, ic2] = matrixP1final_1(v_total_2D,t_total_2D,fnum,node_aimant,edge_aimant,node_ensemble,edge_ensemble);

[A1, M1]= matrixP1vect(v_total_2D,t_total_2D);

%==========================================================================
%                        resolution du systeme lineaire
%==========================================================================
A1 = A1(ic2,ic2);

sol = A1\M;


%--------------------------------------------------------------------------
%       reconstituion de la solution sur les noeuds interieurs
%--------------------------------------------------------------------------

u = zeros(nn,1);
u(ic2) = sol;
u(ibint) = 0;  

%---------------------------------------------------------------------------
%                           calcul du gradient
%---------------------------------------------------------------------------

B = gradient(u,v_total_2D,t_total_2D);

%---------------------------------------------------------------------------
%                  recuperation de la solution sur la cuve
%---------------------------------------------------------------------------

j = 1;
for k = 1:nt
        
        if fnum(k) == 2
            t_cuve(j,:) = t_total_2D(k,1:3);
            j = j + 1;
        end
end
numeros_noeuds_cuve = unique(sort(t_cuve(:)));
v_cuve = v_total_2D(numeros_noeuds_cuve,:);

u_cuve = u(numeros_noeuds_cuve,:);
B_cuve = B(numeros_noeuds_cuve,:);

%affichage norme B pour verification
% figure();
% patch('faces',t_cuve(:,1:3),'vertices',v_total_2D(:,1:2),'facecolor','white','edgecolor','blue');
% axis equal;  
% colorbar;
% title('maillage cuve réarrangé');
%  

B_cuve_X = B_cuve(:,1);
B_cuve_Y = B_cuve(:,2);

B_total_X = B(:,1);
B_total_Y = B(:,2);

X_cuve_2D = v_cuve(:,1);
Y_cuve_2D = v_cuve(:,2);

X_total_2D = v_total_2D(:,1);
Y_total_2D = v_total_2D(:,2);

%==========================================================================
%                        affichage de la solution
%==========================================================================
 
% figure();
% quiver(X_total_2D,Y_total_2D,B_total_X,B_total_Y);
% xlim([0 10]);
% ylim([0 10]);
% title('lignes de champ magnetique sur le maillage de resolution');
 
 
M2 = zeros(size(v_total_2D,1),1);
M2(ic2) = M;
 
%  %  affichage norme B pour verification
%  figure();
%  patch('faces',t_total_2D(:,1:3),'vertices',v_total_2D(:,1:2),'FaceVertexCData',sqrt(B_cuve_X.^2+B_cuve_Y.^2),'facecolor','interp','edgecolor','black');
%  axis equal;  
%  colorbar;
%  title('norme de B');

% figure();
% hold on;
% quiver(X_cuve_2D,Y_cuve_2D,B_cuve_X,B_cuve_Y);
% xlim([4.5 5.5]);
% ylim([4.5 5.5]);
% % x = linspace(0,10,100);
% % y1 = afficher_aimant_sup(x);
% % y2 = afficher_aimant_inf(x);
% % plot(y1,x);
% % plot(y2,x);
% title('lignes de champ magnetique dans la cuve');
% xlabel('X'); ylabel('Y');
% hold off

end
 
