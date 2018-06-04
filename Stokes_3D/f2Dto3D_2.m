function [fx, fy, fz, fxgd, fygd,fzgd,fxy] = f2Dto3D_2(v,t,nv1,nv2,nbquad,X,Y,Z,Nx,Ny,Nz,aimant_centre)
%==========================================================================
%     a partir d'un maillage donnee par les parametres f2Dto3D82 renvoit 
%   le terme source de l'equation de navier-stokes
%
%   fx, fy, fz sont les forces dans des vecteurs colonnes 1*size(v,1)
%   fxgd, fygd, fzgd sont les forces dans des matrices Nx*Ny*Ny pour
%   affichage
%==========================================================================
%--------------------------------------------------------------------------
%   Construction du maillage triangulaire 2D pour z=0 a partir du maillage 3D
%   par des cubes
%--------------------------------------------------------------------------


close all;
addpath('../champ_magnetique');

%--------- extraction d'une couche 2d et triangulation --------------------
lz = 0;                             % niveau z=0
indz = find(v(:,3) == lz);
xz = v(indz,1); yz = v(indz,2); vz = [xz, yz];
tz = delaunay(xz,yz);               % triangulation de Delaunay des points du plan 

figure();
triplot(tz,xz,yz);   % affichage
title('couche 2D extraite du maillage 3D de la cuve')
%_________________________________________________________________________%
%         on voit que les point dans v ne sont pas dans l'ordre           %
%_________________________________________________________________________%

%--------------------------------------------------------------------------
% Calcul du champ magnetique et de la force avec le maillage 2D defini par
% (tz,vz)
%--------------------------------------------------------------------------

%coordonnees des differents maillage
node_ensemble = [0,0 ; 10 0; 10 10; 0 10]; %  maillage total
node_select = [4.5 4.5;5.5 4.5; 5.5 5.5; 4.5 5.5];        %  mailage cuve cote 4

%---------- configuration aimant ------------------------------------------
if aimant_centre == true
    node_aimant = [4.75 3.5; 5.25 3.5; 5.25 4.5; 4.75 4.5];        %  maillage aimant
end
if aimant_centre == false
    node_aimant = [5.1 3.5;5.6 3.5;5.6 4.5;5.1 4.5];        %  maillage aimant
end

edge_aimant = [(1:size(node_aimant,1))',[(2:size(node_aimant,1))'; 1]];       %liste des aretes interieures
edge_ensemble = [1,2; 2,3; 3,4; 4,1];                               %liste des aretes exterieures

edge_total = [edge_aimant; edge_ensemble+size(node_aimant,1)];        %liste de toutes les aretes mises dans l'ordre
node_total = [node_aimant; node_ensemble];                      %liste de tous les noeuds mises dans l'ordre


%generation maillage total de resolution du probleme 2D
L_total = 10;
H_total = 10;
x_2D = linspace(0,L_total,2*10*(Nx-1) + 1)';
y_2D = linspace(0,H_total,2*10*(Ny-1) + 1 )';

[X_total_2D,Y_total_2D] = meshgrid(x_2D,y_2D);
v_total_2D_x = X_total_2D';
v_total_2D_x = v_total_2D_x(:);
v_total_2D_y = Y_total_2D';
v_total_2D_y = v_total_2D_y(:);
v_total_2D = [v_total_2D_x(:) v_total_2D_y(:)];

t_total_2D = delaunay(v_total_2D(:,1),v_total_2D(:,2)); %generation maillage triangulaire total

 figure();
 triplot(t_total_2D ,X_total_2D,Y_total_2D);
 xlim([4.5 5.5]);
 ylim([4.5 5.5]);
 title('maillage de resolution 2D total généré')


%u = zeros(size(X_total_2D(:),1),1);
% 
%  figure();
%  patch('faces',t_total_2D(:,1:3),'vertices',v_total_2D(:,1:2),'facecolor','white');
%  title('avec path')
%  xlim([4.5 5.5]);
%  ylim([4.5 5.5]);

%--------------------------------------------------------------------------
%               calcul du champ magnétique sur la couche z=0
%--------------------------------------------------------------------------

[B_X,B_Y] = champ_magnetique_fct_Tancrede_1(v_total_2D,t_total_2D,node_aimant,edge_aimant,node_ensemble,edge_ensemble,aimant_centre);


f0 = B_Y; %calcul du second menbre sur le plan z=0

%==========================================================================
%                 placement du second menbre dans le maillage 3D
%==========================================================================

fx = zeros(size(v,1),1);
fy = zeros(size(v,1),1);
fz = zeros(size(v,1),1);

% Pour Stokes3D la force magnetique doit etre sous la forme d'une matrice
% colonne nb_noeuds * 1

%----------------- marche ?  ------------------------------------------
%   l'idee est de recuperer les coordonnees au-dessus de chaque noeuds pour
%   assigner les meme valeurs que la couche z=0
% for k=1:size(f0,1)
%     iz = find(abs((v(:,1) - vz(k,1)))<0.01 & abs((v(:,2) - vz(k,2)))<0.01);
%     if size(iz,1) ~= 21
%         lk = 1;
%     end
%     fz(iz,:) = f0(k);                                     
% end
%--------------------- test table de connectivite ------------------------
% test_connectivite = true;
% for ti = 1:size(t,1)
%     a = (fz(t(ti,1))==fz(t(ti,13))) && (fz(t(ti,1))==fz(t(ti,5))) && (fz(t(ti,2))==fz(t(ti,14))) && (fz(t(ti,2))==fz(t(ti,6))) && ...
%         (fz(t(ti,3))==fz(t(ti,15))) && (fz(t(ti,3))==fz(t(ti,7))) && (fz(t(ti,4))==fz(t(ti,16))) && (fz(t(ti,4))==fz(t(ti,8))) && ...
%         (fz(t(ti,9))==fz(t(ti,22))) && (fz(t(ti,9))==fz(t(ti,17))) && (fz(t(ti,10))==fz(t(ti,23))) && (fz(t(ti,10))==fz(t(ti,18))) && ...
%         (fz(t(ti,11))==fz(t(ti,24))) && (fz(t(ti,11))==fz(t(ti,19))) && (fz(t(ti,12))==fz(t(ti,25))) && (fz(t(ti,12))==fz(t(ti,20))) && ...
%         (fz(t(ti,21))==fz(t(ti,27))) && (fz(t(ti,21))==fz(t(ti,26)));
%     if a == false
%         test_connectivite = false
%     end
% end
%   les test montre qu'il n'y pas de probleme avec les tables de
%   connictivite
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

%----------- placement des valuers de f0 dans une matrice 3D --------------
%   placement des valuers sur une matrice Nx*Ny*Nz
%   d'abord placement de f0 sur une matrice Nx*Ny comme les valeurs sont
%   les meme sur chaque plan du maillage selon z
fxy = zeros(2*Ny-1,2*Nx-1);
for yi = 1:2*Ny-1       % selon l'axe y
    for xj = 1:2*Nx-1         % selon l'axe x
        m = (2*Nx-1)*(yi-1)+xj ;
        fxy(yi,xj) = f0(m);
    end
end
fxgd=zeros(2*Ny-1,2*Nx-1,2*Nz-1);
fygd=zeros(2*Ny-1,2*Nx-1,2*Nz-1);
fzgd=zeros(2*Ny-1,2*Nx-1,2*Nz-1);

% comme la valeur selon z est l seule non-nulle on ne touche pas aux autres
% composante ce qui serait des calculs inutiles
for p=1:(2*Nz-1)
    fzgd(:,:,p)=fxy;
end
%--------------------------------------------------------------------------

%---constitution de fz a partir de fzgd en respectant la structure de v --%
hx=1/(2*(Nx-1)); hy=1/(2*(Ny-1)); hz=1/(2*(Nz-1));
%fzgd=zeros(2*Ny-1,2*Nx-1,2*Nz-1);
for m = 1:size(fz,1)
    
    ind = round(v(m,:)./[hx, hy, hz])+1;
    fz(m,:) = fzgd(ind(2),ind(1),ind(3));
end
%-------------------------------------------------------------------------%


%--------------- placement des valeurs dans fz à partir de fzgd ----------%
% for zi = 1:(2*Nz-1)
%     for xi = 1:(2*Nx-1)
%         fz((xi -1) * (2*Nx-1) + 1 + (zi -1) * ((2*Nz - 1) ^ 2):xi * (2*Nx - 1) + (zi-1) * ((2*Nz - 1) ^ 2)) = fzgd(xi,:,zi)';
%     end
% end
%--------------------------------------------------------------------------

%----------- mise ne forme matrice 3D pour affichage avec quiver ----------
%  ainsi onverifie le bon placement des donees
%  cela necessite le rearrangement des donnees dans une matrice 3D
 hx = 1/(2*(Nx-1)); hy = 1/(2*(Ny-1)); hz = 1/(2*(Nz-1));
 for m = 1:size(fz,1)
     ind = round(v(m,:)./[hx, hy, hz])+1;
     fzgd(ind(2),ind(1),ind(3)) = fz(m);
 end
%--------------------------------------------------------------------------



















% Affichage de fx, juste pour vesrifier
% Construction du maillage du cube [0,L]x[|O,P]x[0,H]
% x_2D = linspace(0,L,2*Nx-1);
% y_2D = linspace(0,P,2*Ny-1);
% zl = linspace(0,H,2*Nz-1);
% [X,Y,Z]=meshgrid(x_2D,y_2D,zl);
% 
% hx=L/(2*(Nx-1)); hy=P/(2*(Ny-1)); hz=H/(2*(Nz-1));
% fxgd=zeros(2*Ny-1,2*Nx-1,2*Nz-1);
% 
% for m=1:size(v,1)
%     ind=round(v(m,:)./[hx, hy, hz])+1;
%     fxgd(ind(2),ind(1),ind(3))=fx(m);
% end
% figure()
% xslice = [0,L/2,L]; yslice = P/2; zslice = [0,H/2,H];
% slice(X,Y,Z,fxgd,xslice,yslice,zslice,'linear')
% xlabel('X'); ylabel('Y'); zlabel('Z');
end
