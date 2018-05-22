function [fx, fy, fz, fxgd, fygd,fzgd] = f2Dto3D_2(v,t,nv1,nv2,nbquad,X,Y,Z)
%==========================================================================
%     a partir d'un maillage donnee par les parametres f2Dto3D82 renvoit 
%   le terme source de l'equation de navier-stokes
%
%   fx, fy, fz sont les forces dans des vecteurs colonnes 1*size(v,1)
%   fxgd, fygd, fzgd sont les forces dans des matrices Nx*Ny*Ny pour
%   affichage
%==========================================================================



close all;
addpath '/lib_ALICE/mexfiles/bin'
addpath '/lib'
addpath('../champ_magnetique');

%   parametre de la cuve
global L P H
L = 1;  % largeur
P = 1;  % profondeur
H = 1;  % hauteur
nu = 1; % viscosite

% Nombres de points dans les directions x, y et z pour la cuve
Nx = 11;
Ny = 11;
Nz = 11;


%generation du maillage 3D
[v,t,nv1,nv2,nbquad,X,Y,Z]= maillage_3D(Nx,Ny,Nz,L,P,H);

fprintf('*/ Maillage du cube %dx%dx%d (%f s)\n',Nx,Ny,Nz); 
fprintf('     nombre d''elements (cubes) = %d\n', size(t,1));
fprintf('     nombre de d.o.f pour la vitesse = %d\n', nv2);
fprintf('     nombre de d.o.f pour la pression = %d\n', nv1);

% Construction du maillage triangulaire 2D pour z=0 a partir du maillage 3D
% par des cubes.
lz = 0;                % niveau z=0
indz = find(v(:,3) == lz);
xz = v(indz,1); yz = v(indz,2); vz = [xz, yz];
tz = delaunay(xz,yz);  % triangulation de Delaunay des points du plan 

 %figure();
 %triplot(tz,xz,yz);   % affichage
 %title('maillage cuve initial')

%==========================================================================
% Calcul du champ magnetique et de la force avec le maillage 2D defini par
% (tz,vz)
%==========================================================================

%coordonnees des differents maillage
node_ensemble = [0,0 ; 10 0; 10 10; 0 10]; %  maillage total
node_select = [4.5 4.5;5.5 4.5; 5.5 5.5; 4.5 5.5];        %  mailage cuve cote 4
node_aimant = [4.75 3.5; 5.25 3.5; 5.25 4.5; 4.75 4.5];        %  maillage aimant

edge_aimant = [(1:size(node_aimant,1))',[(2:size(node_aimant,1))'; 1]];       %liste des aretes interieures
edge_ensemble = [1,2; 2,3; 3,4; 4,1];                               %liste des aretes exterieures

edge_total = [edge_aimant; edge_ensemble+size(node_aimant,1)];        %liste de toutes les aretes mises dans l'ordre
node_total = [node_aimant; node_ensemble];                      %liste de tous les noeuds mises dans l'ordre


%generation maillage total de resolution du probleme 2D
L_total = 10;
H_total = 10;
x_2D = linspace(0,L_total,2*10*(Nx-1) + 1)';
y_2D = linspace(0,H_total,2*10*(Ny-1)+1 )';

[X_total_2D,Y_total_2D] = meshgrid(x_2D,y_2D);
v_total_2D_x = X_total_2D';
v_total_2D_x = v_total_2D_x(:);
v_total_2D_y = Y_total_2D';
v_total_2D_y = v_total_2D_y(:);
v_total_2D = [v_total_2D_x(:) v_total_2D_y(:)];

t_total_2D = delaunay(v_total_2D(:,1),v_total_2D(:,2)); %genration maillage triangulaire total

 %figure();
 %triplot(t_total_2D ,X_total_2D,Y_total_2D);
 %xlim([4.5 5.5]);
 %ylim([4.5 5.5]);
 %title('maillage de resolution 2D total')


%u = zeros(size(X_total_2D(:),1),1);

%  figure();
%  patch('faces',t_total_2D(:,1:3),'vertices',v_total_2D(:,1:2),'facecolor','white');
%  title('avec path')
%  xlim([4.5 5.5]);
%  ylim([4.5 5.5]);

[B_X,B_Y] = champ_magnetique_fct_Tancrede_1(v_total_2D,t_total_2D,node_aimant,edge_aimant,node_ensemble,edge_ensemble);

f0 = B_Y; %force magnetique pour la couche du maillage a z =0

% Construction des forces 3D : fx(x,y,z)=fx(x,y,z=0)
fx = zeros(size(v,1),1);
fy = zeros(size(v,1),1);
fz = zeros(size(v,1),1);

%==========================================================================
%                 placement du second menbre dans le maillage 3D
%==========================================================================

% placement des valeurs sur une matrice colonne
for k=1:size(vz,1) 
    iz = find(vz(k,1) == v(:,1) & vz(k,2) == v(:,2)); %on recupere tout les noeuds 'au-dessus'
    fz(iz,:)=f0(k);                                     %les autres composantes son nulles ici
end
 %on verifie avec les tables de connectivité : les points sont biens placé
 
 
% placement des valuers sur une matrice Nx*Ny*Nz
% d'abord placement de f0 sur une matrice Nx*Ny
fxy=zeros(2*Ny-1,2*Nx-1);
for yi = 1:2*Ny-1       % selon l'axe y
    for xj = 1:2*Nx-1         % selon l'axe x
        m = (2*Nx-1)*(yi-1)+xj ;
        fxy(xj,yi) = f0(m);
    end
end

fxgd=zeros(2*Ny-1,2*Nx-1,2*Nz-1);
fygd=zeros(2*Ny-1,2*Nx-1,2*Nz-1);
fzgd=zeros(2*Ny-1,2*Nx-1,2*Nz-1);
% comme la valeur selon z est l seule non-nulle on ne touche pas aux autres
% composante ce qui serait des calculs inutiles
for m=1:(2*Nz-1)
    fzgd(:,:,m)=fxy;
end
    
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
