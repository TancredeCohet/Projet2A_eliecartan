%==========================================================================
%                               extraction du millage 2D a partir du
%                               maillage 3D
%==========================================================================

addpath './lib_ALICE/mexfiles/bin'
addpath './lib'
addpath('../champ_magnetique');

%   parametre de la cuve
global L P H
L = 12;  % largeur
P = 12;  % profondeur
H = 12;  % hauteur
nu = 1; % viscosite

% Nombres de points dans les directions x, y et z pour la cuve
Nx = 4*2;
Ny = 4*2;
Nz = 4*2;

%on veut un aimant de taille 0.5 * 1
taille_aimant = [1,1]; %aimant de taille 0.5 * 1 * inf

%generation du maillage 3D
[v,t,nv1,nv2,nbquad,X,Y,Z]= maillage_3D(Nx,Ny,Nz,L,P,H);

fprintf('*/ Maillage du cube %dx%dx%d (%f s)\n',Nx,Ny,Nz); 
fprintf('     nombre d''elements (cubes) = %d\n', size(t,1));
fprintf('     nombre de d.o.f pour la vitesse = %d\n', nv2);
fprintf('     nombre de d.o.f pour la pression = %d\n', nv1);

% Construction du maillage triangulaire 2D pour z=0 a  partir du maillage 3D
% par des cubes.
lz = 0;                % niveau z=0
indz = find(v(:,3) == lz);
xz = v(indz,1); yz = v(indz,2); vz = [xz, yz];
tz = delaunay(xz,yz);  % triangulation de Delaunay des points du plan 

triplot(tz,xz,yz);   % affichage
title('maillage cuve')

%==========================================================================
% Calcul du champ magnetique et de la force avec le maillage 2D defini par
% (tz,vz)
%==========================================================================

%coordonnées des differents maillage
node_ensemble = [0,0 ; 10 0; 10 10; 0 10]; %  maillage total
node_select = [3 4; 7 4; 7 8; 3 8];        %  mailage cuve cote 4
node_aimant = [5 2; 7 2; 7 5; 5 5];        %  maillage aimant

edge_aimant = [(1:size(node_aimant,1))',[(2:size(node_aimant,1))'; 1]];       %liste des aretes interieures
edge_ensemble = [1,2; 2,3; 3,4; 4,1];                               %liste des aretes exterieures

edge_total = [edge_aimant; edge_ensemble+size(node_aimant,1)];        %liste de toutes les aretes mises dans l'ordre
node_total = [node_aimant; node_ensemble];                      %liste de tous les noeuds mises dans l'ordre


%genberation maillage total de resolution du probleme 2D
x_2D = linspace(0,L,3*Nx + 1)';
y_2D = linspace(0,P,3*Ny + 1)';
[X_total_2D,Y_total_2D] = meshgrid(x_2D,y_2D);

%on remet dans le bon ordre les coordonnées des nouds pour le vecteur v
for i = 0:size(X_total_2D,1)-1
    v_total_2D_x(:,i + 1) = X_total_2D(size(X_total_2D,1) - i,:)';
    v_total_2D_y(:,i + 1) = Y_total_2D(size(X_total_2D,1) - i,:)';
end
v_total_2D = [v_total_2D_x(:) v_total_2D_y(:)];

t_total_2D = delaunay(v_total_2D(1,:),v_total_2D(2,:)); %genration maillage triangulaire total


triplot(t_total_2D ,X_total_2D,Y_total_2D); 
title('maillage de resolution 2D total')


u = zeros(size(X_total_2D(:),1),1);
figure()
patch('faces',t_total_2D(:,1:3),'vertices',v_total_2D(:,1:2),'FaceVertexCData',u,'facecolor','interp','edgecolor','black');
title('avec path')
figure()
[B_X,B_Y] = champ_magnetique_fct_Tancrede(v_total_2D,t_total_2D,node_aimant,edge_aimant,node_total,edge_total);

f0 = B_Y; %force magnetique pour la couche du maillage a z =0

% Construction des forces 3D : fx(x,y,z)=fx(x,y,z=0)
fx = zeros(size(v,1),1);
fy = zeros(size(v,1),1);
fz = zeros(size(v,1),1);

for k=1:size(v,1) 
    iz = find(vz(:,1) == v(k,1) & vz(:,2) == v(k,2)); %on recupere tout les noeuds 'au-dessus'
    fx(k)=f0(iz);                                     %les autres composantes son nulles ici
end

    
% Affichage de fx, juste pour vesrifier
% Construction du maillage du cube [0,L]x[|O,P]x[0,H]
x_2D = linspace(0,L,2*Nx-1);
y_2D = linspace(0,P,2*Ny-1);
zl = linspace(0,H,2*Nz-1);
[X,Y,Z]=meshgrid(x_2D,y_2D,zl);

hx=L/(2*(Nx-1)); hy=P/(2*(Ny-1)); hz=H/(2*(Nz-1));
fxgd=zeros(2*Ny-1,2*Nx-1,2*Nz-1);

for m=1:size(v,1)
    ind=round(v(m,:)./[hx, hy, hz])+1;
    fxgd(ind(2),ind(1),ind(3))=fx(m);
end
figure()
xslice = [0,L/2,L]; yslice = P/2; zslice = [0,H/2,H];
slice(X,Y,Z,fxgd,xslice,yslice,zslice,'linear')
xlabel('X'); ylabel('Y'); zlabel('Z');
