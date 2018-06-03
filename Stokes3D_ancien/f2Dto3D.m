%==========================================================================
%                               extraction du millage 2D a partir du
%                               maillage 3D
%==========================================================================

[v,t,nv1,nv2,nbquad,X,Y,Z]= maillage_3D(Nx,Ny,Nz,L,P,H);

fprintf('*/ Maillage du cube %dx%dx%d (%f s)\n',Nx,Ny,Nz); 
fprintf('     nombre d''elements (cubes) = %d\n', size(t,1));
fprintf('     nombre de d.o.f pour la vitesse = %d\n', nv2);
fprintf('     nombre de d.o.f pour la pression = %d\n', nv1);

% Construction du maillage triangulaire 2D pour z=0 à partir du maillage 3D
% par des cubes.
lz=0;                % niveau z=0
indz = find(v(:,3)==lz);
xz=v(indz,1); yz=v(indz,2); vz=[xz, yz];
tz=delaunay(xz,yz);  % triangulation de Delaunay des points du plan 
triplot(tz,xz,yz);   % affichage

% Calcul du champ magnétique et de la force avec le maillage 2D défini par
% (tz,vz)
f0 = xz+yz;  % c'est un exemple scalaire, juste pour illustrer. A remplacer 
             % par le vrai calcul de la force/champ magnétique

% Construction des forces 3D : fx(x,y,z)=fx(x,y,z=0)
fx=zeros(size(v,1),1);
for k=1:size(v,1)
    iz=find(vz(:,1)==v(k,1) & vz(:,2)==v(k,2)); 
    fx(k)=f0(iz);
end
    
% Affichage de fx, juste pour vérifier
% Construction du maillage du cube [0,L]x[|O,P]x[0,H]
xl = linspace(0,L,2*Nx-1);
yl = linspace(0,P,2*Ny-1);
zl = linspace(0,H,2*Nz-1);
[X,Y,Z]=meshgrid(xl,yl,zl);

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
