%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations de Stokes 3D dans le cube [0,Lx]x[0,Ly]x[0,Lz]
% Elements Finis Q2/Q1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

addpath './lib'

tstart0=tic;

global Lx Ly Lz
Lx=1;  % largeur
Ly=1;  % profondeur
Lz=1;  % hauteur
nu=1; % viscosité

% Nombres de points dans les directions x, y et z
Nx=10;
Ny=10;
Nz=10;

compare = true; % comparaison avec une solution exacte
displaysol = true;  % affichage de la solution approchée (vitesse)

fprintf('********************************************\n');
fprintf(' Stokes 3D dans le cube [0,Lx]x[0,Ly]x[0,Lz]   \n');
fprintf('********************************************\n');
%--------------------------------------------------------------------------
% Maillage (grid mesh) du cube [0,Lx]x[0,Ly]x[0,Lz]
%--------------------------------------------------------------------------
tstart=tic;
fprintf('*/ Maillage du cube %dx%dx%d : ',Nx,Ny,Nz); 

[v,t,nv1,nv2,nbquad,X,Y,Z] = mesh_3D(Lx,Ly,Lz,Nx,Ny,Nz);
fprintf('%g s\n',toc(tstart));
fprintf('     nombre d''elements (cubes) = %d\n', size(t,1));
fprintf('     nombre de noeuds pour la vitesse = %d\n', nv2);
fprintf('     nombre de noeuds pour la pression = %d\n', nv1);

% Recherche des points du bord pour le traitement des conditions aux 
% limites de Dirichlet
ind_bd=find(v(:,1)==0|v(:,1)==Lx|v(:,2)==0|v(:,2)==Ly|v(:,3)==0|v(:,3)==Lz);
ind_bd=unique(ind_bd);

%--------------------------------------------------------------------------
% Construction des matrices de Stokes 
%   A = [ a  b
%         b' 0]
%   Masse : matrice de masse
%--------------------------------------------------------------------------
tstart1=tic;
fprintf('*/ Assemblage des matrices : ');
[DxDx,DyDy,DzDz,DxDy,DxDz,DyDz,DxW,DyW,DzW,Masse]=matrix_3D(v,t);

A=[nu*(2*DxDx+DyDy+DzDz),              nu*DxDy',              nu*DxDz',             -DxW;
                  nu*DxDy, nu*(DxDx+2*DyDy+DzDz),              nu*DyDz',             -DyW;
                  nu*DxDz,               nu*DyDz, nu*(DxDx+DyDy+2*DzDz),             -DzW;
                    -DxW',                 -DyW',                 -DzW', sparse(nv1,nv1)];
fprintf('%g s\n',toc(tstart1)); 

% Construction du second membre (RHS)
x=v(:,1); y=v(:,2); z=v(:,3);  % coordonnées

% terme source f sans pression
[fx,fy,fz]=source_f(nu,x,y,z);
F=[Masse*fx; Masse*fy; Masse*fz; sparse(nv1,1)];

% Traitement de la pression à moyenne nulle par multiplicateur de Lagrange.
[A,F,nv1]=tmtpression(Nx,Ny,Nz,nv1,nv2,t,A,F);

% Traitement des CL
ibT=[ind_bd;ind_bd+nv2;ind_bd+2*nv2];
icT=(1:3*nv2+nv1);
icT(ibT)=[];
F(ibT)=[];
A=A(icT,icT);

fprintf('     taille %d x %d, nnz = %d\n',size(A,1),size(A,2), nnz(A)); 
fprintf('     temps total de preparation du systeme lineaire : %g s\n',toc(tstart0));

% -------------------------------------------------------------------------
% Resolution du systeme lineaire Ax=F
% -------------------------------------------------------------------------
[sol,residu,iter]=linsolver(A,F,nv1,nv2-length(ind_bd));
%[sol,iter]=uzawa(A,F,nv1,nv2-length(ind_bd));

UP=zeros(3*nv2+nv1-1,1);
UP(icT)=sol;
% Vitesse
Ux=UP(1:nv2); Uy=UP(nv2+1:2*nv2); Uz=UP(2*nv2+1:3*nv2);
% Pression
Pr=UP((3*nv2+1):end);

fprintf('---------------------------------\n');
fprintf('Temps total de calcul : %g s\n',toc(tstart0));
fprintf('---------------------------------\n');
fprintf('*/ Sauvegarde (vitesse/pression)\n');
%--------------------------------------------------------------------------
% Construction du format VTK et sauvegarde
%--------------------------------------------------------------------------
% Interpolation de la pression sur le maillage Q2
nt1= (Nx-1)*(Ny-1)*(Nz-1); % nb de cubes Q1
Pr2 = build_Q2Pressure(t(1:nt1,:),Pr,nv2);

% Construction du maillage Q1-isoQ2
tQ1isoQ2 = meshQ1isoQ2(t); 
[X2,Y2,Z2,Uxgd,Uygd,Uzgd,Pgd]=gridformat(Lx,Ly,Lz,Nx,Ny,Nz,v,Ux,Uy,Uz,Pr2);

% Sauvegarde de la vitesse et de la pression
filename_vtk='UP3D.vtk';
fprintf('   - au format VTK (''VTK STRUCTURED_GRID'') : %s\n',filename_vtk);
writevtk3D(X2,Y2,Z2,Uxgd,Uygd,Uzgd,Pgd,filename_vtk);

filename_mat='UP3D.mat';
fprintf('   - au format MATLAB : %s\n',filename_mat);
save(filename_mat,'Uxgd','Uygd','Uzgd','Pgd','X2','Y2','Z2');

%--------------------------------------------------------------------------
% Affichage
%--------------------------------------------------------------------------
U_mod=sqrt(Uxgd.^2 + Uygd.^2 + Uzgd.^2);  % module de la vitesse

if displaysol
    
    % % Affichage Q1 du module de la vitesse
    % U_mod=sqrt(Ux.^2 + Uy.^2+ Uz.^2);
    % R=reshape(full(U_mod(1:Nx*Ny*Nz)),Nx,Ny,Nz);
    % xslice = [0,Lx/2,Lx]; yslice = Ly/2; zslice = [0,Lz/2,Lz];
    % figure;
    % slice(X,Y,Z,R,xslice,yslice,zslice,'nearest')

    % Affichage Q2 du module de la vitesse
    U_mod=sqrt(Uxgd.^2 + Uygd.^2 + Uzgd.^2);
    figure()
    xslice = [0,Lx/2,Lx]; yslice = Ly/2; zslice = [0,Lz/2,Lz];
    slice(X2,Y2,Z2,U_mod,xslice,yslice,zslice,'linear')

    colorbar; title('Solution approchée');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal

    figure()
    quiver3(X2,Y2,Z2,Uxgd,Uygd,Uzgd)
    title('champ de vitesse')
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal

end
%---------------------------------------------------------------
% Comparaison de la solution approchée avec une solution exacte
%---------------------------------------------------------------
if compare
    % Solution exacte 3D
    ux_exact=(x .^ 3 .* (Lx-x) .^ 3 .* y .^ 2 .* (Ly-y) .^ 3 .* z .^ 3 .* (Lz-z) .^ 3) ./ 0.243e3 - (x .^ 3 .* (Lx-x) .^ 3 .* y .^ 3 .* (Ly-y) .^ 2 .* z .^ 3 .* (Lz-z) .^ 3) ./ 0.243e3 - (x .^ 3 .* (Lx-x) .^ 3 .* y .^ 3 .* (Ly-y) .^ 3 .* z .^ 2 .* (Lz-z) .^ 3) ./ 0.243e3 + (x .^ 3 .* (Lx-x) .^ 3 .* y .^ 3 .* (Ly-y) .^ 3 .* z .^ 3 .* (Lz-z) .^ 2) ./ 0.243e3;
    uy_exact=(x .^ 3 .* (Lx-x) .^ 3 .* y .^ 3 .* (Ly-y) .^ 3 .* z .^ 2 .* (Lz-z) .^ 3) ./ 0.243e3 - (x .^ 3 .* (Lx-x) .^ 3 .* y .^ 3 .* (Ly-y) .^ 3 .* z .^ 3 .* (Lz-z) .^ 2) ./ 0.243e3 - (x .^ 2 .* (Lx-x) .^ 3 .* y .^ 3 .* (Ly-y) .^ 3 .* z .^ 3 .* (Lz-z) .^ 3) ./ 0.243e3 + (x .^ 3 .* (Lx-x) .^ 2 .* y .^ 3 .* (Ly-y) .^ 3 .* z .^ 3 .* (Lz-z) .^ 3) ./ 0.243e3;
    uz_exact=(x .^ 2 .* (Lx-x) .^ 3 .* y .^ 3 .* (Ly-y) .^ 3 .* z .^ 3 .* (Lz-z) .^ 3) ./ 0.243e3 - (x .^ 3 .* (Lx-x) .^ 2 .* y .^ 3 .* (Ly-y) .^ 3 .* z .^ 3 .* (Lz-z) .^ 3) ./ 0.243e3 - (x .^ 3 .* (Lx-x) .^ 3 .* y .^ 2 .* (Ly-y) .^ 3 .* z .^ 3 .* (Lz-z) .^ 3) ./ 0.243e3 + (x .^ 3 .* (Lx-x) .^ 3 .* y .^ 3 .* (Ly-y) .^ 2 .* z .^ 3 .* (Lz-z) .^ 3) ./ 0.243e3;
    ux_exact=1e4*ux_exact; uy_exact=1e4*uy_exact; uz_exact=1e4*uz_exact;
    p_exact=zeros(nv2,1);
    
    [~,~,~,Uxgd_e,Uygd_e,Uzgd_e,Pgd_e]=gridformat(Lx,Ly,Lz,Nx,Ny,Nz,v,ux_exact,uy_exact,uz_exact,p_exact);
    U_mod_exact=sqrt(Uxgd_e.^2 + Uygd_e.^2 + Uzgd_e.^2);
    erreur = max(abs(U_mod_exact(:)-U_mod(:)))/max(abs(U_mod_exact(:)));
    fprintf('*/ Comparaison avec la solution exacte :\n');
    fprintf('   erreur relative = %6.3e\n',erreur);
    
    % Affichage 
    if displaysol
        xslice = [0,Lx/2,Lx]; yslice = Ly/2; zslice = [0,Lz/2,Lz];
        figure()
        slice(X2,Y2,Z2,U_mod_exact,xslice,yslice,zslice,'linear')
        colorbar; title('Solution exacte');
        xlabel('X'); ylabel('Y'); zlabel('Z');
        axis equal
    end
end
