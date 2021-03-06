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
[fx, fy, fz, fxgd, fygd,fzgd,fxy] = f2Dto3D_2(v,t,nv1,nv2,nbquad,X,Y,Z,Nx,Ny,Nz); 
%   fx, fy, fz sont les forces dans des vecteurs colonnes 1*size(v,1)
%   fxgd, fygd, fzgd sont les forces dans des matrices Nx*Ny*Ny pour
%   affichage


% terme source f sans pression

F=[fx; fy; fz; sparse(nv1,1)];

% Traitement de la pression à moyenne nulle par multiplicateur de Lagrange.
[A,F,nv1] = tmtpression(Nx,Ny,Nz,nv1,nv2,t,A,F);

% Traitement des CL
ibT=[ind_bd;ind_bd+nv2;ind_bd+2*nv2];
icT=(1:3*nv2+nv1);
icT(ibT)=[];
F(ibT)=[];
A = A(icT,icT);

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
[X2,Y2,Z2,Uxgd,Uygd,Uzgd,Pgd] = gridformat(Lx,Ly,Lz,Nx,Ny,Nz,v,Ux,Uy,Uz,Pr2);

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
   
    
%   affichage du second membre de navier Stokes
    figure();
    quiver3(X2,Y2,Z2,fxgd,fygd,fzgd);
    title('champ de vitesse');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal;

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

