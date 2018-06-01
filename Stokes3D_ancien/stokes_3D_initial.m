%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations de Stokes 3D dans le cube [0,L]x[0,P]x[0,H].
% Conditions limites nulles pour la vitesse sur le bord
% Elements Finis Q2/Q1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

addpath './lib_ALICE/mexfiles/bin'
addpath './lib'
tstart0=tic;
ClearSparseMatrix

global L P H
L=1;  % largeur
P=1;  % profondeur
H=1;  % hauteur
nu=1; % viscosité

% Nombres de points dans les directions x, y et z
Nx=10;
Ny=10;
Nz=10;

F1=1000;       % amplitude de la force pour 1 tourbillon
F2=5000;       % amplitude de la force pour 2 tourbillons

compare = false; % comparaison avec une solution exacte

fprintf('********************************************\n');
fprintf(' Stokes 3D dans le cube [0,L]x[0,P]x[0,H]   \n');
fprintf('********************************************\n');
%------------------------------------------------
% Maillage (grid mesh) du cube [0,L]x[0,P]x[0,H]
%------------------------------------------------
tstart=tic;
[v,t,nv1,nv2,nbquad,X,Y,Z]=maillage_3D(Nx,Ny,Nz,L,P,H);
fprintf('*/ Maillage du cube %dx%dx%d (%f s)\n',Nx,Ny,Nz,toc(tstart)); 
fprintf('     nombre d''elements (cubes) = %d\n', size(t,1));
fprintf('     nombre de d.o.f pour la vitesse = %d\n', nv2);
fprintf('     nombre de d.o.f pour la pression = %d\n', nv1);

% Recherche des points du bord pour le traitement des conditions aux 
% limites de Dirichlet
ind_bd=find(v(:,1)==0|v(:,1)==L|v(:,2)==0|v(:,2)==P|v(:,3)==0|v(:,3)==H);
ind_bd=unique(ind_bd);

%----------------------------------------------------------
% Construction des matrices de Stokes (format creux ALICE)
%   Asp = [ a  b
%           b' 0]
%   Msp mass matrix
%---------------------------------------------------------
tstart=tic;
fprintf('*/ Construction des matrices de Stokes\n');
[Asp,Msp]=matrice_3D(v,t,nv1,nv2,nu);
fprintf(' (%f s)\n',toc(tstart)); 

% Construction du second membre (RHS)

x=v(:,1); y=v(:,2); z=v(:,3);  % coordonnées

% Terme source f
arg=[L,P,H,nu,F1,F2];
[fx,fy,fz]=source_f(x,y,z,arg);

Masse=SparseMatrixToMatlab(Msp);
F=[Masse*fx; Masse*fy; Masse*fz; sparse(nv1,1)];
DeleteSparseMatrix(Msp);

% Traitement des conditions limites de Dirichlet : 
% ne peut traiter que le cas homogène u=0 sur le bord
tstart=tic;
ind=[ind_bd;ind_bd+nv2;ind_bd+2*nv2];
F=RemoveRowsColsSparseMatrix(Asp,ind,F);
fprintf('     traitement des CL (%f s)\n',toc(tstart)); 

tstart=tic;
A=SparseMatrixToMatlab(Asp);
fprintf('     passage du format ALICE->MATLAB (%f s)\n',toc(tstart));
DeleteSparseMatrix(Asp);

% Traitement de la pression 
% tmtopt = 1 -> multiplicateur de Lagrange (moyenne nulle)
%          2 -> pression fixee en un point
tmtopt = 1; 
[A,F,nv1]=tmtpression(Nx,Ny,Nz,nv1,nv2,t,A,F,tmtopt);
fprintf('     taille %d x %d, nnz = %d\n',size(A,1),size(A,2), nnz(A)); 

fprintf('     temps total de preparation du systeme lineaire : %4.2fs\n',toc(tstart0));

% -----------------------------------
% Resolution du systeme lineaire Ax=F
% -----------------------------------
[sol,residu,iter]=linsolver(A,F,nv1,nv2);
%[sol,iter]=uzawa(A,F,nv1,nv2);

% Vitesse
Ux=sol(1:nv2); Uy=sol(nv2+1:2*nv2); Uz=sol(2*nv2+1:3*nv2);
% Pression
Pr=sol((3*nv2+1):(end-1));

fprintf('*/ Sauvegarde (vitesse/pression)\n');
%----------------------------------------------------------------------
% Construction du format VTK et sauvegarde
%----------------------------------------------------------------------
% Interpolation de la pression sur le maillage Q2
nt1= (Nx-1)*(Ny-1)*(Nz-1); % nb de cubes Q1
Pr2 = build_Q2Pressure(t(1:nt1,:),Pr,nv2);

% Construction du maillage Q1-isoQ2
tQ1isoQ2 = meshQ1isoQ2(t); 

% Sauvegarde de la vitesse et de la pression
filename_vtk='UP3D.vtk';
fprintf('   - au format VTK (''VTK UnstructuredGrid Files'') : %s\n',filename_vtk);
write_vtu_unstructured(v,tQ1isoQ2,Ux,Uy,Uz,Pr2,filename_vtk);
%write_vtu_unstructured(v,t,Ux,Uy,Uz,Pr2,filename);

[X2,Y2,Z2,Uxgd,Uygd,Uzgd,Pgd]=gridformat(L,P,H,Nx,Ny,Nz,v,Ux,Uy,Uz,Pr2);
filename_mat='UP3D.mat';
fprintf('   - au format MATLAB : %s\n',filename_mat);
save(filename_mat,'Uxgd','Uygd','Uzgd','Pgd','X2','Y2','Z2','F1','F2');

%-----------
% Affichage
%-----------
% module de la vitesse
U_mod=sqrt(Uxgd.^2 + Uygd.^2 + Uzgd.^2);
figure(1);
xslice = [0,L/2,L]; yslice = P/2; zslice = [0,H/2,H];
slice(X2,Y2,Z2,U_mod,xslice,yslice,zslice,'linear')
colorbar; title('Solution approchée');
xlabel('X'); ylabel('Y'); zlabel('Z');

%---------------------------------------------------------------
% Comparaison de la solution approchée avec une solution exacte
%---------------------------------------------------------------
if compare
    % Solution exacte 3D
    ux_exact=(x .^ 3 .* (1 - x) .^ 3 .* y .^ 2 .* (1 - y) .^ 3 .* z .^ 3 .* (1 - z) .^ 3) ./ 0.243e3 - (x .^ 3 .* (1 - x) .^ 3 .* y .^ 3 .* (1 - y) .^ 2 .* z .^ 3 .* (1 - z) .^ 3) ./ 0.243e3 - (x .^ 3 .* (1 - x) .^ 3 .* y .^ 3 .* (1 - y) .^ 3 .* z .^ 2 .* (1 - z) .^ 3) ./ 0.243e3 + (x .^ 3 .* (1 - x) .^ 3 .* y .^ 3 .* (1 - y) .^ 3 .* z .^ 3 .* (1 - z) .^ 2) ./ 0.243e3;
    uy_exact=(x .^ 3 .* (1 - x) .^ 3 .* y .^ 3 .* (1 - y) .^ 3 .* z .^ 2 .* (1 - z) .^ 3) ./ 0.243e3 - (x .^ 3 .* (1 - x) .^ 3 .* y .^ 3 .* (1 - y) .^ 3 .* z .^ 3 .* (1 - z) .^ 2) ./ 0.243e3 - (x .^ 2 .* (1 - x) .^ 3 .* y .^ 3 .* (1 - y) .^ 3 .* z .^ 3 .* (1 - z) .^ 3) ./ 0.243e3 + (x .^ 3 .* (1 - x) .^ 2 .* y .^ 3 .* (1 - y) .^ 3 .* z .^ 3 .* (1 - z) .^ 3) ./ 0.243e3;
    uz_exact=(x .^ 2 .* (1 - x) .^ 3 .* y .^ 3 .* (1 - y) .^ 3 .* z .^ 3 .* (1 - z) .^ 3) ./ 0.243e3 - (x .^ 3 .* (1 - x) .^ 2 .* y .^ 3 .* (1 - y) .^ 3 .* z .^ 3 .* (1 - z) .^ 3) ./ 0.243e3 - (x .^ 3 .* (1 - x) .^ 3 .* y .^ 2 .* (1 - y) .^ 3 .* z .^ 3 .* (1 - z) .^ 3) ./ 0.243e3 + (x .^ 3 .* (1 - x) .^ 3 .* y .^ 3 .* (1 - y) .^ 2 .* z .^ 3 .* (1 - z) .^ 3) ./ 0.243e3;
    ux_exact=1e4*ux_exact; uy_exact=1e4*uy_exact; uz_exact=1e4*uz_exact;
    p_exact=zeros(nv1,1);

    sol_exact=[ux_exact;uy_exact;uz_exact;p_exact;0];
    sol_mod=sqrt(ux_exact.^2 + uy_exact.^2 + uz_exact.^2);

    R=reshape(sol_mod(1:Nx*Ny*Nz),Nx,Ny,Nz);
    R=permute(R,[2,1,3]);

    xslice = [0,L/2,L]; yslice = P/2; zslice = [0,H/2,H];
    %xslice = x(1:Nx); yslice = y(1:Ny); zslice = Z(1:Nz);
    figure(2)
    slice(X,Y,Z,R,xslice,yslice,zslice,'nearest')
    colorbar; title('Solution exacte');
    xlabel('X'); ylabel('Y'); zlabel('Z');

    U_mod=sqrt(Ux.^2 + Uy.^2 + Uz.^2);
    erreur = max(abs(sol_mod-U_mod))/max(abs(sol_mod));
    fprintf('*/ Comparaison avec la solution exacte :\n');
    fprintf('   erreur relative = %6.3e\n',erreur);
end