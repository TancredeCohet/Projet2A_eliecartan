%=========================================================================% 
%       Resulution finale du probleme avec superposition des deux
%       tourbillons
%=========================================================================%

%   resolution dans la configuration aimant centre
%[X2,Y2,Z2,Uxgd_centre,Uygd_centre,Uzgd_centre,Pgd_centre] = stokes_3D(true);

%   resolution dans la configuration aimant de cote
[X2,Y2,Z2,Uxgd_cote,Uygd_cote,Uzgd_cote,Pgd_cote] = stokes_3D(false);

%   transposition du champ pour l'aimant de coté et somme des deux champ de
%   vitesse
[Uxgd_cote_transp,Uygd_cote_transp,Uzgd_cote_transp] = transposition_champ(Uxgd_cote,Uygd_cote,Uzgd_cote);
% Sauvegarde de la vitesse et de la pression
filename1_vtk='transpose.vtk';
fprintf('   - au format VTK (''VTK STRUCTURED_GRID'') : %s\n',filename1_vtk);
writevtk3D(X2,Y2,Z2,Uxgd_cote_transp,Uygd_cote_transp,Uzgd_cote_transp,Pgd_cote,filename1_vtk);

%   somme des 2 champs
Uxgd_total = Uxgd_centre + Uxgd_cote_transp;
Uygd_total = Uygd_centre + Uygd_cote_transp;
Uzgd_total = Uzgd_centre + Uzgd_cote_transp;
Pgd_total = Pgd_centre + Pgd_cote_transp;

filename2_vtk='total.vtk';
fprintf('   - au format VTK (''VTK STRUCTURED_GRID'') : %s\n',filename2_vtk);
writevtk3D(X2,Y2,Z2,Uxgd_total,Uygd_total,Uzgd_total,Pgd_total,filename2_vtk);

