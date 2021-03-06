%=========================================================================% 
%       Resulution finale du probleme avec superposition des deux
%       tourbillons
%=========================================================================%

%   rapport R entre les forces magnetique
% R = f_1_touribillon / f_2_tourbillon
% on fait donc f_1 = R/R_t * f_2 o� R_t = norm(f_1)/norm(f_2);
%R = 0.005;
for R = [0.69 7.87 2 3 4 5 7 8]
%   resolution dans la configuration aimant centre
[X2,Y2,Z2,Uxgd_centre,Uygd_centre,Uzgd_centre,Pgd_centre] = stokes_3D(true,R);


%   resolution dans la configuration aimant de cote
[X2,Y2,Z2,Uxgd_cote,Uygd_cote,Uzgd_cote,Pgd_cote] = stokes_3D(false,R);

%   transposition du champ pour l'aimant de cot� et somme des deux champ de
%   vitesse
[Uxgd_cote_transp,Uygd_cote_transp,Uzgd_cote_transp] = transposition_champ(Uxgd_cote,Uygd_cote,Uzgd_cote);
% Sauvegarde de la vitesse et de la pression
% filename1_vtk='transpose.vtk';
% fprintf('   - au format VTK (''VTK STRUCTURED_GRID'') : %s\n',filename1_vtk);
% writevtk3D(X2,Y2,Z2,Uxgd_cote_transp,Uygd_cote_transp,Uzgd_cote_transp,Pgd_cote,filename1_vtk);

%   somme des 2 champs

% prise en compte du rapport R
R_i = norm(Uxgd_centre(:)) / norm(Uxgd_cote_transp(:));
Uxgd_cote_transp = (R / R_i) .* Uxgd_cote_transp;
Uxgd_total = Uxgd_centre + Uxgd_cote_transp;
Uygd_total = Uygd_centre + Uygd_cote_transp;
Uzgd_total = Uzgd_centre + Uzgd_cote_transp;
Pgd_total = Pgd_centre + Pgd_cote;

filename2_vtk=strcat('champ_total_R_',num2str(R),'.vtk');
fprintf('   - au format VTK (''VTK STRUCTURED_GRID'') : %s\n',filename2_vtk);
writevtk3D(X2,Y2,Z2,Uxgd_total,Uygd_total,Uzgd_total,Pgd_total,filename2_vtk);
end
