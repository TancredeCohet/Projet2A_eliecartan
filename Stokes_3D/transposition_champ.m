function [Ux_transp,Uy_transp,Uz_transp] = transposition_champ(Ux,Uy,Uz)
%-------------------------------------------------------------------------%
%           transpose le champ par couche de telles maniere 
%               a pouvoir sommer les champs de vitesses
%-------------------------------------------------------------------------%
Ux_transp = zeros(size(Ux));
Uy_transp = zeros(size(Uy));
Uz_transp = zeros(size(Uz));
%------------ pour Uz -----------------------------------------------------
for z = 1:size(Uz,3)
    for x = 1:size(Uz,1)
        Uz_transp(:,size(Uz,1) + 1 - x,z) = Uz(x,:,z)'; 
    end
end
%------------ pour Ux -----------------------------------------------------
for z = 1:size(Ux,3)
    for x = 1:size(Uz,1)
        Ux_transp(:,size(Uz,1) + 1 - x,z) = - Uy(x,:,z)'; 
    end
end
%------------ pour Uy -----------------------------------------------------
for z = 1:size(Uy,3)
    for x = 1:size(Uz,1)
        Uy_transp(:,size(Uz,1) + 1 - x,z) = Ux(x,:,z)'; 
    end
end