function [X,Y,Z,Uxgd,Uygd,Uzgd,Pgd]=gridformat(L,P,H,Nx,Ny,Nz,v,Ux,Uy,Uz,Pr)
%----------------------------------------------------------------------
% Construction du format meshgrid/MATLAB 
%----------------------------------------------------------------------
    % Construction du maillage du cube [0,L]x[|O,P]x[0,H]
    xl = linspace(0,L,2*Nx-1);
    yl = linspace(0,P,2*Ny-1);
    zl = linspace(0,H,2*Nz-1);
    [X,Y,Z]=meshgrid(xl,yl,zl);

    hx=L/(2*(Nx-1)); hy=P/(2*(Ny-1)); hz=H/(2*(Nz-1));
    Uxgd=zeros(2*Ny-1,2*Nx-1,2*Nz-1);
    Uygd=zeros(2*Ny-1,2*Nx-1,2*Nz-1);
    Uzgd=zeros(2*Ny-1,2*Nx-1,2*Nz-1);
    Pgd=zeros(2*Ny-1,2*Nx-1,2*Nz-1);

    for m=1:size(Ux,1)

        ind=round(v(m,:)./[hx, hy, hz])+1;
        Uxgd(ind(2),ind(1),ind(3))=Ux(m);
        Uygd(ind(2),ind(1),ind(3))=Uy(m);
        Uzgd(ind(2),ind(1),ind(3))=Uz(m);
        Pgd(ind(2),ind(1),ind(3))=Pr(m);
    end
    
end
