function [A,F,nv1]=tmtpression(Nx,Ny,Nz,nv1,nv2,t,A,F)
%
% Traitement de la pression à moyenne nulle par multiplicateur de Lagrange.
%
    global Lx Ly Lz

    un=ones(nv1,1);
    w=zeros(nv1,1);
    nbquad=(Nx-1)*(Ny-1)*(Nz-1);
    volume=Lx*Ly*Lz/nbquad;

    for k=1:nbquad
        pt=t(k,1:8);
        w(pt)=w(pt)+volume/8;
    end
    C=[sparse(3*nv2,1);w/(w'*un)];

    A=[A C;
       C' 0];
    F=[F;0];
    nv1=nv1+1;

end