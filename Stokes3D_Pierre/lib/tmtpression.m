function [A,F,nv1]=tmtpression(Nx,Ny,Nz,nv1,nv2,t,A,F,tmtopt);
%
% Traitement la pression :
% tmtopt = 1 -> par multiplicateur de Lagrange : pression à moyenne nulle
%          2 -> pression fixée en un point
%
global H L P

if tmtopt==1
    % -1- Par multiplicateur de lagrange : pression à moyenne nulle
    un=ones(nv1,1);
    w=zeros(nv1,1);
    nbquad=(Nx-1)*(Ny-1)*(Nz-1);
    volume=L*H*P/nbquad;
    %(1/(Nx-1))*(1/(Ny-1))*(1/(Nz-1));
    for k=1:nbquad
        pt=t(k,1:8);
        w(pt)=w(pt)+volume/8;
    end
    C=[sparse(3*nv2,1);w/(w'*un)];

    A=[A C;
       C' 0];
    F=[F;0];
    nv1=nv1+1;
else
    % -2- Pression fixée en un point
    A(3*nv2+1,:)=[];
    A(:,3*nv2+1)=[];
    F(3*nv2+1)=[];
    nv1=nv1-1;
end