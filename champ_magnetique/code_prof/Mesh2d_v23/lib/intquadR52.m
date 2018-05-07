function F = intquadR52(f,t1,v1,t2,v2)
% Calcul de \int_\Omega f(x)l_j(x)dx où l_j fonctions de base P2 sur chaque 
% triangle par la formule de quadrature R52 exact P5.
%--------------------------------------------------------------------------
lambda = [1/3, 1/3, 1/3;
            1,   0,   0;
            0,   1,   0;
            0,   0,   1;
          1/2, 1/2,   0;
            0, 1/2, 1/2;
          1/2,   0, 1/2;
          5/7, 1/7, 1/7;
          1/7, 5/7, 1/7;
          1/7, 1/7, 5/7];
      
weight = [81/320, 1/90, 1/90, 1/90, 16/225, 16/225, 16/225, 2401/14400, 2401/14400, 2401/14400];

%--------------------------------------------------------------------------

fbP2 = [lambda(:,1:3).*(2*lambda(:,1:3)-1), ...
        4*lambda(:,2).*lambda(:,3), ...
        4*lambda(:,1).*lambda(:,3), ...
        4*lambda(:,1).*lambda(:,2)]';
    
W = (ones(6,1)*weight).*fbP2;

nv2=size(v2,1);
nt1=size(t1,1);
F=zeros(nv2,1);
for k=1:nt1

    i1=t1(k,1:3);
    p1=v1(i1,1:2);
    aire=t1(k,5);
    
    nx = lambda*p1; % coordonnées des points d'integration (10x2)
    
    i2=[i1, t2(k+3*nt1,1:3)];
    F(i2)=F(i2)+aire*W*f(nx(:,1), nx(:,2));
    
end

end

