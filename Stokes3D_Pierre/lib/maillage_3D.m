function [v,t,nv1,nv2,nbquad,X,Y,Z]=maillage_3D(Nx,Ny,Nz,L,P,H)

% Construction du maillage du cube [0,L]x[|O,P]x[0,H]
xl = linspace(0,L,Nx);
yl = linspace(0,P,Ny);
zl = linspace(0,H,Nz);
[X,Y,Z]=meshgrid(xl,yl,zl);

nv1=Nx*Ny*Nz;                   % nb de sommets du maillage Q1
nv2=(2*Nx-1)*(2*Ny-1)*(2*Nz-1); % nb de sommets du maillage Q2
v = zeros(nv2,3);               % coordonnées des sommets du maillage Q2

nbquad=(Nx-1)*(Ny-1)*(Nz-1);    % nb de quadrangles du maillage Q1
t=zeros(nbquad,27);             % connectivité du maillage Q2

% Coordonnées des sommets du maillage Q1
k=1;
for m=1:Nz
    for j=1:Ny
        for i=1:Nx
            v(k,1)=X(j,i,m);
            v(k,2)=Y(j,i,m);
            v(k,3)=Z(j,i,m);      
            k=k+1;
        end
    end
end

% indice pour se deplacer dans le tableau v
ind_point=k;
% pour le passage à un plan en Z supérieur
num_plan=1; 
% indice pour se deplacer dans les ind_points premier points de v
k=1;

for i=1:nbquad
    
    % on se trouve au bout de la ligne des points en x
    if mod(k,Nx)==0 
        k=k+1;
    end
    
   % on arrive au bout du premier plan en Z 
    if mod(k,Nx*(Ny-1)+1+(num_plan-1)*Nx*Ny)==0 
        k=k+Nx;
        num_plan=num_plan+1;
    end

    % numerotation des huits sommets du quadrangle
    t(i,1)=k;
    t(i,2)=k+1;
    t(i,3)=k+Nx+1;
    t(i,4)=k+Nx;
    t(i,5)=k+Nx*Ny;
    t(i,6)=k+1+Nx*Ny;
    t(i,7)=k+Nx+1+Nx*Ny;
    t(i,8)=k+Nx+Nx*Ny;
    
    k=k+1;
    
    % numerotation des autres points milieux 
    %  traitement du point 9
    if k<=Nx
        v(ind_point,1)=(v(t(i,1),1)+v(t(i,2),1))/2;
        v(ind_point,2)=(v(t(i,1),2)+v(t(i,2),2))/2;
        v(ind_point,3)=(v(t(i,1),3)+v(t(i,2),3))/2;
        t(i,9)=ind_point;
        ind_point=ind_point+1;
    else
        % connectivité : pour le premier plan on va en 11, sinon en 17
        if (num_plan>1)
            t(i,9)=t(i-((Nx-1)*(Ny-1)),17);
        else
            t(i,9)=t(i-(Nx-1),11);
        end
    end
    
    % traitement des points 10, 11, 12, 21
    if (num_plan>1)
        %point 10
        t(i,10)=t(i-((Nx-1)*(Ny-1)),18);
        %point 11
        t(i,11)=t(i-((Nx-1)*(Ny-1)),19);
        %point 12
        t(i,12)=t(i-((Nx-1)*(Ny-1)),20);
        %point 21
        t(i,21)=t(i-((Nx-1)*(Ny-1)),26);
    else
        %point 10
        v(ind_point,1)=(v(t(i,2),1)+v(t(i,3),1))/2;
        v(ind_point,2)=(v(t(i,2),2)+v(t(i,3),2))/2;
        v(ind_point,3)=(v(t(i,2),3)+v(t(i,3),3))/2;
        t(i,10)=ind_point;
        ind_point=ind_point+1;
        %point 11
        v(ind_point,1)=(v(t(i,3),1)+v(t(i,4),1))/2;
        v(ind_point,2)=(v(t(i,3),2)+v(t(i,4),2))/2;
        v(ind_point,3)=(v(t(i,3),3)+v(t(i,4),3))/2;
        t(i,11)=ind_point;
        ind_point=ind_point+1;
        %point 12
        if mod(k-1,Nx)==1
            v(ind_point,1)=(v(t(i,4),1)+v(t(i,1),1))/2;
            v(ind_point,2)=(v(t(i,4),2)+v(t(i,1),2))/2;
            v(ind_point,3)=(v(t(i,4),3)+v(t(i,1),3))/2;
            t(i,12)=ind_point;
            ind_point=ind_point+1;
        else
            t(i,12)=t(i-1,10);
        end
        %point 21
        v(ind_point,1)=(v(t(i,1),1)+v(t(i,2),1))/2;
        v(ind_point,2)=(v(t(i,3),2)+v(t(i,2),2))/2;
        v(ind_point,3)=(v(t(i,1),3)+v(t(i,2),3))/2;
        t(i,21)=ind_point;
        ind_point=ind_point+1;
    end
    
    %point 13
    if (k-1)==(1+(num_plan-1)*Nx*Ny)
        v(ind_point,1)=(v(t(i,1),1)+v(t(i,5),1))/2;
        v(ind_point,2)=(v(t(i,1),2)+v(t(i,5),2))/2;
        v(ind_point,3)=(v(t(i,1),3)+v(t(i,5),3))/2;
        t(i,13)=ind_point;
        ind_point=ind_point+1;
    else
        if ((k-1)>(1+(num_plan-1)*Nx*Ny) && (k-1)<(Nx+(num_plan-1)*Nx*Ny))
            t(i,13)=t(i-1,14);
        else
            t(i,13)=t(i-(Nx-1),16);
        end
    end
    
    %point 14
    if ((k-1)>=(1+(num_plan-1)*Nx*Ny) && (k-1)<(Nx+(num_plan-1)*Nx*Ny))
        v(ind_point,1)=(v(t(i,2),1)+v(t(i,6),1))/2;
        v(ind_point,2)=(v(t(i,2),2)+v(t(i,6),2))/2;
        v(ind_point,3)=(v(t(i,2),3)+v(t(i,6),3))/2;
        t(i,14)=ind_point;
        ind_point=ind_point+1;
    else
        t(i,14)=t(i-(Nx-1),15);
    end
    
    %point 15
    v(ind_point,1)=(v(t(i,3),1)+v(t(i,7),1))/2;
    v(ind_point,2)=(v(t(i,3),2)+v(t(i,7),2))/2;
    v(ind_point,3)=(v(t(i,3),3)+v(t(i,7),3))/2;
    t(i,15)=ind_point;
    ind_point=ind_point+1;
    
    % traitement des points 16, 20, 25
    if mod((k-1)-(num_plan-1)*Nx*Ny,Nx)==1
        %point 16
        v(ind_point,1)=(v(t(i,4),1)+v(t(i,8),1))/2;
        v(ind_point,2)=(v(t(i,4),2)+v(t(i,8),2))/2;
        v(ind_point,3)=(v(t(i,4),3)+v(t(i,8),3))/2;
        t(i,16)=ind_point;
        ind_point=ind_point+1;
        %point 20
        v(ind_point,1)=(v(t(i,8),1)+v(t(i,5),1))/2;
        v(ind_point,2)=(v(t(i,8),2)+v(t(i,5),2))/2;
        v(ind_point,3)=(v(t(i,8),3)+v(t(i,5),3))/2;
        t(i,20)=ind_point;
        ind_point=ind_point+1;
        %point 25
        v(ind_point,1)=(v(t(i,8),1)+v(t(i,5),1))/2;
        v(ind_point,2)=(v(t(i,8),2)+v(t(i,5),2))/2;
        v(ind_point,3)=(v(t(i,1),3)+v(t(i,5),3))/2;
        t(i,25)=ind_point;
        ind_point=ind_point+1;
    else
        %point 16
        t(i,16)=t(i-1,15);
        %point 20
        t(i,20)=t(i-1,18);
        %point 25
        t(i,25)=t(i-1,23);
    end
    
    %point 17
    if (k-1)<=(Nx+(num_plan-1)*Nx*Ny)
        v(ind_point,1)=(v(t(i,5),1)+v(t(i,6),1))/2;
        v(ind_point,2)=(v(t(i,5),2)+v(t(i,6),2))/2;
        v(ind_point,3)=(v(t(i,5),3)+v(t(i,6),3))/2;
        t(i,17)=ind_point;
        ind_point=ind_point+1;
    else
        t(i,17)=t(i-(Nx-1),19);
    end
    
    %point 18
    v(ind_point,1)=(v(t(i,6),1)+v(t(i,7),1))/2;
    v(ind_point,2)=(v(t(i,6),2)+v(t(i,7),2))/2;
    v(ind_point,3)=(v(t(i,6),3)+v(t(i,7),3))/2;
    t(i,18)=ind_point;
    ind_point=ind_point+1;
    
    %point 19
    v(ind_point,1)=(v(t(i,7),1)+v(t(i,8),1))/2;
    v(ind_point,2)=(v(t(i,7),2)+v(t(i,8),2))/2;
    v(ind_point,3)=(v(t(i,7),3)+v(t(i,8),3))/2;
    t(i,19)=ind_point;
    ind_point=ind_point+1;
    
    %point 22
    if (k-1-(num_plan-1)*Nx*Ny)<Nx
        v(ind_point,1)=(v(t(i,1),1)+v(t(i,2),1))/2;
        v(ind_point,2)=(v(t(i,1),2)+v(t(i,5),2))/2;
        v(ind_point,3)=(v(t(i,1),3)+v(t(i,5),3))/2;
        t(i,22)=ind_point;
        ind_point=ind_point+1;
    else
        t(i,22)=t(i-(Nx-1),24);
    end
    
    %point 23
    v(ind_point,1)=(v(t(i,6),1)+v(t(i,2),1))/2;
    v(ind_point,2)=(v(t(i,2),2)+v(t(i,3),2))/2;
    v(ind_point,3)=(v(t(i,6),3)+v(t(i,2),3))/2;
    t(i,23)=ind_point;
    ind_point=ind_point+1;
    
    %point 24
    v(ind_point,1)=(v(t(i,3),1)+v(t(i,4),1))/2;
    v(ind_point,2)=(v(t(i,3),2)+v(t(i,4),2))/2;
    v(ind_point,3)=(v(t(i,4),3)+v(t(i,8),3))/2;
    t(i,24)=ind_point;
    ind_point=ind_point+1;
    
    %point 26
    v(ind_point,1)=(v(t(i,5),1)+v(t(i,6),1))/2;
    v(ind_point,2)=(v(t(i,5),2)+v(t(i,8),2))/2;
    v(ind_point,3)=(v(t(i,5),3)+v(t(i,8),3))/2;
    t(i,26)=ind_point;
    ind_point=ind_point+1;
    
    %point 27
    v(ind_point,1)=(v(t(i,1),1)+v(t(i,2),1))/2;
    v(ind_point,2)=(v(t(i,1),2)+v(t(i,4),2))/2;
    v(ind_point,3)=(v(t(i,1),3)+v(t(i,5),3))/2;
    t(i,27)=ind_point;
    ind_point=ind_point+1;

end
end