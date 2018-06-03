function [sol,iter]=uzawa(A,F,nv1,nv2)
% ----------------------------------------------
% ----------------------------------------------
% Resolution du systeme de Stokes par UZAWA
% ----------------------------------------------
    fprintf('*/ Résolution du système linéaire par UZAWA\n'); 

    % Parametres pour les iterations d'UZAWA
    r=2000; rho=2000;
    tol=1e-6;
    maxit=150;
    
    % Parametres pour MINRES
    tolMINRES=1e-10; itermaxMINRES=500; 

    % Matrice d'Uzawa

    F1=F(1:3*nv2); F2=F(3*nv2+1:end);
    A1=A(1:3*nv2,1:3*nv2);
    B=A(1:3*nv2,3*nv2+1:end);
    Ar=A1+r*B*(B');

    M1=diag(diag(Ar));
   
    % Initialisation de la pression
    p = zeros(nv1,1);

    % calcul de la vitesse : u = Ar\(F1-B*(p-r*F2)) 
    [u,flag,resM,iterM] = minres(Ar,F1-B*(p-r*F2),tolMINRES,itermaxMINRES,M1);
    fprintf('     residu relatif = %6.3e, nombre d''iteration = %d\n',resM,iterM);
    
    k=0;
    res=1;
    while res>tol && k<maxit

        % calcul de la pression
        p = p+rho*(B'*u-F2);

        % calcul de la vitesse : u = Ar\(F1-B*(p-r*F2)) 
        [v,flag,resM,iterM] = minres(Ar,F1-B*(p-r*F2),tolMINRES,itermaxMINRES,M1,[],u);
        fprintf('     residu relatif = %6.3e, nombre d''iteration = %d,',resM,iterM);
                
        res=norm(u-v,2);
        fprintf(' residu UZAWA  = %6.3e\n',res);
        u=v;
        k=k+1;
    end

    fprintf('nb iterations : %d | erreur : %6.3e\n', k, res);
    res1=max(abs(A1*u+B*p-F1));
    res2=max(abs(B'*u-F2));
    fprintf('Résidus : max|A*u+B*p-F1|=%6.3e\n',res1);
    fprintf('              max|B''u-F2|=%6.3e\n',res2);
    
    sol=[u;p];
    iter=k;
    
end
