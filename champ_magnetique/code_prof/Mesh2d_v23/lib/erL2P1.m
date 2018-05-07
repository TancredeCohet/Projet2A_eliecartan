function [erreurL2]=erL2P1(solexact,sol,t,v)
%--------------------------------------------------------------------------
% Calculs de la norme L2 de l'erreur
%   sol : solution approchée (P1)
%   solexact : fonction donnant la solution exacte
%--------------------------------------------------------------------------
    
    % Données pour la formule de quadrature à 7 points de Gauss (exacte P5)
    %
    a1 = (9-2*sqrt(15))/21;  % a1 = 0.059715871789770;
    b1 = (6+sqrt(15))/21;    % b1 = 0.470142064105115;
    a2 = (9+2*sqrt(15))/21;  % a2 = 0.797426985353087;
    b2 = (6-sqrt(15))/21;    % b2 = 0.101286507323456;
    % coordonnées barycentriques des points d'integration
    CBI = [1/3,1/3,1/3;
           a1, b1, b1;
           b1, a1, b1;
           b1, b1, a1;
           a2, b2, b2;
           b2, a2, b2;
           b2, b2, a2];
    %   Poids associes
    w1 = 9/40;                   % w1 = 0.225;
    w2 = (155+sqrt(15))/1200;    % w2 = 0.132394152788506;
    w3 = (155-sqrt(15))/1200;    % w3 = 0.125939180544827;
    w = [w1; w2; w2; w2; w3; w3; w3];
    %----------------------------------------------------------------------
    
    nt=size(t,1);
    sol=sol(:);      % vecteur colonne
    
    % Formule de quadrature 
    FQ=0;
    for k=1:nt
       
        %  Coordonnees des points d'integration dans le triangle k
        vertices=[v(t(k,1),1:2); v(t(k,2),1:2); v(t(k,3),1:2)];
        INTpts = CBI*vertices;
        
        % solution exacte aux points d'integration
        %s1=solexact(INTpts(:,1), INTpts(:,2)); s1=s1(:);
        s1=CBI*solexact(vertices(:,1),vertices(:,2));     
        
        % solution approchée (P1) aux points d'integration
        ndk=t(k,1:3);
        s2=CBI*sol(ndk);
        
        area=t(k,5);
        FQ=FQ+area*w'*((s1-s2).^2);
    end
    erreurL2=sqrt(FQ);
    
%     %----------------------------------------------------------------------
%     % Formule de quadrature des points milieux (exacte P2)
%     stm=solexact(v2(:,1), v2(:,2));
%     S2=0;
%     for k=1:nt1
%         %  Coordonnees des 3 milieux des aretes du triangle k
%         area=t1(k,5);
%         n4=t2(k+3*nt1,1); n5=t2(k+3*nt1,2); n6=t2(k+3*nt1,3);   
%         S2=S2+area/3 * sum((stm([n4,n5,n6])-sol([n4,n5,n6])).^2);
%     end
%     erS2=sqrt(S2);
    
end

