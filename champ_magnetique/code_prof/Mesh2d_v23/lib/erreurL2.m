function [erU,erS]=erreurL2(solexact,PUx,PUy,Stream,t1,t2,v1,v2)

    global L H
    
    quadra = 1; % quadra=1 -> formule de quadrature de Gauss à 7 points
                %       =2 -> formule de quadrature des points milieux
                
    %----------------------------------------------------------------------
    % Données pour la formule de quadrature de Gauss à 7 points
    %
    a1 = (9-2*sqrt(15))/21;  % a1 = 0.059715871789770;
    b1 = (6+sqrt(15))/21;    % b1 = 0.470142064105115;
    a2 = (9+2*sqrt(15))/21;  % a2 = 0.797426985353087;
    b2 = (6-sqrt(15))/21;    % b2 = 0.101286507323456;
    MT7 = [1/3,1/3,1/3;
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
    
    nt1=size(t1,1);
    
    switch quadra,
        case 1,

            % Formule de quadrature à 7 points de Gauss (exacte P5)
            S1=0; S2=0;
            for k=1:nt1
                %  Coordonnees des 7 points de Gauss dans le triangle k
                area=t1(k,5);
                pts=[v1(t1(k,1),1:2); v1(t1(k,2),1:2); v1(t1(k,3),1:2)];
                M7pts = MT7*pts;
                [Ux,Uy,stm]=solexact(M7pts(:,1), M7pts(:,2),L,H);
                nn = [t2(k,1); t2(k+nt1,1); t2(k+2*nt1,1);
                      t2(k+3*nt1,1); t2(k+3*nt1,2); t2(k+3*nt1,3)]; 
                PUxq=PUx(nn)'*baseP2(MT7(:,1)',MT7(:,2)',MT7(:,3)');    
                PUyq=PUy(nn)'*baseP2(MT7(:,1)',MT7(:,2)',MT7(:,3)');
                S1=S1+area*w'*((Ux-PUxq').^2+(Uy-PUyq').^2);
                ss=Stream(nn)'*baseP2(MT7(:,1)',MT7(:,2)',MT7(:,3)');
                S2=S2+area*w'*((stm-ss').^2);
                
            end
            
        case 2,

            %----------------------------------------------------------------------
            % Formule de quadrature des points milieux (exacte P2)
            [Ux,Uy,stm]=solexact(v2(:,1), v2(:,2),L,H);
            S1=0; S2=0;
            %MT3=[0, 0.5, 0.5; 0.5, 0, 0.5; 0.5, 0.5, 0];
            for k=1:nt1
                %  Coordonnees des 3 milieux des aretes du triangle k
                area=t1(k,5);
                n4=t2(k+3*nt1,1); n5=t2(k+3*nt1,2); n6=t2(k+3*nt1,3);  
                S1=S1+area/3 * sum((Ux([n4,n5,n6])-PUx([n4,n5,n6])).^2+(Uy([n4,n5,n6])-PUy([n4,n5,n6])).^2); 
                S2=S2+area/3 * sum((stm([n4,n5,n6])-Stream([n4,n5,n6])).^2);
            end
            
    end
    erU=sqrt(S1);
    erS=sqrt(S2);
    
end

function [z]=baseP2(l1,l2,l3)
% Fonctions de base P2
z = [l1.*(2*l1-1); l2.*(2*l2-1);l3.*(2*l3-1);
     4*l2.*l3; 4*l1.*l3;4*l1.*l2];

end
