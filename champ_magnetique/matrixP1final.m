function [M, nn, ibint, ic2] = matrixP1final(v,t,fnum,node1,edge1,node2,edge2)
%==========================================================================
%         Calcul (exact) des matrices P1 de masse et de rigidite
%==========================================================================
%                 Fonctions de base P1 : v_j, j=1,...,nv
%
%                          Matrice de masse P1/P1 
%                       /
%                   M = | grad(v_i) dxdy; i = 1,...,nv
%                       /
%                       O
%
%                       Matrice de rigidite P1/P1 :
%                        /  
%                    A = | grad(v_i) . grad(v_j) dxdy
%                        / 
%                        O
%
%==========================================================================

%--------------------------------------------------------------------------
%                             initialision
%--------------------------------------------------------------------------
    nt = size(t,1);                     %nombre de triangle
    nn = size(v,1);                     %nombre de noeuds
    M = zeros(nn,1);                    %initialisation matrice 
    %M = ones(nn,1);   
    v1 = [v findedge(v,node1,edge1)];  %noeuds du bord de la face 1-interieur
    v2 = [v findedge(v,node2,edge2)];  %noeuds du bord de la face 2-exterieur
    
    %matrices elementaire des gradients dans le triangle de reference
    %grad_elem = [1 0 -1;
     %           -1 1 0];
    grad_elem = [-1 1 0;
                 -1 0 1];
    A_val = zeros(nt*9,1);
    
    ind_ligne = 1;              %initialisation premiere ligne

%--------------------------------------------------------------------------   
%                        boucle
%--------------------------------------------------------------------------    
% le principe de l'algorithme est de calculer les valeurs des gradients sur
% chaque triangles et de les placer dans les matrices ensuite. 

    for k = 1:nt
        Tk = t(k,1:3);                                  %on recupere les noeuds definisant le triangle k  [n1 n2 n3]
%--------------------------------------------------------------------------   
%                        calcul du terme source M
%--------------------------------------------------------------------------
        
        if fnum(k) == 1                                 %on teste si le triangle est sur la surface de maillage int�rieure (le terme source du champ B
                                                          %est nul partout en dehors de l'aimant
          l = v1(Tk,3);
          vN1 = v1(Tk,1:3);                       %coordon�es des noeuds du triangle
%          if sommet(v1(Tk,1:2),l)                 %test pour sommets
%              a = 'Coucou'
         
          if size(l(ismember(l,0))) == [1 1]    %on teste si le triangle est le bord de la face 1 
            
           
            z = vN1(:,3)~=0;
            Nb = vN1(z,1:2);                        %on recupere les coordonn�es des noeuds du triangles qui sont sur le bord 
            NNb = Tk(z);                            %on recupere les numeros des triangles sur le bord
            
          
            Ni = vN1(vN1(:,3)==0,1:2);              %noeuds interieurs
            norme_e = sqrt((Nb(1,2) - Nb(2,2))^2 + (Nb(1,1) - Nb(2,1))^2);  %on calcule la norme de l'arete
            n_normal = [Nb(2,2)-Nb(1,2); Nb(1,1)-Nb(2,1)];  %vecteur normal a l'arete
            n_test = [Ni(1,1)-Nb(1,1); Ni(1,1)-Nb(1,2)];
            if n_normal'*n_test < 0                         %on teste si le vecteur normal que l'on a calcul� est dans le bon sens
                n_normal = -n_normal;                       %dans ce cas on le prend dans l'autre sens pour qu'il soit sortant
            end
            M_inter =  (norme_e / 2).*[0 ; 1]'*n_normal;      % [0;1] correspond au vecteur ey
            M(NNb,:) =  M(NNb,:) + M_inter;                              %on place la valeur corespondant a l'arette dans la matrice du terme source 
          end
        end %aux indices des noeuds sur le bord
    end
%
%--------------------------------------------------------------------------
%                             matrice de rigidit� 
%--------------------------------------------------------------------------                     

        pk = v(Tk,1:2);     %3*2 coordonnees des noeuds triangle k [xn1 yn1 ; xn2 yn2 ; xn3 yn3]
        airek=(v(t(k,2),1)-v(t(k,1),1))*(v(t(k,3),2)-v(t(k,1),2))-(v(t(k,3),1)-v(t(k,1),1))*(v(t(k,2),2)-v(t(k,1),2));     %aire du triangle k

        %calcul matrice des indice pour placer les valeurs dans la matrice
        %sparse
        indice=Tk'*ones(1,3);                           %3*3 chaque colonne [n1 n1 n1 ;n2 n2 n2;n3 n3 n3] noeuds du triangle i
        ligne(ind_ligne:ind_ligne+8)=indice(:);         %on stocke sous la forme de colonne
        indice=indice';                                 
        colonne(ind_ligne:ind_ligne+8)=indice(:);       %[n1;n1;n1;n2;n2;n2;n3;n3;n3]
       
        %implementation de A
        %A_inter matrice elementaire 3*3 des gradients dans le triangle k
        terme_jacobien = (airek^2)/4 .* ([pk(3,2)-pk(2,1) pk(1,1)-pk(3,1) ; pk(2,1)-pk(2,2) pk(2,1)-pk(1,1)]^2);
        A_inter = zeros(3,3);
        A_inter(1,1) = (terme_jacobien*grad_elem(:,1))'*grad_elem(:,1);
        A_inter(2,2) = (terme_jacobien*grad_elem(:,2))'*grad_elem(:,2);
        A_inter(3,3) = (terme_jacobien*grad_elem(:,3))'*grad_elem(:,3);
        A_inter(1,2) = (terme_jacobien*grad_elem(:,1))'*grad_elem(:,2);
        A_inter(1,3) = (terme_jacobien*grad_elem(:,1))'*grad_elem(:,3);
        A_inter(2,3) = (terme_jacobien*grad_elem(:,2))'*grad_elem(:,3);
        A_inter(2,1) = A_inter(1,2); 
        A_inter(3,1) = A_inter(1,3); 
        A_inter(3,2) = A_inter(2,3); 
       
        A_inter = A_inter(:);
        A_val(ind_ligne:ind_ligne + 8) = A_inter;
          
        ind_ligne = ind_ligne + 9;      %it�ration
%}
%    end
        
%    A = sparse(ligne,colonne,A_val);        %on place les valeurs calcules avec les valeurs des noeuds
   
%--------------------------------------------------------------------------
%                      On ne retient que les noeuds interieurs
%--------------------------------------------------------------------------
    ic2 = (1:nn)';  
    ibint = v2(:,3)~= 0;                    %noeuds du bord de la face 2_exterieure
    ic2(ibint) = [];
%   A = A(ic2,ic2);                         %on ne garde que les noeuds int�rieurs

%M = ones(nn,1); 
    M(ibint) = [];                          %suppression des elements sur le bord
     
end