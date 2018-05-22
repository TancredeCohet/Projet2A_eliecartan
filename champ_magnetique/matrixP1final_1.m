function [M, nn, ibint, ic2] = matrixP1final(vz,tz,fnum,node1,edge1,node2,edge2)
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
    nt = size(tz,1);                     %nombre de triangle
    nn = size(vz,1);                     %nombre de noeuds
    M = zeros(nn,1);                    %initialisation matrice 
    %M = ones(nn,1);   
    v1 = [vz findedge(vz,node1,edge1,0.001)];  %noeuds du bord de la face 1 = aimant
    
    %les noeuds sur les sommet du maillage ont un probleme de label, on
    %reassigne les labels ici meme si n'est pas tres elegant...
      v1(find(v1(:,1)==node1(1,1)& v1(:,2)==node1(1,2)),3) = 1;
      v1(find(v1(:,1)==node1(2,1)& v1(:,2)==node1(2,2)),3) = 1;
      v1(find(v1(:,1)==node1(3,1)& v1(:,2)==node1(3,2)),3) = 3;
      v1(find(v1(:,1)==node1(4,1)& v1(:,2)==node1(4,2)),3) = 3;
    v2 = [vz findedge(vz,node2,edge2)];  %noeuds du bord de la face 3 = maillage total
    
    
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

% triangle_a_cheval = 0;
    for k = 1:nt
        Tk = tz(k,1:3);                                  %on recupere les noeuds definisant le triangle k  [n1 n2 n3]
%--------------------------------------------------------------------------   
%                        calcul du terme source M
%--------------------------------------------------------------------------
        
        if fnum(k) == 1                                 %on teste si le triangle est sur le maillage de l'aimant
            
%             if k == 28938
%                  on_est_sur_le_triabgle_du_bord = 1;
%             end                                        
            l = v1(Tk,3);                           %l regroupe les valeurs des labels
            vN1 = v1(Tk,1:3);                       %coordonees des noeuds du triangle
        
            if (size(l(ismember(l,0))) == [1 1]) | (size(l(ismember(l,0))) == [0 1])
                
                compter = sum(l == l');
                z = vN1(:,3) ~= l(compter == 1);
%                 if z == [0;0;0]                   %on repere les triangle
%                 a cheval
%                     triangle_a_cheval =1 +  triangle_a_cheval ;
%                 end
                    
                if size(z(ismember(z,0))) == [1 1]          %on elimine le cas d'un triangle bizarre a cheval sur 2 bords
                    Nb = vN1(z,1:2);                        %on recupere les coordonnees des noeuds du triangles qui sont sur le bord 
                    NNb = Tk(z);                            %on recupere les numeros des triangles sur le bord
                    Ni = vN1(vN1(:,3)==l(compter==1),1:2);              %noeuds interieurs
                    norme_e = sqrt((Nb(1,2) - Nb(2,2))^2 + (Nb(1,1) - Nb(2,1))^2);  %on calcule la norme de l'arete
                    n_normal = [Nb(2,2)-Nb(1,2); Nb(1,1)-Nb(2,1)];  %vecteur normal a l'arete
                    n_test = [Ni(1,1)-Nb(1,1); Ni(1,2)-Nb(1,2)];
                    if n_normal'*n_test < 0                         %on teste si le vecteur normal que l'on a calcule est dans le bon sens
                        n_normal = -n_normal;                       %dans ce cas on le prend dans l'autre sens pour qu'il soit sortant
                    end
                    M_inter =  (norme_e / 2).*[0 ; 1]'*n_normal;      % [0;1] correspond au vecteur ey
                    M(NNb,:) =  M(NNb,:) + M_inter;                              %on place la valeur corespondant a l'arette dans la matrice du terme source 
                end                                                    %aux indices des noeuds sur le bord
            end           
        end
    end
   
%
%--------------------------------------------------------------------------
%                             matrice de rigidite
%--------------------------------------------------------------------------                     

        pk = vz(Tk,1:2);     %3*2 coordonnees des noeuds triangle k [xn1 yn1 ; xn2 yn2 ; xn3 yn3]
        %calcul matrice des indice pour placer les valeurs dans la matrice
        %sparse
        indice=Tk'*ones(1,3);                           %3*3 chaque colonne [n1 n1 n1 ;n2 n2 n2;n3 n3 n3] noeuds du triangle i
        ligne(ind_ligne:ind_ligne+8) = indice(:);         %on stocke sous la forme de colonne
        indice=indice';                                 
        colonne(ind_ligne:ind_ligne+8) = indice(:);       %[n1;n1;n1;n2;n2;n2;n3;n3;n3]
       
        %implementation de A
        %A_inter matrice elementaire 3*3 des gradients dans le triangle k
%         terme_jacobien = (airek^2)/4 .* ([pk(3,2)-pk(2,1) pk(1,1)-pk(3,1) ; pk(2,1)-pk(2,2) pk(2,1)-pk(1,1)]^2);
%         A_inter = zeros(3,3);
%         A_inter(1,1) = (terme_jacobien*grad_elem(:,1))'*grad_elem(:,1);
%         A_inter(2,2) = (terme_jacobien*grad_elem(:,2))'*grad_elem(:,2);
%         A_inter(3,3) = (terme_jacobien*grad_elem(:,3))'*grad_elem(:,3);
%         A_inter(1,2) = (terme_jacobien*grad_elem(:,1))'*grad_elem(:,2);
%         A_inter(1,3) = (terme_jacobien*grad_elem(:,1))'*grad_elem(:,3);
%         A_inter(2,3) = (terme_jacobien*grad_elem(:,2))'*grad_elem(:,3);
%         A_inter(2,1) = A_inter(1,2); 
%         A_inter(3,1) = A_inter(1,3); 
%         A_inter(3,2) = A_inter(2,3); 
%        
%         A_inter = A_inter(:);
%         A_val(ind_ligne:ind_ligne + 8) = A_inter;
%           
%         ind_ligne = ind_ligne + 9;      %iteration
%}
%    end
        
%    A = sparse(ligne,colonne,A_val);        %on place les valeurs calcules avec les valeurs des noeuds
   
%--------------------------------------------------------------------------
%                      On ne retient que les noeuds interieurs
%--------------------------------------------------------------------------
    ic2 = (1:nn)';  
    ibint = v2(:,3)~= 0;                    %noeuds du bord de la face 2_exterieure
    ic2(ibint) = [];                        %liste des noeuds interieur au maillage total
%   A = A(ic2,ic2);                         %on ne garde que les noeuds interieurs

%M = ones(nn,1); 

%affichage du second membre
%close all;
 figure();
 patch('faces',tz(:,1:3),'vertices',vz(:,1:2),'FaceVertexCData',M,'facecolor','interp','edgecolor','black');
 axis equal;
 title('second membre');
 xlim([4.5 6]);
 ylim([2 6]);
 colorbar;
 M(ibint) = [];                          %suppression des elements sur le bord
end