%=========================================================================
%                       Calcul du gradient du potentiel U
%=========================================================================
function [B] = gradient(u,v,t);
%           u matrice colonne des potentiels u
%           v matrice colonne des coordon�es des noeuds
%           t matrice des triangles
%-------------------------------------------------------------------------
%           On utilise la formulation variationelle du probl�me ce qui
%           revient � resoudre le syst�me lin�aire M * B = C * U
%
%           Avec les fonctions de base P1 : v_j, j=1,...,nv
%
%                   M est la matrice de masse  
%                           /
%                   M_i_j = | v_i*v_j dxdy; i = 1,...,nv
%                           /
%                           O
%                            /  
%                    C_i_j = | grad(v_i) . v_j dxdy
%                            / 
%                            O
%           Pour claculer la matrice C on r�alise une boucle for en
%           calculant sur chaque triangle les contributions � la matrice C
%--------------------------------------------------------------------------

 nt = size(t,1);                     %nombre de triangle
 nn = size(v,1);                     %nombre de noeuds

 
 grad_elem_ref = [-1 1 0;            %gradient de phi1, phi2, phi3 sur le triangle de r�f�rence
                  -1 0 1];
 
 ind_ligne = 1;              %initialisation premiere ligne
 
 M_elem = [ 2 1 1;           %matrice de masse elementaire
	        1 2 1;
	        1 1 2]/12;
 M_elem = M_elem(:);
 
 M_val = zeros(nt*9,1);

 %-------------------------------------------------------------------------
 %                      boucle sur les triangles
 %-------------------------------------------------------------------------
 for k = 1:nt
     
    Tk = t(k,1:3);                                  %on recupere les noeuds definisant le triangle k  [n1 n2 n3]
    vk = v(Tk,1:2);                                 %coordonn�es des noeuds du triangle k   
    airek = t(k,5);                                 %aire du triangle k

    %----------------------------------------------------------------------
    %on recupere les indices des neuds pour ensuite placer les 9
    %contributions du triangle k dans la matrice C
    indice=Tk'*ones(1,3);                           %3*3 chaque colonne [n1 n1 n1 ;n2 n2 n2;n3 n3 n3] noeuds du triangle i
                                                    % [n1;n2;n3;n1;n2;n3;n1;n2;n3]
    ligne(ind_ligne:ind_ligne+8)=indice(:);         %on stocke sous la forme de colonne
    indice=indice';                                 
    colonne(ind_ligne:ind_ligne+8)=indice(:);       %[n1;n1;n1;n2;n2;n2;n3;n3;n3]
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    %                       Calcul matrice C
    %----------------------------------------------------------------------
    
    %Jacobien
    %Jk = [vk(3,2)-vk(1,2) vk(1,1)-vk(3,1); vk(1,2)-vk(2,2) vk(2,1)-v(1,1)];  !! ERREUR : vk au lieu de v !!!
    Jk = [vk(3,2)-vk(1,2) vk(1,1)-vk(3,1); vk(1,2)-vk(2,2) vk(2,1)-vk(1,1)];
    %Jk = Jk/(2*airek);
    
    %grad_elem = Jk * grad_elem_ref;  ERREUR : c'est la transposée de Jk
    grad_elem = Jk' * grad_elem_ref;
    
    C_elem_x1 = (1/6) * [ grad_elem(1,1)*ones(1,3);
                          grad_elem(1,2)*ones(1,3);
                          grad_elem(1,3)*ones(1,3)];
    
    C_elem_x1 = C_elem_x1(:);
    C_val_x1(ind_ligne:ind_ligne + 8) = C_elem_x1;
          
    C_elem_x2 = (1/6) * [ grad_elem(2,1)*ones(1,3);
                          grad_elem(2,2)*ones(1,3);
                          grad_elem(2,3)*ones(1,3)];

    C_elem_x2 = C_elem_x2(:);
    C_val_x2(ind_ligne:ind_ligne + 8) = C_elem_x2;
    
    %----------------------------------------------------------------------
    %                       Calcul matrice M
    %----------------------------------------------------------------------
    
    M_val(ind_ligne:ind_ligne+8) = airek*M_elem;
    
    ind_ligne = ind_ligne + 9;      %it�ration
 end
    
 
 C_x1 = sparse(ligne,colonne,C_val_x1);                       %on place les valeurs avec les nums des noeuds
 C_x2 = sparse(ligne,colonne,C_val_x2);
 
 M = sparse(ligne,colonne,M_val);
 
 %C_x1 = full(C_x1);
 %C_x2 = full(C_x2);

 
 %-------------------------------------------------------------------------
 %                          Resolution du systeme lineaire 
 %-------------------------------------------------------------------------
 
 B_x1 = M\(C_x1'*u);
 B_x2 = M\(C_x2'*u);
 B = [B_x1 B_x2];
end
 