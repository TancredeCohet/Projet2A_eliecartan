function [A, M] = matrixP1vect(v,t)
% -------------------------------------------------------------------------
% Calcul (exact) des matrices P1 de masse et de rigidite
%-------------------------------------------------------------------------
% Fonctions de base P1 : v_j, j=1,...,nv
%
% Matrice de masse P1/P1 :
% -----------------------
%     /
% M = | (v_i)(v_j) dxdy; i,j=1,...,nv
%     /
%      O
%
% Matrice de rigidite P1/P1 :
% --------------------------
%     /  
% A = | grad(v_i) . grad(v_j) dxdy
%     / 
%      O
%
%-------------------------------------------------------------------------
% Version "vectorisée" (bcp plus rapide: facteur > 10...) de matrixP1.m : 
% on travaille avec des vecteurs à  la place des matrices (construction des 
% matrices tout à  la fin avec A=sparse(i,j,val)
% -------------------------------------------------------------------------
     nt=size(t,1); %nombre de triangle
    
   % Matrice de masse elementaire P1/P1 (exacte)
   M_elem = [ 2 1 1;
	           1 2 1;
	           1 1 2]/12;
    M_elem=M_elem(:);                % stockage colonne par colonne
    
    % Matrice de rigidite elementaire P1/P1 (exacte)
    DxDx_elem = [ 1 -1  0;
		         -1  1  0;
		          0  0  0]/2;
    DxDx_elem=DxDx_elem(:);    % stockage colonne par colonne
    
    DyDy_elem = [ 1  0 -1;
		          0  0  0;
		         -1  0  1]/2;
    DyDy_elem=DyDy_elem(:);     % stockage colonne par colonne, met tout sous la forme d'une seule colonne
    
    DxDy_elem = [ 2 -1 -1;
		         -1  0  1;
		         -1  1  0]/2;
    DxDy_elem = DxDy_elem(:);  % stockage colonne par colonne, met tout sous la forme d'une seule colonne
  
    fprintf('Assemblage des matrices (1...20)='); 
    cmpt = 1; 
    
    ntt = floor(nt/20)*20;
    
    Mval = zeros(nt*9,1);
    Aval = zeros(nt*9,1);
    
    ind_ligne = 1;    
    
    for k = 1:nt %boucle sur chaque triangle
        if mod((k*20),ntt)==0 
            fprintf('%d',cmpt); cmpt=cmpt+1; 
        end
      
      i = t(k,1:3);         %1*3 avec noeuds ieme triangle
      p = v(i,1:2);         %3*2 coordonnees des noeuds du ième triangle  [xn1 yn1 ; xn2 yn2 ; xn3 yn3]
      aire = t(k,5);        %aire du triangle i

     
      indice = i'*ones(1,3);                    %matrice 3*3 chaque colonne [n1;n2;n3] noeuds du ième triangle
      ligne(ind_ligne:ind_ligne+8) = indice(:); %[n1 n2 n3 n1 n2 n3 n1 n2 n3]
      indice = indice';
      colonne(ind_ligne:ind_ligne+8) = indice(:); %[n1 n1 n1 n2 n2 n2 n3 n3 n3]
        
      % Matrice de masse P1/P1
      Mval(ind_ligne:ind_ligne+8) = aire*M_elem;

      % Matrice de rigidite
      Ak = ((p(3,2)-p(1,2))^2 + (p(3,1)-p(1,1))^2) * DxDx_elem ...
	      +((p(2,2)-p(1,2))^2 + (p(2,1)-p(1,1))^2) * DyDy_elem ...   
          -((p(3,2)-p(1,2))*(p(2,2)-p(1,2)) + ...
            (p(3,1)-p(1,1))*(p(2,1)-p(1,1)) ) * DxDy_elem;
	  
      Aval(ind_ligne:ind_ligne+8) = Ak/(2*aire);   
       
      ind_ligne = ind_ligne+9;     
      
    end

    A=sparse(ligne,colonne,Aval); %on place les valeurs avec les nums des noeuds
    M=sparse(ligne,colonne,Mval);
fprintf('\n');
end
  
