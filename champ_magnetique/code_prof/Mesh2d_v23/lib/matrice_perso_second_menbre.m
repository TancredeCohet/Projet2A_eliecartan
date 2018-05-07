function M = matrixP1vectperso(v,t)
% -------------------------------------------------------------------------
% Calcul (exact) des matrices P1 de masse et de rigidite
%-------------------------------------------------------------------------
% Fonctions de base P1 : v_j, j=1,...,nv
%
% Matrice de masse P1/P1 :
% -----------------------
%     /
% M = | grad(v_i) dxdy; i = 1,...,nv
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
    nt = size(t,1); %nombre de triangle
    
    %Matrice de masse elementaire P1/P1 (exacte)
    M_elem = [ -1 1 0;
	           -1 0 1];
    %M_elem = M_elem(:);                % stockage colonne par colonne
  
    fprintf('Assemblage des matrices (1...20)='); 
 
    Mvaleurs=zeros(nt*6,1);
    %Aval = zeros(nt*9,1);
    ind_ligne = 1;
    for k = 1:nt            %boucle sur chaque triangle, on calcule les valeurs pour chaque noeuds
        
      
      i = t(k,1:3);         %1*3 avec noeuds triangle i  [n1 n2 n3]
      p = v(i,1:2);         %3*2 coordonnees des noeuds triangle i [xn1 yn1 ; xn2 yn2 ; xn3 yn3]
      airek = t(k,5);       %aire du triangle k

      indice = i'*ones(1,2);                    %matrice 3*3 chaque colonne [n1 n1 ;n2 n2 ;n3 n3] noeuds du triangle i
      ligne(ind_ligne:ind_ligne + 5) = indice(:); %[n1 n2 n3 n1 n2 n3]
      indice = indice';
      colonne(ind_ligne:ind_ligne + 5) = indice(:); %[n1 n1 n2 n2 n3 n3]
        
      % Matrice de masse P1/P1
      M_inter = airek * [p(3,2)-p(2,1) p(1,1)-p(3,1) ; p(2,1)-p(2,2) p(2,1)-p(1,1)] * M_elem;
      Mvaleurs(ind_ligne:ind_ligne + 5) = M_inter(:); %mise sous la forme d'une colonne
  
      ind_ligne = ind_ligne + 6; %itération 
      
    end
    M=Mvaleurs;

    fprintf('\n');
    size_M=size(M);
    size_M=size_M(1);
    mat_0011=sparse(floor(size_M/2),size_M);
    cmptr=0;
    for i=1:1:floor(size_M/2);
        cmptr=cmptr+2;
        mat_0011(i,cmptr-1:cmptr)=[1 1];
    end
    M=mat_0011*sparse(M);
    
end
