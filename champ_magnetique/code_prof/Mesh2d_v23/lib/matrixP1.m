function [A,M] = matrixP1(v,t)
%--------------------------------------------------------------------------
% Calcul (exact) des matrices P1 de masse et de rigidite
%
% Fonctions de base P1 : v_j, j=1,...,nv1
%
% Matrice de masse P1/P1 :
% -------------------------------------------------------------------------
%     /
% M = | (v_i)(v_j) dxdy; i,j=1,...,nv1
%     /
%      O
%
% Matrice de rigidite P1/P1 :
% -------------------------------------------------------------------------
%     /  
% A = | grad(v_i) . grad(v_j) dxdy
%     / 
%      O
%--------------------------------------------------------------------------

    nv=size(v,1); nt=size(t,1);
    M=sparse(nv,nv); 
    A=M;
    
    % Matrice de masse elementaire P1/P1 (exacte)
    M_elem = [ 2 1 1;
	           1 2 1;
	           1 1 2]/12;

    % Matrice de rigidite elementaire P1/P1 (exacte)
    DxDx_elem = [ 1 -1  0;
		         -1  1  0;
		          0  0  0]/2;
    
    DyDy_elem = [ 1  0 -1;
		          0  0  0;
		         -1  0  1]/2;
    
    DxDy_elem = [ 2 -1 -1;
		         -1  0  1;
		         -1  1  0]/2;
     
    cmpt=1; fprintf('Assemblage des matrices (1...20)='); ntt=floor(nt/20)*20;

    for k=1:nt

      if mod((k*20),ntt)==0 
          fprintf('%d',cmpt); cmpt=cmpt+1; 
      end
      i = t(k,1:3);
      p = v(i,1:2);
      aire = t(k,5);

      % Matrice de masse P1/P1
      M(i,i)=M(i,i)+aire*M_elem;

      % Matrice de rigidite DxDx
      Ak = ((p(3,2)-p(1,2))^2 + (p(3,1)-p(1,1))^2) * DxDx_elem ...
	      +((p(2,2)-p(1,2))^2 + (p(2,1)-p(1,1))^2) * DyDy_elem ...   
          -((p(3,2)-p(1,2))*(p(2,2)-p(1,2)) + ...
            (p(3,1)-p(1,1))*(p(2,1)-p(1,1)) ) * DxDy_elem;
	      
      A(i,i) = A(i,i) + Ak/(2*aire);
            
    end

    fprintf('\n');
end
