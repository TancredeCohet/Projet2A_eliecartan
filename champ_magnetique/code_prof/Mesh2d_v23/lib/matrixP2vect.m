function [AA,MM]=matrixP2vect(v1,t1)
%
% Calcul de la matrice de rigidité AA en P2
% et de la matrice de masse MM en P2
%
% Version "vectorisée" (bcp plus rapide: facteur > 10...) de matrixP2.m : 
% on travaille avec des vecteurs à  la place des matrices (construction des 
% matrices tout à  la fin avec A=sparse(i,j,val)
%
    nv1=size(v1,1);
    nt1=size(t1,1);

    %Matrice de masse elementaire P2/P2
    matMM_elem=[ 6 -1 -1 -4  0  0;
                -1  6 -1  0 -4  0;
                -1 -1  6  0  0 -4;
                -4  0  0 32 16 16;
                 0 -4  0 16 32 16;
                 0  0 -4 16 16 32]/180;
     
    MM_elem = reshape(matMM_elem',1,36);         % stockage ligne par ligne

    %Matrices de rigidité élémentaires P2/P2
    matDxDx_elem=[ 3  1  0  0  0 -4;
                   1  3  0  0  0 -4;
                   0  0  0  0  0  0;
                   0  0  0  8 -8  0;
                   0  0  0 -8  8  0;
                  -4 -4  0  0  0  8]/6;
    
    DxDx_elem = reshape(matDxDx_elem',1,36);    % stockage ligne par ligne

    matDyDy_elem=[ 3  0  1  0 -4  0;
                   0  0  0  0  0  0;
                   1  0  3  0 -4  0;
                   0  0  0  8  0 -8;
                  -4  0 -4  0  8  0;
                   0  0  0 -8  0  8]/6;

    DyDy_elem = reshape(matDyDy_elem',1,36);      % stockage ligne par ligne

    matDxDy_elem=[ 3  0  1  0 -4  0;
                   1  0 -1  4  0 -4;
                   0  0  0  0  0  0;
                   0  0  4  4 -4 -4;
                   0  0 -4 -4  4  4;
                  -4  0  0 -4  4  4]/6;

    DxDy_elem = reshape(matDxDy_elem',1,36);              % stockage ligne par ligne
    DxDy_elem_transp = reshape(matDxDy_elem,1,36);  % stockage colonne par colonne

    cmpt=1;
    ntt=floor(nt1/20)*20;
    fprintf('(1...20)=');

    indMMi=zeros(nt1*36,1);
    indMMj=zeros(nt1*36,1);
    MMval=zeros(nt1*36,1);
    AAval=zeros(nt1*36,1);

    ind_ligne=1;


    for k=1:nt1
        if(mod((k*20),ntt)==0)
            fprintf(strcat(num2str(cmpt),' '));
            cmpt=cmpt+1;
        end
        i1=t1(k,1:3);
        i2=[i1 t1(k,6:8)+nv1];
        p=v1(i1,1:2);
        aire=t1(k,5);

        %Vecteur de masse P2/P2
        MMval(ind_ligne:ind_ligne+35)=aire*MM_elem';

        indice=ones(6,1)*i2;
        indMMi(ind_ligne:ind_ligne+35)=indice(:);
        indice=indice';
        indMMj(ind_ligne:ind_ligne+35)=indice(:);

        %Vecteur de rigidite
        DxDxk= (p(3,2)-p(1,2))^2*DxDx_elem...
           -(p(3,2)-p(1,2))*(p(2,2)-p(1,2))*(DxDy_elem+DxDy_elem_transp)...
           +(p(2,2)-p(1,2))^2*DyDy_elem;
       
        DyDyk= (p(3,1)-p(1,1))^2*DxDx_elem...
           -(p(3,1)-p(1,1))*(p(2,1)-p(1,1))*(DxDy_elem+DxDy_elem_transp)...
           +(p(2,1)-p(1,1))^2*DyDy_elem;
       
        AAval(ind_ligne:ind_ligne+35)=(DxDxk+DyDyk)/(2*aire);
 
        ind_ligne=ind_ligne+36;

    end

    % Construction des matrices
    MM=sparse(indMMi,indMMj,MMval);
    AA=sparse(indMMi,indMMj,AAval);

    fprintf('\n');
end