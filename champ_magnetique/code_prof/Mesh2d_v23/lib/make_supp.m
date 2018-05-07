function [isupp,tsupp]=make_supp(t)
%
% Les tableaux isupp et tsupp1 permettent de donner la liste des 
% triangles ayant un noeud donné comme sommet.
% isupp(i) : adresse dans le vecteur tsupp du premier triangle ayant 
% le noeud i comme sommet. 
% La liste complète des triangles ayant le noeud i comme sommet, 
% est donnée par :
%
%    tsupp(k) pour k = isupp(i),..., isupp(i+1) − 1.
%
    nt=size(t,1);
    I=(1:nt)';
    ss=[t(:,1) I;t(:,2) I;t(:,3) I];
    [i,j]=sort(ss(:,1),'ascend');
    ss=ss(j,:);
    i=1:(3*nt-1);
    isupp=find(ss(i,1)~=ss(i+1,1));
    isupp=[0 isupp' 3*nt];
    tsupp=ss(:,2);

end
