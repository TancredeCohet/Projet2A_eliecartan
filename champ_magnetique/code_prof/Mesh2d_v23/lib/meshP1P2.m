function outmesh=meshP1P2(node,edge,hdata,options);
%--------------------------------------------------------------------------
% MAILLAGES P1 et P2 (cf. doc_maillageP1P2.pdf)
%--------------------------------------------------------------------------
%
% * Maillage P1 :
%
%   v1 : coordonné́es des noeuds v1=[x, y, label]
%        label : référence des noeuds :
%                label=0 -> noeud du domaine intérieur
%                label=i -> noeud sur le bord défini par edge(i), 
%                pour i=1,...,nbedge.
%
%   e1 : aretes e1=[n1, n2, label, t1, t2] 
%        - (n1,n2) sont les deux noeuds de l'arete
%        - (t1,t2) sont les deux triangles de part et d'autre de l'arete 
%          (t2=0, si l'arete est sur le bord).
%
%   t1 : triangles t1=[n1,n2,n3, label,aire, e1,e2,e3]
%        - (n1,n2,n3) trois sommets des triangles
%        - (e1,e2,e3) trois aretes  des triangles
%
%
% * Maillage P2 :
%
%   v2,e2 : coordonnées des noeuds et aretes (meme structure que v1,e1)
%
%   t2    : triangles obtenus en coupant en 4 les triangles du maillage P1.
%           t2 = [n1, m2, m3, label,aire, e1,e2,e3;
%                 n2, m1, m3, label,aire, e1,e2,e3;
%                 n3, m2, m1, label,aire, e1,e2,e3;
%                 m1, m2, m3, label,aire, e1,e2,e3]
%
%           où m1,m2,m3 sont les points milieux des aretes des triangles
%           P1; m1,m2,m3 sont sur les aretes opposées aux noeud n1,n2,n3
%           respectivement.
%
% * Supports (P1 et P2)
%
% isupp,tsupp : permet de donner la liste des triangles ayant un noeud 
%                 donné comme sommet (support).
%      isupp(i) : adresse dans le vecteur tsupp du premier triangle ayant 
%                 le noeud i comme sommet.
%                 La liste complète des triangles ayant le noeud i comme 
%                 sommet, est donnée par :
%
%                 tsupp(k) pour k = isupp(i) + 1, ... , isupp(i + 1).
%
%--------------------------------------------------------------------------

[p,t] = mesh2d(node,edge,hdata,options);

%--------------------------------------------------------------------------
% Construction étendue du maillage P1
%--------------------------------------------------------------------------
% Référence des noeuds :
%  0 -> noeud du domaine intérieur
%  i -> noeud sur le bord défini par edge(i) i=1,...,nbedge
v1=[p, findedge(p,node,edge)];

% Construction d'autres données du maillage:

%    Aretes du maillage P1
t=[t ones(size(t,1),1)];
e1=make_edges(v1,t);

%    Calcul des aires des triangles
t1=area_triangle(v1,t);

%    Recuperation des (trois) aretes de chaque triangle.
t1=edges_in_triangle(e1,t1);

%    Construction du support noeud -> triangle 
[isupp1,tsupp1]=make_supp(t1);

%--------------------------------------------------------------------------
% Construction du maillage P2 (à partir du maillage P1)
%--------------------------------------------------------------------------
[v2,t2]=P1toP2(v1,t1,e1);
e2=make_edges(v2,t2);
t2=area_triangle(v2,t2);
t2=edges_in_triangle(e2,t2);
[isupp2,tsupp2]=make_supp(t2);

% nv1, nv2 : nb de noeuds des maillage P1 et P2
nv1=size(p,1);  nv2=size(v2,1); 

%--------------------------------------------------------------------------
% Indices des noeuds du bord
%--------------------------------------------------------------------------
ib2=find(v2(:,3)~=0);  %noeuds du bord
ic2=(1:nv2)';             
ic2(ib2)=[];           %suppression dans ic2 de tous les points du bord 
                              %ic2 : noeuds interieurs

outmesh=struct('v1',v1,'t1',t1,'e1',e1,'isupp1',isupp1,'tsupp1',tsupp1,...
    'v2',v2,'t2',t2,'e2',e2,'isupp2',isupp2,'tsupp2',tsupp2,...
    'ib2',ib2,'ic2',ic2);

end