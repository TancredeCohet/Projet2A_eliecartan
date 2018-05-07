function outmesh=meshP1(node,edge,hdata,options)
%--------------------------------------------------------------------------
% MAILLAGE P1 
%--------------------------------------------------------------------------
%
%   v1 : coordonnées des noeuds v1=[x, y, label]
%        label : réference des noeuds :
%                label=0 -> noeud du domaine intérieur
%                label=i -> noeud sur le bord défini par edge(i), 
%                pour i=1,...,nbedge.
%
%   e1 : aretes e1=[n1, n2, label, t1, t2] 
%        - (n1,n2) sont les deux noeuds de l'arete
%        - (t1,t2) sont les deux triangles de part et d'autre de l'arete 
%          (t2=0, si l'arete est sur le bord).
%
%   t1 : triangles t1=[n1,n2,n3, label,area, e1,e2,e3]
%        - (n1,n2,n3) trois sommets des triangles
%        - label: rÃ©fÃ©rences du domaine
%        - area: aires des triangles 
%        - (e1,e2,e3) trois aretes  des triangles
%
% * Supports des noeuds
%
%   isupp1,tsupp1 : permet de donner la liste des triangles ayant un noeud 
%                   donné comme sommet (support).
%      isupp1 (i) : adresse dans le vecteur tsupp du premier triangle ayant 
%                   le noeud i comme sommet.
%                   La liste complÃ¨te des triangles ayant le noeud i comme 
%                   sommet, est donnÃ©e par :
%
%                   tsupp(k) pour k = isupp(i) + 1, ... , isupp(i + 1).
%
%--------------------------------------------------------------------------

[p,t] = mesh2d(node,edge,hdata,options);

%--------------------------------------------------------------------------
% Construction attendue du maillage P1
% %--------------------------------------------------------------------------
% Référence des noeuds :
%  0 -> noeud du domaine intérieur
%  i -> noeud sur le bord défini par edge(i) i=1,...,nbedge
v1 = [p, findedge(p,node,edge)];

% Construction d'autres donnÃ©es du maillage:

%    Aretes du maillage P1
t = [t ones(size(t,1),1)];
e1 = make_edges(v1,t);

%    Calcul des aires des triangles
t1 = area_triangle(v1,t);

%    Recuperation des (trois) aretes de chaque triangle.
t1 = edges_in_triangle(e1,t1);

%    Construction du support noeud -> triangle 
[isupp1,tsupp1] = make_supp(t1);

% nv1 : nb de noeuds du maillage P1
nv1 = size(p,1);

%--------------------------------------------------------------------------
% Indices des noeuds du bord
%--------------------------------------------------------------------------
ib1 = v1(:,3)~= 0;       %noeuds du bord
ic1 = (1:nv1)';             
ic1(ib1) = [];           %suppression dans ic1 de tous les points du bord 
                         %ic1 : noeuds interieurs

outmesh = struct('v1',v1,'t1',t1,'e1',e1,'isupp1',isupp1,'tsupp1',tsupp1,...
    'ib1',ib1,'ic1',ic1);

end