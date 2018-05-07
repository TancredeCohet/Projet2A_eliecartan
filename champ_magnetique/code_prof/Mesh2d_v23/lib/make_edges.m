function e = make_edges(v,t)
%--------------------------------------------------------------------------
% Cette fonction permet de définir les aretes des triangles à partir 
% des points v et des triangles t.
%
% e = [n1 n2 label T_gauche T_droit]
%
% n1,n2 : numero des noeuds/extrémités des aretes
% label : références des aretes; 
%         aretes internes : label=0
%         aretes du bord : label = min( ref(n1), ref(n2) )
% 
% T_gauche, T_droit : Triangles à gauche et à droite de chaque arete dans
%                     le sens de parcours n1 -> n2;
%                     T_droit=0 pour une arete sur le bord
% 
%--------------------------------------------------------------------------

    % on prend les points composant chaque segment du triangle
    e=[t(:,[1 2 4]);t(:,[2 3 4]);t(:,[3 1 4])];
    nt=size(t,1);n=(1:nt);n=[n n n]';
    e=[e n n*0];

    % Tri des extrémités des aretes : le noeud de plus petit indice est 
    % dans la premiere colonne de e.
    k=find(e(:,1)>e(:,2));
    tmp=e(k,1);e(k,1)=e(k,2);e(k,2)=tmp;
    
    % Permutation des colonnes 4 et 5 correspondantes
    e(k,5)=e(k,4); e(k,4)=0;

    % Tri décroissant de la deuxième puis de la première colonne de e 
    [~,j]=sort(e(:,2),'descend'); e=e(j,:);
    [~,j]=sort(e(:,1),'descend'); e=e(j,:);

    %----------------------------------------------------------------------
    % Suppression des aretes comptées deux fois.
    % ---------------------------------------------------------------------
    % Recherche des sommets tels que (n1,n2)(j)==(n1,n2)(j+1) :
    % on recherche les aretes qui sont mises 2 fois, i.e. les aretes qui
    % appartiennent à 2 triangles à la fois
    ne=size(e,1);
    i=1:(ne-1);
    j=find((e(i,1)==e(i+1,1)) & (e(i,2)==e(i+1,2)));

    % Eliminations des doublons :
    % On commence par "concaténer" les numéros des triangles.
    jj=find(e(j+1,4)>e(j,4));
    e(j(jj),4)=e(j(jj)+1,4);
    jj=find(e(j+1,5)>e(j,5));
    e(j(jj),5)=e(j(jj)+1,5);
    
    % On élimine les aretes définies 2 fois
    if ~isempty(j)
        e(j+1,:)=[];
    end

    % Si Tgauche=0, on permute les deux noeuds pour avoir Tdroit=0
    J=find(e(:,4)==0);
    tmp=e(J,1);
    e(J,1)=e(J,2);
    e(J,2)=tmp;
    
    e(J,4)=e(J,5);
    e(J,5)=0;
    
    %----------------------------------------------------------------------
    % Label (reference) des aretes
    %----------------------------------------------------------------------
    % Cas des aretes internes -> label=0
    J=(e(:,4)~=0) & (e(:,5)~=0);
    e(J,3)=0;
    
    % Cas des aretes sur le bord: -> label=min(ref(n1),ref(n2))
    J=find(e(:,5)==0);
    n1=e(J,1); n2=e(J,2);
    e(J,3)=min([v(n1,3) v(n2,3)],[],2);
 
end
