H = 1.05;         % hauteur du domaine
L = 1.5;            % longueur du domaine
pas_h = 0.1;   % pas du maillage "interieur"

%--------------------------------------------------------------------------
% Maillage
%----------------------------------------------------------------------------
% Bord exterieur du domaine (rectangle) 
node = [0 0; L 0 ; L H ; 0 H ; L/4 H/4 ; L/4 3*H/4 ; 3*L/4 3*H/4 ; 3*L/4 H/4];
edge = [1 2 ; 2 3 ; 3 4 ; 4 1; 5 6 ; 6 7 ;7 6 ;6 5];

fprintf('-----------------------------------------------------------------\n');
fprintf('Maillage par ''mesh2d'' v2.3 (Copyright (C) 2007 Darren Engwirda)\n');
fprintf('-----------------------------------------------------------------\n');
tstart=tic; %lance chrono temps de calcul

%implementation longueur des arete
hdata = [];
hdata.hmax  = pas_h; 

%definition des options du maillage
optionsmesh = [];
optionsmesh.output = false;



[v1,t1] = mesh2d(node,edge,hdata,optionsmesh);


fprintf('Maillage P1 : nombre de noeuds    (nv1) = %d\n',nv1);
fprintf('              nombre de triangles (nt1) = %d\n',nt1);
fprintf('Temps de construction du maillage: %f s\n',toc(tstart));

%affichage du maillage
patch('faces',t1(:,1:3),'vertices',v1(:,1:2),'facecolor','white','edgecolor','blue');
axis equal; axis on
drawnow 
