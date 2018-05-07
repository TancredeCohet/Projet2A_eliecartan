%==========================================================================

close all;

% Path pour mesh2d (maillage)
addpath 'Mesh2d_v23/' 
addpath 'lib'

tfull=tic;
%--------------------------------------------------------------------------
% Définition de la géométrie et parametres du maillage
%--------------------------------------------------------------------------
% node    : coordonnées [x;y] des noeuds du bord du domaine
% edge    : aretes du bord définies par les numeros des noeuds [n1, n2]
% pas_h   : pas de discretisation du maillage
% hdata   : taille(s) du maillage
% options : options pour le maillage (dhmax,verbose...)

%------------------------------------------------------------------
% Domaine rectangle [0,L]x[0,H]
%------------------------------------------------------------------
L=1; H=1.01;               

node = [0 0; L 0; L H; 0 H];
edge = [1 2; 2 3; 3 4; 4 1];

%node = [0 0; L 0; L H; 0 H; 0 H/2.5];
%edge = [1 2; 2 3; 3 4; 4 5; 5 1];


% Taille du maillage
pas_h=[0.08, 0.04, 0.02];
erreurmax=[];
erreurL2=[];
sizemesh=[];

hdata = [];
optionsmesh=[];
optionsmesh.output=false;  % pas de sortie maillage

for k=1:length(pas_h)

    %--------------------------------------------------------------------------
    % Maillage P1
    %--------------------------------------------------------------------------
    %   v : coordonné́es des noeuds v=[x, y, label]
    %       label : référence des noeuds :
    %               label=0 -> noeud du domaine intérieur
    %               label=i -> noeud sur le bord défini par edge(i), 
    %               pour i=1,...,nbedge.
    %
    %   t : triangles t=[n1,n2,n3, label,aire]
    %        - (n1,n2,n3) trois sommets des triangles
    %        - label = 1
    %        - aire : aire des triangles
    %--------------------------------------------------------------------------
    fprintf('-----------------------------------------------------------------\n');
    fprintf('Maillage par ''mesh2d'' v2.3 (Copyright (C) 2007 Darren Engwirda)\n');
    fprintf('-----------------------------------------------------------------\n');

    hdata.hmax  = pas_h(k);   

    outmesh=meshP1P2(node,edge,hdata,optionsmesh);
    
    v1=outmesh.v1; t1=outmesh.t1; e1=outmesh.e1;
    v2=outmesh.v2; t2=outmesh.t2; e2=outmesh.e2;

    ibd2=outmesh.ib2; iin2=outmesh.ic2;

    % nv1, nt1, ne1 : nb de noeuds, de triangles et d'aretes du maillage P1
    % nv2, nt2, ne2 : nb de noeuds, de triangles et d'aretes du maillage P2
    nv1=size(v1,1); nt1=size(t1,1); ne1=size(e1,1);
    nv2=size(v2,1); nt2=size(t2,1); ne2=size(e2,1);

    fprintf('Maillage P1 : nombre de noeuds    (nv1) = %d\n',nv1);
    fprintf('              nombre de triangles (nt1) = %d\n',nt1);
    fprintf('Maillage P2 : nombre de noeuds    (nv2) = %d\n',nv2);
    fprintf('              nombre de triangles (nt2) = %d\n',nt2);

    % Taille du maillage
    sm=sqrt(max(abs(t1(:,5))));
    sizemesh=[sizemesh,sm];

    fprintf('Maillage P1 : nombre de noeuds    (nv) = %d\n',nv1);
    fprintf('              nombre de triangles (nt) = %d\n',nt1);
    fprintf('taille du maillage hmax = %e\n',sizemesh(k));
    
%     % Affichage du maillage
%     patch('faces',t1(:,1:3),'vertices',v1(:,1:2),'facecolor','white','edgecolor','blue');
%     axis equal; 

    %--------------------------------------------------------------------------
    % Construction des matrices de masse M et de rigidité A   
    %--------------------------------------------------------------------------
    [A,M]=matrixP2vect(v1,t1);

    fprintf('Taille de la matrice de rigidité : %d x %d (nnz = %d)\n',size(A), nnz(A));

    % Terme source
    % ------------
    %f = @(x,y) 10*ones(size(x)); 
    f = @(x,y) 4*pi^2*(1/L^2+1/H^2)*sin(2*pi*x/L).*sin(2*pi*y/H); 
    %f = @(x,y) 2*ones(size(x));
    %b = M*f(v2(:,1),v2(:,2));
    b = intquadR52(f,t1,v1,t2,v2);
    
    % Conditions limites
    % ------------------
    g = @(x,y) zeros(size(x));
    %g = @(x,y) 0.5*(x.*(L-x)+y.*(H-y));
    ubd = g(v2(ibd2,1),v2(ibd2,2));

    % Traitement des CL
    b=b-A(:,ibd2)*ubd;

    A=A(iin2,iin2);
    b(ibd2)=[];

    % Résolution du système linéaire
    sol=A\b;
    
    % Reconstruction de la solution
    u=zeros(nv2,1);
    u(iin2)=sol;
    u(ibd2)=ubd;

    %--------------------------------------------------------------------------
    % Affichage solution
    %--------------------------------------------------------------------------
    patch('faces',t2(:,1:3),'vertices',v2(:,1:2),'FaceVertexCData',u,...
          'facecolor','interp','edgecolor','black');
    axis equal; 
    colorbar;
    drawnow

    %--------------------------------------------------------------------------
    %  solution exacte
    %--------------------------------------------------------------------------
    %solexact = @(x,y) 0.5*(x.*(L-x)+y.*(H-y));
    solexact = @(x,y) sin(2*pi*x/L).*sin(2*pi*y/H); 
    uexact = solexact(v2(:,1),v2(:,2));

    % erreur max
    erm=max(abs(u-uexact));
    erreurmax=[erreurmax, erm];
    fprintf('erreur max : %e\n',erm);

    % erreur L2
    erL2=erL2P1(solexact,u,t2,v2);
    erreurL2=[erreurL2, erL2];
    fprintf('erreur L2 : %e\n',erL2);
end

% Determination de l'ordre de convergence
loglog(sizemesh,erreurL2,'x-')
xlabel('h'); ylabel('erreur L2');
hold on
loglog(sizemesh,sizemesh.^3,'r-')
legend('erreur L^2','h^3');
[r,m,b]=regression(log(sizemesh),log(erreurL2));
fprintf('\nordre de convergence O(h^m) avec m=%g\n',m);