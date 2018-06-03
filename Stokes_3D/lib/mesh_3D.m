function [vnodes,t,nb_vertex,nb_dof,nb_cubes,X,Y,Z] = mesh_3D(Lx,Ly,Lz,Nx,Ny,Nz)
%**************************************************************************
% Build 3D mesh of the rectangular cuboid  [0,Lx] x [0,Ly] x [0,Lz]
% with respectively NX, Ny and Nz discrete points in each corresponding
% direction.
%
%    vnodes : vector containing the mesh points for finite elements Q1(pressure)
%             and Q2(velocity)
%         t : connectivity table, indices vector to access to the points of 
%             vnodes belonging to an elementary cube : vnodes(t(i,k)) is the 
%             k-th point of the cube i (k in {1,...,17}, i in {1...nb_cubes})
% nb_vertex : number of degre of freedom for the pressure
%    nb_dof : number of degres of freedom for one component of the velocity
%  nb_cubes : number of elementary "cubes"
%         X : 3-tensors for the mesh along the X axis
%         Y : 3-tensors for the mesh along the Y axis
%         Z : 3-tensors for the mesh along the Z axis
%**************************************************************************

% mesh grid
xl=linspace(0,Lx,Nx);
yl=linspace(0,Ly,Ny);
zl=linspace(0,Lz,Nz);
[X,Y,Z]=ndgrid(xl,yl,zl);

% mesh counting
nb_vertex=Nx*Ny*Nz;
nb_edges=(Nx-1)*Ny*Nz+Nx*(Ny-1)*Nz+Nx*Ny*(Nz-1);
nb_faces=Nx*(Ny-1)*(Nz-1)+(Nx-1)*Ny*(Nz-1)+(Nx-1)*(Ny-1)*Nz;
nb_cubes=(Nx-1)*(Ny-1)*(Nz-1);
nb_dof = nb_vertex + nb_edges + nb_faces + nb_cubes;

% the first nb_vertex components of vnodes contain the mesh points for the
% pressure, i.e. all the vertex of the mesh, unfolded along the matlab convention :
% X is the most variant component :
% thus we have X(1,1,1)->x(2,1,1)->...->X(Nx,2,1)->...->X(Nx,NY,2)
vnodes=[X(:),Y(:),Z(:)];
% indexes
[ix,iy,iz]=ndgrid(1:Nx-1,1:Ny-1,1:Nz-1);
first_point = [ix(:),iy(:)-1,iz(:)-1] *[1;Nx;Nx*Ny];
% store the indices for the vertex of the cubes (Q1 mesh)
t=[first_point ,...
	first_point + 1 ,...
	first_point + Nx + 1 ,...
	first_point + Nx, ...
	first_point + Nx*Ny,...
	first_point + Nx*Ny + 1,...
	first_point + Nx*Ny + Nx + 1,...
	first_point + Nx*Ny + Nx];

% to store the  coordinates
coords=cell(27,1);
indices=cell(27,1);

%point 9 : middle of 1 and 2
coords{9}  = 0.5*( vnodes(t(:, 1), :) + vnodes(t(:, 2), :) );

%point 12 : middle of 1 and 4
coords{12} = 0.5*( vnodes(t(:, 1), :) + vnodes(t(:, 4), :) );

%point 13 : middle of 1 and 5
coords{13} = 0.5*( vnodes(t(:, 1), :) + vnodes(t(:, 5), :) );

% we add the three computed points to vnodes
vnodes = [vnodes; coords{9}; coords{12}; coords{13} ];

% add indices for point 9
int_cubes = (1:nb_cubes)';
shift_9  = nb_vertex;
indices{9} = shift_9 + int_cubes;

% add indices for point 12
shift_12 = shift_9 + nb_cubes;
indices{12} = shift_12 + int_cubes;

% add indices for point 13
shift_13 = shift_9 + 2 * nb_cubes;
indices{13} = shift_13 + int_cubes;

% temp function to compute the shift index
shift_ind=@(x,y,z) x + (Nx-1)*y + (Nx-1)*(Ny-1)*z;

% t[11,(i,j,k)]=t[9,(i,j+1,k)]
indices{11} = shift_9  + int_cubes + shift_ind(0, 1 ,0);

% t[10,(i,j,k)]=t[12,(i+1,j,k)]
indices{10} = shift_12 + int_cubes + shift_ind(1, 0, 0);

% t[14,(i,j,k)]=t[13,(i+1,j,k)]
indices{14} = shift_13 + int_cubes + shift_ind(1, 0, 0);

% t[15,(i,j,k)]=t[13,(i+1,j+1,k)]
indices{15} = shift_13 + int_cubes + shift_ind(1, 1, 0);

% t[16,(i,j,k)]=t[13,(i,j+1,k)]
indices{16} = shift_13 + int_cubes + shift_ind(0, 1, 0);

% t[17,(i,j,k)]=t[9,(i,j,k+1)]
indices{17} = shift_9  + int_cubes + shift_ind(0, 0, 1);

% t[18,(i,j,k)]=t[12,(i+1,j,k+1)]
indices{18} = shift_12 + int_cubes + shift_ind(1, 0, 1);

% t[19,(i,j,k)]=t[9,(i,j+1,k+1)]
indices{19} = shift_9  + int_cubes + shift_ind(0, 1, 1);

% t[20,(i,j,k)]=t[12,(i,j,k+1)]
indices{20} = shift_12 + int_cubes + shift_ind(0, 0, 1);

%append to t the computed indices
t = [t , [indices{9:20}]];

%build points 21,22,25 on the right values of 9, 12 and 13
coords{21} = vnodes(t(:, 9), :) + vnodes(t(:, 12), :) - vnodes(t(:, 1), :);
coords{22} = vnodes(t(:, 9), :) + vnodes(t(:, 13), :) - vnodes(t(:, 1), :);
coords{25} = vnodes(t(:, 12), :) + vnodes(t(:, 13), :) - vnodes(t(:, 1), :);
vnodes = [vnodes ; coords{21}; coords{22}; coords{25}];

% add indices for point 21
shift_21 = shift_9 + 3 * nb_cubes;
indices{21} = shift_21 + int_cubes;

% add indices for  point 22
shift_22 = shift_21 + nb_cubes;
indices{22} = shift_22 + int_cubes;

% add indices for  point 25
shift_25 = shift_22 + nb_cubes;
indices{25} = shift_25 + int_cubes;

% t[23,(i,j,k)]=t[25,(i+1,j,k)]
indices{23} = shift_25 + int_cubes + shift_ind(1, 0, 0);

% t[24,(i,j,k)]=t[22,(i,j+1,k)]
indices{24} = shift_22 + int_cubes + shift_ind(0, 1, 0);

% t[26,(i,j,k)]=t[21,(i,j,k+1)]
indices{26} = shift_21 + int_cubes + shift_ind(0, 0, 1);

%append to t the computed indices
t = [t, [indices{21:26}]];

% point 27 = point 9 + point 12 + point 13 - 3 * point 1 + point 1
coords{27}  = vnodes(t(:, 9), :)+vnodes(t(:, 12), :)+vnodes(t(:, 13),:)...
	-2*vnodes(t(:, 1), :);
shift_27    = shift_21 + 3*nb_cubes;
indices{27} = shift_27 + int_cubes;

vnodes=[vnodes; coords{27}];
t=[t, indices{27}];

% Now we need corrections for limits case

% a inline function to flatten the 3D-arrays
% flat(x) ~ x(:)
flat  =  @(x) reshape(x,numel(x),1,1);

% Face Nx-1 :
[ix,iy,iz]  = ndgrid(Nx-1,1:Ny-1,1:Nz-1);
false_ind = flat(shift_ind(ix,iy-1,iz-1));
nb_nx_1 = size(false_ind,1);
int_false = (1:nb_nx_1)';

%for point 10
shift_10_correc = shift_27 + nb_cubes;
coords{10} = 0.5*(vnodes(t(false_ind, 2), :) + vnodes(t(false_ind, 3), :));
t(false_ind, 10) = shift_10_correc + int_false;

%for point 14
shift_14_correc = shift_10_correc + nb_nx_1;
coords{14} = 0.5*(vnodes(t(false_ind, 2), :) + vnodes(t(false_ind, 6), :));
t(false_ind, 14) = shift_14_correc + int_false;

% for point 23
shift_23_correc = shift_14_correc + nb_nx_1 ;
coords{23} = 0.25*(vnodes(t(false_ind, 2), :)+vnodes(t(false_ind, 3), :) + ...
	vnodes(t(false_ind, 6), :)+vnodes(t(false_ind, 7), :));
t(false_ind, 23) = shift_23_correc + int_false;

%append to vnodes the the corrected vectors
vnodes = [vnodes; coords{10}; coords{14};  coords{23} ];

%for point 15
[ix,iy,iz]  = ndgrid(Nx-1,1:Ny-2,1:Nz-1);
false_ind_p = flat(shift_ind(ix,iy-1,iz-1));
t(false_ind_p, 15)= t(false_ind_p + shift_ind(0,1,0), 14);

%for point 18
[ix,iy,iz]  = ndgrid(Nx-1,1:Ny-1,1:Nz-2);
false_ind_p = flat(shift_ind(ix,iy-1,iz-1));
t(false_ind_p, 18)= t(false_ind_p + shift_ind(0,0,1), 10);

% Face Ny-1 :
[ix,iy,iz]  = ndgrid(1:Nx-1,Ny-1,1:Nz-1);
false_ind = flat(shift_ind(ix,iy-1,iz-1));
nb_ny_1 = size(false_ind,1);
int_false = (1:nb_ny_1)';

%for point 16
shift_16_correc = shift_10_correc + 3 * nb_nx_1;
coords{16} = 0.5*(vnodes(t(false_ind, 4), :) + vnodes(t(false_ind, 8), :));
t(false_ind, 16) = shift_16_correc + int_false;

%for point 11
shift_11_correc = shift_16_correc + nb_ny_1;
coords{11} = 0.5*(vnodes(t(false_ind, 3), :)+vnodes(t(false_ind, 4), :));
t(false_ind, 11)= shift_11_correc + int_false;

% for point 24
shift_24_correc = shift_11_correc + nb_ny_1;
coords{24} = 0.25*(vnodes(t(false_ind, 3), :)+vnodes(t(false_ind, 4), :) + ...
	vnodes(t(false_ind, 7), :)+vnodes(t(false_ind, 8), :));
t(false_ind, 24)= shift_24_correc + int_false;

vnodes = [vnodes;  coords{16}; coords{11}; coords{24}];

%for point 15
[ix,iy,iz]  = ndgrid(1:Nx-2,Ny-1,1:Nz-1);
false_ind_p = flat(shift_ind(ix,iy-1,iz-1));
t(false_ind_p, 15) = t(false_ind_p + shift_ind(1,0,0), 16);

%for point 19
[ix,iy,iz]  = ndgrid(1:Nx-1,Ny-1,1:Nz-2);
false_ind_p = flat(shift_ind(ix,iy-1,iz-1));
t(false_ind_p, 19) = t(false_ind_p + shift_ind(0,0,1), 11);

% Face Nz-1:
[ix,iy,iz]  = ndgrid(1:Nx-1,1:Ny-1,Nz-1);
false_ind = flat(shift_ind(ix,iy-1,iz-1));
nb_nz_1 = size(false_ind,1);
int_false = (1:nb_nz_1)';

%for point 17
shift_17_correc = shift_16_correc + 3 * nb_ny_1;
coords{17} = 0.5*(vnodes(t(false_ind, 5), :) + vnodes(t(false_ind, 6), :));
t(false_ind, 17) = shift_17_correc + int_false;

%for point 20
shift_20_correc = shift_17_correc + nb_nz_1;
coords{20} = 0.5*(vnodes(t(false_ind, 5), :)+vnodes(t(false_ind, 8), :));
t(false_ind, 20) = shift_20_correc + int_false;

% for point 26
shift_26_correc = shift_20_correc + nb_nz_1;
coords{26} = 0.25*(vnodes(t(false_ind, 5), :) + vnodes(t(false_ind, 6), :) + ...
	vnodes(t(false_ind, 7), :) + vnodes(t(false_ind, 8), :));
t(false_ind, 26) = shift_26_correc + int_false;

vnodes =[vnodes;  coords{17}; coords{20}; coords{26}];

%for point 18
[ix,iy,iz]  = ndgrid(1:Nx-2,1:Ny-1,Nz-1);
false_ind_p = flat(shift_ind(ix, iy-1, iz-1));
t(false_ind_p, 18)= t(false_ind_p + shift_ind(1,0,0), 20);

%for point 19
[ix,iy,iz]  = ndgrid(1:Nx-1,1:Ny-2,Nz-1);
false_ind_p = flat(shift_ind(ix, iy-1, iz-1));
t(false_ind_p, 19)= t(false_ind_p + shift_ind(0,1,0), 17);

% for corrections it remains the last three edges

%edge (Ny-1)x(Nx-1)
[ix,iy,iz]  = ndgrid(Nx-1,Ny-1,1:Nz-1);
false_ind = flat(shift_ind(ix,iy-1,iz-1));
nb_nx_ny = size(false_ind,1);
int_false = (1:nb_nx_ny)';

%for point 14
shift_15_correc_p = shift_17_correc + 3 * nb_nz_1;
coords{15} = 0.5*(vnodes(t(false_ind, 3), :) + vnodes(t(false_ind, 7),:));
t(false_ind, 15) = shift_15_correc_p + int_false;

vnodes = [vnodes;  coords{15}  ];

%edge (Ny-1)x(Nz-1)
[ix,iy,iz]  = ndgrid(1:Nx-1,Ny-1,Nz-1);
false_ind = flat(shift_ind(ix,iy-1,iz-1));
nb_ny_nz = size(false_ind,1);
int_false =(1:nb_ny_nz)';

%for point 19
shift_19_correc_p = shift_15_correc_p + nb_nx_ny;
coords{19} = 0.5*(vnodes(t(false_ind, 7), :) + vnodes(t(false_ind, 8), :));
t(false_ind, 19) = shift_19_correc_p + int_false;

vnodes = [vnodes; coords{19}];

%edge (Nx-1)x(Nz-1)
[ix,iy,iz]  = ndgrid(Nx-1,1:Ny-1,Nz-1);
false_ind = flat(shift_ind(ix,iy-1,iz-1));
nb_nx_nz = size(false_ind,1);
int_false = (1:nb_nx_nz)';

%for point 18
shift_18_correc_p = shift_19_correc_p + nb_ny_nz;
coords{18} = 0.5*(vnodes(t(false_ind, 6), :) + vnodes(t(false_ind, 7), :));
t(false_ind,18)= shift_18_correc_p + int_false;

vnodes = [vnodes ; coords{18} ];

end
