function [] = write_vtu_unstructured(v,t,Ux,Uy,Uz,Pr,filename)
%***********************************************************************
% Copyright ©2013 Inria - Université de Lorraine
% This software is under GNU General Public Licence
% (http://www.gnu.org/licenses/gpl-3.0.txt)
%***********************************************************************
% VTU Save a 3-D vector array in VTK format.
%  VTU(v,t,Ux,Uy,Uz,pression,filename) save the
%  velocity field(Ux, Uy,Uz) and pressure(pression) in
%  the file filename
% INPUT :
% v : matrix nb_dof x 3 representing the vectors
% t : indices matrices (nb_dof x 8) currently coming from
%     buildQ2_Mesh
% Ux : x-component of the velocity field
% Uy : y-...............................
% Uz : z-...............................
% Pr : pression field
% filename : name of the file to save VTK data

Ux=full(Ux);
Uy=full(Uy);
Uz=full(Uz);
% 
fid = fopen(filename,'w');
sx = size(v,1);
stx = size(t,1);

fprintf(fid,'<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n');
fprintf(fid,'<UnstructuredGrid>\n<Piece NumberOfPoints="%d" NumberOfCells="%d">\n' ,sx,stx);

fprintf(fid,'<PointData Scalars="pression" Vectors="vector_vitesse">\n');

fprintf(fid,'<DataArray type="Float32" Name="vector_vitesse" NumberOfComponents="3" format="ascii">\n');

for i=1:sx ,
	fprintf(fid,'%e %e %e\n',Ux(i),Uy(i),Uz(i)  )    ;
end

fprintf(fid,'</DataArray>');
fprintf(fid,'<DataArray type="Float32" Name="pression" format="ascii">\n');
for i=1:sx ,
	fprintf(fid,'%e\n',Pr(i)  )    ;
end

fprintf(fid,'</DataArray>\n');
fprintf(fid,'</PointData>\n');
fprintf(fid,'<CellData>\n');
fprintf(fid,'</CellData>\n');

fprintf(fid,'<Points>\n');
fprintf(fid,'<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">\n');

for i=1:sx ,
	fprintf(fid,'%f %f %f\n',v(i,1),v(i,2),v(i,3));
end
fprintf(fid,'</DataArray>\n</Points>\n');
fprintf(fid,'<Cells>\n<DataArray type="Int32" Name="connectivity" format="ascii">\n');

for i=1:stx ,
	fprintf(fid,'%d %d %d %d %d %d %d %d\n',t(i,1)-1,t(i,2)-1,t(i,3)-1,t(i,4)-1,t(i,5)-1,t(i,6)-1,t(i,7)-1,t(i,8)-1);
end


fprintf(fid,'</DataArray>\n');
fprintf(fid,'<DataArray type="Int32" Name="offsets" format="ascii">\n');

offset=8;
for i=1:stx ,
	fprintf(fid,'%d ',offset);
	offset=offset+8;
end
fprintf(fid,'\n</DataArray>\n');
fprintf(fid,'<DataArray type="UInt8" Name="types" format="ascii">\n');

for i=1:stx ,
	fprintf(fid,'12 ');
end
fprintf(fid,'\n </DataArray>\n');
fprintf(fid,'</Cells>\n');
fprintf(fid,'</Piece>\n');
fprintf(fid,'</UnstructuredGrid>\n');
fprintf(fid,'</VTKFile>\n');
fclose(fid);

% pvdFile = sprintf('%s/vtk/.fluid.pvd',dirPVD);
% fid2 = fopen(pvdFile, 'a');
% fprintf(fid2,'<DataSet timestep="%f" group="" part="0" file="%s"/>\n',temps,filename);
% fclose(fid2);


end
