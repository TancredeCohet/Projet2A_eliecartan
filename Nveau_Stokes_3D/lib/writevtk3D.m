function [] = writevtk3D(X,Y,Z,Ux,Uy,Uz,Pr,filename)
    
    %  writevtk3D saves a 3D vector/scalar array in VTK structured format.

    [Nxt,Nyt,Nzt]=size(X);
    fid = fopen(filename,'w');
    
    fprintf(fid,'# vtk DataFile Version 3.0\n# Velocity/pressure data\nASCII\n');
    fprintf(fid,'DATASET STRUCTURED_GRID\nDIMENSIONS %d %d %d ',Nxt,Nyt,Nzt);
    fprintf(fid,'POINTS %d DOUBLE\n',Nxt*Nyt*Nzt);
    for k=1:Nzt 
        for j=1:Nyt
            for i=1:Nxt
                fprintf(fid,'%e %e %e\n',X(i,j,k),Y(i,j,k),Z(i,j,k));
            end
        end
    end 
    fprintf(fid,'\nPOINT_DATA %d\n',Nxt*Nyt*Nzt); 
    fprintf(fid,'VECTORS velocity_vectors DOUBLE\n');
    for k=1:Nzt
        for j=1:Nyt
            for i=1:Nxt
                fprintf(fid,'%e %e %e\n',Ux(i,j,k),Uy(i,j,k),Uz(i,j,k));
            end
        end
    end
    
    fprintf(fid,'SCALARS pressure DOUBLE\n');
    fprintf(fid,'LOOKUP_TABLE p\n');
    for k=1:Nzt
        for j=1:Nyt
            for i=1:Nxt
                fprintf(fid,'%e\n',Pr(i,j,k));
            end
        end
    end
   
    fclose(fid);
end
