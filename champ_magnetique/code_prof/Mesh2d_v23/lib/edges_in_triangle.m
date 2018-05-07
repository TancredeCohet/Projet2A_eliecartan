function t2=edges_in_triangle(e,t)
[ne,n]=size(e);
[nt,n]=size(t);
tt=[e(:,4) (1:ne)';e(:,5) (1:ne)']; %matrice ne*2 triangle gauches et droites numero à droite)


%pause
[n,i]=sort(tt(:,1),'ascend'); %trie les triangles par ordre croissant 
tt=tt(i,:); %tris dans l'ordre croissant
i=find(tt(:,1)==0);

tt(i,:)=[];

tt=reshape(tt(:,2),3,nt)';
%reorganisation
ttt=[t(:,2)+t(:,3) t(:,3)+t(:,1) t(:,1)+t(:,2)];
t4=reshape(e(tt,1)+e(tt,2),nt,3);
I=find(t4(:,2)==ttt(:,1));
tmp=tt(I,1);tt(I,1)=tt(I,2);tt(I,2)=tmp;
tmp=t4(I,1);t4(I,1)=t4(I,2);t4(I,2)=tmp;
I=find(t4(:,3)==ttt(:,1));
tmp=tt(I,1);tt(I,1)=tt(I,3);tt(I,3)=tmp;
tmp=t4(I,1);t4(I,1)=t4(I,3);t4(I,3)=tmp;
I=find(t4(:,3)==ttt(:,2));
tmp=tt(I,2);tt(I,2)=tt(I,3);tt(I,3)=tmp;
tmp=t4(I,2);t4(I,2)=t4(I,3);t4(I,3)=tmp;

t2=[t(:,1:5) tt]; %noeud1 n2 n3 where ire cote1 c2 c3

end
