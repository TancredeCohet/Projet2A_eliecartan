function [v2,t2]=P1toP2(v,t,e)
%
% Construction du maillage P2 à partir d'un maillage P1 qui est défini par 
% v (coordonnéees des noeuds), t (triangles/connectivité), e (aretes)
%
    [nv,n]=size(v);
    [nt,n]=size(t);
    x2=(v(e(:,1),1)+v(e(:,2),1))/2;
    y2=(v(e(:,1),2)+v(e(:,2),2))/2;
    w2=e(:,3);
    v2=[v;[x2 y2 w2]];
    t2=[t(:,1) t(:,8)+nv t(:,7)+nv t(:,4)
        t(:,2) t(:,6)+nv t(:,8)+nv t(:,4)
        t(:,3) t(:,7)+nv t(:,6)+nv t(:,4)
        t(:,6)+nv t(:,7)+nv t(:,8)+nv t(:,4)];

end
