function t2=area_triangle(v,t);
    t2=t(:,1:3)';
    nt=size(t,1);
    x=reshape(v(t2,1),3,nt);
    y=reshape(v(t2,2),3,nt);
    x(2,:)=x(2,:)-x(1,:);
    x(3,:)=x(3,:)-x(1,:);
    y(2,:)=y(2,:)-y(1,:);
    y(3,:)=y(3,:)-y(1,:);
    A=(x(2,:).*y(3,:)-x(3,:).*y(2,:))/2;
    t2=[t(:,1:4) A'];
end
