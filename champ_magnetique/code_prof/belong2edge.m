function [bool]=belong2edge(v,t,k);
    xGt=1/3*sum(v(t(k,:),1));%calcul du centre de gravité du triangle k
    yGt=1/3*sum(v(t(k,:),2));
    d1=[abs(xGt-1),1];%distance à la droite x=-1
    d2=[abs(xGt+1),2];%distance à la droite x=1
    d3=[abs(yGt-1),3];%distance à la droite y=-1
    d4=[abs(yGt-1),4];%distance à la droite y=1
    d=min(d1,d2,d3,d4);
    a=find([d1 d2 d3 d4]==d);
    if a==1;
        x=-1;%projection du centre de gravité sur le côté x=-1
        y=yGt;
        cpt=0;
        for i=1:3;
            a=(v(t(k,mod(i+1,3),2))-v(t(k,mod(i,3),2)))/(v(t(k,mod(i+1,3),1)-v(t(k,mod(i,3),1))));
            b=(v(t(k,mod(i,3),1)*v(t(k,mod(i+1,3),2))-v(t(k,mod(i+1,3),1)*v(t(k,mod(i,3),2)))/(v(t(k,mod(i+1,3),1)-v(t(k,mod(i,3),1))))));
            if y==a*x+b;
                cpt=cpt+1;
            end
        end
        if cpt>=1;
            bool=true
        else
            bool=false;
        end
    elseif a==2;
        x=1;
        y=yGt;
        cpt=0;
        for i=1:3;
            a=(v(t(k,mod(i+1,3),2))-v(t(k,mod(i,3),2)))/(v(t(k,mod(i+1,3),1)-v(t(k,mod(i,3),1))));
            b=(v(t(k,mod(i,3),1)*v(t(k,mod(i+1,3),2))-v(t(k,mod(i+1,3),1)*v(t(k,mod(i,3),2)))/(v(t(k,mod(i+1,3),1)-v(t(k,mod(i,3),1))))));
            if y==a*x+b;
                cpt=cpt+1;
            end
        end
        if cpt>=1;
            bool=true
        else
            bool=false;
        end
    elseif a==3;
        cpt=0;
        for i=1:3;
            a=(v(t(k,mod(i+1,3),2))-v(t(k,mod(i,3),2)))/(v(t(k,mod(i+1,3),1)-v(t(k,mod(i,3),1))));
            b=(v(t(k,mod(i,3),1)*v(t(k,mod(i+1,3),2))-v(t(k,mod(i+1,3),1)*v(t(k,mod(i,3),2)))/(v(t(k,mod(i+1,3),1)-v(t(k,mod(i,3),1))))));
            if y==a*x+b;
                cpt=cpt+1;
            end
        end
        if cpt>=1;
            bool=true
        else
            bool=false;
        end
    else
        cpt=0;
        for i=1:3;
            a=(v(t(k,mod(i+1,3),2))-v(t(k,mod(i,3),2)))/(v(t(k,mod(i+1,3),1)-v(t(k,mod(i,3),1))));
            b=(v(t(k,mod(i,3),1)*v(t(k,mod(i+1,3),2))-v(t(k,mod(i+1,3),1)*v(t(k,mod(i,3),2)))/(v(t(k,mod(i+1,3),1)-v(t(k,mod(i,3),1))))));
            if y==a*x+b;
                cpt=cpt+1;
            end
        end
        if cpt>=1;
            bool=true
        else
            bool=false;
        end
    end
end

        
        
        
           
        
    