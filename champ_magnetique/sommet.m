function [test] = sommet(vk1,l)
%test si le triangle consid�r� est un sommet du maillage int�rieur

test = false;
    for i = 1:3
        if ismember([-1 -1],vk1(i,:),'rows') | ismember([1 -1],vk1(i,:),'rows') |ismember([1 1],vk1(i,:),'rows') | ismember([-1 1],vk1(i,:),'rows')
            l(i) = [];
            if size(l(ismember(l,0)))==1
                test = true;
            end
        end
    end

end