function tQ1 = meshQ1isoQ2(t)
    nbv = size(t,1);
    tQ1 = zeros(8*nbv,8);
    k1=0;
    for k=1:nbv
        sk = t(k,1:27);
        k1=k1+1;
        tQ1(k1,:)=[sk(1),sk(9),sk(21),sk(12),sk(13),sk(22),sk(27),sk(25)];
        k1=k1+1;
        tQ1(k1,:)=[sk(9),sk(2),sk(10),sk(21),sk(22),sk(14),sk(23),sk(27)];
        k1=k1+1;
        tQ1(k1,:)=[sk(21),sk(10),sk(3),sk(11),sk(27),sk(23),sk(15),sk(24)];
        k1=k1+1;
        tQ1(k1,:)=[sk(12),sk(21),sk(11),sk(4),sk(25),sk(27),sk(24),sk(16)];
        k1=k1+1;
        tQ1(k1,:)=[sk(13),sk(22),sk(27),sk(25),sk(5),sk(17),sk(26),sk(20)];
        k1=k1+1;
        tQ1(k1,:)=[sk(22),sk(14),sk(23),sk(27),sk(17),sk(6),sk(18),sk(26)];
        k1=k1+1;
        tQ1(k1,:)=[sk(27),sk(23),sk(15),sk(24),sk(26),sk(18),sk(7),sk(19)];
        k1=k1+1;
        tQ1(k1,:)=[sk(25),sk(27),sk(24),sk(16),sk(20),sk(26),sk(19),sk(8)];
    end
    
end

