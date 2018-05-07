function [sol,residu,iter]=linsolver(A,F,nv1,nv2)
        
        fprintf('*/ Resolution du systeme lineaire');
        % choix d'un preconditionneur/solveur
        % precond=0 : MINRES sans preconditionneur
        % precond=1 : MINRES avec preconditionneur diagonal
        % precond=2 : MINRES avec preconditionneur ILU
        % precond=3 : GMRES  avec preconditionneur ILU
        precond =1;
        
        tstart1=tic;
        
        if precond==0     % aucun preconditionnement
            fprintf(' (MINRES sans preconditionneur)\n');
            M1=[]; M2=[];
        elseif precond==1 % preconditionnement diagonal
            fprintf(' (MINRES avec preconditionneur diagonal)\n');
            B=A(1:3*nv2,3*nv2+1:end-1);
            diagA=diag(A(1:3*nv2,1:3*nv2));
            M1=diag([diagA; diag(B'*diag(sparse(1./diagA))*B); 1]); 
            M2=[];
        elseif precond==2  % preconditionnement ILU du Laplacien pour MINRES
                           % le preconditionnement de MINRES (MATLAB) doit
                           % etre symetrique defini positif.
            
            fprintf(' (MINRES avec preconditionneur ILU)\n');
            
            setup.type='nofill';
            [IL,IU]=ilu(A(1:3*nv2,1:3*nv2),setup);
            fprintf('     temps de factorisation ilu : %f s\n',toc(tstart1));
            
            tstart2=tic;
            Msp1=CreateSparseMatrix(3*nv2+nv1,3*nv2+nv1);
            indice1=[1:3*nv2];
            AddMatElem(Msp1,indice1,indice1,IL);
            % On rajoute l'identite sur les d.o.f de la pression
            indice2=3*nv2+[1:nv1];
            AddMatElem(Msp1,indice2,indice2,speye(nv1,nv1));
            M1=SparseMatrixToMatlab(Msp1);
            DeleteSparseMatrix(Msp1);
            clear IL;
            %
            Msp2=CreateSparseMatrix(3*nv2+nv1,3*nv2+nv1);
            AddMatElem(Msp2,indice1,indice1,IU);
            AddMatElem(Msp2,indice2,indice2,speye(nv1,nv1));
            M2=SparseMatrixToMatlab(Msp2);
            DeleteSparseMatrix(Msp2);
            clear IU;
            fprintf('     passage du format ALICE->MATLAB (%f s)\n',toc(tstart2));
            fprintf('     temps total de construction du preconditionnement: %f s\n',toc(tstart1));
            
            %La version directe ("equivalente") suivante ne fonctionne pas ! 
            %M1 = speye(size(A));
            %M1(1:3*nv2,1:3*nv2)=IL; %!!! Ne fonctionne pas !!! -> Out of memory 
            %clear IL;
            %M2=M1;
            %M2(1:3*nv2,1:3*nv2)=IU;
            %clear IU;
            
        else  % preconditionnement ILU de la matrice de Stokes penalisee 
              % par epsilon*pression pour GMRES  
            
            fprintf(' (GMRES  avec preconditionneur ILU)\n');
            
            setup.type='nofill';
            Asp=CreateSparseMatrix(3*nv2+nv1,3*nv2+nv1);
            indice1=[1:3*nv2+nv1];
            AddMatElem(Asp,indice1,indice1,A);
            indice2=3*nv2+[1:nv1];
            epsil=1e-6;
            AddMatElem(Asp,indice2,indice2,epsil*speye(nv1,nv1));
            Aepsil=SparseMatrixToMatlab(Asp);
            DeleteSparseMatrix(Asp);
            fprintf('     passage du format ALICE->MATLAB (%f s)\n',toc(tstart1));
            
            tstart2=tic;
            setup.type='nofill';
            [M1,M2]=ilu(Aepsil,setup);
            clear Aepsil;
            fprintf('     temps de factorisation ilu : %f s\n',toc(tstart2));
            fprintf('     temps total de construction du preconditionnement:%f s\n',toc(tstart1));
        end    

        fprintf('     resolution du systeme lineaire...\n');
   
        tol=1e-10; itermax=3000; 
        if precond<=2,  % solveur MINRES (preconditionneur s.d.p)
           [sol,flag,residu,iter] = minres(A,F,tol,itermax,M1,M2);
           fprintf('     residu relatif = %6.3e, nombre d''iteration = %d\n',residu,iter);
        else            % solveur GMRES
           restart=60;
           [sol,flag,residu,iter] = gmres(A,F,restart,tol,itermax,M1,M2);
           fprintf('     residu relatif = %6.3e, nombres d''iterations = %d (outer), %d (inner)\n',residu,iter(1),iter(2));
        end

        fprintf('     temps total de resolution du systeme lineaire : %f s\n',toc(tstart1));

end

