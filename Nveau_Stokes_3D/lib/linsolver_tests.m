function [sol,residu,iter]=linsolver(A,F,nv1,nv2)

    fprintf('Solveur du système linéaire:\n');
    fprintf('   0  UMFPACK(MATLAB)\n');
    fprintf('   1  MINRES\n');
    fprintf('   2  SIMPLE\n');
    fprintf('   3  GCR-SIMPLE\n');
    fprintf('   4  BICGSTAB-SIMPLE\n');
    fprintf('   5  Complément de Schur(MINRES/PCG)\n');
    solver=input('(MINRES par defaut) : ');
    if isempty(solver), solver=1; end;

    if solver==0,       % solveur direct UMFPACK

        disp('Resolution par UMFPACK...');
        tstart = tic;
        spparms('spumoni',1);
        sol=A\F;
    %     r=symamd(A); 
    %     spparms('autommd',0); spparms('autoamd',0); 
    %     sol=A(r,r)\F(r);
    %     sol(r)=sol;
        spparms('spumoni',0);
        fprintf('   temps de resolution : %f s\n',toc(tstart));

    elseif  solver==1   % solveur iteratif MINRES

        disp('Resolution par MINRES (MATLAB)...');
        % choix d'un preconditionnement
        fprintf('preconditionneur:\n');
        fprintf('   0  aucun\n');
        fprintf('   1  diagonal\n');
        fprintf('   2  bloc\n');
        fprintf('   3  ilu\n');
        %fprintf('   4  geometric multigrid\n');
        precond=input('(preconditionnement diagonal par defaut) : ');
        if isempty(precond), precond=1; end;

        tstart = tic;
        if precond==0,        % aucun preconditionnement
            M1=[]; M2=[];
        elseif precond==1,    % preconditionnement diagonal
            B=A(1:3*nv2,3*nv2+1:end-1);
            diagA=diag(A(1:3*nv2,1:3*nv2));
            M1=diag([diagA; diag(B'*diag(sparse(1./diagA))*B); 1]); 
            M2=[];
        elseif precond==2,    % preconditionnement bloc
    %        tstart=tic;
            disp('preconditionnement bloc');
    %        B=A(1:3*nv2,3*nv2+1:end);
    %        D=diag(sparse(1./diag(A(1:3*nv2,1:3*nv2))));
    %         M1=[A(1:3*nv2,1:3*nv2) sparse(3*nv2,nv1+1);
    %             B'                 -B'*D*B]; 
    %         M1(3*nv2+nv1+1,3*nv2+nv1+1)=1;
    %        M1=[];
    %        M2=[speye(3*nv2,3*nv2)  D*B;
    %            sparse(nv1+1,3*nv2) speye(nv1+1,nv1+1)];
    %        fprintf('   temps de construction du preconditionnement:%f\n',toc(tstart));
            M2='m_block3D'; M1 = [];
            mparams = struct('AA',A,'F',F,'nv2',nv2);
        else

            % Iterative solve (ILU preconditioner)

            %tstart=tic;
            %setup.type='nofill';   % to save storage
            %[M1,M2] = ilu(A,setup);
            %fprintf('  time to compute ILU preconditioner is: %d\n',toc(tstart));

            tstart=tic;
            setup.type='nofill';
            [IL,IU]=ilu(A(1:3*nv2,1:3*nv2),setup);
            fprintf('   temps de factorisation ilu:%f\n',toc(tstart));
            M1 = speye(size(A));
            M1(1:3*nv2,1:3*nv2)=IL;
            clear IL;
            M2=M1;
            M2(1:3*nv2,1:3*nv2)=IU;
            clear IU;
            mparams=struct();
            fprintf('   temps de construction du preconditionnement:%f\n',toc(tstart));
        end    

        if precond <2,
            tol=1e-10; itermax=3000;
            [sol,flag,residu,iter] = minres(A,F,tol,itermax,M1,M2);
        else
            %tol=1e-10; itermax=500;
            % % zero initial guess
            %x0=zeros(size(F));
            %[sol,flag,residu,niter] = gmres(A,F,[],tol,itermax,M1,M2,x0,mparams);
            [sol,flag,residu,niter]=gmres(A,F,60,[],[],M1,M2);
            iter=niter(2);
        end
        fprintf('  residu relatif = %6.3e, nombre d''iteration = %d\n',residu,iter);
        fprintf('  temps de resolution : %f s\n',toc(tstart));

        % tstart = tic;
        % disp('Resolution par MINRES...');
        % sol=zeros(length(A(:,1)),1);
        % [sol,tol,iter] = minresiecn(A,F,sol,1e-10);
        % fprintf('  residu relatif = %6.3e, nombre d''iteration = %d\n',tol/norm(F),iter);
        % fprintf('  temps de resolution : %f s\n',toc(tstart));

    elseif solver==2  % Solveur iteratif SIMPLE 

        disp('Resolution par methode SIMPLE');
        tstart = tic;
        tol=1e-6; itermax=4000;
        [sol,residu,iter]=simple(A,F,nv1,nv2,tol,itermax);
        fprintf('  residu relatif = %6.3e, nombre d''iteration = %d\n',residu,iter);
        fprintf('  temps de resolution : %f s\n',toc(tstart));

    elseif solver==3   % solveur iteratif CGR-SIMPLE

        disp('Resolution par methode GCR-SIMPLE');
        tstart = tic;
        restart=0; tol=1e-10; itermax=500; 
        x0=zeros(size(F));
        [sol,residu,iter] = pgcr(A,F,3*nv2,restart,tol,itermax,x0);
        fprintf('  residu relatif = %6.3e, nombre d''iteration = %d\n',residu,iter);
        fprintf('  temps de resolution : %f s\n',toc(tstart));

    elseif solver ==4 % solveur BICGSTAB-SIMPLE

        tstart=tic;
        disp('Resolution par methode BICGSTAB-SIMPLE');
        M2='m_simple'; M1 = [];
        mparams = struct('AA',A,'F',F,'dim2',3*nv2);
        tol=1e-10; itermax=500;
        x0=zeros(size(F)); % zero initial guess
        [sol,flag,residu,iter] = bicgstab(A,F,tol,itermax,M1,M2,x0,mparams); 
        fprintf('  residu relatif = %6.3e, nombre d''iteration = %d\n',residu,iter);
        fprintf('  temps de resolution : %f s\n',toc(tstart));

    else   % solveur iteratif MINRES pour le complément de Schur

        disp('Resolution par complement de SCHUR (MINRES/PCG)');
        tstart = tic;
        disp('Resolution de la pression...');
        pp=zeros(nv1,1);
        tol = 1e-10;
        [pp,residu,iter]=schur_Stokes(A,F,3*nv2,nv1,pp,tol);
        fprintf('  residu relatif de la pression = %6.3e, nombre d''iteration = %d\n',residu,iter);
        fprintf('  temps de resolution pour la pression: %f s\n',toc(tstart));

        tstart = tic;
        disp('Resolution de la vitesse (pcg)...');
        %  Resolution en vitesse : A*uu = f- B*pp
        tol=1e-10; itermax=1000;
        [uu,flag,residu,iter] = pcg(A(1:3*nv2,1:3*nv2),F(1:3*nv2)-A(1:3*nv2,3*nv2+1:end)*pp,tol,itermax);   % X = PCG(A,F,TOL,MAXIT,M1,M2,X0)
        fprintf('  residu relatif de la vitesse = %6.3e, nombre d''iteration = %d\n',residu,iter);
        fprintf('  temps de resolution pour la vitesse (pcg): %f s\n',toc(tstart));
        sol =[uu;pp];
    end
end

