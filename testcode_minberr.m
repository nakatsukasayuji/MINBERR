%% Simple test codes for MINBERR for PSD Ax=b
n = 1000; 
A = randn(n); A = A'*A; normA = norm(A); 
b = randn(n,1); 

maxiter = 200; % try run fixed #iter
reorth = 1; 
tol = 0; 
[xsol,iter,berrhist] = minberr(A,b,maxiter,reorth,tol,normA);
[xsol_cg] = pcg(A,b,tol,maxiter);

disp(['backward error minberr: ',num2str(norm(A*xsol-b)/norm(xsol)/normA), ',  pcg: ',num2str(norm(A*xsol_cg-b)/norm(xsol_cg)/normA)])

%% nonsymmetric Ax=b

n = 1000; 
A = randn(n); normA = norm(A); 
b = randn(n,1); 

maxiter = 200; % try run fixed #iter
tol = 0; 
[xsol,iter,berrhist] = minberrNE(A,b,maxiter,tol,normA)
[xsol_lsqr] = lsqr(A,b,tol,maxiter);

disp(['backward error minberr: ',num2str(norm(A*xsol-b)/norm(xsol)/normA), ',  lsqr: ',num2str(norm(A*xsol_cg-b)/norm(xsol_cg)/normA)])
