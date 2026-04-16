# MINBERR
MINBERR is a minium-backward error Krylov solver that converges universally. 
Given a symmetric positive semidefinite (PSD) linear system Ax=b, after k iterations (using k matrix-vector multiplications), MINBERR computes an approximate solution x_k in the Krylov subspace 
span(b,Ab,...,A^{k-1}b) that minimizes the *backward* error, the smallest perturbation Delta A in A such that (A+Delta A)x_k=b. 


# MINBERR-NE 

Typical usage: 
```
 [x,berrhistory] = minberr(A,b,maxiter,tol)
``` 

# Reference
This project is based on the paper "Towards Universal Convergence of Backward Error  in Linear System Solvers" by Michal Derezinski, Yuji Nakatsukasa, and Elizaveta Rebrova, arXiv April 2026

[keys]( url)
