# MINBERR
MINBERR is a minium-backward error Krylov solver that converges universally. 
Given a symmetric positive semidefinite (PSD) linear system $Ax=b$, after $k$ iterations (using $k$ matrix-vector multiplications), MINBERR computes an approximate solution $x_k$ in the Krylov subspace 
$Span(b,Ab,...,A^{k-1}b)$ that minimizes the *backward* error, the smallest perturbation $\Delta A$ such that $(A+\Delta A)x_k=b$. 
It can be shown (see reference below) that, regardless of $A$,$b$, (i.e., no matter how large the system, or how ill-conditioned the problem), the backward error after $k$ steps is bounded by $3/k^2$. 

# MINBERR-NE 

Typical usage: 
```
 [x,berrhistory] = minberr(A,b,maxiter,tol); % if A is PSD
 [x,berrhistory] = minberrNE(A,b,maxiter,tol); % if A is nonsymmetric
``` 

# Reference
This project is based on the paper "Towards Universal Convergence of Backward Error  in Linear System Solvers" by Michal Derezinski, Yuji Nakatsukasa, and Elizaveta Rebrova, arXiv April 2026

[keys]( url)
