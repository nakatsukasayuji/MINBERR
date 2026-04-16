function [xsol,iter,berrhist] = minberr(A,b,maxiter,reorth,tol,normA)
% Solves Ax=b for a positive semidefinite nxn A.
% MINBERR finds the minimizer of the backward error in the rylov subspace
% Span(b,A*b,A^2*b,...) using bidiagonalization.
% Universal convergence holds: the backward error after k iterations is
% bounded by 3/k^2, no matter what A,b are given.
%
% Based on paper "Towards Universal Convergence of Backward Error in
% Linear System Solvers" by Michal Derezinski, Yuji Nakatsukasa, and
% Elizaveta Rebrova, arXiv April 2026.

if nargin<4, reorth = 0; end
if nargin<5, tol = 0; end
if nargin<6, normA = normest(A,0.1); end
% normA is needed to compute berr =||Ax-b||/norm(A)/norm(x)

delta = 1e-10; % failure probability
numinvit = @(k,delta) 100; % works in most cases
%numinvit = @(k,delta) 2.23*log(k/delta^2); % if 1/2 relative berr is sufficient

n = size(A,1);
Q = b/norm(b); qend = Q;
qnew = A*qend;
y = [1];
Tmain = zeros(1,n);
Toff = zeros(1,n);
berrhist = [];
for ii = 1:maxiter
    if ii==1
        coeff = qend'*qnew; Tmain(ii) = coeff(end);
        qnew = qnew - Q*(coeff);
    else
        qnew = qnew - Q(:,end-1:end)*(coeff);
    end
    if reorth
        qnew = qnew - Q*(Q'*qnew); % if stability desired, reorthogonalize (bit expensive but recommended)
    end
    Toff(ii) = norm(qnew);
    Q = [Q qnew/Toff(ii)];
    qend = Q(:,end);
    qnew = A*qend;
    coeff = Q(:,end-1:end)'*qnew; Tmain(ii+1) = coeff(end);

    %T = spdiags([[Toff(1:ii) 0]' Tmain(1:ii+1)' [0 Toff(1:ii)]'],[-1 0 1],ii+1,ii+1);% below is faster
    if ii==1
        T = spdiags([[Toff(1:ii) 0]' Tmain(1:ii+1)' [0 Toff(1:ii)]'],[-1 0 1],ii+1,ii+1);
    else
        T = blkdiag(T,Tmain(ii+1));
        T(end,end-1) = Toff(ii); T(end-1,end) = Toff(ii);
    end

    % we could do chol(T'*T-eps^2*speye(size(T,2))) to test convergence fast,
    % but for stability we compute berr every step.
    % (usually less expensive than reorth, which is recommended anyway)
    [berr,y] = minberrfromLanczos(T,numinvit(ii,delta));
    berrhist = [berrhist berr/normA];

    if berr/normA < tol, % converged
        disp(['converged: BERR reduced to ',num2str(berr)]);
        break
    end

end
iter = ii;
xsol = Q(:,1:numel(y))*y*(norm(b)/(T(1,1:2)*y(1:2)));
end


function [berr,y] = minberrfromLanczos(T,maxit)
% input: T s.t. Lanczos AQ=QT, initial y for inv-power
% output backward error ||Ax-b||/(||A||*||x||)
% to compute soln x,
% xsol = Q(:,1:numel(x))*x*(normb/(T(1,1:2)*x(1:2)));
TT = T(2:end,1:end-1);
y = randn(size(TT,1),1);
TTT = TT'; % keyboard
for ip = 1:maxit % inv power iteration
    ypre = y;
    y = TTT\y;
    y = TT\y;
    y = y/norm(y);
    if subspace(y,ypre)<1e-8,
        break
    end
end
berr = norm(TT*y);
end
