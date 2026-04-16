function [xsol,iter,berrhist,V,y] = minberrNE(A,b,maxiter,tol,normA)
% Solves Ax=b for a general nxn A. 
% MINBERRNE finds the minimizer of the backward error in the normal Krylov subspace
% Span(A'b,(A'*A)*A'b,(A'*A)^2*A'b,...) using bidiagonalization. 
% Unlike MINBERR for PSD systems, *universal* convergence is not established, 
% but it is observed to do so like 1/k up to a log(cond(A)) factor. 
% Based on paper "Towards Universal Convergence of Backward Error in 
% Linear System Solvers" by Michal Derezinski, Yuji Nakatsukasa, and 
% Elizaveta Rebrova, arXiv April 2026.

if nargin<3, maxiter = 100; end % maximum #iterations to run
if nargin<4, tol = 0; end       % if berr<tol, stop iteration
if nargin<5, normA = normest(A,0.1); end
% normA is needed to compute berr =||Ax-b||/norm(A)/norm(x)

tau = tol^2;  % dqds tolerance
delta = 1e-5; % failure probability
numinvit = @(k,delta) 100; % works in most cases
%numinvit = @(k,delta) 2.23*log(k/delta^2); % if 1/2 relative berr is sufficient

if nargout < 3
    computeberr = 0; % don't compute berr at each step, bit faster
else
    computeberr = 1; 
end

bet = norm(b); % start Golub-Kahan bidiagonalization
u = b/bet; U = [u];
vv = A'*u;
alp = norm(vv);
V = vv/alp(end);

vnow = V(:,end);
unext = A*vnow - alp(end)*U(:,end);
bet = [bet norm(unext)];
unext = unext/norm(unext);
berrhist = [];
for ii = 1:maxiter
    U = [U unext];
    vnext = A'*unext - bet(end)*vnow;
    alp = [alp norm(vnext)];
    V = [V vnext/norm(vnext)];
    vnow = V(:,end);
    unext = A*vnow - alp(end)*U(:,end);
    bet = [bet norm(unext)];
    unext = unext/norm(unext);

    if computeberr==0 % dqds to test for convergence        
        if ii == 1
            q = bet(2:end).^2;
            e = alp(2:end).^2;
            d = q(1)-tau;
        else
            qnew = d + alp(end)^2;
            temp = bet(end)^2/qnew;
            d = d*temp - tau;
        end
            if d<0 | ii == maxiter % converged or terminate
                iter = ii;
                % disp(['converged ii=',num2str(ii)]),
                B = spdiags([bet(2:end)' [0 alp(2:end)]' ],[0 1],ii+1,ii+1);
                Bt = B';                
                y = invit_forminsvd(B,numinvit(size(B,2),delta));

                scale = norm(b)^2/((b'*A)*(V*y));
                xsol = V*scale*y;     return 
            end

    %{
    if computeberr==0 % dqds to test for convergence
        q = bet(2:end).^2;
        e = alp(2:end).^2;
        tau = tol^2;
        d = q(1)-tau;
        for ip = 1:numel(q)-1
            qnew(ip) = d + e(ip);
            temp = q(ip+1)/qnew(ip);
            enew(ip) = e(ip)*temp;
            d = d*temp - tau;
            if d<0 | ii == maxiter % converged or terminate
                iter = ii;
                % disp(['converged ii=',num2str(ii)]),
                B = spdiags([bet(2:end)' [0 alp(2:end)]' ],[0 1],ii+1,ii+1);
                Bt = B';                
                y = invit_forminsvd(B,numinvit(size(B,2),delta));

                scale = norm(b)^2/((b'*A)*(V*y));
                xsol = V*scale*y;     return 
            end
        end
    %}
    else % compute berr every step
        B = spdiags([bet(2:end)' [0 alp(2:end)]'],[0 1],ii+1,ii+1);
        Bt = B';
        y = invit_forminsvd(B,numinvit(size(B,2),delta));
        berrnow = norm(B*y)/normA;
        if berrnow<tol, % converged
        %    disp(['converged: BERR reduced to ',num2str(berrnow), ' at iter=',num2str(ii)]);     
        break
        end
        berrhist = [berrhist berrnow];
    end
end

scale = norm(b)^2/((b'*A)*(V*y));
xsol = V*scale*y;

iter = ii;
end

function y = invit_forminsvd(B,maxit) 
% inverse iteration to find smallest right singular vector of B
y = randn(size(B,2),1);
Bt = B'; 
for ip = 1:maxit
    ypre = y;
    y = Bt\y;
    y = B\y;
    y = y/norm(y);
    if subspace(y,ypre) < 1e-8 
        % this angular tolerance for invit (implies O(1e-15) convergence in singval)
        % disp(['invit converged iter=',num2str(ip)]), conv = 1;
        break, end
    end
end