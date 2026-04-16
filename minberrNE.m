function [xsol,berrhist,y] = minberrNE(A,b,maxiter,tol,normA)
% finds solution with minimum backward error in range(Q) for Ax=b, in the
% subspace Span(A'b,(A'*A)*A'b,(A'*A)^2*A'b,...) using bidiagonalization. 
% NOTE: berrhist is the norm(A*x-b)/norm(x), not berr = norm(A*x-b)/(norm(A)*norm(x))

if nargin<3, maxiter = 100; end
if nargin<4, tol = 0; end
if nargin<5, normA = normest(A,0.1); end

bet = norm(b);
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

    
    if nargout==1    % dqds 
    %     B = spdiags([bet(2:end)' [0 alp(2:end)]' ],[0 1],ii+1,ii+1);
    q = bet(2:end).^2; 
    e = alp(2:end).^2;
    tau = tol^2;
    d = q(1)-tau;
    for ip = 1:numel(q)-1
        qnew(ip) = d + e(ip);
        temp = q(ip+1)/qnew(ip);
        enew(ip) = e(ip)*temp;
        d = d*temp - tau;
        if d<0, % converged
            conv = 1;
            disp(['conv ii=',num2str(ii)]), 
    B = spdiags([bet(2:end)' [0 alp(2:end)]' ],[0 1],ii+1,ii+1);
    Bt = B';
    y = randn(size(B,2),1);

    for ip = 1:100
        ypre = y;
        y = Bt\y;
        y = B\y;
        y = y/norm(y);
        if subspace(y,ypre)<1e-8,
            % enable if desired
            % disp(['invit converged iter=',num2str(ip)]), conv = 1;
            break, end
        end
            
            scale = norm(b)^2/((b'*A)*(V*y));
            xsol = V*scale*y;

            return
%            keyboard
            break
        end
    end

    else % compute berr every step
    B = spdiags([bet(2:end)' [0 alp(2:end)]' ],[0 1],ii+1,ii+1);
    Bt = B';
    y = randn(size(B,2),1);

    for ip = 1:100
        ypre = y;
        y = Bt\y;
        y = B\y;
        y = y/norm(y);
%        if subspace(y,ypre)<1e-8,
        if 1-abs(y'*ypre)<1e-8,
            % enable if desired
            % disp(['invit converged iter=',num2str(ip)]), conv = 1;
            break, end
    end
    berrnow = norm(B*y)/normA; 
    if berrnow<tol, % converged
            disp(['converged: BERR reduced to ',num2str(berrnow), ' at iter=',num2str(ii)]);     break
    end
    berrhist = [berrhist berrnow];
    end

end

if ~exist('y','var')
        B = spdiags([bet(2:end)' [0 alp(2:end)]' ],[0 1],ii+1,ii+1);
    Bt = B';
    y = randn(size(B,2),1);

    for ip = 1:100
        ypre = y;
        y = Bt\y;
        y = B\y;
        y = y/norm(y);
        %if subspace(y,ypre)<1e-8,
        if 1-abs(y'*ypre)<1e-8,
            % enable if desired
            % disp(['invit converged iter=',num2str(ip)]), conv = 1;
            break, end
    end
end
scale = norm(b)^2/((b'*A)*(V*y));
xsol = V*scale*y;

%return
        %B = spdiags([bet(2:end)' [0 alp(2:end)]' ],[0 1],ii+1,ii+1);
 