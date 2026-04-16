function [xsol,berrhist,y] = minberr(A,b,maxiter,reorth,tol)
% Solves Ax=b for a positive semidefinite A. 
% MINBERR finds the minimizer of the backward error in the Krylov subspace.
% It converges *universally* at least as 1/k^2 after k iterations. 

if nargin<4, reorth = 0; end
if nargin<5, tol = 1e-16; end

n = size(A,1);
Q = b/norm(b); qend = Q; 
qnew = A*qend; 
y = [1];
Tmain = zeros(1,n);
Toff = zeros(1,n);
berrhist = [];
for ii = 1:maxiter
    %disp(ii)
    if ii==1
    coeff = qend'*qnew; Tmain(ii) = coeff(end); 
    qnew = qnew - Q*(coeff); 
    else
    qnew = qnew - Q(:,end-1:end)*(coeff); 
    end
    if reorth
    qnew = qnew - Q*(Q'*qnew); % if desired, reorthogonalize    
    end
    Toff(ii) = norm(qnew);     
    Q = [Q qnew/Toff(ii)];
    qend = Q(:,end);         
    qnew = A*qend; 
    coeff = Q(:,end-1:end)'*qnew; Tmain(ii+1) = coeff(end);     
    
    T = diag([Tmain(1:ii+1)]) + diag(Toff(1:ii),1) + diag(Toff(1:ii),-1);

    %[berr,y] = minberrfromLanczos_dqds(T);
    
   
    if ii>1
    %[berr,y] = minberrfromLanczos(T,[y;0]);
    [berr,y] = minberrfromLanczos(T);
    else % 1st iter
    [berr,y] = minberrfromLanczos(T,1);        
    xsol = Q(:,1:numel(y))*y*(norm(b)/(T(1,1)*y(1)));    
    end    
    

    berrhist = [berrhist berr];

    if berr<tol, % converged
            disp(['converged: BERR reduced to ',num2str(berr)]);     break
    end

end
    xsol = Q(:,1:numel(y))*y*(norm(b)/(T(1,1:2)*y(1:2)));    
end

function [berr,y] = minberrfromLanczos(T,y)
% input: T s.t. Lanczos AQ=QT, initial y for inv-power
% output backward error ||Ax-b||/(||A||*||x||)
% to compute soln x, 
% xsol = Q(:,1:numel(x))*x*(normb/(T(1,1:2)*x(1:2)));
    TT = T(2:end,1:end-1);     
    if nargin < 2, 
    y = randn(size(TT,1),1);
    end
    TTT = TT'; % keyboard
    for ip = 1:100 % inv power iteration
    ypre = y; 
    y = TTT\y; 
    y = TT\y; 
    y = y/norm(y);     
    if isnan(y)|isinf(y), keyboard, end
    if subspace(y,ypre)<1e-8, 
        % enable if desired:
            % disp(['invit converged iter=',num2str(ip)]); 
            break
    end        
    end
    berr = norm(TT*y);       
end


function [berr,y] = minberrfromLanczos_dqds(T,y)
% input: T s.t. Lanczos AQ=QT, initial y for inv-power
% output backward error ||Ax-b||/(||A||*||x||)
% to compute soln x, 
% xsol = Q(:,1:numel(x))*x*(normb/(T(1,1:2)*x(1:2)));
    TT = T(2:end,1:end-1);     
    if nargin < 2, 
    y = randn(size(TT,1),1);
    end
    TTT = TT'; % keyboard
    for ip = 1:100 % inv power iteration
    ypre = y; 
    y = TTT\y; 
    y = TT\y; 
    y = y/norm(y);     
    if isnan(y)|isinf(y), keyboard, end
    if subspace(y,ypre)<1e-8, 
        % enable if desired:
            % disp(['invit converged iter=',num2str(ip)]); 
            break
    end        
    end
    berr = norm(TT*y);       
    [berr min(svd(TT))]
    keyboard
end
