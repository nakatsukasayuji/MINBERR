function [xsol,y] = minberrfromQ_faster(A,T,Q,b,yini)
% finds solution with minimum backward error in range(Q) for Ax=b. 
% Assumes Lanczos A*Q = Q*T + (final column)
% faster version, uses that T is tridiag. This is for *nonsymmetric* Ax=b.
% Input: Lanczos decomp A*Q=Q*T
    Tori = T; 
    if isa(A, 'function_handle'), bA = A(b)'; else bA = b'*A; end
    
    btil = bA*Q;
    Tori(1,1) = Tori(1,1) - norm(btil)^2/(norm(b)^2); % can probably be made more efficient. 

    Tsparse = sparse(Tori);        

    k = size(Tori,1);
    conv = 0;
    if nargin < 5, y = randn(k,1); else y = yini; end
    for ii = 1:100 % inverse iteration
        xpre = y;
        y = Tsparse\y;
        y = y/norm(y);
        if subspace(y,xpre)<1e-8,
            disp(['converged',num2str(ii)]), conv = 1;
            break, end
        %disp([subspace(x,xpre) ii])
    end

    if conv == 0  % do RQI if it hasn't converged
        iinow = ii;
        I = speye(k);
        Torinrm = max(diag(Tori));
        for ii = iinow:iinow+10 % RQ iteration
            lamnow = y'*Tsparse*y - 1e-10*Torinrm;
            xpre = y;
            y = (Tsparse-lamnow*I)\y;
            y = y/norm(y);
            %disp([subspace(x,xpre) ii])
        end
    end
    
    alp = norm(b)^2/(bA*(Q*y));
    xsol = Q*(y*alp);

    return
end
