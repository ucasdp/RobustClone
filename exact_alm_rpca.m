function [A1 E1 iter] = exact_alm_rpca(D, lambda, tol, maxIter)


[m n] = size(D);

if nargin < 2
    lambda = 1 / sqrt(max(m,n));
end

if nargin < 3
    tol = 1e-7;
elseif tol == -1
    tol = 1e-7;
end


if nargin < 4
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end


% initialize
Y = sign(D);
norm_two = svds(Y, 1);
%norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;


A1 = zeros( m, n);
E1 = zeros( m, n);
normD = norm(D, 'fro');
tolProj = 1e-6 * normD;
total_svd = 0;
%mu = 10/norm_two; 
mu =0.5/norm_two; 
r = 6;          

iter = 0;
converged = false;
stopCon = 1;
sv = 5;
lenS = sv;
while ~converged       
    
    iter = iter + 1;
    
    % solve the primal problem by alternative projection
    primal_converged = false;
    primal_iter = 0;
    sv = sv + round(n * 0.1);
    while primal_converged == false
        
        h = D - A1 + (1/mu)*Y;
        E2 = max( h - lambda/mu,0) + min( h + lambda/mu,0); 
        
        if choosvd(n, sv) == 1
            [U S V] = svds(D - E2 + (1/mu)*Y, sv);
            %[U S V] = lansvd(D - E2 + (1/mu)*Y, sv, 'L');
        else
            [U S V] = svd(D - E2 + (1/mu)*Y, 'econ');
        end
        diagS = diag(S);
        lenS = length(find(diagS > 1/mu));
        if lenS < sv
            sv = min(lenS + 1, n);
        else
            sv = min(lenS + round(0.05*n), n);
        end
        A2 = U(:,1:lenS)*diag(diagS(1:lenS)-1/mu)*V(:,1:lenS)';    
        
        if norm(A1 - A2, 'fro') < tolProj && norm(E1 - E2, 'fro') < tolProj
            primal_converged = true;
        end
        A1 = A2;
        E1 = E2;
        primal_iter = primal_iter + 1;
        total_svd = total_svd + 1;
               
    end
        
    Z = D - A1 - E1;        
    Y = Y + mu*Z;
    mu = r * mu;
    
    
    %% stop condition    
    stopCon = norm(Z, 'fro') / normD;
    if stopCon < tol
        converged = true;
    end    
    
    if ~converged && iter >= maxIter
        disp('Reach maximum iterations ') ;
        converged = 1 ;       
    end
end

if nargin == 5
    fclose(fid);
end

