%Unscented Tranformation
%{
    Input:  
            - ffun: function
            - x: state vector, may have dimension nx1 or nxN
            - Px: covariance error matrix, may have dimension nxnx1 or nxnxN
    Output:
            - y: propagated mean
            - Py: propagated covariance
            - X: Sigma points
            - Y: propagated sigma points

   OBS: It is assumed that the state is partitioned as x = [x_q x_e]' and
        the function input is the state vector x, input vector u, and step 
        time dt.
   
%}

function [y, Py, Pxy] = UT_q(fun,x,Px)

 n = size(Px,1);  % state vector dimension
    nsp = 2*n;          
    W = (1/nsp)*ones(1, nsp); % weights
    
    ns = size(fun(x(:,1)),1);
    
    % Setermination of sigma ponts (SP)
    P_root = sqrt(n).*chol(Px(:,:,1))';
    P_root = [P_root -P_root];
    x_m =  + repmat(x(:,1),1,nsp);
        
    % Propagation of SP  
        
    for k=1:nsp
       X(:,k) = o_sum(x_m(:,k), P_root(:,k));   
       Y(:,k) = fun(X(:,k)); 
    end
    % mean
    y(:,1) = weighted_mean(Y,W);
        
    if size(Y(:,1),1)>=4
        % Covariance
        Py = zeros(ns-1,ns-1);
        % Cross-Covariance
        Pxy = zeros(size(Px,1), size(Y(:,1),1)-1);
    else
        % Covariance
        Py = zeros(ns,ns);
        % Cross-Covariance
        Pxy = zeros(size(Px,1), size(Y(:,1),1));
    end
        
    for k=1:nsp
        error_y = o_minus(Y(:,k), y(:,1));
%         keyboard
        Py(:,:,1) = Py(:,:,1) + W(k)*error_y*error_y';
        error_x = o_minus(X(:,k), x(:,1));
        Pxy = Pxy + W(k)*error_x*error_y';
     end
        
    
end