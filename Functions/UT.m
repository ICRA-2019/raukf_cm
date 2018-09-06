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
%}

function [y, Py] = UT(fun,x,Px)

   n = size(Px,1);  % state vector dimension
   nsp = 2*n;          
   W = (1/nsp)*ones(1, nsp); % weights
    
   ns = size(fun(x(:,1)),1);
    
   Py = zeros(ns-1,ns-1,size(x,2));
   

   % Setermination of sigma ponts (SP)
   P_root = sqrt(n).*chol(Px)';
   X =  [P_root -P_root] + repmat(x,1,nsp);
        
   % Propagation of SP  
   for k=1:nsp  
    Y(:,k) = fun(X(:,k)); 
   end
   % mean
   y = weighted_mean(Y,W);
   
   % covariance
   for k=1:nsp
    error = o_minus(Y(:,k), y);
    Py = Py + W(k)*error*error';     
   end
   
end
    
