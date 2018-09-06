%--------------------------------------------------------------------
% Antonio C. B. Chiella, UFMG, Jul 2017 
%
% Compute the forecast (prediction) step of curbature-UKF 
%
%--------------------------------------------------------------------
%{
    Output: 
           x_kk1 = mean
           Pxx_kk1 = covariance matrix
    Entrada:
           ffun = process function x_k = f(x_k-1,w_k-1, u_k-1,k-1)  
           x_k1k1 = prior state estimation
           u_k1 = measured input
           P_k1k1 = prior covariance matrix 
           Qk1 = multiplicative process noise covariance
           Qk2 = additive process noise covariance
           dt = sampled time

OBS: It is assumed that the state is partitioned as x = [x_q x_e]'.
%}
function [x_kk1, Pxx_kk1_v] = forecast(ffun, x_k1k1, u_k1, Pxx_k1k1_v, Qk1, Qk2, dt)
    
% Augmented state vector
n_q = size(Qk1,2);
x_aug = zeros(n_q,1);
x_k1k1_aug = [x_k1k1; x_aug];
    
% Augmented Covariance matrix
Pxx_k1k1_v_aug = blkdiag(Pxx_k1k1_v, Qk1);
    
% Unscented transform
[x_kk1, Pxx_kk1_v] = UT_q(ffun, x_k1k1_aug, Pxx_k1k1_v_aug, u_k1, dt);
    
Pxx_kk1_v = Pxx_kk1_v + Qk2; %covariance
   
end


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

function [y,Py] = UT_q(fun,x,Px, u, dt)

    n = size(Px,1);  % state vector dimension
    nsp = 2*n;          
    W = (1/nsp)*ones(1, nsp); % weights
   
    ns = size(fun(x(:,1), u, dt),1);
    Py = zeros(ns-1,ns-1,size(x,2));
    
    for i = 1:size(x,2)

        % Setermination of sigma ponts (SP)
        P_root = sqrt(n).*chol(Px(:,:,i))';
        P_root = [P_root -P_root];
        x_m =  repmat(x(:,i),1,nsp);
        
        % Propagation of SP  
        for k=1:nsp
         X(:,k) = o_sum(x_m(:,k), P_root(:,k));  
         Y(:,k) = fun(X(:,k),u,dt); 
        end
 
        % mean
        y(:,i) = weighted_mean(Y,W);
        
        % covariance
        for k=1:nsp
            error = o_minus(Y(:,k), y(:,i));
            Py(:,:,i) = Py(:,:,i) + W(k)*error*error';     
        end
    end
    
   
end

