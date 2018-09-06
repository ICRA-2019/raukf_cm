%--------------------------------------------------------------------
% Antonio C. B. Chiella, UFMG, Jul 2017 
%
% Compute the data-assimilation (update) step of curbature-UKF 
%
%--------------------------------------------------------------------
%{
    Output: 
           x_kk = mean
           Pxx_kk = covariance matrix
           v_k = innovation
           Pyy_kk1 = covariace matrix of innovation
    Entrada:
           hfun = obsservation function y_k = h(x_k,u_k) 
           x_kk1 = prior state estimation
           yk = measurement
           Pxx_kk1 = prior covariance matrix 
           Rk = measurement noise covariance
           
%}
function [x_kk, Pxx_kk_v, v_k, Pyy_kk1] = data_assimilation(x_kk1, Pxx_kk1_v, hfun, y_k, Rk)

[y_kk1, Pyy_kk1, Pxy_kk1] = UT_q(hfun, x_kk1, Pxx_kk1_v);

Pyy_kk1 = Pyy_kk1 + Rk;


% Data-assimilation step

v_k = o_minus(y_k,y_kk1); % innovation
K_k = Pxy_kk1/(Pyy_kk1); % Kalman gain
x_kk = o_sum(x_kk1, K_k*v_k);
Pxx_kk_v = Pxx_kk1_v - K_k*Pyy_kk1*K_k';

end

%{
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
        
    
    % Covariance
    Py = zeros(ns-1,ns-1,size(x,2));
    % Cross-Covariance
    Pxy = zeros(size(X,1)-1, size(Y(:,1),1)-1);
    for k=1:nsp
        error_y = o_minus(Y(:,k), y(:,1));
        Py(:,:,1) = Py(:,:,1) + W(k)*error_y*error_y';
        error_x = o_minus(X(:,k), x(:,1));
        Pxy = Pxy + W(k)*error_x*error_y';
     end
        
    
end
%}