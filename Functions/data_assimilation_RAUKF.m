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
function [x_kk, Pxx_kk_v, v_k, Pyy_kk1, R_adp] = data_assimilation_RAUKF(x_kk1, Pxx_kk1_v, hfun, y_k, Rk, v_k_w, N_w)

[y_kk1, Pyy_kk1, Pxy_kk1] = UT_q(hfun, x_kk1, Pxx_kk1_v);

v_k = o_minus(y_k,y_kk1); % innovation
% Adaptive step
v_k_w = [v_k_w v_k];
k = size(v_k_w,2);
n_y = size(Rk,2);
R_adp = Rk;
W_i = eye(n_y,n_y);
if k>=N_w,
  
    m = median(v_k_w(:,k-N_w+1:k),2);
    m_abs = abs(v_k_w(:,k-N_w+1:k) - repmat(m,1,N_w));
    MAD_i = median(m_abs,2);
    Pyy = zeros(n_y,n_y);
    for i=k-N_w+1:k
        sd = 3*1.4826*MAD_i;
        w_ii(:,i) = min(ones(n_y,1),sd./m_abs(:,i-k+N_w));
        W_i = diag(w_ii(:,i));     
        Pyy = Pyy + (1/(N_w))*(W_i*v_k_w(:,i))*(W_i*v_k_w(:,i))';
    end

    R_adp = diag(diag(Pyy - Pyy_kk1));
    R_adp = max(R_adp, Rk);
        
end

% Pyy_kk1 = Pyy_kk1 + inv(W_i)*R_adp;
Pyy_kk1 = Pyy_kk1 + R_adp;



% Data-assimilation step
K_k = Pxy_kk1/(Pyy_kk1); % Kalman gain
x_kk = o_sum(x_kk1, K_k*W_i*v_k);
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