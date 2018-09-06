%--------------------------------------------------------------------------
% Robust Adaptive Unscented Kalman Filter for Attitude  Estimation
% represented by unit quaternion
%
% Antonio C. B. Chiella, Bruno O. S. Teixeira, and Guilherme A. S. Pereira
% Federal University of Minas Gerais
% 
%--------------------------------------------------------------------------

%% =========================== Initial setup ==============================
clear all % Clear all variables
clc % Clear the workSpace

addpath('./Functions') % Add "Functions" folder


% Load the data
load data.mat

% Load the filtering configuration
filter_configuration

%% ============================ Filtering =================================

h = waitbar(0,'Wait!!!'); % Progress bar
for k=2:k_f
    
    % sampling time
    dt = t(k) - t(k-1);
    
    % compute the measurement
    [q_m(:,k), R_k(:,:,k)] = UT(@comp_quat,[a_m(:,k);B_m(:,k)],Ram);
    
    % ====================== Quaternion UKF ===============================
    % forecast
    [x_kk_1(:,k), Pxx_kk_1(:,:,k)] = forecast(@ffun, x_kk_1(:,k-1), ...
                                            w_m(:,k-1), Pxx_kk_1(:,:,k-1),...
                                            Q1, Q2, dt);
    % data-assimilation
    [x_kk_1(:,k), Pxx_kk_1(:,:,k), v_k_1(:,k), Pyy_kk1_1(:,:,k)] = ...
    data_assimilation(x_kk_1(:,k), Pxx_kk_1(:,:,k), @hfun, q_m(:,k), R_k(:,:,k));
     
    %============================ RAUKF ===================================
     % forecast
    [x_kk_2(:,k), Pxx_kk_2(:,:,k)] = forecast(@ffun, x_kk_2(:,k-1), ...
                                              w_m(:,k-1), Pxx_kk_2(:,:,k-1),...
                                              Q1, Q2, dt);
    % data-assimilation
    [x_kk_2(:,k), Pxx_kk_2(:,:,k), v_k_2(:,k), Pyy_kk1_2(:,:,k), R_adp(:,:,k)]...
     = data_assimilation_RAUKF(x_kk_2(:,k), Pxx_kk_2(:,:,k), ...
     @hfun, q_m(:,k), R_k(:,:,k), v_k_2, N);
    %======================================================================
    
    % Progress bar
    progres = round((k / k_f)*100);
    waitbar(k / k_f,h,sprintf('%d%% Running...',progres))
end
delete(h) 

%% ============================ Graphic ===================================

% Convert the estimated quaternion to Euler angles
for k = 1:k_f
    [x_kk_1_euler(:,k), Pxx_kk_1_euler(:,:,k), Pxy] = UT_q(@quat2euler,x_kk_1(1:4,k), Pxx_kk_1(1:3,1:3,k));
    [x_kk_2_euler(:,k), Pxx_kk_2_euler(:,:,k), Pxy] = UT_q(@quat2euler,x_kk_2(1:4,k), Pxx_kk_2(1:3,1:3,k));     
end
 %%
 figure(1)
 hold on
 xlabel('t(s)')
 ylabel('\phi [rad]')
 plot(t_r, deg2rad(phi_gt), 'k')
 plot(t, x_kk_1_euler(1,:), 'b')
 plot(t, x_kk_2_euler(1,:), 'r')
 legend('truth', 'UKF', 'RAUKF')
 
 figure(2)
 hold on
 xlabel('t(s)')
 ylabel('\theta [rad]')
 plot(t_r, deg2rad(theta_gt), 'k')
 plot(t, x_kk_1_euler(2,:), 'b')
 plot(t, x_kk_2_euler(2,:), 'r')
 legend('truth', 'UKF', 'RAUKF')
 
 figure(3)
 hold on
 xlabel('t(s)')
 ylabel('\psi [rad]')
 plot(t_r, deg2rad(psi_gt), 'k')
 plot(t, x_kk_1_euler(3,:), 'b')
 plot(t, x_kk_2_euler(3,:), 'r')
 legend('truth', 'UKF', 'RAUKF')
 