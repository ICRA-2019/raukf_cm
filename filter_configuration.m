%--------------------------------------------------------------------
% Initial configuration of UKF_AE 
%--------------------------------------------------------------------

% covariance matrix of multiplicative noise (noise of the gyrometer)
Q1 = diag([0.008^2, 0.0065^2, 0.0086^2]); 
          
% Covariance matrix of additive process noise   
% diag([quaternion_3x1(only for stability) gyro_bias_3x1])
Q2 = diag([0.00000000000000000001, 0.00000000000000000001, 0.00000000000000000001, 0.000000001, 0.000000001, 0.000000001]); 

% Covariance matrix of measurement diag([acc_3x1 mag_3x1])
Ram = diag([0.0361^2; 0.0455^2; 0.0330^2;0.0011^2; 0.00098^2; 0.00098^2]);
% Compute the first measurement
[q_m(:,1), R_k(:,:,1)] = UT(@comp_quat,[a_m(:,1);B_m(:,1)],Ram);

% Initial state
x_00 = [q_m(:,1); 0; 0; 0]; % qw,qx,qy,qz, b_wx, b_wy, b_wz, 

% Initial state covariance matrix (rotation vector 3x3 and bias 3x3)
Pxx_00_v = diag([0.01, 0.01, 0.01, 0.00001, 0.00001, 0.00001]);

% ========================== Quaternion UKF ===============================
x_kk_1 = x_00;
Pxx_kk_1 = Pxx_00_v;
Pyy_kk1_1 = R_k;
v_k_1 = zeros(3,1);

% ================================= RAUKF =================================
x_kk_2 = x_00;
Pxx_kk_2 = Pxx_00_v;
Pyy_kk1_2 = R_k;
v_k_2 = zeros(3,1);
R_adp = R_k;
N = 20;
sigma = 3;

