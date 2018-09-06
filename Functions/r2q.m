%--------------------------------------------------------------------
% Mapping rotation vector parametrization, q_rho, to its unit quaternion 
% representation. The quaternion is represented by q = [q0 q_x q_y q_z]^T.
%
% Antonio C. B. Chiella, UFMG, Jul 2017 
%--------------------------------------------------------------------
function q = r2q(q_rho)
    
    q_rho = real(q_rho);
    theta = norm(q_rho);
    
    
    if theta > 1e-100
        q(1,1) = cos(theta*0.5);
        q(2:4,1) = sin(theta*0.5)*q_rho/norm(q_rho);    
    else
        q(2:4,1) = sin(theta*0.5)*ones(3,1);
        q(1,1) = cos(theta*0.5);
%         q(1,1) = 1;
%         q(2:4,1) = q_rho/2;
    end
    
    q = real(q);
% q = euler2quat([q_rho(1);q_rho(2); q_rho(3)]);
        
end