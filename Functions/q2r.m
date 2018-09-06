%--------------------------------------------------------------------
% Mapping a unit quaternion q to its rotation vector parametrization
% q_rho. The quaternion is represented by q = [q0 q_x q_y q_z]^T.
%
% Antonio C. B. Chiella, UFMG, Jul 2017 
%--------------------------------------------------------------------
function q_rho = q2r(q)
%{
To resolve the ambiguity between q and -q and thus ensure that the shortest path between
both quaternions is considered, the argument of arccos is checked to be non-negative.
If it is negative, q1 and q2 are separated by more than 180 degrees along an axis q23/|q23|. Inverting
either q1 or q2, brings them closer together.
%}
   q = real(q);
   if q(1)>=0
        theta = 2*acos(q(1));
   else
       theta = -2*acos(-q(1));
   end
  
    
    if norm(q(2:4))>1e-100
        q_rho = theta*q(2:4)/norm(q(2:4));
    else
%         q_rho = theta*q(2:4);
        q_rho = 2*q(2:4);
    end
    
    q_rho = real(q_rho);

% [psi, theta, phi] = quat2angle(q, 'ZYX');
% q_rho = [phi, theta, psi];
end
