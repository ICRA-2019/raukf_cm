function q = comp_quat(x)
a = x(1:3);
b = x(4:6);

%=====teste=====
% e_m = (norm(a)-9.81)/9.81;

%=======

a = a/norm(a);
b = b/norm(b);

if (a(3)>=0)
    q_a = [sqrt((a(3)+1)/2), -a(2)/sqrt(2*(a(3)+1)), a(1)/sqrt(2*(a(3)+1)), 0];
else
    q_a = [-a(2)/sqrt(2*(1-a(3))), sqrt((1-a(3))/2), 0, a(1)/sqrt(2*(1-a(3)))];
end



% alpha = -10*e_m + 2;
% if(alpha>1)
%     alpha = 1;
% end
% if(alpha<0)
%     alpha = 0;
% end
% q_i = [1 0 0 0];  
% 
% q_b = weighted_mean([q_i',q_a'], [1-alpha alpha])';

b = quat2rotm(q_a)'*b;

T = b(1)^2 + b(2)^2;
if (b(1)>=0)
    q_m = [sqrt(T+b(1)*sqrt(T))/sqrt(2*T), 0, 0, b(2)/(sqrt(2)*sqrt(T+b(1)*sqrt(T)))];
else
    q_m = [b(2)/(sqrt(2)*sqrt(T-b(1)*sqrt(T))), 0, 0, sqrt(T-b(1)*sqrt(T))/sqrt(2*T)];
end
% q= q_m';

 q = quatconj(quatmultiply(q_a, q_m))';
%q = quatmultiply(q_a, q_m)';
end

