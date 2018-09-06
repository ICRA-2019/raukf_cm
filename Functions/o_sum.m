%--------------------------------------------------------------------
% Given a vector x1 = [x1_q x1_e]' and x2 = [x2_v x2_e]', where x1_q 
% is a unit quaternion vector, x2_v is a rotation vector, and x1_e 
% and x2_e are elements in Euclidian space, the function returns
% the sum s = (x1 + x2).
% NOTE: Quaternion follows Hamilton definition qw + iqx + jqy + kqz, where
% the vector representation is given by q = [qw qx qy qz]'.
%
% Antonio C. B. Chiella, UFMG, Jul 2017 
%--------------------------------------------------------------------
function s = o_sum(x1,x2)
    
    
    if length(x1)>4
        x1_q(1,:) = x1(1:4);
        x2_q(1,:) = r2q(x2(1:3)); % transform x2_v to x2_q
  
        s_q(:,1) = real(quatmultiply(x2_q,x1_q)); % rotation vector parameterization
        s_e(:,1) = x1(5:end) + x2(4:end);
        s = [s_q; s_e];
    end
    
    if length(x1) == 4 % only unity quaternion
        x1_q(1,:) = x1(1:4);
        x2_q(1,:) = r2q(x2(1:3)); % transform x2_v to x2_q
  
        s_q(:,1) = real(quatmultiply(x2_q,x1_q)); % rotation vector parameterization
        s = s_q;    
    end
    
    if length(x1)<4 % only elements in Euclidian space
         s_e(:,1) = x1(1:end) + x2(1:end);
         s = s_e;    
    end
   
end