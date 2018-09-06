%--------------------------------------------------------------------
% Given a vector x1 = [x1_q x1_e]' and x2 = [x2_q x2_e]', where x1_q 
% and x2_q are unit quaternion vectors and x1_e x2_e are elements in
% Euclidian space, the function returns the subtraction s = (x1 - x2).
% NOTE: Quaternion follows Hamilton definition qw + iqx + jqy + kqz, where
% the vector representation is given by q = [qw qx qy qz]'.
%
% Antonio C. B. Chiella, UFMG, Jul 2017 
%--------------------------------------------------------------------
function s = o_minus(x1,x2)
    
    
    
    
    if length(x1)>4
        x1_q(1,:) = x1(1:4);
        x2_q(1,:) = x2(1:4);
        s_r(:,1) = q2r(real(quatmultiply(x1_q,quatinv(x2_q)))); % rotation vector parameterization
        s_e(:,1) = x1(5:end) - x2(5:end);
        s = [s_r; s_e];
    end
    
    if length(x1)==4
        x1_q(1,:) = x1(1:4);
        x2_q(1,:) = x2(1:4);
        s_r(:,1) = q2r(real(quatmultiply(x1_q,quatinv(x2_q)))); % rotation vector parameterization
        s = s_r;
    end
    
    if length(x1)<4
        s_e(:,1) = x1(1:end) - x2(1:end);
        s = s_e;  
    end
        
end