%--------------------------------------------------------------------
% Given a set of vectors xi = [xi_q xi_e]' and weights wi where xi_q  
% is a unit  quaternion vector and xi_e is elements in Euclidian space, 
% the function  returns the weighted mean.
% NOTE: Quaternion follows Hamilton definition qw + iqx + jqy + kqz, where
% the vector representation is given by q = [qw qx qy qz]'.
%
% Antonio C. B. Chiella, UFMG, Jul 2017 
%--------------------------------------------------------------------
function wm = weighted_mean(xi, wi)

N = 5; % Optimization iterations
n_i = length(wi);

if size(xi,1)>4
    xi_q = xi(1:4,:);
    xi_e = xi(5:end,:);
    q_m = xi(1:4,1); % inicialização
    e_i = zeros(3,n_i);
    e_s = ones(3,1);
    % keyboard
    while norm(e_s)>10^-6
        for i = 1:n_i
    %         keyboard
            e_i(:,i) = wi(i)*q2r(quatmultiply(xi_q(:,i)',quatinv(q_m')));
        end   
        e_s = sum(e_i,2);
        q_m = quatmultiply(r2q(e_s)',q_m')';    
    end


    e_m = zeros(size(xi,1)-4,1);
    for i = 1:n_i
        e_m = e_m + wi(i)*xi_e(:,i);
    end
    
    wm = [q_m; e_m];
end
    
if size(xi,1)==4
    xi_q = xi(1:4,:);
    q_m = xi(1:4,1); % inicialização
    e_i = zeros(3,n_i);
    e_s = zeros(3,1);
    % keyboard
    while norm(e_s)>10^-6
        for i = 1:n_i
    %         keyboard
            e_i(:,i) = wi(i)*q2r(quatmultiply(xi_q(:,i)',quatinv(q_m')));
        end   
        e_s = sum(e_i,2);
        q_m = quatmultiply(r2q(e_s)',q_m')';    
    end
   
    wm = q_m;
    
end
    
if size(xi,1)<4
    xi_e = xi(1:end,:);
    e_m = zeros(size(xi,1),1);
    for i = 1:n_i
        e_m = e_m + wi(i)*xi_e(:,i);
    end
    wm = e_m;
end



    

        
end