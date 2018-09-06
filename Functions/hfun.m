% Observation map
%{
    Input:  
            - x: state vector, may have dimension nx1 
    Output:
            - y: estimated measurement       
%}
function y = hfun(x)
 e_w = x(1); e_x = x(2); e_y = x(3); e_z = x(4);

 y = [e_w; e_x; e_y; e_z];

end