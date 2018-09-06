% Process function
%{
    Input:  
            - x: state vector at time k-1
            - u: input vector (gyro measurements)
            - dt: sampling time
    Output:
            - x: state vector at time k       
%}
function x_k = ffun(x, u, dt)
  
  % Discount the respective bias and noise from input measurement
  x_noise = x(8:10,1); 
  x_bias = x(5:7,1);  
  u = u - x_noise(1:3,1) - x_bias;

  p = u(1); q = u(2); r = u(3);
  
  Omega =[ 0, -p, -q, -r;
           p, 0, r, -q;
           q, -r, 0, p;
           r, q, -p, 0];
       
  s = (dt/2)*norm(u);
  if (s > 0)
     sin_s = sin(s)/s;
  else 
     sin_s = 1;
  end
  M = cos(s).*eye(4) + (0.5*dt*sin_s)*Omega;

  x_k(1:4,1) = M*x(1:4,1);
    
  %Bias propagation
  x_k(5) = x(5);
  x_k(6) = x(6); 
  x_k(7) = x(7);
 
end