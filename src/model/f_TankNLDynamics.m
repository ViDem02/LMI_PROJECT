function dx = f_TankNLDynamics(t, x, u, p)

% Nonlinear dynamics, to be given to ODE45. 

h1 = x(1); 
h2 = x(2);

s1 = sqrt(2 * p.g * h1);
s2 = sqrt(2 * p.g * h2);

dx = zeros(2, 1);
dx(1) = ( u(1)*p.q0 - u(2)*p.a2*s1 - p.alpha1*s1 ) / p.A1;
dx(2) = ( u(2)*p.a2*s1 - u(3)*p.a3*s2 - p.alpha2*s2 ) / p.A2;

end