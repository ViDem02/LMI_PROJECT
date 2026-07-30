function [A, B, C, D, E, F] = f_TankStateSpace(xeq, ueq, p)
% Linearized state-space matrices

h1 = xeq(1);  h2 = xeq(2);
u2 = ueq(2);  u3 = ueq(3);

phi1 = sqrt(2*p.g*h1);       phi2 = sqrt(2*p.g*h2);
beta1 = p.g / phi1;          beta2 = p.g / phi2;   % d(phi)/dh at the eq. point

% -- A = df/dx  (lower-triangular: tank 1 does not depend on h2)
A = [ -(u2*p.a2)*beta1/p.A1 ,            0                ;
       (u2*p.a2)*beta1/p.A2 ,  -(u3*p.a3)*beta2/p.A2      ];

% -- B = df/du
B = [  p.q0/p.A1 ,  -p.a2*phi1/p.A1 ,        0            ;
          0      ,   p.a2*phi1/p.A2 ,  -p.a3*phi2/p.A2     ];

% -- E = df/dw ,  w = [alpha1; alpha2]
E = [ -phi1/p.A1 ,      0       ;
          0      ,  -phi2/p.A2   ];

% -- Performance output z = [h1; h2] = x
C = eye(2);
D = zeros(2, 3);
F = zeros(2, 2);

end