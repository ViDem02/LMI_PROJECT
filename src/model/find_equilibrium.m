function [xeq, ueq] = find_equilibrium(cfg)
% Steady-state inputs holding (h1,h2) fixed: 
% solve h1dot = h2dot = 0 for u1,u3.

p = cfg.plant;
h1 = cfg.eq.h1;
h2 = cfg.eq.h2;
u2 = cfg.eq.u2;

phi1 = sqrt(2*p.g*h1);
phi2 = sqrt(2*p.g*h2);

u1 = u2 * p.a2 * phi1 / p.q0;           % from h1dot = 0
u3 = u2 * p.a2 * phi1 / (p.a3 * phi2);  % from h2dot = 0




xeq = [h1; h2];
ueq = [u1; u2; u3];

end