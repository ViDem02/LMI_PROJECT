function [xeq, ueq] = find_equilibrium(cfg)

p = cfg.plant;
h1 = cfg.eq.h1;
h2 = cfg.eq.h2;
u2 = cfg.eq.u2;

phi1 = sqrt(2*p.g*h1);
phi2 = sqrt(2*p.g*h2);

u1 = u2 * p.a2 * phi1 / p.q0;
u3 = u2 * p.a2 * phi1 / (p.a3 * phi2);

xeq = [h1; h2];
ueq = [u1; u2; u3];

end