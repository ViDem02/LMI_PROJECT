function res = run_nl(uFun, cfg, t_sim, dt)
% Integrates the nonlinear system under a given ctrl law. 


p = cfg.plant;
eq = cfg.eq;

xeq = [eq.h1; eq.h2];
x0 = xeq(:) + cfg.design.dx0(:);
tspan = 0:dt:t_sim;

[t, x] = ode45(@(t,x) f_TankNLDynamics(t, x, uFun(t,x), p), tspan, x0, ...
               odeset('RelTol',1e-6,'AbsTol',1e-8));

N = numel(t);
m = numel(uFun(t(1), x(1,:)'));
u = zeros(N, m);
for i = 1:N
    u(i,:) = reshape(uFun(t(i), x(i,:)'), 1, []);
end

res.t = t;  res.x = x;  res.u = u; 
res.dx = x - xeq(:)';

end