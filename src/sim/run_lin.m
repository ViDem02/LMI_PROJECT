function res = run_lin(lin, cfg, t_sim, dt)
% Integrate the linearised model in deviation coordinates.


A = lin.A; xeq = lin.xeq(:);  ueq = lin.ueq(:);

tspan = 0:dt:t_sim;
dx0   = cfg.design.x0(:);


[t, dx] = ode45(@(t,dx) A*dx, tspan, dx0, ...
                odeset('RelTol',1e-6,'AbsTol',1e-8));

res.t  = t;
res.dx = dx;
res.x  = xeq' + dx;
res.u = repmat(ueq', numel(t), 1);

end
