function res = run_observer_lin(lin, K, L, xhat0, cfg, t_sim, dt)
% Integrate the coupled true+observer dynamics on the linear plant.
% Please see file run_observer for more insightful comments.

A = lin.A;  B = lin.B;  C = lin.C;
xeq = lin.xeq(:);  ueq = lin.ueq(:);

x0 = xeq + cfg.design.dx0(:);

tspan = 0:dt:t_sim;
[t, z] = ode45(@(t,z) obs_rhs_lin(z, A, B, C, K, L, xeq, ueq), ...
               tspan, [x0(:); xhat0(:)], odeset('RelTol',1e-8,'AbsTol',1e-10));

x    = z(:, 1:2);
xhat = z(:, 3:4);

res.t    = t;
res.x    = x;
res.xhat = xhat;
res.u    = ueq' + (xhat - xeq') * K';
res.e    = x - xhat;
res.r    = (x - xeq') * C' - (xhat - xeq') * C';   % output residual y - yhat

end

function dz = obs_rhs_lin(z, A, B, C, K, L, xeq, ueq)
    x    = z(1:2);   xhat = z(3:4);
    u    = ueq + K*(xhat - xeq);         % controller uses the estimate
    dx   = A*(x - xeq) + B*(u - ueq);    % LINEAR plant
    dxh  = xhat - xeq;
    r    = C*(x - xeq) - C*dxh;          % innovation
    dxhat= A*dxh + B*(u - ueq) + L*r;
    dz   = [dx; dxhat];
end
