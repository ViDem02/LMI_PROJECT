function res = run_observer(lin, K, L, xhat0, cfg, t_sim, dt)


A = lin.A;  B = lin.B;  C = lin.C;
xeq = lin.xeq(:);  ueq = lin.ueq(:);



p  = cfg.plant;

x0 = xeq + cfg.design.dx0(:);

tspan = 0:dt:t_sim;
[t, z] = ode45(@(t,z) obs_rhs(t, z, A, B, C, K, L, xeq, ueq, p), ...
               tspan, [x0(:); xhat0(:)], odeset('RelTol',1e-6,'AbsTol',1e-8));



% We simulate two coupled dynamical systems at once:
%   - the TRUE plant state      x    
%   - the observer's estimate   xhat 
%
% ode45 integrates a SINGLE state vector, so the two are stacked into one
% augmented vector of length 4, plant on top and estimate below:
%
%       z = [ x ; xhat ]  in R^4
%
% Packing is required because the two systems are mutually coupled at every
% instant and therefore cannot be integrated separately:
%
%   - the plant needs the estimate, through the 
%     control law u = ueq + K(xhat-xeq);
%
%   - the observer needs the true state, through 
%     the residual r = C(x - xhat).
%
% Neither trajectory is available ahead of the other, so a single 
% solver must advance on a common time step.
%
% Convention:
%   - initial condition : z0 = [x0 ; xhat0] in absolute coordinates
%   - unpack (below)     : x = z(1:2);  xhat = z(3:4);
%   - repack (below)     : dz = [dx ; dxhat];
% ----------------------------------------

x    = z(:, 1:2);
xhat = z(:, 3:4);

res.t    = t;
res.x    = x;
res.xhat = xhat;
res.u    = ueq' + (xhat - xeq') * K';
res.e    = x - xhat;
res.r    = (x - xeq') * C' - (xhat - xeq') * C';   % output residual y - yhat

end

function dz = obs_rhs(t, z, A, B, C, K, L, xeq, ueq, p)
x    = z(1:2);
xhat = z(3:4);

u   = ueq + K * (xhat - xeq);         % controller uses the estimate
dx  = f_TankNLDynamics(t, x, u, p);       % true (nonlinear) plant

dxh = xhat - xeq;
r   = C*(x - xeq) - C*dxh;            % residual = measured - predicted output
% this is innovation y - yhat, follows from C(x-xhat)

dxhat = A*dxh + B*(u - ueq) + L*r;    
%       predict             + correct

dz = [dx; dxhat];
end
