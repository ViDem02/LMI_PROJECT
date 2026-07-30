%% decompose_h2.m - which channel drives the linear/nonlinear h2 gap?
% Run this AFTER q1_lin.m so `res`, `xeq`, `ueq`, `cfg` are in scope,

startup_project;
cfg = load_params;
p   = cfg.plant;
[xeq, ueq] = find_equilibrium(cfg);

h1b = xeq(1); h2b = xeq(2);
phi1b = sqrt(2*p.g*h1b);   phi2b = sqrt(2*p.g*h2b);  
beta1 = p.g/phi1b;          beta2 = p.g/phi2b;       
u2 = ueq(2); u3 = ueq(3);

% reuse the large-deviation trajectories from q1_lin (res(2)) 
R    = res(2);
t    = R.nl.t;
xnl  = R.nl.x;         
xlin = R.lin.x;        

% helper: outflow velocity, true vs linearised (tangent at op. point) 
phi_true = @(h)  sqrt(2*p.g*max(h,0));
phi_lin  = @(h,hb,phib,betab)  phib + betab*(h - hb);   % 1st-order tangent

% tank-2 inflow (+ u2 a2 phi1) and outflow (- u3 a3 phi2), per model 
% nonlinear: phi evaluated with true sqrt on the nonlinear trajectory
in_nl   =  u2*p.a2 .* phi_true(xnl(:,1)) / p.A2;
out_nl  = -u3*p.a3 .* phi_true(xnl(:,2)) / p.A2;
h2dot_nl = in_nl + out_nl;

% linear: phi evaluated with the tangent on the linear trajectory
in_lin  =  u2*p.a2 .* phi_lin(xlin(:,1), h1b, phi1b, beta1) / p.A2;
out_lin = -u3*p.a3 .* phi_lin(xlin(:,2), h2b, phi2b, beta2) / p.A2;
h2dot_lin = in_lin + out_lin;

% the two competing ERROR channels (linear minus nonlinear) 
err_in  = in_lin  - in_nl;     % inflow-channel error  
err_out = out_lin - out_nl;    % outflow-channel error 
err_tot = h2dot_lin - h2dot_nl;

%% Figure: contributions over time
figDec = figure('Name','h2 flow decomposition'); 

 hold on; grid on
plot(t, err_in,  'LineWidth',1.4);
plot(t, err_out, 'LineWidth',1.4);
plot(t, err_tot, 'k','LineWidth',1.6);
yline(0,':','Color',[.5 .5 .5]);
xlabel('time [s]'); ylabel('linear $-$ nonlinear [m/s]','Interpreter','latex');
legend('inflow-channel error','outflow-channel error','net error', ...
       'Interpreter','latex','Location','best');
title('Channel analysis on error');


%% Numeric summary: integrated (net) contribution of each channel
fprintf('\nIntegrated over the run (trapz):\n');
fprintf('  inflow-channel error  : %+.3e  (linear over/under-fills via phi_1)\n', trapz(t,err_in));
fprintf('  outflow-channel error : %+.3e  (linear over/under-drains via phi_2)\n', trapz(t,err_out));
fprintf('  net h2dot error       : %+.3e  (sign => which model''s h2 is faster)\n', trapz(t,err_tot));

%% Save everything

style_fig(figDec, 'Size',[10 10])
snapshot_run(cfg, struct(), figDec, 'decompose_h2')