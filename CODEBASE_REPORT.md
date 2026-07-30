# Codebase Report: LMI_PROJECT

## Directory Tree

```
LMI_PROJECT/
├── .gitattributes
├── config/
│   ├── design.json
│   ├── equilibrium.json
│   ├── observer_gains.json
│   └── plant_params.json
├── decompose_h2.m
├── find_u2.m
├── q1_lin.m
├── q2_feedback.asv
├── q2_feedback.m
├── q3_observer.m
├── report/
│   └── after-figureplacement/
│       ├── figures/
│       │   ├── find_u2.png
│       │   ├── q1_cmp.png
│       │   ├── q1_decomp.png
│       │   ├── q1_err.png
│       │   ├── q2_control.png
│       │   ├── q2_state_long.png
│       │   ├── q2_state_short.png
│       │   ├── q3_fb.png
│       │   └── q3_track.png
│       ├── references.bib
│       ├── report_fp.aux
│       ├── report_fp.bbl
│       ├── report_fp.blg
│       ├── report_fp.fdb_latexmk
│       ├── report_fp.fls
│       ├── report_fp.log
│       ├── report_fp.out
│       ├── report_fp.pdf
│       ├── report_fp.synctex.gz
│       └── report_fp.tex
├── selected_results/
│   ├── fig_01.png
│   ├── fig_02.png
│   ├── find_u2.png
│   └── h2 flow decomposition.png
├── src/
│   ├── io/
│   │   ├── load_params.m
│   │   ├── snapshot_run.m
│   │   └── style_fig.m
│   ├── model/
│   │   ├── f_TankNLDynamics.m
│   │   ├── f_TankStateSpace.m
│   │   └── find_equilibrium.m
│   └── sim/
│       ├── make_input.m
│       ├── run_lin.m
│       ├── run_nl.m
│       ├── run_observer.m
│       └── run_observer_lin.m
└── startup_project.m
```

## File Contents

# .gitattributes

```
# Auto detect text files and perform LF normalization
* text=auto
```

# config/design.json

```json
{
    "dx0": [-0.05, -0.05],
    "fake_zero" : 1e-3,
    "verbose": 0,
    "alpha_b": 2,
    "alpha_c": 0.0115,
    "alpha_obs" : 2,
    "alpha_d": 0.006,
    "gain_letter": "c"
}
```

# config/equilibrium.json

```json
{
  "h1": 0.2,
  "h2": 0.2,
  "u2": 0.3200
}
```

# config/observer_gains.json

```json
{
  "a": [[-694.7117, -53.2899],
        [-139.7823, -192.2493],
        [27.9687, 5.0172]],
  "b": [[-3532.1, -34.5],
        [-840.1, -872.3],
        [16.6, 0.8]],
  "c": [[-2.7426, -1.4209],
        [2.5719, -0.4021],
        [3.0223, 3.2663]]
}
```

# config/plant_params.json

```json
{
  "a2": 5.1e-5,
  "A1": 0.04,
  "a3": 5.0e-5,
  "A2": 0.04,
  "q0": 4.9e-5,
  "g": 9.81,
  "alpha1": 0,
  "alpha2": 0
}
```

# decompose_h2.m

```matlab
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
```

# find_u2.m

```matlab
% Script used to find the best u2 to allow the best margin 
% of maneuver for u1 and u3.

clearvars;

startup_project();
cfg = load_params;


space = 0:0.01:0.485;

array = zeros( width(space),3 );


i = 1;
for u2 = space

% configuration of the system
p = cfg.plant;
h1 = 0.2; h2 = 0.2;
phi1 = sqrt(2*p.g*h1);  phi2 = sqrt(2*p.g*h2);
u1   = u2 * p.a2*phi1 / p.q0;         % from h1dot = 0
u3   = u2 * p.a2*phi1 / (p.a3*phi2);  % from h2dot = 0
ubar = [u1; u2; u3];

% margins eval - considering that ctrl u will be capped between 0 and 1
lowMargin = ubar;
upperMargin = 1 - ubar; 
margin = min(lowMargin, upperMargin); 
[worst, iworst] = min(margin);
feasible = all(ubar >= 0) && all(ubar <= 1);


fprintf('ueq = [%.4f %.4f %.4f]\n', ubar);
fprintf("margins: %10.4f  %10.4f  %10.4f \n\n", margin(1), margin(2), margin(3))

array(i,:) = margin';
i = i + 1;
% how much space will each control have
end


T = array2table(array, ...
    'VariableNames', {'u1_margin','u2_margin','u3_margin'});

figUs = figure(1);
plot(space, T.u1_margin, 'LineWidth', 1.5);
hold on;
plot(space, T.u2_margin, 'LineWidth', 1.5);
plot(space, T.u3_margin, 'LineWidth', 1.5);
hold off;

grid on;
xlabel('u_2');
ylabel('Margin');
legend('u_1','u_2','u_3','Location','best');
title('Control margins vs. u_2');

%% Save everything


style_fig(figUs, 'Size',[8 8], 'MaxXTicks',6, 'MaxYTicks',6)
snapshot_run(cfg, struct('table', T), figUs, 'find_u2');
```

# q1_lin.m

```matlab
%% Init

clearvars;
startup_project;
cfg = load_params;

%% Linearization

p = cfg.plant;
eq = cfg.eq;

[xeq, ueq] = find_equilibrium(cfg);
[A, B, C, D, E, F] = f_TankStateSpace(xeq, ueq, p);

disp('xeq ='), disp(xeq), disp('ueq ='), disp(ueq)

fprintf("\n\nMatrices\n")

disp('A ='), disp(A), disp('B ='), disp(B), disp('E ='), disp(E)
disp('C ='), disp(C), disp('D ='), disp(D), disp('F ='), disp(F)


lambda = eig(A);
fprintf("\n\nEWs analysis\n")
fprintf('\topen-loop eigenvalues : %.10e , %.10e\n', lambda(1), lambda(2));
% both of them have negative real parts. 

fprintf('\tslowest time constant : %.0f s\n', 1/min(abs(real(lambda))));

% This means that the system is stable but quite slow. 

%% Simulation


% create the input
uFun = @(t, x) ueq;

% model for run_lin
lin  = struct('A', A, 'B', B, 'E', E, 'xeq', xeq, 'ueq', ueq);   

% run non linear simulation
t_sim = 2500;
dt    = 0.5;


% small vs large dev coordinates
cases = { [-0.02; -0.02], '\delta x small = -0.02m';
          [-0.18; -0.18], '\delta x large = -0.18m' };

res = struct('nl', {}, 'lin', {}, 'x0', {}, 'label', {});
for i = 1:size(cases, 1)
    dx0 = cases{i, 1};
    sim_cfg = cfg;
    sim_cfg.design.dx0 = dx0;

    % same scenario on each model, same time grid
    rnl  = run_nl(uFun, sim_cfg, t_sim, dt);    % nonlinear
    rlin = run_lin(lin, sim_cfg, t_sim, dt);    % linear 

    res(i).nl    = rnl;
    res(i).lin   = rlin;
    res(i).x0    = xeq + dx0;
    res(i).label = cases{i, 2};
end



%% Comparison figure between Linear and Non Linear

cNL  = [0      0.4470 0.7410];   % nonlinear (blue)
cLIN = [0.8500 0.3250 0.0980];   % linear    (red)
figCmp = figure(1); clf(figCmp); figCmp.Name = 'Q1: linear vs nonlinear';
tiledlayout(1, 4, 'TileSpacing', 'compact');
for i = 1:numel(res)
    for k = 1:2   % tank index
        nexttile((i-1)*2 + k); hold on; grid on; grid minor
        plot(res(i).nl.t,  res(i).nl.x(:,k),  'Color', cNL,  'LineWidth', 1.3)
        plot(res(i).lin.t, res(i).lin.x(:,k), 'Color', cLIN, 'LineWidth', 2.3, ...
             'LineStyle', '--')
        yline(xeq(k), ':', 'Color', [.5 .5 .5], 'LineWidth', 2);
        if i>1
            xlabel('time [s]'); 
        end
        ylabel(sprintf('h_%d [m]', k));
        title(sprintf('h_%d vs %s', k, res(i).label));
        if i == 1 && k == 1
            legend('nonlinear', 'linear', 'equilibrium', 'Location', 'best');
        end
    end
end

%% Model mismatch - small on LEFT (blue), large on RIGHT (red)

figErr = figure(2); clf(figErr); figErr.Name = 'Q1: linearization error';
cL = [0 .447 .741];    % blue - small
cR = [.85 .325 .098];  % red  - large


errS = vecnorm(res(1).nl.x - res(1).lin.x, 2, 2);
errL = vecnorm(res(2).nl.x - res(2).lin.x, 2, 2);

yyaxis left
hS = plot(res(1).nl.t, errS, '-', 'LineWidth', 1.6, ...
          'DisplayName', res(1).label);
hS.Color = cL;
ax = gca; ax.YColor = cL;
ylabel('|x_{nl}-x_{lin}| [m] - small', 'Color', cL);
grid on; grid minor

yyaxis right
hL = plot(res(2).nl.t, errL, '-', 'LineWidth', 1.6, ...
          'DisplayName', res(2).label);
hL.Color = cR; hL.LineStyle = '-';
ax.YColor = cR;
ylabel('|x_{nl}-x_{lin}| [m] - large', 'Color', cR);

%ax.YAxis(2).Exponent = -2;        % force ×10^{-2} on the RIGHT axis

xlabel('time [s]');
title('Lin error');
legend('Location', 'best');

%% Save everything


style_fig(figCmp, 'Size',[30 8], 'MaxXTicks',5, 'MaxYTicks',5, 'FontSize',11);
style_fig(figErr, 'Size',[8 8], 'MaxXTicks',5, 'MaxYTicks',5);
snapshot_run(cfg, struct('lin', lin), [figCmp, figErr], 'q1');
```

# q2_feedback.asv

```
%% Init

clearvars;
startup_project;
cfg = load_params;



%% Start from linearized system

[xeq, ueq] = find_equilibrium(cfg);
[A, B]     = f_TankStateSpace(xeq, ueq, cfg.plant);

fprintf("\n\nControllability analysis\n");
fprintf('\tcontrollability rank : %d / %d\n', rank(ctrb(A,B)), size(A,1));
fprintf('\topen-loop EWs      : %.3e , %.3e\n', eig(A));
% negative - we can move forward.


%% Create an arbitrary stabilizing full state feedback


%   We solve LMI:  W > 0,  He(A W + B X) < 0,
%   returning K = X W^-1 so that A + B K is Hurwitz. (Theory eq. 24-25.)

n = size(A,1); m = size(B,2); fz = cfg.design.fake_zero;
opts = sdpsettings('solver', "mosek", 'verbose', cfg.design.verbose);

yalmip('clear');
W = sdpvar(n,n);
X = sdpvar(m,n);
constr = [ W >= fz*eye(n);
           (A*W + B*X) + (A*W + B*X)' <= -fz*eye(n) ];
optimize(constr, [], opts);

[primalfeas, dualfeas] = check(constr);
feas = all([primalfeas; dualfeas] > -fz);

fprintf('\n\nArbitrary stabiliser\n')
if feas
    fprintf('\tconstraints satisfied.\n');
else
    fprintf('\tconstraints NOT satisfied.\n');
end

Ka = value(X) / value(W);

%% Full state feedback with guaranteed convergence rate 


%   K = FSF_CONVRATE(A, B, alpha, cfg) solves  W >= I,
%   He(A W + B X) <= -2 alpha W, returning K = X W^-1. (Theory eq. 26.)


n = size(A,1); m = size(B,2); fz = cfg.design.fake_zero;
opts = sdpsettings('solver', "mosek", 'verbose', cfg.design.verbose);

alpha_b = cfg.design.alpha_b;

yalmip('clear');
W = sdpvar(n,n);
X = sdpvar(m,n);
constr = [ W >= eye(n);
           (A*W + B*X) + (A*W + B*X)' <= -2*alpha_b*W ];
optimize(constr, [], opts);

[primalfeas, dualfeas] = check(constr);
feas = all([primalfeas; dualfeas] > -fz);

fprintf('\n\nConv-rate stabiliser\n')
if feas
    fprintf('\tconstraints satisfied.\n');
else
    fprintf('\tconstraints NOT satisfied.\n');
end

Kb = value(X) / value(W);


%% Guaranteed convergence rate within control bounds


% The following solves theory eq. (27):
%
% min kappa
% s.t.  W > I
%       He(A W + B X) <= -2 alpha W
%       [ kappa I   X' ;
%         X      kappa I ] >= 0


n = size(A,1); m = size(B,2); fz = cfg.design.fake_zero;
opts = sdpsettings('solver', "mosek", 'verbose', cfg.design.verbose);
alpha_c = cfg.design.alpha_c;

yalmip('clear');
W   = sdpvar(n,n);
X   = sdpvar(m,n);
kap = sdpvar(1,1);

constr = [ W >= eye(n);
           (A*W + B*X) + (A*W + B*X)' <= -2*alpha_c*W;
           [kap*eye(n)  X';
            X           kap*eye(m)] >= fz*eye(n+m) ];
optimize(constr, kap, opts);


[primalfeas, dualfeas] = check(constr);
feas = all([primalfeas; dualfeas] > -fz);

fprintf('\n\nFeasible stabiliser\n')
if feas
    fprintf('\tconstraints satisfied.\n');
else
    fprintf('\tconstraints NOT satisfied.\n');
end
fprintf('\tgain-size bound: |K| < %.4g\n\n', value(kap));


Kc     = value(X) / value(W);



fprintf('\nClosed-loop EWs:\n');
fprintf('\t(a) %.3e , %.3e\n', eig(A+B*Ka));
fprintf('\t(b) %.3e , %.3e\n', eig(A+B*Kb));
fprintf('\t(c) %.3e , %.3e\n', eig(A+B*Kc));


fprintf('\nGains:\n');
disp("Ka = "), disp(Ka), disp("Kb = "), disp(Kb), disp("Kc = ") , disp(Kc)
disp("kappa = "), disp(kap)




%% Simulate all four cases on the nonlinear model
t_sim = 3000;  dt = 0.1;
viz_cfg = cfg;

r_ol = run_nl(make_input([], xeq, ueq), viz_cfg, t_sim, dt);
r_a  = run_nl(make_input(Ka, xeq, ueq), viz_cfg, t_sim, dt);
r_b  = run_nl(make_input(Kb, xeq, ueq), viz_cfg, t_sim, dt);
r_c  = run_nl(make_input(Kc, xeq, ueq), viz_cfg, t_sim, dt);

runs   = {r_ol, r_a, r_b, r_c};
labels = {'open loop', '(a) arbitrary', '(b) conv-rate', '(c) feasible'};
cset   = {[.5 .5 .5], [0 .447 .741], [.85 .325 .098], [.466 .674 .188]};  % ol/a/b/c
uCol   = {[0 .447 .741], [.85 .325 .098], [.929 .694 .125]};              % u1/u2/u3



%% State response comparison  (figure 3)

figSt = figure(3); clf(figSt); figSt.Name = 'Q2: state response';
for k = 1:2
    subplot(2,1,k); hold on; grid on; grid minor
    for j = 1:numel(runs)
        plot(runs{j}.t, runs{j}.x(:,k), 'Color', cset{j}, 'LineWidth', 1.4);
    end
    yline(xeq(k), ':', 'Color', [.3 .3 .3], 'HandleVisibility', 'off');
    xlabel('time [s]'); ylabel(sprintf('h_%d [m]', k));
    if k == 1
        legend(labels, 'Location', 'best'); title('Closed loop vs open loop');
    end
    xlim([0 3000])
end


%% Control signals in each loop  (figure 4)
f4 = {r_a, r_b, r_c};  f4lab = {'(a) arbitrary',  '(b) conv-rate', '(c) feasible'};
xspan = {[0 25], [0 11], [0 600]};        % per-subplot time window
figU = figure(4); clf(figU); figU.Name = 'Q2: control signals';
for s = 1:3
    subplot(3,1,s)
    yyaxis left;  plot(f4{s}.t, f4{s}.u(:,1),   '-', 'LineWidth', 1.3); ylabel('u_1 [-]')
    yyaxis right; plot(f4{s}.t, f4{s}.u(:,2:3), '-', 'LineWidth', 1.3); ylabel('u_2, u_3 [-]')
    %yline(0, '--', 'Color', [.6 .6 .6]); yline(1, '--', 'Color', [.6 .6 .6]);
    grid on; title(f4lab{s}); %legend('u_1', 'u_2', 'u_3', 'Location', 'best')

    xlim(xspan{s});                        % shared x-axis: set once per subplot
    if s == 3                              % custom y-ranges only on the third
        yyaxis left;  ylim([0.63 0.9]);   % u1 axis
        yyaxis right; %ylim([0.25 0.40]);   % u2, u3 axis
    end
end
xlabel('time [s]');




%% Comparative table
% ---- setup ----
Tsim   = 3000;                       % long enough for the slow (c, open-loop) modes
x0abs  = xeq + cfg.design.dx0;       
sp     = xeq;                        % regulation setpoint (absolute levels)
band   = 0.02;                       % 2% settling band, as a fraction of the step size
tgrid  = (0:0.1:Tsim).';             % uniform grid
odeopt = odeset('RelTol',1e-8,'AbsTol',1e-10);
p      = cfg.plant;
tmin  = 1.4;                          % transient floor [s]: ignore any settling
                                     % detected in the first tmin seconds, so the
                                     % initial overshoot of the aggressive gains
                                     % is not mistaken for true settling.

% ---- controllers to compare (label, gain) ----
ctrls = { 'open-loop', zeros(size(Ka)); ...
          'K_a',       Ka; ...
          'K_b',       Kb; ...
          'K_c',       Kc };

% ---- pre-allocate one array per table column ----
label        = strings(4,1);   % controller name
settlingTime = zeros(4,1);     % time to enter & stay in the 2% band [s]
inputMax     = zeros(4,1);     % max over all u components, whole run
inputMin     = zeros(4,1);     % min over all u components, whole run

for i = 1:4
    label(i) = ctrls{i,1};  K = ctrls{i,2};

    rhs   = @(t,x) f_TankNLDynamics(t, x, ueq + K*(x - xeq), p);
    [t,X] = ode45(rhs, tgrid, x0abs, odeopt);
    U     = (ueq + K*(X.' - xeq)).';

    tol = band .* abs(x0abs - sp);           % 2% band, per state

    ts_perTank = zeros(2,1);
    for k = 1:2
        outside = abs(X(:,k) - sp(k)) > tol(k);   % true when outside the band
        outside(t < tmin) = true;                 % <-- discard first tmin seconds
        last_out = find(outside, 1, 'last');
        ts_perTank(k) = t(last_out);              % last exit from the band
    end
    settlingTime(i) = max(ts_perTank);            % slower tank governs

    inputMax(i) = max(U(:));   inputMin(i) = min(U(:));
end

% ---- table ----
T = table(label, settlingTime, inputMax, inputMin, ...
    'VariableNames', {'controller','settlingTime_s','inputMax','inputMin'});


fprintf("\n\nComparison among gains\n\n")
disp(T)

% ---- LaTeX table body (.3cm row spacing) ----
fprintf('\n\n%%  LaTeX table body \n');
for i = 1:4
    fprintf('%s & %.0f & %.3f & %.3f \\\\[.3cm]\n', ...
        label(i), settlingTime(i), inputMax(i), inputMin(i));
end



%% Save everything

style_fig(figSt, "Size", [8 13], "FontSize",11)
style_fig(figU, "Size", [8 8])
snapshot_run(cfg, struct('Ka', Ka, 'Kb', Kb, 'Kc', Kc, 'T', T), [figSt, figU], 'q2');
```

# q2_feedback.m

```matlab
%% Init

clearvars;
startup_project;
cfg = load_params;



%% Start from linearized system

[xeq, ueq] = find_equilibrium(cfg);
[A, B]     = f_TankStateSpace(xeq, ueq, cfg.plant);

fprintf("\n\nControllability analysis\n");
fprintf('\tcontrollability rank : %d / %d\n', rank(ctrb(A,B)), size(A,1));
fprintf('\topen-loop EWs      : %.3e , %.3e\n', eig(A));
% negative - we can move forward.


%% Create an arbitrary stabilizing full state feedback


%   We solve LMI:  W > 0,  He(A W + B X) < 0,
%   returning K = X W^-1 so that A + B K is Hurwitz. (Theory eq. 24-25.)

n = size(A,1); m = size(B,2); fz = cfg.design.fake_zero;
opts = sdpsettings('solver', "mosek", 'verbose', cfg.design.verbose);

yalmip('clear');
W = sdpvar(n,n);
X = sdpvar(m,n);
constr = [ W >= fz*eye(n);
           (A*W + B*X) + (A*W + B*X)' <= -fz*eye(n) ];
optimize(constr, [], opts);

[primalfeas, dualfeas] = check(constr);
feas = all([primalfeas; dualfeas] > -fz);

fprintf('\n\nArbitrary stabiliser\n')
if feas
    fprintf('\tconstraints satisfied.\n');
else
    fprintf('\tconstraints NOT satisfied.\n');
end

Ka = value(X) / value(W);

%% Full state feedback with guaranteed convergence rate 


%   K = FSF_CONVRATE(A, B, alpha, cfg) solves  W >= I,
%   He(A W + B X) <= -2 alpha W, returning K = X W^-1. (Theory eq. 26.)


n = size(A,1); m = size(B,2); fz = cfg.design.fake_zero;
opts = sdpsettings('solver', "mosek", 'verbose', cfg.design.verbose);

alpha_b = cfg.design.alpha_b;

yalmip('clear');
W = sdpvar(n,n);
X = sdpvar(m,n);
constr = [ W >= eye(n);
           (A*W + B*X) + (A*W + B*X)' <= -2*alpha_b*W ];
optimize(constr, [], opts);

[primalfeas, dualfeas] = check(constr);
feas = all([primalfeas; dualfeas] > -fz);

fprintf('\n\nConv-rate stabiliser\n')
if feas
    fprintf('\tconstraints satisfied.\n');
else
    fprintf('\tconstraints NOT satisfied.\n');
end

Kb = value(X) / value(W);


%% Guaranteed convergence rate within control bounds


% The following solves theory eq. (27):
%
% min kappa
% s.t.  W > I
%       He(A W + B X) <= -2 alpha W
%       [ kappa I   X' ;
%         X      kappa I ] >= 0


n = size(A,1); m = size(B,2); fz = cfg.design.fake_zero;
opts = sdpsettings('solver', "mosek", 'verbose', cfg.design.verbose);
alpha_c = cfg.design.alpha_c;

yalmip('clear');
W   = sdpvar(n,n);
X   = sdpvar(m,n);
kap = sdpvar(1,1);

constr = [ W >= eye(n);
           (A*W + B*X) + (A*W + B*X)' <= -2*alpha_c*W;
           [kap*eye(n)  X';
            X           kap*eye(m)] >= fz*eye(n+m) ];
optimize(constr, kap, opts);


[primalfeas, dualfeas] = check(constr);
feas = all([primalfeas; dualfeas] > -fz);

fprintf('\n\nFeasible stabiliser\n')
if feas
    fprintf('\tconstraints satisfied.\n');
else
    fprintf('\tconstraints NOT satisfied.\n');
end
fprintf('\tgain-size bound: |K| < %.4g\n\n', value(kap));


Kc     = value(X) / value(W);



fprintf('\nClosed-loop EWs:\n');
fprintf('\t(a) %.3e , %.3e\n', eig(A+B*Ka));
fprintf('\t(b) %.3e , %.3e\n', eig(A+B*Kb));
fprintf('\t(c) %.3e , %.3e\n', eig(A+B*Kc));


fprintf('\nGains:\n');
disp("Ka = "), disp(Ka), disp("Kb = "), disp(Kb), disp("Kc = ") , disp(Kc)
disp("kappa = "), disp(kap)




%% Simulate all four cases on the nonlinear model
t_sim = 3000;  dt = 0.1;
viz_cfg = cfg;

r_ol = run_nl(make_input([], xeq, ueq), viz_cfg, t_sim, dt);
r_a  = run_nl(make_input(Ka, xeq, ueq), viz_cfg, t_sim, dt);
r_b  = run_nl(make_input(Kb, xeq, ueq), viz_cfg, t_sim, dt);
r_c  = run_nl(make_input(Kc, xeq, ueq), viz_cfg, t_sim, dt);

runs   = {r_ol, r_a, r_b, r_c};
labels = {'open loop', '(a) arbitrary', '(b) conv-rate', '(c) feasible'};
cset   = {[.5 .5 .5], [0 .447 .741], [.85 .325 .098], [.466 .674 .188]};  % ol/a/b/c
uCol   = {[0 .447 .741], [.85 .325 .098], [.929 .694 .125]};              % u1/u2/u3



%% State response comparison  (figure 3)

figSt = figure(3); clf(figSt); figSt.Name = 'Q2: state response';
for k = 1:2
    subplot(2,1,k); hold on; grid on; grid minor
    for j = 1:numel(runs)
        plot(runs{j}.t, runs{j}.x(:,k), 'Color', cset{j}, 'LineWidth', 1.4);
    end
    yline(xeq(k), ':', 'Color', [.3 .3 .3], 'HandleVisibility', 'off');
    xlabel('time [s]'); ylabel(sprintf('h_%d [m]', k));
    if k == 1
        legend(labels, 'Location', 'best'); title('Closed loop vs open loop');
    end
    xlim([0 15])
end


%% Control signals in each loop  (figure 4)
f4 = {r_a, r_b, r_c};  f4lab = {'(a) arbitrary',  '(b) conv-rate', '(c) feasible'};
xspan = {[0 25], [0 5], [0 600]};        % per-subplot time window
figU = figure(4); clf(figU); figU.Name = 'Q2: control signals';
for s = 1:3
    subplot(3,1,s)
    yyaxis left;  plot(f4{s}.t, f4{s}.u(:,1),   '-', 'LineWidth', 1.3); ylabel('u_1 [-]')
    yyaxis right; plot(f4{s}.t, f4{s}.u(:,2:3), '-', 'LineWidth', 1.3); ylabel('u_2, u_3 [-]')
    %yline(0, '--', 'Color', [.6 .6 .6]); yline(1, '--', 'Color', [.6 .6 .6]);
    grid on; title(f4lab{s}); %legend('u_1', 'u_2', 'u_3', 'Location', 'best')

    xlim(xspan{s});                        % shared x-axis: set once per subplot
    if s == 3                              % custom y-ranges only on the third
        yyaxis left;  ylim([0.63 0.9]);   % u1 axis
        yyaxis right; %ylim([0.25 0.40]);   % u2, u3 axis
    end
end
xlabel('time [s]');




%% Comparative table
% ---- setup ----
Tsim   = 3000;                       % long enough for the slow (c, open-loop) modes
x0abs  = xeq + cfg.design.dx0;       
sp     = xeq;                        % regulation setpoint (absolute levels)
band   = 0.02;                       % 2% settling band, as a fraction of the step size
tgrid  = (0:0.1:Tsim).';             % uniform grid
odeopt = odeset('RelTol',1e-8,'AbsTol',1e-10);
p      = cfg.plant;
tmin  = 1.4;                          % transient floor [s]: ignore any settling
                                     % detected in the first tmin seconds, so the
                                     % initial overshoot of the aggressive gains
                                     % is not mistaken for true settling.

% ---- controllers to compare (label, gain) ----
ctrls = { 'open-loop', zeros(size(Ka)); ...
          'K_a',       Ka; ...
          'K_b',       Kb; ...
          'K_c',       Kc };

% ---- pre-allocate one array per table column ----
label        = strings(4,1);   % controller name
settlingTime = zeros(4,1);     % time to enter & stay in the 2% band [s]
inputMax     = zeros(4,1);     % max over all u components, whole run
inputMin     = zeros(4,1);     % min over all u components, whole run

for i = 1:4
    label(i) = ctrls{i,1};  K = ctrls{i,2};

    rhs   = @(t,x) f_TankNLDynamics(t, x, ueq + K*(x - xeq), p);
    [t,X] = ode45(rhs, tgrid, x0abs, odeopt);
    U     = (ueq + K*(X.' - xeq)).';

    tol = band .* abs(x0abs - sp);           % 2% band, per state

    ts_perTank = zeros(2,1);
    for k = 1:2
        outside = abs(X(:,k) - sp(k)) > tol(k);   % true when outside the band
        outside(t < tmin) = true;                 % <-- discard first tmin seconds
        last_out = find(outside, 1, 'last');
        ts_perTank(k) = t(last_out);              % last exit from the band
    end
    settlingTime(i) = max(ts_perTank);            % slower tank governs

    inputMax(i) = max(U(:));   inputMin(i) = min(U(:));
end

% ---- table ----
T = table(label, settlingTime, inputMax, inputMin, ...
    'VariableNames', {'controller','settlingTime_s','inputMax','inputMin'});


fprintf("\n\nComparison among gains\n\n")
disp(T)

% ---- LaTeX table body (.3cm row spacing) ----
fprintf('\n\n%%  LaTeX table body \n');
for i = 1:4
    fprintf('%s & %.0f & %.3f & %.3f \\\\[.3cm]\n', ...
        label(i), settlingTime(i), inputMax(i), inputMin(i));
end



%% Save everything

style_fig(figSt, "Size", [8 13], "FontSize",11)
style_fig(figU, "Size", [8 8], "MaxXTicks", 6)
snapshot_run(cfg, struct('Ka', Ka, 'Kb', Kb, 'Kc', Kc, 'T', T), [figSt, figU], 'q2');
```

# q3_observer.m

```matlab
%% Init

clearvars;
startup_project;

cfg = load_params;
[xeq, ueq] = find_equilibrium(cfg);
[A, B]     = f_TankStateSpace(xeq, ueq, cfg.plant);

%% Select gain

gain_letter = cfg.design.gain_letter;
Kgain = cfg.observer_gains.(gain_letter);


%% Analysis of observability

C1 = [1 0];   % observe only h1
C2 = [0 1];   % observe only h2
n  = size(A,1);


r1 = rank(obsv(A, C1));
r2 = rank(obsv(A, C2));

fprintf("\nObservability analysis\n")
fprintf("\tRank of C1 is %d, therefore h2 cannot " + ...
    "be observed through h1\n", r1);
fprintf("\tRank of C2 is %d, therefore h1 can be observed " + ...
    "through h2\n", r2);

C = C2;

%% Finding the observer gain

% We find L so that A - L C is Hurwitz, while maintaining convergence 
% rate >= alpha, by solving the standard observer LMI:
%
%       P > 0,   He(P A - Y C) <= -2 alpha P,   L = P^-1 Y.
%
%   This is the dual of the state-feedback problem. Note C is the measured
%   output map (here [0 1], i.e. h2), not the performance output z.

n = size(A,1); q = size(C,1); fz = cfg.design.fake_zero;
opts = sdpsettings('solver', 'mosek', 'verbose', cfg.design.verbose);

alpha_obs = cfg.design.alpha_obs;

yalmip('clear');
P = sdpvar(n,n);
Y = sdpvar(n,q);
constr = [ P >= eye(n);
           (P*A - Y*C) + (P*A - Y*C)' <= -2*alpha_obs*P ];
optimize(constr, [], opts);

[primalfeas, dualfeas] = check(constr);
feas = all([primalfeas; dualfeas] > -fz);

disp('Observer gain')
if feas
    fprintf('     constraints satisfied.\n');
else
    fprintf('     constraints NOT satisfied.\n');
end

L = value(P) \ value(Y);
disp('L ='), disp(L)      

%% EWs analysis
fprintf('\n\nEW analysis\n')
fprintf('\n\ncontroller EWs (A+B*K) : %.3e , %.3e\n', eig(A + B*Kgain));
fprintf('observer   EWs (A-L*C) : %.3e , %.3e\n', eig(A - L*C));

%% Simulate the observer on the nonlinear model
t_sim = 400;  dt = 0.1;    

r_fs = run_nl(make_input(Kgain, xeq, ueq), cfg, t_sim, dt);

x0  = xeq + cfg.design.dx0(:);
xhat0 = xeq;
xhat0(2) = x0(2);
lin  = struct('A', A, 'B', B, 'C', C, 'xeq', xeq, 'ueq', ueq);   % model bundle
r_ob = run_observer(lin, Kgain, L, xhat0, cfg, t_sim, dt);


%% Observer estimate vs actual state
cTrue = [0      0.4470 0.7410];   % actual state (blue)
cEst  = [0.8500 0.3250 0.0980];   % estimate     (red)
cObs  = [0.4660 0.6740 0.1880];   % observer-based loop (green)
figEst = figure(5); clf(figEst); figEst.Name = 'Q3: observer estimates';
for k = 1:2
    subplot(2,1,k); hold on; grid on; grid minor
    plot(r_ob.t, r_ob.x(:,k),    'Color', cTrue, 'LineWidth', 1.5);
    plot(r_ob.t, r_ob.xhat(:,k), 'Color', cEst,  'LineWidth', 1.3, 'LineStyle', '--');
    yline(xeq(k), ':', 'Color', [.3 .3 .3], 'HandleVisibility', 'off');
    ylabel(sprintf('h_%d [m]', k));
    if k == 1
        legend('actual x', 'estimate xhat', 'Location', 'best');
        title('Observer tracking (only h_2 is measured)');
    end
    xlim([0 400])
end
xlabel('time [s]');


%% Full-state vs observer-based closed loop
figCmp = figure(6); clf(figCmp); figCmp.Name = 'Q3: full-state vs observer';
for k = 1:2
    subplot(2,1,k); hold on; grid on; grid minor
    plot(r_fs.t, r_fs.x(:,k), 'Color', cTrue, 'LineWidth', 1.5);
    plot(r_ob.t, r_ob.x(:,k), 'Color', cObs,  'LineWidth', 1.3, 'LineStyle', '--');
    yline(xeq(k), ':', 'Color', [.3 .3 .3], 'HandleVisibility', 'off');
    ylabel(sprintf('h_%d [m]', k));
    if k == 1
        legend('full-state feedback', 'observer-based', 'Location', 'best');
        title('Close loop behaviour');
    end
    xlim([0 400])
end
xlabel('time [s]');


%% Performance table

t_tab = 3000;  dt = 0.1;              % long horizon so the slow Kc loop settles

r_fs_t = run_nl(make_input(Kgain, xeq, ueq), cfg, t_tab, dt);   % full-state
x0 = xeq + cfg.design.dx0(:);  xhat0 = xeq;  xhat0(2) = x0(2);
lin = struct('A',A,'B',B,'C',C,'xeq',xeq,'ueq',ueq);
r_ob_t = run_observer(lin, Kgain, L, xhat0, cfg, t_tab, dt);    % observer-based

sp = xeq;  band = 0.02;  tmin = 1.4;  step = abs(cfg.design.dx0(:));
loops = {'full-state', r_fs_t; 'observer-based', r_ob_t};
label = strings(2,1); ts1 = zeros(2,1); ts2 = zeros(2,1);
uMax = zeros(2,1); uMin = zeros(2,1);
for i = 1:2
    label(i) = loops{i,1};  R = loops{i,2};
    t = R.t; X = R.x; U = R.u;  tol = band.*step;  ts = zeros(2,1);
    for k = 1:2
        outside = abs(X(:,k)-sp(k)) > tol(k);
        outside(t < tmin) = true;
        ts(k) = t(find(outside,1,'last'));
    end
    ts1(i)=ts(1); ts2(i)=ts(2);  uMax(i)=max(U(:)); uMin(i)=min(U(:));
end
T3 = table(label, ts1, ts2, uMax, uMin, ...
    'VariableNames',{'loop','ts_h1_s','ts_h2_s','uMax','uMin'});

fprintf('\n\nPerformance table\n')
disp(T3)
fprintf('\n%% LaTeX body\n');
for i=1:2
    fprintf('%s & %.0f & %.0f & %.3f & %.3f \\\\[.3cm]\n', ...
        label(i), ts1(i), ts2(i), uMax(i), uMin(i));
end

% estimation convergence on the nonlinear plant: when is |h1 - h1hat| < 1 mm and stays?
e1  = abs(r_ob_t.x(:,1) - r_ob_t.xhat(:,1));
thr = 1e-3;
idx = find(e1 >= thr, 1, 'last');
fprintf('\n\nMoment of convergence on non-linear plant\n')
fprintf('\th1 estimate within %.0e m after t = %.1f s\n', thr, r_ob_t.t(idx));


%% Observer sanity check on the linear plant to see correction dynamics

r_ob_lin = run_observer_lin(lin, Kgain, L, xhat0, cfg, t_sim, dt);

e   = vecnorm(r_ob_lin.x - r_ob_lin.xhat, 2, 2);
idx = find(e >= 1e-3, 1, 'last');
fprintf('\n\nMoment of convergence on linear plant\n')
if isempty(idx)
    fprintf('\tEstimate within 1e-03 for the whole run.\n');
else
    fprintf('\tEstimation error stays below 1e-4 after t = %.2f s\n', ...
        r_ob_lin.t(min(idx+1,numel(r_ob_lin.t))));
end

%% Design K_d 
m = size(B,2);
alpha_d =  0.0095;  %cfg.design.alpha_d;

yalmip('clear');
W = sdpvar(n,n);  X = sdpvar(m,n);  kap = sdpvar(1,1);
constr = [ W >= eye(n);
           (A*W + B*X) + (A*W + B*X)' <= -2*alpha_d*W;
           [kap*eye(n)  X';
            X           kap*eye(m)] >= fz*eye(n+m) ];
optimize(constr, kap, opts);
Kd = value(X) / value(W);

fprintf('\n\nK_d design (alpha_d = %.4g)\n', alpha_d);
fprintf('\t|K| < %.4g,  EWs %.3e , %.3e\n', value(kap), eig(A+B*Kd));
disp("Kd = ")
disp(Kd)

%% extended statistics

r_ob_d = run_observer(lin, Kd, L, xhat0, cfg, t_tab, dt);   % observer-based, K_d

loops = {'full-state Kc',      r_fs_t; ...
         'observer-based Kc',  r_ob_t; ...
         'observer-based Kd',  r_ob_d};
nL = size(loops,1);
label = strings(nL,1); ts1 = zeros(nL,1); ts2 = zeros(nL,1);
uMax = zeros(nL,1); uMin = zeros(nL,1);
for i = 1:nL
    label(i) = loops{i,1};  R = loops{i,2};
    t = R.t; X = R.x; U = R.u;  tol = band.*step;  ts = zeros(2,1);
    for k = 1:2
        outside = abs(X(:,k)-sp(k)) > tol(k);
        outside(t < tmin) = true;
        ts(k) = t(find(outside,1,'last'));
    end
    ts1(i)=ts(1); ts2(i)=ts(2);  uMax(i)=max(U(:)); uMin(i)=min(U(:));
end
T3 = table(label, ts1, ts2, uMax, uMin, ...
    'VariableNames',{'loop','ts_h1_s','ts_h2_s','uMax','uMin'});
fprintf('\n\nPerformance table\n'); disp(T3)


%% Save everything

style_fig(figCmp, "Size",[10 12],"FontSize",10)
style_fig(figEst, "Size",[10 12],"FontSize",10)
snapshot_run(cfg, struct('L', L, 'Kgain', Kgain, 'T3', T3), [figCmp, figEst], 'q3')
```

# report/after-figureplacement/figures/find_u2.png

*(binary file, content omitted)*

# report/after-figureplacement/figures/q1_cmp.png

*(binary file, content omitted)*

# report/after-figureplacement/figures/q1_decomp.png

*(binary file, content omitted)*

# report/after-figureplacement/figures/q1_err.png

*(binary file, content omitted)*

# report/after-figureplacement/figures/q2_control.png

*(binary file, content omitted)*

# report/after-figureplacement/figures/q2_state_long.png

*(binary file, content omitted)*

# report/after-figureplacement/figures/q2_state_short.png

*(binary file, content omitted)*

# report/after-figureplacement/figures/q3_fb.png

*(binary file, content omitted)*

# report/after-figureplacement/figures/q3_track.png

*(binary file, content omitted)*

# report/after-figureplacement/references.bib

```
@book{hespanha,
  author    = {Hespanha, Jo{\~a}o P.},
  title     = {Linear Systems Theory},
  publisher = {Princeton University Press},
  address   = {Princeton, NJ},
  year      = {2009},
  edition   = {1st}
}

@misc{zaccarian_lmi,
  author       = {Zaccarian, Luca and Sferlazza, Alessandro},
  title        = {Analysis and Synthesis of Feedback Systems:
                  Quadratic Functions and {LMIs}},
  howpublished = {Lecture notes},
  note         = {Automatic Control course, University of Trento},
}
```

# report/after-figureplacement/report_fp.aux

```
\relax 
\providecommand*\new@tpo@label[2]{}
\providecommand \babel@aux [2]{\global \let \babel@toc \@gobbletwo }
\@nameuse{bbl@beforestart}
\providecommand\hyper@newdestlabel[2]{}
\providecommand\HyField@AuxAddToFields[1]{}
\providecommand\HyField@AuxAddToCoFields[2]{}
\providecommand\BKM@entry[2]{}
\BKM@entry{id=1,dest={73656374696F6E2E31},srcline={182}}{5C3337365C3337375C303030495C3030306E5C303030745C303030725C3030306F5C303030645C303030755C303030635C303030745C303030695C3030306F5C3030306E}
\BKM@entry{id=2,dest={73656374696F6E2E32},srcline={192}}{5C3337365C3337375C3030304D5C3030306F5C303030645C303030655C3030306C5C3030306C5C303030695C3030306E5C303030675C3030305C3034305C303030615C3030306E5C303030645C3030305C3034305C3030304C5C303030695C3030306E5C303030655C303030615C303030725C303030695C303030735C303030615C303030745C303030695C3030306F5C3030306E}
\BKM@entry{id=3,dest={73756273656374696F6E2E322E31},srcline={195}}{5C3337365C3337375C3030304F5C303030705C303030655C303030725C303030615C303030745C303030695C3030306E5C303030675C3030305C3034305C303030705C3030306F5C303030695C3030306E5C30303074}
\BKM@entry{id=4,dest={73756273656374696F6E2E322E32},srcline={229}}{5C3337365C3337375C3030304C5C303030695C3030306E5C303030655C303030615C303030725C303030695C3030307A5C303030655C303030645C3030305C3034305C3030306D5C3030306F5C303030645C303030655C3030306C}
\babel@aux{english}{}
\@writefile{toc}{\contentsline {section}{\numberline {1}Introduction}{1}{section.1}\protected@file@percent }
\newlabel{sec:intro}{{1}{1}{Introduction}{section.1}{}}
\@writefile{toc}{\contentsline {section}{\numberline {2}Modelling and Linearisation}{1}{section.2}\protected@file@percent }
\newlabel{sec:q1}{{2}{1}{Modelling and Linearisation}{section.2}{}}
\@writefile{toc}{\contentsline {subsection}{\numberline {2.1}Operating point}{1}{subsection.2.1}\protected@file@percent }
\newlabel{eq:q1_balance}{{1}{1}{Operating point}{equation.1}{}}
\@writefile{lof}{\contentsline {figure}{\numberline {1}{\ignorespaces \small  \add@extra@listi {sml}Actuator margins $\min (\bar  u_i,\,1-\bar  u_i)$ against the free parameter $u_2$, with $\bar  u_1,\bar  u_3$ set by the steady-state balance~\eqref  {eq:q1_balance}. }}{1}{figure.caption.1}\protected@file@percent }
\providecommand*\caption@xref[2]{\@setref\relax\@undefined{#1}}
\newlabel{fig:q1_margins}{{1}{1}{\small Actuator margins $\min (\bar u_i,\,1-\bar u_i)$ against the free parameter $u_2$, with $\bar u_1,\bar u_3$ set by the steady-state balance~\eqref {eq:q1_balance}}{figure.caption.1}{}}
\newlabel{eq:q1_eqpoint}{{2}{1}{Operating point}{equation.2}{}}
\@writefile{toc}{\contentsline {subsection}{\numberline {2.2}Linearized model}{1}{subsection.2.2}\protected@file@percent }
\newlabel{eq:q1_nl}{{3}{1}{Linearized model}{equation.3}{}}
\newlabel{eq:q1_jac}{{4}{1}{Linearized model}{equation.4}{}}
\newlabel{eq:q1_beta}{{5}{1}{Linearized model}{equation.5}{}}
\BKM@entry{id=5,dest={73756273656374696F6E2E322E33},srcline={299}}{5C3337365C3337375C3030304F5C303030705C303030655C3030306E5C3030302D5C3030306C5C3030306F5C3030306F5C303030705C3030305C3034305C303030645C303030795C3030306E5C303030615C3030306D5C303030695C303030635C30303073}
\BKM@entry{id=6,dest={73756273656374696F6E2E322E34},srcline={317}}{5C3337365C3337375C3030304C5C303030695C3030306E5C303030655C303030615C303030725C3030305C3034305C303030765C303030655C303030725C303030735C303030755C303030735C3030305C3034305C3030306E5C3030306F5C3030306E5C3030306C5C303030695C3030306E5C303030655C303030615C303030725C3030305C3034305C303030725C303030655C303030735C303030705C3030306F5C3030306E5C303030735C30303065}
\@writefile{lof}{\contentsline {figure}{\numberline {2}{\ignorespaces Free response of the linearised and nonlinear models under the constant equilibrium input, for different initial deviations. Far from the operating point the nonlinear response fills faster, as the linearisation overestimates the $\sqrt  {h}$ outflow.}}{2}{figure.caption.2}\protected@file@percent }
\newlabel{fig:q1_cmp}{{2}{2}{Free response of the linearised and nonlinear models under the constant equilibrium input, for different initial deviations. Far from the operating point the nonlinear response fills faster, as the linearisation overestimates the $\sqrt {h}$ outflow}{figure.caption.2}{}}
\newlabel{eq:q1_entries}{{6}{2}{Linearized model}{equation.6}{}}
\newlabel{eq:q1_A}{{7}{2}{Linearized model}{equation.7}{}}
\newlabel{eq:q1_B}{{8}{2}{Linearized model}{equation.8}{}}
\@writefile{toc}{\contentsline {subsection}{\numberline {2.3}Open-loop dynamics}{2}{subsection.2.3}\protected@file@percent }
\newlabel{eq:q1_sameEW}{{9}{2}{Open-loop dynamics}{equation.9}{}}
\@writefile{toc}{\contentsline {subsection}{\numberline {2.4}Linear versus nonlinear response}{2}{subsection.2.4}\protected@file@percent }
\BKM@entry{id=7,dest={73656374696F6E2E33},srcline={386}}{5C3337365C3337375C303030535C303030745C303030615C303030745C303030655C3030302D5C303030465C303030655C303030655C303030645C303030625C303030615C303030635C3030306B5C3030305C3034305C303030445C303030655C303030735C303030695C303030675C3030306E5C3030305C3034305C303030765C303030695C303030615C3030305C3034305C3030304C5C3030304D5C303030495C30303073}
\BKM@entry{id=8,dest={73756273656374696F6E2E332E31},srcline={391}}{5C3337365C3337375C303030465C303030725C3030306F5C3030306D5C3030305C3034305C3030304C5C303030795C303030615C303030705C303030755C3030306E5C3030306F5C303030765C3030305C3034305C303030695C3030306E5C303030655C303030715C303030755C303030615C3030306C5C303030695C303030745C303030795C3030305C3034305C303030745C3030306F5C3030305C3034305C3030304C5C3030304D5C30303049}
\citation{zaccarian_lmi}
\BKM@entry{id=9,dest={73756273656374696F6E2E332E32},srcline={425}}{5C3337365C3337375C303030415C303030725C303030625C303030695C303030745C303030725C303030615C303030725C303030795C3030305C3034305C303030735C303030745C303030615C303030625C303030695C3030306C5C303030695C3030307A5C303030655C30303072}
\newlabel{fig:q1_decomp}{{3a}{3}{Error of $\dot h_2$ split into inflow ($+u_2a_2\varphi _1$) and outflow ($-u_3a_3\varphi _2$) channels}{figure.caption.3}{}}
\newlabel{sub@fig:q1_decomp}{{a}{3}{Error of $\dot h_2$ split into inflow ($+u_2a_2\varphi _1$) and outflow ($-u_3a_3\varphi _2$) channels}{figure.caption.3}{}}
\newlabel{fig:q1_err}{{3b}{3}{Norm of the mismatch $\|\bm {x}_{\mathrm {nl}}-\bm {x}_{\mathrm {lin}}\|$ for the two initial deviations (separate axes)}{figure.caption.3}{}}
\newlabel{sub@fig:q1_err}{{b}{3}{Norm of the mismatch $\|\bm {x}_{\mathrm {nl}}-\bm {x}_{\mathrm {lin}}\|$ for the two initial deviations (separate axes)}{figure.caption.3}{}}
\@writefile{lof}{\contentsline {figure}{\numberline {3}{\ignorespaces Linearisation error. Channel decomposition~(a) and overall norm~(b); the peak scales by roughly two orders of magnitude between the small and large deviations.}}{3}{figure.caption.3}\protected@file@percent }
\newlabel{fig:q1_error}{{3}{3}{Linearisation error. Channel decomposition~(a) and overall norm~(b); the peak scales by roughly two orders of magnitude between the small and large deviations}{figure.caption.3}{}}
\@writefile{toc}{\contentsline {section}{\numberline {3}State-Feedback Design via LMIs}{3}{section.3}\protected@file@percent }
\newlabel{sec:q2}{{3}{3}{State-Feedback Design via LMIs}{section.3}{}}
\@writefile{toc}{\contentsline {subsection}{\numberline {3.1}From Lyapunov inequality to LMI}{3}{subsection.3.1}\protected@file@percent }
\newlabel{eq:q2_lyap}{{11}{3}{From Lyapunov inequality to LMI}{equation.11}{}}
\newlabel{eq:q2_baselmi}{{12}{3}{From Lyapunov inequality to LMI}{equation.12}{}}
\newlabel{eq:q2_recover}{{13}{3}{From Lyapunov inequality to LMI}{equation.13}{}}
\BKM@entry{id=10,dest={73756273656374696F6E2E332E33},srcline={439}}{5C3337365C3337375C303030475C303030755C303030615C303030725C303030615C3030306E5C303030745C303030655C303030655C303030645C3030305C3034305C303030635C3030306F5C3030306E5C303030765C303030655C303030725C303030675C303030655C3030306E5C303030635C303030655C3030305C3034305C303030725C303030615C303030745C30303065}
\citation{zaccarian_lmi}
\BKM@entry{id=11,dest={73756273656374696F6E2E332E34},srcline={450}}{5C3337365C3337375C3030304D5C303030695C3030306E5C303030695C3030306D5C303030755C3030306D5C3030302D5C303030675C303030615C303030695C3030306E5C3030305C3034305C303030645C303030655C303030735C303030695C303030675C3030306E5C3030305C3034305C303030775C303030695C303030745C303030685C3030305C3034305C303030615C303030635C303030745C303030755C303030615C303030745C3030306F5C303030725C3030305C3034305C303030615C303030755C303030745C303030685C3030306F5C303030725C303030695C303030745C30303079}
\citation{zaccarian_lmi}
\BKM@entry{id=12,dest={73756273656374696F6E2E332E35},srcline={469}}{5C3337365C3337375C303030435C3030306F5C3030306D5C303030705C303030615C303030725C303030695C303030735C3030306F5C3030306E}
\BKM@entry{id=13,dest={73656374696F6E2E34},srcline={520}}{5C3337365C3337375C3030304F5C303030625C303030735C303030655C303030725C303030765C303030655C303030725C3030305C3034305C303030445C303030655C303030735C303030695C303030675C3030306E}
\BKM@entry{id=14,dest={73756273656374696F6E2E342E31},srcline={523}}{5C3337365C3337375C303030575C303030685C303030695C303030635C303030685C3030305C3034305C303030735C303030655C3030306E5C303030735C3030306F5C303030725C3030305C3034305C303030665C303030615C303030695C3030306C5C303030755C303030725C303030655C3030305C3034305C303030695C303030735C3030305C3034305C303030725C303030655C303030635C3030306F5C303030765C303030655C303030725C303030615C303030625C3030306C5C30303065}
\@writefile{toc}{\contentsline {subsection}{\numberline {3.2}Arbitrary stabilizer}{4}{subsection.3.2}\protected@file@percent }
\@writefile{lof}{\contentsline {figure}{\numberline {4}{\ignorespaces Control signals for the three designs. NB: differing vertical scales, $u_1$ on the left axis and $u_2,u_3$ on the right. Different time scales capture the transient. }}{4}{figure.caption.4}\protected@file@percent }
\newlabel{fig:q2_control}{{4}{4}{Control signals for the three designs. NB: differing vertical scales, $u_1$ on the left axis and $u_2,u_3$ on the right. Different time scales capture the transient}{figure.caption.4}{}}
\@writefile{toc}{\contentsline {subsection}{\numberline {3.3}Guaranteed convergence rate}{4}{subsection.3.3}\protected@file@percent }
\newlabel{eq:q2_lmi_b}{{14}{4}{Guaranteed convergence rate}{equation.14}{}}
\@writefile{toc}{\contentsline {subsection}{\numberline {3.4}Minimum-gain design with actuator authority}{4}{subsection.3.4}\protected@file@percent }
\newlabel{eq:q2_lmi_c}{{15}{4}{Minimum-gain design with actuator authority}{equation.15}{}}
\@writefile{toc}{\contentsline {subsection}{\numberline {3.5}Comparison}{4}{subsection.3.5}\protected@file@percent }
\@writefile{lot}{\contentsline {table}{\numberline {1}{\ignorespaces Closed-loop comparison on the nonlinear plant, from $\delta \bm  {x}_0=-0.05$~m. Settling time $t_s$ is the time to enter and remain within $2\%$ of the initial step size; a $1.4$~s floor is applied to avoid the overshoot being mistaken with settling time. The commands of $K_a$ and $K_b$ leave the admissible range $[0,1]$, so their settling times are not physically realisable; only $K_c$ respects the actuator bounds.$u_{\text  {max}}$ and $u_{\text  {min}}$ are considered across all three controls. }}{4}{table.caption.5}\protected@file@percent }
\newlabel{tab:q2_cmp}{{1}{4}{Closed-loop comparison on the nonlinear plant, from $\delta \bm {x}_0=-0.05$~m. Settling time $t_s$ is the time to enter and remain within $2\%$ of the initial step size; a $1.4$~s floor is applied to avoid the overshoot being mistaken with settling time. The commands of $K_a$ and $K_b$ leave the admissible range $[0,1]$, so their settling times are not physically realisable; only $K_c$ respects the actuator bounds.$u_{\text {max}}$ and $u_{\text {min}}$ are considered across all three controls}{table.caption.5}{}}
\BKM@entry{id=15,dest={73756273656374696F6E2E342E32},srcline={529}}{5C3337365C3337375C3030304F5C303030625C303030735C303030655C303030725C303030765C303030655C303030725C3030305C3034305C303030675C303030615C303030695C3030306E}
\citation{zaccarian_lmi}
\BKM@entry{id=16,dest={73756273656374696F6E2E342E33},srcline={555}}{5C3337365C3337375C303030455C303030735C303030745C303030695C3030306D5C303030615C303030745C303030695C3030306F5C3030306E}
\newlabel{fig:q2_state_short}{{5a}{5}{Short timescale ($0$--$15$~s): $K_a$, $K_b$ settle quickly but overshoot; $K_c$ rises gently}{figure.caption.6}{}}
\newlabel{sub@fig:q2_state_short}{{a}{5}{Short timescale ($0$--$15$~s): $K_a$, $K_b$ settle quickly but overshoot; $K_c$ rises gently}{figure.caption.6}{}}
\newlabel{fig:q2_state_long}{{5b}{5}{Full timescale ($0$--$3000$~s): $K_c$ settles in $\approx 385$~s against the open-loop $\approx 2840$~s}{figure.caption.6}{}}
\newlabel{sub@fig:q2_state_long}{{b}{5}{Full timescale ($0$--$3000$~s): $K_c$ settles in $\approx 385$~s against the open-loop $\approx 2840$~s}{figure.caption.6}{}}
\@writefile{lof}{\contentsline {figure}{\numberline {5}{\ignorespaces Closed-loop state response on the nonlinear plant from $\delta \bm  {x}_0=-0.05$~m, shown at two timescales. The apparent instantaneous settling of $K_a$ and $K_b$ in~(b) is obtained only with the inadmissible commands of Figure~\ref {fig:q2_control}.}}{5}{figure.caption.6}\protected@file@percent }
\newlabel{fig:q2_state}{{5}{5}{Closed-loop state response on the nonlinear plant from $\delta \bm {x}_0=-0.05$~m, shown at two timescales. The apparent instantaneous settling of $K_a$ and $K_b$ in~(b) is obtained only with the inadmissible commands of Figure~\ref {fig:q2_control}}{figure.caption.6}{}}
\@writefile{toc}{\contentsline {section}{\numberline {4}Observer Design}{5}{section.4}\protected@file@percent }
\newlabel{sec:q3}{{4}{5}{Observer Design}{section.4}{}}
\@writefile{toc}{\contentsline {subsection}{\numberline {4.1}Which sensor failure is recoverable}{5}{subsection.4.1}\protected@file@percent }
\@writefile{toc}{\contentsline {subsection}{\numberline {4.2}Observer gain}{5}{subsection.4.2}\protected@file@percent }
\newlabel{eq:q3_obs}{{16}{5}{Observer gain}{equation.16}{}}
\newlabel{eq:q3_lmi}{{17}{5}{Observer gain}{equation.17}{}}
\@writefile{toc}{\contentsline {subsection}{\numberline {4.3}Estimation}{5}{subsection.4.3}\protected@file@percent }
\BKM@entry{id=17,dest={73656374696F6E2E35},srcline={586}}{5C3337365C3337375C3030304F5C303030625C303030735C303030655C303030725C303030765C303030655C303030725C3030302D5C303030625C303030615C303030735C303030655C303030645C3030305C3034305C303030665C303030655C303030655C303030645C303030625C303030615C303030635C3030306B}
\bibstyle{abbrvnat}
\bibdata{references}
\bibcite{zaccarian_lmi}{{1}{}{{Zaccarian and Sferlazza}}{{}}}
\newlabel{fig:q3_track}{{6a}{6}{Observer estimate against the true state (only $h_2$ measured). The $h_2$ estimate is exact from the start; the $h_1$ estimate starts at $0.2$~m, dips near $t=1.3$~s, then converges}{figure.caption.7}{}}
\newlabel{sub@fig:q3_track}{{a}{6}{Observer estimate against the true state (only $h_2$ measured). The $h_2$ estimate is exact from the start; the $h_1$ estimate starts at $0.2$~m, dips near $t=1.3$~s, then converges}{figure.caption.7}{}}
\newlabel{fig:q3_fb}{{6b}{6}{Closed loop with true-state against observer-based feedback, both using $K_c$. Observer-based is slightly faster on $h_1$, then matches as the estimate becomes exact}{figure.caption.7}{}}
\newlabel{sub@fig:q3_fb}{{b}{6}{Closed loop with true-state against observer-based feedback, both using $K_c$. Observer-based is slightly faster on $h_1$, then matches as the estimate becomes exact}{figure.caption.7}{}}
\@writefile{lof}{\contentsline {figure}{\numberline {6}{\ignorespaces Observer behaviour: state estimation~(a) and the resulting observer-based closed loop~(b).}}{6}{figure.caption.7}\protected@file@percent }
\newlabel{fig:q3}{{6}{6}{Observer behaviour: state estimation~(a) and the resulting observer-based closed loop~(b)}{figure.caption.7}{}}
\@writefile{toc}{\contentsline {section}{\numberline {5}Observer-based feedback}{6}{section.5}\protected@file@percent }
\newlabel{eq:q3_fb}{{18}{6}{Observer-based feedback}{equation.18}{}}
\@writefile{lot}{\contentsline {table}{\numberline {2}{\ignorespaces Full-state against observer-based feedback, both using $K_c$, from $\delta \bm  {x}_0=-0.05$~m (same $1.4$~s settling floor as Table~\ref {tab:q2_cmp})}}{6}{table.caption.8}\protected@file@percent }
\newlabel{tab:q3_cmp}{{2}{6}{Full-state against observer-based feedback, both using $K_c$, from $\delta \bm {x}_0=-0.05$~m (same $1.4$~s settling floor as Table~\ref {tab:q2_cmp})}{table.caption.8}{}}
\gdef \@abspage@last{6}
```

# report/after-figureplacement/report_fp.bbl

```
\begin{thebibliography}{1}
\providecommand{\natexlab}[1]{#1}
\providecommand{\url}[1]{\texttt{#1}}
\expandafter\ifx\csname urlstyle\endcsname\relax
  \providecommand{\doi}[1]{doi: #1}\else
  \providecommand{\doi}{doi: \begingroup \urlstyle{rm}\Url}\fi

\bibitem[Zaccarian and Sferlazza()]{zaccarian_lmi}
L.~Zaccarian and A.~Sferlazza.
\newblock Analysis and synthesis of feedback systems: Quadratic functions and {LMIs}.
\newblock Lecture notes.
\newblock Automatic Control course, University of Trento.

\end{thebibliography}
```

# report/after-figureplacement/report_fp.blg

```
This is BibTeX, Version 0.99d (TeX Live 2025)
Capacity: max_strings=200000, hash_size=200000, hash_prime=170003
The top-level auxiliary file: report_fp.aux
The style file: abbrvnat.bst
Database file #1: references.bib
Warning--empty year in zaccarian_lmi
You've used 1 entry,
            2773 wiz_defined-function locations,
            597 strings with 4743 characters,
and the built_in function-call counts, 400 in all, are:
= -- 33
> -- 21
< -- 1
+ -- 7
- -- 6
* -- 23
:= -- 65
add.period$ -- 4
call.type$ -- 1
change.case$ -- 5
chr.to.int$ -- 1
cite$ -- 3
duplicate$ -- 15
empty$ -- 33
format.name$ -- 9
if$ -- 80
int.to.chr$ -- 1
int.to.str$ -- 1
missing$ -- 0
newline$ -- 14
num.names$ -- 4
pop$ -- 17
preamble$ -- 1
purify$ -- 4
quote$ -- 0
skip$ -- 13
stack$ -- 0
substring$ -- 2
swap$ -- 1
text.length$ -- 0
text.prefix$ -- 0
top$ -- 0
type$ -- 11
warning$ -- 1
while$ -- 3
width$ -- 0
write$ -- 20
(There was 1 warning)
```

# report/after-figureplacement/report_fp.fdb_latexmk

```
# Fdb version 4
["bibtex report_fp"] 1785446013.27528 "report_fp.aux" "report_fp.bbl" "report_fp" 1785446014.331 0
  "./references.bib" 1785441267.66428 513 322299608c31c4c4682038b525a470bc ""
  "/usr/local/texlive/2025basic/texmf-dist/bibtex/bst/natbib/abbrvnat.bst" 1233624470 26532 66bed5cc59495faaa3f9ae9af9fa5e42 ""
  "report_fp.aux" 1785446014.23401 16602 d412ef4121943406e30a8bf1b8d49c04 "pdflatex"
  (generated)
  "report_fp.bbl"
  "report_fp.blg"
  (rewritten before read)
["pdflatex"] 1785446013.29181 "/Users/vi/Documents/Projects/LMI_PROJECT/report/after-figureplacement/report_fp.tex" "report_fp.pdf" "report_fp" 1785446014.33112 0
  "/Users/vi/Documents/Projects/LMI_PROJECT/report/after-figureplacement/report_fp.tex" 1785446011.93699 33004 3c54951de9aba1410ba8fcd1c30f3e06 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/enc/dvips/base/8r.enc" 1165713224 4850 80dc9bab7f31fb78a000ccfed0e27cab ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvb8r.tfm" 1136768653 4484 b828043cbd581d289d955903c1339981 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvb8t.tfm" 1136768653 6628 34c39492c0adc454c1c199922bba8363 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8r.tfm" 1136768653 4712 9ef4d7d106579d4b136e1529e1a4533c ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8t.tfm" 1136768653 7040 b2bd27e2bfe6f6948cbc3239cae7444f ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvro8r.tfm" 1136768653 4964 f223217e5e1f85fa3742fb0480aba9e8 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvro8t.tfm" 1136768653 7292 d408122871055bbc541d7280c154e210 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/jknappen/ec/ecrm1000.tfm" 1136768653 3584 adb004a0c8e7c46ee66cad73671f37b4 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmbsy5.tfm" 1246382020 1120 1e8878807317373affa7f7bba4cf2f6a ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmbsy6.tfm" 1246382020 1124 14ccf5552bc7f77ca02a8a402bea8bfb ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmbsy7.tfm" 1246382020 1120 7f9f170e8aa57527ad6c49feafd45d54 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmbsy8.tfm" 1246382020 1120 200be8b775682cdf80acad4be5ef57e4 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmbsy9.tfm" 1246382020 1112 cbc11b646ccc26599775160605aaee3a ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmex7.tfm" 1246382020 1004 54797486969f23fa377b128694d548df ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmex8.tfm" 1246382020 988 bdf658c3bfc2d96d3c8b02cfc1c94c20 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmex9.tfm" 1246382020 996 a18840b13b499c08ac2de96a99eda4bc ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmmib5.tfm" 1246382020 1496 c79f6914c6d39ffb3759967363d1be79 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmmib6.tfm" 1246382020 1516 a3bf6a5e7ec4401b1f52092dfaaed242 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmmib7.tfm" 1246382020 1508 6e807ff901c35a5f1fde0ca275533df8 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmmib8.tfm" 1246382020 1528 dab402b9d3774ca98baa037071cee7ae ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmmib9.tfm" 1246382020 1528 159d57adcba064aab4277245c826577d ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbsy10.tfm" 1136768653 1116 4e6ba9d7914baa6482fd69f67d126380 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbx10.tfm" 1136768653 1328 c834bbb027764024c09d3d2bf908b5f0 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbx12.tfm" 1136768653 1324 c910af8c371558dc20f2d7822f66fe64 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbx5.tfm" 1136768653 1332 f817c21a1ba54560425663374f1b651a ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbx6.tfm" 1136768653 1344 8a0be4fe4d376203000810ad4dc81558 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbx7.tfm" 1136768653 1336 3125ccb448c1a09074e3aa4a9832f130 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbx8.tfm" 1136768653 1332 1fde11373e221473104d6cc5993f046e ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbx9.tfm" 1136768653 1328 5442e22a7072966dbaf88ca900acf3f0 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmex10.tfm" 1136768653 992 662f679a0b3d2d53c1b94050fdaa3f50 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmmi12.tfm" 1136768653 1524 4414a8315f39513458b80dfc63bff03a ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmmi6.tfm" 1136768653 1512 f21f83efb36853c0b70002322c1ab3ad ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmmi8.tfm" 1136768653 1520 eccf95517727cb11801f4f1aee3a21b4 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmmi9.tfm" 1136768653 1524 d89e2d087a9828407a196f428428ef4a ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmmib10.tfm" 1136768653 1524 554068197b70979a55370e6c6495f441 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmr12.tfm" 1136768653 1288 655e228510b4c2a1abe905c368440826 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmr6.tfm" 1136768653 1300 b62933e007d01cfd073f79b963c01526 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmr8.tfm" 1136768653 1292 21c1c5bfeaebccffdb478fd231a0997d ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmr9.tfm" 1136768653 1292 6b21b9c2c7bebb38aa2273f7ca0fb3af ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmss10.tfm" 1136768653 1316 b636689f1933f24d1294acdf6041daaa ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmss12.tfm" 1136768653 1324 37b971caf729d7edd9cbb9f9b0ea76eb ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmss8.tfm" 1136768653 1296 d77f431d10d47c8ea2cc18cf45346274 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmss9.tfm" 1136768653 1320 49357c421c0d469f88b867dd0c3d10e8 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmsy10.tfm" 1136768653 1124 6c73e740cf17375f03eec0ee63599741 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmsy6.tfm" 1136768653 1116 933a60c408fc0a863a92debe84b2d294 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmsy8.tfm" 1136768653 1120 8b7d695260f3cff42e636090a8002094 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmsy9.tfm" 1136768653 1116 25a7bf822c58caf309a702ef79f4afbb ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmbx10.pfb" 1248133631 34811 78b52f49e893bcba91bd7581cdc144c0 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmex10.pfb" 1248133631 30251 6afa5cb1d0204815a708a080681d4674 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi10.pfb" 1248133631 36299 5f9df58c2139e7edcf37c8fca4bd384d ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi12.pfb" 1248133631 36741 fa121aac0049305630cf160b86157ee4 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi6.pfb" 1248133631 37166 8ab3487cbe3ab49ebce74c29ea2418db ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi7.pfb" 1248133631 36281 c355509802a035cadc5f15869451dcee ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi8.pfb" 1248133631 35469 70d41d2b9ea31d5d813066df7c99281c ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi9.pfb" 1248133631 36094 798f80770b3b148ceedd006d487db67c ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmib10.pfb" 1248133631 36912 b448ef9ad9d7228ec3c6e71005136d55 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr10.pfb" 1248133631 35752 024fb6c41858982481f6968b5fc26508 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr5.pfb" 1248133631 31809 8670ca339bf94e56da1fc21c80635e2a ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr6.pfb" 1248133631 32734 69e00a6b65cedb993666e42eedb3d48f ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr7.pfb" 1248133631 32762 224316ccc9ad3ca0423a14971cfa7fc1 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr8.pfb" 1248133631 32726 0a1aea6fcd6468ee2cf64d891f5c43c8 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr9.pfb" 1248133631 33993 9b89b85fd2d9df0482bd47194d1d3bf3 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmss10.pfb" 1248133631 24457 5cbb7bdf209d5d1ce9892a9b80a307cc ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmsy10.pfb" 1248133631 32569 5e5ddc8df908dea60932f3c484a54c0d ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmsy7.pfb" 1248133631 32716 08e384dc442464e7285e891af9f45947 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmsy8.pfb" 1248133631 32626 4f5c1b83753b1dd3a97d1b399a005b4b ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmsy9.pfb" 1248133631 32442 c975af247b6702f7ca0c299af3616b80 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cmextra/cmmib7.pfb" 1248133631 39386 0bb525d495ceab21fe797301a586e0d7 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cmextra/cmmib9.pfb" 1248133631 36510 24bd6e95867d9b5c2691bf9144a22cda ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/urw/helvetic/uhvb8a.pfb" 1136849748 35941 f27169cc74234d5bd5e4cca5abafaabb ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/urw/helvetic/uhvr8a.pfb" 1136849748 44648 23115b2a545ebfe2c526c3ca99db8b95 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/type1/urw/helvetic/uhvro8a.pfb" 1136849748 48169 b4fd9d908b9ee8c65d4305ad39071c5e ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/vf/adobe/helvetic/phvb8t.vf" 1136768653 2340 0efed6a948c3c37d870e4e7ddb85c7c3 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/vf/adobe/helvetic/phvr8t.vf" 1136768653 2344 44ff28c9ef2fc97180cd884f900fee71 ""
  "/usr/local/texlive/2025basic/texmf-dist/fonts/vf/adobe/helvetic/phvro8t.vf" 1136768653 2348 717f0e01beff9c30c0df4945e471886f ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/context/base/mkii/supp-pdf.mkii" 1461363279 71627 94eb9990bed73c364d7f53f960cc8c5b ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/atbegshi/atbegshi.sty" 1575674566 24708 5584a51a7101caf7e6bbf1fc27d8f7b1 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/babel-english/english.ldf" 1496785618 7008 9ff5fdcc865b01beca2b0fe4a46231d4 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/babel/babel.sty" 1739572459 143388 b008e1666520ff43d0878397f2926242 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/babel/locale/en/babel-en.ini" 1661803479 3966 caeee5a9e5771d4446aa1ca9015ba1b2 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/babel/locale/en/babel-english.tex" 1498512262 336 ed676b5e7dfd862bc78d634f6a973f37 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/babel/txtbabel.def" 1735765002 6945 a248d839b1f26b388440c973966647b5 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/bigintcalc/bigintcalc.sty" 1576625341 40635 c40361e206be584d448876bba8a64a3b ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/bitset/bitset.sty" 1576016050 33961 6b5c75130e435b2bfdb9f480a09a39f9 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/gettitlestring/gettitlestring.sty" 1576625223 8371 9d55b8bd010bc717624922fb3477d92e ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/iftex/iftex.sty" 1734129479 7984 7dbb9280f03c0a315425f1b4f35d43ee ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/iftex/ifvtex.sty" 1572645307 1057 525c2192b5febbd8c1f662c9468335bb ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/infwarerr/infwarerr.sty" 1575499628 8356 7bbb2c2373aa810be568c29e333da8ed ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/intcalc/intcalc.sty" 1576625065 31769 002a487f55041f8e805cfbf6385ffd97 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/kvdefinekeys/kvdefinekeys.sty" 1576878844 5412 d5a2436094cd7be85769db90f29250a6 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/ltxcmds/ltxcmds.sty" 1701727651 17865 1a9bd36b4f98178fa551aca822290953 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/pdfescape/pdfescape.sty" 1576015897 19007 15924f7228aca6c6d184b115f4baa231 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/pdftexcmds/pdftexcmds.sty" 1593379760 20089 80423eac55aa175305d35b49e04fe23b ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/stringenc/stringenc.sty" 1575152242 21514 b7557edcee22835ef6b03ede1802dad4 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/generic/uniquecounter/uniquecounter.sty" 1576624663 7008 f92eaa0a3872ed622bbf538217cd2ab7 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsbsy.sty" 1717359999 2222 2166a1f7827be30ddc30434e5efcee1b ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsgen.sty" 1717359999 4173 d22509bc0c91281d991b2de7c88720dd ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsmath.sty" 1730928152 88370 c780f23aea0ece6add91e09b44dca2cd ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsopn.sty" 1717359999 4474 23ca1d3a79a57b405388059456d0a8df ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amstext.sty" 1717359999 2444 71618ea5f2377e33b04fb97afdd0eac2 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/atveryend/atveryend.sty" 1728505250 1695 be6b4d13b33db697fd3fd30b24716c1a ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/base/atbegshi-ltx.sty" 1738182759 2963 d8ec5a1b4e0a106c5c737900202763e4 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/base/atveryend-ltx.sty" 1738182759 2378 14b657ee5031da98cf91648f19642694 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/base/fontenc.sty" 1738182759 5275 0d62fb62162c7ab056e941ef18c5076d ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/base/inputenc.sty" 1738182759 5048 0270515b828149155424600fd2d58ac5 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/bookmark/bkm-pdftex.def" 1702241854 8818 aa5157b46368efebf023abff55611467 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/bookmark/bookmark.sty" 1702241854 18245 97e6be180cf07bb6f7008cfdaaecfce5 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/booktabs/booktabs.sty" 1579038678 6078 f1cb470c9199e7110a27851508ed7a5c ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/caption-koma.sto" 1696191071 6333 57259249669cec150aa89a6e45395d3f ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/caption.sty" 1696191071 56128 c2ccf1a29d78c33bc553880402e4fb9a ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/caption3.sty" 1696191071 72619 ee90b6612147680fd73c3b1406a74245 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/subcaption.sty" 1690576852 12494 0c0cdb824278a4d51cefeb2e79901315 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/enumitem/enumitem.sty" 1738874546 52272 63d293bc0d496619edb57585740861a2 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/epstopdf-pkg/epstopdf-base.sty" 1579991033 13886 d1306dcf79a944f6988e688c1785f9ce ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/etoolbox/etoolbox.sty" 1739306980 46850 d87daedc2abdc653769a6f1067849fe0 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/geometry/geometry.sty" 1578002852 41601 9cf6c5257b1bc7af01a58859749dd37a ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-cfg/color.cfg" 1459978653 1213 620bba36b25224fa9b7e1ccb4ecb76fd ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-cfg/graphics.cfg" 1465944070 1224 978390e9c2234eab29404bc21b268d1e ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-def/pdftex.def" 1713382759 19440 9da9dcbb27470349a580fca7372d454b ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/dvipsnam.def" 1717359999 5009 1ca6c92de20f17acac654d4e0598c200 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/graphics.sty" 1730496337 18363 dee506cb8d56825d8a4d020f5d5f8704 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/graphicx.sty" 1717359999 8010 6f2ad8c2b2ffbd607af6475441c7b5e4 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/keyval.sty" 1717359999 2671 70891d50dac933918b827d326687c6e8 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/mathcolor.ltx" 1667332637 2885 9c645d672ae17285bba324998918efd8 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/trig.sty" 1717359999 4023 2c9f39712cf7b43d3eb93a8bbd5c8f67 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/grfext/grfext.sty" 1575499774 7133 b94bbacbee6e4fdccdc7f810b2aec370 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/hycolor/hycolor.sty" 1580250785 17914 4c28a13fc3d975e6e81c9bea1d697276 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/hpdftex.def" 1730838014 48154 82da9991b9f0390b3a9d3af6c8618af4 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/hyperref.sty" 1730838014 222112 c22dbd2288f89f7ba942ac22f7d00f11 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/nameref.sty" 1705871765 11026 182c63f139a71afd30a28e5f1ed2cd1c ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/pd1enc.def" 1730838014 14249 ff700eb13ce975a424b2dd99b1a83044 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/puenc.def" 1730838014 117112 7533bff456301d32e6d6356fad15f543 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrartcl.cls" 1729800159 244773 292e8a87db6633ccdc51a3b380022269 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrbase.sty" 1729800159 100966 88ffa900fd38336daf0681d35d342ea9 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrkbase.sty" 1729800159 22043 fd1e842c83a1e0d0711eaa2cec7f7705 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlayer-scrpage.sty" 1729800159 60545 b309f0e93fdc770348e7c62115822a88 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlayer.sty" 1729800159 87414 5835b180c5ceb8e3c1d388575ab68ae3 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlfile-hook.sty" 1729800159 11185 08107e8d26d093ccd4c424c2b74809f6 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlfile.sty" 1729800159 3328 17a5a2d4f4e9d388803c10ac9fffe9d3 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlogo.sty" 1729800159 2162 e219c1ddf641a7cd0ee0103af3ac7f3d ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrsize10pt.clo" 1729800159 6190 ec1ef4257e686a475ebcb771c2149a92 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/tocbasic.sty" 1729800159 111413 13a9e9f758789eb6259d46ebbe3c3c6d ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/typearea.sty" 1729800159 58618 afeffc7e809c19b544ab28b3b61ad0a7 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/kvoptions/kvoptions.sty" 1655478651 22555 6d8e155cfef6d82c3d5c742fea7c992e ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/kvsetkeys/kvsetkeys.sty" 1665067230 13815 760b0c02f691ea230f5359c4e1de23a7 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/l3backend/l3backend-pdftex.def" 1716410060 29785 9f93ab201fe5dd053afcc6c1bcf7d266 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/latexconfig/epstopdf-sys.cfg" 1279039959 678 4792914a8f45be57bb98413425e4c7af ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/mathtools/mathtools.sty" 1728349409 62711 3b0ed1d68b598a75677c10d16646d5e5 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/mathtools/mhsetup.sty" 1616101747 5582 a43dedf8e5ec418356f1e9dfe5d29fc3 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/microtype-pdftex.def" 1739394495 49650 26a5e891c8da4553198575ba0517c0e5 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/microtype.cfg" 1739394495 27015 bd167d0154f271c424b157d8894ae4a4 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/microtype.sty" 1739394495 102775 6624742dafeb6f262a13657f9f77f048 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/mt-cmr.cfg" 1739394495 22906 b0be544692bb405a84d147a96af5d777 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/natbib/natbib.sty" 1291685959 45456 1c8843383c0bd05870c45fa0ebea6cc2 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/psnfss/helvet.sty" 1586716065 1499 de0ad166b701b820e03588a29bb30798 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/psnfss/t1phv.fd" 1586716065 1483 47067fbe7c3ffed1ede7aaa7b8549d7a ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/ragged2e/ragged2e.sty" 1687466457 10085 3dececaeed388e31285c366e124c0c4e ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/refcount/refcount.sty" 1576624809 9878 9e94e8fa600d95f9c7731bb21dfb67a4 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/rerunfilecheck/rerunfilecheck.sty" 1657483315 9714 ba3194bd52c8499b3f1e3eb91d409670 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/array.sty" 1730496337 14552 27664839421e418b87f56fa4c6f66b1a ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/bm.sty" 1717359999 13236 1aba485b47e3c79a47dc0906cc971820 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/calc.sty" 1717359999 10214 61188260d324e94bc2f66825d7d3fdf4 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/multicol.sty" 1730496337 32879 f4cfa2226c995da5b04a89ee47abf182 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/tabularx.sty" 1717359999 7243 e5dac1240636811edb77568b81818372 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/url/url.sty" 1388531844 12796 8edb7d69a20b857904dd0ea757c14ec9 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/wrapfig/wrapfig.sty" 1137111090 26220 3701aebf80ccdef248c0c20dd062fea9 ""
  "/usr/local/texlive/2025basic/texmf-dist/tex/latex/xcolor/xcolor.sty" 1727642399 55384 b454dec21c2d9f45ec0b793f0995b992 ""
  "/usr/local/texlive/2025basic/texmf-dist/web2c/texmf.cnf" 1739380943 42148 61becc7c670cd061bb319c643c27fdd4 ""
  "/usr/local/texlive/2025basic/texmf-var/fonts/map/pdftex/updmap/pdftex.map" 1767117959 365213 3ebf252b02690d51af809a89a8c4f06f ""
  "/usr/local/texlive/2025basic/texmf-var/web2c/pdftex/pdflatex.fmt" 1752831409 2663826 3324249ba99bbe56829524fcf0d7692a ""
  "/usr/local/texlive/2025basic/texmf.cnf" 1741448314 628 81d1ada0f601fb1ffa62f1323192d946 ""
  "figures/find_u2.png" 1785441939.02641 46810 51b72161e706318317c0e661e0991605 ""
  "figures/q1_cmp.png" 1785443119.968 135155 51bca834469111386a7ace14b4e82d3c ""
  "figures/q1_decomp.png" 1785442631.48394 59575 8d569ed7e41980af8886d789db543c3f ""
  "figures/q1_err.png" 1785442451.78592 87442 a76f45152a369382e50e51be3a43a8b4 ""
  "figures/q2_control.png" 1785444485.00701 86524 7f96361ea78aee3bd132c181ac9beb49 ""
  "figures/q2_state_long.png" 1785444108.14416 114415 fdfcda3a1f0018cbfc084878f06a70a2 ""
  "figures/q2_state_short.png" 1785444816.41225 114490 b031409d3953343dd6523dcc8e76e67c ""
  "figures/q3_fb.png" 1785445507.94115 89068 708c4371c468689585e41a2d8302dbb1 ""
  "figures/q3_track.png" 1785445508.67331 90568 df28380e2e5d0bdf893f2bcc86cfc8c7 ""
  "report_fp.aux" 1785446014.23401 16602 d412ef4121943406e30a8bf1b8d49c04 "pdflatex"
  "report_fp.bbl" 1785446013.28952 527 d313001c3a51764d3baccde679617f56 "bibtex report_fp"
  "report_fp.out" 1785446014.10073 0 d41d8cd98f00b204e9800998ecf8427e "pdflatex"
  "report_fp.tex" 1785446011.93699 33004 3c54951de9aba1410ba8fcd1c30f3e06 ""
  (generated)
  "report_fp.aux"
  "report_fp.log"
  "report_fp.out"
  "report_fp.pdf"
  (rewritten before read)
```

# report/after-figureplacement/report_fp.fls

```
PWD /Users/vi/Documents/Projects/LMI_PROJECT/report/after-figureplacement
INPUT /usr/local/texlive/2025basic/texmf.cnf
INPUT /usr/local/texlive/2025basic/texmf-dist/web2c/texmf.cnf
INPUT /usr/local/texlive/2025basic/texmf-var/web2c/pdftex/pdflatex.fmt
INPUT /Users/vi/Documents/Projects/LMI_PROJECT/report/after-figureplacement/report_fp.tex
OUTPUT report_fp.log
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrartcl.cls
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrartcl.cls
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrkbase.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrkbase.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrbase.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrbase.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlfile.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlfile.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlfile-hook.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlfile-hook.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlogo.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlogo.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/keyval.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/keyval.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/tocbasic.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/tocbasic.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/typearea.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrsize10pt.clo
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrsize10pt.clo
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrsize10pt.clo
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/typearea.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/base/fontenc.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/base/fontenc.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/jknappen/ec/ecrm1000.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/base/inputenc.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/base/inputenc.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/babel/babel.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/babel/babel.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/babel/txtbabel.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/babel-english/english.ldf
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/babel-english/english.ldf
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/babel-english/english.ldf
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/babel/locale/en/babel-english.tex
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/babel/locale/en/babel-english.tex
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/babel/locale/en/babel-english.tex
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/babel/locale/en/babel-en.ini
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/geometry/geometry.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/geometry/geometry.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/iftex/ifvtex.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/iftex/ifvtex.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/iftex/iftex.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/iftex/iftex.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/multicol.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/multicol.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/graphicx.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/graphicx.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/graphics.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/graphics.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/trig.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/trig.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-cfg/graphics.cfg
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-cfg/graphics.cfg
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-cfg/graphics.cfg
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-def/pdftex.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-def/pdftex.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-def/pdftex.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/subcaption.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/subcaption.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/caption.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/caption.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/caption3.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/caption3.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/caption-koma.sto
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/caption-koma.sto
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/caption-koma.sto
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/wrapfig/wrapfig.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/wrapfig/wrapfig.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/mathtools/mathtools.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/mathtools/mathtools.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/calc.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/calc.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/mathtools/mhsetup.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/mathtools/mhsetup.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsmath.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsmath.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsopn.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amstext.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amstext.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsgen.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsgen.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsbsy.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsbsy.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsopn.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/bm.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/bm.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/booktabs/booktabs.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/booktabs/booktabs.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/tabularx.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/tabularx.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/array.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/array.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/ragged2e/ragged2e.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/ragged2e/ragged2e.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/psnfss/helvet.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/psnfss/helvet.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/microtype.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/microtype.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/etoolbox/etoolbox.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/etoolbox/etoolbox.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/microtype-pdftex.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/microtype-pdftex.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/microtype-pdftex.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/microtype.cfg
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/microtype.cfg
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/microtype.cfg
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/enumitem/enumitem.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/enumitem/enumitem.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/xcolor/xcolor.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/xcolor/xcolor.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-cfg/color.cfg
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-cfg/color.cfg
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-cfg/color.cfg
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/mathcolor.ltx
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/mathcolor.ltx
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/mathcolor.ltx
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/dvipsnam.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/dvipsnam.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/dvipsnam.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/hyperref.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/hyperref.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/kvsetkeys/kvsetkeys.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/kvsetkeys/kvsetkeys.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/kvdefinekeys/kvdefinekeys.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/kvdefinekeys/kvdefinekeys.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/pdfescape/pdfescape.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/pdfescape/pdfescape.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/ltxcmds/ltxcmds.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/ltxcmds/ltxcmds.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/pdftexcmds/pdftexcmds.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/pdftexcmds/pdftexcmds.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/infwarerr/infwarerr.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/infwarerr/infwarerr.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hycolor/hycolor.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hycolor/hycolor.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/nameref.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/nameref.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/refcount/refcount.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/refcount/refcount.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/gettitlestring/gettitlestring.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/gettitlestring/gettitlestring.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/kvoptions/kvoptions.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/kvoptions/kvoptions.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/stringenc/stringenc.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/stringenc/stringenc.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/pd1enc.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/pd1enc.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/pd1enc.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/intcalc/intcalc.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/intcalc/intcalc.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/puenc.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/puenc.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/puenc.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/url/url.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/url/url.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/bitset/bitset.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/bitset/bitset.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/bigintcalc/bigintcalc.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/bigintcalc/bigintcalc.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/atbegshi/atbegshi.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/base/atbegshi-ltx.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/base/atbegshi-ltx.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/hpdftex.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/hpdftex.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/hpdftex.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/atveryend/atveryend.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/base/atveryend-ltx.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/base/atveryend-ltx.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/rerunfilecheck/rerunfilecheck.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/rerunfilecheck/rerunfilecheck.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/uniquecounter/uniquecounter.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/generic/uniquecounter/uniquecounter.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/natbib/natbib.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/natbib/natbib.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlayer-scrpage.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlayer-scrpage.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlayer.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlayer.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/psnfss/t1phv.fd
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/psnfss/t1phv.fd
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/psnfss/t1phv.fd
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8t.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/l3backend/l3backend-pdftex.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/l3backend/l3backend-pdftex.def
INPUT ./report_fp.aux
INPUT ./report_fp.aux
INPUT report_fp.aux
OUTPUT report_fp.aux
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/context/base/mkii/supp-pdf.mkii
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/context/base/mkii/supp-pdf.mkii
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/context/base/mkii/supp-pdf.mkii
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/epstopdf-pkg/epstopdf-base.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/epstopdf-pkg/epstopdf-base.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/grfext/grfext.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/grfext/grfext.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/latexconfig/epstopdf-sys.cfg
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/latexconfig/epstopdf-sys.cfg
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/latexconfig/epstopdf-sys.cfg
INPUT ./report_fp.out
INPUT ./report_fp.out
INPUT report_fp.out
INPUT report_fp.out
INPUT ./report_fp.out
INPUT ./report_fp.out
OUTPUT report_fp.out
OUTPUT report_fp.pdf
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/bookmark/bookmark.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/bookmark/bookmark.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/bookmark/bookmark.sty
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/bookmark/bkm-pdftex.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/bookmark/bkm-pdftex.def
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/bookmark/bkm-pdftex.def
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8t.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvb8t.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8t.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmr9.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/mt-cmr.cfg
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/mt-cmr.cfg
INPUT /usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/mt-cmr.cfg
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmr6.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmmi9.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmmi6.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmsy9.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmsy6.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmex9.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmex7.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmex7.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbx9.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbx6.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbx5.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmmib9.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmmib6.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmmib5.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmbsy9.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmbsy6.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmbsy5.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvb8t.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvb8t.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmex7.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbx10.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbx7.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmmib10.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmmib7.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbsy10.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmbsy7.tfm
INPUT ./figures/find_u2.png
INPUT ./figures/find_u2.png
INPUT ./figures/find_u2.png
INPUT ./figures/find_u2.png
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8t.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8t.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvro8t.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/vf/adobe/helvetic/phvb8t.vf
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvb8r.tfm
INPUT /usr/local/texlive/2025basic/texmf-var/fonts/map/pdftex/updmap/pdftex.map
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/enc/dvips/base/8r.enc
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/vf/adobe/helvetic/phvr8t.vf
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8r.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/vf/adobe/helvetic/phvb8t.vf
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvb8r.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/vf/adobe/helvetic/phvr8t.vf
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8r.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/vf/adobe/helvetic/phvb8t.vf
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvb8r.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/vf/adobe/helvetic/phvro8t.vf
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvro8r.tfm
INPUT ./figures/q1_cmp.png
INPUT ./figures/q1_cmp.png
INPUT ./figures/q1_cmp.png
INPUT ./figures/q1_cmp.png
INPUT ./figures/q1_decomp.png
INPUT ./figures/q1_decomp.png
INPUT ./figures/q1_decomp.png
INPUT ./figures/q1_decomp.png
INPUT ./figures/q1_err.png
INPUT ./figures/q1_err.png
INPUT ./figures/q1_err.png
INPUT ./figures/q1_err.png
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmss10.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmss8.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmss8.tfm
INPUT ./figures/q2_control.png
INPUT ./figures/q2_control.png
INPUT ./figures/q2_control.png
INPUT ./figures/q2_control.png
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/vf/adobe/helvetic/phvr8t.vf
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8r.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/vf/adobe/helvetic/phvr8t.vf
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8r.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmss9.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmss8.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8t.tfm
INPUT ./figures/q2_state_short.png
INPUT ./figures/q2_state_short.png
INPUT ./figures/q2_state_short.png
INPUT ./figures/q2_state_short.png
INPUT ./figures/q2_state_long.png
INPUT ./figures/q2_state_long.png
INPUT ./figures/q2_state_long.png
INPUT ./figures/q2_state_long.png
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/vf/adobe/helvetic/phvr8t.vf
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8r.tfm
INPUT ./figures/q3_track.png
INPUT ./figures/q3_track.png
INPUT ./figures/q3_track.png
INPUT ./figures/q3_track.png
INPUT ./figures/q3_fb.png
INPUT ./figures/q3_fb.png
INPUT ./figures/q3_fb.png
INPUT ./figures/q3_fb.png
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8t.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmr8.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmmi8.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmsy8.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmex8.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmss8.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbx8.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmmib8.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/amsfonts/cmextra/cmbsy8.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8t.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmr12.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmmi12.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmsy10.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmex10.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmss12.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbx12.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmmib10.tfm
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/public/cm/cmbsy10.tfm
INPUT ./report_fp.bbl
INPUT ./report_fp.bbl
INPUT report_fp.bbl
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/vf/adobe/helvetic/phvr8t.vf
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/tfm/adobe/helvetic/phvr8r.tfm
INPUT report_fp.aux
INPUT ./report_fp.out
INPUT ./report_fp.out
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmbx10.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmex10.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi10.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi12.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi6.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi7.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi8.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi9.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmib10.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cmextra/cmmib7.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cmextra/cmmib9.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr10.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr5.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr6.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr7.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr8.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr9.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmss10.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmsy10.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmsy7.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmsy8.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmsy9.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/urw/helvetic/uhvb8a.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/urw/helvetic/uhvr8a.pfb
INPUT /usr/local/texlive/2025basic/texmf-dist/fonts/type1/urw/helvetic/uhvro8a.pfb
```

# report/after-figureplacement/report_fp.log

```
This is pdfTeX, Version 3.141592653-2.6-1.40.27 (TeX Live 2025) (preloaded format=pdflatex 2025.7.18)  30 JUL 2026 23:13
entering extended mode
 restricted \write18 enabled.
 file:line:error style messages enabled.
 %&-line parsing enabled.
**/Users/vi/Documents/Projects/LMI_PROJECT/report/after-figureplacement/report_fp.tex
(/Users/vi/Documents/Projects/LMI_PROJECT/report/after-figureplacement/report_fp.tex
LaTeX2e <2024-11-01> patch level 2
L3 programming layer <2025-01-18>
(/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrartcl.cls
Document Class: scrartcl 2024/10/24 v3.43 KOMA-Script document class (article)
(/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrkbase.sty
Package: scrkbase 2024/10/24 v3.43 KOMA-Script package (KOMA-Script-dependent basics and keyval usage)
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrbase.sty
Package: scrbase 2024/10/24 v3.43 KOMA-Script package (KOMA-Script-independent basics and keyval usage)
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlfile.sty
Package: scrlfile 2024/10/24 v3.43 KOMA-Script package (file load hooks)
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlfile-hook.sty
Package: scrlfile-hook 2024/10/24 v3.43 KOMA-Script package (using LaTeX hooks)
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlogo.sty
Package: scrlogo 2024/10/24 v3.43 KOMA-Script package (logo)
))) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/keyval.sty
Package: keyval 2022/05/29 v1.15 key=value parser (DPC)
\KV@toks@=\toks17
)
Applying: [2021/05/01] Usage of raw or classic option list on input line 254.
Already applied: [0000/00/00] Usage of raw or classic option list on input line 370.
)) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/tocbasic.sty
Package: tocbasic 2024/10/24 v3.43 KOMA-Script package (handling toc-files)
\scr@dte@tocline@numberwidth=\skip49
\scr@dte@tocline@numbox=\box52
)
Package tocbasic Info: omitting babel extension for `toc'
(tocbasic)             because of feature `nobabel' available
(tocbasic)             for `toc' on input line 135.
Package scrartcl Info: You've used standard option `10pt'.
(scrartcl)             This is correct!
(scrartcl)             Internally I'm using `fontsize=10pt'.
(scrartcl)             If you'd like to set the option with \KOMAoptions,
(scrartcl)             you'd have to use `fontsize=10pt' there
(scrartcl)             instead of `10pt', too.
Class scrartcl Info: File `scrsize10pt.clo' used to setup font sizes on input line 2466.
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrsize10pt.clo
File: scrsize10pt.clo 2024/10/24 v3.43 KOMA-Script font size class option (10pt)
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/typearea.sty
Package: typearea 2024/10/24 v3.43 KOMA-Script package (type area)
\ta@bcor=\skip50
\ta@div=\count196
Package typearea Info: You've used standard option `a4paper'.
(typearea)             This is correct!
(typearea)             Internally I'm using `paper=a4'.
(typearea)             If you'd like to set the option with \KOMAoptions,
(typearea)             you'd have to use `paper=a4' there
(typearea)             instead of `a4paper', too.
\ta@hblk=\skip51
\ta@vblk=\skip52
\ta@temp=\skip53
\footheight=\skip54
Package typearea Info: These are the values describing the layout:
(typearea)             DIV  = 8
(typearea)             BCOR = 0.0pt
(typearea)             \paperwidth      = 597.50793pt
(typearea)              \textwidth      = 373.44246pt
(typearea)              DIV departure   = -4%
(typearea)              \evensidemargin = 39.76274pt
(typearea)              \oddsidemargin  = 39.76274pt
(typearea)             \paperheight     = 845.04694pt
(typearea)              \textheight     = 538.0pt
(typearea)              \topmargin      = 0.36087pt
(typearea)              \headheight     = 15.0pt
(typearea)              \headsep        = 18.0pt
(typearea)              \topskip        = 10.0pt
(typearea)              \footskip       = 42.0pt
(typearea)              \baselineskip   = 12.0pt
(typearea)              on input line 1802.
)
\c@part=\count197
\c@section=\count198
\c@subsection=\count199
\c@subsubsection=\count266
\c@paragraph=\count267
\c@subparagraph=\count268
\scr@dte@section@maxnumwidth=\skip55
Class scrartcl Info: using compatibility default `runin=bysign'
(scrartcl)           for `\section on input line 5107.
Class scrartcl Info: using compatibility default `afterindent=bysign'
(scrartcl)           for `\section on input line 5107.
\scr@dte@part@maxnumwidth=\skip56
Class scrartcl Info: using compatibility default `afterindent=false'
(scrartcl)           for `\part on input line 5115.
\scr@dte@subsection@maxnumwidth=\skip57
Class scrartcl Info: using compatibility default `runin=bysign'
(scrartcl)           for `\subsection on input line 5125.
Class scrartcl Info: using compatibility default `afterindent=bysign'
(scrartcl)           for `\subsection on input line 5125.
\scr@dte@subsubsection@maxnumwidth=\skip58
Class scrartcl Info: using compatibility default `runin=bysign'
(scrartcl)           for `\subsubsection on input line 5135.
Class scrartcl Info: using compatibility default `afterindent=bysign'
(scrartcl)           for `\subsubsection on input line 5135.
\scr@dte@paragraph@maxnumwidth=\skip59
Class scrartcl Info: using compatibility default `runin=bysign'
(scrartcl)           for `\paragraph on input line 5146.
Class scrartcl Info: using compatibility default `afterindent=bysign'
(scrartcl)           for `\paragraph on input line 5146.
\scr@dte@subparagraph@maxnumwidth=\skip60
Class scrartcl Info: using compatibility default `runin=bysign'
(scrartcl)           for `\subparagraph on input line 5156.
Class scrartcl Info: using compatibility default `afterindent=bysign'
(scrartcl)           for `\subparagraph on input line 5156.
\abovecaptionskip=\skip61
\belowcaptionskip=\skip62
\c@pti@nb@sid@b@x=\box53
Package tocbasic Info: omitting babel extension for `lof'
(tocbasic)             because of feature `nobabel' available
(tocbasic)             for `lof' on input line 6337.
\scr@dte@figure@maxnumwidth=\skip63
\c@figure=\count269
Package tocbasic Info: omitting babel extension for `lot'
(tocbasic)             because of feature `nobabel' available
(tocbasic)             for `lot' on input line 6353.
\scr@dte@table@maxnumwidth=\skip64
\c@table=\count270
Class scrartcl Info: Redefining `\numberline' on input line 6523.
\bibindent=\dimen141
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/base/fontenc.sty
Package: fontenc 2021/04/29 v2.0v Standard LaTeX package
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/base/inputenc.sty
Package: inputenc 2024/02/08 v1.3d Input encoding file
\inpenc@prehook=\toks18
\inpenc@posthook=\toks19
) (/usr/local/texlive/2025basic/texmf-dist/tex/generic/babel/babel.sty
Package: babel 2025/02/14 v25.4 The multilingual framework for pdfLaTeX, LuaLaTeX and XeLaTeX
\babel@savecnt=\count271
\U@D=\dimen142
\l@unhyphenated=\language23
 (/usr/local/texlive/2025basic/texmf-dist/tex/generic/babel/txtbabel.def)
\bbl@readstream=\read2
\bbl@dirlevel=\count272
 (/usr/local/texlive/2025basic/texmf-dist/tex/generic/babel-english/english.ldf
Language: english 2017/06/06 v3.3r English support from the babel system
Package babel Info: Hyphen rules for 'canadian' set to \l@english
(babel)             (\language0). Reported on input line 102.
Package babel Info: Hyphen rules for 'australian' set to \l@ukenglish
(babel)             (\language7). Reported on input line 105.
Package babel Info: Hyphen rules for 'newzealand' set to \l@ukenglish
(babel)             (\language7). Reported on input line 108.
)) (/usr/local/texlive/2025basic/texmf-dist/tex/generic/babel/locale/en/babel-english.tex
Package babel Info: Importing font and identification data for english
(babel)             from babel-en.ini. Reported on input line 11.
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/geometry/geometry.sty
Package: geometry 2020/01/02 v5.9 Page Geometry
 (/usr/local/texlive/2025basic/texmf-dist/tex/generic/iftex/ifvtex.sty
Package: ifvtex 2019/10/25 v1.7 ifvtex legacy package. Use iftex instead.
 (/usr/local/texlive/2025basic/texmf-dist/tex/generic/iftex/iftex.sty
Package: iftex 2024/12/12 v1.0g TeX engine tests
))
\Gm@cnth=\count273
\Gm@cntv=\count274
\c@Gm@tempcnt=\count275
\Gm@bindingoffset=\dimen143
\Gm@wd@mp=\dimen144
\Gm@odd@mp=\dimen145
\Gm@even@mp=\dimen146
\Gm@layoutwidth=\dimen147
\Gm@layoutheight=\dimen148
\Gm@layouthoffset=\dimen149
\Gm@layoutvoffset=\dimen150
\Gm@dimlist=\toks20
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/multicol.sty
Package: multicol 2024/09/14 v1.9i multicolumn formatting (FMi)
\c@tracingmulticols=\count276
\mult@box=\box54
\multicol@leftmargin=\dimen151
\c@unbalance=\count277
\c@collectmore=\count278
\doublecol@number=\count279
\multicoltolerance=\count280
\multicolpretolerance=\count281
\full@width=\dimen152
\page@free=\dimen153
\premulticols=\dimen154
\postmulticols=\dimen155
\multicolsep=\skip65
\multicolbaselineskip=\skip66
\partial@page=\box55
\last@line=\box56
\mc@boxedresult=\box57
\maxbalancingoverflow=\dimen156
\mult@rightbox=\box58
\mult@grightbox=\box59
\mult@firstbox=\box60
\mult@gfirstbox=\box61
\@tempa=\box62
\@tempa=\box63
\@tempa=\box64
\@tempa=\box65
\@tempa=\box66
\@tempa=\box67
\@tempa=\box68
\@tempa=\box69
\@tempa=\box70
\@tempa=\box71
\@tempa=\box72
\@tempa=\box73
\@tempa=\box74
\@tempa=\box75
\@tempa=\box76
\@tempa=\box77
\@tempa=\box78
\@tempa=\box79
\@tempa=\box80
\@tempa=\box81
\@tempa=\box82
\@tempa=\box83
\@tempa=\box84
\@tempa=\box85
\@tempa=\box86
\@tempa=\box87
\@tempa=\box88
\@tempa=\box89
\@tempa=\box90
\@tempa=\box91
\@tempa=\box92
\@tempa=\box93
\@tempa=\box94
\@tempa=\box95
\@tempa=\box96
\@tempa=\box97
\c@minrows=\count282
\c@columnbadness=\count283
\c@finalcolumnbadness=\count284
\last@try=\dimen157
\multicolovershoot=\dimen158
\multicolundershoot=\dimen159
\mult@nat@firstbox=\box98
\colbreak@box=\box99
\mc@col@check@num=\count285
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/graphicx.sty
Package: graphicx 2021/09/16 v1.2d Enhanced LaTeX Graphics (DPC,SPQR)
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/graphics.sty
Package: graphics 2024/08/06 v1.4g Standard LaTeX Graphics (DPC,SPQR)
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/trig.sty
Package: trig 2023/12/02 v1.11 sin cos tan (DPC)
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-cfg/graphics.cfg
File: graphics.cfg 2016/06/04 v1.11 sample graphics configuration
)
Package graphics Info: Driver file: pdftex.def on input line 106.
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-def/pdftex.def
File: pdftex.def 2024/04/13 v1.2c Graphics/color driver for pdftex
))
\Gin@req@height=\dimen160
\Gin@req@width=\dimen161
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/subcaption.sty
Package: subcaption 2023/07/28 v1.6b Sub-captions (AR)
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/caption.sty
Package: caption 2023/08/05 v3.6o Customizing captions (AR)
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/caption3.sty
Package: caption3 2023/07/31 v2.4d caption3 kernel (AR)
\caption@tempdima=\dimen162
\captionmargin=\dimen163
\caption@leftmargin=\dimen164
\caption@rightmargin=\dimen165
\caption@width=\dimen166
\caption@indent=\dimen167
\caption@parindent=\dimen168
\caption@hangindent=\dimen169
Package caption Info: KOMA-Script document class detected.
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/caption/caption-koma.sto
File: caption-koma.sto 2023/09/08 v2.0e Adaption of the caption package to the KOMA-Script document classes (AR)
))
\c@caption@flags=\count286
\c@continuedfloat=\count287
Package caption Info: KOMA-Script scrextend package detected.
\caption@addmargin@hsize=\dimen170
\caption@addmargin@linewidth=\dimen171
)
Package caption Info: New subtype `subfigure' on input line 238.
\c@subfigure=\count288
Package caption Info: New subtype `subtable' on input line 238.
\c@subtable=\count289
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/wrapfig/wrapfig.sty
\wrapoverhang=\dimen172
\WF@size=\dimen173
\c@WF@wrappedlines=\count290
\WF@box=\box100
\WF@everypar=\toks21
Package: wrapfig 2003/01/31  v 3.6
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/mathtools/mathtools.sty
Package: mathtools 2024/10/04 v1.31 mathematical typesetting tools
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/calc.sty
Package: calc 2023/07/08 v4.3 Infix arithmetic (KKT,FJ)
\calc@Acount=\count291
\calc@Bcount=\count292
\calc@Adimen=\dimen174
\calc@Bdimen=\dimen175
\calc@Askip=\skip67
\calc@Bskip=\skip68
LaTeX Info: Redefining \setlength on input line 80.
LaTeX Info: Redefining \addtolength on input line 81.
\calc@Ccount=\count293
\calc@Cskip=\skip69
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/mathtools/mhsetup.sty
Package: mhsetup 2021/03/18 v1.4 programming setup (MH)
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsmath.sty
Package: amsmath 2024/11/05 v2.17t AMS math features
\@mathmargin=\skip70

For additional information on amsmath, use the `?' option.
(/usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amstext.sty
Package: amstext 2021/08/26 v2.01 AMS text
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsgen.sty
File: amsgen.sty 1999/11/30 v2.0 generic functions
\@emptytoks=\toks22
\ex@=\dimen176
)) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsbsy.sty
Package: amsbsy 1999/11/29 v1.2d Bold Symbols
\pmbraise@=\dimen177
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/amsmath/amsopn.sty
Package: amsopn 2022/04/08 v2.04 operator names
)
\inf@bad=\count294
LaTeX Info: Redefining \frac on input line 233.
\uproot@=\count295
\leftroot@=\count296
LaTeX Info: Redefining \overline on input line 398.
LaTeX Info: Redefining \colon on input line 409.
\classnum@=\count297
\DOTSCASE@=\count298
LaTeX Info: Redefining \ldots on input line 495.
LaTeX Info: Redefining \dots on input line 498.
LaTeX Info: Redefining \cdots on input line 619.
\Mathstrutbox@=\box101
\strutbox@=\box102
LaTeX Info: Redefining \big on input line 721.
LaTeX Info: Redefining \Big on input line 722.
LaTeX Info: Redefining \bigg on input line 723.
LaTeX Info: Redefining \Bigg on input line 724.
\big@size=\dimen178
LaTeX Font Info:    Redeclaring font encoding OML on input line 742.
LaTeX Font Info:    Redeclaring font encoding OMS on input line 743.
\macc@depth=\count299
LaTeX Info: Redefining \bmod on input line 904.
LaTeX Info: Redefining \pmod on input line 909.
LaTeX Info: Redefining \smash on input line 939.
LaTeX Info: Redefining \relbar on input line 969.
LaTeX Info: Redefining \Relbar on input line 970.
\c@MaxMatrixCols=\count300
\dotsspace@=\muskip17
\c@parentequation=\count301
\dspbrk@lvl=\count302
\tag@help=\toks23
\row@=\count303
\column@=\count304
\maxfields@=\count305
\andhelp@=\toks24
\eqnshift@=\dimen179
\alignsep@=\dimen180
\tagshift@=\dimen181
\tagwidth@=\dimen182
\totwidth@=\dimen183
\lineht@=\dimen184
\@envbody=\toks25
\multlinegap=\skip71
\multlinetaggap=\skip72
\mathdisplay@stack=\toks26
LaTeX Info: Redefining \[ on input line 2953.
LaTeX Info: Redefining \] on input line 2954.
)
\g_MT_multlinerow_int=\count306
\l_MT_multwidth_dim=\dimen185
\origjot=\skip73
\l_MT_shortvdotswithinadjustabove_dim=\dimen186
\l_MT_shortvdotswithinadjustbelow_dim=\dimen187
\l_MT_above_intertext_sep=\dimen188
\l_MT_below_intertext_sep=\dimen189
\l_MT_above_shortintertext_sep=\dimen190
\l_MT_below_shortintertext_sep=\dimen191
\xmathstrut@box=\box103
\xmathstrut@dim=\dimen192
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/bm.sty
Package: bm 2023/12/19 v1.2f Bold Symbol Support (DPC/FMi)
\symboldoperators=\mathgroup4
\symboldletters=\mathgroup5
\symboldsymbols=\mathgroup6
Package bm Info: No bold for \OMX/cmex/m/n, using \pmb.
LaTeX Font Info:    Redeclaring math alphabet \mathbf on input line 149.
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/booktabs/booktabs.sty
Package: booktabs 2020/01/12 v1.61803398 Publication quality tables
\heavyrulewidth=\dimen193
\lightrulewidth=\dimen194
\cmidrulewidth=\dimen195
\belowrulesep=\dimen196
\belowbottomsep=\dimen197
\aboverulesep=\dimen198
\abovetopsep=\dimen199
\cmidrulesep=\dimen256
\cmidrulekern=\dimen257
\defaultaddspace=\dimen258
\@cmidla=\count307
\@cmidlb=\count308
\@aboverulesep=\dimen259
\@belowrulesep=\dimen260
\@thisruleclass=\count309
\@lastruleclass=\count310
\@thisrulewidth=\dimen261
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/tabularx.sty
Package: tabularx 2023/12/11 v2.12a `tabularx' package (DPC)
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/tools/array.sty
Package: array 2024/10/17 v2.6g Tabular extension package (FMi)
\col@sep=\dimen262
\ar@mcellbox=\box104
\extrarowheight=\dimen263
\NC@list=\toks27
\extratabsurround=\skip74
\backup@length=\skip75
\ar@cellbox=\box105
)
\TX@col@width=\dimen264
\TX@old@table=\dimen265
\TX@old@col=\dimen266
\TX@target=\dimen267
\TX@delta=\dimen268
\TX@cols=\count311
\TX@ftn=\toks28
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/ragged2e/ragged2e.sty
Package: ragged2e 2023/06/22 v3.6 ragged2e Package
\CenteringLeftskip=\skip76
\RaggedLeftLeftskip=\skip77
\RaggedRightLeftskip=\skip78
\CenteringRightskip=\skip79
\RaggedLeftRightskip=\skip80
\RaggedRightRightskip=\skip81
\CenteringParfillskip=\skip82
\RaggedLeftParfillskip=\skip83
\RaggedRightParfillskip=\skip84
\JustifyingParfillskip=\skip85
\CenteringParindent=\skip86
\RaggedLeftParindent=\skip87
\RaggedRightParindent=\skip88
\JustifyingParindent=\skip89
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/psnfss/helvet.sty
Package: helvet 2020/03/25 PSNFSS-v9.3 (WaS) 
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/microtype.sty
Package: microtype 2025/02/11 v3.2a Micro-typographical refinements (RS)
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/etoolbox/etoolbox.sty
Package: etoolbox 2025/02/11 v2.5l e-TeX tools for LaTeX (JAW)
\etb@tempcnta=\count312
)
\MT@toks=\toks29
\MT@tempbox=\box106
\MT@count=\count313
LaTeX Info: Redefining \noprotrusionifhmode on input line 1087.
LaTeX Info: Redefining \leftprotrusion on input line 1088.
\MT@prot@toks=\toks30
LaTeX Info: Redefining \rightprotrusion on input line 1107.
LaTeX Info: Redefining \textls on input line 1449.
\MT@outer@kern=\dimen269
LaTeX Info: Redefining \microtypecontext on input line 2053.
LaTeX Info: Redefining \textmicrotypecontext on input line 2070.
\MT@listname@count=\count314
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/microtype-pdftex.def
File: microtype-pdftex.def 2025/02/11 v3.2a Definitions specific to pdftex (RS)
LaTeX Info: Redefining \lsstyle on input line 944.
LaTeX Info: Redefining \lslig on input line 944.
\MT@outer@space=\skip90
)
Package microtype Info: Loading configuration file microtype.cfg.
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/microtype.cfg
File: microtype.cfg 2025/02/11 v3.2a microtype main configuration file (RS)
)
LaTeX Info: Redefining \microtypesetup on input line 3065.
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/enumitem/enumitem.sty
Package: enumitem 2025/02/06 v3.11 Customized lists
\labelindent=\skip91
\enit@outerparindent=\dimen270
\enit@toks=\toks31
\enit@inbox=\box107
\enit@count@id=\count315
\enitdp@description=\count316
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/xcolor/xcolor.sty
Package: xcolor 2024/09/29 v3.02 LaTeX color extensions (UK)
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics-cfg/color.cfg
File: color.cfg 2016/01/02 v1.6 sample color configuration
)
Package xcolor Info: Driver file: pdftex.def on input line 274.
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/mathcolor.ltx)
Package xcolor Info: Model `cmy' substituted by `cmy0' on input line 1349.
Package xcolor Info: Model `hsb' substituted by `rgb' on input line 1353.
Package xcolor Info: Model `RGB' extended on input line 1365.
Package xcolor Info: Model `HTML' substituted by `rgb' on input line 1367.
Package xcolor Info: Model `Hsb' substituted by `hsb' on input line 1368.
Package xcolor Info: Model `tHsb' substituted by `hsb' on input line 1369.
Package xcolor Info: Model `HSB' substituted by `hsb' on input line 1370.
Package xcolor Info: Model `Gray' substituted by `gray' on input line 1371.
Package xcolor Info: Model `wave' substituted by `hsb' on input line 1372.
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/graphics/dvipsnam.def
File: dvipsnam.def 2016/06/17 v3.0m Driver-dependent file (DPC,SPQR)
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/hyperref.sty
Package: hyperref 2024-11-05 v7.01l Hypertext links for LaTeX
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/kvsetkeys/kvsetkeys.sty
Package: kvsetkeys 2022-10-05 v1.19 Key value parser (HO)
) (/usr/local/texlive/2025basic/texmf-dist/tex/generic/kvdefinekeys/kvdefinekeys.sty
Package: kvdefinekeys 2019-12-19 v1.6 Define keys (HO)
) (/usr/local/texlive/2025basic/texmf-dist/tex/generic/pdfescape/pdfescape.sty
Package: pdfescape 2019/12/09 v1.15 Implements pdfTeX's escape features (HO)
 (/usr/local/texlive/2025basic/texmf-dist/tex/generic/ltxcmds/ltxcmds.sty
Package: ltxcmds 2023-12-04 v1.26 LaTeX kernel commands for general use (HO)
) (/usr/local/texlive/2025basic/texmf-dist/tex/generic/pdftexcmds/pdftexcmds.sty
Package: pdftexcmds 2020-06-27 v0.33 Utility functions of pdfTeX for LuaTeX (HO)
 (/usr/local/texlive/2025basic/texmf-dist/tex/generic/infwarerr/infwarerr.sty
Package: infwarerr 2019/12/03 v1.5 Providing info/warning/error messages (HO)
)
Package pdftexcmds Info: \pdf@primitive is available.
Package pdftexcmds Info: \pdf@ifprimitive is available.
Package pdftexcmds Info: \pdfdraftmode found.
)) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/hycolor/hycolor.sty
Package: hycolor 2020-01-27 v1.10 Color options for hyperref/bookmark (HO)
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/nameref.sty
Package: nameref 2023-11-26 v2.56 Cross-referencing by name of section
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/refcount/refcount.sty
Package: refcount 2019/12/15 v3.6 Data extraction from label references (HO)
) (/usr/local/texlive/2025basic/texmf-dist/tex/generic/gettitlestring/gettitlestring.sty
Package: gettitlestring 2019/12/15 v1.6 Cleanup title references (HO)
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/kvoptions/kvoptions.sty
Package: kvoptions 2022-06-15 v3.15 Key value format for package options (HO)
))
\c@section@level=\count317
) (/usr/local/texlive/2025basic/texmf-dist/tex/generic/stringenc/stringenc.sty
Package: stringenc 2019/11/29 v1.12 Convert strings between diff. encodings (HO)
)
\@linkdim=\dimen271
\Hy@linkcounter=\count318
\Hy@pagecounter=\count319
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/pd1enc.def
File: pd1enc.def 2024-11-05 v7.01l Hyperref: PDFDocEncoding definition (HO)
Now handling font encoding PD1 ...
... no UTF-8 mapping file for font encoding PD1
) (/usr/local/texlive/2025basic/texmf-dist/tex/generic/intcalc/intcalc.sty
Package: intcalc 2019/12/15 v1.3 Expandable calculations with integers (HO)
)
\Hy@SavedSpaceFactor=\count320
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/puenc.def
File: puenc.def 2024-11-05 v7.01l Hyperref: PDF Unicode definition (HO)
Now handling font encoding PU ...
... no UTF-8 mapping file for font encoding PU
)
Package hyperref Info: Hyper figures OFF on input line 4157.
Package hyperref Info: Link nesting OFF on input line 4162.
Package hyperref Info: Hyper index ON on input line 4165.
Package hyperref Info: Plain pages OFF on input line 4172.
Package hyperref Info: Backreferencing OFF on input line 4177.
Package hyperref Info: Implicit mode ON; LaTeX internals redefined.
Package hyperref Info: Bookmarks ON on input line 4424.
\c@Hy@tempcnt=\count321
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/url/url.sty
\Urlmuskip=\muskip18
Package: url 2013/09/16  ver 3.4  Verb mode for urls, etc.
)
LaTeX Info: Redefining \url on input line 4763.
\XeTeXLinkMargin=\dimen272
 (/usr/local/texlive/2025basic/texmf-dist/tex/generic/bitset/bitset.sty
Package: bitset 2019/12/09 v1.3 Handle bit-vector datatype (HO)
 (/usr/local/texlive/2025basic/texmf-dist/tex/generic/bigintcalc/bigintcalc.sty
Package: bigintcalc 2019/12/15 v1.5 Expandable calculations on big integers (HO)
))
\Fld@menulength=\count322
\Field@Width=\dimen273
\Fld@charsize=\dimen274
Package hyperref Info: Hyper figures OFF on input line 6042.
Package hyperref Info: Link nesting OFF on input line 6047.
Package hyperref Info: Hyper index ON on input line 6050.
Package hyperref Info: backreferencing OFF on input line 6057.
Package hyperref Info: Link coloring OFF on input line 6062.
Package hyperref Info: Link coloring with OCG OFF on input line 6067.
Package hyperref Info: PDF/A mode OFF on input line 6072.
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/base/atbegshi-ltx.sty
Package: atbegshi-ltx 2021/01/10 v1.0c Emulation of the original atbegshi
package with kernel methods
)
\Hy@abspage=\count323
\c@Item=\count324
\c@Hfootnote=\count325
)
Package hyperref Info: Driver (autodetected): hpdftex.
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/hyperref/hpdftex.def
File: hpdftex.def 2024-11-05 v7.01l Hyperref driver for pdfTeX
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/base/atveryend-ltx.sty
Package: atveryend-ltx 2020/08/19 v1.0a Emulation of the original atveryend package
with kernel methods
)
\Fld@listcount=\count326
\c@bookmark@seq@number=\count327
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/rerunfilecheck/rerunfilecheck.sty
Package: rerunfilecheck 2022-07-10 v1.10 Rerun checks for auxiliary files (HO)
 (/usr/local/texlive/2025basic/texmf-dist/tex/generic/uniquecounter/uniquecounter.sty
Package: uniquecounter 2019/12/15 v1.4 Provide unlimited unique counter (HO)
)
Package uniquecounter Info: New unique counter `rerunfilecheck' on input line 285.
)
\Hy@SectionHShift=\skip92
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/natbib/natbib.sty
Package: natbib 2010/09/13 8.31b (PWD, AO)
\bibhang=\skip93
\bibsep=\skip94
LaTeX Info: Redefining \cite on input line 694.
\c@NAT@ctr=\count328
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlayer-scrpage.sty
Package: scrlayer-scrpage 2024/10/24 v3.43 KOMA-Script package (end user interface for scrlayer)
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/koma-script/scrlayer.sty
Package: scrlayer 2024/10/24 v3.43 KOMA-Script package (defining layers and page styles)
Package scrlayer Info: patching LaTeX kernel macro \pagestyle on input line 2182.
)
Package scrlayer-scrpage Info: Makeing stand-alone element `pagehead' from
(scrlayer-scrpage)             alias to `pageheadfoot' on input line 1106.
)
LaTeX Font Info:    Trying to load font information for T1+phv on input line 164.
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/psnfss/t1phv.fd
File: t1phv.fd 2020/03/25 scalable font definitions for T1/phv.
)
LaTeX Font Info:    Font shape `T1/phv/m/n' will be
(Font)              scaled to size 9.49997pt on input line 164.
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/l3backend/l3backend-pdftex.def
File: l3backend-pdftex.def 2024-05-08 L3 backend support: PDF output (pdfTeX)
\l__color_backend_stack_int=\count329
\l__pdf_internal_box=\box108
) (./report_fp.aux)
\openout1 = `report_fp.aux'.

LaTeX Font Info:    Checking defaults for OML/cmm/m/it on input line 164.
LaTeX Font Info:    ... okay on input line 164.
LaTeX Font Info:    Checking defaults for OMS/cmsy/m/n on input line 164.
LaTeX Font Info:    ... okay on input line 164.
LaTeX Font Info:    Checking defaults for OT1/cmr/m/n on input line 164.
LaTeX Font Info:    ... okay on input line 164.
LaTeX Font Info:    Checking defaults for T1/cmr/m/n on input line 164.
LaTeX Font Info:    ... okay on input line 164.
LaTeX Font Info:    Checking defaults for TS1/cmr/m/n on input line 164.
LaTeX Font Info:    ... okay on input line 164.
LaTeX Font Info:    Checking defaults for OMX/cmex/m/n on input line 164.
LaTeX Font Info:    ... okay on input line 164.
LaTeX Font Info:    Checking defaults for U/cmr/m/n on input line 164.
LaTeX Font Info:    ... okay on input line 164.
LaTeX Font Info:    Checking defaults for PD1/pdf/m/n on input line 164.
LaTeX Font Info:    ... okay on input line 164.
LaTeX Font Info:    Checking defaults for PU/pdf/m/n on input line 164.
LaTeX Font Info:    ... okay on input line 164.
Package scrbase Info: activating english \abstractname on input line 164.
Package scrbase Info: activating english \partname on input line 164.
Package scrbase Info: activating english \figurename on input line 164.
Package scrbase Info: activating english \tablename on input line 164.
Package scrbase Info: activating english \contentsname on input line 164.
Package scrbase Info: activating english \listfigurename on input line 164.
Package scrbase Info: activating english \listtablename on input line 164.
Package scrbase Info: activating english \refname on input line 164.
Package scrbase Info: activating english \indexname on input line 164.

*geometry* driver: auto-detecting
*geometry* detected driver: pdftex
*geometry* verbose mode - [ preamble ] result:
* driver: pdftex
* paper: a4paper
* layout: <same size as paper>
* layoutoffset:(h,v)=(0.0pt,0.0pt)
* modes: includefoot 
* h-part:(L,W,R)=(56.9055pt, 483.69687pt, 56.9055pt)
* v-part:(T,H,B)=(28.45274pt, 788.14136pt, 28.45274pt)
* \paperwidth=597.50787pt
* \paperheight=845.04684pt
* \textwidth=483.69687pt
* \textheight=746.14136pt
* \oddsidemargin=-15.36449pt
* \evensidemargin=-15.36449pt
* \topmargin=-76.81725pt
* \headheight=15.0pt
* \headsep=18.0pt
* \topskip=10.0pt
* \footskip=42.67912pt
* \marginparwidth=74.68849pt
* \marginparsep=12.8401pt
* \columnsep=22.76219pt
* \skip\footins=9.0pt plus 4.0pt minus 2.0pt
* \hoffset=0.0pt
* \voffset=0.0pt
* \mag=1000
* \@twocolumnfalse
* \@twosidefalse
* \@mparswitchfalse
* \@reversemarginfalse
* (1in=72.27pt=25.4mm, 1cm=28.453pt)

(/usr/local/texlive/2025basic/texmf-dist/tex/context/base/mkii/supp-pdf.mkii
[Loading MPS to PDF converter (version 2006.09.02).]
\scratchcounter=\count330
\scratchdimen=\dimen275
\scratchbox=\box109
\nofMPsegments=\count331
\nofMParguments=\count332
\everyMPshowfont=\toks32
\MPscratchCnt=\count333
\MPscratchDim=\dimen276
\MPnumerator=\count334
\makeMPintoPDFobject=\count335
\everyMPtoPDFconversion=\toks33
) (/usr/local/texlive/2025basic/texmf-dist/tex/latex/epstopdf-pkg/epstopdf-base.sty
Package: epstopdf-base 2020-01-24 v2.11 Base part for package epstopdf
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/grfext/grfext.sty
Package: grfext 2019/12/03 v1.3 Manage graphics extensions (HO)
)
Package epstopdf-base Info: Redefining graphics rule for `.eps' on input line 485.
Package grfext Info: Graphics extension search list:
(grfext)             [.pdf,.png,.jpg,.mps,.jpeg,.jbig2,.jb2,.PDF,.PNG,.JPG,.JPEG,.JBIG2,.JB2,.eps]
(grfext)             \AppendGraphicsExtensions on input line 504.
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/latexconfig/epstopdf-sys.cfg
File: epstopdf-sys.cfg 2010/07/13 v1.3 Configuration of (r)epstopdf for TeX Live
))
Package caption Info: Begin \AtBeginDocument code.
Package caption Info: hyperref package is loaded.
Package caption Info: wrapfig package is loaded.
Package caption Info: End \AtBeginDocument code.
LaTeX Info: Redefining \microtypecontext on input line 164.
Package microtype Info: Applying patch `item' on input line 164.
Package microtype Info: Applying patch `toc' on input line 164.
Package microtype Info: Applying patch `eqnum' on input line 164.
Package microtype Info: Applying patch `footnote' on input line 164.
Package microtype Info: Applying patch `verbatim' on input line 164.
LaTeX Info: Redefining \microtypesetup on input line 164.
Package microtype Info: Generating PDF output.
Package microtype Info: Character protrusion enabled (level 2).
Package microtype Info: Using default protrusion set `alltext'.
Package microtype Info: No font expansion.
Package microtype Info: No adjustment of tracking.
Package microtype Info: No adjustment of interword spacing.
Package microtype Info: No adjustment of character kerning.
Package microtype Info: Loading generic protrusion settings for font family
(microtype)             `phv' (encoding: T1).
(microtype)             For optimal results, create family-specific settings.
(microtype)             See the microtype manual for details.
Package hyperref Info: Link coloring OFF on input line 164.
 (./report_fp.out) (./report_fp.out)
\@outlinefile=\write3
\openout3 = `report_fp.out'.

Class scrartcl Info: loading recommended package `bookmark'.
(scrartcl)           Using `bookmark' together with `hyperref' is recommended,
(scrartcl)           because of handling of possible bookmark level gaps.
(scrartcl)           You can avoid loading `bookmark' with KOMA-Script option
(scrartcl)           `bookmarkpackage=false' before \begin{document} and
(scrartcl)           you can avoid this message adding:
(scrartcl)             \usepackage{bookmark}
(scrartcl)           before \begin{document} on input line 164.
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/bookmark/bookmark.sty
Package: bookmark 2023-12-10 v1.31 PDF bookmarks (HO)
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/bookmark/bkm-pdftex.def
File: bkm-pdftex.def 2023-12-10 v1.31 bookmark driver for pdfTeX and luaTeX (HO)
\BKM@id=\count336
))
LaTeX Font Info:    Font shape `T1/phv/m/n' will be
(Font)              scaled to size 13.67995pt on input line 171.
LaTeX Font Info:    Font shape `T1/phv/b/n' will be
(Font)              scaled to size 13.67995pt on input line 171.
LaTeX Font Info:    Font shape `T1/phv/m/n' will be
(Font)              scaled to size 8.54997pt on input line 172.
 (/usr/local/texlive/2025basic/texmf-dist/tex/latex/microtype/mt-cmr.cfg
File: mt-cmr.cfg 2013/05/19 v2.2 microtype config. file: Computer Modern Roman (RS)
)
LaTeX Font Info:    Font shape `T1/phv/b/n' will be
(Font)              scaled to size 9.49997pt on input line 182.
LaTeX Font Info:    Font shape `T1/phv/b/n' will be
(Font)              scaled to size 11.39996pt on input line 182.
<figures/find_u2.png, id=4, 209.1012pt x 220.6644pt>
File: figures/find_u2.png Graphic file (type png)
<use figures/find_u2.png>
Package pdftex.def Info: figures/find_u2.png  used on input line 211.
(pdftex.def)             Requested size: 170.71652pt x 180.16144pt.

Underfull \hbox (badness 2799) in paragraph at lines 212--212
[][][]|\T1/phv/m/n/9 Ac-tu-a-tor mar-gins $[]\OT1/cmr/m/n/9 ([]\OML/cmm/m/it/9 ;  \OT1/cmr/m/n/9 1 \OMS/cmsy/m/n/9 ^^@ []\OT1/cmr/m/n/9 )$
 []


Underfull \hbox (badness 1028) in paragraph at lines 212--212
[]\OML/cmm/m/it/9 ; []$ \T1/phv/m/n/9 set by the steady-state bal-
 []

LaTeX Font Info:    Font shape `T1/phv/m/n' will be
(Font)              scaled to size 6.64998pt on input line 224.
LaTeX Font Info:    Font shape `T1/phv/m/n' will be
(Font)              scaled to size 4.74998pt on input line 224.


LaTeX Font Info:    Font shape `T1/phv/m/sl' will be
(Font)              scaled to size 9.49997pt on input line 262.
[1

{/usr/local/texlive/2025basic/texmf-var/fonts/map/pdftex/updmap/pdftex.map}{/usr/local/texlive/2025basic/texmf-dist/fonts/enc/dvips/base/8r.enc} <./figures/find_u2.png (PNG copy)>]
<figures/q1_cmp.png, id=43, 732.5769pt x 208.1376pt>
File: figures/q1_cmp.png Graphic file (type png)
<use figures/q1_cmp.png>
Package pdftex.def Info: figures/q1_cmp.png  used on input line 324.
(pdftex.def)             Requested size: 455.24408pt x 129.3425pt.


[2 <./figures/q1_cmp.png (PNG copy)>]
<figures/q1_decomp.png, id=58, 199.2243pt x 221.8689pt>
File: figures/q1_decomp.png Graphic file (type png)
<use figures/q1_decomp.png>
Package pdftex.def Info: figures/q1_decomp.png  used on input line 363.
(pdftex.def)             Requested size: 142.26378pt x 158.4353pt.

Underfull \hbox (badness 3333) in paragraph at lines 365--365
[][][]|\T1/phv/m/n/9 Er-ror of $[]$ split into in-flow ($\OT1/cmr/m/n/9 +\OML/cmm/m/it/9 u[]a[]'[]$\T1/phv/m/n/9 ) and out-flow
 []

<figures/q1_err.png, id=59, 225.4824pt x 221.8689pt>
File: figures/q1_err.png Graphic file (type png)
<use figures/q1_err.png>
Package pdftex.def Info: figures/q1_err.png  used on input line 371.
(pdftex.def)             Requested size: 170.71652pt x 167.98903pt.
Package microtype Info: Loading generic protrusion settings for font family
(microtype)             `cmss' (encoding: OT1).
(microtype)             For optimal results, create family-specific settings.
(microtype)             See the microtype manual for details.
<figures/q2_control.png, id=65, 224.2779pt x 212.7147pt>
File: figures/q2_control.png Graphic file (type png)
<use figures/q2_control.png>
Package pdftex.def Info: figures/q2_control.png  used on input line 428.
(pdftex.def)             Requested size: 170.71652pt x 161.92125pt.


[3 <./figures/q1_decomp.png (PNG copy)> <./figures/q1_err.png (PNG copy)>]
LaTeX Font Info:    Font shape `T1/phv/m/n' will be
(Font)              scaled to size 5.69998pt on input line 480.
<figures/q2_state_short.png, id=89, 296.5479pt x 340.3917pt>
File: figures/q2_state_short.png Graphic file (type png)
<use figures/q2_state_short.png>
Package pdftex.def Info: figures/q2_state_short.png  used on input line 497.
(pdftex.def)             Requested size: 173.51375pt x 199.16928pt.

Overfull \hbox (2.79723pt too wide) in paragraph at lines 497--499
 [][] 
 []

<figures/q2_state_long.png, id=90, 289.08pt x 344.487pt>
File: figures/q2_state_long.png Graphic file (type png)
<use figures/q2_state_long.png>
Package pdftex.def Info: figures/q2_state_long.png  used on input line 505.
(pdftex.def)             Requested size: 167.13278pt x 199.16928pt.


[4 <./figures/q2_control.png (PNG copy)>]

[5 <./figures/q2_state_short.png (PNG copy)> <./figures/q2_state_long.png (PNG copy)>]
<figures/q3_track.png, id=120, 257.5221pt x 320.397pt>
File: figures/q3_track.png Graphic file (type png)
<use figures/q3_track.png>
Package pdftex.def Info: figures/q3_track.png  used on input line 565.
(pdftex.def)             Requested size: 160.08643pt x 199.16928pt.
<figures/q3_fb.png, id=121, 265.4718pt x 318.2289pt>
File: figures/q3_fb.png Graphic file (type png)
<use figures/q3_fb.png>
Package pdftex.def Info: figures/q3_fb.png  used on input line 574.
(pdftex.def)             Requested size: 166.1544pt x 199.16928pt.
LaTeX Font Info:    Font shape `T1/phv/m/n' will be
(Font)              scaled to size 7.59998pt on input line 587.
LaTeX Font Info:    Font shape `T1/phv/m/n' will be
(Font)              scaled to size 11.39996pt on input line 610.
 (./report_fp.bbl)

[6 <./figures/q3_track.png (PNG copy)> <./figures/q3_fb.png (PNG copy)>] (./report_fp.aux)
 ***********
LaTeX2e <2024-11-01> patch level 2
L3 programming layer <2025-01-18>
 ***********
Package rerunfilecheck Info: File `report_fp.out' has not changed.
(rerunfilecheck)             Checksum: D41D8CD98F00B204E9800998ECF8427E;0.
 ) 
Here is how much of TeX's memory you used:
 21621 strings out of 474662
 405124 string characters out of 5751109
 1041088 words of memory out of 5000000
 44330 multiletter control sequences out of 15000+600000
 605184 words of font info for 105 fonts, out of 8000000 for 9000
 319 hyphenation exceptions out of 8191
 108i,20n,107p,10974b,678s stack positions out of 10000i,1000n,20000p,200000b,200000s
</usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmbx10.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmex10.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi10.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi12.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi6.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi7.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi8.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmi9.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmmib10.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cmextra/cmmib7.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cmextra/cmmib9.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr10.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr5.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr6.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr7.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr8.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmr9.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmss10.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmsy10.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmsy7.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmsy8.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/public/amsfonts/cm/cmsy9.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/urw/helvetic/uhvb8a.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/urw/helvetic/uhvr8a.pfb></usr/local/texlive/2025basic/texmf-dist/fonts/type1/urw/helvetic/uhvro8a.pfb>
Output written on report_fp.pdf (6 pages, 1100677 bytes).
PDF statistics:
 312 PDF objects out of 1000 (max. 8388607)
 239 compressed objects within 3 object streams
 52 named destinations out of 1000 (max. 500000)
 17078 words of extra memory for PDF output out of 17280 (max. 10000000)
```

# report/after-figureplacement/report_fp.out

```

```

# report/after-figureplacement/report_fp.pdf

*(binary file, content omitted)*

# report/after-figureplacement/report_fp.synctex.gz

```
      Ä½[¯^9Òvß¿Âö¨¿æù ['bØHâÀ004=jmE}¤hä¿çy¸E¾ïÖîoìÌèÝ¬kÅbÅÿõ¿þø¿½û/ÿæ÷éóß~}­¿û÷¿þþÇ×úõÿù3Ê~øûþß~üãw¿~ùüÃúôÛÿùîGüøÿÃ¿ÿ¯ÿéùÿÓÿøïþ·>½ûý·O_~xûÓw¾ÿéÃû?PðóÛßäþ×~ÿ/ï¾¼ÍëþøüéûñíÏ? øç÷QÆÿõíç?²ä¾ÿÛÏ_øóßòÿ?þöËÛï?ÿøéÃï_~À?o?}ùñçùñçÏ'KûçY~ê»ùüå'K÷çYNýçøóO~Y%ß?þöÛÇoüWàûÛûßéßzûûã?ÿðñÝ?þþöç[þò}ùíÇ9°ÔêÏóçÿÏ;­~ÿÍü­²ÕVÔüþîí§woGQ3¨?üôÛ¯_ÞýzùìçÂëé'fß6hÞ¿ûõÝ'ÿõí_ßý|üÿÈìÛÆËÙ¯_~{÷Så÷måÂïûw¿¾ÿùÃçÇÎÿåç¿5¶ß6V®b¢wàz­f0úYçÝo¿¼ûòéíÇØ¬ß6vª¼~â_~ú;þ¸o53·ÙsÆÊß~ûùó¿üñó?þ6êyÎ(i&çüqï9Ceføydø9¦1üòéÃûÙs¦Êìûzß%Äés&ÆcïßÿV:¹Có¹åÇ·¿{òÃç?þzþ¿ü9C¤r\°{Îä2±³?ûeâ÷=gð¬³}ÎHù¿Ñ/ðÏê¿£Ï*¿¼ýòx¾úkdùÁr°Á8©ØçA¸ÇÏï¾üñûÈð9#åí/É³þ;²{ÎØáÏ/#»çLìÈí9ccàö×ÏãÜa34n¿ý>Êæ30Nýë/#ç¿Âáýòö¯ÛásÃ!8ýñóÛOã´á?Þ~úôvl÷¡ðéíû÷ïþfÞµ#Ãçß?ÿúÓçÏ?<¾ûùïïFÍuÏ¿|øñÓoôqû¯åsÃ;6Þ_ë?FÏ¿ó{ÎØX}}}ýsÈ»_ÿøåÃw¿´cúònÕoÎFfÆ¿*ÇÁqøÞ?å`ÑÈþR½gÿ§<¬¿ýýtýí/C/ûçGtë§Opê±ÂÿsâÏý×ÈòyÉÇ¿ã3?üúîäÛÿYÛð©¬1XÞ}+ó®ÿ>/"A/ÿøËß>×ÇUç%Ç.<ûÏí·ëüúÓÿ þÓ§þkdúQôøcLÿì3~B×ß]Õ3<kÆy÷Ó¿ýQVÿ#ÃçEõïß}ùòáËÏï>#êùõýôçÈþyÃé·âeî¿FÏNhï~~L3µ~úýoìºY
ÏêýBg»þ;. >gêéòýq/>gÞùãÓÏüß(Ô·¶ÔòáäùÏÈîÛÆJg÷¾¶[ÿ9²}Î)+jo¿üõÝûÏ¾Iþ¹9Gx*ñ9cäðïï>ýãÝ¯%|CøÇ¯\|þññÝ§?GæÏsþøõÃÿõÇ»b~Þ}ºþ52×ö+:êÃ_ÏÆ5î%äßþãÝ'þýûÛ÷ã$þì:rå=ò|ÎÜszÿ_ôïÿþVY=gÞùÙþõí¡Nýô«ûWÚúa3êíÕCLß6¸~äú×/Çpøåã?|þã÷ß)õ¿ð¯Êô9#ìÝï¿üFßÿþñ}ÿcÚCJÏiï?ýDFfÏÊÿ£¸Ô¤üüqÕ0åBýöÇºù£þ¹ì»ðËÛOÛq³êYÛ5áÇ_ä ÈÏnC°ùåû9"¯ÿøÇ2D
ßý·ïýðÓß~);ÁJ©ïþó¯¾¼Ößýóú	çkõÝÿÑþ»²ôåõwÿ6F¿ûõwoô+ÖEL
¯¼õA%^[²	V¿r)ç­yUN}÷x%8ÿ}­^á¿ß½¸<«Rkü<Þ·ï^~÷üïÍ7Tù½
6æÊêã¨U4Q¥Âþà$°Ñdëó@ CÎÞ'$¸!É+«;úÕ¿ä_
ïõªáK¼òJ½2(¶¦ö|ÇÝõiìÊC«Ä¢ßÌ
ãTu	èÁûÕ§oNýés.«ÇfÊM³ÈÇ'|ÚÀ®m÷æiÕÃSe;oÉ¶f¤&+¬?)±ñÎ>Q	OäMk°kk>
ë´:ûµMMûär¹Â­w&e¦åAàÆù¸lïëô¿u£ÞûØÇ5üÒê¸lö{hïh¦¹_,"ßëÿN}2qçq
¿\;hïh»â]ñrAS b²2.mä·3:ü4öã0ü8M0ß¸ý;*Hü=xîû`ð£ôqYºÆæµ6ëÕâÍÖtnaè¿79)h+öÒwàÍ	MÔø©¦»-Kíú%Y?N$1º®¬wT»òïÚæ9¼?.;y©Ê
&UùÔ¦ùïÛª¬-|c=«ò¥t;©ò{ýBñf	»¥4³wS/¾å§
ì¦ÏÕ=A¢?I!¸i1ñ:|!ÕéjÇæâoÚäÛ\öëR{Ljï©@óß¬àPÎ­¦ÚO!c?Ð¢>JÖabÎÚ}®£n¸Æ¥ì¼S:¼Ò	îA¤@ÄÖÆiø
¾¡§À.Ü·~\Aï«dRêå/Öåg-/·ü^nåÿÞýw_O°Uß®a-Ñ:Àû" ËmEm><)MåFio²N3Íhç¾qä$øëb4F>èÍà\ø ©å!|æò³ùÒÍÝÿòèÛÔÉFJ2DéNmÌ«d²9µÉ¾ÚÃØudã­ÑÑðÑË_`41è`*°§`Ñ¶¿?c¹³ªØÙÒM6Ú<v4:8k-´ÎÍ\´V	Í¨¾.&¥%~RÖ/ø£UV; ÂGÁÇ9nøÁØøÁfï´øè¸èÀb ÀÇ7¹èD³¡<cßÚÁïÅàÎ>Ñ#WÔ;ñaS^ñÑn0÷¢^c¡¢J[Ág.oø.t¥Ó¥LðíkùÛap\ñQLÎaQ5l´RO´P~Â÷Â°üVyñ­Rf,×\3¨LE2?Åå9_øW>Z»l´ß¢âä¼¬bP÷v>­D3:³8NW9\)¥Ý«l­ÛáFÌN¿ò¨ÂåsªhEv;.ªx->±×	Ò·âßÃd$¨þÄ&¥ ³Ù?ÅÐX<§­I\§­­)ÂÙ¡;aÖb6v(Ð+&6yñ¹ÛY` ®GWCÃÃz6BÔb¯DAa`ÃBLD#ûJä*QÌø¸ÔÎâã=ñÒ+ð)j#lBL0×&ß´#qFi¸ÞÄµ3e×ÖÛÀdfg£ò)«MKhjR`RvIm¥ðÉÁÔ!¤½lB¸·kÅvnc¤ðÓ9ØAáFYA`ì,|ÿfLÏ_e8#øêMm0yðhh*£8¸]S÷Ö;Qµð-wdT{X_AæmøÍ·è²;ÆeP$Ì[öa;º0í£lì0D,F ÷ºfÃåÐ9¢	ÑèwMËYÎAN3ð­-`¸à%È¸·µßÕÖB²l*­¤)®Âö¾Ðð__ÁñIç´lb7SÝ3 ©ÅïÏbÌlÖvìuq5mè ti+~³,À°ê¡½éÑfØÝM×ªcTþPf¯ÛÖ&ØÖ -Çàé*B­;ý_oxþþQhÃù®Å¨®ÙâR+zë  ÜäÕ­k¿FôêêC]ëuqÅÇ¬aÍ¦ïOI'ã1¹`¨ÍØp»1ÀZhø;PZ7£kå0[Ç$ð1B×ZàÛd¡×^àcLÁbEï2è¬¾7Þ1|Ë#c-eÇ/Ä ð9w^à'Ød­Ekg<ðìµì}-VÙI|x9	þ!E¤&ÊÒÂÍåé4²<
J56Uùè.$½C0PÕ6v[váÅ¼¶¶Â±'FãßØ(øäálÐÐVgÝ$øðèKF8Ñt@§ð7ºàh¢»äRr½£Ð¸Æ¦V¬n&+ðÄ>q¸Íø.-§Î¦æ(ñÙ)Óïsçóf]¾¶¹îJjx³Qÿñ²«:8$|ñÀá/Aªaxç^yåBhiÑ@pL<CQîQîàè\$ÉÞ~¤+;ãgéÞ6A6ÌÛë4vÆ»R¹#çòõTC¼ E÷êFGÒWô4?¢#wÌ	AqB[t$+Ë
óA»0,Ê]»¢0 RHI[F4´±ðÝ¡Äv *-ðT,ð7!DÉßìÓBN¸VlN}U#r. ~ä¢
¢º[	NE¹©Dt¤Gª-;ÒSÏ¸ç¾STÆå|À%æÎ{¹ókø¡ËÛ'F-ÉÀ([QK]?Å#T­MÙ´ÁcÄ*ª`XÁo*N0`Sæâ¨ Ú~®áaÐÑ5cnHs=8SÁé¸pJ.»FfRYd[õâ®6!ïdíVH´µb3á/tBºFªà$A&Gíwµ¡W}Jyn_ÍìEbì\A<nªa¬ÉRHFKteËs slÚlÚÃ&>Ë ÉÊDD®Þv¤ÉBx&ÔÆpÄd·#UÊ«uÊ»4p)bHb XØÍ5ÍÕ\ Lv#B(»"¤`¬Ü0Hñ7&÷FV½¹,+QÞÞÁÑukÝ~|Æó,æoz±Ö¨EÄ½üëQa£v?´KàÃ4)¸Ø¢8ÔCÃïN\Ç@hè$úS0»´(ð°4VsPÏ,2¼IÁÙs¤0t o6¥z ÆÀÊ8MxÄ9Ön¹^äëz.SÞÉ#Óßxü¥}Õ\g^4BQ¿J~¸M°äjÑ³0#:D©"\óf¸%ûI j:]ÑåSx0¯>À: ²
»#V°&ùEÀ2>-A¾£¾Ù Nì]ö¦.pYNÒJvõ] ³^ëDðuÄ'¼v¦JOÓ¢ü'å¥y¡±6jO;@£(Bûñ4ÕH¶F¸ÎÄC¼Fm±á¦T1#¶ÐzpÛ¶è«]¥Þñ}[©øu_©¿¬	¢è:s¨;¼<ä9ãÃ¢Âéw1KF´10úx `#|hßÊEÌ,¡bôfgûúá
>Ì8KÁ¥ÃeV[tZê ÅBÌs=a²@<>Ç£!3Ø÷òßÑwª¿¬^9b¢Ht79§Êz8}<Êé.b`ºAvC÷>ºÑÛkÂ°µ½Rÿ.ÛøÏp6ï!F}1PJ1rHrÌ94
uçA]ü7A#µâBf+YûÃÉg&ã?z`WúÄl ÓÆÂ¤Y'6ðÃä®áÎ@9j;òö=Z½ gnJ£u3Gãµ ú©·#çTNí´uØQWã2C£.ñÐ¤oIÍ¾2Ù§FI0ÿÒây(¬65Üß°C!¨Ë¾ãÃr±âs[;DyÞ9Aøf¿Õ¤­¢o7)J.Ë.tW×²I1QÃ×[4iÜiR¸è¸mRSæ¥«kËtõõ° Ü<Mê¯cÓh:mGõ}Ø*Ç`Mà&Ñn$!L
.#©BO"¹¨àÂ<¤;GêAoH0¡xê3Ç£q1îÑû¡Np4\¨WÑ
 Çì/86¹_s0zW öÙp]]..-aàbãhW+ oÂÏ··JVÑV!úÍí ÐuÒ¶cXñ4Tey¥jÌQ²¹,÷&®e,(¡øÑ©Jy«ÖaÄ¡X~Cg3 ¿ã4á¼[ûz ×-ß^+ÇlCÑì³ü¬;`RN&tö0~ª®;j`T²ÈzÇS5ëïÁOÕ´6Kw$ãÑJÇxq-]ÎÑAH¡y/Ôè$<³pkEhb¹¦n¬J*Gø³ðï20:±aQK<é4´=Jx5Ñ$	â	³1rr²³xØ¨(è*½Ýµ¯ÛÄûº^04nE1ELZ4U¶qì²fö²×9ä­(Öñ(fÖûÞbá¨%½ t
ã ÛsC FX=«ç1Ók	Ý5&.rÆ}	c ÃM¶Zæ³àÈÈ6ëíÇ#VU</¡)óëvÜ[Ì[1fùñøîàsÜ~<&òàWJh-æOL#£ãöhÚö¼uÇÔ²äØÍæÄq,í©ÏA/NöÅÄôÁÜ0ßÐkåÛÔºZþþ,ÇHsVâ#:Õ&|qHYâG´6®<ÅoßVåç9Øc¿|q|;Ærÿöà?¯K N_B®,Ñç+#Ø´9P£OZÎPMSUË?åû¾³4¦+YäÚË]tWå1sñZTÃcÛ°
S5®aî2\_ÉrLpK¢îa
>¢ÙéÀ:ÄÍ)
4X®sO«^Ü5Ò¦F¯8Uy÷5ÉP¡°áÞ^êâ¾ Ï½Á°ÊöûÛ6Ä1T 8ÂJfh.O9¹64¾I4³F«f§³ |³ßhP CÇ9qÜhQàÁZÙ¢ÛºM
wÞ©MÚwA³Zì7nÛ¤&ã¼MjËzêÜ¦B|´ÒIø°T½çÍn¨j8´3!".Yú7`Â.Ô"pú3z[cRpÈu7 o£u3GÇ£µv Ã,aÂr£á~?KÁam40"ÓAr´ *ÿàyÈY
 Ýª²Â½ÞñêB¶!<¼úÙúõ4ïb5+ ò kÚÅ+Ú{·%ÁÓÍnA¢ç(÷Kù
¼°íÛ¹ÕháNx^,FÈ\#Oº$§>mVBµ8©MÑøÝ¯:|'8e®¡·Ù}¦m/ç¼2Í8XÌrw"¼ãÔ	&A-È»HzG(Ç@{k'²>|¹ü]|C#{[æÅ¿e«h÷=êm³p
Ye· ÄØM9ì	]È)ê,	a»4Ï)5BÃu×Mçx!4JK<øºj
u0ò _Zà	ú\üñJØ è¢Ïfæel8ÂÏÙ/8fà¦	^E¨ÇFRÈaW¯øä¾Ï5ODEUZr	<çª¾5ÇYûzUÎÈ60Ö¹m£·;½$¤+ä´ßbb)-DEã»¬æ.;~Ñ¥0¼É£·³æ$G}åêhI÷±7§«Á|´ ­¦LZ9,bÙbO
ñÌs`_V®åêø.ï)ç.&"+ÛQ§/ËD3!JÛN²wµì$Xk°fÐgy¢P"P)5>là£çqÓÁÖ& ·Sï,¼ORqµ,g¯¤,\x#ä%:üXB½\ÏGåÇb¹6[Ë6åçwxË«Ï¹Âr§¢àqGk'^o9^ÑO$ôïh&Ið'Tk(­V[Æðj=þ÷ÚªÃÆ7ïMø´È²Bº/(½Æ¨_bT!Ê~<mÅåø¦ÀñÜ.ü*k^y°:S -ñ53®ø;üZ¤¸_6`Irûî`KùòHovàJs-x<s­¾óÇ!·âþ}eÛg]ç!<|4·Z¦Æ¸ Rß¾e:aìxïO0ü¤ü h¥³Ìé\4TMg
Çð+³¬pGVÇìÆ²Z-Ý2Ê¬`>xEBqQñ¤ËODÜ¡ü)ÝqýÔÀÑ«µmhè{aº¡ýßït¯Î¸ìDo|:§Û6*ïHCb'âq¶ã qÞµ`Þ wÓ('Gæ/Å§bpÂ¢n?õc­g´£¢CÓà÷ti]9\Úà]ÚàÉ[M?Gò+âaFÊÇ§¯2V'2¼Z.¹6ÉO#ÇG?¸0Ê]Ãõf¸ibôÌssÓ.%8DV/íÒ®®¥]BpN»ÔwÆz²Þ¨¼µK'ÑòÂ.qÏ-]â	^K£]Ø¥¸#+7s²KÐ4àÓ3Ìà¶]z*úã]Ø%x"p=ÔSíÒþ»Ä¸ñ0à¶éÀÿ0ú]âz'7<%ó¸Ã\óXXÛ¥ö©wìÒü©³Ý ·íNûÔÁîÌ*wÀæõSæÑ*+YÉÞ3Zº¹p;Û&Ç=ÖäÊIÊi¼Ñgä-1Þåñ+çäè:Þà¦#¦a).pç& 'Æ[Ù5ôE­Kà%')§ÅmâÌ ôá	ÔêÕZàQ=e¼×¯ñÜÍ3PÄg]hâxx,£a±/ a`¦È;£G½§ú½¸^÷ß®ÚÊ{[/íP,8K'xí3¢¨,vN*2×^ËFØÁ$Ëy¿gz»ºæÙAð¼P¤¸£HL:R)_'ã
)1Ú`Â¸LÆ ll(bÁÔÑøÛ°Ïj1«AMFä³¬]j½t?Ój4 åÄKU±x×³==Hôj,ñ8¥sÍµØm*Aw,4»òVX+´Ê¢MG¶<}­ëÈÌìyÌõcr
GB²\äÕÅzi½â¯ÀjJàÕëºl¹ªYSþWCVOÅrÅªÃÂùð[èßWü¢4ï4n+ó+Ú©Sïn ¶3vðÇ
¨9§.ýS×.ÝÍ] Ãëµ»Ô£ÒpYüìÒAÏq¹Jø-	É0?FNp¦e6¤`Ðµ7ø"»ÜÎ!:&ÓLcPÑÛz®.z×"bjÖjÈd¸á»Hè3ñ§ 1Ûq>Öö¯þ¸áñZµzbÃY ®føÃÀ
jµ 0Í«Z=o äzµzÜÙFÊ[ÁÑÂ?à²à#ËÛ§vâXÄ/áµ%¡æµ+¦EÌÛ­a,ììØî
dM(Iø*èÍð+ä#+Ðm¡3á(Æ`ó¥÷N¼AÄgÔâAwê2N>ëê|ÐkVÉö ¼ßZ.Îf&ÿÊë$¼ÜèRC<gÙ
æ)öÌ4sOYfujnð§¼Ñ	xcÌíqÅ;|lÛJÈaïu#t°!ÌX!1óò±Aøf¿Õ¼Ê!w³!JzÕ»ºHûNJtâ ¸Ý·4LÞu"ð§Qåà!æ¶¯ÇY£YÏüY¼ÈE3Ïþg3¢év>Vôn>ÊÜùð¼¤ÎÝîäÌv7*,íÊ0uóh%2©èDX®dÐÛÚÔhy«f¶6 ÔW6Þ3f Y åa,!·®bB¸
¼säû9fàøBrYt7<ab©í:p@¬-çºü;RÔÊÎv«Ff0ÂÐ;<¥|dx»êÖL8ÖJÇ SS:ÌLÙùæptzÜ»]
v #²eÒÎ,;<lXT®Þ8·D;ëZ
v^`±¼´Ñø¼#çè®³xÐZêâu°ó¦çµÔVA5±äJ
]ÌÓ"¿ÙÁo5YÓ	äfÅ"®lW×²É¢a¢l²¸Ód<¥aì¡a\-,ÃIbàôqÁã°P<<ô|·Ph[®¤î³Ó<ÐÔ	»ËsâàÈK¨<¢6qì£÷¼×12Ï[QÉÖ¥Q!Û ç
¸Åy¦¢¶[;Ù	«cùMÃë¦éæº¶Ø²OÐbêpGáåx±¼¿uhÓÙaYº.%}hõÉ¦ÏsZ¯xù±^vwÕ¨ 3 ãëøÎú¹²éNÄÈTÞ&ôðpBh+õKàMª3r],Lo<np®ÉRWZ^·ò³Zpu¶ã¼N±ñº}Õ]Ô<ïe9Onù ªaªVhÑTM¬dL.æJG\Éà·RvÒA;0Ð¢ ¥Ì±£ã|©¯àf)yÁÖ#wÇDñª t°(Z§ÝcbaÒj'	1~¸¸­6Ñdoàòãò"/YÝøy-eÑ*é8ï&öL¹EU¼VÿnûÕ]Çû[3¡QL&¥Õ÷VY¥³(Ü-ÿºà8º!ÛBÞzñðë^F±KNôVn«­ü¬Ûyº¦º³gR*åÛdÙ¤'v,=Ø;&Û
C5¥µme¿xbØ/òôãÌñKf²é
»r5§­gºrN·òvÅs3qð
 Ü,Á_ktÇÑLTIräEìÀpAF/«
&ty»ªàßÃ&É1yOñÚ- )1ÍÑBÌMa3®é£¡s<JÂ	Ù	o¸5Ø&xkÛ¯6°6Ç£çoè£gæ8ÖÃ*!^
fVrmÃ'ä¤Ôz-
qwÎ²Ö½{2ñVJ.
[«®9QÜ,ÃêJ)1ñ@fsl½ÁÄ¦g²Ì´E¼a>À4&á9z¸¨²à§Uò:s°dÈ¤²>+¡æë=Jm	ç¥H#	
ìNnKyIèOºùÛAäx@`ñ!f ÝÖÖ[Ùj6eH¼(c¥2ðñÊl«Ê;ÁæÄo&"Úr[¨(¦q>n1®U	ÀrYj-WûÄ jûgDî¨^Ë-rPLå,©ø¹Qµ°×Á]à
YÇL{mpAq3¦-AÞóJÿ0µQÂ÷dõ©'2^O)Wµ®ÒS¥c4!"ê`&S?qìÑr!p&	ÉLøµù¬Bè_jKýdÊUù
Ìæ­Æ¶ì3MT£(]×ñÕÆmU<=Ä;0:^Û¶Jvö7fB:µ¼1¾©ÝyÇk&Ä A ¶­ÒÄL8cÅÑÃ3áÎ½§ö6ôÁÇ°ÔÊæÝ¬E^÷,ò´¸qY@àÚð{#«å¼å¼8\Õ Ì|[e#<c(Aü	Q«A¾>qð>ÇG?ÝÏÇ7$Gcøz1)Àùbpä`ß;B¢Sc9¸!¤0«oð
Á4*nB:»heÞò«nÁ¢NYULCvýÌ+>þfÂòr×´62bºAÇ))îÇ8Êq5¨^ûï
²¼év ëÎ§ã1½³¼.	paNÁT¬¼!ÊSEÁ?ñÊ* UÊ6'+åá!T.FH [×fVz¦3pþK5SpOYàaè+	ÛnðB LvÉÙï§±O%\.;ýs­ÑÁË°yÿ­çx¹,[ùeÐ(2wÌdÂ]¤¤ZRÀw`ÒqI2Q,êà£/É/XaV9VEfVÌ÷ü¢O3óq
« 1ê8aã`¹@¨C,C/«Ý+ëvI¢Ïëvé5ÂðÊmE©åmÝ®<Ë=S,ËÊ}¬dÐ¤ Ëà½ti*:5JG±<8K#>p¨­s¼hYèÒZOåÌÅÜ³XÙ»Óv,¿FLF~¦^óâa'Ì|Ñ±2Eó
BzXAû¦Füá@hþHçlæHt¦.Þ4-sÅ$.B~Aöû<_5bÂ%ÄÜòÒõI&øOJ¹­¨¹$¨³³>ª]_.GÑjp@"O.K@c×ª ÏpÃðÌMhÌbÙ®FLbà¬¡·±KRQ^^åbÜZ~ÊËÛ`®W¦.¤F,&vYñ/Ø!çU»Ý¹9-ÁÇ¦â¬dã²¯,ÈB= Z!&è_ÇíÄR%bñ´E«âe#ÞTb
B|É;pøKÓÅlÐ\Õ«LøìêC·?õFnôÏBøjÝ Þâ!¨Fø°4;-RÙ¶nù¬íØåÍXåðOpdlçlØ)ßÕAÇ'ÉßR²ÿ
åí¢±¹&è¬Ô/>ÁT¶x%@óiN#«
Ì4a»J¾+ÖNË£ÓÂäûeyºq°ÙZ­$!<ËÍà\^"ùï¥ìZ¥ûÎÓÇwßùÂ1¥
eÿ]®ÙNËâ|IÙÐÅ%ÿ{/¯³væãyq`ï1)ø6k3e½6;ÎÚ\À¢<.rYQ
C¯`5Í¾cBôe3nbx
3±c¾$f;ÝHÍ}`ÇôÌº  Úc¹"F×&Æ±}`/FcJ:ÈÖLÜ¹çµQÓ¯	qî(*K×T/ LY ïåfµQ*rO2<làÕ KÂMª¦ËRÔëùÂÁa¸aN·ëÛºÄþõë ÀÃH¼*²aÈ÷®è=UÓö<í.ÃM¹NUrÅlÏmïäÒÖWÂÄB;ô÷ÃÐ,	ÊF·Õ ÄLÖ!ÍûÂhÃYõg
Í­ÕBÙÁ.e;RßM"Â>5k ;QªTMÞòÚK-çö.[3 ÙÂI8XI®vçnîÒHrÞ6ôÚ(s½¼IÄBY÷r¯Ö½oËû´¡;¦Õ4;Ä(¼×¸cÜc&ãØ½+m}/Ò$LXyçh¸1UÄ%ëÕÐ&Ãßís*ïÎë ÛWÙ|÷eèçøX tú¸W¾ù6Í*<¿$¾i$i-À3M}L'	sV
Çð>)»ýªè_ö¿flõ¢ªgÒ	eÑ@ô9	;ºÓ1qìNÇãètT«ñÄ¯Nv\¦æûÞ·ÑähMRCî."Ö8÷`¹¢qÜ[a5H°=ùa£æçX!Mù9«¨0\5/¿ß`t!Îò®½*¯1±³ñ¶uÅ5;nî¾Ù ïîÜãÉê;M¦@áKíSPE¹!×îxábÊåÐÆåÅTÎ`	F}9M¢kxÃípBy#22Tò c-?ñùÈßÃ¥Ë_G®Ðb¬]úÛ°wØÚÄ`Î wÎuí
ÛE^§hØâgb3Yº/½x`W.½©¸UÞµ?ãár±áÄ;Ïhá´äÜÊïpÎÌÏFÎg/ñ,qt²þZªm¾¦Q@¹°#(1&¤ ·ÅÔÚ9ÞàX´A*©Íó¸ÃÜØÀ³|«F0>ÑdÏ®æFà@¹t5S`¢WÂnÑlÀO<ã!¾¬î|·kã½*ÌÕÆ4¨°ÌRz·¡y¤º¢[ðwî]Ä¾óëÜ .ÀÜèº¹<C~"ÇSI\[kOírýIòCzî1;/ìtT£^êEôJqeZPDú¶YRð 
^Rð5K­RaèÛdÊÏ´mêãÕæÜ-Dö1}ÚàBwãàyTZxiÏdQ³)YÎõ.êV-ÿxs´ßù¿Xâ£¸²K@Ä^^ñ_)Y
1÷j³¬ìÙw¦ãcïKíyÜÂê(	yp/ò0ÿQvòË¸áÍ9 fhÚchÎá%Ë&fDÇa+÷Ö»³(|[Î/:Á)· ðÁÚì¥<,ÅàjîçJKà 	1Åò­ÛðLÂÂåáq°#Æ-!÷&cl¸êL_«²ÝWlU©Ådº
X[Á)ÐßáJ½×ô³­üc\:Ù#\1í+ÙZvagxzÊ3 ËyÃÖ¨¹úcf©z¸4ÖÑçìÕg¯[åcü
v%QhÕ0A~°¼q6Õ£ù¶+Ñî;ytÄ©Þ\|t+9rv1éÅÄ;)3!Cåª;B6ÌßIÂË²|/@ ÊÙ¤äeU.ÃÓÔj[ÏÏº3!­¯w×ºo¼ÂÊnv0iÊ
R	cI¦á¶ôOÌBãx¦'Æ¡Æa ,§æÑs!C;×²¡ÉcÙãÈ®#¬Ô½}.;UôVÖ1røäò¢<37Wèåç×Í³ý¦ÃM¾M\²c?±ÃÄ¬ùP´d«¹õ~e×É È~¨­
b¹ú:ÃÛäA«¹Æ}JQ=37èý×hf¸çó­3¿nÛÍÜÆ;$µû ].°xÙ@pB£	é»¼ña?ËájË´ "ÈáÑ»8ívcz-/;ù¦Ë®éÁ>3#¿ôÏH12r§á;]%ÁÜ0xé¡ax*ß.*÷ëGê ÇªàeUpý¹i¾ÕÏûÞ*IBH
guÏ¯ÊEÓÀíRæÀ ^¦HR÷y+sR¸fÔÓW¡Ã¹­)»!3IïéÎÖ f¦ØDãSB¸Á&V¥\ØÄrê¹¼ãpÓ&2)1§¡ý}¸»ï=z¶òðB1¢¼ì DyI0M»>ñÇ<¬Ã f+|*E¦òl¸[/ëÅôå²/Y[®ó¶
vAÁÓz¬ºRhÃ NRáQ
Ãâ@¹´ d¦¥¶P|d½[%#ÑpYi¾ßÖ b*­ªKêãSìBÕ*>%§\(W­sä
Q©Ëgl=HçæhÞdL	Ùr'üÇ
~O-É6Ð5TRño¹>ßÓ!¯è^ÞæÍ·m,úIé*5ñb,¡0/ràsÅZr,7ôÊÏr£ÅÏJ®390{Þå#ÎUàÈ«tÜJ?a'
3Z£{ùïRê©!e9ý8ç1êð
±À×dÄæã,þ¸Ë^ª4ðñ&Jô¡'r]r|IÇ§º0oÖåÅ²b&<1¯Ê
Ê¹=bûò¶yÃ>3\qÈYçý\Î{å[½*«^f®­}ûYÛÊS/Ø=ÊË¥f¥£ ÷wNe /oæf(cLOlüû¹Ì_±ï5~­íIå©Ï_ËEãO _ jWÀ¢ñ;è/Þ¼çð!Ê£r
ctÀ±¹³°%Õ\ùx
óðny×¾O÷ÄÚÔ%ÅÖn_6éù¦vYy`ÚÍ=°<ÀkÇEq-¿æ#Sæhøå§ÁÅ¯´0ýÛzfÛÏmoÆ_óª[¥QñÙØ:dúp:à]Ñe3à5Ó¯7b§ÙüÍLÀ1×`ö
|´îP0ÞÂ&Ù/±ÇÄ±¦$3óÍº{h¸.'_ ËúNìÞz¾Ü"bI	fJ[ð¶:ïBÀË¼VÝ²ôoåsÚA¾4~l³Â'l<ì _O ¢©\¶f
xÑw|³º^ÞÉÆËÁù`Ü	ú·òp¢ßï³Â'lÈN8²&
0zKU'4h<>}ycøip¢+Ñ	ýn'Ô
° :¡D'Ì :a´®E'tÐ_V½zßàIVû6öll£ÅÁåGiµ»<O±ÚÑhµ_¶9íÄmÓ¬ø*Í8WòRî3å>õÓ±¿Vìdù|»ÇðÍ.ú¿8-§LýåÓ¤J*{ëöª{Èd­¹¬_9^	«¹0<Éô<¾øp¿àó©^þ²RD>ÔFÏcâsÖ~¬áäø:.w&üÌ«ô\¨¬årÃÅÞßhNª={À8ZÑöíu§F|
c+aû­èÓn<P¢BÈ:Jù>£=xÀºøðTà'ï¥üÌÃ^g¦­·AðÉ0cF+Ï|IJÛD3'Í »f¶ÓóÄ °Yx.pYÎJ
î¾qÝYRdîBi	fz^L =åP¨¿³YRÄ{#IÙÜ8ÍE#ò}=ð nÍé4ÞWc<ù¼8Ü tE¿í®i#ÚÊá¯¯0ªâÛAS 6J<£¶Bi¥WÑ
Ï¨­0N!zÍmtéö*ë-°¸éNù~BL.íÑuNÊÄýqÞW!Öú°J*YÓèËáîåiÞl ç^D¶ã}âò3¬ àä0¤;%1Gm(7|w¸¼õµ.åGÑòË³a' `Ðð¨ò
ÍLå^©¸&<ñyÓâ@YâÉ[*Ñ±wØ×3¨¥LøvD)¯mªoØ½r¬Ø|÷ËÀ*¼ãKcT
ÑyYÍ73zVx'Ì¹,p]	u×.1ôº®¤+»(êQðÒ]Cî\¢Ôúí§C àyñ¼"Ïë\åM\ýòx±Ldmô3ÍlzÅæ%§ß
8ÆáûÙHl;Ä}x1(s¦Þ(ÈÌ©!Cø¿Pý¥B_*§±FÌáôl³ÓW]FT$ËçvåÈéÛ#9+<WÚ¡Ü@ydÞ÷Ä\å#ó:æM*Ðwe+­õ|ÐÞ/(0	Ù5 wätå(gX5j¥YÝ®{Ìy1Çù¹Ø}0³*~õ6X`µvÀZ¹4Í=ÂEÓ{>ßMÏ¾Ï37}`Jb·hú )pGNøúGMA¨m	æpy
|ÖD4=òõêQ^l 7¬Ñg#±ïY£ãcY`­ÑtË­k´@8û¥x Ä¦«ØRÁ2$©?ô8AtÀ
àwÄÄp±Õ-÷sßXñÿh[×ÊABPy.YèÞ ¸­{ðúù¢[zµöàÎò>Ç<¤x%cF"x¦"Í^Ò±wØÂKbÞHz4`ßöÓCòÙÆKb5ÄIxIÆdyï%uà%¡¿:/i ÜöSï§7^Ò ïÞK%XÔanØ¥ö7²Ø÷ìqØ(a. [vi8Ø¥B5é1ÆtSÜaËkjsø ¸3-âoëV³T¹­ÄÕ*É¼î1Ïw»Ä,Åw5ï¤ùIô>´\ÔÚwjuLçªÂ\kqYR-lwI%»°Mð6±ÈCàV,¾îÈøÄJX6}à[nÑôp9&2MZ4=c 3èà´»9 îÔcyísÕô\<<Ñ©!aKþ:Ñô\1ÈiÑô)
ô¢u:àñ½[5=4¤¤WÌ;àsÞF4}bbtxä£6½Ø nØÀ?ÏFbß³Dc}FWxÝ²kÄÁ.ªoFçÕ§8Îi/pvÿg·Õñtì³3Û{³¶÷:yê¤æeHËÉ:¹¾Vy:çÅÝ2U>-ÒUxV®¢^Ê_~ìå¯]jç:ÙYN±D©s9Ì/>ÍÕpk¦(õºÄIvw
ì4$³Íü6ðÊkLÉõà pÅ6Ååó°1ejçõÏJ\ùµT|ò/4Rb)Íü9Ú{.Þò¥(mºîHÍ7½øØ uñ¹½ÊÈè­Tu\\"Â--á×6èMcTð;Ô!7¢¬SLL~Âè§¡	(bFB&¯Èßo(_@M¬ÀÇ±UÛßP¾~FÃÊÆÎ}¸z­h|-Ê·\ð^PàSÝØâ/à®=8êÁ±Y\G3ëeÝ¡ÜÎ^(/ó'*Þ4ÌõÝ£åàf4mýjù(>ßîÞ8rÙÇÌd°n¯Êï*óR`ÙWXì<¼ÝÕßl ûÈÉ±ÙYwµàËRyÔ­âÔÙâF]2'²gL'-à ¨:´òPÊ6×äfª`­ø¤lKÞQþbOx¡7ßj2Í£ùÜu¾ÑdÌW¯Ë, lW×²ÉhxïK4YÜi2tSMn1;³É
)}|Élc§E'ãPá+4A7Ä>¢²¡Ä[s>;º
×IT4%
	ñÁ-Ø:
ÜeÞWÂ0hFRÂÕÁÞ L#ÎÉðao1·e½p=0'2ÔKB+a*Ïl®#×÷ÙèKÓåU¤vFÁ
S¤¯³»øa¹>@GREmÇ«ÅGp£¹v ÁÁëÄ« ¶rYZkyÉHë_	ÛîúL8Ú°Ñ£æBÃÝ±¯:ÏiòD¯S3ÇoGjyõ£S%ðXöé´ç-u,W)¢ÄGG«Aü_N?[+ð3¯hY¯Ötõ¸ßÅb%/ë§1s0µ °&yAá³.ÄÑVÌ(XñN\ÌFR$fL3ZRðyÈ~AÄ/ ñgaR¾1=×axhyÇïñ=LAQr-.êeÊeùb¦`Òaåu¸ËmAÁ¼µ|APÄr3K@2pAõ øx¢6Äê0PôëQ×¡@7@¨só7+¤éü
QzÑpþfQÈ%øH+´ ð¢PBt@<.n!®ç°!îñÄ©ÈsþWüÇ~?ãÜü°K:¦+àÁÒqécIùòH_O°s<:åû
ð|ùn T
ÏUô²g
xªLñ Æ£ÏúÔF<nË=îÐûÊã/
zÞ1#Òb·a0ý7b?î°§½¢ ØÊ%qõÐïª¯ØòÐ¤©b/IbãÝ {g\zª²à<èM5§ zg¹ÞSizÕ_ct¨:"úÓð±Bvvâ×¶q.¢5_Öà<$¸Ã
®%¦Ø²&'Óñ¦Ô»~ÉçÄM×l&Î-gOÄ¡½:!¥Âã@\òHlÂã@\ñÆòû<NÄ3j!ÀøxnW1Må+?.£©áÝØXØ	Ê¼Áåx§ßæñÎ#õSäÂ¡éH_-(<g¬¤È¼GÆ¨t¦È¼O²463 dP×ÆF KcÃeï4¼×|?^Ýz´?
ûq-M.gÑ{ê+v36öÕØ|lØ¡¤Î°¶çC®ÎKdóQÙÁD}lÅ#;;¨iÿj¦.²©±f(t¨WP©·LQ§pÔ}Ssó:<
2f)]ÇmS3ðMÍÀã¶©xÌ¦¦ó¸cj:¥©À·MQÉû8Ø¤njø^AÔi´ó ³3R^ÿR¼`)Ì)ÏÍe%Y17BIIÁU¤qÂ×FPðâì¥_# å^>µ45]>¯è sD¿qÂ^ö§¢?îÐ±1\X±7c7õWì6Þì¥±vÊåèOÅæC]alâ`éøÀÝ+iné[>!124yPÖþ5Õà\äëfJ%-
Îc"¯
N«ÿ¶Áa:êrí®ÁÙcNYp,Jip:;§ó§ó¸cp:ap{§ñX¾cp4÷Ü ÓjÌ3K-Q«·ïôÖÎÝï?û¬¶ïOòô}á.OÝöÃFoä3qØOªy-Ñ§}a  TFsRZååEÜÐËOüéÁÏ¬¼åØÊéÞò¢é\Î¼ÍÇYö$&ì\A¡O (0@ùêD<ì ¯iynQTÀY ûhÑkÏtñÅcìµFÎY~æ\Å).BÏ{YO)(ß 46ÊÀ½	³¨<x>¶½¨÷(}ZT¹%£=-s\tDbf÷¸ H·íÍT3À=ÿ¨+úÏoI
6_³åùä¢`&br"Y
5uAçÕÇs$·5ª:6¯p<Õ«d¼ÌàTKi+¤i)ÍeÑ°¶Â¨BÏY8Þ6ýT	³ký§¯ +Í÷B8á?îðûZZyf
~Ø%]â¯/kiKÊ÷DúZÁ·äìTç¨&j³¶+éXüñº#s/ï	Á¹ãÊve^iBIn: Q$#k¤83àµ>r¬ÖÎ&úãýêe&iC;r>à±ú 2
]® =òÎÅ'zãÓ§Â¶"3v©Ø'áð#iÍã
;!Ûøâ¯D§d­ &;¿)+fÛÏI(rO|®»º°WÑY/åK<£íÊÇ
ËµßYùèê3áP¾ð:)u|M9l>¼`ñtìÇ¶T=>ê¸¥^Qãi×Mõ
»jÞ-_(Ø|fÞTÓS>AÛ;;óñqîUÅl3^'ÂRÚKAzEÖúAE;¸£U°Ä/Z5ªMø¶ÚäÌt±A¨
3¤§\'µ1Ì®`´´Y¹Ê2¯7Í6OÆÜáÕJbRÚaÆCóÂcÜpâ×^­øKíá¸+e7øÖhÈ¡>ÆÒ)âfáLùRÈB¯Õe¹ÙÁçäèÌÃ+ÐY;°6ØßÈfÆvåÑh+¶ÌzRWEgI:°Ï#Å#¨U«M©1¾F|ÂÐ¢Vs-Ùö
æeHru+ÓûZ7ÇÆ«ê6^ûìù<¨xòð@ÁPDÆ¬¡JfÏK~µ¢?ë|ë¬â³©K*|à3îµ»/ y7dg]L¹À:b;XUI¦¹fÿ4ý±¿Íý4ýã+e!¯Ô¿
º«Ä«þÍCõúù¾iä-ï_8v{ÅÜâ)µçS"å$÷}µçÝ¬íyãyzm°ÎM¾ ÒæÕÝrÉÓ±ÏjûrÉ4ôå.Ïñ]¬'ôúpMÇH:íÔ7`¾y.OÕQ®é¸g¼£ûßÅüù©ý)º^ :M]
7IÛ¨òSÆÇ%ðÑû©	NÑÒ¹õ*ÊzÄòFhùÎ7ØÐÏù´^[Stî®RÞI	Oñ@¿-ã¸<¬ñþ¼ìï%>QIú®}+ôÚ¾_<o*^Õ(ÙËO|îYÑ¨Íå¼|2)Dg|¦Dç¹7Ï<k6|­yE~Ú Ä~<ËÒ³Ö4DvG+F[\õ©	ëèîîþ¾và9ü·WS2Í7ú¾©®Ê$1.Ú´/ÕIÃ-¯~nËKûÞ©q³w^l b4ßïÆÿ®hø¶/3Á4|.`ÜµNN\X^äAËyò)`îQ¹#º3U<µ °<Ïl#Âð
6³ +¸Y*(fÅh6PhÌFç$-à¦ª½\Ñ -ò®'Åç±»Îô£}ú(¼ïOüþU¸(aE78-x?nÊwø³sBåHºPayB/«ó;OéïÕW'ß½íÜ,|¯'H6´ä7ËSw.J3×Ð¸k·Ìà=`©í»Vz{aÆÏ^_k{.´½ª,¹£íNY¦Ú~-ßáÏÚ~BB\©fªÐ:ÅAÑ?ÝÕ<Ñ7®fÞqåçJ$ÙJÛ×òt#¸â!Må%ýÊ}K·¬Mû¦Ç~VöÚékeïÐ²W-%5Ð¿CµBZu®#ñÿùßrÊþ»ÿO÷o1^ÿß¢Öº×.Z¦Ïx¸ûª_×¦ÑèHpþÛj*~QÎ<NùÞÜâ:Uù=såÁ®åxÄy<é¿ñ$Ø¯lx|yÛ{MxÀvHÈ|`õ/ùWÃ{f½·ÎOt×§±+ÅÌµßÌ
ãTu	èÁûÕ§oNýés.«o,¹,U\Â~Ju-y:ò-ÙÖË,ëOJl<oG=I	OäMk°kk>
ë´:5ÖTÑx¯Vøxë?c	§4ÎÇeãLÿMõÞÇ>®áVÀe³ßC{¿Dû3ÍýbY¼Za¸Ù+ÿúdþâÑý_Á¯#A ×CáÚû%Ú®øFW¼\ÐÀ¦øtÒÖÝwOÆ{§?ýä]&«k6Äßg±á¾Fµ³TËÒ5î4¯°Y¯.o°Ù»¿ ÊO4À{i;ðæ&j|TCÓÝ¥výÃ¬']WÖ;ª]ùJwmóÞ¼Tåª|jÓü÷mUæ£Í#{\®q'U>`s¯_(Þ,a·fvànjòÅ·ü´ÝÔã¹º'H´Ðâ'É17M#&ÞB/Ä³:]íØ\¼RàM|;ãË~]joIí=hþû{@9·Ú5>¾øxúTd(9X>ÏÌûUQÕðAÖûl|ÇÕw§\øF\êIWYfdHåD{Õ
 º4m¨è!ü*ªa÷»ÂâXÇ§ {=7Q3b`¢ò°÷¿
ªzNÉçU¸äÒ¹ÜC5VØßWôz\X®ï4´ó8Åû³öh«>»¤i×Ü'/EÿZÍZÈ
O<°~áhÑ¸o©eW6´þA¿N³DTk5ÁäÈ#.«9_ì²ßtÃ®æCÉú^~â?²å:ãg>uÚ.
ºâ4v3¤8xLÄ\%ü¸æðZµ:yÈ)!¾E©"£Û8´upN2>ãóUv<ÐshãÄO3S`_¦&¤@Ø¡¦à:ì$)ï=1Ç¸ 4åñw½%Â·®Ä78æfË²?µçon!#A£Ù\(÷ê¼$GÝVÆ&ÂHQµ©ôÿà{YÁYÏ1»AK§ªx9GÃGy:Ò2ó×æ«ßÚ$GË´^I|ËiÏv2:&TÒ²Ë}Ã·íh|à±;9a÷x7dñq±$íÙ\×QrÌ¼EóÀÏÓ`y8ZùUV¹dxkóUV1gXlRèGcJ¯£WèùNIï¹I9XÄÕ£Õ9ñºÛNFð
qÁÑ)¦[4óL®+ûªÏX³3 &»D¾ûèfaÈÅE±óV7ôf'ñc/¯v2GR/?NÙ«n'Ñì]Â¡såYåÈD9Ô7(ærÍ¼F*Iæ©Ó©o¿
¬á2ýêÐ É¹»!ö.ÉÑðn¦M#SÏ¥º´A'´L­¢X¦¿dO»Àd¹JÂ*>îjÆâ<7^ <¯"/ |\/jr«æ|ïÛ¥­pB:Yf
Gl¨j"dr.¸»YaÎÙb*Ô~TÞ@È¤½T+cJX¡eU°ªÒMUÖ£WWUyÏ>¢ÅnpÄ@	ÅÐÎÍÜÎ½Fe(T_M;³øêé'DUå)Fg1ÃÌi§ÁÊB¸.yìëéã`áuys&´9Ç¾	yäP'i§,È¦¾Âô'fX\ 2,ÍamvÆ)©áLÕÈ$B[áCTä42í° xËôüRm`: x½ßm~aeókLÑl~`Bÿ4e+6ïdÕCZnMä#s9fx«£íå§ÿÑÆI|qË4¡ Ú&QÅÃðbxMÎÙùa<Yï$2ü>-½KÂ}Î$HÈç÷Ôb´ª^áqLÝqzû°Äc´Ökpm0(7ßÁM÷¾¹] wPH	Þ;¯ó OùG nÍ×u	oÜ+¹èfÓ?]2x·HO`´õ¡%ð&Õyo¢gÐ2p éRËrEIô)ÑPð'³s7ãPB`µÒ½<sjd0Áæ(È.jM"|q#ÉV¦ºéJÆ«A3Ô_i¦dßÔÊ¤$WLp[BoafÍ¢F¨¹ãc33á<ÉXÜ )ØpDÝXM Rc áÅM4>b_bM+Â¥ó)yb Æ²  Þ¥è¥vèãÒÏî2
<{.µYIÈ¥kævo22åÃ³(<´Ô¢CxE^Ó®'
®±B^¼â	óðaoyEÉSà¡¢.>Â\5þN62òUàd;#æa!#ëaïéÍ(/etY4¬Êyvùèh¡<°®¶pý÷cHi;vÌ)bÂ¢Ff14n6NMlyñCÉ.G§Är[]¡geb4Ð0ÿ³û¡bÎJ£¥mfÏg3 eÅÍ<f_MùÎLÏ
t¸p|M¹¡×:ó¡
6|I®+ø®2sßYÒZùÃ¦ü¬Æk¦S¢:UÀdiúðYà"çùóçÃÜØ%­¤LË¡ÀÊ¹2LÕÃ£J%XiÕów#ËÌW`¥Ô0kpÛ·RÏà>YÆ[I-mj0ÑÄ¦Yã;ÞÊì¸µ%´àhIg¨Kc2st¼E Ìcy';-¾¡\ûo¯´¦CyzhÃyt ¼áhR«\OÖ=&àcàP×6Åw /Í·åuê0XsG¬À
S¤cyµrÐ¶yÐðy+æ{¼!äÛL#aù}ðûùtF°RÑOUÚmÓu¢e7%±×4Vj ñ¼á¥F	c'D8Ýbô/ïI)d²Ôé
ü8¾¾"Ü=Rseó/ã ü}==R¦AÌº¦LÕNnÔ%®ÿ>ðõ¨È/ êëâ9c¦æòuõµ]t*à¡S´%ÂnD~«|¿¯,yÎ=¼Ó[Å´ë#/ÍKJYé=aäS ^è¼zf2£åCD>îôAÀ«>HÂúÀ$/Üý-}°6zøS³>lëZéßÌÖÑ%¡à¶>p[æ} 8øÌ¸].-F_÷¿Lº-lèNÑL8ÖÞ¦4RoªEôæöÑÑTàÞè±üýYnrIV?ãGxÞÖåÏzû7rb@ù¤íymoQ¿±ßìàcËüæ}4¥/dÔr+S(Ý¾<ô
bjÚ'_àÉoùzûlÏ5Õâ³\ñ7øûQèÃ*©¸]¯/4¬è^>AËú|5íM±AßÕ}MzD¾qïÝÀPyÚ]ý
YC¿X©¯'rÍ
¯¹Ö¹-À-îÊdgé«oð¿¢gä¼ñ³¡#
ä¢ÏÒmW0sÄuð8)Ua_ÌåæT/ñBAÐ
ïòfl¹Ì¬ÈCË­i-w"!kù,d¶ïOòÄ²ÄSû¯´¥&ý"VY¹6]ÓOÞZ¡Ál¼¯àKfA¡5ßs×Á¶ 4
pÌì«|æº*åSÑwèbTh
§40jl
3æF~*ZEwAªï¹¶ïN­Ôµp[ÙtÙºmk xõ|]W0ï{Ìct+4yíJàòñèÎ¼ÒXÊK²
P)ÑñÍ$AáNÁ,gÜLbÊ­E# í¸Ãüéø]ÒÆ_Ä!×OlÍâtÀq¸4Ç·æ.@hHþMC I&
I%fæ£E£±AåcçýhU½¢eÁ ¸m	Ð/¾dLÁÇö!^S×ÚÃh+:ÒÞXdF.J£àù{rìÏ .ñäÓ«X®ô>¦ä	Å¶:5ÌÖ,íÛ3×¾Õ	ø±¯R¾'ð^ xfÐ?äÛÆÝRda¾`W¬¬3èÞzÀw:¨AÆ5A¶õ(øÊyäáyß¶ówõEJ#Ô¼lýV`ý¾hùýUmÈ\¦¦°,]½ëDþ<MvúXÌ=è2·,sºÂ¼»;,]¿M¼õªìÃsMYAé.ßÃÔÛ	±¨g«záÎ8ññ\ß­8×ºÅâasxÝB8fÙB4 yÊ(øb¯áÀÀg.Üj£G=Å=uSð`
×ï½hEBDÕí6©ÑY_}Ø£ÌâÎòÃ"O$&%?¬!4Ãq"<ì ¿Ü,ÆbÉº:Ï¯Ì^c5jðÈN¿eRf@9oJOW¦£þôÇºt'¾ Ð9ªÄWµw4ôÓ¬è2vAç#Sypv^À×%kÖ­<©ìãédMÌ;àódvv'C`þ`LnðU¿ ÃG3Îá ¸]«ËÆz×ÞXkàêá´^yc!4PhØów>0*ëÎXÜqÆS9¦ÍÎwg MáZyÊ£/ÖqÖ¾Xk¬¢G[ÛpçEÉ¡fÊcJ&% w-ÿ¸ÃcÊËÃ÷ãRÈÍÿ	K¨?%­úBG¯ <eñÂ>ªÇíêLjOçËÛºVKÆjÅÕq9d:àÎ±|[ÊNCæ/£<µ	o®a|ìØÚ¦10úÒÀp#åÅåkÔb4
Û£É÷¥MX&G¡ªd2±Ô0n.0¡n×JZå¯)òÞÁ<å6«ci¼Îlµjë¯jü4ìÇ¶GÎ8.ÛÜà´©¾a*ybËçÛÃ¾Ä;gzG¹D»*'®pGe=âeâªu{ôl4f®û¥üªç§\æ&jÖ)òà¬pG+=¸-×ÞqoÚpt-_¾°í³zÇ§<lÀ"æDÜÐ§â?îð¥î©Ë!#<MÏp#BÃo6üÀ¿Xß
[]½8ËÓ0oëÚ¤Vs§êi&1æ¦ý-¹cíÊÄïêZx°ÍÇrþ4^:àÎxÙ+>Ø0øÚbè§É³Å^Êè3ÀÑky]:9ÛxÔ¬Ôàõi;£	Aæ[ÑÔ?ä¦O\'Ö+bã¶°â$/)2ßY7VÛæýiØ;l9Ä`rÜeç¦yïØO1ï£Õ®w¿UèP×p[]-0e]}6ïP3[
ê~)?å	æÝ£)Q®§Û
i¶\@Z*ä î-ónygQ±òAPM©ß8lð?»V[ðÏZÝðOÑ}©Å¶b~bØÝyi#Ñ§cä@ñ³y·½¢ÔrÞGe¸ÐÊßåp·Ð?Å/Oø|oWÍøZuÍ¸  ´VI±'(x0Î9)«ÖåÅ=)¬6\IR
ØÏÀã'Âr¡2X	ðz¹H#Xy«²
äãßÚIåU	I²Å(Ìp@À,I
Ø)¾b*)2OîùEBû.:ÊðãJ*	Hi¼¯ £K
æ».w\f C¬+ðñÔ÷z¦¡i°æ¾+ÌàQÙ£±ßVyHÓÑ8äú
G/V§þ_U¨NI¦ãè!ªçöy(÷ñi=ü]å]9(s¥xÜQ'öà
XÒ¡7¹ýt9±¤|yO¤¯á[Õ)ßW@Ê¹èNT
L69PÛf
xo/«^Ä_^¤Ñ9J4ãolAÂ×ÎÞ;=®dî/Ù/±ÇÇmÍÁkØt9¹l§cety×÷ÂãD\òËùDaNìzÃ	
ÅWUâf}kÏâØ$	x
ÓcË*ÖBc	¯¸Fp.Âòèq8Ö®Ð÷|xãÏÜ_ùÔ«`}ª|{3>à;àaøZ÷ª´ü¨ðÅã¼«êúã}rkÂczhU·~âòÎRðúbEèWâUCï| D¼v¯BÇ9×ÕÇ#ã6gÉ¼î1çÂ¾u+ætb¸!wÀæ¼cudþR¶£qÞÆvç¥ÅtlD¶|ª}G×L[vÂ{K46ÜjÌ\1©!Ö ¡|Ül§Êc÷Ê/8BÑ^¶É0W6
Á= #Ûðrz±l,Kåã¡.&|ê§#BÅÔ?XÜoC jC<óñ;I1ªÝLs÷oßQÅüÈöâ&t?ýE[.@/æ	Î*\$=ÛûÓDNQ¢|¹WN>ÂÞWô¥½Y¯>à;àaxÜ [¹õCÙ­È7ú?
ýq>Íq0FåÆpAWeùÖS¢w>
 #¶­Sæ±é8úyÜc^ç É¼úyÜaÞæÎüå¢©Jc;9 s@Åç8úÞ[¢àª|«%êp"Ö9 öýRù0\p¢
s@Lý<,÷jÏ+àåT­ÖØXÊ§ÍsÝw q2"iøZÜ7VÊ:nYúã}2DÉ<d0.è÷
ÏÃR1óã°» n
»@g*£Óåqù­dÞ wó¢KaÅ<b^Â4*wÀ=æÎâÂ`\Úñ®Á(À==·6ås¦\
Æ¸i0(¬ád%Fù³ó{qÁ¹e0Òà×\à-@óf&±ýE2dÖ
ð¤5ÛØsØ².ÕBçÎBã!×:§¬!tFÃÂË¶*<wWë8)df½®ãÄò´9óVNÄ[ÊeT!XoÊÌ¬
3SÉÉ$&%´Hð:êÏWóXdÉHÏU:æ³ÁiV &P^²_bOmnxùoÖWìÞyÁ cnw^gtbwÝFä¸[Ý6ð:äËuà©-Pj$¨@'Á¹6W+ô=ÍGãÌç$ú* (|[3±ä ú-/6ýnä?ûq}¯JfÎm7¢lö./Øug»wIìòNm°ûÑ¶Lóç¯ë(WÀuÈÌÃRÑ9ËPôLi9×ýR~J]»Øí]n.(Eg}\E¹n¬¢òð¼WQêw÷Ö*IAL¼q ¥nWY? ^îA/6 4÷¡6ÃØºJÖ}`ä¨Ú*bôi!ë`åà¥®Uøþ¸C#ó8Tch[¤Üþ\ÐÐÛ²ê~sµãôå ê;Ê&Ó5óÊg*ç«àÛÅí;ÌëæÅÊçµeêÊçù©Â-Ø!ç»ÍàÌ
Sàåfnc¶îYD ñØV¸¬|ßa3ó!]+1Ýqîi¬ÄìbìÉùËÜórz±ÜÓ÷ÜÀ§~z,	ûÆÇk:~Ö5ÂäÐÚÎ[Ý>.(f@ïg@ÉïoóE
©?RpôVÚä©aË\^Ì~eVf7gXAHZ¢|ÙÍ¹®{Jô¥£o¦VøñúÎá(ºFxÀC'|Ø¾VÇ8QqA1u¦ºÈàmò©/÷ À¬ÜS´yiÍFîë»Ï_ÚNXà^~­ ¼Ø nÍ7ÿ4ôÇú4=ÓJ±Zã~zÑÔ5)skzÒåèvä	ÝÞø4U
ýæ×FYy,3÷vE¬ÓA	L£Y0¯{Ì½ö:mÌ_ÊOå/qøT9=âS¾cÓáAéTVlx¶ßrºÑLYÅéòD¬ÓSù.¶ÍÓÓGèÂ0=é9óä`iÓÓf°´GdªuÜÈmª|0Ó6À×
à)­¤¬3àÖXý§¡?îÐ§¡M³ìIC{DÂÐ6Ì¥ AFß¸9úh«Ö¡MPª¥¢yÜcÎ¤sæØkö¥eîmÖ©A?­6áP#^­¸íàÖÐ._ÛÚ¬JÏdÛ¡}Á¹5´Íà]FÉM@u¦1½ýE2O;àI«g·±çE¯ÄgÊi±zÖåyÊêYc4®µ¥O¢Ëû"óê#ï_¥0,7ÚéáaimsÖøéØgµí¬ñ,O;k<ÈSAæy[¸lo§P®ÞàsyÏk±Þ·¤À| "Ï;~r1ù-øa^Ðña.'/ïôõkÅ|®S¾¯ Z£; Qðj¯ È%kòUüàµs|àE-/×æSå»as{Ôö\ÎJ#Üá|bÒåL>ÇlÈèÄî¡ %NC%àµãzZ-:Çóñ`Ù.Ð÷|Sèç¸^tèËh,ÍpVøDéàk¤·rÄGMëËÿ$ìÇöä¡s<-£GÝjôí7òÕè£çiÝoÒ«`î¸ÕÒL³ÏuÕ¤]Îðû\C¬v<rÛïçº_ÊOñQ¹Ø?E®FpÙ4äAÔ
ÿ n
Y,×8®FâÊürBF+Tà-Àñ^E /÷ Àf\6>.øà>µ<ßã&aåg
Ûx£N j;zë}iA1ªÑ«X¶ô3ë]¥MZüåSô«\^ÌÆ-Ö³­ã«$Ñ>ÀÜ'Ó=¾´i0µÂb;wcçU8^UÂàë	°àhI! U¦Õb¡)j%|¹Õl|ªrW>Oäé3óðµK$ù¬àÆìðOÂ~ÜaOsI.Iå-µ²K.Øð2qc.)Ø|ÇFwìa.ÉçÃ".4p¯¥ãÀµ&"çª2^Ä:är+{ùkÝ/å§D8IºKD®¼WQ+ÃÝÑgä´w"sÉE\ÑSÃ\{wQà[fæ¦á¶ýE2=ÏlÀêÛØ³/,±ìÊ©îò<Å©nF§úe1J³J  cPÂæµÑ¥æh%)¥ÊÍÌï:èÐ(àW¨l$lpL\2)2W8|ZPÀ×Ïnº^Ì ²Æ.³ ò"jÌ
ñ\NÊ¼á6Íyþd
@±)Emaâyï,)¬Ï)JHÔ&/g%)ÐªìLmÑµ-rpLG[ðì¹k\@hølÞ¢=ÓÑ{/ÛÃóù`;¸NËHÌg)k@ZÑ´¸À
ÆÁJQ²í-)àf{>¬$Ú£ÒÉ1Ï]µG
|lnèí6XûJLN·ë6OÇ>«íë6Ó÷u.O[·ñýF+SÊiæµÛßÏN¢OwÄâaÚrt	ÊaqÊ¥ÃZ^ñ3<ØlDyä»õ<Ç>'¾ÄÍ³*ù£Äå5xQw(ÉH+ntEÉIó4*ó
Í<ETÃ
xØ*+&
å[×b4
LÅÌm,+Ï§h«PYaÊI]±øòèAiY9àT,«PY1Ó[Qpé·ìÃ

oÊÅ¢rLuÁå¼`5*+g,fÇE³4iñå|mÚÅ»ï°ïÊ
Ó T4I
ÏÎó¢ïËVë¼¿.)f0ªèÃÐXeî,Ô'*Â[ÔAë}ìåL¬0sMAR$C7ÏJ
¨ìsæäYÁ¾ z$I
æ	QN.àë¸ª£<äfªä)9f ¦
øf²­0nèË¶ÂßFN´âøxZê6 ÚÞèítÔµ@¨V>
f»<vnZù Ñg+ÏÛðÉ¸¡}¢RSgé¤ç«¼Höy';ËÝ#¿²KLÌ3_#;+Yæe¬·*Kk)o þb0SCyF øPgpò{éa§uwU1'åý¥ÐFû¯DöAÕòÛ`-D-Úäaoù®ÎÇ Uù@Ã¾ß_ö²;0@Ã*NôÌÆ¸åaØb&ÑD^õ·]^ÈJæ/É7Þ|d  &Ô`³Ô|[~ÔÜ/×pcWÂiìä«þRsËBØ8o¼Úµj¾}(	9Ç°%<e²$XÒVTx£OABT®:njÄ´S:91	ËöÑè?'Õf·#@ÀGøcñ¢è/Ohûõ¹Îh1ÿé@<!
Èi]BqQ«IMnõqf§V#Es3GnðÅqWõå¥]ö)b"íå¥qöÈO±§[//Ý»lñ=ÑkyçÇÚ­¼ÚxqV]ÙCAµÝ&' è43raB4VÎ,ãðJØ-µ	ø¢F­
a_#smójAÈk£exCèÑ¥¼'$ÑKì[BIåÀ¥±i/UàS÷#Ü\föÞr|2Ènbü³ï§ÝÓÐâe« Ö°ÅúmúO3g"D([d¶ýÀì
Lò/DáUn{I :ú`°hG½h¦ÏMåH_Ù:¡ÿ³Lsàü¨}¬ÃwQ¸h.Ã÷ò^¬âÆÿHØìH`
39¬aJq;è» ]²$äA|låaï#½ëc£2a!£ËÐ·_Ý,×L8V&
mÍ¨nÒøeyã°=	t¨½gj¸Þ¼`¾b.ç%Æ»^~ÊË#sºWË'¡mûÄ³af2T¢´É;2¨43=!ßvÑh+ÁÇ\àzøcL§²vüàn:¹q`(xãhèöÁX7ÀÃP«bµÚ~µåã¢1IB[?®Ë8Ã»ùÆØdÜ'dtÜÌ¯â>½-!S$«±0-Ãìv|f[²©·_JöZ+E¼E¤Îñ¥DQ[¥CpÎhY&®¹°%Ì\Y´J¦º	0|ÀlU§iu£*:1oG¾^TXzèýo4^0/¹Ë±ó¤¹(±¥Ä{1e0¯,î£eñwF±%L\t´BDÿÞFÑ¼»ä.Ø/¸A²xóøpêËáxn2iH)®ªÏC³ntë1Ìäç ÃÜèÿ'îKvôÈut÷ý
ä"
T¿qÖ6»Z\àâ¾ûý>FhQú3]uÐ].4S¤DgB~äñÒÚè+wUÍ \5~Ï0ãD÷YÏô\óWtäv¦N
úëíèLôúyÓ8øiÌi÷¿jï#L+Â¦ÚfIßîá
Øîì&hÙÕ~Ü2LûàÞ¿>Ê³ËwìÆF}Í$ ã{;=b*2ð%;ÉtÜ«YÁW2lð%LÅ_d]·XÁc{¤&·é__ºÕêÁ@ã©YsDíø`d ú¾ÈÚ¤«.gBôPÀU¡gò®z&|ATQrìM½o¬	Ýù;	Û'!³¾ßa¿8P¡Äå
ÇÍß¹©ñC¼3JñÈ¢Û2­ú4±.ÂÒ	!_$LO'|W Ì[:8l
ý,{ÓvÄfè(0MngÅP¥ºÕV¦5\{î&äx.AwÕ^¼¨A~;À;G,c^|m8bÅO]íÜ¦
Nìì|O­	i(ñaÔ°òÎtSb k=°Y¹ªo'@çø¼AÆ	´\"À1EöA0ey!R|JÛib&ðaõ©ÒH
`¥v5÷Ê
`
JÃ¢$:©òXÙ]GOçVÂÞ+<6ëU8¶¤C}çY	WÀ~KB{äåÁëÓ&ÃÃÐÑ¯Ûü=ÒKû~W
]äì!Ë÷T`è{Å'U:²8ÏëfXñD@ªßðÅ¡Äý½0©®¿5Æy¡ðv4ñcBeÉ
lºé#Væ}fe ÕëÖÑ¾«¢¯ÃàÊ¦Æ»ö¡"è¼ßP@¤E^©(
:;B&h VD1ÖiV×ÈFþ¤Ìä.ÜwºpôlæPó©NéB­Â±±õ\¬é|.Ù÷æ^VNÍ]Ø£-¨:§ru^P4Æ
ía°!7ÿ`)¨60Ö/Ò
îÕú=@eâ¸þýÆ°-x×Ü¿k)sGÌS¡ã.ÇÂd®¯ø.áu¬eü§bàã@°Ógn¯cÊ ;·*Ù»VaHCÃå°>NÖT9%db^°z
\Ñ¬^"^©
k#ß>³2ÔÁ¾p	NøsBAP'{2¾atÙý½`¼ª­/ôõ¯ÝÛÜÙez¤Øñýíð½-;$¿ã»èÌ«É*íìâN¢bbë f\Ò Ë5o £[7¬<ÃùÝ£6æÔ#dKr¹0jFd3yÆEµFféØ\îkq.± ²e±¬ ?ØÛ.ªH°Àµgû@¹¬è 	+exª#Ö2DrÉÆ°¤'BW°x®ªóûa+a(øzH]`ÀÒ¬)ÂL*z
½àcÜ9õó»0-É;käBZ!gG7âuÜ" ní§×i)z/ÓsA5¶9=Ê¾vBê7Þz]ëÀ«È'ÐOY÷è¿Nè÷ü×©ó·)bÓ§
4t]^[mÄ(|å1}b^´=Ø 5ðÊág¹<+Ã8-ñÚåT,9Í|¼fRyºc^h£Ä
óxÍGÖ||:1ÿ¦»1p¢Ô2u#Ã	Ý%õÙg½ìS76ö%F7c#A¤Ùí +¡kN=!`£ïÓ¨7ûGá¦=.XpÔ<û6uj¥EôÊ7ó´f±è¡*d¦WËg ^/zl90ß,z&«¢CÆiíòìNªÖ^uôE?HOÏB;IVkx½è¿þë®=æÑøÅE?Ð¿²è±ù@w^-ðréÀ¨nÖz]b>§ëtT1Oç"	7Ì1¢<ÅÖÌ'Àkæø/§§¸ú¦»ñ+u².½sÝ.~¼^ôâÁËÕZôla:Ã~Âyµèi/éª½ðB¾I£jÃx:
hG!O«ý¬o©y3|ÿª\ºèÞÊÀ~¯Î½Ù n.ùÌ}O¢&EQfÂ»CÕCÇUd*9÷d·<YXÈJæ1áZZhd¼Ê)?Ñ5uÑÑ½leHÏLJÍ#öÀëvlyFK¿ÞW 
^¿3<
=/É=IW=c7ä6S¼fOó@QGÓ»Y7	kB×	#ï©yh¸ò¦ËÄÑIÔýaEÂå¹"éP/keOÈÈ[r$³öÕ $F^°§«7AÎÌW\È×aû
º3ÍÂAÃº»¶æ=yÕp·ì	¶2_y;s¶Nü¬;æÍ$¤'C-UÒoT¼%>òYOÑ"l=R¤¢{
8L!³©
y¸læ=Ç,çÑzNB_äÌZ«ÒÇjZ2rç¯3KûDÞ¸Bó;JOï^õÌàg¼MGÂ¹Y«^8®­ôæÉKM³ÛdÀwûGæµÒ-»#Ópð¾æFnES¶É±qûþq3ï%wî	}«Ím1ZÍò¨YâCf3y?AÆòQÂíý$»fÐèëiÕxôÕP®jòÑW÷ÔrJçFl?ÖNÖ¾Gçä¬uaÇóÎÙ1¢Ó³®ì0YøfL}G§r·ÃLöYÎ?öå0#¨¥ã¢"ÄÂËr:w"LôÑ
f WV+G¾
î\>u9mÔø#!Oßåþx­
E,åã
Hýª.^Ü¢A*\BºÉuBF½Ðñ¡Q´z8GObtçtanC#¤²Ü(#ó¢=²¼à®|$Ìl3%SÁÆÍx;À/gGÈ[yv?ÝÊÈs~#§J°#×Çw¼ %}È÷eÝ!ú»ÿ?i®¾+ f$ÆgoÙÙ|å=Çx¿ûï`^ÎËÿ]>) jPðÜ­Fm<j ÿèÔSA&L¤×!V
$]O-ØðÇ	þb¼Ìüúj¼ð¹Ìf¼eíÆ|° ³êÎ	ðº;±ºCöÖÆ
)Lë¿0  i3­çQ@PåèS1ïOÃ4Gñ9
&¹õ]À'Ì=¸ÇT~ë hVf®Ø^ø³¼qT!ó**Â÷<'ÌÊëÂ¦ ZÐ÷·äp¿%fäGÝÚFÐ}I{°-Í7Wª¦~¥Oô ^?¥£_M¼n½ßßßÿ}}wIt¥mßß×ï©së$x1HLÑ×â¼7ºÄÚÄKÍc¥ZQøØ»¸½.øQêæO÷È*×&8}çtßh@g­VÜnÞ ¯»oÛßLýÓA|ÏÉÃÕµ\iÆ¯Õa'¾uRÕOªMW]suZ©PËRHUÊWeÈF!qµòfÔ»àõacåí tVê`°ºÉðIéSnrÚ *J©º:Ê8Q®Õ«y¼ï
ðIu-¾R+Ó#ÀpQ9d]jvÏ
Î3o ¦¨ç¡Öb·
À";ÄIx»YðI+¥¶k F¦Ü¬t&«ªªã£ç£ª3^Wk1Bâfø¢Zjåõªî6f2|'«)$ÃÚ4·ÞNÆ×¸Ño
÷âý´?W¥\²[åO¼b	ÓµX9x\QâJx]FL9É
èJñ%ÔM­2#[åM­``21h¼Î¯f]ë1xrq½RÔ
?h
:vc"j
O7£¤-ïÙ­øÈ(Ó~Ü[xßßÛnù»ÜÇ
¶ýO}Pg©íS÷×i:dø¸Qz¤6~0îJXñðM_©{'ñ.®ùïôPã¯pKq×ÍÑþëÞc³b'ðÛ¬é£ÖçØÄ{ÊoUéß78CÃCh.ó	ït³6ÈËbNúñ6Ú/îy£aHrqDjFQ81l(èvGõ}­ç |ROtetÕM=+õ#¬¥ÖëèJ×³fN.Ý2')jJPõT ú{îGÛ÷Uô_'ôçÕ©½k­©_Î*ÐÐÑÁmßÑU÷~²{/äÄ²¼Ø½2 u±e)»W3Ë3µÆç]`±a;MÔ+9Æ=l·Ý+Ò¹¡]ªé­7@`¦ ¤KOJæl£ä^jAâÖR%lÝÚà£1[(ha²ëV#3<Füûh+CëÌL^©iÔ*)uÞ	)Í3uÃaì·³O5Ûí¬½²y:ôiÃ)ÂÚo`ßÅ«xYøÃ*õiWèvC-ð"ö46¾¶Ó6¢rõ±ÿødÊ_'¥ÏÅåÇ/`* ;_»A¹|åFÆÇ²®$Eô¡´êiti²Ê¯îLÆ«hHNv·FÔQ£øYã%ÀÈ©úPs¤jÂ+xåJ(¯Æ¶ðÎQo:à]]Õ}äØ&7A)äm;¡{¦üè¸áoD¸§Ë©­ºVðI­<
ØæZ½u÷Æ9Ý>ÑCX³¥]Ä ½ºãdGaÊ«Eyíêi°%&³"\²41O¡xqrÖnðé¸nleúä¸©j-w*ÑaºdôãõJt	Mz(ô`{$¤ë$'ñªÃú×yWçÕïäÔ¶¼wôofÍx;ZQ®5çãÈò¥TuEc>:-jzUø`>9¯ëÈ·â£ ±wÊQöë¦òIXùÖé0½/bF/ÏW6^7Yá
è'ªó¦B·úâÑÛ{~)³þ{È°²í`ÓwÞFç¨¿3Á;aùôâÿí{º]}Jª~WjU×qùß|*«Ã	±àó}.Êpàøn«ZÅ¯?³Â3ðXÜP> r/AQDL^(
î"V
ºgîå+Eå»å¤Ë`ôNËð«+cÉt¬ÒÔ­³®ÓÂ2êßQÀâ«qà¼Á´Õ¬ø~Ó¸
|³­»DD#Öx;:«JÓ/À»0î©ÝæÂ¸Ä­gqì¼}N¥ªJe@QRí{I÷õÞÖ4£ÇB3·W7M¤I`I ^ytØ¨ÑW%3ÊC"g7JÿÎØñÖï÷÷À>K[Æ¶ñúRd<'°Ü&2Xè¾Qq	&ªZ07Å÷üNón\ªúðêw)(tBzÉÈæ´ÂZõJÝ¾ÐG²$lRXmt,dÉnªÊ¯¶«ZèvóÕSUùZ\üWBgù8§SbW2«"ÀJf¾A9D© 	#u¬jÌ10À¶&ÄîÕ·¼sÌÉQ7m`úôâ¦Y}ìF`zê:Î@SÍÑÍ¬Úðv4 \zìN¾0Ó85BÏy(õlS2æÁòùy¬ÊÂMôT"@åK.ÞËã»Ó 2Peü9EØJBn!\Cä®äì¯ÝlË§þh¼¡·2"ÏêÄ¥(ÚgexV²ê)v²Ôee®©´ö=¬;]Ìc|¬£õÊ®`ëU¡ØñÙ 4ì*vs±§v4/3©XfUÍQÝ*<dæQ&ÓÊ±wÌÔï}£è¢<6bÆ?EñMg@=h%ò¸Æ©.8Ì01MÇÆA3r&,¼<òr	é!í`À8¢¿ùòS2©K=ÏGÆ,2|9µv7ÃR¡Ñ Æo1»¢"ô«;ÖIÇ D&òùHi=%ðãõHqd&>·?qäg²U½©9ä
GHÍ`­cÂ1óEWµec_!ÊVÂ°qÐF2ÚÙÎ1rìÈ]E¬3VyýÎ°¯Ü6ïÉ0Bß?îÊô¶ÚÈõ¤ñáÒV7­%ÃÒtx`¾D$0WñjDæoû~ã·pýû·Fr:¹R$2*^@]®rñ/ßK
|ï¢ø0þ,È¯ø$µW¬õêsÐm£eùø¸ ¡¹½ÏS·&e(cÂ¸!ØZBÀyyz]5újXT6¼£\µÁwS±ã;¥MRTìW²JÝ|"«|XÂYÈj©N4?Èl#³&Çkó[è9³ÒÇõPÃY¦`ëÀ'6~`&lu-L©îkQÐõ%®®"ê¢	GÂÅMÔ%Å¾W°!ï\QÎÐñ¾S0
ÃjB>H¡úÐKlèMhëý¹Ädx=l4ab¥´Î³ÑÆÝ8Ì½Ú ¡éÓ#äi³k_OW]Yº¡g icáÃ\'tô8yEÄãÌçµ8
Ç dFqÞT%ºBÔ#Ç.*VÂð}'<ÇãÕ×¿/ôa>ñ¸üFoedÃ¬ävbãð¡´ºaQ¹[É¬g¼ÌqSõ]20IM¢¬w¤>ñþ¿L¥µE!Ì#ÅøÐaaÆúÇ"£Ð§ä<ïÂö5¡*kªHrë2Þ	!³®mÀè0Yæ-rlC¡ç«ÙT/<D[;T¥Ib©¡Ú¦;:CóÓë±}×GBÆwD¢©'md8Góµ×PÃêt|SÃ<ò2$
=LÍ#0­±rÎîs	JcÙT!ßtàÔ+©òÝRÐÌ×§Õµ(>EÝ
)[6ÓÚ³É®ÖS ub;ÜLkFÎÆBV&¤Ù	ün?ã´ü8dÉud½n+á
Ø
7V6W0LDN§vÞßè­ðhxåÔ÷ÓÑhû÷÷w¨Ó|K¥ðy(â«ÂÏôÖáñäk¥Ó¢Â/|à5¾¤¥fH	 KìR'úzÙW
_­¼PWÚìï"8z%Ú
&ëUe"ê
«ÈgnC9yê(^)²	°ôØAÁ¢4 F;ß^9~A×	íÅú_)t/ÑÀcÇy<+_è"¡«o÷C¹Øèì¾çüf<xüÄÐWÂMuÔPz
*
TØ¾0±	UOj9e:µUÃGjì"(j,½&×·¨°+ È%7>nÉÒÅN-Ù0µµßBÉO÷Â¼sÆR"F0óï¼×è¹æ}ÙÛP¤ü+6K3H7¸FÆÀ©Ldí{HTjÕw4ªì ìr*»«²dU*»¬62PÆõ©ÎÂ8äÍ§òhÞÅ`t=éD_9V.TSU5W/éÍG³¢LãÛ>8ÒÅ.¡¹&Ç`aìMÛ{M"HÌÆÜ÷½æ¢{×¥ì2r´¶#$º.!óÉr÷ZÎadò¼è»º*0nËuH:uQ}×rßV°cªaÅ±HÎM¸X»é§ª®4^rÒux;ÀGgùÂm827]LóO]	§åµ«.øÊUïG83fÃ mR	0wèEFú0ºì¤:¹sæBÎ:2æÞjÅ2¤=q+örL~ÿýéÇÜ
ÑÀ=U÷¯Çè?Íì+köÛ	ÜzÑÓ»ÚE=Â¥-«°Òa£LÜL:ô'ãèí)¥KG[b|®|dfàÌEW.w}T®D?j¸è
Û&,³rÜ1  ÀÓXÌé¼ÕÀÒpÑlÇW.IAÜö«\Íï;ØÓ¬NÝÈ3Ìp8oG02S´kBÇTöUóFÈ+Q±×6xîIgþ¹
aÄW½õNjÃR¹ÐÍY£­JûÊÌa{ªR§¹Ñ[áÒANÖï°.-e|ûÞÔëxÈ¶bmd®F5{%þs`%ë¥e&¾§@\jý	4ëïÌ!i*¦P®ðVk_;æíä+Ý¬òh?É* 1ZÀô-ÅÏfH	þÔ\j3mG&¤+Åm Å_Û»Ô*HUx 2÷½·ÈaSypËù\yÆðñQ]Éè§Ð2è(@Êùàuèáfm:L1Ó­æXt7n:°zfÓØ4OlbÇ¢* 1£3#tQ´ ÿO,ßX;=ìtÅ1|Â  PB#³C¯µf*Ls$ìaøê&L¾VÂÞ\±ë '»KÌ²^õè2Wh4Y"d.ïñèä}ßÔ1óê%Õcã`Âó´WW¥H =-µxvñªJM««èHºÕôtÑk
ÏLrì½ÔÆ®ãf2×±»ÆæÆçb^¯dÀ)§Ê3Ä%ÔÌ¨	±¼.îûÑËsÓÍÔ¦F{ÃÒås¨
pÕ°é®ÈôiíÙw²qÜíqeºäLV^ö,{Ü©¡÷ï©Ù|ïÐN³÷ïm£¯	§W+Ö3ö=¡LVdî¡*²·=¸iÏ,eËC»ãqØ±hº8÷c'yç ¶v£´Þ<tkÝãqyAñÒ*û¾~h>î-ÕTé';Þý@ïîÌKnááãé&vÔ©¦±­'Æ6×ãX Ç§Ör¡¸Ô¢òÑ$ý½çí:üu*§²ÔÈÀ§Åé"èî^kmÄ¹g«C[äÎINWPù<>[9F¾Èg±Ü'm+!u<
9&Ø|â[±V%C¹²©#ã@0*ªàÇà#ÈM½¦X4!§8s?§-:D+JúB­Ñçìô"î;5Y«Ú ífDÐÝ¾`ü×cwñØÆÆirÚ1ü?F±­éÆCë£ÂjÎÕKÒFÕ¦ .qÓ++îÄQà´!WtÛôJ	tlÞôJám`JjÞð±!{ù¯.i¼MÉÌæ:É¦V*c£3iª.sºø©ËÖ.ÅÐcÑ]ÊÇy¢ºK"ysh°hy}©×#_Rk·mYTUø'Ý¥T%¡oêØ÷æVÔfoÌ§íKîµLgCGï{í}%ÞÙ0Ø}[)v#QÔ`%Kæd8eLÏÜ k;M|¯¾COe$UCÎñ(ìPåiN^ñ»n/æ
tË0óBÞ,Ç¦
|ù¥kÎPIÕ£9B[ëïSÑ>oÚ-#ÿ¢qøäyM¨ª-/6 :)AS ¨
¡X=ú2÷Aá¡¹Òp¬°8Éõx$¬L6ôH3ú5V1 <¢n òÒX
`87·ï2Ýw{Ö3zátZ³Üý¬=ãVGÝ³0A°PîY·»ÂbcBõ¬êðvÝáÊÔï
À$·¯,rYí8	ØHÚ²é2º¾iÁ,ãI×>ÐÕø4@qÌÉyÞ¯ôx¿Y»®dXS9i ¦o( FDÕQ¨êV»ç5"EhSHÔ³§uU<smæ=~¤ç=1t¹|}À,êÇ:`!Mõ¢öÑ0kÊqQ÷mE.ïWuÇñ¼z·SÜf#ÆK»ùæ¿ËsNIQîÄÉ]W·3Þu]ûóúÌþÕÿßBÃ5â÷ÔBÄ<>tÉýM¾ßuüµ~BôséEuâtí½òñ|:gL/÷ØÑ¯æüu>{S
î©P.éôõ:^è¯ëxbf.o{¼^w>àøê$íè§þmÐgÿ6Î;ñMmxmF1ýöýÆO%Cû	ê{fÂ"Â¿rè8Or>öil-vªè_jé}J%N·¼+ï4w`ÙUu¥1úÁÇëîÃpÿÙ°mÿãSå®ßù2ÈÔ0çBÔ&ÆxÌ¦x;èÐÈD5zVâðóÄá<:­/Gçý P«ùóaü_Du»ü4.çm5Æcwàí ºÊÛ©Eè°=Ó¯«G^ÑV
Ï¯L#¹R0ËÊi
ÇØnCá`Äº¡X'FjÆææp¼jßvtè ÃM}÷âßc¿Nì¶ÓÏòéCiÞø£U¸ùÂ^(<¡0ä°»Å¥Ü}?á¯jÎ
Õ+éA¥å
}È<iJ_Ü|ºÚºùèkåf£{}¡fSOþv}ZèÏÇ¤YKèÜÆsXVð°í§^ú{ì¿G}?Ût3ÛÛ4Ó_>íAÿíÏï'üu¶ßP5!TëjÐsÊçrüºªy£TÍýD_ýRÍv³}_!w<´¨¼°ôdr?Ò{Ñ~è¤¿Ç~ìmÐ÷}@7½Í2ý¥å¼MàoJp·iU[s;ôñOþúïjvÿõÿÜý7Öëÿõ]JxÃ×§ôS5¬Ñàþ{, Õ|~~o§Á½9ÿºÜ Ï\"ÿô	c0IÑÐì²)MRü^oóÉF7.RáÓ~aüÇðþf¹¯ÎO
×_óPªç(ÿX?t»¢fÔð}Bðs×ôCÇÜcA¿8rÙ w©¸)ýâºóuäWuÛ3Ú³h¬XcúH}qÞÈÞ¼aÏÞüÖ-9lË¶ÒÐl¤kT®
>xÂo[dÁÎùØvÎr ñ?Ú©5ö×þèu
Üvûgh?·hÿ¤»ß·w'/GåiLÖÏêÿþ\	
¸_
 ýÜ¢>¿o>ßJÚ~ø>«ãg¿ý¢=6§j6ÔÏf±ÙáþL|g4,Í×=î²¯]°u^=(~la«vÿ@q¡¾(öVàË
MøZM]÷º.mè·Öùq#©Õõd}¢:}ßÍãSßüÞÛAÞNåÓSùMëÏ¯§²õÐí:_÷¸ËT¾`ë¨?(~la¯&ÍªÀ½ÉÝò¯ìå<^ûB6³øKõ>#Þj?×éôcëçÝ>ôÉï3þØëvöv½÷Z¾íãaP®½ö´ É>^­>üÍóV3¹vÛiO±£A©'Ù'®e¬$>û=\Ãäko 1åÉû¶Êò?òõF{C!Ç)ÉjNã´*e
µt2§ 5d7õrX³¿Bql7zÁ¹&ÜÊ&Ìt
çr;h:ô
I®.Ôàþ3¨æH~8æUS¬ÎèSq¯sö
]&êÏË&åª¶ø·	?:ø[PÞ&=*ßÔÓ?Ø41cVË°è×8nðÇÝzÍ½ß(e­SÎ|È|WÑ Æø£RÖ0cG&×ZAdéÛXýûBK/täÖ2zýûhÆ÷oM®j3dç=ï¦É÷)Oe¶ï¹&_§2ÛZsQÜ/;'°v)5²Ô)²l¶=/dQîØ6à6v,Ñ_½¶H±ï®ÎYSùÐ$%[²¤EÇúôFe¿óeLÅSAÁÔUsÔ¸0ð²ó£²ïê{ÁBxÉ&ù?cÓÇ©XNÈÚð¦¥|R=ÈHÉ%ºSÇWf]@)+YÍwªão2Ëû~¾Z[é ÿÙEþT?æbôÄ±ÈÔ
^hI·ÌJäÚëÐ¶$Âkðíè7æN
Ñx9F|¤Áó u6u6Ì'À'Ìsa(¤óè¬ákIÅ|¼f+ûù[QüXÝ_ÌO\ÕC¬@u 3ô«Ut¿Äü8)T\ÒU)¥0GTW0.ªo4C¸õubª"ºðN±Ià[ö7	¡ãpKl{ÈÛá{ôÊ}+v|OBÝom#Ä7&¶U®hÅy|ÍÑH!hvé&_³æÜ¿¿æÌ¨Ìÿ1qn¢ËÔÀx6k/1Ç!v½ñýgÇ>ÑÕ|Å/)òú.Ém|Ôß¦¤úêÞ×ïØ®án£®fPC¬ËÓÆ®ä¬ý®lüéQl°½©È#æÈ-]o
ÌÍum9màC:;c~KË_"¹ã?3!;F-»=W³Ý gæóÿbr{Ì§ÚrLÿ1dóU¬§Z½/,i½zþgÃï3Ð¼§¢ÆÊkÅã@0óÀ?'HlÊç ù¦X?<JÿìYçu}Ò¬·CµÞ¯Éø;©1ð`	uÒ®F=
6iªXÎ2Õñ§É¶õkNmX7¥Þ±[ûásÖh kØ
D;2Xo«éÇÇ­cÄ§LbB5ö
=ý¹
Ð!ú¡ÞF7Á.Ï#ßR$ÙäWÂj$áã"&3Z"dbËë0dÆ}*º
4¹ÓVSpzºâøvw·ê#m´' Äi	~Ä&ÎªéùÜ&tÄ¶8	Â.ë<ytR¾1eÓäÉÙacLºÖp×ÕbÐµfît)§­ÖD'J7ÕêíhÕ­¿^uïwLõxJ&lè%#3a3ÐØ£ç,_£&dKÆ7=29o*ª0È£c\³zæÃ \9z7Á0ºøcð.\îÂ1`aóá¸`jûË@ZTv ÈËºYICâÑÒe=C$> "ø¬â0ª,oBÐñ
|°øùFNcÙR¢¥¡«|Ká¯Ë¦pX½o² 1K# Iì°UûÝÚÑùÌLÝ·à¡+p[àçTUoÞgÀ7uH¸Õåé ØzÏ`ç³ßÚË5ôÚÊ£µVÜÎcÈêjdÊj°Cöñ
»
Fê¡þ@Øòh%ÙÝÆÓÆOèßoüÄµ]­þnK2)*>eúQ¡
¬¹¥8Aw
?ö. ýÄ0~MRu^ÄHYãûÍ	*+ìfë+ËH}¶ÑYe.!ÄMaâiÅÑf}fîllkapòGF+!s©kU
#vMLÁÄv'Âj¯uëÇºt"DyÖ1.ÊJè K2²Ó°0ý&@o 	r;è~vî
ü­ëà³ÄGÎyýtússGÙ¤&c¨!Ú³FËWsU¹tBæ}%dÀÕ#áÉY¿03ÊzôGÂÊ­=mÚXyT]ÎmäBÍ¥zÍyGè_>q¤S¥jBG¼bxeêÕ¦2âè®È76p
`6MrÃ`"ÆÞazöCô#võ
àzqð¸Ë·ïM2UÌF;>MËàr+ÄèE
 'ôÔIWV'>TW
]ªqCAÅß¸
)E®¯L5çÂ«óÑ¢.+,h
f¾cðUEéÙ´\à:« ·C"úVeÚÒ05°©7+à{óB=æLC´¤Ö9Z[MeiäJ=ý¶éÂ#1Õ`
½/Æ»ú.ñ)Û÷c¸·bÍöÀ)ÏtOvþ{fNúå»DZÑ£åX»C£tD&6F"§åBg%¥Ò¥~oxçÉe%±Þ T~D
ýÍ1'xC"Q3q/kuæHÎÕÑx ÔðH(nËjÍÌD9Û
XÿÌ\£8j]¢x''·+ I85£ëÀh×!Úc`ÿ@OÔ£Îp×Åá¡;
o6ÓÈ2e;¸Â­5, g¨¬êÀ4ÆO¬gÆÈð°Ð:Æ94GB[*KhBFMaÈ!
f¦ÖVmð¼ÁqFa
=ßE\Cî;;qrÞµr¤ûA-ºË¥¹{°Â¸£¹´íK­bÇ7ýXy«Sýè
uy]/²èùèeÞÔã|ôìañDY	àéÝßÌm­'C±2>Ø±DìÌ1¡ÛÐþÊqìO1Æ7yz¤¦SÊRægTè­©Á´
#Ã»î<­Ð)² ²­{2l^nÈ²ÌSÆÈ¼-Q?md÷÷~Ë]MNéÄ»9ìê?K×¹._ZªÇµ9ÒÕÒ]Cëh5ÓD²Á^nuQÑÊ3ÊSQÉBü¦qÓ¹!¨E;V:¦lÚÀÇ¡q×6¥S Y3 ódtå`B&ý;UAÔå|åªÉ-`506PëfF³µåXûÎ,uA÷Em8Rg»bº7^£ýXwX3#¾/ÖùJèx^£ÏØY»OU¶¸rÄÕèÄ1xY5Çs?övÔ3ÐIÂrLYÓ<maá}ÀTPòúÍÈTæ°ÞÌiW+}Ptå8KP]HR¹ú¬jíwÎ%ªµ§Ð`ÊEè;ÅëÙ\!FØÊãÖ«60 Y
z=oÈ]E57¾k­[»íªR>3m`¼2qCï&à@Y¿'F
LU}§TéßÞ ¹ßKîólú	yø ¾'¦ØÕüäÜÎ
FÌ1õ9!V>:5JPøÔ.©ßl Ì;b5'Q7L
¦)h£þÐ¬ Ý
yMÙ±I -@ÃÀ\#!mªKwWèäZ¨®âë®¯@À¨Â}`QÍ*s(zzÜ©ç7ìc5¹ê4E­ÌºVX÷Ø¾#Øb¡.¼²=4E]8Ó{Uz6* }\×
tTüË<P­æY^¼¡À0Ùy&6h¯;ðvL\ûuÀ»0n¨ÉçF¹TæBc2fÞè%É\Ã°Û°¾Ýï29òPzÁß)oq7¡aBqðI=5Ù"MéÏ>n^Ä`5©¦ëòkqW"0Î8³`1Çrûo.¼, PÎè91POô_ôq+åã	ú¶**x»©¾é¾}Rh	âVÐâë'ù½C2ñâ¡ÓýÑç,³t}~á$çèWI
o£9?ÖïX]-ßãUÇ?/
ïNNx!ª-ClV÷² ^Å*ñ¦ãMÕ]KÚ-¬ ûâµÏÝx÷ëÄ30âp.Ì«×[rÓa,LÓ¿7üÂ9<'Àæ»ØÍ¶ökÈ¿Èò^j$<ás!ÏÇqû²2-(
©^ÉVp`ô}åTLE²s­*£GÓÁ[S¾nÝ'2;nTjqÊ§t²ô/¦UÃm<Æ÷ýÀîCKe¸|Ç»î-UÕû­	ÄLÙÆ3o`O¦õÄ(|õù÷¬È0¨+{èF¾öf
-í[æ{ØtlÔüÿ¬ !Z¤@';¶Üq^\ Æ^Ü+iëJ Úé>¢ 2Ú½ñv
¢,Ò¥jæð	sÞWªiykD©.é&ÁÎ@²jx'À]êVå;
&U¤gb´7ý5 4	óÂ2ñÃ¦¿¸cÓ]Å,D±3ùg±z(ð+Zz'l(P6}ÄUà &ÀË»D¦>´V3ÏCÑ¨ÑI~ÝHWôË? ¯XÌbq,Pf²é¯ø¤IÙ¸í¶¿qYtñyøXªDSË«³¼tÿÂþ\ÖIµÈºôRÖí'Y·AF¬0ÄSõK§VñHö|7¥Kå~VÇÖÜ¶lfVHa«ØaK/4;æMb?+ÕîXØV·£M k'À­saì%ç¡Ýñz¡#6õîF|¬A¤¿%×¤ÒïtaEÀÇW<K:jxÎ ?bJ`9æzìP)txÜº®Lx©å}û×	[éy@)¡&?£½ý¦Ç	ã¡íT=v¬8j¬º$%½©Æk{Ì`º7^ë{
.ÛzPøF}Ï_àSÈÃ.óI¸ÃþM6
ûSI(ØiÐ+Ix@ppïª)¹Ñ­bÞ&ÔÑKg ^ojë>ðÚ]ojLì(+óø9½óR|lÇ­IÁwvªI37Õ À]jf£«TÆ2¦CocmÓ^æõ]pÄÌ³
Õ\KwX_ÌÆÖ 47ñrÂ×]_büRÉáG·_¾ÐPâ2Ó±ÑM,­AÒN%oÜf2À'
(t¨LÛÉ±$6a >c^ÅCOGg:R^; ¾à+ëdx >ÄÍrNÜ³­Í%'¨^O	ðº¹þº\ÞMo¡kÝ2´¼º|
¼º¼ÿÂþ\F»ðHYÆÐK¹Gdäa1hmØ¶âjcFãÊÍÍ¶5gúv§$.¿ÞFÖ=²ÔzôuÑ.§õ}o:Õ%ìç8wøzÜÅÄï
Êý3°*td«´ïoïÿ¾¾G*%ïïë÷ytùhlPcÔÛÚ¾'½Øíû·)¾xïU¡ãû](5Sh¡ô']Ûæõ3K¤+óÚ-P¹1+j~afáh>¶UÃõ;ÝÍ¨÷öAÂHÇ
¶aZöj¢ëµbmê5AÕ!ú;î9÷	iezq(îe: ½ÌöÈôÆ¾Ke®Cz(¬íçµ
¯ýUyÿìt¿ï¯GÎB²óúlºËbxkLÇËµ¥dë]UïÜÁROmËL"¨ û½+lìÞí
w íyÉ¯XEÇ´hZ¡?°´óÀ«\ë5+>v3eCAß­ä7 >Ä`ÂfÅ*RÊ<¯Ea2Ê¨IÁPÍ¡ån{¥ªÑl\fZÂÂAR@%¿.ßVgÊÅu3#ç+Eày|Ür·_õcB¢wúÅsÑvBK
êô}Ýi¹¯#J/b0îëvw%¢ÃZÅîÄ!E>Ð¡ÕË»HznCßC½¯úd?î¼¢ý:·÷ØnÂÏNM¬ú|]
)
ÿêÊºáãk¢ÏÆï°ÉYâþø"Þ;Âçx-ðÆ{FÖ¤nì§ð<úqÀ>_O9Kº0õµ7ÏC, d¾(¥ÍósI#,¾^¡@û±~w6ÔözCÛðB¡ãRÝîÏbûLÝ§ú®EmµryPÃçV¶#>¯ÀÇ-,£5¦Qòê¥RÐ4@£ ø¡ôÚ|xqèXå*¥ ;µì9`5`ÇÂ¶#Æxâ«ëÜwGnËØm{èHDØ}4]ØiÔ®G?YZ_`O®
=²Z¾ZA·$S µÄU½~Ó¯Áx@L0VÈÊg¡ Æúî<%xÓsDw L#!´FÍ vó)>ä)µágè÷)éÏúïòy;ñùOpD·tY2óÐ+ó<}&
lÆ@µ±sjñTrMI©££Þ8Ng¿|ÚG4^7vg3 9R$¬iuÄfDè6Ì;àæÁÂçTóoº¥1F§>N¹?ZLÓËê
[Ìø!ÖÅ,¡!sÃGe;±h^õ¸Øï@(W;00ªÂû)ùGMocPê×+à)v³çñSe¸·ÿ ¬0dÖ}ºÐî
´kØ¦#·¡4k[YºAøÃè%8m¦ù_Zù ÎàwÜ+?É¤h]G'éªÀÍü4¨4^f1¡ÑR.1$±Qd3o"GgãÜ?ûÒ®wànÀéDHÂØ$.§­
å¥Ì®qU¼$6Sq >é9&X=yË¼ÊÄùx½\þ!ýå=­ÂØMå2@,âc¹(YÖhàª#Í·ÍÞÈ¿];o¡ÙìÝâ×9hµàæÐ´)Yüà ¯$KãqÇJ2v«l4ð«¬.SÜ Ým7`×mô¯®UöÖðbÉ`Ú´f{g"Ã¹Wp4byåm6=0:;À°zPÛ;;>®"þKØ]Æµ¡vì¨>É9Ó
HØØ
e[tIÅR¢ÕHØù7ÂWÆjVÅJLwÔ")Ùòö Fìw ÉK:TuSAÁrãt#Þ*~öÒ\Å£#'q[´|gBç«]àx6WòJÁ·§èj7¬.=jbõã ØbÊüh%|?_Aúïòy;ñùOpD?¿þëXÁ (>=~âx4èÆ×P:ýðB§'¢3`X
yFB(ù |Æ<ñA×Ã³c>3;¹\¬ÖÙ{Q°ÝåíÒl0</»U&Êò0®v0Ð¿=O£Á h%yç×³
à¬ýØ(ìÏ4 "{¹°QÐôJ: 
h0ÍGÏíÄYÝæxa0o8µÍT@	<näÍ ¼P{±ºÌ×9æPòö0bXîòz @2´¦ÀÆÏGÉºðÚ`\&m¨éi0 ²vãK6Sq >é¹1ÊçèæmW;nx½\þ!ýåÂÇ¥Î¸Ýré ÏËñ¹\V¬QÀTGo½ñÁ°Cü:­LÊ¼ bq×±Ì¤þ_Ë[ISD}dVå{K==Q³¦]0[À®Û´Á0Zs2ªµ¾«¢R®<©¢ºÁ 7®òÊ`ÿ%ì.c!
rLOí[ ©b^­ ¾äàÆCË ÅK*""w³ÝNÅIU,èó4»Ì3y^NÏÃp!Ñ#o05ÀÉc20â¨_Úfècþ¨Éf R®%¸U^P	¾ÐõôóLÓ@5y:4Ì µµ;y5e´¥éýþ¸õJ)_r'hôéªüß
¯÷=ÍÁåêÆ	e|å8HéxVB#Ã®æäáÛB«o2ßÈÃbõYû¯F³Ã¸`ZÛS-¬±!³O%ñ=Ü×ÃR'*­+aà[É4Öhï v.8îô+ad^¿VuÈ·3ÿ¡Ñ©sélò¦ª%JL)UâÛ~2;BuÑÍ9ôÞ°¿üõr Êb×SÁ)$ª§ÍØE^XÖVXgéÿW?yY¬í·ÝLÐò!<1­¡ íHúÕ`+¹ Îu^Ûj<äA×2:QÝýæ´Q!ÖGëaÐÙ&Ûèªc°#Hç>
sè-±õ«Ûïiaäø":ç1Ä.ÿeÉ¬cx,l;( T£ÆpÜÝÚn¿1ï1ÜÝ~bÍQôn5º°Ã2ÙSÎôð[Çð^GØJÄÅRc$8æny.yØºTôgÜPÏbùÏ2#	*BÇ·tE/ùFHß
	¾±ò¤2ær,1TÌ£Å£Ã¥¯« Cà&ß	ÑeQ|ª*Ý2éo§½¹cíJô\»½Âc±ÅyÞ­%bðL@z¦7ë8#(ó<:`ã¦*Ý»¬q4=¼øÌ¯bðÄxã9êÞÊåÙxcÃ0wìÈ=þ2D2¨zy!sF%ÃyéTïÊIp#kßÚî>ëwáiÕRôkoêÔ7|ïfc#«ÐÀ=¼ºz`®ë¥A®ºÐ} Voêa-ú\Þtì»ÉBWdô	Mèø6úSGÙ`meBE(Òq`,ý]jØ2æ¨ÈCTXü4tCñÔ²Î Q"º:2±¢"ÄôÄ(ÄSØÔÌ»ªyöà%HÁ¾D¦&é¹îóGãdwêa¯é¿^®êÛÞ9b¾Xñ¡mwØtòá¤£¿a  hEqÓAÔåªåÐAá¬BÑDfÌ±¿sÑµâÚ·rÌ^«"ÉvèÀ±L/µ
£*5Ò1r3ÈÜË³jC/Ñ3GÈzõ£¸æÅ¸(¿£ï}â=Ò£çGU¡÷püM¾¯+`+ø+½æTÅ{ÁÏ`.\7úe"NAÙÚ÷&ºm ß·þÏ³È8{¼!\ð²Kª,ø)aL¿@ÂúÙÏª
ßÒÒË¾(hPA8+\e·¡È|¡2ÕõíèûrKwp(íÎï(Â¥2âX× Ý¯Ö×zB¿ÈâW¸Vþ$g¥ð0ARÞô^ nîª¦à±
/ÑI`
à¦nçný¸'jÅm¦Bµ
ÅÞyev÷,/ó,¯ £§3ûã_¯xDQ¼F_(Ò(n©&¿3EW#µ%`oêdjÌ¿'YÆPHÐµYhd¥d>d÷w&á	»°uÀSþpªâ1þFÂªÇÚ\	{E,Æ{ÒRLn^8
pÔc³®;.¯Îk¼¡çr<Vy#³sóÄÊ1Áv1¦92\´8>¬ÝÖªôî@5KJ^ u
3W¥øvwØûM0£Äk
O L°nìÏrE¤­!ùx*
(ñàÔNµòÙI ç­ks°ÍaÐBèmy÷Z7@¯õîä£HðÚbs­;¯*ýÔWo'@§!ÌÖvx
Î2µ!pË ú^Ö3ÈA-M&g£!@½tjxpApïò¦ªÅðxÄK¬P°¯-t!¬1c;0§½a>ß¡ç«¼VÂ>g¼sä	E(º*0í¥w7côiã¼³ÑÝí]6
¼oÞó¿êë½Äµ8Îp~ò®&êÀféGñQBÛjy?ºBL
ÐïíV@´Û7K*ºñfñ}<§Ýk!\]}66ÏÈKqß6¶¡¾Iµ¼ïè­p¬7ïñý.+²$ö<Ü²
¾ó3»ÿdØ³Â3ÜJù;Ôd1l?«*_
ÅçwÐíÂ¯O"¾¶HL¼¥	såU·|lÌîªæhÐÌ$O',AZ4¡N#½wc46mp´xæ0uBª¾+@¨°´ó«`ôaìÈc! ó!ëÉa§l>ÎH[±[kÂh A?òbH
KÌôkÜô.ïGÌTÔÊ×5~æØyÕdj9>.­@_¾b×99!ô°lkQ×ªØðýPÐUÅBm<VáæÝ@ÊÍ¨;v·ó=kwK¸&iBÝYp.3H¯&düfÎUAOÜ!ï^ýTÇÁãim(BfH-qcÂ¬aÖ>EèéñÝq=ôaô\Q½æ/¯ºu½
GõÜKP²f´Öt^WâTr¯FÏHó0Dh3þÐé<A®íNEÓÏ¼òôÍ½¡·Â[ÔöþýçýýÚ®ðéïË)µâêkAáWä2Z¾ó)(Ê¬Y¬/CÛAQ@Nt,PPÛ«
SÛÞNÆKâÜ· A´iå{Z,RrüÒ5Ï
l!3lÔØÙ¬wV±¢ÏmÑ lVû.7O¬( ¥ÓúW µ½  L þRe@pïÊ@ºå°x¥[îs3WSÐ¼aÈEè¦Wô2ZµQ ãY1µâPx³³i9Éò¼Z)UÞõj
:äó I0èâ]¤ Ìï]6µÂwÍ¦Û±0å¢¨.ÆQSÔD
¸*
ô:ifõq.Ô¡-þÅrÐÖé aÍ#G°hm2ñXèÎþã1ÿA{á*T£^ÍHß5;ÊöýÆÑöïßEbn,²râÒgÎ¬õ»1kÏZ.YsÉK;])^_Y_¦àù7³4&.v¾©C>. oFyâ<:é×1ûÂ¾áÆjFgÃ\ãñ´ïÈÑ_ÕèëÑ_ÕCMÀw«2×À;)c%Þax NèE^J¯dwÔà÷dØ*ïyW2æÙ	\º{2?1Ml ­åu02í¬%a#Lé48ÈðÿY4GcRªæHvf²;U>0<tt0äÈ¼>GBÐIÒVE=ý1Ï%¦Ì+ ¬	yM]]',øf²&,Vã±ÄR¸'­½æLaìG£®Ì© ¹>Ýf.¼à#Ê[=rØà|`0EUc1?Mroè?1Q©Pbó
ïBß	ùì£núGYO'ýHçë6¬ëÆ ôîPVfOÂ¦t,þÙh&d°jº2é½×u,ÌHæ7U ñ¢¬8VúoVÓ)¢/+½^odVBo°T¿Î×>²4RKxÈDYÌ^TÔsHÇÅÈÐztÒ0ó|4çÆþY½t0Ó¥_®CåêqäÆ³®në<öHcç¾çTRÌqþN
ÑwnôV³7¾su¹Ò÷!Eö¡ìþYÂ»ÌSÎØUýËzñÍ¹¹1PÜrÛåEwQtÌ*åS-§ò6tÐPà
é
Èô
0Ðç¤	¹ìE×^	Åg0#!tádKÕ9ÜLhw!o?iÂÂT%Jäaaö^Ï	KKôs9VØ5kBGÃ!MÃù¹Ú¥â
"9>²n¼ãNló³0®ÚéÆÐJèyÕþ$Ìâ¹Àtª\,LÝyÚË¦ÁÛiÀÍÀéÛ ÎÓìÛ ?dÀ¨ÝêÌã_c
&i½]¾HaWFUÝô_¤Á1Õì"´_gÓ#ØÈ  ºô#³&'³U ¿>"
ÍñÿÔd¶zj7"zg}ÝÈ©5"OØü2"XKreo'ÀØ±¬²^§¿1ÒÓ¯ç¸"htN<­ÓZa-Ocl§ÎÆäÃ ±tôe.Ioìà%L­¼»bn³Så +3æiÖ¡«æjqÉèæ¾àmfiÂÃZç¾_xïÅÌòhk­ÛÌ¬c9ëÌ6öif½÷¿¶ê ×§çøKu Ñ¦ôn(òâ¢±ÁÆëL_Ü³DKu%Ë<È§²'«¼Á£ÝÈºèxð?+
µ,»¼³4lØ:4·¦kMÞðaGÆTÄ ]ÞÓÏEw,H¬µQl472tÂ®~	<*[àAH%ÅóvGÔÓUýÙ®ÇrÖéJÙ<yük 3¤|;zû+´(íg(<r¼Í3cüiV`íV	%²ÒÉ)8%ÊE¿Â.|®ýÜÃÎÂÎLÑ°{ê6:åÎOªb·¼¦Ñ!ðV5«±ùDB&Ãù».É<Ý´6ÞcVBµÂd²ëÊ`ì(º
pÌqlÌ<	¼¨	é;$7#û^Ã.Ä#Ç¢	ùB=&,-"ÐSÒ¸£À2w%\[a%É\7ç=Ô;ÐÁ	¾;ú¶Ç÷MÚ¦ÂÜ+>6ûësÁ¯ÆÛÉ_}ý­Æç(T]YË*E áY!ÿ¦ #M^S Fd¦)øö³æ¨)óª)"ï¬ý¦&åS!3ùY¨¶
ÍÅ3ÏðÜa´A·cÐs(¡äÇ Ë3}±ïxÇ;tËØ>zÌ-õrqÝà^Ç£¾}¿ñi^6µïßdWI;§¢aÃÁ1WÃÏè¢þÎ«bvùò7'¦zý=òDX[S¨5+>gcn)ìXØTx½RÐO §)[ÎkVÀta|\ 
6õè)DÕ0ÆôX0b'ãm2Õ1¨ïé¤©ñ² z¾\ EäïËn«6øÎ Û\|ÔI ØFÈ"ã;`9^ÈÞöàÆÑm0và§ÚÑÞÝð c¸öv0Cäqzÿþvø~ÃU²/î[éK¹U>Ø±+kdÌÉ(Áªzo?!ôê$a PÏ|eÐq8õ»_jGhéÈõ§jèÙê&ß¨yÿ­ªÅíM®íÉ&È³ÕC§ñ\ª+ÇÌ4Jyîs	úÚ	¡³tÕ°®LôL«ä4×P~ltGB¬0&·ÐÞB)a	Ô,ÕNÁgUçì:Uþ»dø)é	ë<ñÀÉË}=¡åå··9Þ¹B`ù
a4÷GÂXé`uU°6ã)|"¶Ò~agô^Ü<VXàªg§«´áKQâzìÔîgîSî.¾!qÕ»õéuwÑG
ûÞ±»Æ¦°®¡CNÛD½ü6_ù2ÕëT¤Á|£÷½î¡ªïL#³©}oÂ?ò-OÅòR-øN&i¬­bQ9¹)^Øe4¬º°sCÈG+&òBVs¤äS-*Ý¶®eD÷Üv¥¢kéCì¦§ÛKCÌGÈs£ûÕÊp
±£ÄwéTwK'mñ
[9H ã¡nÆ^¬11MxR×aÐ=]&ß¤*ýã0<2±èy`!%b*7¡3vf=NEC×1Mèx ãÖÈá¡q;Â[2¯m%Â´wzäyµaóDX	rå¸V®Úã8ÒFõÌº¸ÎZEvôq5¦5,Ý:¤'j7+È_¾£zúH=¿©rÍÒÒ­Õa+Ý*§âäø¾ÛÌ7ÿ9ÊYñ~éðR¯8ìö½I=úYì¾ó9tËäUn¡>Éò	ù]7Ü
ß½ÙßÁ<ÙPÈÑ5ÊZ§Iqäôd¥ðh¶ì¾+}måJf¥ ÒrÕI®+ã VÖ¦©ðÅÊ¦ÊiuÒgÌMÁÜcqÃdÏ_ Gª0Ñ
äËPÆÐ|À7]Byñ(@§Öbôdcì*÷:kt{.k
.[7¬"\4ÛU¢}t_1¦»î
 GdÊø¸äX¯÷LA1nyÐóðhÌÀU{iÓ&¾Ú´é;gÊUËô74t=é¤Ì>4²Èx5²û»Ì;Æa.á';ßÈ$dBýW:ËhCÜÅõ`ÎÊXlÐAÐ@KôAçtMÜjr9åá¡î1¬d¸íJ%&¤3È¢Û(6NÖðAì¦°ï5[òÎ·åôÚÜUeÊÌY;µì[aQ:=Ç
s'WÝ2ÚsÿJNX!óCX»1_ÙáC©'Bg (Bh0Ló­-ÎwBËÎiB¦¹£T;¢[°pm<ÛÐ_<Cú Z}· SÐèT L§ú\RÅß·Ù3Þ"E±éD[÷<SÔÑ0kB¨Ðëu;Á3=·3(Bo¸W°È{èî`ú\òý«æè$ÂºdFU´^ÂhB,¦y"¯UáÐs7éXb`
N·i|ä³êµQböp|ßmÌÁåòôï{èuõ;}^ ¤¡_[
=2$FéßïJAT3OÂò~É(|Ì'
põ[ÛÄç­}Çd§T_øc¿.¾Í'CU0iÁ/Ð¡Oµ>ke ÏÔ$´ðV|Ë±åÝ¸<ó¯b¥k
hö7.ÂóAµ"ð;%ÚJ'0dMÁ=Hfé
ÀgøÍ*AÔóºç9;@cÅ¸Æ[§Xeïß·j
LRPþ1h
ê=ôCí÷,äílàslß!¹½¼¿*1ç¦¾_¾ïH§1A=|U©ôÈÉézV:o¦CC>µX{O~6 ´UõÐÑM»Øiåv
>!ç±¢`hXI_Ý ·è¢	îpræÿóq	ÌÎ_lçr-Ée1Ô-}U"IÌÒzO¯òÀ­Oü_'üKí¼å×	üv oèuØ®·ß¾P¥=añªkidÎ#	%VÀ\ÄÏãc~ÞÌÍû£z¬G|_ôï±ÑØ¢µØ½Æ'¿+=Ç|Å'_@½GU|òb¶-Â¥{
Ëh#wÐGëÓiåÓ§N.aJá3R5Ñó]|NQ°!³¡¨¥2_¢H<¤rE²ÉXNfð<«zWÀëðÙhì/LL»¶ôÉDØ >'ÂÐWDæªÖ}
ë¤òLD[t¼ÚBôë:Ô=&ñømeÅ§÷°tø°îÎwÿÑ©)Î3Ø¿Lf±ÖÓÑ^»x­]>SÒ-ãg
Ô|W
²et;EÁèÏ.éuâøÀß¸¡`ØÍ¼)#Ô\3õ6G¨A÷lðT¥È(9a·Ñ"CcÄ:âÔ{ôýextìÍÞ!§54ýâ ôèwÁtlð¢Ñ-ÍXFio-£
úm(ïÁ0v³yñ1¤Ñ§8DáãLóF¹DdbtFÑ¤ÌoíUØòTÇílÉ°õG!÷÷é#çëw:ôéû·V­îÖªj5Âê±Uµ(¦zs"cÌ*ÖÃBÇbCGö¡ePÚñ+!3`		ñ=4aàÁ8ï!/dÚ¯P3Ü×sS[	ÔGÂB@JÅ. ßðÎ±	ü·Ð
I`yÑ-b · 0ÀX4zHìÓÔÅµ¾ZÞª9%#h§È±z8¨5¸ã$pX¸q·ød;¥°N;×[½Â;Çdk2»0ÔÕKìÑÌk­_m(b
b¿®µ¦¥íÝcÑ-¹¦ö-!a¶ T«ÌUHRFÝÐ{kµº	ù.É1þÚM­Ç.è«;.uB'ä1sG#¼}Ö«ÈÓK7é'ÂÈ§­ÆjÂn¢5B³Û(
W÷«`DÇîËKöÞ
|Ç¦ïwÙ1IìæÁÔ}ãxÉïÀo]a°á ÙA	
Å¦Á2 äAFëZúdÇ, ¹²~·¾u4Å³û:wË³Þ÷ptTíï ßZ¶r:$^%4e[Ö»÷Y×åæÙHqÆuÖÝÉè*ô0Ð½ÖÛIß&ÉÐ7üH=
ÑÌhL92)ròE2db÷"¬0FÌÝ´Ç½tSÔ°,C,«DÞTïlsOÂÇcy»«ç#=XJzäù&ÑÉAÓ¾ñü[øAb53PÉ©
ôV`NMÈ«þk³ß2 ÜÞ­LðF#!£Ç¤OÑã Ú	¡C5:oë:dl1É®¬-k3!m8Ê3­0õ#Âh½eÒ<¤cË*ôý[Ç\Èæ{%IcÜõ	ÇÀz"â
SëIß]ÈCf¨½¦ZjN'ïân²çýf&£"dpy¥¾ â¬¥Ë,Éq/ÊÇiÁÀàOCQ]§w_Ê³zzÚxN-¨Ñ¸±g.+`»Â vÔVÊÑO6SÌFoèÝÔb:}ÿy§ß ­£öý}ÿýË@O(Eobý^éó=}ÿÖêÃ÷ ©êúÀ`x¯~åöti¥À<ÚýûÛá{ÛÎùwy>FÑWüi¹2ö×Êa%QS]!ô+ÔLR*÷ä+ä{EÕIúÁa +Íªk·¦÷£{mbúZE_ èbr$ìº%J½óRófÓÄ¤
»êf^ÍÍ-ÿ¸WI_B}¾rx ]m0Ún%!ë#óÏ*~¦ã¥1¾4/W%6w%(+0Uÿ´×\Þ1ú2Ã¯=ÝøPªè$«OÍ}ºE¨GðÛ¬éZ£;ÊoUéß
y/×(6 ÖQ¦iß ûc³(
¤¼¦àÂHO
#1¼øy¬¬J¨°gØ¶ªd]ÝÐ}ê»:&ïÖR)FRT
 ×Y¶ÔÂW
lNV
9TMÝ á
 Ö¿xV¬eðîyiä#ùl(2V¯ßÔ!â<Ø}aÓÍÖ±çÃÐß@¿¡KÆaè¨Q;µÓé&7øÔ\Ú
}9
W£{£ôïÒøþóþÎH¼ªYð3¦Î[ñ³gÞ¢ñVqÁ/ïï§ú¼¾ß|Då]¾W\U]zé3 Âg8ä¤ëiùÙÇ¬*¤ÂæëRa)ÚÏÑÐe3Ã)Fo ÔhVÌ3%F«¢(ÆËñÐ
àK°i8ì»p½ËZËà)
£
)ÀP[Ó´êÌìÃ°èùY6 *´öV@AWt*¥ÙP@-â¸RÌ³6Mý=ÉSãÂÉ-jÕÊIpÉä5î"ÔÝWÐLÊC÷U= ¸ÏÛ¨)¾ÿÉCKK!BÿÚÉG§s];ÍÐ01üR9'ÅÚÈ$¿âÃ¿0@J];92».cT¯ÃåvÎ
¬_æù¸_WcZD
f·ChózvÅ·F>ÉÑF©ùø)ÏwøEÒ³¥{{	ÚúÐsl:ß`MW2sÔ³¥üöYþ}ÒÏQ@6Ê
@g"¶
Ð);õJMdlMceæ?
àI¾}6X¢Eª ÀpÉÑ
ÌøÈûFn*
Aà*ZºfyÒcl\>+ã ªºòE ¤xR¶êrtáßEQ¡®EÝQT¼hGK PyÊô¥þpÏ@fº?ÚD*lcÜEÍ¥¯³¢óuô»à¡è,MÎ¨QStÜx`aÉð&õÕ­¯³}QtÂÄúoCAáwf-ñ´áÓ¢Î¬¾[+þ6ñÀ_¾b$¨(¨s(
zu§(|º2*@c,gE20vE=lc4T?î+"A·`îyEè^Q®(ø2a)!ép¥`þÃ3KQÀHÉ|o¥),ïÓ(W¥Åã74Rx¡(`0r×¢Àà#ðv4Vt:¤û¶¢Ø?¿PA¡+¾X°¡l µÐeM2ìM9GÅ¬( RHV¨6OÑ/»àã^ÒmwÔ Éé;èDà)ùhSg>-©	Â~óÇßå9ÇBg'ó
KäÐï=?¯Ïüç_ýøLÀø}÷Çú¡cðòÓõý®ã¯-ð¢KÜ(ª7¤kï5/#v¥4w"/JØÑ¯æüu>{S
î©PúÚ¥¯×ñB]Ç3³pyÛãýóºC
Û¯NÒ~êß}öï¿/¼èE$w¼6£ ß<LñýÆö(o©×ï-ÉÊ§PÝ®EáÓMé~Å&äÂg6{Fu¿ ÕÒâ:É¾×jûÆ\
ã¡¸b]Ãýqî?Û ÞþÿüT¹kfyÒ_SÃ¢¸Ù64Æc6íÀÛé/2£¦ 8ü<q8NkÄËÑy? Ôjþ|ÿWÑ1Ç*}v|%`ÝõAÔAÜ·H?)èJN-BÕìé³»®Z¾,fJEá
1ê0ôåcOQ8Ï
GuC±N.Õ9Ì9Íáx9Õ¾íè.Ð@úî7Å¿Ç~Øm§åÓÒ2þ¼ñG«pó½Px$Áßóõk±Ý÷þªæÜP½TjYÞÐÇ
ÝISúâæûÔÕÖÍw@_+7Ýë5zò·ëÓnf-¡sWÏaYÁgÀv¶zéï±WbüõýlÐÍloÓLùd¶3§B*z¶?¿ð×Ù~CÕxR­sªA_Î)Øg)Ç¯«7úAÕìÐOTùµÐ/Õl7Û÷õBpÇCÊKOö'÷#Ý°íNú{ì×ÉÞ}?Ùt3ÙÛ,Ó_XÐdScÄ.Uµ½ªùMþä¯ÿN<Áû¯ÿçÿë¿±^ÿ¯Dè.~g^z¾#çvr²(¦:jX£3Áý÷Xåó{û<
æýá®ßW\"ÿô<×¤èÈF×Åï¿~¿7Áùd#úA©Dã§(±ÖhÉïoûJà<ñÔpý5¥:yòõCÇ¸+jF
)Ð'?wM?tùÇ=Æxðóù¹:rÙ w©¸)ýâºóuäWuÛ3Ú³h¬XcÆ@ûâ$¼½yÃ½ù5¬[rÐC¤^î©Æ'¹5çÚàç~Û"tÎÇ¶sÿÑNý¬±¿öðG¯kà¶Û?Cû¹Eû'Ýý¾ý¼;ax9*ÿKc²¶xVÿwðçJPÀýRøíçíôùÅP|ÛÐàðùVÒöÃ÷Y?c¼ø-ìåì±Á<U³¤~4Í÷gr&·§æ¿¶_÷¸Ë¾vÁÖyõ ø±­ÚýÅUÆlÿ¢ nØ[	</74UâWj5uÝëº´¡ßrXçÇ¤V×õêô}7O}ówxly;;LOå{6­?¿ÊÆÒ®Sùñu»Lå¶úâÇöjÒ¬
ÜËüÐ-ÿ:À^Îãµ¸/Ôh3¿TÉ¸é3bá­æðxNO9¶~ÞMàCü>ãí¸ngoéÙ{O õçÛ>åÚkOûø	ìãÕêÓN3¥0ÆW;.pùZØøÍd}H²A½ÜbÏ
Ò>ßØ-)Hûü­ág[è º²É¨Ã.lrfxÄ´b×â$sÏòYy%û÷ë;]ÎßÉ÷f2S!´¤Ð;ÿîüUÁæ_Yé}êòíá·´Z5ìÇñïÅ¯®ó²Õûêë8µóIjmH3åïãÅÐ\<£í@/$ùµ>2ð:¥²×+ø¨üÌ8èÅ°ôQÖkäí\ËÙ|AãÂ'µô´ÁàNÈ×!×ýÏ/Ëöçæ!öy´£Ýp?¯a 82re4À=ãÓÞð=+.c¸éÂæ\yÖî¦Ì
ø¿Sª$N(n÷.8¹÷mñ§²;yÏlþèq»ûcE¢Gæë´¡÷Fña :ãj@q«Ä§±qT	ÃrïQ
$y5.´z~x|7½ÜìRS; 0¢}°Ð1z`ûÒÂûÜ	é3Êø+GÏØÆô²]9ò2ZOO{4GLdÜ,ÅÑf[³õÇÒ	],Ìs¯	±ø¶óXÀµhB½N>	xnCrº5Ç¾,Ç ÙGoÌ\d FÄªWµ2ÊwìÁ} r°Éh£ÌXa/·èÅé[V,£(ÉFØ/ÑONWZc¿à
*Ã4 J}à¾	Z+Ç}M²;=+ç!ºÉ 9:¾Û¬qð97­½V
 Ë
ªÊyfc4éÔsB 4GìX¯ÒÇ>ÇíÊÎÞQM¡ëaBÊq
}ðé Ç
¦Ðq:Mó,¢&elT}oÄæ$ ló>
ÜÄÛQI7a¿£5!Öq
ùHÈBJ¢«µºïDX	&Eè%4o&Ó¤4ÎüÞÓTq(Âªú +BWø$.(Â·|LA'AÂÔXN×cé+Ür=±ì ß ©êúÂ +S;Gª´KBùb 1U«ø¤Vº»ê:Ã<xðéª¯Ä=âÔû}
­+à°¸0]`ä|"a2&ümèc{ñlÏô}iìPðÓùù)jþTeÇÌ¼PiU¨MÌB±rÄö

,
 zy.¦ øîq­¤iàeãÒªAx³C`ºcpì«i}°.­oöMcóq"»¶i¯½Ï±ÐØM+ï7Y
ü¢ëdV*\%µõ
`äÊÜ,
 )° q¬ïíÌèáö7w^Cd §ü0-q_¿¾¤WûE;üÚ£,ouÄøÎya4Þ|ýÙ	_ßØî`aÅO¢Ä$ãóE4¼±{3îN#jNÆ°AÂ¨bãÕ)Ëv_¹µÅ3ÇÁÊL"xIÌ¾ãÆø[¼åWdÌþ:}_dÉ8§já\ õ=Àt-&«ïÉT,¢ £k­å
b!Ùã¦ÂÏA·²dÌ®pk[a«ÉØéûÒ9Îx
rÕ(ì÷L}j³|(ãZyg%4j¬c¢(ê»c*¦²Ì;a[üVíH'õã~´¦¤"và®À gvkM(Át/&Ü4õÕVaÞÕ~Ô¶ù*qk¼2÷/H)öËF×ó½ï+ ÓáW"hÑTD^]$Ê³ÚL}Èú;Ãa¹ª¿×£D
~nùÌt*Þó¿/ma"¥'Y­ò s©wïä¬úd¯ý=`¼ãá@/y¹jË00¢Sd9¹B=Õºæ\Ò¥UJî}äLa°ÕXÆMðÉÎ±S%å][&+Ùò}·| djÌesûX>¼±tìþ9D6sý#M;×?·]
2¦³HÌÃ(/"'^ø1iåñßãsa||yÂøøù}¨kÉÌrJJ.ýôDLH÷+¼ß}ôÑ£ÂMAÅn­cªÝ¯Zì´MÆTy°3y,°²^QbÐß1¨¨ATÅHè/&TØ×:2 hô]ôP{­fGgvÐ51mÔøXØ Å5@<^$>ë÷Lë.[ýS	ZúNÛáC×êÌÛ©ó×Nf¼ëañSèÖ áéNU¯QôÙå;Cw(£kµÓ°2ºÝ¾ù&*Ì½².àF¶[À4p%DÒö>¥ï óýÉ=íÐË`íÈý[ è:Q_©T¬x11CÔíy1E©¬Ì'tñ9.ýêW¢ÂpFn))¶G$´(nÞ
¼/%Þ==u8Ýÿà)[â©I÷FÏÏ°{gF_KKQÎqWÈ}/5©=t>XÆï+eh^¢ì{Ñ °×¨¾µLÃ³¼YôbÞ-U§ ªâ½|§V}Pß¡VVYêëwÆcNIU¶J´Æ®Ã×,EyoÒ©5ÉåÂø^+C©¨É2È²5K¤*ÏÛ(fsSß+xÔáÉ8Çú{BdÅ!É³3ªÌKÇ@
ê:
$&[_Ùu2³èÆîä-#Æ%ÓÒ+uÞ-ZîXZÐÍ}öjµ¦nµÞWð9&¡½ÏðPrgÛ>cù00úQaÜÄå3OyF¹~æAÛ}W¨_Ò¡BÕa®¤¼A·¨è¿ÃOÆ¥^¨H61ZÈjÓàÚÉx£µãDM<w[Ú+yÒzìB'C4Ì\Ï,vê;¬D§¿c@Vß±ùFUx9
cVÁ`M±r]ÉÞöàÆ
ªÔ4¥Þ×ï©t&a£ûé%ìàQP#H·KØfè
Ó|´ïáGe1àíp=IúZÉ ó ÊcÕ6ªÍªÍÌõ^µ=M'±}oìV0sW íóûòêOÄÍ.5aÊ8Ã7Úè0´rc;$@åÑËTe:·ÄÝQbjæDP3v¤Nô¶ÿÜ°Ñ¡µB#2õÏrÎÍ¯·W^ò<êE$Iæµ~Gu$pßR¶$ñej_dÀµ9è%O²ÔÉäøÃ« 2u(
ºo¥¤\+å@7i¡Ó¨bx¾Ê+Cí2Ðr%ÃÎ¨0jkü`¦­¤»°í<u!Óí1§ÀJÆC°Ê;´}%1±
Æ¯mPèTm×ïm­®Å,ßwk¸PbôR©ÏÝuìÆn@¿2<êL7Ú¾!ºÁ®µ-µU¼Ùùö¼
7¨}f¤bW^-ÆÚßóªôîêÊ÷ÔX*ueiiSÎñÀOV1qfÕòÞð{­ÇeL¼õ{Àjej0õ¯±Õw¬5îlkµ¢$½u§jEx§Ø¥¥äÿ÷n+æHºæyßEx@h¿É;kHú8¡uïcþO_&yDV5«3ª"Ü?³W{LÉ,Zý½TQÊ&ÉÝ¬µÏ"oiÈÅ{Ñ°RÁz+]av%ÝÕû`«*5|¸sÞ³a5lÎÙàt7°1ØaûÙn°¶3°jÀ³èý{?~ËfYz{6§i[G_¦É´M^5/î9±J¿÷Ùë¯úÎMduîo¢úLÍèÂ÷âK,Ç,K3ôºMGÛ²E/óz­ÜË dîìØï«uÈµñýáï6±ãû·àþàI;¢âÁ±Ë÷¸è½ûwË
NM¨®ÙK+B¹íta	7Äe®FàÉ«Aøñ26¶nfþ{@ø:;ic°Oo»4);^â­C4û¶ynN|4¥xÄ¡ÿñ|'av3IGËÐa²MøÊÜ-YÀåì\aaÆ99Y¯QovXuÑ+øæ;L>vï²ã,';s~B&¢`> ÛA½qÎÄuë^qÛ®whGr·£%l×)âAÒ,)rG@ë[³vöe®´»M¢)Ã©Ë­eÓjÐ£Ä²³4î:LdÈÊ¡r¨«õDHW7®nÂeÆ¯ÍÁñ:÷.:Å¾D¼§h×	ú¤hSÓæÈúür'ìe³LöñèSoßã|áÿþL¸×ÅÚÆ/ß9X²ú{&ìz÷trB9(óûÇåûNá¼xÅ_PóÌ^/Ã-¹ÿÉ¥oâ¿±ÓåÅ'¬ÍÜÓiãþÀ/ªM=©ã¤Ë49í¢n!±ñÏáè-MT´0ºIeêÆØÌÙvè#àðmJÆ#òåG,ø©Tõ"Tø¥DtÅwÂ1TZcÍlm­À qàÍd¸sh¦¯z@ÄK3ùÚÙùâVÐ½²¹ã¥7ú|2ÒËÈ±·8*Ùì­Dì"X%7{`pbÿR¬õãLU(ïùá¦åÛ÷í'¥ùëENÊfï¸ïÏw°µß;ÉCP¡©fÅ_#\¼þ.3¾MÔ·thÙBÐ·WÚðÖ2Â1IºÂæom£§42I3+ÎÕó2!U~¡xnP4¢bS4+¥MYÈ(?Dç$ú¥Ê]N#pÉ³îNÙh9ükDVjV;AZ¼îQü û´æ1¦¬²ô·N[î¿Áþd,ð9!õ·¹D6ÔµDÝ-w;@Ë-L¸¾¶ù4û®wsú);?6AËK$Éw/j2FÚ'¿<§È#¹ì,·gãûÇøî9ySÙÔ\
éÍL,lË©4×ïXÁ!+8Ðf8ÃªåjÙ¨JU?]Y¢7ÈÞËk%{: ðÛ£±÷aàÆÙ­hÖU¨S¥Å¿µø ©;d0÷¬øl? EmêqYõX©ov=:Påh³îfÌ¢;s² ·E{¯|Åµµ=Ô¡zu>.ô¢H`{è nÐZ0Ý=+Gï¸6w2ÿÜ 5ý
5aª÷Mu=»HÅ¸Þ¿{ eY'ö^¸ÀÇ	#ù¨ 2ï¼Úcµö¹çyp±M­i%u³)á7«ó´ÁÛ¨&<PQfê}¶ç"m¶äúqùÞù+ïwªêã^9O\¬:¦³o^»L¦l@ã"./³m{VÅu=yÒ×ÂQ6Åû´eUæÜuÜÚéóB0×t `8í¡²×# &`wsÕK28¸ùÚÛáãBÓX¶ä±f¢ìÁk¶fIQtj&°EÇMXÈÌ¥+°@µN&hyN&ª·N©¸ùR@B¡É ºvkq)5"zÍ±ïÍ9·ÿo)®ò¡§xÄÎõ_#ðË@{Øû÷(*1Æ÷ïÏ¤n¡RÒÆÏÍ+aúp¦§ÒïWõûwü¯ù³è´ÁUÌ»oÑ¿TÞdëò7-éP~ÌsR[ô6ÌúÂþý+~BÚ7´¥|G8Ý¤¼oE·n«4Âc¡áy©7Q}9!Zâ¨	´ªÛò¯nt®dyò©	ã5s5ûp
HV/w'>ÒcýT«l®]M.}×(£Å¶åUÊÂ ÄÎS*&¾,äìÔ«ñä0¾ãMÝRðnZ dö$½mñ¼½x¾À
o¼Ü}nß³RfOlÂóV,çYùÊYu{öâ{AñÔ6bç¥ÑjSª·ò ÏbrÁY"[ÚçãF)J	ÛùÓÖ=Å&YMÑ@dMJnãN)b Ü^¿·ê$ü¤Y-¯hr{;Úõc8å©:W½·û{=xÓü[=*Ñn;Þvr·E½&Ðµåbr9ºCoÿÏ6
|¯CQôä¦%(Bö-D©"âû ]¸ÂÕPñª[PÆ¾N+BÅäk­ª¬lÞw êA|­gG åspU¹Êï±Ê¿eµÞÂõ¬ÌAöáCTíä_/vVy¾töñ¢ö~üý;¡·Ê²¯R&ÜôïO³´êx=èÅqÜ¹!AK÷ää{"ÜQÉeQMÒã~¿FÆèãõ\ùV"/`ã­Ö4ãUãÁ7 @f1{3AÉh×fvgt­H	C|^Eððº-XÄ^[àK­µ´wpÖ¡Ãé³Óg5½üCtÂ§"Q¼ÍÝTñ%Q/Usßþ5Ñ j|9Ô5¥WÈð9¬e
]æ3w´¯wöÇÉ~d¯UÃ¼êt]'á©ëM¾ËãÚ¡®gÁ«Lµº¸;Çu(Ëï?2_ÚVéLÚSüHÊÄUÀa#¢.GÕ?=«ø¶à½ÅW¯"¸@á®ïaïyDÑÜÚ;©X«uÌ×ÐÔo°$bÈÛºÁòQ~å<fa¶éçÜPWBS×6ñ'å
 ONg(egKnµ×´ùÔCD¶ØyO´ñ½(úÓgîKV7Â2*ÖfÜ
7NgÖÇiÍ"¶ðÕìÍÖÙmxð£RK)Z"kj{gø´2öÐd,Ë ¬pnµÄkÅô|ÎKþsQæú¿e®ÿo)®ëOñ4ù[Êûå%¿¬ó²wâqJgþ@*ôï}bÊ ¡#yb­Ú¹OJ95wY, ô,R%l-ÇrÏAW·Rß[ªàUr¢GLß+Ëu¡äªp»i?][«
Â©½kû=EQ-éÉþ!`p SÄÓHÉº©Ú[&¢fÞÊO¹XEñtFM½Þ-ÛËÞ!sí~(j¬Ëªk)ÙdNuä2|QG$JÞeSøeø"÷K»­×:
M:ÊçzÎxDOPEaTm:ôì[Vp·<úÕït:·ÎÛãOsgïËâàDÌÎïßû¼õ<uPü©XGlå_x#ev~é7C3ïüD£;¿¨9í	ÿþ¥Ó³{Ü²K"GÂò ®èªYÞø¨ë&Vv¨6h§2Ùv$ÒyàÃ4àÁ3^#d-!5BD+é5"Æ±!U<IlàùýC!
w=u©à¾taæq·â+ÉÃ3Bã4ÕtJ¾y ×v¸'öî·ÕC¾×´ýûÃ?|Óöïß:BTÙÐæÅïíÙsìßk{gTÍÚBDÞÎo[hqçÀâ	 Á
!ÛhhmöÖ
!ch¢ØZ!
XìBøñ"àN¶T³TøûàÅug~sw§ßgÀoÑG2F_üN³o§©ÂÒcvWaä{ÍZths$U:,ñôÕUËV¶Ì´3õ
{¿nhQVXî°Z¸³äþ3ücÑ==zLsöt.7	ë+Di/É:®>Ê¡(IÈ:Åæ,ùI1´X^¥A¢Q]0«ýà½ÑE5Òp¨\8EºW÷[%êÁ c{£ëh¯Ãm±òºÊcd|b±vR@æ)Ï².G6:Üt( -Èá·yerÒ@Ò-À"JÀ£/!ãx*eîCw¬¾,û·pÊÂol:4Nn0×Võ\túOdMVÝ*TÂÈ¹L]é£bÁ®b¡¸pÈ¼¸~?AÞ®èsÍe*hÂ1_jâãÅþá"Ä0z@!êß»pãt$èïhÑWý]ö®,ü?bÎ:<åD§Áaû«^ÇyEì¸j3k Á¨}õq<Þè,ý;Õlûã5»~ÿþúÛIQüRÙ¹¸ù½Õ$©ïG»_ø8´Kºß¥p¿ê¶p¯«´²Ø,ªwþ*ã4?ë"¦Â§Fÿô²»NöRpäxz0EJ×\d«Okä5õ+nýoI(ÁßÙ=«4Ô;ûöió+ÅºP?ÎTÂ§7ßîÛOJÓg¹Qqü/Wüóï±·ªºêïíà=Ëpÿ©kö¡ä_äûä¦9üïñoyà{áñr7ùñ/¶A¤à&%^u8p£0ß¸oioÜ¬¿ó4Ä³¢!ºgMjîJK4Üa2¿Ùöû¦7oïF0°½Ï¼¡¢ÔµÃÃ<r7Q2rúûÈýùÞFr¹çÂ:M>î¤+ªw£	øOQ¶bï#õãÂþótESþøÏ'æ¹42G(
=^Hr:÷ð{»/oÿÞ^Q!qÂ93é9<|,MýjìÔÂXÖ½¤/0T²Çkâ{­&P^_q.¦ÅÌ5¡X	]ÍÏIósÆîïÀx5§*ÿB[h§$ð.êÚ¦â­#ÜËg¾Þ»	îØ»í]#JS6|£Zì£ ³·
;¡íºdqÏF×{''¥X|¸uõâTX¹f«¨ëWª
ð®¾p´#¤Í³ìcö	óÚÉDJÀCìu#9ò½§ª& írÕ}*»6?ÊâC>T0DÑÚîSÎ3JRy´0$Á7õå´'Ò÷éÑC#2iÛ]ù¸°ÿÝtþ¼±Ïð3qÞÔ]ü²¡0;µ¬ÉÎdolåÚÁ&~µâràÕæÍ¬
?4fp/íLàpdèõH´ñr*Ð;xMì»·@ß'ö{:Üºa(ùÎC³³"·t{ô`'áãFèC;Ï.È¬íyÜ	R­Â·;éóBà~ËäKEæû!ï¯	ñP¿Ð[D$(ouY!G¬YEDlæUøf|3
i3©PxÓP.çPE3âuã¡U$=$k¨7zch¯²Ú/ö{:shvÍþó¡ýdø£8¼lÝò¤¬OS¨¿Èiâ>nÞ©²ñeÉîP9À"ù;
ó­wÂ2ÌOÝ1Gó^4ÎÞÒa·`&°Ñã¡²;aUÔïoClXÚàÚÊ0Gê¹sTªö9þë0¶ef~iÝ¨Ûç°®ßÆRV1¬·»¡5z_á!ÔßTÆá
ÑWÝæ[_5éûÓIx!¯éô_óàÇÑcÝºyº¿Jýì=ìnk·J¶¸Ù$o×RõW©iÖËðÆyZcª6é½¼JMY³oÂaª£ ËkÐe	ËÙ}Ò\sê°$<'ðLtÌJð¬ËPüÏ)5-¹qÐìâÀÅq~®,MpÈ_1ßsôM3.Hl¡vq~&Ï¡ÅÕñ@à»Ñ`Cójsí¶õÁ¢J1¡^º¬S3-4ÑâÇ>RÌV:(.Y}hx­2þ3MNtn¢'±orsþ]Ó©
ýií¤ù²¡ÑKqM2íÀ/N®@^Guu&á©Î)ScÌ8ôtkuFæíGºuî&3}#óNÆ¬cn ö£zø;ÓlÃ¿N sDZ×Á!­k¾E)ÊNX6Ñî@ÈxÕ«:+L'r¨·¾ L_Ä}ÇjW`ÅñtÒí¨-÷FgÅ5I[òfsZU0ãã­êæ²}øjBe¡P´ÆA°ì¹¾ö%QVe±Õ)b3W×)V[]S /±^·§è¹ó­àhV»g%Ë|»¾eåm(H,´j®Å.¾Ù¸Æ«p·:;p'kÉu-.Zëýê;ñ['3ØÇòÉËþ0¿ó¥<Ö­B<dC÷Ø`<WK¨Ro°1
	kêY@:ìãò½óÂUMó,Ùd÷òBn¡lT°b°
óÖl÷GnXºàBv/E1Òl½¶äjKnÉáÁÁàòmÇ5ûvQý¬sæINòÈ;©\í5G¼][Fýøî1ñ
Ä¿
þ3v`ËíF¥'vo>ÐC«äÊu{r8õÌÍcFt·AÇCT|Uý¦¶tÂçXÆ²
a®(ç9á®&æ©à·¸¨"W\à>ºÊZKFï@ÂGFU×	¾0îTWôòV§(
±¾ÞZÏIÍeªêò½"óx ª&»&ÌåPv^÷°\l¯Bí}¾Kbxr(:vý·æxéâÕEÉ¶CÜ/ÚE`n3Û,d×)è
zbÒ²[ì´Ñló¡Ë½¬!z]â^I\[£nî­40`]îEÓÜ	Ç'ÐÁ¯4úJFV,ü¡j?ì}ÒÎï!ÏypjÞWððI·¦¿ø¦Ú+©5kAD9
e.FØl83kOEóË\g»IîÓ@}ïüÞ½º{ËFjQ	ì5²±¥æ)°kØÈpM L4gûÒÎ¨bàÁ¤8ì]¶-ïñõâ²a<û×ídE¨¶xR·¢xB
G£¾r/®90¨6
¼¸õ{å]&Û}40ñV¯ÔkQ3NWBÔ­Fôë)VÍí^ùeÄírM
Èõd×q)
V*¤µ ¬rü¸ÐçRaDXÝC#xÍâtÔè-BQC]æÎJKô×Êµûµ¤ÇÈÿÊ@»6gíGñ¹×[KØÃIv¶jô.¬[[eù©~Êoï9(Ó@Ç{­#º<ÕÑ@©·[=Ç6rß;á,ßªä´|/ÐSy£3=ìãÄ®EìÍó{FT¥¼Å]È»íâKÜs*sìP¸ÕÙ-*GÂeñ5Æ©`ò¬9ùÑ{.)ÆaùöÄ¤µ÷:0­i¸;PóK÷ÁµLÖuÀ !× S$ìPR­2¥ÉÜ(ë±ÁRTX`¯õb¤XJvä³¸4Ö@b.º¼EÀÉ¦¹ÝNfbÍæÚÜÎI{]H|R5!á[ÜQDoYÉþÝ¯â
{å>.ôb)¦ÝtÂ§"`îõ$À&ÞyÙíàÛù¤\yÂå"{ÑT¦È% #5®6|Ã)q1ÉzD± ¼ªéBxjzÑ	\¥ý¡"Ñ«*-zxlïJFÇ¼öÉlãQ¹uÞhÄ+P4´SCñ7l«½§T^ENâÚfÖp³ùîêxÅïìCßq[óUÇg«¯_q!Ï»nðÐjç÷¡·-Ö,8'Ça¬³8|Å§Ç+"C<q¨;oªiläçÿZBLþ«¥µªÂ£o³zhB¦¨l;6Íù×Ä9ª^&5bðöu×¥"hTAÄ¯Lñk'åP`ý\ÚÈuEzfk¯9Ùçð
(ÃJö÷¢Êõ
¦¾U	dã¡òDL5î^yúÁÖ]N$Ä*"ëØê)~\èËñiîS5Az`Bðq#µ2!oïY+¼,ÐæÀ?y¬û¼§<;}.¢¢QÊvè@ÈD:^æSÉx%Py~ìÉÍo÷W¼ÒWÀäÛ#ks¸/äFôDOØ¢ú`Â¸Oïy¦â¯ t=®ÛyüyË½1ZÒà=é¡9FeE={¨s¡ÁtÂÏB3§ÆqÞ%½/x_hð@e'û¼èOÌd>;Y"~pvXÉ¸}
¢ßX¶H;#ðÒ-·êÛ+QÃC2;ljÕFvùXì¸ÿÑ<k
´²g(1(àÇ¾ÜZ8ÞÏê]³`^ëà×Û=ÙÆÆÄkÈü¸Æy¨«Ë8Ö!C(³C«àÚÖÑkadÏ>°Ì^så$Ws¨¼¨äðaÝtìEÁ»µÝo¦úAÅdTèEQaèÖE¶ø)eupAä|ú¢Ù×ó0¥ûÞê ûÁÈõ­ýkcYß%«ÊWÂGÛt ´ðÖY<ÛK¼ ìeðìJ±º^öVJatìª8·Q¢f&9zjõ\Ý)
wàN8I®oef.bRyÐBL>.:ûË\,OñÙ¿÷cÜ>Ów¿±ÉÎ£qñºÇ=cÿ¾ï5Fùþv9Ðÿ~áèuNP<K 32²¿²ß\Lx_õ¡ÞÆ÷¿»ß¿u!Ýl]?^ßi&pi),½Ñ'Û9fÍÇ«	ò|fì>ém![Í¾YÌ	KÛYúwé?ëïìóòýÇë{"¦; lEeë00'[5LU[xÏ°,º7f
;,çæâþ#:qEÁ¤YEÙË7¬vøFÑ8Yí°1Hp4
|[_
Ì%r:xË1ðÄ9
$Q`«sb¬bA'ÊåÈ¸M¾êXÙueOX¢È±@L¼[Y¶d¼ËÐ@YÌR¹	g¸â1I³¾fæµ§('âíý0SÄNºrG\~×6,ÎTV²p¸¸$§r*cóÿñEñòáíìÒë,ã6Èd(I¦Þê:KÛ"¢(Ù	¨·¢TÚJQ´OT¶à	yÌ±&ÔÈgUºD(-ãúçºè5Í»ê;ÐÿH8qÇ
ªCv 9	PQbDH|u%ìe Íû{ÏCÖõæçw$ã0Î½lIÖEÃÒ¿Á2oýÛð<Â
3Þ+´ÁÒÉöËÄ¢`¢´_rc7!bfp¿ÀÓ lbÛ\RÝ©=GP@Ùp6¯Àv	¬.
ÞHxL¨RXårkm;
ÌµçV¢ú`50áÐ-ï9N`);äX9[×1e¢* §Òâ×QÖ´õQÍ=·|øË×@I-k¯9ÆfëT?q±lVÑ²déÅöÛ]°¤[F»å¹¯ÎÃS{®ÊúïCO})g
MAá¾êVk¢)b~­®Y°ÞSmÀpQËÅpu³ò¨¡¢>ì=ìª:7=¤yB+Y¤5ÖPÖ¿³%fLmß3¬WÙd¼S!DF64rí0©Ió{¾ÃªT%3Î¥ã¹IÄkÇlj/®@'PàrÊ.R»KÓTTÅ-÷ÏÜþ¨q"Ú@.EID²ºÉwy|p|;éØÒ¡ÙwrO­V		~o5ÌJ#öèþÈ}¸0çÈ½ÕpLKö Þ:÷Ö"B®©	 ¬Sä­$*Eù=ãÈqÎ¤$ÂÈNÙöæt²¶§xjí¾	Î·Ñ>ËX¹ÕÊ`sGf¯íØÌB¢±ïUÞPÜ²¸½ØëÃ5Y´\ð [wß S@î=}¡¶wÂQÚLÆæ//vm#æû[x`nßSÓÒÃüÞ¸ê^ëUÏ"M*à7)XÁ;<Òñ
6Ö¦Gc¬FQÑRo¹a\ß(ã°×Ø²ç7ÄÑ¨j[BÏb» ^TêSXL}ìÚïU£ µM¦=E<áV¥È,OwG<Q³ðxÀðHÙì½Ö2ãª©Nez6_¯7 Xq¡êPeÙãÃNpøB$ï§â¶cËá8.=7Å(ÉØ)\ÊÈ½	¹ê¢^ì8]kæ¹lÏ¯µv²¼eQ)F¬³ªÑ)ÊjVW·ÊÓ<Ë*`Í#òµ((î ,ã¸+§ª\Ù+ oá ÕXSÖ,°j $)w¯ãÌQtòh©®ÅÜâN89,yY¯·^>W61çJ¤ÝÙ_É¸å{À{=Fõ°ï±£hìDÒL¹#2¯	¿k;¹khBhp»Á4;ê¡MúÐí°94
lFtKq6{YÞ¬Üò³íXIØâ¾×ï
(z!6g Ñ_|ÂC¿b`¾Ýs,Dvw¢rw/*Ø¤j ¨q©üVT¤U@¬ô×pËÑµÖ³º¨
ÀÚð
$hm;<Û¼½ÇyÆ­¨>5÷­ÈÎ¯Yp\rñEö¨qßëSãòmäËÎÖy\¿+ BÞÇ+°pr²:÷ÅXt¢F<³zr6]ÝuÚxbC71±ÙÔv¤xÎÑãÏû¢Z
+µ_+¢Ø>ìE)6?"^5ß®!\½Ö ÁêùDdS_è[gºC^ï)î£ /DÔñ×G§æ"`²¹i	Liæ÷!7Ýëtª'ßIGeDä*9\'§ã©¢A
ÉùaÂ:KPN©
VØÄQ×&ZV
Ö¨R$ÅuRéakÚnð-Áå(²â Ük`Ó×J½G8aÖ@üQ8î¯.9b}í]Ñu§AD¡mu{V¢]Å¦Iô¬¸ÕXVlECóÊe_8£¾çQÅ³ÜµUdÚÞ
(`­Ì7 k¢ÎQTT»÷»Õ½ÕxX1WU)Â£¿E4In©50òäîÞ³N¦cÆô\ççí¼F7öUÆ¬3îµ^ÔYímr§Md¸K½áfÔó=mn9ÁND.ªo±¦ÌµUÚcvI·e¯Q	¦uË±½÷ºÕ0ó
íOqÏâQ|V'» Þ¿xQ%\x±ye'Îv ÿþ|Bn3ê_öó"q­þ
á½ÔwÙOÅè«J&«ËÓÞ`5¯	ÙUJV*ý F`<Øô!Ë¼HÝ7ÖèÒZ®ÖN¬©]t`)%BTâ[9«nBvÔ9UE­qÄ«Ù`Ã~<f¨gÔpáâ«q8õ«3zHàôBtÄ9â*wÉ=nq\¦Np¸L-ãf_5®Å±mu(-þ÷y"X-Ee//cÇô]ò)nÃ	µgùqA´zþ½¤ZævAðÒÑÕ4²øµpzõG1ÿÐ [Zù;º-áéug7¿µuÏ;;QUgï0	='áÇCÈMÇòK]
÷e$I<Fö@~»>/ÎÜyvi¿B|Ær(\iû4
mµòÒî7BovMëEñ(H[Z\ôB5ÄÁø6¼HP÷J{:ª.¡´mH
ÁXÉÇ]>/¶Ûf>7:6ý-êÔÞè²ÏåiitÂhÙ1 ºnuMø«$Yfu'áiõÌsw^ÊÔ@·»ç­(¼²Y®õéÞð<ç]PÔÍHØ^ûQÉ¯?iFO`Üu¨õfôxÀ\
õ;NÒlÆ=éfTÞ¡|Óbf!<ÍØÊE®¥îhf
2Ù_óIµ¡Ôßg[âµ×ÙeoÃ3©É/¢NòýR:ÍÿJè@øUÚÙ"Ô-"´gøNð×C Ú,o"²ôK­7aÝ^J®æª)ùÙßî¤Ïá+ë»×æ'¤«9¤ØIUC>/¯¦[¤ØÖÄX¥GÖt8»
f½ªÚLa²×f]Ez¦"ü¤6CìµñMú-k%¡è&º7øçÀ³ö1§µ·Ìÿ·jûeoÕ>ZÅÃªtE"c~­ByÍ¯IOèçg?Î¿Ìÿá;!LÂÇðW'×,{Ù'Öi¦yïJç´E4PÙ±uÆ.)âzIXÜ«â¤ÌeoÒìRÀc&dIø1ê­éeãGÈé¸owÒçðò4ÚV1U\ªÿ;^FÇº~h!­¶GÕF>/¯æ·_¤ÕV)­FC/ÒÊ/Òj«ÍV{mäÑµQ¤gÒ*ÂOj3¤Õ÷OEWcH7DXDÌV3©Ud¢IááKé4¿c=¡_>ý(}Â"}zï0	7Â_ÐÌ^öI¿d|ó+ëMê×§¥]gcb]úæ<¯½zVåÕz©âÄÄ±]Ôe¯åüñ%móÈküv'}^_°½pCÌôtV11³Ww­å¦ÈPµÑ¤Ïá«1£òþ »sÊ1d®"=ÓR~R!Fö¼B(ºC¾¼!â"D¶jI­"içØÆ6|)_¢æ×ª'ôò¥³åK\äKÏð&áãFø«ºÙË>	<a!`2öÀnÓnb/ÎØåK|9çäbJ¤8ÁWÅÒbG\	CtÂ·;éóBøJ|ÄE|ôt~tB)rÎtÖu}Ð®æ[m¦(PµÑ¤Ïá«y)±Õfj)çÚL-e«Í"{mDÐµQ¤gò)ÂOj3Å^±§xÙSÅeëÔ^ç=U\¤Î[µÓ"Z¶jIÿÕÏ¸ç9ã~på@zâÇeÁ~?vkËù~«q±×"\±#~ð+äÆþé@}Äÿ\öR[ÈkÀöðfyþÙØ­1+ûr¡p :ê\¶=/GL»Âv5Bï¼ã{6Êý;f\µ,ßäðºì{çOHÒÁ$¦NÂMÓ·ÿÖ"uN¥ìiç~7=ÏhJª,¯ç÷'%ÑËLò~iG#Ó-uX ¯Ý==Q<kÅ«móê5ED: qNò({ZÍoÞ£ûmÄ0¹y4é#aÔei`
Ä _ëhE¸t²÷ç] "Qe{cTúHGhØjBJxãOEý9xBcÍVhSéÅøxÔhim Â}×%¶Ìª·Ë°Ù¯5aiÙ¡§ì% Ý¥±úÅ¬De¬²ÌNìÆ\@Y\l[w`!p	·Ñá*¦[éPy4ñlèÖÉ;}¦ñÇt_NBïKüt»wd±F¼¶¬}9rÅHßézzÃ{­'6oWÞÈKa³OØÑ²xkïpÐVÔ±«ð84â-±bÜUZá'(À¾¹{ã\RÇj;%íXfeg)ö$Ó×aà-#øaïyÈöKÔ(¥>­²I¨Lf³'ï;Ìáe°qYád°WÐô¢èSä$¢ÔOÂÇÐ³ÂÁK0þVµæ` ­>
¿o{â%QÅªì9[E¶8Ù0Ûw`Lok²zº¶£Þ¢Û&>/@)
¯Ãtsb
lã!E4sµ9®dXÂlÌ©òÕ`¬òq¡ÏN(¾$Ý³02B§t«µYÂÙ+Ìànä?ý¤+'¸èÄáÚÎ-ðÎÇ40ZáõxkoZô t<Äý­¨Ñó;h áÂß3G¢íâ
vïnÞÁ!ÛþU&hJ¾÷Ci±µâÑ²ª~X¬YæÃ$Tâ=èÂã´¤òÚhÏJ¹¤Vöíèú	ç±2lËÞ¹<Ã¨fí)Î©½§¸â<{Vv§ÅùtÞßÄy{Ü8Ø§81cÖddù8ijxJU0ÑÊrÝ|qÐâ±Å¸)çMÆ{OÐZ´^4ô¢Ô&>NDðlVq)_¤[Qd^7G
:|dÝÉc;40µ·I?.ô)\1æ0 >!0OÎñy´kL6\ÔÁ85&K@·6XöâfYîÚ¹È¥íán9 ¨\4ðãBÕèzNÂSÏÓªÆh#åTO`4
c/.Á:xÞu©§3ETèUqñf]ía°ò"±æªkÛ­Ôz5Îqk6áQÌ(WºëPæ[¬z"âc
õúb²N´:EÑérB¥¡ìtå>.ô)kEy
Õëe«ãò24×ðMÙ'rNX7Â@È. »ìV|O0ïcö#+"³ßPcÛþéÜÄxÝEí@xàU7 Ïg=ä(Z§{CúïÀp\YíÂ²,ôÇÃiìUt?ìcYhû<ÎÙ!1ÍCÿ+¶$«`£õOjD6VÂs;£ÒÃÖÛãåRK¨ððßS5°àW¥(*7^ÑYÕU(60ÓTU}¦X]j;¦NøTB¸³§î¶ Øï¿Á¸HÃÄmeÉº¸¢CòdðÖ"õQÛrgÓôhõ$r\Ì
8g¹Ã~ý£#¿ç(@¡Iµ¬kïòRa DÆ±ÅÞáéÙêñNÓáñÓý{DÇ&r3oÊ¡;8¾5FÓkkãsÕøCãÈ´M ¢Ï#>®²®ê$<UùXxd	Ñ¸Ï¸àAû©#ô*®ìtQe®õÄk]EÀÜQÔïT=Ú-ÀW«8
íü#æl+40Îw ×òwî³ü}®ùÛú[®s)Â¦­ßV<î9»Ìï]g)Á¦sûnMûoI¨_Úîß{ó^]«^ò×?nVvüHüÅÏW	ûEKÄç, ó¢Åü&Õ%¬Â¯sþñ¯¦iÔýÎ¼tzç·M±YØÞÒ¦ÑrIa¶íçJ2¿òÍC4F^r¿¾?eüóHü	èûÖjÄto½Îò=ÂÜ{Tq|]q°¿ªóÏ+õ½5UçÞ2mqO~½/ö¯ËxKÌl©|ùþý²'Ì~qö[ûvê{ûþõâÒÉ×GTDqýÓøþð§ÂzÔ÷[*T:E6,]îü¥Ái~Êò¿ó[^©£
Â?/UÓ/FYÜ9÷EáÚ9ôÒÏÍøÎõêîîþGïÀ×Zú?0T5*q06+&[LNû`Òo£éD>'´YênÚ	*ï·î½Ó+ñeï|^j6ÿ¼ÿG;MuÆCßàçLæÝèDÍñÖ'ò±ñDïBvj:Í8P³èá	²;û#®ø
hóèY#Fa= ö1d 11§S¸¾jßN¸éJ ÁCyÏâ¿ü>°ûJ¿Ê§JËøÇÃ?k5/dôÂ¿½ºöýÏË÷ÿ®æ<T=ÞPjZ>ÔEÃ¹kJ¿¸ø¾ëjûâ;©_+7ÝëJ¶´äß.T%X÷Fêã½[vòpí·Vú×WbüéõóhÔÃhïÃLùÉhF¶°Eö÷ï7þ}´?T5 ÞQûêÔ/ÇÏÚ¿¬j>ìUsP¢ÊïþRÉN£ý\)OihQùâÒý=õ+îJ8öK#ýkÉï½wúy°Oêa°÷Q¦¿<;ÿ¾þ¦?ªÒ¾ùmü'Ìü¯ðÿ)óõÿmp$C¸äJZ6ñ'­¾+¸`}Ù.çßepuº}þì×Î|}xÊ÷ÇW©nYþ;ûºJÑhñÙ!)þVÀýd#ra	ÍJ4Þ.n6×¬Éß¿ïçOu×?×®Ô0o½üÇþap<5³*ÒO ßOU¿4ù·[,jjýêÈE3O©x¢iýÙÍ
Ì/3U¶sB§cÅõo Ã¿8æKk>´÷Öü5®GrØÚ<4Îfc	Wï}ªnô¾yçióãØ8ûÄÿÎFýYeÿ<ÓßZ]Íþ3¶ïG¶§¹?'_õÊÿPì5~Sÿô÷ ç©ð¶ïG¶Ûç/ºâÛÓÏ]I;vßÏÊø³·ý®­Ì¦-$õû¢YV¸$gr^v úëw[×^´}\½!þ8Òv¿²`gê/
àÎ}Àøå¦rüR-M÷uYz×SØÇÇÃ¤f×{Ò7ÔíûißÚæ_IûÇ±CyÐôP~FÓþû×C ?C#ûóøõÌ»
åmïõ7ÄGÚWfWà¾Éoºå?/´/ÇñÝ/è0©Ëæf-m5ßÀûpzcûçÓ ¾´ÉßOøÇ±_£wÐôè}Ðþ{ß
åÞjÛþø´ì÷]Él%iÒnÕ­ëÕ¦Y!u¬Yz½	ÎýM^Æ]*àpß/ÒÜzC¬ù¡øf(óãH>z"KÍæ¨½¬QÓç3VSy÷­ÎÛi_@¼x_eÊ·2}^¨CF_óÃ¼¢¬¼x{Ùÿ^ókLÎÉwk
ÍyíÊûz^wâþGgocõûû$ìgÃ~äúMä+áêÖDV¢Gnu®=3txF¯.<ðÑs¢¾þ\öRÅ)Tüô ç,Ñg¡<o'û2®ÔÂË Ü¹öpîÔ_/æÜÀ¹ÁÜky:ZËüþ­ÏÂ(Ê¦§üãÂóPW(Êc6;Òîß¥±JóûRÉGã3%w×añ-UV¸6î[ÁnÀa~ÜQx `ÐÎfWÌ¡	^Ðëé
¶«éZO/ãPsHÔy\v N÷c9îôbÀÛz²ð©	Åë~ÂcÆW=àÄ;pÈµdZúsÔ#¡È¬
Àëå4RÜé#ÅmÅi¼ªÇB(Ñ>vC±âR%Æ¥±&Y^ÍN|
éþkåx?©;a¶ùVñ±ÈE47E¦½+~DNlR¶&tÑÕ¿x5¾?¡|+Å_½Ìþpø^ÚªÓÁ6ÒiãÇ;SLñD ¡LEàCYeÁëÏ[ÙAø¸ÆÁ~YZãSpßnÎa£,^%ÔÁ¸O_{('¯dS¼ÅªAjôURòkªÃH¡dïë^óè§[;å¾
>Ôhì=ÛY¦>L<Ä¡PãYÂ¢7¾-ÞWð2âdS©
Ñ4Âª	_«2ø	òÚÒ4fY3³0Ë?¯-úÜiÉ²âU?lÚáWÌÁ>ÚuÿmÌæ¹P	ÒkÜ3YEÜ¯Ì¯ÝÍý¿+°qEptDºáEöËÀøÑ	Å"7îýóMw÷/ègÀÜ§È¯¡øw2R%û¹¾øâµþP+gòMclÜ¼¬´Ëê5êçd
ÄÒé¤4Óxi$j-OÕßVÅA&zXÔß[4ù:Qr¼{)`2®ycj
ý®	Axµ7Â@;¡#dÜ¾äÀP^Ïâlã¯
T
îµV¨,ã¬YÛãN)ÊJàpµêU;áÑõzúbÌÉñ¤ú]?|H5®Ý6x^´ºH éèEÔ¤´¨Û@ðèêqÕñ¾yéÛsü¸Ðg]ú¨7BWñÙ3alµçYÿÝ¢ã=SçåñM:¾Ï¹]}×<^\p}w@+Ô9­­¶Mø¼\ÏF#ðÓäí!ù.ëð!ÃáGª÷¶«ª/1Y7"ü5G.*)Ùuñîµ*·÷:MèIµhÀÉâu+ô·­A©opRßÞGZgl&Ò\97úTÄNÄ¢8õÛ×Ü~¸7+·Ó½çà&âo.ï6-§Ð&½øçZ| ¾Î0ÿF.{©D$É8_8v;²ÏuøDF¨Í£Ëëµ×ô©w{AÁNûáîUH2ç÷ÿÑÎêï¶$â]îé$òô³4s5ß)¿\µRxÖr¡pp½,5ÉG:ýûR­æµù-Ó"wãDïñ
gàÃùÙÓÜ
­1g*×eì×%%Ö­ø5GtÙÒj Ô»Ú²g²Õ­²W$ó9GÜ=áDKx9·Ð@ÐZl§)a6ÏGúHQ&Oj^
XË`ûÜ	ÅÈÀ6åäÊ­1P¸ñÙH4¯æyd«>¡"ø"Ry¯gÀ³Ý[=É,#o)õÇÐSLF6·!\û*É^XÙÈ(ëÏ+>³nYÙýbôyËâå¯Eápó
X¹L0ù6l(óä±³ÄÃü=EÑsøz<u Ñx ¬R¡@'¢L/r9>FË\ovàNê×*ÏD
Ù7ý7öxMÈ<HNö!ÌÈ¦ÝHFTFÆDTñ¬pÖ?Þ7Þcñ×vl±£Ì×yhì!ÇDìx\ÇlÀ}J94í0Úp!àw6s£¬J¼ïþÿ"ÕO[GÂß¨àÃ,<3çFâRÏ>SU×«ªðÔóàæ%½­öTÏÚ<úéÄA¯î>já}»îO"s#\H~Ã
rÖ@ÙgÔ¶-?ç(»N%Uw¨µ¹ù¾)RR£X)]8q ²Lv/j^DlG½*sïcf|=¢ÔØ8!÷*Ñ-õ,­?p¥£Rô<$ªÕ&0à/éæö¥rmî)ÞvàN8Ê½®O¹·»!hrOVÈïºõ£PW	â5¿wþD04«¿sþGýÜ¿ÿ~pBð^§QîÕìÆ.`çñKÔã`6ÈÐÏÜ],[ÛbwÚxeÏÜÀ@ã"¦¿Ñ|¯×æ÷À 
i~øSNTÆ÷o!8{O)£{ç¨R*,GÜãoü¸þ	¹Iøñ"à÷?8	¸#Oøbï£	7fÓfFÅSqÖ³ÁÊàËù9W:1[â8Éa=Í¨_ñÒ\¸¦óûÃpøÌüþ­#¾ªS)¡äð|O	>¸ÛùK-Ã¹¿ZQYRVü¢8¢+ ÞKEóõaÆ½ÎCcY¨yxâÿ«	
ÇèZà±(5©/×²ln®ýÝæ88æ`Pßög[Î·Pºkøá7íUÐìÌ¸úÚlÿõbÍ¬`TÐY^¥ï/eØ¸õoË
&ÃGÆ[Ù`¾ÃD]i«Ù+Æ´ÇKnÅâjÜMXÿ.*xA''£P{rgrO®â`Î/mñ¹/xÚñÇÍk÷~¡f7øº?(¹yhÜ*ZXÅCÒ¥¹j
¯B{A­4}9¡	D@³Ù3²ØYðª{É)¹·F½,rP¨-UÄÑ5Ò}2¤ÏÂ­njÑuRÎ&²÷	-Z»u«¢¹¤ÀÞÞ±áëµ¸Ê-ÁÞµBâåÓ²¾EcêuPyR2«(:36o-çp¥(ÛÅ"¯·Nv\;´Ð@QÁR/z.ª£ÚÜî-ëDü¨#ÞBRTuø¸ÐG"íM\	ýÃa¶vs¡m |(cY%«Ým<)úL±ùã
j<-¯Ç¾m}?ÆS/UE\:]Ez®L[ê¢ÛÁ?Ë²srîSÑÒäÆ­íðZ&P6ø»R]¢	Áûº|Èªyñ´kEczÑ	å ÕKö·9á½¨ÕÆëulø}¿á÷]ÇrÜæ°tt·uø¯»ö@|Î>VQ°Õ¿÷UT¾Ëó;Q!\ÍñÌ\JNj­Ñ\*¿økrí¾(G+BÃ«äd±ì»÷ï|ÆýMýM×l*!ý\RÉÅÅ°dOOÖ²ýJðãB)f×N2áSb#^å½Áë-[©qÎOYIîµTÞ0¦huAWÇVB{+Pd¢mjÉäl>ûIØæ'vxäÀa dý!+é¸^¶|?³
xDO&Æ9Âk©#fqÆÜF'a/óè´uò
,DJ©:æ¦ÓÃf sðås¯C5mo®S¬¸`\²z¥8¥9Î·Úc(Ã{þµSú·~+5NZc0zMQ>Ë.J÷°Y¸Ex¼l*`#¡ÈoÀ$â¨-¿;PdTiÀ3âE-\§{QQÃ:ÇÒe;ÆuÊlû×H¾Í­¹°.%
$XÌkø¦.Fw°Ìüès* %ÌÜCª¬	~­,y®Ñ³L/2(;<âXáz§ND¬ô~øÙ
Gø+V¬ýµµæ¶µJ¤ù½/U8Þö¥ÄcÚùeHpÅ/k¯.Îïïc9SéTY	§øøë÷Ïý;Æaã1QkÊoXÖ÷èu2öf³"ÙUUvMøÞ	­Ñ«RZZ9YÝ<¸Tt;XÎÛøÛ	Òü±-6 já¾OÂáTªÆÞÁU´NQ~MôyÌDÛ´=·OÍw;I}?I¢m¹æ^¶³¼)ßcÕÏO¨ddí°Ü.Òõï²ÜeÖ­í{e[ªÊ¦J°e~ÏÆutz6þØÓÃß;J±J×û×ðË84R$c0ñY6RQESóQ®RÄ´Ûü½rúL±Êvåª>nU?5îí8/£1Sá¹ÇÇJQô±RÍ½ð;}¤16_	À
WÏú'X¬a0öóF*hâVÌI®ã+(QS¬ò:·7 èÖ	Ò	Íqr8¤"Ó6/sco¹>S,¢ò­éãFèYª!pâx)|uÆÆSá	5Ø,TØå=>.ô"&¾¬OMø²ÛÙ©P1nèS¡re¸:=dçLAÊWPþ«ÅF
t\rÖ]~qd|+JúCuû}Á©¶óEÏ]×nÒuîra§ºN×Ö¢&Ñ¯mÔB×ñZÇsÍ1!$jnÞ¨5OT=æ;°Ä(ûCs5[¢FVÂÃ@õ&'SîUÑGèiëÊñ©	_TAÊö0P¥Ñ<çº¸IDw]è¼3Ä¦? e ÊÎøÚW²Ç%"pÖõiGÐCMÄøXRÎOë!©¡ìÛ	S^u¼"ºÊ°öH®Òö¡Hs[v=ÙÂÖ0¯ö\Ø+XÜ»,µ-N°®¿D£!©äZR·ä|á·¢;¬HbÖU3y$'
EâN½ÿÜ¿×M8Æ4}i["ÝÚ
óâ´®× 
Z¦éÖ\5UßÙV
¡Vn0¢v¥¶3Úp¢¤y"B_Ó·`UË²o2o'HUÝâT¦dªzsH*0YÝÏ¼t(µì3âkå~Jí7¦âÂ¶ÐW
Dõ+öZÝ,s¥]ií@¢ÄtÆ(!¢ ¯õ·k%ãë¡¨UÄæÚò[ÄÅç:b çz-ªcE7·¬!ÖCxöï«½§Ù ÝO²Ý¶]Ðnÿ¸Ðg²-O«Xø¸ÆRÉM9õª]¾Ñwµ¤q9¹ØnÎ.µná¦"m±*ö:EÙdöY:R	 õìYÃ¥îuz{c¤ÛQV$ÑlØ«ú¶Wß	¼D6æ´a.êzVÒâ­Ô¢i-uµ²JV³
¸«Y4ç\ÍñÍñÀr$ñ*³äaï|~§ùrG&|>Ë}ø¦\F?Sþþ|]iÓï}¬_¡?a$bµñ·ØõQóÓ¯í7þÊóÀ¤ù«¬îÞk~Yï²ÈÌ¤ ÜhLáãFI¥öQµòBøºE
Â4Ú¹EkJð^ÎöTM×^$@ÈJ<]P5S±ha&T«®Ù$ü¤f²ÛC8ÕÌ{Ì.§GÙÊIGöë?ÍPpÖkDÆR9²ziØBÄPäfbGiÉTÝz0VÅ(JÑ³o!ü¤õJ{ØèZg¨²~]H¯BðÂÐë9ËÛ	åz$É~·þÒÛMý¡¢
ëfË¾n=ÇzPáøÅ:Bþ¸~çÅ¡°Gfqí¥ûxÙ³çØü´®Êr> êë øáéÂù|kHé¥¤Ãï4D¤¾ºE¤(ÙD«[$DZ&ñN#Dû¶í(gG$¶yÞëÆU¹KÏí
Æ üxª±õÖ`mícòw`è-Yç2GP2¸l6se?\jÆW%ê4¹ìßþarÙ¿ëÍ¨~O5Å8£RÊ²73ìx6~sZtùþýÇë»ì·E0G<4jê³ù6n«èÝ¬oÅÎ¥U]µ0ùõM ¾îvöí ^XÆAýÃÒJÃwh÷%6±µeò0lÀú÷Ø°¨ïÈy^îÙ}¶éçld^ð+'êq.¡Æ
8ò³ÍiÛÈ¥©î£Æ·exú2µoÀ û¢£ª:äÄÏÔ)²Ø{Q"
9coO1Æìc¶pNæ¼ÎsêêÃ5+ii ìypfQ¢?¶ð@;íE¥Îí@Ç9hpåcsÜw $Èí·~\èm99sÀ6¶Ï©Xz1Ünþ*+ÅÉss;ëå+.Òàé÷|DÜ[Ñ*4ìCUýË»eì0&¡{>;J\Ò-&¢­ÖµsC	<XLÖ¢ ®Ö)¦Ì0_{­wrU³´V	æÚëUÇi8Hlúxö´¤Xè
OÓ¬n./ÒÁðxOQ4Ï"
Àuyiã]¥µê*Åvàü;p'­ã²xESu_
»¬ÔÑÑÙ{í@AýH¦y#e+¼
uvÂúwü¸¤¢¶C^J¯Ýd#û(CÒÍ½c6Ç.jzBl»r«@sPK#
Xò\ä6qee+«Í»©tç-GÂºb4d²1¸æ(Mç¨RRâyÉÅOl2G)7 l.E+±Èâc½Uä¦¢¡°^sÄ}@9ndÜad¸$:¼)k?ñb[TG{ ÒÜX ]ê b½R)OsrR-Ê3º¹a]»æ(K$
t²±å¾ü#Ö¿DûÞ[MöÜ6;Mè¢Ceµ2EÔ QcW%å'BØ£¨'DÿíìM9Fyä¡|vó{XÁ§øK<lRüØiä6·7m>ÎÁòk°%áÉM#DÙ`kn5F#vJaAüxZo4mo>6ï¢¨½ü.8zÓÏ×c{ãîûÄ´(ô?ÅL¥3½¶.'§AvtÁ^1{Y¯dÄ.{y0nË2Rý-à«ÈÜÙeôóÞÙÿ¼°O%R¬õãLÕ¨§¶/Ô·+îÛ/æDÎµ?Lõiã1q»dÓ
4	CÔÈ­>.Ü3Ímkh#ïÅ{1Ñ¬0×Elã®	[Jÿÿ*Òh¥£k9ºãÒ'rUk½{RÀaDa±-í}JÜðÆÃKj¼gîk>o
xKIh7½¬Ð¼h¿õóFß_÷´Jêï%£¹¥§iWO==I?ééã{OïËÄÐã3'ë0RàAÙHºè5]zÎáÁO¬¹Ö¥óú0Û	×aÖU¿fÁr¿NÃìf8Õ%X×-p¼Âº³Þ"<K&Px÷jC<qó¬+ÓXtÆÔü<¼#Av"Á/ÿèIm4ç-diâ·$U6Þ#ØbtÙ2!LÒì7%DõrÏÃsPP$Ü6Õ<Æçø
LT[`¢òÌØ?Í-iÛ{Æý±²vWß.#LscÚû}#yØ2¯JÏ)àq°çX´d$ý£û[%úÞßÜttïRÆ|sW£FÔä½wKéþ¸¸Õ%;üxâ2ÜYÉv»LW¤EÇ^VùÖøüIf¢;ZÜ;¾2áúo£chÆ¼Ù¶¿8êvî¿æþù¨EÏOQ·¾u'Æ·Q§æ2áD^ájK
Ù;f#ßÔ q²Ð9éÉÃñ&åèMÖÃJnËãnQ´À¯	Üj¶½L\óÙ	Y3e9
Zàq¡Î,<.Ü7à­O£¼åV²ê6ô<öÊb¥Û­[«CÎsm­ÄB¼
Bkøä².å.j±!ªA,xQWÈªgµªá°éVyø¶×=GÍ´á	ûï'£¹6aK32î0aßH_MØ3ã2a}"û%Ýß*Ñ÷þÎ&aÜª{/·1/êSj_#Î$2'xeW$(´´®{Ó0ºØ	W	P"´üW\1O)±H[>JÈjl>eQñ&|Xg&~=lUR/ÍýÃ¶rÃÀ§¶ôa{`|¶;Ã2l¥±ywXï½¬A¤÷LVtÀî÷'Qýé7/ëÇhïP÷b^.-%úÖKµÿ¹«ñs;âp}»ÿ~ÿÛY^£=ÛùoZQv×Bæ·ÿn°BÀtUQÝ[µ;¬Ïÿb*9lw²xN××;NÚ÷£é
ôÒ7¼eT@å´ÄP¦aó§ÇdÖmÀ2""±ÕÀÂ°rÍ1Ê(å)F<d¥x&^^»¢W&_9¥6¦4#
vJí5°ð8>ï3GwK«OBaZéás£zAÈFÃ«Ù(sE÷¼¨kÏ;®§d·¬Rú÷¢Ê#ÿ×È½,®W`e5x«Y¶<µ¸k*öTÔæ*Êîcm6x°ÒT@ÜL0¹åè¹`®VÕãÍ«ÏD5|Þ³B÷"ÔÛYÉ&³Ní=Ç Ê>ÏöæìNq{÷ÍA°FéLÑPÖfëµÊoÿ{@wO½­5úþä02cü÷äq7ÔÕ 0nÌ¼ám/yµ¸Ãàv6'°v;Cìë-7i;®¬â&1sC6mYá¬LîjÅ ÇN6èqíýLÁ%Þâ¦¦¨ÍÕ§|(ªhâÍâcBHu¥Þ¹6ºqb£W`Áú³Z
¬¸%«;0Lù'"LÅpR³"T):©[´KÕû4Öâ¹h Wó)Ø[Df¤öBQEº×AôÉ­ÁÍ¯?.½o9r©Ã"ó¬uÄ.7±¢;Øe9^;;èb­eLSo@neMxbm:7o^o"À8¤ºåèE¥Mpî@5e¹Å)XwàN8KÜ­«SWÍ"Äâ·¡yóÌ¯Mvã{q"´îz¾ÿx2%zrÅsF¯½ÄÓ}¨"}BÐÆÓQd}W¾Ã´ÉyXÆ÷ 
ú»¨´¼}Ï¦b´>¿ÿx¾K/Î¾z2`)Åk¯¿Ú5}55Ìw³³Ä< âÎ;Ð
 _ðh è©Ý¨oÀ}¤(Zô®>»¤sÌRõå¶Eg`@%
ý¾ù7 À¢uQGcb´Í{©\ââg e2@Ytóð9¦8ÑÇQ8CIt*©¹:¦h>£D@ þ>î¿Èçý¢¹~U¤êqkÛïui´n·x$öúã­ÊùzN	GëuvÒºø1×Y¥M	¬
î&×áòù}ôä)R&éÙ/}PØÖ68}±ûSðòÊNÁ³ôÁ-ÄßYÒ@I3þ½©{5k-¨÷×Á]+AZõtr²:âà6ð¥çqë¨¶9	wjw L+¦V
HÄ¤®@íØ¢Ea±Ø0OÂÇ0Ñ:ì°Uîaùq 2ÆëÈ~ÞðøCVÒã¢¤-«	Ä+{tzjf*ÉTçpkáßgC-'2J±èq±­R¬(ÆÞ&¢qâK/"ÆøÖÜ:Ùc!ã¤(ë!ÌÂ(+FÄÞ}Ø]I_k)Í£%D×µËãc¼]Üob¡p[-½äuh²µ
/éx§u DÙ½&ç±·zìÀîiA¥8=Å05´EßÁ©LÅ¯=ºèâÃÞ3xe£¾}ÏX ôïßï®wþ2Ã9ÙÆ/ò<Å¨ÓÎrM¥ßù¥¹h~+s7j
`mEæð¸Õ;8l!È6«:k{OÝz¼ÌoQ#ùQLÕN?pÈ©ø#nZô«	«Ñ³ÒGà+ºdúÊ³2MHI¦IÐy$Á¸úWÂ<[w¨GåÑõ 8T*k|ïop@dÈõ@­km'9;¡=NzØºöRÚë&M'Oõ.ÏÍK&T\ú,Iýxfä®}JÖDøÏg÷r`è³y¾4GpÙÓyu4û¾{ÁìnÙ¥<,¯RÊwîÍüÚêáÜ±Ún0Ü4U7í°LÁRçé5×[nE*Q°dho¹ñ\Ç;¬âò3ï°m@¤QV8öÞ¯;üs~S¦cwà ¸æ_¸-µ1§ö`/Jßl·.'L$T<kë÷Ný/âuVlÀÛýÊÈªJa ñÑêü¡V- C:ÅJ,ºn
èèÒ>ëaÅ¹>÷ý·Z#j»{Úâ^¯Ö)F3âu¤rçÛ]àbÂ-¯Õ]×Ñ@x-Un÷¢û°yÎôb¼«§-2lÖx®d]ErºKqèÚÞé$½[þÇ à]×ëñ"TÜ>G;z|ÃèòÜûYw"»y g¡Ç¯­×Y;iqìpuËQfha¿·÷Ñ³UÙR\KÍIV£å6"bºãÞ¬Ïóï;2ÙÇwÙoÿ¥´ß»à5Ú6áödK¤´<Z§EkRàBuM§y¶ÉP´Ù
äÀ þò§ËïÁM %gÕwÑÞ¤ÎÇ"DÑÇ>S¬²[ÐSÄ9O³u)âíJ7¤±½öÄïw{57T­,mÿ©l£mä¨´é-k0¸ðrÜaS+óññ¾^þ9Rxë÷Õ{.(íÎ ëê´3>ßØ{Ez~Rjd^3A7Í ñÒQ¥È sAÄ/ôQjJ6+:Ñ³sYK=VLß©èÄ'á'smc^ÿ¤_«î6|r¼ôÕmQd·c¥Þ¢­jIiÁw³JÑµwnéD(	ýUeE?"âÜæ»è9±y­S@ÊÍoù"ª²,àQKQbÆr põ4ÃÒ¿Ö°mnTá.ÑüU9UÛ­J±4wå¥£)×³þ
(+»wK÷vªÅ»µ¨cq3":ÛËËsRNÕ>	3¬ú\ EæX»>n±ä{ôk«Óx& vTi¯Ü,
×^\ð«2JA<!øJÒÕs!úëú¹,ùp'uÙ=`ð'º@tÑ¶ËÃ>ÖvT°¼$ce©¯COIÒÔª`Ü´¥3,û±aá¸D9Ã*Ù!vØ\Ã]_°§gñ
ÝlvÏ	â7fÎÒ0±ÐØ+íUIÅ±¿âÎd'äQjÌ»Ð¬ËDGÏ¢\ÊÀÁbXËÐ¸ì#ðô°y/º(Ä³káÇ¾)ªó²PàFúÐHrÙs+Jµ¸c9Ô¡VâP}Ya¸ÜY*÷©-b`¾­c¦fÐò0ÎeXvtª¸xÂJíÚa­ïPmn.Ë¦¥øRIÆ\¶JïLW¶ø¸/ÿQO@'jìºÜÚì\ï/»ÞÔ¢ç)` þ:B8æÂYËÏ¦ï)âÕ²¡ãpfÓdþb2LA=¼¥S
ÔtVÈyé*åÒµ) ·k9ú+°¦´j][hx-
äg³hÉ_û¡sV¾9Ý_ãXÉð­dÜHyàª¢´ã{HÑc%hï)Ûô
þYT/À± íÀpZédÑWZé°¤·csKd¾b÷ë´9²Áò½VXñ¶ócîË~Rñçæ;n~ÿÜ¿3Üs½ÍóD??ø¦0ü±º$p42ÆöÈÍ?^HdùøþqùÞwÎ×
:Â5I¿â$.*~Qé¯øk=ûúýi1â%ÆÖÂ_·X»KùÐdBJÍnÏ¨ý½:»Öv´æÕ¤&8¼t'M T?4ì¸ÀÔ-Dìc¡<)Lºo,o¬q¥¡1ä0Uü÷s{ 	Ô9§"¤:²BàúÚá(Z#Ã¹RÖBH{(tÃÓ ìacÚSl?_i8_^=Ú+ì [}0r/\v1É°6YW~®}E ÇÿÏÿâ$"_Éòyy÷ZqD~ûYþzÈÎÉòlä÷N:	ÁRÕæùÀÕkè
!º>Þhó3Ê+ ·8~ýwßo.6[Ò²ôáömi^î®Ý¨Ê	gÉà½8 1ìî
<À«y.ù@`éB íÍùºÊëÚhNÒ;rÄRÑFK3¬ØÛ{!HÖ`ÏìÐË»0@§¸2ê³½¹¬ÅruÅ³¶Eõw^ÔD²®@sU­N*$c!ª2.P¸óÀc<À1~@Hu#Eâo´xçUª	ÝoÄny¤õ^mîÅ¼4>ôm<ä	nhí:¼ÇC²)¡®ïÎþûñxw¶h¼;[JÔßÅiË*ë¹võ5×d×ìÛå±°ý¤,m9à;¯:d§3 øÈä\r9Ê0ß?÷ï¬çxîFo|Ï93`Éà¬:vÒXqEªë!RH.l0Ûa»ÃÃÈ}\¾÷ä´8·fy«©âSÁ*Ö»
6
oEüÊÞ+ªæ´¬ã¦SÑp¥r²¸·ÅfÿRBÑÛS`bk§P4®²¢U"Ì«×Zñ6aä§¸½9f_>¸VêærÍ±²çäV0ïEÃYÑ{NÊÉ³ ðdÀfîßÊØ#+ ôhÅÁø¥UpÄý¢âìNäú,¦>VgßËdtåd¬
u¢LboE!2:["*\VÁ%.Ï/­âñGÉàèo®EoEñ\
r©ÄZ;¨çøq¡%u5Vu¹&b!ji!«(Ká~'ºM¨84ê¸Æþ|]DR·îÏµÃÑHÒâ8ÒÙ{Iº°#×¶jC \U?`]bfÑØP§¶äØýñæ!µR²*9Ù8¦f# Ï5.ª09sÏi
FdÇ¤SÄS&*EÏQñ·²sÿþ°§yÈVNã
c÷2ÎSv\S2µpéT®À,ë(6©
0d®9r&êp´×[=V¢C¶fmYuÇ½MÐá+òP÷EQ¨%ª*
»õ²T~J0jjÕY5kc5§°Í.°»!á0Ê®îZFb×·Íû,\·F=V>.ô)áKQXQ@Þ=Õ<×0Ñ/ÇDÇM±%'*EÙ!áoù¢,¢ß¨w|¢«MCÑ&)TVm²±þ]²¢oî³L,VôçCä\ó.ev8`°÷<d{;ªÅfj4,g©ÉªwäÃÇ;ÌQ"ñ°£¥,¢¬Icn\ù¥6DXÿÍí;e±É½÷l
·¥Þ²))b6¥«4FQÉU^ár«låÆg]:jBÛÌméùYþw-WÑÍ'¾øm^(oÀæàTÏäª	¼oZõ"«Ø *ýq¡Ïeo+sA§(êú"þ¥ü\ZX[Ì¡(²¨6ìÀ¹&µGlqé·^FY~
f§ª2|S,÷¬R%lã[WËú³¨ÛÕµçÖ;°]
¾æø±gdqÈ¹6) qÓj¼
.Qd¸Ë_zîóBÝ+íy"1SÀöå¥V<)|¤àqÄ¯KeeÃíU,ðº½§T@¡¯óW ¬èSlr½7E+¥'d elJþ6ÞC:÷fÑîÀ¹6Ün,u+¡¡><¾hNÑK3q5C¾¯@,\Ó*ÇJê±3ò*EOXÖým¼q^*ã´¼ù`[r.J^:Îå-Fçòö¤xZÞ8®1éË`«Â]3cì¯#8Y§<aûf2óÔï_¢´Iµñ³Ì¤R?.Isò¿8ÙhT¯ø	ñÚÎ/Â÷½¿=2wIo²,
Áë~HÈ.w@WÐAø¸Æîí? ¡ç!E¥uR"p(|À\A#²ì[ê©ØSq¯¢ìQxXN*{±T4tÌUR¨ÅYB¢¼îNÐK
z¨¡M{ëtæ§¢Þézù+9m:±FDL!Ç@*¿vÔg¢YÜgj{/\^²N}Ïw~0b¨à\Ü®\·Ð\÷û÷U8¼6ïÎ5¾ëÙÂ¡Ê×wiJõ£8,CVâÃ1ë½q(³9ÉËö åÚ!z;ÝÎ}?ÝÎÃ«Zgy&ãU
×)AïÞæ¿ý	ëßwÃV}Çâ³T
«M×lJó^4¬DK´õLÚ äv¾ár8`½Sw búÏßÐóe­x¢ÿùP¯@¢ F¬=êÈº'íÀDpç®Àb¨¢×ÀâxScïÀÕÒ¡qÚÎ8]³Ydt<%³ívé}a¬öáäxª"¡úk¢üÞ+ L¿jòiGõ;Pfm7ðãB)VÌ;æØMò7fj9tÂË§ñ¹ð­'¨k7îDÞ`¶ôzÈÎa331@ô"BÂ^[×©þ0¥:gL7 lBDþ®±á-iR)DY£pMqØ¸¦j¹]Ijº]aCØ?¹(~Ø{xtÂAàHÆ)oí=¶¶Æ5¹þqÅF|Ë&r¬Êfcd~í0éP,1³¹ZyMA±»Ñ*XÓ1l¡6y+b1§¥7Â@$Ý0m9©Øy3	¿·Ð
f¦¸Ñ§Lç¼.UþTnóck9xç<³Í¸	Åh/.g_ïÝ=ÞåyÝ³¶¥X~6közÊÓÄÃ°¡ò6$Þ7¹þ ,ÄvW`Âè>ú*A¡æ+°9²²EÃ5m~©|%*×7pc¨i×ÊË³¡]òì@Y3ÎE.óUÑgR·ôx[_7Y=²ö0Þ¡ÛMæ^\Þ¾¤½¯fËÉúfÐ4@:ËÖkn*O7YbDzC³Íoe9(ÉèÞwYá¯½óìì¡¯
±î«±àÈk!J Ä×Ú´¥Èjé®tMQ¸¬¢ìOk9¤èMÀpO>RÄûxZ_'O@Øv4ÇS'(µâäWXû¶Ìó
X4á¸VÜù¸NËjfKtÉ±ò(&ÁÞ|Y±Vïøª^ÉæxÕOæé-4ø,ÛÎ×\Äºß¤>ö
ªAÝàí>äviæ 3É¢¬>|#ÆA[ÝqfMú*.
ÐöÕýû÷¾-«xYPü·a6`=ê{çgÝFô÷ïªzO¼1k}àÀñ4¡G/â¸8å»¶­}ßÖÇ½ ,¯bÊw×0Lh{8êLJûí0¬juÌÌuf)´Ãú÷j[ðÓ=¹J,¸'7JQ+fAAÁð6Û_p"{¤ÓÎ'2XÌòzùñ1	SJ"¡* 'uH×E4¦Qï@pà¤zo±}¤(ã+ã¾j>7B4Kps5µyî©¸ì}±k$tQ¹Ö³E«wUÛ]xog`áÖ»
äI{cHñÊØãÏur§jð9GÜÊVF}ÜÚ²;QÖÌ^w ÎùLõW ,:ÊU}kÁ)¦&PLæ+[«Diiîo×[[?.ôb2ÍNQö!>®Ä±Ê&g-+_Íµ
W4uú7Â@`s®ùDóoó¯âXEne¬¢ç²OÏÙ¥ØÙâUQDPØÜì±AH{a±eElãØ²Ï]*£ÜsMªÛÁrvËz²wÂP¥¦2Îµ~u­#ìÑUÉövÐLëÝàxÛÉ½Ì1ópæ>ø÷ï=ÊÓd«²Ël­PPæØ-Ù ×tËMÉÀC¤½RÒÖØ3¨ädUÜ:\Úh'Õ/4/Zóûçþp¶ÆÝä¿ÈqÑ¤Ýàâ¿ÊQ]ÐJlw«·PåÂgÇÉ^Dá`¶vþ¸ÐGD#ÈK
"/;éÚ7¼Ç2&Y
y1ÂË^ËÝ5E®IçÒ.òåEE'«\F,²Ä\<ÔAfEó]~«¡7°¾QÀ,ôý(]N z
}SÆWßíÃ\e<\²uYò¸¥.A5·ã­±ºò²µ­¿VÍCdåÂô,cs,û¸Yz-7yÓ"JÄË\=ßê@0z·ô»©;VÎrÄO¥S#c1!"mË5Gm¥i{sû×©W>DõÃl}ÏÊ[!]¤ÎU#N£`bÜÊðìÀ	@¸ûr¥;¡­cÿ¥V2ÏÕ]ÌoUzâyú÷kü«i¾9è2¸Ö
¤ùsç{íÏ¾¿>óã?ÇÄÛ-íµÆþ\Iðþ±_Ð¢º·ïOÿ<ú¾5ÁÃ¢ñ Ý[¯³üBø,ÖFFº4â`UçWê{kªÎ½eJt±ôëe|±]Æ[bfKåãÌ÷ï¿:Hû­};õ½}»®ØÜBM¾>¢bà§ßþ$î?öïÿ.¨t§×Ë¿ð$=9Í/êobÅÙøñö8õ^ª¦?^jm
FÛ¼qÚ9¦ìRÝÿýÆõêîîþGïÀGûß?TYQ:Û¾}V,¦YELãm4ÈÇá8Çù­ê¦ Rø~KáÞ;½_öÎç fóÏ»ñ´	Æmò·«¿G'j·N<e²Sè8²õMjÖÊÚÁýQ"Ô:ìwîÁ~iÂyÛÜO(¢jÖbCªss:+áË¡öí{®<÷¼(þkÉï»¯ô«|ú¡´<ü³VxøBF_(<<g|¾ÿyù~ãßÕªgÒJMËºh8wMéßw]m_|'õkåæ {ýBÉüÛåªëÃH]q¼wËN¾£ýÖJÿZòJ?½~ízí}é/?í­¥¢Gûû÷ÿ>Úªï¨}LuêcJvh)Ç_W5öª9¨?Qå÷L©d§Ñ~.Ï§4´¨|qéÁþúw%Eû¥þµä÷ÁÞ;ý<Ø'õ0Øû(Ó_ú
î³þ¦ÇüQöUÌoíoþüg"ðÐü¯øÿÇÿùý?ÿ÷ÿ_ÜõØuçwþ4{_øjçÁAò$,P#ylK (@ÿúê^NU÷¡d$@Í=UÕku-½T½«áÝÛÿúñýßþã»·¿ýþëþöþï|ûÅïÿòþkúã÷øÝ¾øÿ·ýÍ¾ýðþï?||ûÕ7ßxóÍw?Ñ¿~õõ{¼ýõ	þÓ7?üúÏþë«_Y_ý7ÛfÉe<¨ÆCæiÀàf?ß+B1¼ÝàüÿXd ñù¾}åüpöÿË[¥*ß D¤ô
hMF°&>¯'Á~ç3
!(Ù hütïúÀ=ù§áýÌzo	´+®óTª
ë,)?t³¡f´ð[zàqÕõÍÀ_<b|&^olé,»Ô]ÂxIuÝAz9ò­¶­Zlãh¬_ØbÇÉ)_Ä'òf4OØu4_uJ[ù°44K¨9×^]á§¯#pxp#6<þOõ¹Î~»_F]ÃþÚãí÷ýòójãæ¬ü?Íìñì^¬à× ë¥ðÚãm÷ùÆT<,h°ù|ëé{®Ï,&üöÂø»(«é7ÔïÉ²Xh¸7É³øºÆzíI¾ºP|¹IïátÌ¦¾P 7ì¥À
MÕøVMCw»-mê%Hþ8Ôêº½£Ú}_ññnl~NÙOËI^²riV>¹Iþ¾ÍÊ|ÝÀJV¾|]ã
V>`rÖ/_.a·Fp79ùb[~ØÀnò±¬î-ZpñÚ197#DÙ/Ä®rL~^1ðfL>¿à§å¼.¹·Ã4÷$þ÷p(å¨]ýï+hò¿¥×§oÂ©ÏÝ¤2"Câmw´~Òô¬bøL\ló;p52øI[<á@®!<¹¨g¤âÅ±]ãuIc7,
ÅK%G®7 ÊK(ÎSäXÆº4pñ.'§³mQÃ§}¿Ôrè;ænãFsàÖìÿI¸/Drdäèo¶0`!3î¹u´À~ÓÐW,N´xÒ'ýaË	F¿øHèH¼+&¦=g§ãd3Ms
		ã¿2gküiSPC-ÏÏ¨E¶*"]r.ÖÙk2n4ÊâAâ@økÅSìµÝÚ\ÃqÇê¹õWjq¡¸hBTB`8ºÖg#Ö%¼3­NÜSExÏ^öùÄñ<]}÷x]ô÷âÈõïmy[e§>á­AÍµeúYßdq¹ Ð®¬îxù²ê,²çìª±Ö8cMRtÖ"`àÔ&,f"Ô¨wÐªµÛu\-= ¥m
njª;;¿-1 01ï¨¥|GÄjÙxMÈy^W¿ÛÀ{DYjlì>)ÃäIÀãØë1®ÉbàÐÛõ­x#¦øBê~È3XbÐÍER$4R%v5ÖEN!â}2%¼
ÏGò±-ÞB+-÷uøIÖuìVEÄq¿¯£òmË­³ÜÅö8¾ßËïxËoên|¨¼éxcx¯x¨ðÃEÕÐ¶ ¼64PGX1ùwqÕSVS
ÈQÍXß¤zk6,á
eÞµì²×­°HejÔ ¤¾I·O±q»uWå$AA·cáFéÐ£!ÖkBX¨n	±¢¯Eõ!  °×SI¬ì#ÒÊ÷T+NËïÁfVÒ©¸´åV%&¾7«øuÜfX\ðÍÈ§9qVñú[77#o{Üò¬Å#4ÓÚbÈn±~éAV%KIÑ-&)|­.Ía`tQaHça¨´ï59´ìßÛ·Û¦ðSDÞ1£ñÉ"Ãû¢þýnó½/ah£ðÕ÷³|~âu9$aÀ(|<CÎÎ*üh~`	üpûÐ55^³êôåoñì!è¢ URSdX*iAAÚÅ-Õ:~} ÇUZQ®$¾x©ò6dâÙ ©ÕCx­¦úÅ ">·US ìó(èJA*
HüÖp(»h.ýN).ê ùÌÌ#m_à÷ðjjJ-Ð¾aÖØû'Æú
ø,Uæ17®¾b½ÔR#¢¸wß'$O$×ü¶KóH8gNq& û^H>æO
rR!¯e	FLÊôêïÙ½A0ZRÿ4dóB,GÌ.nULë<cD'äÃ7^¢¿iøÝû^Eè<D[8ÇØ÷x°_`ikHMÐðxÇKá¯À3~NEªe/ðüp]+y	ÞúÈ£¨$Ø')'#¿{É
´tWe~2/©l*Âmc ñEü9Ä§4BòË<IÜir±%a41V4aD`W÷È!Ý6 Ryk4 !\[UEG]UµXmôùY÷môx+×²){Nîªn¼GÍPtã½Gb%¹õ @@,ú´`r§ÙêÚ°  FÜ
¸¬géH\åSÌmåI\Øã¤OÅ¬f¨Éyä}µÌJÝvW>xÐgÕ>VNõ1É ²¾"z²U7Ñ×M5ºÄH¦
[Vú+øÄµ	º¯.ÄÐ\ddÂp"põëïË@;½"Sg¯·/	=ÏÅâ,EA÷¡ò<%SÉí=@|åì!pX.ÏÜòdy¯Æ²¯(² $KÒ[ãvÓh1DÍ¯0}Ûq±ÄT¥UqX|ÞØW$õR.q%DL./)Òñyö­îÀ;b¸ôOÑ
êqb.ëÁyÂÜe "RKëüXáÊò	¥uÎ ZàAZç ¦ØsfAÐ'.(°O&Ív^pÑâsòMý Îaê5¬ý(Iýàcée0<2k5 ÙÂÊ¢
8_¥`×G8rÏ 3Øæ©¹íXSË¿"vòqf0\ðYªbµÙ¹<¨SBRËtÊ­kê+ôa. Ç±åË+ùÌ6=è$p£h ÓªGBâÖf7="A3ÙºHvÚÛ|ö{ê³Ó¿øH6Gý
 TUN!yú÷ä÷åmõýY²æ´h Yl5·K¡%E­%#=f§)HçzÞi $í-«:E­ ë ¯$(ùgVÅJrVl&BdHÕ&h
$£ÚP§8YS@O"F³¢pÇ[X=$Ð3&± À`°õðlú®§^l_ë!AJ	D£¤"pÔ§,ªWÉ§eIPxÄ¥6^7L¼ÛÓÍEÜckIR8â²ªÃKföérW Þ5) .áIÊÀ«ðk	*[¬ Çkñ¤ùÂÈæ#÷ÑFvÂT}Ø-°G¢¥<*¥B¢ämd0õqáJÖ#0E£eÔÉÎïh5*#ë±Õw­ Öðâ)É2~mãáõ!º°$+È	·'ÖÂ«Z ë¢£¼ÿª>zmÚµÐø¨8GâÌ8³%Ä+þêõ$[ÒmÇÄ"i.ÌEÈÄ8ïÞ
!e±"$^Åm4hµÎ§q@zì+)H¯]U$_l^T/×««¤Ý 0ø¦×°ÑiMH£Xó»¦8çüG#7À¾Ú@Ò¾$M8^Û	©§°rú*M
Y	{@Ê{ZIÁIH'[ó¬eÅpëâã$ ÂY´ù×ÛgW¢(×¦	iùäB¯4%ÝF`G6NU¢Õº]z¤£IQ§	I
$²P¶Ã[/jt(s8¢è¡ôçsQ£È¾®;OÓ@ã¢¥@×9½ÄE¼:ØØ¤ý¡Q{k¨3¿HC?¯Ö}¯oRÌlÞ<øpad¤À'«Lµ¤+ ûmº>î©.(ª» å ñVS(@«mÅ¦¢.¼,kÌQa*J¾Ðì3á8ço±ê¡f'=½¯4fnÌ¦8p¸0(¤-2R¢ãt=ú¸F=ì«ô|÷ñýÄêÈÕùý¡Q ÿKyQRFÄö¬¿6¸Þ*¿d¯$ù¢ª>á(SS$²u1l@zñ3 ¸Jf¢x: T/N3G|tì½ícFh5nº`GÁ¼hí%qÉ]D#.¹g|o(GkZ!{©wR$:C»N¡À{UùÒòðJv\ïX[qDü¨íèïü5uXxüRÇ;·Î A`HÑP²i-¹]ÿJ¨v²¡Nd®Pßi©Ø¤«©$Ûp43U´!yO0Á*:¼û%#¥îÚgÉ "TC3É¦è=	x3ïJôzÈ~ëmØóû6F6&Dúd$a\ÒNQ|\V¸pÔ-Àt26FÊ¢sI&ùÈrÃIÁÖÊÔçÞ¹ÂÕÆEHóI*]ãå5¥yÈY^5[^u³èÌêÓ{Ì
àIâ* WÚ­H|a÷DO	õÝÖdãú\ºä<¾è#¢:ùE«"x¬ÙLëB·!Ã¼.Rê*÷»4#Âå¤Õ[B<î,z8HQUbn]"¾Å®DÈëÀö§$DQòn=¢'ãæ"º»*}NÙO2.G=Î©.;M#uËúÁ>ëbÍè¡åì·Úd¨>I(Ýþ¼(CB!¿l(ÃÄæ\EwßÕÇtôþTã9¾~çlqÌñ×ï	fluã{Ó(°VÝ(3úÐ wÕSI6ÙâÈc9~­.ySnªfd!Â1!¹ÙÝæ{ÃGl×¢{SÎóT¶jp^=ÕØÑTtd ßß
q6¥E$¶8:×!MÞ×] ÷ØÝ0#³á¢)ç(aÑl'zl!µÉp;
tµÌ­.ÒNÒ}¹çIÉ4>6Í"ë)R1êIEÎ¼\£WÇíØi8qÚ¹QÏi<Ne`Ä/F)T²F¦EJ[»¡@N¨KUðõ	¹!²Æ{k3uCÁ"æ0!8º7Û9s±º¹á^'¶pxùªô+
N´¡:¤ÍqÛÆ×BütUHãT,Kt	 Ð]CþÊxé.V«¦Z=3È;w9NØ$?ÑHÂéðDhÜ$6E¯Ò¯ä 1¼¨
éíÝ,!»Ú7x:¶Ue\D5Z0`£_¶qôºÒ+Us+fÕ,_¡§¹ 2tsÒd¹Êû«+b1éIËÝN2þÐOæ§Þ]Ò÷Ñãý&Äc+WQ`2 rÉqò¸­w3seÖçg+}Î¯Ò"á+¶C×Û£ôþ=T_`¤ïZË'·í{ótS*.ÚQ­%#3vM±¯È²+y/Èº%×É3ÜüÜÕF§
Scº:ænGV±sÁ ÈpB[X8­[R£çë`±n$åÉ*KÄ©
²DÛåHÜRÕ%zÄOq*ûÖÿip«N
ä'ØRáp#ìnJ¤ùæEÙ ÷
@ÉYãL¡×T©eRdRG¤OhcÃÁ";	Òú©æÒ¢®¡ÝÎBEsIGû.ZMÑÍìVÜ[5¥¤ôÛ6Å1zÈNä²åC\,ÄßaÜ6ÒÂ^Ø¦ÒZ$®ÜÕ{ÎZ5j8Ú><{ÈðëDÌO^G$²·rÑ6å
>JÄ#¿w®îXÜR¬Ùzdªß"1Æ¹·ÛæÕëuér	e9eTRsªÖðP¯Ô³3wEX=N¬ÌnÊ JHY5Åã@yd¤}eàËÞKÏ>.`/J$××¬Ûx·GÆ¤Uàx
¿òØ«,ã³áBüè¾rEÙ.^ËP¡.
XêV¤ä°ÕëãnMdFH'zß¿ß2´õ8lÅc.uÓ$#ºzÒdH»XqÚz!ëIÁ&öRY4½²®
[=¬YW9]ÕwreÑ]^ 2^×Õ âzOÒÚî¸8Ö%¢G9Z]"Ij²û¼êÓÝ>JDêz3µñ^°gZÓ-~ÿe©¡!~Óe´ìêEÕ.»_`Û"9¹­É4 S¡Ê~¤¬*ÁsiáÔÍ0GÜÚ2®-dzºè5!Å8w$°hÌJRäÜyN¹²a9d½NÒ !(ûM¤(3Ñ(BÄòIÉ²ÈX/Û
@:C/Ç> löÅ± À»ùÃúÃJÐù^Á'Wü§©WCmUøIt¨HS;ÄT!À]YÌÑÅzKEÚõÂ¾Xr÷¼_Ìª'Õ«ç»b=¢ÀÎÃs-zRqN.!
Ú{Köh5­´GÆÈº¨=³qÒ¬èè­$¬uhT+lîB(âP-¬}'½..fòvMó1wùÈ""çxqÝ
²$
K!|#ÂIlÈhpI!¦ ;
k2\Z@^'EGJ)ÓqGXðÚ:¨Ñ"1n2/NY¢õ¥»í8¶C=ß4ù iÎ$Çw4!M¨¶nkÄÞe0>ÀÊà-Y"$ìJ$Å`Â<Åû¦úy: [#>j ,X5£À¹ÂAß®*}¾{)Jäôí9èª5ä¤ÀÑ0¤ÚàP!ËÛ1e´'w+çô ãUi¶í ;>MF£.JDj=Þ	Ù8ÁZÖäÜ{ÞíÜ÷Ý\pY¯ ÒA8CÜ® Ò)Þð­2IH«c-(áÃGÂC-7qø½@B}·ÓØÇaÜ8v\ey§e}6#kä¦¤Ëó^¹Çnmq¿"F½`ÉÏÂ¡¾9Oj#'·]ÞWòD*;´@%KôÁ$4kE&	%`©áHôcÇVÝeËöÌIÍúÍäDî;£âF!Çý"¾KcX¥\ð1}þÈÁî( QáçP2^BôïÊ#l½BjÏ·<+6"urÃBÑñLHÓ£§eTÑ2G¸Ç
ÆÙ6ãI§ò]M|1¼¨d~³1,.fPr&Àÿv?ºßmÀ®õù {ØS><×¤O
\£ñXòñ´txüÙ"â:SQØþàË
àÎÚ×¾Xx½c¸?"ØS×Ä¦©Á|!û'Î4»»z\Åþó¥Üý°ópý±w#Aö¯;NË:tÇábñI¤è»gó^´[o3"ÞqKÄÖ[¦n4ry¿9ôvmy·ÁþÌb4¶!O"Ç[6þuNÆ¹c§Èd´G¶)Odjü?¨s%WáÇÃ
â`ýÂ6OÑ®ÎÁ¼SÁ(SV9ÜÂåq
OAÇønWÏåÉ'-¤¥P¡¶bHm1wåõÔIÀmøÅhìp Ish g8`xå p3Á¸bwa"äE2UÊµìÐ|<êBLvBÌor¸!æ°ßÁæbÛz¤CN¼S[
1<<)|wBv:Ã,ò¤¼wèPªñµbà¶s¸QX¹¥]uPwåpKÀ»à-bh0
7#ñ+1ÜrÝCÅãÉî]=j¸±íæÎÃí-®-øÌ'ÄàYð¶AùÕì¸´)>C²vnð¢aØ¸ìø®Õqo2PEÇ;0Q!þbxàöùÎT´¸½KcqÝäùLô³bÙÑÓ:ll½hjÑCkÕxF@k4Yª­§~ã½F.n¼JÂë`âÂ*¾ÓÔÏßqPhL'¤K¤$YöQd¶Y¶±XMg+nú,Cm­ÐúJÎlk¤Åã\jl $ÎÍ=o °#DÜ=Vdøe3dã3Üã¹ÄFÐw¦È^Æ\qñAæñrG[CÎhÂZïj$ÉKó^³ê<®ix²)TÔdfh
dH*^O9Iâg$ÁphEH«·v}ÀER\Æ:ddDEÏTÏÔVºÄc¹áûAXiÕ×Ðú3@4Åjª_ô[ò	Q(7ÇKQÙOFIDTOU"fÎ ÙöN6Oµ¤BPc×Iõé@'1#VÈ>Ä¢Qí_qãHÕ¤ÈrÀâö,7²RKDJßNv~¯àp²|Õwþ$ÃWVCâ,âlÝ<êTöÑ$Êp²%Dô Ò§¬MÏ¨ý¸«´h*"îÝjjÄí Uß-^mÂfQ%"¯2
íJ$·í5a±GõánÅ¡Ç½Â¸W Dt1aù2êÄUîÍg.ÅE­î¡ê~à°Ç©ûnÒ	±.Ú!ºnÚs5^A³¥kÄ+;2üÓnö"dÇmº.¸¡9d©ÃÓ/]U(x»]GæI*#!Kñî°\Ú_ÛÍaÛ¹R+è>ÂET]bE'&ÆÁ$¹'mÊ|`éëÑ&ÀmFóÆ æå&Í#ö½<%GæÀn ¼AH«¥LK¬ïÝ	.	»o+	îùÀëW	[Ä
½Uð(lµïç÷Y+ð-\ ø4d¸ ®ð+ÞpQ$>âôâ.Ä	£X e8iAÕ«[PT4lì'(
ueÑ*\±¦¥«?Y*ém³¨#dvyQN4Ùhäç®(Ð	ì,)
<¥ÈuÑÜàq¾§:«5vBhó(fÑXqGàù/b«÷áöd.ÝNëpï¡¡´ï$ªO/3}:¿¶«SU¤îëdd¾@+²ñÛ
²ÔÈrBÞ¬ZQ<ÇTÅÄZ'ZXNFQGrlZa±Y]ò°3)Ü¤
´¨ØJÞ4|yròÜ×V¢"Üy§·ýbüñ¨Íæý HÁ.ÕBÞ@Yp­ü( ^"Y¬Ñ 0n»qDQä Ñd¢zd(Þ4Þ!Ô¯2I`B8³'D2ªgìª@ßêqóËL,'û ±<¯©NÈr¹lyÚÏ»W° ûQ5âåmvE÷¡ââ\Z4¢¯åZb<ãzT~4JNÑÃ§¼äTmÓÃWx!y[ÅíWm}ÊN¨Ë¬65vY(	%`hìYLòr¾©±KzÕÐ[daÇìzãð¸±tád<§dÙ!7\dµË®ñLjµõf,òúxU5ô§«u§­gñDÚ´hBr¢cñ´ì
4¸Æ¨É%ËÀ ¬M=v*ê®ï4ÓcÙÖX©Ó¼¨±"Np£]#sï«Áúc¨îv!¾·9«6øBÄ"7µäNTn¼dâ¨ßÅbñ§´ Ì3víq
ñía6çÐzøqÜ<<x·j#CZ%.Ã(
@ìÓ÷_`£Ýs¹Ôª«òdØáIÿPN_tUTyî~·´H^cÓ..K%w8ô4©3tN¯rlã¸ÕbÔs7ë$ð¡,%a ë5»í>yì7ô£<}G8oïh¨8`%æ<Yô½Â¯9{àR_%7KlEòp4FLîüHCÂ/%Y$ÒàÔîBÖ«`'.TEØ¸¨4È
ÒSv2¤`qº6j:yên×HZÄmÙ:%ÕR·¤%Þ¡(B}BäYâõ/	±cB|¸ø\È§	=bAoj$£ÃuvÂ®D°G½`
×ªG­É=Ú')H\iÒ)cW~;ë¸{ûN5¥â,iªqÒå8õY7Ù+½7ÐA]¹²'$´ )"DäÑ4-IH`Äm1Ùè8OÐ0eñ¤­,ª`¡$×Ö$²Ul²0¼75Û6bÿ4!Ïä¶È®Tr#g#´ëßûE½ÂHSC\(eÃ®¸&$±nkÜ2¡Ç¾³wzæ¸m+qKÄdÖÏý°%·¸y
¶<Ñ»PÄ0Æñ£Å.eqÏAÃ;Y7H9a*.;â´+AnKU@'VÙ!õ*.ØJtxniSÅ."w(BqP"¼ÛÀTÈ8ÚôÖÀsºvËâ
"¥6Äó¥Ó!½÷a ¡1$#-@îY¤] Û
ØãsFx·t
Z<¸Ç®ªJçÐfª¶
Î¯(ªònéÓLF¯yÊfÄ³¶[¦Âý"W¨úPùÅDÔ%âVq;û>´E*£f³¿¼ÙÂÏ²¶çÙ¿ë PåE­¬êyÐ±~7y¬Yº}¹àF+=	ä UMÁA=â(©É´%¹©K¼ÛÀG¹òm5¤p{Aá.u.Þ¬4"a-ºVKyHÌpþ !Wq_yrw2¿yMH5×xë8´­lJÁÎq\CÜüVUáb!.£îªªxÇ´ðòÈçË;s¢DàÎxt¶)ÑeL¯7¥ÁÚIn<eÊÞhF:
¾¯¯ Þ0·±R[,°Z®ìpõ'êqìAÁv]ÿJB	X*æDþh¼ýÐ¹"'­±ÈÃÔÐ»Í²kÇ÷Çó;BÓcFàWäjßÎD«*§âýäùv2ÆG]3Äs,)ÌEuåám-*'WÖò5	ÀÙgud¼ëUSvÑj-Çê» E6*°4À½3Wn¼pf¥ðf£tävé}p1u  êìñ&¦ ©$NªwÐe=ôNz:¹syãd~NöhÅãXÉô%ÔìlòñQÕèòø¨Ñð ¢N£Úq¬©EtÌN®YEÉíC§
â½ãU¤$+þÈ»!#··AHewdHeX5ºÀk
Q É]RE¦nKD
¹¡H°møð	}!¶ëøvÂÑøË¢Ædji_c"ç¤Î3Û XÖÀ¹A®*Ó´ÖBVUìNH
ÓI«¡}¶HÃ^ÏçccÏ·%!Ic?iÝgYEè°Sl;Bx^#©½)Ø½ßzÜc2ºó|_ÑïTdtFÞÝH^YñÓ"ê ¤MSçZ¸.GÆX¾5\	fÂ» Qd¼¶ÓÒ¹W U:$l'¶ÀÝÆnuÒ4 h@Gl@xÃÑÎÈRf)ÉC<
Uýd?m;WBäÃón&¼Û	Q«ª`á
îTßd#"8¯&$ã ¸Üâ­¸ñ­q9ÇíJÇM²ÝHfnRsÜB¶pÐ¼-CC)Bh&ÜUyáûå-¦ÜÊzÀ)örOôV9,ÞYßqÃÛÚ÷¦É2Q­­ØÊî*ÉÁð2rÆ5Ù\;g8ÜgªI çñC0.à ¢³$Ì3o¹èV2Üâ*
^3 I%¼Ñ-±ÖVQbè>G8ìÂ7ADSî6ð^b@ ¶NÀÝÐª¢)Cbïm#ÏÑ¸í¶ÞKDLÀn£´ª²á4Û6"Yov+B<OpL)ðÉø×&âÐ*á£DúêËÄ÷
PÌ¬.û1bÅ.öOÄÃé?A¸å^tsºÄ¿ë§Ã.,<B²¶9ßÚndI}Àzðw3û
ÞK¤nãÖ¦}
hUãìöû\³WAoí­'QÁ{Æ°5º
ÐªÂ>wqz8ZáÙ­¡L½f´FÈé^~Z@Û	ðQ"rºF§8tÜæPÄv&É¸±d/¦ÿ«ÓR×ed®f;É¨àÝ¸`ÝêpëG
®1ä5æÀ0<ñ[G6ÓQsø&]rÙk¡À±ªÑâ±^i«ªp iò-`Vq ÙbìülÚ8¬A(K³G¥èðJ.ïV¶­`Ø­î_	[ùâ{B^
7}<¿ã½¬6åýù~ þüL/B^|·ÈÃÎZTTLîA4¹¸Ûºm@gÎp¯ ýÏÕöù¼ï
qû¤¦"²YëZk­R§ÎùêÙ&),.<<n0§âÉÉF2/dFbÿcyú3P1vìàfÌ5
¤SÆ&ãà Aóç3ÀD2eÐLÐÊXqïÈq,QãÝ»ùyEqåv¢ éÈ´g¢ÔñÇY£ñ \]ç(ê×rÎDÝ¼ÆåMDºV{m4ü£Â+ 
À§àÂM¯N@>b:ÉúÆÝ­^$u½b¨S>ìA÷ ®pl¬Ç9¹êÎ30 w;@+5è¶ òoE°¦§µ]Ï0ÚÈu¤28ÐejGåt ±ËËtæ¬HÃÛùUï]×OZJõ2[t.)`c@ááÎbâ4è~¸=qó}Q÷m8
hý³Øy2zâ&ÀK&6­7nO"C{x8«s|vYôÄa²ùDuì8¤ZLNªvu1q
t¿àgÈ^OÇkÔ0ÍÏ)\p¸7Zð;=
Ðºp<¼nò ×I7ÊµâÐ©vz!ÑA ¸;
|ä
¥ÅÀ#G*i5
tJ~ÀÎ Äv3ØB®ÆéaÄyfFã¡g!Ã½Fh#@«)ã~Æ	p#NÐ½(%Ü ìK&
ïq.ÆRÓ;ÉÈD'<^³âÄMò³Æ
³oêO§ñêÏ¾ÔnÒ­Â+ 
À§àqaÂFR ºiû)ª¢""lº¸jU2­[S Z\¸Û:EFú6«z "Þùò#Ûþ5ý@öI3J¨µ;bc#Á¤õºÖndÈíÍh=è~¸¥¼Üd\´rN "kÚ`_0ù77!èùºXtgèUÕ
ºß nwg¨\ÑhÓúÜwçDXt§k(Ù®mtwè~¸¥¡Ü¤¡TÝ·z]

!)º´×ÝP S¶)À-
á&
!ë~0õïÛ
,uËä¡F²x#ò^÷>w~)ÞOb[ôy
â^I²"£^ð·="¯ñµjh½P¢7ô¥D÷Do^i > ¤¢ÇV@pøDãª dõEÄºRu·¼Â-aË Üí "**=pEÍ±Ï°XäënM>DÆòDêZE%Mµ>.ûtâmïy¹ðaºß nÉ5?iÙ¶.¦;OÎrÍObZtwi1CäªîhÐýp»;CLî1½éÎÓ²;]LËît«»£@÷À-1í'i¬ê¾
º]LK.ru7èX
pKLûIËYÁ½JDÛ³Âï¼!Ã£îMÀÞIÊÞ¬AÜ8	B ¥}S ¯%y+è¸¡/pp«ð
Hðét,] K@À².E³
p·4.gåÀ
 "$x2Kô×Û#`%°µ±	`N­k¸`àétÜ§á:7
aÚö û
àÄ
 më¸óä,±Â$Ew 9©êÝo ·»3°èÎÀî,»Ó°ìN¦º;
t¿ÜÀa³ªîÛ »Ñ°¤èÂTwCN¥ ·p°¬ûÀÔ¿o7½&1})>N²Xôy
úc;èÔ3]]ÎÃ~£KÔ²9Ú#ò«T¤Ã@´9wÂÿ8]y\C£Ï¨E¶*à8d±èüìa OG>(GÄÀ<^Í7Ø«'nÿxòØß#¶éû95ðÛ4psQNÉ¸ç¨ÊI&àâóøþÐ¨Pø9´±yZðÎ´:³C|Í¨ÚRÁ%ìñý,	WU]·©]#³4B5 ³ä-X>_	}'tOrõ àEmAÄ>ðc½D¶{§Û Ù`§'v]b¸¸§pYOUUÅ½FU!W;a¥.ñ¾½$D.¼¨ÛÁA@Ièodçî6ð^"Ò¦àígÜ+@Æfy>Í ÃÃs 6Aä~^Ôø0Æê~8e7r)¨ÎÓÞÉ,ÂJ5rwx/¾Æämí#3 md×a
ÄydZs#©aÉ"¸Â Å}¢¨YËa¿¶LóuÚN  iSw;@+ßpZ9r}]ð²ìsU0ÅïI¾â1bT$qVvs¥à½Dës%ªÎi@«
eÍ4»²¶àÕ¬Ë§ÅJô8éÀµzETH1oÅÞ³â¶¨&ìzé®£ùDxÓ¤SÀËóÍ©<ÏJÍ5Ýt%Üt	^ï<¾w=#²7ªï¥"g¶úk*×æêKáÓªÅSù½D»Õ4üþÎlýTcÖÇ%!BgÝ-"ùxê¯î+éÖadÄ¡®f2Mß®»kÜÓXZ1ß·¨>H¾#ï»çxØ×ùÆo´;àø3'?ø³ÃBæöVÍä­ùDC¨Î^ÂxxÉðl
¨¸ñ+/kôsB/Û4óÍ"M·åj´öl_ääûp´Ï÷ËÏïy©Øø\T¨ø¼ÇìMÀ·±N³÷éèEp±{¡»^7xYy`zp!ÿ¸ÜÝ7¾Å3^!WÄÔöã3®ã#_(ïðX}]Vå+`YÕB«;¢pY,oT_É­ëÑeªK²á¸CåM~§!³Èy'«ñ$)!nd5ÌÛ|ÍÞöjóÜcÍàÂ:"gTQ.".Vë!.H¥h4a&ÓÙ7xIÐhôXÅ±h¾¢bó¥EÐ"Û
MØµª^0	_ IPýf@fîtUOBRv[UBbèUç2Û©¸¹s',9â>nJE¸XuäËTeSúp!,C@YIè7Vy×RK©""Æ²á¶-V"ì¶h­¡Òúô;Ñúé@Çqp¢Êb¯â1ûøþéøN#¼°þò¼¨r
vÀdïîé{-Ç!Ü.ù:U}·È)iqÁAÄaOÿ%G,·h5Æ½bB$îåòæ`ÖEåäj¦@#¼OT bì)Ë¢*öû¢Ó5a=ê!´µØæ¹;)¿+d¿à]+â{i =¶àé`ÁjhXYñhÅ-¦ûÃkj-"\°LiípA<ÏÀüòçiÔÒ´¯uÁZ7YjF*\·«:È¼îí:óa'¿müv	|èQÁ¢qA*G¯¡¼`F<Ùn)Í}Í vô£;¶ÐëhªÉÝUìWéåm<Ðo·qW¥Ü­ñ~yÛ´é¥LÚÑwãÛ ×ñýtàE¸·yà5¤4ð«?ñ,õãdº Ê)óC_Çü:_kNdXJ|<«x'Ù6 ÕÓ§ð	ý:"ÙîÄum]ri¶÷.°é~ÚL÷6ç>Øÿ=«-#5iHc©c¤²S6¶1Æ¸pÓ
¼d'l	Ò%=M JxÜ°Ö³s¿¨Õüü4þ¿N¢Ë[ÂÓ 
#­»>ã2+ðrÊ0nÕ"tìy¤V-é¼7ü}Z8ÐxËI6¥¦@T7²á5d.ÇlxN°ÜdµÝÚPà¢½k¥øóÝ4ý,ñæÄ½ÀÅTtÃà!CÌVw~ÿvó}/ÍªWÒJ-ËzI¥¼³^¨|¯¶T¾zÛ¸YØ^/hÙ4Ý3Eôid
½tq	ÞÜ¾¥W¼ãç¬¯¹}@ÜÞØLyÛñZiÅí×ï;|Éí'T1ÄJòTÞä)<ÏÎñå¦æ¾15;ôS^Vú¢­¸}Ý!WehQy`if¿¾¥ÛÖ¢}3H?¯xÉìmÒ×Ì> fo\¦¿´sÓ~Pcì@©ÖÍ|8OMþøêWdSÿ^ýîï?üôñ]ïÞþ×ï?üøöß½ýí÷_ÿô·÷ÿøãÛ/>|ÿ÷_Ó¿ÿÃïþôÅüû¿ýëoþóí÷?|ÿáãÛ¯¾ùøþÃo¾{ü>üõ«¯ßäí¯Oð¾ùá×_ýôéÕ¯¬¯¾øþÇ_ýíÏ}ÿîÕo¾ÿéïüÁW¿rAÿòã×¾ûáãO{÷êÆ¥5«E 
```

# report/after-figureplacement/report_fp.tex

```
% ======================================================================
%  Personal LaTeX Template
%  Victor Gabriele De Matteis
%
%  Usage:
%    1. Copy this file and rename it for your project.
%    2. Fill in the title block below.
%    3. Replace \bibliography{references} with your own .bib file.
%    4. Compile: pdflatex -> bibtex -> pdflatex x2
%
%  Citation style: numeric [1], compact (abbrvnat).
%  To switch to author-year, replace:
%    \usepackage[numbers,square,sort&compress]{natbib}
%    \bibliographystyle{abbrvnat}
%  with:
%    \usepackage[round,authoryear]{natbib}
%    \bibliographystyle{apalike}
% ======================================================================

\documentclass[10pt,a4paper]{scrartcl}

% --------------------------------------------------
% Encoding and language
% --------------------------------------------------
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[english]{babel}

% --------------------------------------------------
% Page layout
% --------------------------------------------------
\usepackage[
top=10mm,
bottom=10mm,
left=20mm,
right=20mm,
includefoot
]{geometry}
\usepackage{multicol}

\setlength{\columnsep}{8mm}      % distance between columns
\setlength{\columnseprule}{0pt}  % no vertical rule

% --------------------------------------------------
% Mathematics and graphics
% --------------------------------------------------
\usepackage{graphicx}
\usepackage{subcaption}
\usepackage{wrapfig}
\usepackage{mathtools} % loads amsmath automatically
%\usepackage{amssymb}
\usepackage{bm}

% Let figures pack densely with text instead of drifting to float pages.
\renewcommand{\topfraction}{0.92}
\renewcommand{\bottomfraction}{0.85}
\renewcommand{\textfraction}{0.06}
\renewcommand{\floatpagefraction}{0.75}
\setcounter{topnumber}{3}
\setcounter{totalnumber}{4}

% --------------------------------------------------
% Tables
% --------------------------------------------------
\usepackage{booktabs}
\usepackage{tabularx}
\usepackage{array}
\usepackage{ragged2e}

% A ragged-right fixed-width column for use in tabularx:
%   \begin{tabularx}{\textwidth}{@{}L{3cm} X@{}}
\newcolumntype{L}[1]{>{\RaggedRight\arraybackslash}p{#1}}

% --------------------------------------------------
% Typography
% --------------------------------------------------
\usepackage[scaled=0.95]{helvet}
\usepackage[expansion=false]{microtype}   % expansion=false avoids font errors
\renewcommand{\familydefault}{\sfdefault} % sans-serif throughout

% --------------------------------------------------
% Lists
% --------------------------------------------------
\usepackage{enumitem}
\setlist[itemize]{
leftmargin=*,
topsep=2pt,
parsep=0pt,
partopsep=0pt,
itemsep=3pt
}
\setlist[description]{
leftmargin=0pt,
topsep=2pt,
parsep=0pt,
itemsep=3pt
}

% --------------------------------------------------
% Colour and hyperlinks
% --------------------------------------------------
\usepackage[dvipsnames]{xcolor}
\usepackage[hidelinks]{hyperref}

\definecolor{accent}{RGB}{25,55,140}   % change this to restyle everything

% \emph produces upright accent-coloured text instead of italic.
% Use \textit{} when you genuinely want italics.
\renewcommand{\emph}[1]{\textit{#1}}

% --------------------------------------------------
% Citations
% --------------------------------------------------
\usepackage[numbers,square,sort&compress]{natbib}

% --------------------------------------------------
% Header and footer
% --------------------------------------------------
\usepackage[automark,headsepline=false,footsepline=false]{scrlayer-scrpage}
\clearpairofpagestyles

%  Edit the three footer fields below:
\ifoot{Victor Gabriele De Matteis}          % left:   your name / student ID
\cfoot{\pagemark}            % centre: page number (always keep this)
\ofoot{AC Course Project 2025/2026}    % right:  course name or document label

\setlength{\footskip}{15mm}  % pushes footer toward the page edge
\pagestyle{scrheadings}

% --------------------------------------------------
% Spacing
% --------------------------------------------------
%\setlength{\parindent}{0pt}
%\setlength{\parskip}{0.45em}

% --------------------------------------------------
% Section formatting
% --------------------------------------------------
\RedeclareSectionCommand[
beforeskip=1.4ex plus 0.2ex minus 0.1ex,
afterskip=0.6ex
]{section}
\RedeclareSectionCommand[
beforeskip=1.0ex plus 0.2ex minus 0.1ex,
afterskip=0.4ex
]{subsection}
\RedeclareSectionCommand[
beforeskip=0.7ex plus 0.2ex minus 0.1ex,
afterskip=0.3ex
]{subsubsection}

\setkomafont{section}{\large\bfseries\color{accent}}
\setkomafont{subsection}{\normalsize\bfseries\color{accent}}
\setkomafont{subsubsection}{\normalsize\bfseries}

% A section command that adds a faint rule beneath the title.
% Use \RuleSection{Title} instead of \section{Title} where you want it.
\newcommand{\RuleSection}[1]{%
\section{#1}%
\vspace{0.2em}{\color{accent!35}\hrule height 0.4pt}\vspace{0.3em}%
}

% ======================================================================
\begin{document}
% ======================================================================

% --------------------------------------------------
% Title block
% --------------------------------------------------
\begin{center}
{\Large\bfseries\color{accent} Automatic Control Course Project}\\[0.25cm]
{\small Università di Trento --- Automatic Control --- A.Y.\ 2025/2026 $\qquad\qquad$}
{\small Victor Gabriele De Matteis, ID. 265757}\\[0.2cm]
\end{center}


% ======================================================================
%  BODY
% ======================================================================


\RuleSection{Introduction}
\label{sec:intro}

All the constructions in this report are carried out on a linear model obtained by linearizing the dynamics about an equilibrium $(\bar x, \bar u)$, so all coordinates are expressed in deviation variables 

$$\delta x = x - \bar x, \qquad \delta u = u - \bar u, \qquad \delta\omega = \omega - \bar\omega,$$

therefore, every gain, estimate and response that follows is understood relative to this operating point. 
Absolute quantities are recovered as  $x = \bar x + \delta x$ and $u = \bar u + \delta u$.

\RuleSection{Modelling and Linearisation}
\label{sec:q1}

\subsection{Operating point}

The steady state conditions $\dot h_1 = \dot h_2=0$ provide two algebraic equations in five unknowns, so the admissible equilibria form a three-parameter family. The operating point therefore is a design choice rather than a property fixed by the plant. 

We linearize about $\bar h_1 = \bar h_2 = 0.2$ m, the maximum admissible level and therefore the most demanding steady state the closed loop is required to hold. 
Fixing the two levels leaves us with one parameter: $u_2$, upon which both $u_1$ and $u_3$ depend:

\begin{equation}
    \bar u_1=u_2\,\frac{a_2\varphi_1}{q_0},\qquad
    \bar u_3=u_2\,\frac{a_2\varphi_1}{a_3\varphi_2},\qquad
    \varphi_i=\sqrt{2g\,\bar h_i}.
    \label{eq:q1_balance}
\end{equation}

\begin{wrapfigure}{r}{6.4cm}
    \centering
    \includegraphics[width=6cm]{figures/find_u2}
    \caption{\small Actuator margins $\min(\bar u_i,\,1-\bar u_i)$ against the free parameter $u_2$, with $\bar u_1,\bar u_3$ set by the steady-state balance~\eqref{eq:q1_balance}. }
    \label{fig:q1_margins}
\end{wrapfigure}

We resolve this last freedom from the margin sweep of Figure \ref{fig:q1_margins}. Controller margins will become especially important in the next sections, where all inputs will be bounded also during the transient dynamics. The most balanced choice is considered to be the point where these curves cross, it's the best compromise among controller margins. This occurs at $u_2 \approx 0.32$, which we adopt. 

For the operating point, the leakages are set to their nominal healthy-plant value $\bm{\bar \omega} = 0$  for two reasons. First: with the prescribed values for leakages, the equilibrium $\bar h_1 = 0.2$ is not feasible given capping of controllers, as the following calculation demonstrates: 

$$\bar u_1 = \frac{[u_2 a_2 + \alpha_1(h_1 = 0.2)]\,\varphi_1}{q_0} = 1.06 > 1.$$

Furthermore, $\omega$ leaks are framed as disturbances, not as operating conditions, therefore it's not logical to find an equilibrium considering any different value for $\bm{\bar \omega}$, even if mathematically possible. The resulting equilibrium is the following: 
\begin{equation}
    \bar{\bm{x}}=\begin{bmatrix}0.200\\0.200\end{bmatrix}\!\text{m},\qquad
    \bar{\bm{u}}=\begin{bmatrix}0.660\\0.320\\0.326\end{bmatrix}.
    \label{eq:q1_eqpoint}
\end{equation}

\subsection{Linearized model}

Collecting the two balances, the plant follows the following ODEs: 

\begin{equation}
    \dot h_1 = \frac{1}{A_1}\big[\,u_1 q_0 - u_2 a_2\,\varphi_1
                 - \alpha_1\varphi_1\,\big],
    \qquad
    \dot h_2 = \frac{1}{A_2}\big[\,u_2 a_2\,\varphi_1 - u_3 a_3\,\varphi_2
                 - \alpha_2\varphi_2\,\big],
    \label{eq:q1_nl}
\end{equation}

with outflow velocity $\varphi_i \equiv \varphi(h_i) = \sqrt{2 \ g \ h_i}$. Writing $\bm f = [\dot h_1, \dot h_2]^\top$ as a function of $(\bm x, \bm u, \bm \omega)$, the linearized model is obtained from its Jacobians at the operating point, evaluated at our chosen equilibrium:

\begin{equation}
    A=\left.\frac{\partial\bm{f}}{\partial\bm{x}}\right|_\star,\quad
    B=\left.\frac{\partial\bm{f}}{\partial\bm{u}}\right|_\star,\quad
    E=\left.\frac{\partial\bm{f}}{\partial\bm{\omega}}\right|_\star,
    \label{eq:q1_jac}
\end{equation}

The only nonlinearity is $\varphi_i$, whose derivative at the operating level is: 


\begin{equation}
    \beta_i:=\frac{\mathrm{d}\varphi_i}{\mathrm{d}h_i}
    =\frac{g}{\varphi_i}=\frac{g}{\Big({2g\,\bar h_i}\Big)^{1/2}} .
    \label{eq:q1_beta}
\end{equation}

Since the inputs and disturbances enter \eqref{eq:q1_nl} linearly, the entries of $B$ and $E$ are read directly from the coefficients of $\bm u$ and $\bm \omega$, whereas $A$ alone carries the derivative $\beta_i$ through the chain rule. For $\dot h_1$ and analogously for $\dot h_2$, the following holds: 

\begin{equation}
    \frac{\partial\dot h_1}{\partial h_1}
    =-\frac{(u_2a_2+\alpha_1)\,\beta_1}{A_1}
    \;\xrightarrow{\ \bar{\bm{\omega}}=0\ }\;
    -\frac{u_2a_2\,\beta_1}{A_1},\qquad
    \frac{\partial\dot h_1}{\partial u_1}=\frac{q_0}{A_1},\qquad
    \frac{\partial\dot h_1}{\partial\alpha_1}=-\frac{\varphi_1}{A_1},
    \label{eq:q1_entries}
\end{equation}

Evaluating each matrix at the operating point yields the following matrices: 

\begin{equation}
    A=\begin{bmatrix}
        -\dfrac{u_2a_2\beta_1}{A_1} & 0\\[.2cm]
        \phantom{-}\dfrac{u_2a_2\beta_1}{A_2} & -\dfrac{u_3a_3\beta_2}{A_2}
      \end{bmatrix},\qquad
    B=\begin{bmatrix}
        \dfrac{q_0}{A_1} & -\dfrac{a_2\varphi_1}{A_1} & 0\\[.2cm]
        0 & \dfrac{a_2\varphi_1}{A_2} & -\dfrac{a_3\varphi_2}{A_2}
      \end{bmatrix},
    \label{eq:q1_A}
\end{equation}
\begin{equation}
    E=\begin{bmatrix}
        -\dfrac{\varphi_1}{A_1} & 0\\[.2cm]
        0 & -\dfrac{\varphi_2}{A_2}
      \end{bmatrix},\qquad
    C=I_2,\qquad D=0_{2\times3},\qquad F=0_{2\times2}.
    \label{eq:q1_B}
\end{equation}

The numerical values of $A$, $B$, $C$, $D$, $E$ and $F$ are reported in the accompanying code.

Three features are worth noting. First $A$ it's lower-triangular, thus tank 1 does not depend on $h_2$, reflecting the mono-directional cascade in which tank 1 feeds 2 but not conversely. Second, the disturbance-channel gain $E$ is numerically large relative to $A$ and $B$: this offset is driven by the magnitude of leakages $\alpha_i = \mathcal O(10^{-5})$, so the effective disturbance $E \ \delta \bm \omega $ remains small. 
Third, $C=I_2$, and $D=F=0$ because the performance output is the state itself. The matrix $C$ will change when considering one faulty sensor. 

\subsection{Open-loop dynamics}

As the matrix $A$ is lower-triangular, its eigenvalues are its diagonal entries  $\lambda_1=-\frac{u_2a_2\beta_1}{A_1},\lambda_2=-\frac{u_3a_3\beta_2}{A_2},$
which are both real and negative: $\lambda_{1,2}\approx-2.02\times10^{-3}$. This means that the operating point is an asymptotically stable node. 

The two eigenvalues coincide for a structural reason: substituting the balance equation $\bar u_3a_3=\bar u_2a_2\varphi_1/\varphi_2$ from~\eqref{eq:q1_balance}
into $\lambda_2$, and substituting $A_1, A_2 \to A; \quad h_1, h_2 \to h; \quad \varphi_1, \varphi_2 \to \varphi$ 
since they're equal, we get the expression 

\begin{equation}
    \lambda_1=\lambda_2=-\frac{\bar u_2\,a_2\,\beta}{A},
    \label{eq:q1_sameEW}
\end{equation}

which is the same for both $\lambda_1, \lambda_2$. This implies that the eigenspace is one-dimensional: $A$ is defective. This does not affect stability but determines the mild sigmoidal shape in the rise of $h_2$. 

Finally, the eigenvalue magnitude  corresponds to a dominant time constant of $1/|\lambda| \approx 495$ s: this means that the plant is stable but structurally slow, a limitation that bounds the improvement that any input-limited controller could achieve.

\subsection{Linear versus nonlinear response}

Figure \ref{fig:q1_cmp} compares the linearized and nonlinear dynamics under the constant equilibrium input, from two initial conditions: a small deviation close to operating point and a large one, far from the operating point. Both models converge to the equilibrium, although for a large deviation the separation among the two trajectories is significant. 


\begin{figure}[tb]
    \centering
    \includegraphics[width=16cm]{figures/q1_cmp}
    \caption{Free response of the linearised and nonlinear 
    models under the constant equilibrium input, for different initial deviations. Far from the operating point the
    nonlinear response fills faster, as the linearisation overestimates the
    $\sqrt{h}$ outflow.}
    \label{fig:q1_cmp}
\end{figure}


The mechanism is more clear for $h_1$: the outflow velocity $\varphi(h) = \sqrt{2gh}$ is concave, so its tangent at the operating level lies above the true curve. The linearization therefore overestimates the outflow at low levels, since it enters negatively in $\dot h_1$. This overestimation makes the dynamics of the linear model lag behind. For $h_2$ it's more complicated, since the nonlinear term enters both positively and negatively, specifically, $\varphi_1$ enters positively and behaves as a inflow component, while $\varphi_2$ enters negatively and behaves as an outflow component. What matters is the impact that those non linear terms have in the overall linearization error, so we decompose their contribution to the overall linearization error. We start by restating the linear and non linear dynamics of $h_2$:

\begin{equation}
\dot h_2^{\text{lin}}
=\frac{1}{A_2}\Big[\,u_2 a_2\,\varphi_{\text{lin}}(h_1)
                   -u_3 a_3\,\varphi_{\text{lin}}(h_2)\,\Big],
\qquad
\dot h_2^{\text{nl}}
=\frac{1}{A_2}\Big[\,u_2 a_2\,\varphi_{\text{true}}(h_1)
                   -u_3 a_3\,\varphi_{\text{true}}(h_2)\,\Big],
\end{equation}
with $\varphi_{\text{true}}(h_i)=\sqrt{2g\,h_i}$. Subtracting term by term gives the two-channel error

\[
\begin{aligned}
\dot h_2^{\text{lin}} - \dot h_2^{\text{nl}}
&= \underbrace{\frac{u_2 a_2}{A_2}\,\big[\varphi_{\text{lin}}(h_1)-\varphi_{\text{true}}(h_1)\big]}_{\text{inflow-channel error, }\;e_{\text{in}}}
 \;-\; \underbrace{\frac{u_3 a_3}{A_2}\,\big[\varphi_{\text{lin}}(h_2)-\varphi_{\text{true}}(h_2)\big]}_{\text{outflow-channel error, }\;e_{\text{out}}} .
\end{aligned}
\]

Figure~\ref{fig:q1_decomp} shows this decomposition: the net error $\dot h_2^{\text{lin}}-\dot h_2^{\text{nl}}$ is an order of magnitude smaller than either errors, and it is negative. A negative net error means that the linear model fills $h_2$ more slowly than the nonlinear one, which is why the linear response of $h_2$ lags as well.




\begin{figure}[tb]
    \centering
    \begin{subfigure}[t]{0.49\linewidth}
        \centering
        \includegraphics[width=5cm]{figures/q1_decomp}
        \caption{Error of $\dot h_2$ split into inflow ($+u_2a_2\varphi_1$)
        and outflow ($-u_3a_3\varphi_2$) channels.}
        \label{fig:q1_decomp}
    \end{subfigure}
    \hfill
    \begin{subfigure}[t]{0.49\linewidth}
        \centering
        \includegraphics[width=6cm]{figures/q1_err}
        \caption{Norm of the mismatch
        $\|\bm{x}_{\mathrm{nl}}-\bm{x}_{\mathrm{lin}}\|$ for the two initial
        deviations (separate axes).}
        \label{fig:q1_err}
    \end{subfigure}
    \caption{Linearisation error. Channel decomposition~(a) and overall
    norm~(b); the peak scales by roughly two orders of magnitude between the
    small and large deviations.}
    \label{fig:q1_error}
\end{figure}

Figure \ref{fig:q1_err} quantifies the overall linearization error through $\| \bm x_{\text{nl}} - \bm x_{\text{lin}} \|$: the peak error goes from  $\approx2.5\times10^{-4}$~m for the small deviation to $\approx2\times 10^{-2}$~m for the large one -- roughly two orders of magnitude --
confirming that the linearization is accurate near $0.2$~m and degrades progressively as the state moves away from it.

\RuleSection{State-Feedback Design via LMIs}
\label{sec:q2}

With the plant linearized, we design a static state feedback $\delta \bm u = K \ \delta \bm x$ that makes the fluid level go back to the operating point. When applied to the plant, the command is: $\bm u = \bm{\bar u} + K(\bm{x} - \bm{\bar x})$. The feed-forward term $\bm{\bar u}$ sustains the equilibrium while $ K(\bm{x} - \bm{\bar x}) $ supplies the correction. The three gains obtained below rest on a common construction which is presented in this first section.

\subsection{From Lyapunov inequality to LMI}

The pair $(A,B)$ is controllable as $\operatorname{rank}[B \ AB] = 2 = n$, therefore its closed-loop eigenvalues can be assigned arbitrarily by state feedback. It's worth noting that the open loop is already asymptotically stable: state-feedback here it's not implemented to secure stability but to overcome the structurally slow response.  


In deviation coordinates the closed loop is $\quad \mathsf{d/dt} [\delta\bm{x}]
=(A+BK)\delta\bm{x}$, which is asymptotically stable if and only if $A+BK$
is Hurwitz, equivalently if there exists $P=P^\top\succ0$ such that
\begin{equation}
    \operatorname{He}\!\big(P(A+BK)\big)\prec0,\qquad
    \operatorname{He}(M):=M+M^\top.
    \label{eq:q2_lyap}
\end{equation}


Condition \eqref{eq:q2_lyap} is non linear because of the product $PBK$. Following product $PBK$, and is therefore not a linear matrix inequality. Following \citep{zaccarian_lmi}, we can get to a linear LMI by the change of variables $W := P^{-1} \succ 0$: pre and post multiplying \eqref{eq:q2_lyap} by $W$ and setting $X := KW$ gives 


\begin{equation}
    W \succ 0, \qquad \operatorname{He}(AW + BX) \prec 0
    \label{eq:q2_baselmi}
\end{equation}

now linear in $(W,X)$. Any feasible pair recovers the controller and its
Lyapunov certificate as
\begin{equation}
    K=XW^{-1},\qquad P=W^{-1}.
    \label{eq:q2_recover}
\end{equation}

Controllability guarantees that \ref{eq:q2_lyap} has at least one solution, so a stabilizing $K$ always exists. 
%It's worth noting that there are no further requirements for the response 


\subsection{Arbitrary stabilizer}
\begin{wrapfigure}{r}{6.4cm}
    \centering
    \includegraphics[width=6cm]{figures/q2_control}
    \caption{Control signals for the three designs. NB: differing vertical
    scales, $u_1$ on the left axis and $u_2,u_3$ on the right. Different 
    time scales capture the transient. }
    \label{fig:q2_control}
\end{wrapfigure}
The fist design solves the aforementioned problem in the previous section, with no objective nor constraints. At there there's no further optimization, the solver returns any arbitrarily feasible point. The found gain is $K_a$, which is available in the accompanying code. 

The following gain places the closed-loop eigenvalues at $\lambda_{1,2}(A+BK_a)=-0.500$. Simulated on the nonlinear plant, this gain drives the commands far outside the admissible range, as reported in (Table \ref{tab:q2_cmp}).


\subsection{Guaranteed convergence rate}
To command the speed of response we impose a certified convergence rate $\alpha$, by replacing~\eqref{eq:q2_baselmi} by \citep[eq.~(26)]{zaccarian_lmi}
\begin{equation}
    W\succ I,\qquad \operatorname{He}(AW+BX)\preceq-2\alpha W ,
    \label{eq:q2_lmi_b}
\end{equation}

where $\alpha$ is the exponential decay of the state norm $|\delta\bm{x}(t)|\le M\,e^{-\alpha t}|\delta\bm{x}(0)|$. The restrained normalization $W \succ I$ does not restrict feasibility and will be used in the next section.

Choosing $\alpha = \alpha_b = 2$ yields $K_b$, which is available in the accompanying code. It yields eigenvalues at $\lambda_{1,2}(A+BK_b)=-2.207$, which are smaller than the required $-2$ value imposed by the condition. The response is faster but the gain is an order of magnitude larger, thus the controls are even more out of bound. 

\subsection{Minimum-gain design with actuator authority}

The third design fixed the rate and minimizes a bound on the gain, so that the commands remain within the actuator range. 


Following
\citep[eq.~(27)]{zaccarian_lmi},
\begin{equation}
    \begin{aligned}
        \min_{W,X,\kappa}\quad & \kappa \qquad
        \text{s.t.}\quad & W\succ I,\qquad
            \operatorname{He}(AW+BX)\preceq-2\alpha W, \qquad
        & \begin{bmatrix} \kappa I & X^\top\\ X & \kappa I \end{bmatrix}\succ0 .
    \end{aligned}
    \label{eq:q2_lmi_c}
\end{equation}

A Schur complement on the last block gives $X^\top X \prec \kappa^2 \ I$, hence $\| X \| \le \kappa$; moreover, we have the condition that $W \succ 1$. These two conditions bound the gain $\| K \| = \|XW^-1\| \le \|X\| \le \kappa$. The rate has been manually tuned to $\alpha_c = 0.0115$, which makes sure that no control exceeds the bound $[0,1]$ on the nonlinear plant. The solution of the constrained optimization is reported in the accompanying code, and has eigenvalues at $\lambda_{1,2}(A+BK_c)=-0.0115=-\alpha_c$. Although its performances are far from the ones given by $Ka$ and $Kb$, one could appreciate through Table \ref{tab:q2_cmp} that the settling time is just $\approx 13\%$ of the uncontrolled settling time, which represents a sevenfold improvement on the open loop.

\subsection{Comparison}

Table \ref{tab:q2_cmp} has been produced in order to provide a quantitative evaluation of each gain's performance. 
Additionally, Figures \ref{fig:q2_state_short} and \ref{fig:q2_state_long} show the behaviour of the system with and without state-feedback: the first figure highlights the overshoot, while the second shows settling times of the slower dynamics. Finally, Figure \ref{fig:q2_control} shows the control trajectories using the three different gains. The best performance is provided by the second controller which employs $K_b$, while the controller which could be implemented in reality is the third one, although, as will be shown in the next sections, it's not robust to sensor failure. Using the third controller allows to improve the open loop response by approximately seven times, while still respecting actuator feasibility. 

\begin{table}[htbp]
    \centering
    \small
    \renewcommand{\arraystretch}{1.1}
    \begin{tabular}{l r r r}
        \toprule
        Controller & $t_s$ [s] & $u_{\text{max}}$ & $u_{\text{min}}$ \\
        \midrule
        Open loop            & 2840 & 0.660 & \phantom{-}0.320 \\
        $K_a$ (arbitrary)    &    8 & 38.06 & --1.323         \\
        $K_b$ (rate)         &    2 & 178.99& --4.896          \\
        $K_c$ (min-gain)     &  385 & 0.868 & \phantom{-}0.012  \\
        \bottomrule
    \end{tabular}
    \caption{Closed-loop comparison on the nonlinear plant, from $\delta\bm{x}_0=-0.05$~m. Settling time $t_s$ is the time to enter and remain within $2\%$ of the initial step size; a $1.4$~s floor is applied to avoid the overshoot being mistaken with settling time. The commands of $K_a$ and $K_b$ leave the admissible range $[0,1]$, so their settling times are not physically realisable; only $K_c$ respects the actuator bounds.$u_{\text{max}}$ and $u_{\text{min}}$ are considered across all three controls. }
    \label{tab:q2_cmp}
\end{table}


\begin{figure}[tb]
    \centering
    \begin{subfigure}[t]{6cm}
        \centering
        \includegraphics[height=7cm]{figures/q2_state_short}
        \caption{Short timescale ($0$--$15$~s): $K_a$, $K_b$ settle quickly
        but overshoot; $K_c$ rises gently.}
        \label{fig:q2_state_short}
    \end{subfigure}
    \hspace{2cm}
    \begin{subfigure}[t]{6cm}
        \centering
        \includegraphics[height=7cm]{figures/q2_state_long}
        \caption{Full timescale ($0$--$3000$~s): $K_c$ settles in
        $\approx385$~s against the open-loop $\approx2840$~s.}
        \label{fig:q2_state_long}
    \end{subfigure}

    \caption{Closed-loop state response on the nonlinear plant from
    $\delta\bm{x}_0=-0.05$~m, shown at two timescales. The apparent
    instantaneous settling of $K_a$ and $K_b$ in~(b) is obtained only with
    the inadmissible commands of Figure~\ref{fig:q2_control}.}
    \label{fig:q2_state}
\end{figure}



\RuleSection{Observer Design}
\label{sec:q3}

\subsection{Which sensor failure is recoverable}
If the sensor for $h_1$ fails, while sensor for $h_2$ remains functional, the matrix $C$ becomes $C_2 = [0\quad 1]$. In the other case, we'll have matrix $C_1 = [1\quad 0]$. For each case we can check whether the full state can still be recovered: i.e. whether the pair $(A,C_i)$ is observable. 

The observability matrix has rank $2$ for $C_2$ and rank $1$ for $C_1$. This tells us that a failure in sensor $h_1$ is recoverable while a failure in $h_2$ is not recoverable. The reason is the cascade structure of $A$. The matrix is lower-triangular, so $h_1$ does not depend on $h_2$ while $h_2$ does depend on $h_1$.  From now on, $\tilde C = [0 \quad 1]$.


\subsection{Observer gain}

The observer mirrors the model by working with a state estimate $\widehat{\delta \bm x}$. During its update, two component are summed together: the dynamics of the estimated state and a correction term, which depends on the difference among measured output and predicted output: 

\begin{equation}
    \mathsf{d/dt}\Big[\widehat{\delta\bm{x}}\Big]
    = A\,\widehat{\delta\bm{x}} + B\,\delta\bm{u}
      + L\,\big(\delta\bm{y}- \widehat y\big)
    = A\,\widehat{\delta\bm{x}} + B\,\delta\bm{u}
      + L\,\big(\delta\bm{y}-C\,\widehat{\delta\bm{x}}\big),
    \label{eq:q3_obs}
\end{equation}

The correction term is interesting at it is used as a proxy for the current "goodness" of our estimation. The matrix $L$ tells not only the sensibility that the result has on the state, it also tells the direction in which the estimate should be corrected. 

Moreover, the estimation error $\bm{e}=\delta\bm{x}-\widehat{\delta\bm{x}}$ follows $\dot{\bm{e}}=(A-LC)\bm{e}$, therefore, we want $A-LC$ to have stable eigenvalues with a chosen rate, in oder to have asymptotic stability towards zero. 

This is the dual of the state-feedback problem, because the estimation-error matrix $A-LC$ has the same eigenvalues as its transpose $(A-LC)^\top = A^\top - C^\top L^\top$. Assigning those eigenvalues by the choice of $L$ is therefore the same task as choosing a state-feedback gain for the pair $(A^\top,C^\top)$, which was exactly the problem of Section~\ref{sec:q2} with $A\to A^\top$ and $B\to C^\top$. Therefore, we solve the dual LMI:
\begin{equation}
    P\succ0,\qquad \operatorname{He}(PA-YC)\preceq-2\alpha_{\mathrm{obs}}P,
    \qquad L=P^{-1}Y .
    \label{eq:q3_lmi}
\end{equation}

With $\alpha_{\mathrm{obs}}=2$ the observer eigenvalues come out at $-3.254$ (both equal), well below the required $-2$. The estimate therefore converges much faster than the controller, whose eigenvalues are at $-0.0115$. It's worth noting that in this work we take the difference $A-LC$, meanwhile in \cite{zaccarian_lmi} the summation $A+LC$ is taken into consideration. Once the sign of the elements of $L$ are adjusted, the two formulations are equivalent.

\subsection{Estimation}

The observer is tested on the nonlinear plant with a sensor failure on sensor $h_1$. The initial deviation coordinates are $\delta \bm x = [-0.05, -0.05]$: $h_2$ is initialized at the correct value while $h_1$ is initialized at the equilibrium value $h_1(0) := h_1^{\text{eq}} = 0.2$.

Figure \ref{fig:q3_track} shows the result. The gain $\alpha_2$ dictates the kink in the trajectory at $t = 1.3s$: the trajectory then reaches within $1$ mm of true $h_1$ after $164.3$ s, as the tanks approach the operating point. The lag of the observer model is due to the fact that the observer is based on the linearized version of the model: if the observer is run on the linear state, the correction happens at $1.30$ s. 

\begin{figure}[tb]
    \centering
    \begin{subfigure}[t]{0.49\linewidth}
        \centering
        \includegraphics[height=7cm]{figures/q3_track}
        \caption{Observer estimate against the true state (only $h_2$
        measured). The $h_2$ estimate is exact from the start; the $h_1$
        estimate starts at $0.2$~m, dips near $t=1.3$~s, then converges.}
        \label{fig:q3_track}
    \end{subfigure}
    \hfill
    \begin{subfigure}[t]{0.49\linewidth}
        \centering
        \includegraphics[height=7cm]{figures/q3_fb}
        \caption{Closed loop with true-state against observer-based feedback,
        both using $K_c$. Observer-based is slightly faster on $h_1$, then
        matches as the estimate becomes exact.}
        \label{fig:q3_fb}
    \end{subfigure}
    \caption{Observer behaviour: state estimation~(a) and the resulting
    observer-based closed loop~(b).}
    \label{fig:q3}
\end{figure}


\section{Observer-based feedback}
With no $h_1$ sensor, the controller must use its estimate. By the separation principle the design splits in two. In the coordinates $(\delta\bm{x},\bm{e})$, where $\bm{e}=\delta\bm{x}-\widehat{\delta\bm{x}}$ is the estimation error, the closed-loop matrix is block-triangular, so its eigenvalues are those of the controller $A+BK_c$ together with those of the observer $A-LC$. The controller gives $\lambda_{1,2}=-0.0115$ and the observer gives $\lambda_{1,2}=-3.254$, so the closed loop has eigenvalues \footnotesize $\{-0.0115,\,-0.0115,\,-3.254,\,-3.254\}$, \normalsize that is, the union of the two sets: the system is therefore stable. The controller will now use $K_c$ from section \ref{sec:q2}.  The control law is:

\begin{equation}
    \bm{u}=\bar{\bm{u}}+K_c\big(\widehat{\delta\bm{x}}\,\big)
          =\bar{\bm{u}}+K_c(\hat{\bm{x}}-\bar{\bm{x}}).
    \label{eq:q3_fb}
\end{equation}

Figure \ref{fig:q3_fb} compares the observer-based loop with the full-state loop of section \ref{sec:q2}. Both reach operating point, although the observer based response is faster on $h_1$. The reason follows from Figure \ref{fig:q3_track}: $h_1$ is underestimated, driving the controller to fill more aggressively the first tank. 

Table \ref{tab:q3_cmp} makes the comparison precise: the observer-based loop settles faster on both tanks, although at the cost of actuators saturation. The smallest command drops to $-0.094$, which is outside of feasible range. A new gain $K_d$ is computed, in order to get a feasible gain compatible with the observation model distance from the real model. The manual tuning settles on $\alpha_d=0.0095$, a little below $\alpha_c$, which gives $K_d$. Its numeric value is available in the accompanying code. The associated eigenvalues are $\lambda_{1,2}=-0.0095$. 

The results of the simulation are stated in Table~\ref{tab:q3_cmp}: the controls are now feasible ($u_{\min}=0.004$), at the cost of a slower
response. Settling rises from $192$ to $255$~s on $h_1$ and from $367$ to $453$~s on $h_2$, about a third slower than the infeasible observer-based $K_c$.



\begin{table}[htbp]
    \centering
    \small
    \renewcommand{\arraystretch}{1.1}
    \begin{tabular}{l r r r r}
        \toprule
        Loop & \large $t_{s,h_1}$ \normalsize[s] & \large $t_{s,h_2}$ \normalsize [s] & $u_{\text{max}}$ & $u_{\text{min}}$ \\
        \midrule
        Full-state feedback ($K_c$)   & 303 & 385 & 0.868 & \phantom{$-$}0.012 \\
        Observer-based ($K_c$)        & 192 & 367 & 0.965 & --0.094         \\
        Observer-based ($K_d$)        & 255 & 453 & 0.891 & \phantom{$-$}0.004 \\
        \bottomrule
    \end{tabular}
    \caption{Full-state against observer-based feedback, both using $K_c$,
    from $\delta\bm{x}_0=-0.05$~m (same $1.4$~s settling floor as
    Table~\ref{tab:q2_cmp})}
    \label{tab:q3_cmp}
\end{table}

% ====================================================================== 
\bibliographystyle{abbrvnat}
\bibliography{references}   % <- replace with your .bib filename (no extension)

\end{document}
```

# selected_results/fig_01.png

*(binary file, content omitted)*

# selected_results/fig_02.png

*(binary file, content omitted)*

# selected_results/find_u2.png

*(binary file, content omitted)*

# selected_results/h2 flow decomposition.png

*(binary file, content omitted)*

# src/io/load_params.m

```matlab
function cfg = load_params()

root = fileparts(fileparts(fileparts(mfilename('fullpath')))); % src/io -> root
configDir = fullfile(root, 'config');

cfg.plant  = read_json(fullfile(configDir, 'plant_params.json'));
cfg.eq  = read_json(fullfile(configDir, 'equilibrium.json'));
cfg.design  = read_json(fullfile(configDir, 'design.json'));
cfg.observer_gains = read_json(fullfile(configDir, 'observer_gains.json'));

cfg.meta.root      = root;
cfg.meta.configDir = configDir;
cfg.meta.timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS'); %#ok<TNOW1,DATST>

end

% ------------------------------------------------------------------------
function s = read_json(f)
if ~isfile(f)
    error('load_params:missingFile', 'Config file not found: %s', f);
end
s = jsondecode(fileread(f));
end
```

# src/io/snapshot_run.m

```matlab
function outdir = snapshot_run(cfg, data, figs, tag)
% Record a run so results are traceable to their parameters.


stamp = datestr(now, 'yyyy-mm-dd_HHMMSS'); %#ok<TNOW1,DATST>

outdir = fullfile(cfg.meta.root, 'results', [stamp tag]);
if ~exist(outdir, 'dir'), mkdir(outdir); end

% Copy the exact configuration used.
jsons = {'plant_params.json', 'design.json', 'equilibrium.json', 'observer_gains.json'};
for i = 1:numel(jsons)
    src = fullfile(cfg.meta.configDir, jsons{i});
    if isfile(src), copyfile(src, outdir); end
end

% Save the data payload.
if ~isempty(fieldnames(data))
    save(fullfile(outdir, 'run.mat'), '-struct', 'data');
end

% Save figures.
for i = 1:numel(figs)
    base = fullfile(outdir, sprintf('fig_%02d', figs(i).Number));
    savefig(figs(i), [base '.fig']);
    try
        exportgraphics(figs(i), [base '.png'], 'Resolution', 300);
    catch
    end
end

% Human-readable stamp.
fid = fopen(fullfile(outdir, 'run_info.txt'), 'w');
fprintf(fid, 'Run recorded: %s\n', stamp);
fclose(fid);

fprintf('Run recorded in %s\n', outdir);

end
```

# src/io/style_fig.m

```matlab
function style_fig(fig, opts)
%STYLE_FIG  Make a figure legible at print size.
%
%   style_fig(fig)
%   style_fig(fig, Name, Value, ...)
%
% Name-value options (MATLAB plot style):
%   'Size'       [W H] in centimetres           (default [8 8])
%   'FontSize'   axis/label/title font, points  (default 9)
%   'MaxXTicks'  max number of x ticks          (default 5)
%   'MaxYTicks'  max number of y ticks          (default 5)
%
% Examples:
%   style_fig(figUs)                              % 8x8 cm, 9 pt, <=5 ticks
%   style_fig(figCmp, 'Size',[16 16])             % 2x2 panel figure
%   style_fig(figU,   'Size',[8 24])              % tall 3-stack, text beside
%   style_fig(figErr, 'MaxXTicks',4, 'MaxYTicks',3)
%
% snapshot_run exports with exportgraphics at the figure's on-screen size,
% so a figure set to [8 8] here becomes an 8x8 cm PNG whose text is 9 pt.
% Included in LaTeX at width=8cm (no scaling) the text stays 9 pt.

    arguments
        fig
        opts.Size      (1,2) double {mustBePositive} = [8 8]
        opts.FontSize  (1,1) double {mustBePositive} = 9
        opts.MaxXTicks (1,1) double {mustBeInteger, mustBePositive} = 5
        opts.MaxYTicks (1,1) double {mustBeInteger, mustBePositive} = 5
    end

    W  = opts.Size(1);   H  = opts.Size(2);
    fs = opts.FontSize;

    % Fix the physical size so exportgraphics writes W x H cm.
    set(fig, 'Color', 'w', ...
             'Units', 'centimeters', 'Position', [3 3 W H], ...
             'PaperUnits', 'centimeters', 'PaperSize', [W H], ...
             'PaperPositionMode', 'manual', 'PaperPosition', [0 0 W H]);

    ax = findall(fig, 'Type', 'axes');
    for k = 1:numel(ax)
        a = ax(k);
        set(a, 'FontName', 'Helvetica', 'FontSize', fs, ...
               'LineWidth', 0.5, 'TickLength', [0.02 0.02]);
        thin_ticks(a, 'X', opts.MaxXTicks);
        thin_ticks(a, 'Y', opts.MaxYTicks);
        if ~isempty(a.Title),  a.Title.FontSize  = fs; a.Title.FontWeight = 'bold'; end
        if ~isempty(a.XLabel), a.XLabel.FontSize = fs; end
        if ~isempty(a.YLabel), a.YLabel.FontSize = fs; end
    end

    set(findall(fig, 'Type', 'legend'),   'FontSize', max(fs-1, 7));
    set(findall(fig, 'Type', 'colorbar'), 'FontSize', fs);
end

function thin_ticks(a, axname, nmax)
% Keep at most nmax ticks, preserving MATLAB's "nice" tick values.
    t = a.([axname 'Tick']);
    if numel(t) > nmax
        step = ceil(numel(t) / nmax);
        a.([axname 'Tick']) = t(1:step:end);
    end
end
```

# src/model/f_TankNLDynamics.m

```matlab
function dx = f_TankNLDynamics(t, x, u, p)

% Nonlinear dynamics, to be given to ODE45. 

h1 = x(1); 
h2 = x(2);

s1 = sqrt(2 * p.g * h1);
s2 = sqrt(2 * p.g * h2);

dx = zeros(2, 1);
dx(1) = ( u(1)*p.q0 - u(2)*p.a2*s1 ) / p.A1;
dx(2) = ( u(2)*p.a2*s1 - u(3)*p.a3*s2 ) / p.A2;

end
```

# src/model/f_TankStateSpace.m

```matlab
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
```

# src/model/find_equilibrium.m

```matlab
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
```

# src/sim/make_input.m

```matlab
function uFun = make_input(K, xeq, ueq)

% Build a state-feedback law u(t,x) = ueq + K*(x - xeq).

% Control law:  u = ueq + K*(x - xeq)
%
% K is designed on the linearised model in deviation coordinates
%   dx = x - xeq,  du = u - ueq,   d(dx)/dt = A*dx + B*du,
% where the LMI returns K such that du = K*dx stabilises it. Writing that
% relation back in absolute variables gives u = ueq + K*(x - xeq).

xeq = xeq(:); ueq = ueq(:);
if isempty(K)
    uFun = @(t, x) ueq;
else
    uFun = @(t, x) ueq + K * (x(:) - xeq);
end

end
```

# src/sim/run_lin.m

```matlab
function res = run_lin(lin, cfg, t_sim, dt)
% Integrate the linearised model in deviation coordinates.


A = lin.A; xeq = lin.xeq(:);  ueq = lin.ueq(:);

tspan = 0:dt:t_sim;
dx0   = cfg.design.dx0(:);


[t, dx] = ode45(@(t,dx) A*dx, tspan, dx0, ...
                odeset('RelTol',1e-6,'AbsTol',1e-8));

res.t  = t;
res.dx = dx;
res.x  = xeq' + dx;
res.u = repmat(ueq', numel(t), 1);

end
```

# src/sim/run_nl.m

```matlab
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
```

# src/sim/run_observer.m

```matlab
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
```

# src/sim/run_observer_lin.m

```matlab
function res = run_observer_lin(lin, K, L, xhat0, cfg, t_sim, dt)
% Integrate the coupled true+observer dynamics on the linear plant.

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
```

# startup_project.m

```matlab
function startup_project()

    root = fileparts(mfilename('fullpath'));
    
    addpath(genpath(fullfile(root, 'src')));
    
    fprintf('Project ready. Root: %s\n', root);

end
```
