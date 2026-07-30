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