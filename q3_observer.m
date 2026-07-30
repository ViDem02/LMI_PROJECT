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