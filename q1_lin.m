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