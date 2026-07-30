
%% Init
clear all

startup_project;
cfg = load_params;

%% Linearization

p = cfg.plant; 
eq = cfg.eq;

xeq = [eq.h1; eq.h2];
ueq = [eq.u1; eq.u2; eq.u3];

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


lambda = eig(A);
fprintf('open-loop eigenvalues : %.3e , %.3e\n', lambda(1), lambda(2));
% both of them have negative real parts. 

fprintf('slowest time constant : %.0f s\n', 1/min(abs(real(lambda))));


%% Simulation


% create the input
uFun = @(t, x) ueq;

% model for run_lin
lin  = struct('A', A, 'B', B, 'E', E, 'xeq', xeq, 'ueq', ueq);   

% run non linear simulation
t_sim = 2500;
dt    = 0.5;


%initial height close to eq VS initial height far from eq
cases = { [-0.02; -0.02], 'close (h^{init} = 0.18 m)';
[-0.18; -0.18], 'far (h^{init} = 0.02 m)' };

res = struct('nl', {}, 'lin', {}, 'x0', {}, 'label', {});
for i = 1:size(cases, 1)
    dx0 = cases{i, 1};
    sim_cfg = cfg;
    sim_cfg.design.x0 = dx0;    

    % same scenario on each model, same time grid
    rnl  = run_nl(uFun, sim_cfg, t_sim, dt);    % nonlinear
    rlin = run_lin(lin, sim_cfg, t_sim, dt);    % linear 

    res(i).nl    = rnl;
    res(i).lin   = rlin;
    res(i).x0    = xeq + dx0;
    res(i).label = cases{i, 2};
end



%% Comparison figure

cNL  = [0      0.4470 0.7410];   % nonlinear (blue)
cLIN = [0.8500 0.3250 0.0980];   % linear    (red)
figCmp = figure(1); clf(figCmp); figCmp.Name = 'Q1: linear vs nonlinear';
tiledlayout(2, 2, 'TileSpacing', 'compact');
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
        title(sprintf('h_%d - %s', k, res(i).label));
        if i == 1 && k == 1
            legend('nonlinear', 'linear', 'equilibrium', 'Location', 'best');
        end
    end
end