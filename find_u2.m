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