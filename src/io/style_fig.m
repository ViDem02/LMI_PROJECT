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