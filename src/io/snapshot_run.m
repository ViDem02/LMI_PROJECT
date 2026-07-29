function outdir = snapshot_run(cfg, data, figs, tag)
% Record a run so results are traceable to their parameters.


stamp = datestr(now, 'yyyy-mm-dd_HHMMSS'); %#ok<TNOW1,DATST>

outdir = fullfile(cfg.meta.root, 'results', [stamp tag]);
if ~exist(outdir, 'dir'), mkdir(outdir); end

% Copy the exact configuration used.
jsons = {'plant_params.json'};
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