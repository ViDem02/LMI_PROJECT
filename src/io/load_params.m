function cfg = load_params()

root = fileparts(fileparts(fileparts(mfilename('fullpath')))); % src/io -> root
configDir = fullfile(root, 'config');

cfg.plant  = read_json(fullfile(configDir, 'plant_params.json'));
cfg.eq  = read_json(fullfile(configDir, 'equilibrium.json'));
cfg.design  = read_json(fullfile(configDir, 'design.json'));

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