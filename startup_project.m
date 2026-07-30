function startup_project()

    root = fileparts(mfilename('fullpath'));
    
    addpath(genpath(fullfile(root, 'src')));
    
    fprintf('Project ready. Root: %s\n', root);

end
