function opts = initOpts(opts)
    %% ================== File info ==========================
    % Author            : Tiep Vu (http://www.personal.psu.edu/thv102/)
    % Time created      : Wed Jan 27 23:57:31 2016
    % Last modified     : Wed Jan 27 23:57:32 2016
    % Description       : 
    %     INPUT:
    %
    %     OUTPUT: 
    %
    %% ================== end File info ==========================
    %% show, max_iter, checkgrad
    if ~isfield(opts, 'show_progress')
        opts.show_progress = 0;
    end 
    %%
    if ~isfield(opts, 'checkgrad')
        opts.checkgrad = 0;
    end 
    %%
    if ~isfield(opts, 'max_iter')
        opts.max_iter = 100;
    end 
end