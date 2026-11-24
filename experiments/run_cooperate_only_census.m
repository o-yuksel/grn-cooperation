function run_cooperate_only_census(stage, varargin)
% RUN_COOPERATE_ONLY_CENSUS  How often does one mutation end cooperation in evolved cooperators that
% do not retaliate? Runs the per-mutation census (its snapshots, mutation distribution and code) on
% the wirings that cooperate with their clone but fail the membership test.
%
% Produces : SI S3.3, "Cooperators that do not retaliate lose cooperation less often than members"
%            (input: the rows that cooperate_only_losses('full') reads and prints)
% Inputs   : publication/<architecture>_100reps.mat   (runs of the four architectures; written by
%              run_signaling_architectures)
% Outputs  : publication/cooperate_only_census/full_cooperate_only.mat   (one row per scored wiring, with the
%              mutations drawn for it and which end its cooperation; written by cooperate_only_census)
% Run      : run('startup.m'); run_cooperate_only_census('full')   (all 100 replicates of each
%              architecture; 6.8 min on 2 workers); run_cooperate_only_census(
%              'full', reps) takes one replicate list per architecture ([] for all 100)
% Seed     : rng(20261016 + 1000 m + r + 500000) for architecture m and replicate r, set in
%            cooperate_only_census on the default generator of a pool worker (Threefry)
%
% Scored: readout against the clone at least 1e-6 and ||x*|| (size of the settled expression state)
% in [0.5, 14], the census's rule. 2 workers at one thread, below-normal priority. Run startup.m first.

    mainData = grn_data_file('publication');
    outDir = fullfile(mainData, 'cooperate_only_census');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    set_priority(feature('getpid'));

    cfg = struct('dataDir', mainData);
    switch stage
        case 'full'
            cfg.group = 'cooperate_only'; cfg.modes = 1:4; cfg.reps = {[], [], [], []};   % all 100
            if ~isempty(varargin), cfg.reps = varargin{1}; end
            cfg.outFile = fullfile(outDir, 'full_cooperate_only.mat');
        otherwise
            error('unknown stage %s', stage);
    end
    start_pool();
    fprintf('[%s] stage %s start\n', datestr(now, 'HH:MM:SS'), stage);
    cooperate_only_census(cfg);
    fprintf('[%s] stage %s done\n', datestr(now, 'HH:MM:SS'), stage);
end


function start_pool()
    p = gcp('nocreate');
    if isempty(p), p = parpool('Processes', 2); end
    spmd
        maxNumCompThreads(1);
        set_priority(feature('getpid'));
    end
    if ispc
        fprintf('pool: %d workers, 1 thread each, BelowNormal\n', p.NumWorkers);
    else
        fprintf('pool: %d workers, 1 thread each\n', p.NumWorkers);
    end
end


function set_priority(pid)  % lowers the worker's priority through PowerShell on Windows; no-op elsewhere
    if ~ispc, return; end
    system(sprintf(['powershell -NoProfile -Command "(Get-Process -Id %d).PriorityClass = ' ...
        '''BelowNormal''"'], pid));
end
