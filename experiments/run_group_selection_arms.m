function run_group_selection_arms(nReps, nWorkers)
% RUN_GROUP_SELECTION_ARMS  Runs the three group-selection arms (group selection, the random-split
% control, undivided), 185 populations of 1,000 each for 500 generations, then their analysis.
%
% Produces : the runs and the analysis behind SI S2.7, Table S2 and Fig. 3 (group selection);
%            Methods, "Game and evolution"
% Inputs   : nReps      populations per arm (default 185, the paper's runs)
%            nWorkers   pool size (default 2)
% Outputs  : the three files of run_group_selection (phase 'main') and those of
%            group_selection_arms (tag 'N1000'; with nReps other than 185,
%            'N1000_<nReps>reps')
% Run      : run_group_selection_arms   (one pool of 2 workers); run_group_selection_arms(3) runs
%            3 populations per arm
% Seed     : 770001-770185 (group selection), 780001-780185 (random-split control), 790001-790185
%            (undivided), set in run_group_selection

    if nargin < 1 || isempty(nReps), nReps = 185; end
    if nargin < 2 || isempty(nWorkers), nWorkers = 2; end
    tAll = tic;
    o1 = run_group_selection('fission', 'main', nReps, 500, nWorkers);
    o2 = run_group_selection('random_split', 'main', nReps, 500, nWorkers);
    o3 = run_group_selection('control', 'main', nReps, 500, nWorkers);
    fprintf('TOTAL simulation wall %.1f min\n', toc(tAll) / 60);
    delete(gcp('nocreate'));
    tA = tic;
    if nReps == 185
        group_selection_arms();
    else
        % group_selection_arms names the random-split arm by its stored spelling (its stored analysis)
        group_selection_arms(struct('fission', o1.file, stored_value('groupSelectionMode', 'random_split', 'stored'), o2.file, ...
            'control', o3.file), ...
            sprintf('N1000_%dreps', nReps));
    end
    fprintf('analysis wall %.1f min; TOTAL %.1f min\n', toc(tA) / 60, toc(tAll) / 60);
end
