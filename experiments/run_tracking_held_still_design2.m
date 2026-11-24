function out = run_tracking_held_still_design2(cfg)
% RUN_TRACKING_HELD_STILL_DESIGN2  The population test of the load with the target held still in
% design 2 (output self-activation a = 0.20), with 30 generations without mutation at the end.
%
% Produces : the design-2 runs of SI S5.1 with 30 generations without mutation (read by
%            tracking_held_still_checks and run_tracking_generations_without_mutation)
% Inputs   : those of run_tracking_held_still: publication/tracking_runs/
%              tracking_runs_N1000.mat (written by run_tracking_task) and
%              publication/tracking_single_genotypes/lowering_settings.mat (written by run_tracking_single_genotypes)
% Outputs  : publication/tracking_held_still_30gen/design2/main.mat and main.txt (written by
%              run_tracking_held_still)
% Run      : run_tracking_held_still_design2(struct('tolerance', 0.0110))   (about 46 min on 2 workers)
% Seed     : those of run_tracking_held_still (replicate r continues with rng(20261001 + r, 'twister'))
%
% Calls run_tracking_held_still unchanged with cfg.designIndex = 3, design 2's index among the four
% designs. The tolerance is that of the runner's check that the control's mean fitness reproduces
% the stored one; run_tracking_held_still_design2(struct('pilot', true)), a short run, prints it.

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    cfg.designIndex = 3;                                                  % design 2, a = 0.20
    if ~isfield(cfg, 'outDir')
        cfg.outDir = fullfile(grn_data_file('publication'), ...
            'tracking_held_still_30gen/design2');
    end
    out = run_tracking_held_still(cfg);
end
