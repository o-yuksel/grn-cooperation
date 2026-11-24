function run_fixed_opponent_census(nWorkers, dataDir)
% RUN_FIXED_OPPONENT_CENSUS  For every member of the per-mutation census: its worst-case mutant, and
% how each census mutation moves its readout against that mutant held fixed (the input of Fig. 4A).
%
% Produces : Fig. 4A (input: the census that sample_census_members reads, whose members
%            tipping_input_flag scores for the data fig4_tipping_input draws)
% Inputs   : publication/rho_by_architecture.mat   (the per-mutation census; per_mutation_census)
%            publication/retaliation_margin_vs_rho.mat   (written by retaliation_margin_vs_rho)
%            publication/<architecture>_100reps.mat   (runs of the four architectures;
%              run_signaling_architectures)
% Outputs  : publication/retaliation_fixed_opponent/retaliation_fixed_opponent_census.mat and .txt
%              (the census and its printed report), and ..._census_ck<m>.mat (one checkpoint per
%              architecture m; a rerun resumes from them)
% Run      : run_fixed_opponent_census   (the census on 4 workers, the publication folder of the
%              data root; about 1 h, stored report 57.0 min); run_fixed_opponent_census(nWorkers,
%              dataDir) sets both
% Seed     : none here; measure_retaliation_fixed_opponent uses rng(20260928 + k, 'twister'), one
%            fixed k per draw set, and rng(0) for the census's own half-loss bootstrap
%
% Opens a pool of nWorkers if none is open; the client and every worker run one thread each.

    if nargin < 1 || isempty(nWorkers), nWorkers = 4; end
    if nargin < 2 || isempty(dataDir), dataDir = grn_data_file('publication'); end
    maxNumCompThreads(1);
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('local', nWorkers); end
    wait(parfevalOnAll(pool, @maxNumCompThreads, 1, 1));
    f = parfevalOnAll(pool, @maxNumCompThreads, 1); wait(f); thr = fetchOutputs(f);
    fprintf('RUNNER: pool %d workers; worker maxNumCompThreads %s; client %d\n', ...
        pool.NumWorkers, mat2str(thr(:)'), maxNumCompThreads);
    cfg = struct('fullCensus', true, 'nWorkers', pool.NumWorkers, 'dataDir', dataDir, 'pilot', false);
    measure_retaliation_fixed_opponent(cfg);
end
