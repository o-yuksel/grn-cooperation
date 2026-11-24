function T = group_share_of_members(fpath)
% GROUP_SHARE_OF_MEMBERS  The share of members in the last generation of a group-selection run,
% and the share of cooperators that do not retaliate.
%
% Implements : the share of members of Methods, "Success, share of members and statistics", for
%              the group arms (Table S2, SI S2.7); called by group_selection_tests
% Inputs     : fpath   a run file of a group arm (cfg, res, finalCoop): written by
%                      run_group_selection (N = 1,000), or the N = 400 files (part of the data
%                      archive)
% Outputs    : T   struct, per replicate and pooled: wirings sampled, members (field `nVm`),
%                  wirings meeting the cooperation condition (C1), cooperators that do not
%                  retaliate (C1 without the retaliation condition C2, `nC1noC2`), the share of
%                  members (stored as under its older name) and the share of cooperators that do not retaliate (of all
%                  sampled and of the C1 wirings); also printed as a table
% Seed       : rng(42) (the default generator) once; then per replicate, in order, randperm picks
%              30 of its last-generation individuals; membership_test seeds each wiring itself
%              and restores the stream
%
% The membership test runs with the run's settings, 600 steps, the simulate_pair route, mutation
% size at its default 0.2. Every replicate of the file is sampled.

    d = load(fpath, 'cfg', 'res', 'finalCoop');
    cfg = d.cfg;
    reps = 1:numel(d.res);
    params = struct( ...
        'sensingMode',       cfg.sensingMode, ...
        'hillCoefficient',   cfg.hillCoefficient, ...
        'discreteOutput',    cfg.discreteOutput, ...
        'michaelisConstant', cfg.michaelisConstant, ...
        'degradationRate',   cfg.degradationRate, ...
        'timeStep',          cfg.timeStep, ...
        'developmentSteps',  600, ...
        'backend',           'simulate_pair');

    nSub = 30;
    rng(42);
    nR = numel(reps);
    nIn = zeros(nR, 1); nC1 = zeros(nR, 1); nC1noC2 = zeros(nR, 1); nS = zeros(nR, 1);
    t0 = tic;
    for i = 1:nR
        rep = d.res{reps(i)};
        nAg = size(rep.finalNetworkB, 3);
        for a = randperm(nAg, min(nSub, nAg))
            W = double(rep.finalNetworkM(:,:,a)) .* double(rep.finalNetworkB(:,:,a));
            c = membership_test(W, params);
            nS(i) = nS(i) + 1;
            nIn(i) = nIn(i) + c.is_css;
            nC1(i) = nC1(i) + c.is_weak;
            nC1noC2(i) = nC1noC2(i) + (c.is_weak && ~c.is_css);
        end
    end
    secs = toc(t0);

    T = struct('file', fpath, 'reps', reps(:), 'finalCoop', d.finalCoop(reps(:)), ...
        'nSampled', nS, 'nVm', nIn, 'nC1', nC1, 'nC1noC2', nC1noC2, ...
        'occupancy', sum(nIn) / sum(nS), 'shareC1noC2', sum(nC1noC2) / sum(nS), ...
        'shareC1noC2_givenC1', sum(nC1noC2) / max(sum(nC1), 1), 'seconds', secs);

    [~, name] = fileparts(fpath);
    fprintf('\n%s\n%5s %9s %8s %5s %5s %9s\n', name, 'rep', 'finalCoop', 'sampled', 'mem', 'C1', 'C1notC2');
    for i = 1:nR
        fprintf('%5d %9.3f %8d %5d %5d %9d\n', reps(i), T.finalCoop(i), nS(i), nIn(i), nC1(i), nC1noC2(i));
    end
    fprintf(['pooled: share of members %.4f (%d/%d) | C1 without C2 %.4f of sampled, %.4f of C1 | ' ...
             '%.1f s for %d wirings\n'], T.occupancy, sum(nIn), sum(nS), T.shareC1noC2, ...
             T.shareC1noC2_givenC1, secs, sum(nS));
end
