function S = repair_test_observed_arm(mode, outDir)
% REPAIR_TEST_OBSERVED_ARM  Which failing populations of the census before failure can the repair
% test use? Reports the populations of run_repair_test_observed, why some are skipped, and the
% restarts of the observed arm, with and without the defector that later took over.
%
% Produces : SI S4.3, the failing populations with no member at the restart point (printed block
%            [0], the populations skipped and why)
% Inputs   : mode, 'main' (the only one); outDir (optional, default the repair_test_observed folder
%              under the data root's publication folder)
%            publication/repair_test_observed/repair_test_observed.mat   (run_repair_test_observed)
% Outputs  : publication/repair_test_observed/repair_test_observed.txt   (the printed report: [0] the
%              populations, repair size, starting cooperation)
% Run      : repair_test_observed_arm('main')   (about 1 min)
% Seed     : none (no random draw)
%
% Arms: A observed, with ten copies of the defector; A0 without them.

    if nargin < 1 || isempty(mode), mode = 'main'; end
    assert(strcmp(mode, 'main'), 'repair_test_observed_arm:mode', 'mode must be ''main''');
    if nargin < 2 || isempty(outDir), outDir = fullfile(grn_data_file('publication'), 'repair_test_observed'); end
    L = load(fullfile(outDir, 'repair_test_observed.mat'), 'out'); out = L.out;
    txt = fullfile(outDir, 'repair_test_observed.txt');
    if exist(txt, 'file'), delete(txt); end
    diary(txt);
    S = tabulate_runs(out);
    print_header(out, S);
    diary off;
end


% =========================================================================
function S = tabulate_runs(out)
% Per population i and arm: S.p collapse share, S.g1 cooperation at generation 1, S.seedOK whether
% every run started from the population given.
    ev = out.events; runIdx = find([ev.inRun]);
    S.ev = ev(runIdx); n = numel(runIdx);
    arms = out.cfg.arms;
    for a = 1:numel(arms)
        A = arms{a};
        [S.p.(A), S.g1.(A)] = deal(nan(n, 1));
        S.seedOK.(A) = true(n, 1);
        for i = 1:n
            j = find([out.runs.event] == S.ev(i).num & strcmp({out.runs.arm}, A));
            if isempty(j), continue; end
            r = out.runs(j).r;
            S.p.(A)(i) = mean(r.collapse);
            S.g1.(A)(i) = mean(r.hc(:, 1));
            S.seedOK.(A)(i) = all(r.seedUsed);
        end
    end
    S.arms = arms;
end


% =========================================================================
function print_header(out, S)
    c = out.cfg; ev = out.events;
    fprintf('REPAIR TEST, OBSERVED ARM (MAIN)  seed %d; m = %d; nGen = %d; k = %d; lambda = %.3f; arms %s\n', ...
        c.rngSeed, c.m, c.nGen, c.k, c.lambda, strjoin(c.arms, ' '));
    fprintf('Generators experiments/run_repair_test_observed.m, analysis/s4_losing_cooperation/repair_test_observed_arm.m\n');
    fprintf('Census input %s\n', c.censusFile);
    fprintf('Run time %.1f min for %d jobs on %d workers\n', out.runMin, out.nJobs, c.nWorkers);
    fprintf('\n[0] EVENTS: %d census collapses with a defector that later took over; usable %d; skipped %d\n', ...
        numel(ev), sum([ev.usable]), sum(~[ev.usable]));
    for i = find(~[ev.usable])
        fprintf('  skipped %s: %s (shares at the restart point mem/coop/oth %.3f/%.3f/%.3f)\n', ev(i).id, ev(i).skipReason, ev(i).shares);
    end
    fprintf('  the five run-time events, left out of the main run: %s\n', strjoin(out.pilotIds, ', '));
    fprintf('  in this run: %d events\n', numel(S.ev));
    ok = cellfun(@(A) all(S.seedOK.(A)(isfinite(S.p.(A)))), S.arms);
    fprintf('  seeded start used in every run (gen-1 snapshot == restart population with the defector added): %s\n', ...
        ternary(all(ok), 'YES', 'NO'));
    rp = [S.ev.nReplaced]; md = [S.ev.meanRepairDist];
    fprintf('  repair: cooperators that do not retaliate replaced per event median %d [%d, %d]; mean W distance to the member median %.3f\n', ...
        round(median(rp)), min(rp), max(rp), median(md, 'omitnan'));
    if isfield(S.ev, 'coop1A')
        a = [S.ev.coop1A]; b = [S.ev.coop1B];
        fprintf('  STARTING COOPERATION under the runs'' encounter rule entering gen 1 (no defector added), from the restart point: A (observed) median %.3f [%.3f, %.3f]; B (repaired) median %.3f [%.3f, %.3f]\n', ...
            median(a), min(a), max(a), median(b), min(b), max(b));
        g1 = cellfun(@(A) mean(S.g1.(A), 'omitnan'), {'A', 'A0'});
        fprintf('    recorded gen-1 cooperation in the runs (mean over events): A %.3f | A0 %.3f\n', g1);
    end
end
