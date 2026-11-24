function failed_population_histories()
% FAILED_POPULATION_HISTORIES  Did failed populations ever reach cooperation? Reads the cooperation
% history of every failed population in the four signaling architectures.
%
% Produces : the failed-population rows of Table S3 (number, exceeded 0.90, median peak), SI S2.5,
%            and main text, "Mutual sensing lets gene networks escape the dilemma" (every failed
%            saturating-contrast population first exceeds 0.90)
% Inputs   : publication/<architecture>_100reps.mat   (written by run_signaling_architectures)
% Outputs  : printed only. (1) saturating contrast, failed populations: peak and final cooperation,
%              shares whose cooperation exceeded 0.90, 0.75, 0.50, median generation of the peak,
%              and how many of all 100 exceeded 0.90 and ended below 0.50; (1b), a second table printed after
%              (1), per architecture:
%              failed populations, median and minimum peak, share whose peak exceeded 0.90, median
%              peak generation
% Run      : failed_population_histories   (about 6 min)
% Seed     : none (no random draws)
%
% Failed: mean cooperation in the last generation below 0.95. "Exceeded 0.90": cooperation above
% 0.90 at some generation, the SI's rule ("after exceeding 0.90"; Table S3, "Exceeded 0.90").
d = load(grn_data_file('publication/saturating_contrast_100reps.mat'), 'results');
R = d.results;
n = numel(R);
nsnap = size(R{1}.historyM, 4);
finalCoop = zeros(n, 1); maxCoop = zeros(n, 1); argmaxGen = zeros(n, 1);
for r = 1:n
    c = R{r}.historyCooperation;
    finalCoop(r) = c(end);
    [maxCoop(r), argmaxGen(r)] = max(c);
end
succ = finalCoop >= 0.95;
fail = ~succ;
fprintf('\nreps=%d  snapshots=%d  success(final>=0.95)=%d  fail=%d\n', ...
        n, nsnap, sum(succ), sum(fail));

% (1) The cooperation history of the failed saturating-contrast runs
fprintf('\n(1) cooperation trajectory, FAILED runs (n=%d)\n', sum(fail));
fprintf('    max ever  mean=%.3f median=%.3f min=%.3f | final mean=%.3f median=%.3f\n', ...
        mean(maxCoop(fail)), median(maxCoop(fail)), min(maxCoop(fail)), ...
        mean(finalCoop(fail)), median(finalCoop(fail)));
for thr = [0.9 0.75 0.5]
    fprintf('    exceeded %2.0f%%: %3.0f%% of failed runs\n', ...
            100*thr, 100*mean(maxCoop(fail) > thr));
end
fprintf('    exceeded 90%% then ended <50%%: %d of %d reps\n', ...
        sum(maxCoop > 0.90 & finalCoop < 0.50), n);
fprintf('    median generation of peak, failed runs: %d of %d\n', ...
        round(median(argmaxGen(fail))), numel(R{1}.historyCooperation));

% (1b) The same history, in a second table, for the failed runs of each of the four architectures (Table S3, SI S2.5)
fprintf('\n(1b) cooperation trajectory of FAILED runs, all four architectures\n');
sweep = { 'graded-separate',     'graded_separate_100reps.mat';
          'graded-contrast',     'graded_contrast_100reps.mat';
          'saturating-separate', 'saturating_separate_100reps.mat';
          'saturating-contrast', 'saturating_contrast_100reps.mat' };
fprintf('    %-20s %5s %8s %8s %9s %8s\n', ...
        'architecture', 'nfail', 'medpeak', 'minpeak', 'peak>0.9', 'pkgen');
for a = 1:size(sweep, 1)
    ds = load(grn_data_file(['publication/' sweep{a,2}]), 'results');
    Rs = ds.results; ns = numel(Rs);
    fc = zeros(ns,1); mc = zeros(ns,1); pg = zeros(ns,1);
    for r = 1:ns
        c = Rs{r}.historyCooperation;
        fc(r) = c(end);
        [mc(r), pg(r)] = max(c);
    end
    fl = fc < 0.95;
    fprintf('    %-20s %5d %8.3f %8.3f %8.0f%% %8d\n', sweep{a,1}, sum(fl), ...
            median(mc(fl)), min(mc(fl)), 100*mean(mc(fl) > 0.90), ...
            round(median(pg(fl))));
    clear ds Rs;
end

fprintf('\nDONE\n');
end
