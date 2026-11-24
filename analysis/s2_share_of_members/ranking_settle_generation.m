% RANKING_SETTLE_GENERATION  When does the ranking of the four signaling architectures settle?
% Finds the last generation at which their ranking by mean cooperation, and by success, differs
% from their ranking at generation 500.
%
% Produces : the two generations, and the saturating-contrast mean cooperation over five windows of
%            100 generations, at generation 500 and at its maximum (SI S2.4); printed only
% Inputs   : publication/<architecture>_100reps.mat   (run_signaling_architectures; the runs
%              Fig. 2 draws)
% Outputs  : console lines (no file written)
% Run      : ranking_settle_generation   (not timed)
% Seed     : none (no random draws)
%
% Mean cooperation at a generation: the mean over replicates, in percent; success at a
% generation: the share of replicates at or above 0.95 then.

dataDir = grn_data_file('publication');
modes = {'graded_separate', 'graded_contrast', 'saturating_separate', 'saturating_contrast'};
nG = 500;
M = zeros(4, nG); S = zeros(4, nG);
for m = 1:4
    d = load(fullfile(dataDir, [modes{m} '_100reps.mat']), 'results');
    n = numel(d.results); C = zeros(n, nG);
    for r = 1:n, C(r, :) = d.results{r}.historyCooperation(1:nG); end
    M(m, :) = mean(C, 1) * 100;
    S(m, :) = mean(C >= 0.95, 1) * 100;
end
[~, finM] = sort(M(:, nG)); [~, finS] = sort(S(:, nG));
lastM = 0; lastS = 0;
for g = 1:nG
    [~, o] = sort(M(:, g)); if ~isequal(o, finM), lastM = g; end
    [~, o] = sort(S(:, g)); if ~isequal(o, finS), lastS = g; end
end
fprintf('mean-cooperation ranking last differs from its generation-%d ranking at generation %d\n', nG, lastM);
fprintf('success ranking last differs at generation %d\n', lastS);
% Windows of generations; adjacent windows share their end generation
w = [1 100; 100 200; 200 300; 300 400; 400 500];
for k = 1:size(w, 1)
    fprintf('saturating contrast mean, generations %d-%d: %.1f\n', w(k, 1), w(k, 2), mean(M(4, w(k, 1):w(k, 2))));
end
[mx, ix] = max(M(4, :));
fprintf('saturating contrast: %.1f at generation %d; maximum %.1f at generation %d\n', M(4, nG), nG, mx, ix);
