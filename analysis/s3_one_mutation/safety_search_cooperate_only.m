function safety_search_cooperate_only(d)
% SAFETY_SEARCH_COOPERATE_ONLY  After the safety search that need only cooperate with the clone, do
% the final wirings still retaliate, still respond (alpha), and lose cooperation to any of their
% 1,000 test mutations?
%
% Produces : main text, "Conditional cooperation stays within a mutation's reach of the
%            threshold"; Discussion; SI S3.4; Table S7
% Inputs   : d   the safety-search folder (default publication/safety_search), holding
%              long_flags.mat (safety_search_losses) and search_setup.mat, search_main_starts_1_4.mat,
%              search_main_starts_5_12.mat and long_jobs/long_<k>_C.mat (run_safety_search)
% Outputs  : publication/safety_search/cooperate_only_report.txt   (the lines this run prints: one line per
%              final wiring and the summary counts)
% Run      : safety_search_cooperate_only   (a few seconds)
% Seed     : none (no random draws)
%
% Read from the stored records, per start: alpha at the start and the end of the search,
% cooperation with the default defector at the end (100 encounters under the runs' encounter
% rule), the final wiring's defecting single-weight mutants and membership at the last
% checkpoint, and its share of the 1,000 test mutations that end cooperation (long_flags.mat).
% Alpha is compared with the median final alpha of the 1,000-step searches that must also
% retaliate (long_jobs/long_<k>_S.mat; alpha undefined left out), as run_safety_search('longreport')
% prints it.
    if nargin < 1, d = fullfile(grn_data_file('publication'), 'safety_search'); end
    F = load(fullfile(d, 'long_flags.mat')); s = load(fullfile(d, 'search_setup.mat'));
    R = load_c_runs(d); aS = median_alpha_s(d);

    %% Report
    report = @() report_lines(F, s, R, aS); %#ok<NASGU> called inside evalc
    rep = evalc('report();');
    fprintf('%s', rep);

    %% Save
    fid = fopen(fullfile(d, 'cooperate_only_report.txt'), 'w'); fprintf(fid, '%s', rep); fclose(fid);
end


function report_lines(F, s, R, aS)
% One line per final wiring, then the summary counts; prints the report's lines.
    fprintf('rep | alpha start end | coopDef end | nDef start end | member end | coop lost\n');
    T = zeros(numel(R), 6);
    for j = 1:numel(R)
        r = R{j}; i = find(F.k == r.k & F.grp == 3, 1);
        T(j, :) = [r.ends{2}.alpha, r.ends{2}.coopDefRef, r.fresh.nDefMut(end), r.fresh.member(end), F.lostC(i), r.ends{1}.alpha];
        fprintf('%3d | %6.3f %6.3f | %.2f | %2d %2d | %d | %.3f\n', s.reps(r.k), T(j, 6), T(j, 1), T(j, 2), r.nDef0, T(j, 3), T(j, 4), T(j, 5));
    end
    zero = T(:, 5) == 0; a = T(:, 1);
    fprintf('finals %d; reach zero cooperation loss %d\n', numel(R), sum(zero));
    fprintf('cooperate with the default defector in every encounter: %d of %d (others %s)\n', sum(T(:, 2) == 1), numel(R), mat2str(T(T(:, 2) < 1, 2)'));
    fprintf('alpha median start %.3f, end %.3f (defined %d of %d)\n', median(T(:, 6)), median(a(~isnan(a))), sum(~isnan(a)), numel(a));
    fprintf('among the %d that reach zero: alpha >= %.3f in %d; alpha above its start in %d; undefined %d\n', ...
        sum(zero), aS, sum(zero & a >= aS), sum(zero & a > T(:, 6)), sum(zero & isnan(a)));
    fprintf('members %d of %d (all with no defecting mutant: %d); non-members %d, all with defecting mutants: %d\n', ...
        sum(T(:, 4) == 1), numel(R), all(T(T(:, 4) == 1, 3) == 0), sum(T(:, 4) == 0), all(T(T(:, 4) == 0, 3) > 0));
end


function aS = median_alpha_s(d)
% The median alpha at the end of the 1,000-step searches that must also retaliate.
    L = dir(fullfile(d, 'long_jobs', 'long_*_S.mat')); a = nan(numel(L), 1);
    for i = 1:numel(L)
        x = load(fullfile(L(i).folder, L(i).name)); r = x.r;
        if size(r.trace, 1) == 1000, a(i) = r.ends{2}.alpha; end
    end
    aS = median(a, 'omitnan');
end


function R = load_c_runs(d)
% Records of the search that need only cooperate (`mode` 'C'): those in search_main_starts_1_4.mat and
% search_main_starts_5_12.mat, then every long_jobs/long_*_C.mat.
    R = {};
    for nm = {'search_main_starts_1_4.mat', 'search_main_starts_5_12.mat'}
        y = load(fullfile(d, nm{1}));
        for i = 1:numel(y.R), if strcmp(y.R{i}.mode, 'C'), R{end + 1} = y.R{i}; end, end %#ok<AGROW>
    end
    L = dir(fullfile(d, 'long_jobs', 'long_*_C.mat'));
    for i = 1:numel(L), x = load(fullfile(L(i).folder, L(i).name)); R{end + 1} = x.r; end %#ok<AGROW>
end
