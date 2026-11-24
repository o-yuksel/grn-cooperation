function Pu = report_added_input(J, rs, cfg)
%REPORT_ADDED_INPUT  Prints the drift back after the split, by arm, and returns its series.
%
% Implements : population-median alpha, ||x*||, beta_rms relative to CONTROL at the split, and the best
%              genotype's alpha and ||x*||; used by run_tracking_held_still and run_tracking_moving_target
% Inputs     : J (the arms' runs), rs (replicates), cfg (arms)
% Outputs    : Pu (per arm: ra, rx, rbeta, bestA, bestX)
% Seed       : none

    fprintf('\n---- drift back after the split: population-median alpha, ||x*||, beta_rms as a ratio to CONTROL at the split (median over replicates) ----\n');
    ks = [J(1).track.k]; nK = numel(ks);
    Pu = struct();
    for a = 1:4
        ra = nan(numel(rs), nK); rx = ra; be = ra;
        for q = 1:numel(rs)
            c = J([J.r] == rs(q) & [J.arm] == 1); z = J([J.r] == rs(q) & [J.arm] == a);
            if isempty(z) || z.void, continue; end
            ra(q, :) = [z.track.alpha] / c.track(1).alpha;
            rx(q, :) = [z.track.xn] / c.track(1).xn;
            be(q, :) = [z.track.beta_rms] / c.track(1).beta_rms;
        end
        Pu.(matlab.lang.makeValidName(cfg.arms{a})) = struct('ra', ra, 'rx', rx, 'rbeta', be);
        [bA, bX, kw] = best_series(J, rs, a);
        Pu.(matlab.lang.makeValidName(cfg.arms{a})).bestA = bA;
        Pu.(matlab.lang.makeValidName(cfg.arms{a})).bestX = bX;
        sel = unique([1, find(mod(ks, 30) == 0), nK]);
        fprintf('%-10s gen+%s: alpha %s | ||x*|| %s\n', cfg.arms{a}, mat2str(ks(sel)), ...
            mat2str(median(ra(:, sel), 1, 'omitnan'), 3), mat2str(median(rx(:, sel), 1, 'omitnan'), 3));
        bk = find(isfinite(median(be, 1, 'omitnan')));
        fprintf('           beta_rms / CONTROL split at gen+%s: %s\n', mat2str(ks(bk)), ...
            mat2str(median(be(:, bk), 1, 'omitnan'), 3));
        fprintf('           best genotype alpha, ||x*|| / CONTROL split ThBest at gen+%s: %s | %s\n', mat2str(kw), ...
            mat2str(median(bA, 1, 'omitnan'), 3), mat2str(median(bX, 1, 'omitnan'), 3));
    end
end
