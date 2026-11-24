function tipping_flag_vs_input_size()
% TIPPING_FLAG_VS_INPUT_SIZE  Does the member's own tipping input matter, or would the size of the
% input a mutation adds, |x_j dW_ij|, or of the mutation, |dW_ij|, catch as many losses when it
% flags as many mutations?
%
% Produces : main text, "Conditional cooperation stays within a mutation's reach of the
%            threshold"; Methods, "Mutation tests in the game"; SI S3.3, "The member's own
%            tipping input matters beyond the size of the input a mutation adds"
% Inputs   : publication/tipping_input/tipping_input_7000.mat   (tipping_input_flag, stage 'run')
% Outputs  : publication/tipping_input/tipping_input_baseline.mat   (OUT: per side the recall
%              and precision of each rule with intervals, and the flag share f)
%            prints per side the loss rate, f, recall and precision of each rule with intervals,
%              the paired recall differences, and the same per architecture
% Run      : tipping_flag_vs_input_size   (about 30 s)
% Seed     : ties in size fall in a random order (rng(20261016, 'twister')); the bootstraps use
%            rng(20261016 + 9 + s, 'twister') for side s (1 cooperation, 2 retaliation); both
%            around a saved and restored stream
%
% Sides: cooperation (loss: the mutant no longer cooperates with its clone) and retaliation (loss:
% the mutant's readout against the worst-case mutant H is positive). f is the tipping-input
% flag's share of draws, weighted by the individuals carrying the wiring; each size rule flags
% draws in decreasing size until its weighted share reaches f. Recall is the share of losses
% flagged, precision the share of flagged draws that are losses; a random flag would have recall
% f. Flags are set once on the pooled sample and held fixed in the bootstrap (2,000 resamples of
% members within architecture, the differences paired); per-architecture recall uses them too.

    O = grn_data_file('publication/tipping_input');
    L = load(fullfile(O, 'tipping_input_7000.mat')); S = L.S; R = L.R; K = numel(S(1).dW);
    names = {'graded_separate','graded_contrast','saturating_separate','saturating_contrast'};
    ar = [S.arch]'; w = [S.count]'; n = numel(w);
    get = @(fn) cell2mat(cellfun(@(o) o.(fn)(:)', R, 'UniformOutput', false));
    dl = cell2mat(arrayfun(@(s) s.dW(:)', S, 'UniformOutput', false));
    lost = cell2mat(arrayfun(@(s) s.lost(:)', S, 'UniformOutput', false));
    retL = cell2mat(arrayfun(@(s) (s.rH(:) > 0)', S, 'UniformOutput', false));
    W2 = repmat(w, 1, K);
    fprintf('[DATA] %d members (per architecture %s), %d draws, K = %d\n', n, mat2str(accumarray(ar, 1, [4 1])'), numel(dl), K);
    fprintf('[CENSUS SHARES, this sample] clone-weighted share of ALL draws: lose cooperation %.4f; lose retaliation against H %.4f; either %.4f; both %.4f\n', ...
        wm(lost, W2), wm(retL, W2), wm(lost | retL, W2), wm(lost & retL, W2));

    % Per side: the stored tipping-input flags, the losses, and the size of the input each draw
    % adds, |dW_ij x_j| with x_j at the member's state in that encounter
    sides = {'cooperation', get('predC'), lost, abs(dl .* get('xjC')); ...
             'retaliation vs H', get('predR'), retL, abs(dl .* get('xjR'))};
    for s = 1:2
        P = sides{s, 2}; Y = sides{s, 3}; push = sides{s, 4};
        f = wm(P, W2);
        F = {P, top_share(push, W2, f), top_share(abs(dl), W2, f)};
        lab = {'TIPPING INPUT', 'SIZE-ONLY', 'MUTATION-SIZE-ONLY'};
        fprintf('\n[%s] loss rate %.4f; TIPPING INPUT flag share f = %.4f\n', upper(sides{s, 1}), wm(Y, W2), f);
        [pt, ci] = boot(F, Y, w, ar, 2000, 20261016 + 9 + s);
        for k = 1:3
            fprintf('  %-19s flag share %.4f | recall %.4f [%.4f, %.4f] | precision %.4f [%.4f, %.4f]\n', lab{k}, ...
                wm(F{k}, W2), pt(k, 1), ci(k, 1, 1), ci(k, 1, 2), pt(k, 2), ci(k, 2, 1), ci(k, 2, 2));
        end
        fprintf('  %-19s flag share %.4f | recall %.4f (expected) | precision %.4f (expected)\n', 'RANDOM', f, f, wm(Y, W2));
        fprintf('  TIPPING INPUT - SIZE-ONLY recall %+.4f [%+.4f, %+.4f] | TIPPING INPUT - MUTATION-SIZE-ONLY recall %+.4f [%+.4f, %+.4f]\n', ...
            pt(4, 1), ci(4, 1, 1), ci(4, 1, 2), pt(5, 1), ci(5, 1, 1), ci(5, 1, 2));
        OUT(s) = struct('side', sides{s, 1}, 'pt', pt, 'ci', ci, 'f', f, 'rows', {{'TIPPING INPUT', 'SIZE-ONLY', 'MUTATION-SIZE-ONLY', 'TIPPING INPUT - SIZE-ONLY', 'TIPPING INPUT - MUTATION-SIZE-ONLY', 'RANDOM (flag share)'}}); %#ok<AGROW>
        for m = 1:4
            km = ar == m;
            rc = cellfun(@(Fk) rec(Fk(km, :), Y(km, :), W2(km, :)), F);
            fprintf('    %-20s flag share %.4f | recall TIPPING INPUT %.4f, SIZE-ONLY %.4f, MUTATION-SIZE-ONLY %.4f, RANDOM %.4f\n', ...
                names{m}, wm(P(km, :), W2(km, :)), rc, wm(P(km, :), W2(km, :)));
        end
    end
    save(fullfile(O, 'tipping_input_baseline.mat'), 'OUT');   % row 6: flag share and interval
end


function F = top_share(score, W2, f)
% Flag the highest scores until the weighted share first reaches f; ties fall in a random order
% (a permutation under rng(20261016, 'twister'), then MATLAB's stable sort).
    rs = rng; rng(20261016, 'twister'); p = randperm(numel(score)); rng(rs);
    [~, o] = sort(score(p), 'descend'); idx = p(o);          % stable sort after a random permutation
    cw = cumsum(W2(idx)) / sum(W2(:));
    k = find(cw >= f, 1);
    F = false(size(score)); F(idx(1:k)) = true;
end


function [pt, ci] = boot(F, Y, w, ar, B, seed)
% Rows: recall/precision of F{1..3}; row 4 = F1 - F2, row 5 = F1 - F3 (recall in col 1, precision in col 2);
% row 6 = F1's clone-weighted flag share (the RANDOM rule's expected recall), col 1.
    W = repmat(w, 1, size(Y, 2)); n = numel(w);
    a = cellfun(@(Fk) sum(W .* (Fk & Y), 2), F, 'UniformOutput', false);
    q = cellfun(@(Fk) sum(W .* Fk, 2), F, 'UniformOutput', false);
    y = sum(W .* Y, 2); tot = sum(W, 2);
    st = @(c) stat(c, a, q, y, tot);
    pt = st(ones(n, 1));
    rs = rng; rng(seed, 'twister'); bs = nan(6, 2, B);
    for b = 1:B
        c = zeros(n, 1);
        for m = unique(ar)'
            ix = find(ar == m); c = c + accumarray(ix(randi(numel(ix), numel(ix), 1)), 1, [n 1]);
        end
        bs(:, :, b) = st(c);
    end
    rng(rs); ci = cat(3, prctile(bs, 2.5, 3), prctile(bs, 97.5, 3));
end


function M = stat(c, a, q, y, tot)
    M = zeros(6, 2);
    for k = 1:3, M(k, :) = [c' * a{k} / (c' * y), c' * a{k} / (c' * q{k})]; end
    M(4, :) = M(1, :) - M(2, :); M(5, :) = M(1, :) - M(3, :);
    M(6, :) = [c' * q{1} / (c' * tot), NaN];
end


function r = rec(F, Y, W2), r = sum(W2(:) .* (F(:) & Y(:))) / sum(W2(:) .* Y(:)); end
function p = wm(Z, W2), p = sum(W2(:) .* Z(:)) / sum(W2(:)); end
