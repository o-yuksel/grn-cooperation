function c = coop_clone_batch(Wm)
%COOP_CLONE_BATCH  Cooperation of many wirings with their clones under the runs' encounter rule.
%
% Implements : mutate_with_refusal('coop', W, W, config) for many W at once: 100 encounters of W against its
%              clone, evolve_population's CPU path of simulateInteractionBatch (saturating contrast:
%              contrast sensing, Hill 1, sensor gain 1, mu 0.1, K_M 1, dt 1, 300 steps, stability every 40 at 0.01,
%              snapshot rule as there), W rounded to single as the evolve_population call does; pairs are
%              independent, so batching does not change any pair's trajectory; used by run_safety_search and
%              safety_search_losses
% Inputs     : Wm (4 x 4 x n wirings)
% Outputs    : c (n x 1, the share of the 100 encounters in which W cooperates)
% Seed       : the 100 start-state pairs drawn once under rng(20261022, twister) (rand(4,100)*0.01 focal then
%              opponent, as in evolve_population); the caller's stream is restored

    persistent X0F X0O
    if isempty(X0F)
        s = rng; rng(20261022, 'twister'); X0F = rand(4, 100) * 0.01; X0O = rand(4, 100) * 0.01; rng(s);
    end
    n = size(Wm, 3); W = double(single(Wm)); np = 100 * n;
    Wa = W(:, :, repelem(1:n, 100));
    xf = repmat(X0F, 1, n); xo = repmat(X0O, 1, n);
    fs = xf; os = xo; stF = false(1, np); stO = false(1, np); aIdx = 1:np; xfo = xf; xoo = xo;
    xF = xf; xO = xo; Wf = Wa; last = 0;
    for t = 1:300
        na = numel(aIdx);
        inF = reshape(pagemtimes(Wf, reshape(xf, 4, 1, na)), 4, na);
        inO = reshape(pagemtimes(Wf, reshape(xo, 4, 1, na)), 4, na);
        dO = xo(1, :) + xo(2, :) + 1e-12; dF = xf(1, :) + xf(2, :) + 1e-12;
        sigF = xo(1, :) ./ dO - xo(2, :) ./ dO; sigO = xf(1, :) ./ dF - xf(2, :) ./ dF;
        inF(3, :) = inF(3, :) + sigF; inF(4, :) = inF(4, :) - sigF;
        inO(3, :) = inO(3, :) + sigO; inO(4, :) = inO(4, :) - sigO;
        inF = max(inF, 0); inO = max(inO, 0);
        xf = xf + (inF ./ (1 + inF) - 0.1 * xf); xo = xo + (inO ./ (1 + inO) - 0.1 * xo);
        xf = max(xf, 0); xo = max(xo, 0);
        if t - last >= 40
            mF = max(max(xf, [], 1), 1e-9); mO = max(max(xo, [], 1), 1e-9);
            cF = sum((xf - xfo).^2, 1) ./ mF.^2 < 0.01^2; cO = sum((xo - xoo).^2, 1) ./ mO.^2 < 0.01^2;
            both = cF & cO;
            fs(:, aIdx(both)) = xf(:, both); os(:, aIdx(both)) = xo(:, both); stF(aIdx(both)) = true; stO(aIdx(both)) = true;
            nf = cF & ~both; fs(:, aIdx(nf)) = xf(:, nf); stF(aIdx(nf)) = true;
            no = cO & ~both; os(:, aIdx(no)) = xo(:, no); stO(aIdx(no)) = true;
            if any(both)
                keep = ~both; aIdx = aIdx(keep); Wf = Wf(:, :, keep); xf = xf(:, keep); xo = xo(:, keep);
                if isempty(aIdx), break; end
            end
            xfo = xf; xoo = xo; last = t;
        end
    end
    xF(:, aIdx) = xf; xO(:, aIdx) = xo;
    fs(:, ~stF) = xF(:, ~stF);
    c = mean(reshape(fs(1, :) > fs(2, :), 100, n), 1)';
end
