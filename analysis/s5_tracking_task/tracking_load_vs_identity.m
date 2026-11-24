function out = tracking_load_vs_identity(cfg)
% TRACKING_LOAD_VS_IDENTITY  In the tracking task, how much fitness does one mutation remove from
% each design's generation-350 best genotype, with the target moving and held still, and how close
% does the congruence identity come to it without solving any mutant?
%
% Produces : the measured and predicted per-mutation loads (Fig. 6A; Table S13; Methods,
%            "Non-social tracking task"; main text, "Responsiveness carries a mutation load in a
%            non-social trait"); the .mat is also read by tracking_load_excess,
%            run_tracking_single_genotypes and fig6_load_cells
% Inputs   : publication/tracking_runs/tracking_runs_N<nPop>.mat (genotypes and
%              scoring windows; written by run_tracking_task) and tracking_runs_N<nPop>_alpha.mat
%              in the same folder (written by tracking_responsiveness; read if present)
%            at nPop = 80 only: publication/tracking_runs/tracking_runs_N80_load.mat (part
%              of the data archive), compared with the first replicates (`routing_check`)
%            cfg (optional): nPop (1000; 80 is the N = 80 run), betaDraws (2000), popDraws (400),
%              popCap (25), doPop (true), nBoot (10000), rngSeed (20260917), nWorkers (5),
%              repList (all), inFile, alphaFile, outFile, txtFile, refFile
% Outputs  : publication/tracking_runs/tracking_runs_N<nPop>_load.mat and _load.txt
%              (the console text)
% Run      : tracking_load_vs_identity()   (about 45 min on 5 workers)
% Seed     : draws on replicate r's best genotype in design f (its index in aList): rng(20260917 +
%            1000 f + r, 'twister'); the population sample (`doPop`): rng(20260917 + 500000 +
%            1000 f + r, 'twister'); the bootstraps use their own fixed streams
%
% w is a genotype's fitness exp(-gamma ||p - o||^2) (gamma = 50, p the output, o the optimum),
% averaged over the scoring window (the last 36 generations, each at its own input and optimum).
% L = 1 - E[w(mutant)] / w(parent) over betaDraws draws of the mutation distribution (W33 and W44
% not mutated): the fitness one application removes, not a load per generation. L_identLin
% replaces each draw's output shift by the identity's first-order shift C Phi dW x* (C Phi, the
% outputs' response to each gene's input, SI S1.2) and scores it with the same fitness.


    if nargin < 1 || isempty(cfg), cfg = struct(); end
    dd = fullfile(grn_data_file('publication'),'tracking_runs');
    if ~isfield(cfg,'nPop'),      cfg.nPop      = 1000; end
    if ~isfield(cfg,'inFile'),    cfg.inFile    = fullfile(dd, sprintf('tracking_runs_N%d.mat', cfg.nPop)); end
    if ~isfield(cfg,'alphaFile'), cfg.alphaFile = fullfile(dd, sprintf('tracking_runs_N%d_alpha.mat', cfg.nPop)); end
    if ~isfield(cfg,'outFile'),   cfg.outFile   = fullfile(dd, sprintf('tracking_runs_N%d_load.mat', cfg.nPop)); end
    if ~isfield(cfg,'txtFile'),   cfg.txtFile   = fullfile(dd, sprintf('tracking_runs_N%d_load.txt', cfg.nPop)); end
    if ~isfield(cfg,'refFile'),   cfg.refFile   = fullfile(grn_data_file('publication'),'tracking_runs','tracking_runs_N80_load.mat'); end
    if ~isfield(cfg,'betaDraws'), cfg.betaDraws = 2000; end   % draws per best genotype
    if ~isfield(cfg,'popDraws'),  cfg.popDraws  = 400; end    % draws per individual, population read
    if ~isfield(cfg,'popCap'),    cfg.popCap    = 25; end     % individuals read; >= N: every wiring
    if ~isfield(cfg,'doPop'),     cfg.doPop     = true; end
    if ~isfield(cfg,'nBoot'),     cfg.nBoot     = 10000; end
    if ~isfield(cfg,'rngSeed'),   cfg.rngSeed   = 20260917; end
    if ~isfield(cfg,'nWorkers'),  cfg.nWorkers  = 5; end
    if ~isfield(cfg,'repList'),   cfg.repList   = []; end

    L = load(cfg.inFile); D = L.out;
    P = D.P; aList = D.aList; nA = numel(aList);
    repList = cfg.repList; if isempty(repList), repList = 1:D.nReps; end
    nR = numel(repList);
    % Stored arm names: 'moving' (target moving) and 'fixed' (target held still)
    arms = {'moving','fixed'};

    diary_file(cfg.txtFile, true);
    fprintf('\n=========== TRACKING TASK, MUTATION LOAD, N = %d, %d replicates ===========\n', P.n_pop, nR);
    fprintf('source %s (run finished %s)\n', cfg.inFile, D.finishedAt);
    fprintf('N = %d, %d generations, mutation p = %.2f at sd = %.3f, s = %.1f\n', ...
        P.n_pop, P.n_generations, P.mutation_prob, P.mutation_sd, P.selection_strength);
    fprintf('per-entry p * 16 entries * N = %.1f (nominal; W33 and W44 are not mutated)\n', 16*P.mutation_prob*P.n_pop);
    fprintf('draws: %d per best genotype; population read %s\n', cfg.betaDraws, ...
        popdesc(cfg, P.n_pop));

    start_pool(cfg.nWorkers);

    nJ = numel(arms)*nA*nR;
    jobs = cell(nJ,1); ja = zeros(nJ,1); jf = ja; jr = ja;
    q = 0;
    for ai = 1:numel(arms)
        for f = 1:nA
            for r = 1:nR, q = q+1; ja(q) = ai; jf(q) = f; jr(q) = r; end
        end
    end
    ThB = cell(nJ,1); PopT = cell(nJ,1); PopF = cell(nJ,1);
    WU = cell(nJ,1); WO = cell(nJ,1); WW = cell(nJ,1);
    for q = 1:nJ
        R = D.(arms{ja(q)}); f = jf(q); rr = repList(jr(q));
        ThB{q} = R.ThBest{f,rr}; PopT{q} = R.popTh{f,rr}; PopF{q} = R.popFreq{f,rr};
        WU{q} = R.winU{f}; WO{q} = R.winO{f}; WW{q} = R.winW{f};
    end
    seedBase = cfg.rngSeed; nD = cfg.betaDraws; nP = cfg.popDraws;
    popCap = cfg.popCap; doPop = cfg.doPop; rl = repList; npop = P.n_pop;
    t0 = tic;
    parfor q = 1:nJ
        jobs{q} = one_cell(ThB{q}, PopT{q}, PopF{q}, npop, WU{q}, WO{q}, WW{q}, P, ...
                           seedBase, jf(q), rl(jr(q)), nD, nP, popCap, doPop);
    end
    fprintf('\n%d cells measured in %.1f min\n', nJ, toc(t0)/60);

    mets = {'L','L0','w_par','err_rms','beta_win_rms','beta_win_mean','dLpred', ...
            'beta_ident','L_ident','L_identLin','beta_lin_rms','alpha_op','xnorm_op','L_first', ...
            'Lpop','L0pop','w_parpop','beta_pop_rms','nDistinct','nEval','batchErr'};
    S = struct();
    for ai = 1:numel(arms)
        for m = 1:numel(mets), S.(arms{ai}).(mets{m}) = nan(nA,nR); end
    end
    for q = 1:nJ
        z = jobs{q}; ar = arms{ja(q)}; f = jf(q); r = jr(q);
        for m = 1:numel(mets), S.(ar).(mets{m})(f,r) = z.(mets{m}); end
    end

    out = struct('cfg',cfg,'aList',aList,'nPop',P.n_pop,'nReps',nR,'repList',repList, ...
                 'arms',{arms},'S',S,'sourceFinished',D.finishedAt);

    fprintf('\n---- solver check: steady_batch against tracking.steady, max over cells ----\n');
    for ai = 1:numel(arms)
        fprintf('[%-6s] max|steady_batch - tracking.steady| = %.3g\n', arms{ai}, ...
            max(S.(arms{ai}).batchErr(:)));
    end

    out.routing = routing_check(S, arms, aList, cfg, repList);
    report_cells(S, arms, aList, nR, cfg);
    out.paired = report_paired(S, arms, aList, nR, cfg);
    A = [];
    if exist(cfg.alphaFile,'file'), LA = load(cfg.alphaFile); A = LA.out; end
    out.relation = report_relation(S, A, D, arms, aList, P, repList);
    report_confound(S, arms, aList);

    out.finishedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    save(cfg.outFile, 'out');
    fprintf('\nSaved %s\n', cfg.outFile);
    diary_file(cfg.txtFile, false);
end


% =====================================================================
function s = popdesc(cfg, N)
    if ~cfg.doPop
        s = 'SKIPPED';
    elseif cfg.popCap >= N
        s = sprintf('the whole population (%d individuals) at %d draws each', N, cfg.popDraws);
    else
        s = sprintf('a uniform subsample of %d of %d individuals at %d draws each (budget cut, see the header)', ...
            cfg.popCap, N, cfg.popDraws);
    end
end


% =====================================================================
function rt = routing_check(S, arms, aList, cfg, repList)
%ROUTING_CHECK  At nPop = 80 and 2000 draws only (not in the paper's run): the values of the
%   first replicates against those stored in cfg.refFile, printed as the largest difference.
    rt = struct('applicable', false);
    if cfg.nPop ~= 80 || ~exist(cfg.refFile,'file') || cfg.betaDraws ~= 2000
        fprintf('\n---- CHECK AGAINST THE STORED RUN: not applicable at N = %d, %d draws ----\n', cfg.nPop, cfg.betaDraws);
        return
    end
    LR = load(cfg.refFile); C = LR.out;
    n0 = C.nReps;
    cols = find(repList <= n0);
    if isempty(cols), return; end
    fprintf('\n---- CHECK AGAINST THE STORED RUN: replicates 1..%d against %s ----\n', n0, cfg.refFile);
    fprintf('Same genotype, same window, same seed, same draw count, so these must be 0.\n\n');
    mets = {'L','L0','w_par','err_rms','beta_win_rms','beta_win_mean', ...
            'beta_ident','L_ident','L_identLin','beta_lin_rms','alpha_op','xnorm_op'};
    fprintf('%-8s %-16s %14s\n','arm','metric','max abs diff');
    for ai = 1:numel(arms)
        for m = 1:numel(mets)
            a = S.(arms{ai}).(mets{m})(:,cols);
            b = C.S.(arms{ai}).(mets{m})(:,repList(cols));
            d = max(abs(a(:) - b(:)));
            rt.(arms{ai}).(mets{m}) = d;
            fprintf('%-8s %-16s %14.3g\n', arms{ai}, mets{m}, d);
        end
    end
    rt.applicable = true;
end


% =====================================================================
function z = one_cell(Th, popTh, popFreq, nPop, Uu, Ou, wt, P, seedBase, f, r, nD, nP, popCap, doPop)
%ONE_CELL  One replicate of one design and arm: L, L0 and the identity's predictions on the
%   best genotype, a check of steady_batch, and the population read. No evolution is run.
    z = struct(); z.nEval = size(Uu,2);

    rng(seedBase + 1000*f + r, 'twister');
    TA = draw_mutants(Th, P, nD);
    m = load_of(Th, TA, Uu, Ou, wt, P);
    z.L = m.L; z.L0 = m.L0; z.w_par = m.w_par; z.err_rms = m.err_rms;
    z.beta_win_rms = m.beta_rms; z.beta_win_mean = m.beta_mean;
    z.L_first = P.selection_strength * m.beta_rms^2;
    z.dLpred = m.L0 - m.L;

    % Check: steady_batch against tracking.steady on the first 25 draws (largest difference)
    ns = min(25, nD); Xb = steady_batch(TA(:,:,1:ns), Uu(:,1), P);
    e = 0;
    for k = 1:ns, e = max(e, max(abs(tracking.steady(TA(:,:,k), Uu(:,1), P) - Xb(:,k)))); end
    z.batchErr = e;

    g = identity_beta(Th, Uu, wt, P);
    z.beta_ident = g.beta; z.L_ident = P.selection_strength * g.beta^2;
    z.alpha_op = g.alpha; z.xnorm_op = g.xnorm;

    ll = lin_load(TA, Th, g, wt, P);
    z.L_identLin = ll.L; z.beta_lin_rms = ll.beta_rms;

    % ---- population read: popCap individuals drawn at random without replacement (so a wiring
    % enters in proportion to its carriers), or every distinct wiring when popCap >= N; the
    % loads are averaged over them weighted by frequency ----
    z.Lpop = NaN; z.L0pop = NaN; z.w_parpop = NaN; z.beta_pop_rms = NaN;
    z.nDistinct = size(popTh,1);
    if doPop
        rng(seedBase + 500000 + 1000*f + r, 'twister');
        if popCap >= nPop
            Uw = popTh; fr = popFreq(:)';           % every distinct wiring
        else
            cnt = round(popFreq(:) * nPop);
            ind = repelem((1:size(popTh,1))', cnt); % expand to individuals
            pick = ind(randperm(numel(ind), min(popCap, numel(ind))));
            u = unique(pick, 'stable');
            Uw = popTh(u,:);
            fr = accumarray(arrayfun(@(k) find(u==k,1), pick), 1)';
            fr = fr / sum(fr);
        end
        nd = size(Uw,1);
        accW = 0; accW0 = 0; accP = 0; accB2 = 0;
        for i = 1:nd
            Ti = reshape(Uw(i,:), P.total, P.total);
            Tm = draw_mutants(Ti, P, nP);
            mi = load_of(Ti, Tm, Uu, Ou, wt, P);
            accW  = accW  + fr(i) * mi.Ew;
            accW0 = accW0 + fr(i) * mi.Ew0;
            accP  = accP  + fr(i) * mi.w_par;
            accB2 = accB2 + fr(i) * mi.beta_rms^2;
        end
        z.Lpop = 1 - accW/accP; z.L0pop = 1 - accW0;
        z.w_parpop = accP; z.beta_pop_rms = sqrt(accB2);
    end
end


% =====================================================================
function m = load_of(Th, TA, Uu, Ou, wt, P)
%LOAD_OF  The load of one parent genotype under the draws TA, over the scoring window (points
%   Uu, optima Ou, weights wt): L against the optimum, L0 against the parent's own output, and
%   the RMS and mean output shift per draw.
    M = size(TA,3); s = P.selection_strength; nE = size(Uu,2);
    accW = zeros(M,1); accW0 = zeros(M,1); accP = 0; accE = 0;
    accD2 = zeros(M,1); accD = zeros(M,1);
    for g = 1:nE
        xp = tracking.steady(Th, Uu(:,g), P); pp = xp([3 4]);
        ep = pp - Ou(:,g);
        accP = accP + wt(g) * exp(-s * sum(ep.^2));
        accE = accE + wt(g) * sum(ep.^2);

        X = steady_batch(TA, Uu(:,g), P);
        dOpt = X([3 4],:) - Ou(:,g);
        dPar = X([3 4],:) - pp;
        accW  = accW  + wt(g) * exp(-s * sum(dOpt.^2,1))';
        accW0 = accW0 + wt(g) * exp(-s * sum(dPar.^2,1))';
        accD2 = accD2 + wt(g) * sum(dPar.^2,1)';
        accD  = accD  + wt(g) * sqrt(sum(dPar.^2,1))';
    end
    m.w_par = accP; m.err_rms = sqrt(accE);
    m.Ew = mean(accW); m.Ew0 = mean(accW0);
    m.L = 1 - m.Ew / m.w_par;
    m.L0 = 1 - m.Ew0;
    m.beta_rms = sqrt(mean(accD2));
    m.beta_mean = mean(accD);
end


% =====================================================================
function g = identity_beta(Th, Uu, wt, P)
%IDENTITY_BETA  The identity's prediction of the RMS output shift under the mutation
%   distribution, from C Phi and x* at each window point (no draws; W33 and W44 not mutated),
%   with alpha and ||x*|| averaged over the window; C Phi and x* are kept for lin_load.
    n = P.total; C = [0 0 1 0; 0 0 0 1];
    mutable = true(n); mutable(3,3) = false; mutable(4,4) = false;
    acc = 0; aacc = 0; xacc = 0;
    nE = size(Uu,2); CPs = cell(1,nE); xs = cell(1,nE);
    for g = 1:nE
        x = steady_fp(Th, Uu(:,g), P);
        h = Th*x + Uu(:,g); act = h > 0;
        fp = zeros(n,1); fp(act) = P.K ./ (P.K + h(act)).^2;
        Df = diag(fp);
        Phi = (P.mu*eye(n) - Df*Th) \ Df;
        CP = C*Phi;
        ca = sum(CP.^2, 1)'; xb = x.^2;
        a2 = 0;
        for k = 1:n, a2 = a2 + ca(k) * sum(xb(mutable(k,:))); end
        acc  = acc  + wt(g) * P.mutation_prob * P.mutation_sd^2 * a2;
        aacc = aacc + wt(g) * norm(CP,'fro');
        xacc = xacc + wt(g) * norm(x);
        CPs{g} = CP; xs{g} = x; %#ok<AGROW>
    end
    g = struct('beta', sqrt(acc), 'alpha', aacc, 'xnorm', xacc, 'CP', {CPs}, 'x', {xs});
end


% =====================================================================
function r = lin_load(TA, Th, g, wt, P)
%LIN_LOAD  The load the identity predicts: each draw's first-order shift C Phi dW x*, scored by
%   the task's fitness over the window.
    M = size(TA,3); s = P.selection_strength; nE = numel(g.CP);
    accW = zeros(M,1); accD2 = zeros(M,1);
    for e = 1:nE
        CP = g.CP{e}; xs = g.x{e};
        d2 = zeros(M,1);
        for k = 1:M
            dp = CP * ((TA(:,:,k) - Th) * xs);
            d2(k) = sum(dp.^2);
        end
        accW  = accW  + wt(e) * exp(-s * d2);
        accD2 = accD2 + wt(e) * d2;
    end
    r.L = 1 - mean(accW);
    r.beta_rms = sqrt(mean(accD2));
end


% =====================================================================
function report_cells(S, arms, aList, nR, cfg)
    fprintf('\n\n=================== PER CELL, over %d replicates ===================\n', nR);
    show = {'L','LOAD L = 1 - E[w(mutant)]/w(parent), best genotype, true optimum'; ...
            'L0','L0: the same draws, optimum re-centred on the parent (pure sensitivity)'; ...
            'w_par','w(parent): the unmutated genotype''s own fitness (TASK DIFFICULTY)'; ...
            'err_rms','RMS tracking error of the unmutated genotype over the window'; ...
            'beta_win_rms','beta_rms: RMS readout shift per mutation, at the operating point'; ...
            'beta_win_mean','beta_mean: MEAN readout shift per mutation'; ...
            'beta_ident','beta_ident: the closed form in alpha and x*, no draws'; ...
            'L_ident','L_ident = s*beta_ident^2'; ...
            'L_identLin','L from the identity''s first-order shift, scored by the REAL fitness'; ...
            'alpha_op','alpha_all at the operating point'; ...
            'xnorm_op','||x*|| at the operating point'; ...
            'Lpop','LOAD, population read'; ...
            'L0pop','L0, population read'; ...
            'nDistinct','distinct wirings among the final individuals'};
    for s = 1:size(show,1)
        v0 = S.(arms{1}).(show{s,1});
        if all(isnan(v0(:))), continue; end
        fprintf('\n---- %s ----\n', show{s,2});
        fprintf('%-8s %-6s %12s %12s %28s\n','arm','a','mean','median','bootstrap 95% CI of mean');
        for ai = 1:numel(arms)
            for f = 1:numel(aList)
                v = S.(arms{ai}).(show{s,1})(f,:);
                ci = boot_mean_ci(v, cfg.nBoot, cfg.rngSeed + 13*s + f + 100*ai);
                fprintf('%-8s %-6.2f %12.5g %12.5g    [%10.5g,%10.5g]\n', ...
                    arms{ai}, aList(f), mean(v,'omitnan'), median(v,'omitnan'), ci(1), ci(2));
            end
        end
    end
end


% =====================================================================
function pd = report_paired(S, arms, aList, nR, cfg)
    fprintf('\n\n============ PAIRED DIFFERENCE, TARGET HELD STILL MINUS MOVING (n = %d) ============\n', nR);
    fprintf('Replicate r shares its founder population across arms, so these are within-replicate.\n');
    fprintf('CI = percentile bootstrap of the mean paired difference, %d draws.\n\n', cfg.nBoot);
    mets = {'L','L0','w_par','err_rms','beta_win_rms','beta_win_mean','beta_ident','L_ident', ...
            'L_identLin','beta_lin_rms','L_first','alpha_op','xnorm_op','Lpop','L0pop'};
    pd = struct('metric',{},'a',{},'diffs',{},'meanDiff',{},'medDiff',{},'ci',{},'pSign',{},'pSRank',{},'ratio',{},'nDown',{});
    fprintf('%-15s %-6s %12s %26s %9s %9s %9s %8s\n', ...
        'metric','a','mean diff','bootstrap 95% CI','fix/mov','p(sign)','p(srank)','n down');
    for m = 1:numel(mets)
        mm = mets{m};
        if all(isnan(S.moving.(mm)(:))), continue; end
        for f = 1:numel(aList)
            vMoving = S.moving.(mm)(f,:); vHeld = S.fixed.(mm)(f,:);
            dv = vHeld - vMoving;
            ci = boot_mean_ci(dv, cfg.nBoot, cfg.rngSeed + 17*m + f);
            ps = sign_test(dv); pr = srank_test(dv);
            rr = median(vHeld,'omitnan') / median(vMoving,'omitnan');
            k = numel(pd)+1;
            pd(k).metric = mm; pd(k).a = aList(f); pd(k).diffs = dv;
            pd(k).meanDiff = mean(dv,'omitnan'); pd(k).medDiff = median(dv,'omitnan');
            pd(k).ci = ci; pd(k).pSign = ps; pd(k).pSRank = pr; pd(k).ratio = rr;
            pd(k).nDown = sum(dv < 0);
            fprintf('%-15s %-6.2f %12.5g    [%10.5g,%10.5g] %9.4f %9.4f %9.4f %5d/%d\n', ...
                mm, aList(f), mean(dv,'omitnan'), ci(1), ci(2), rr, ps, pr, sum(dv<0), sum(~isnan(dv)));
        end
        fprintf('\n');
    end

    fprintf('---- SIGN AGREEMENT: does held-still-minus-moving point the same way for L and beta? ----\n');
    fprintf('%-6s %12s %12s %12s %12s\n','a','d L','d beta_rms','d beta_mean','same sign (rms|mean)');
    for f = 1:numel(aList)
        dL  = mean(S.fixed.L(f,:)             - S.moving.L(f,:));
        dBr = mean(S.fixed.beta_win_rms(f,:)  - S.moving.beta_win_rms(f,:));
        dBm = mean(S.fixed.beta_win_mean(f,:) - S.moving.beta_win_mean(f,:));
        fprintf('%-6.2f %12.5g %12.5g %12.5g      %s | %s\n', aList(f), dL, dBr, dBm, ...
            tf(sign(dL)==sign(dBr)), tf(sign(dL)==sign(dBm)));
    end
end


% =====================================================================
function rel = report_relation(S, A, D, arms, aList, P, repList)
    fprintf('\n\n============ THE EIGHT CELLS: L AGAINST beta ============\n');
    s = P.selection_strength;
    nC = numel(arms)*numel(aList);
    lab = cell(nC,1); Lm = nan(nC,1); L0m = Lm; Lpm = Lm; Br = Lm; Bm = Lm;
    Bi = Lm; Li = Lm; Lil = Lm; Bl = Lm; Wp = Lm;
    bSaved = Lm; bOp = Lm; bAlpha = Lm; bAlphaRms = Lm;
    LmCI = nan(nC,2); LilCI = nan(nC,2);
    k = 0;
    for ai = 1:numel(arms)
        for f = 1:numel(aList)
            k = k+1; R = S.(arms{ai});
            lab{k} = sprintf('%-6s a=%.2f', arms{ai}, aList(f));
            Lm(k) = median(R.L(f,:)); L0m(k) = median(R.L0(f,:));
            Lpm(k) = median(R.Lpop(f,:),'omitnan');
            Br(k) = median(R.beta_win_rms(f,:)); Bm(k) = median(R.beta_win_mean(f,:));
            Bi(k) = median(R.beta_ident(f,:)); Li(k) = median(R.L_ident(f,:));
            Lil(k) = median(R.L_identLin(f,:)); Bl(k) = median(R.beta_lin_rms(f,:));
            Wp(k) = median(R.w_par(f,:));
            LmCI(k,:) = boot_mean_ci(R.L(f,:), 10000, 991 + k);
            LilCI(k,:) = boot_mean_ci(R.L_identLin(f,:), 10000, 4441 + k);
            bSaved(k) = median(D.(arms{ai}).beta(f,repList));
            bOp(k)    = median(D.(arms{ai}).betaOp(f,repList));
            if ~isempty(A)
                bAlpha(k)    = median(A.S.(arms{ai}).beta_meas(f,:));
                bAlphaRms(k) = median(A.S.(arms{ai}).beta_measFP_rms(f,:));
            end
        end
    end
    fprintf('\n%-16s %9s %9s %9s %10s %10s %10s %10s %11s %8s %9s %10s\n', ...
        'cell','L','L0','L(pop)','beta_rms','beta_mean','beta_ident','L_ident','L_identLin','w_par','run beta','run betaOp');
    for k = 1:nC
        fprintf('%-16s %9.4f %9.4f %9.4f %10.4f %10.4f %10.4f %10.4f %11.4f %8.4f %9.4f %10.4f\n', ...
            lab{k}, Lm(k), L0m(k), Lpm(k), Br(k), Bm(k), Bi(k), Li(k), Lil(k), Wp(k), bSaved(k), bOp(k));
    end

    % The two L columns are means over replicates with bootstrap intervals; the ratio built from
    % Lm and Lil is a ratio of medians, labeled as one, and the ratio of the two means (the
    % ratio Table S13 reports) is printed beside it.
    fprintf('\n---- for the figure: measured L against identity-predicted L, with CIs ----\n');
    fprintf('%-16s %10s %20s %12s %20s %14s %12s\n', ...
        'cell','L (mean)','L 95% CI','L_identLin','L_identLin 95% CI', ...
        'L/L_idLin (medians)','(ratio of means)');
    for k = 1:nC
        fprintf('%-16s %10.4f  [%8.4f,%8.4f] %12.4f  [%8.4f,%8.4f] %14.4f %12.4f\n', ...
            lab{k}, mean_of(S, arms, aList, k, 'L'), LmCI(k,1), LmCI(k,2), ...
            mean_of(S, arms, aList, k, 'L_identLin'), LilCI(k,1), LilCI(k,2), ...
            Lm(k)/max(Lil(k),eps), ...
            mean_of(S, arms, aList, k, 'L') / ...
                max(mean_of(S, arms, aList, k, 'L_identLin'), eps));
    end

    fprintf('\n---- Spearman of the cell medians (n = %d cells) ----\n', nC);
    pairs = {'L',Lm,'beta_rms (this file, operating point)',Br; ...
             'L',Lm,'beta_mean (this file, operating point)',Bm; ...
             'L',Lm,'saved beta (run_tracking_task, probe input, 300 draws)',bSaved; ...
             'L',Lm,'betaOp (arm operating point)',bOp; ...
             'L',Lm,'alpha-file beta_meas (probe input)',bAlpha; ...
             'L',Lm,'alpha-file beta_measFP_rms (RMS)',bAlphaRms; ...
             'L',Lm,'beta_ident (closed form, no draws)',Bi; ...
             'L',Lm,'L_identLin (identity''s shift, real fitness)',Lil; ...
             'L0',L0m,'beta_rms',Br; ...
             'L0',L0m,'beta_ident',Bi; ...
             'L0',L0m,'L_identLin',Lil; ...
             'L(pop)',Lpm,'beta_rms',Br};
    fprintf('%-8s %-45s %10s %10s\n','y','x','Spearman','p');
    sp = struct('y',{},'x',{},'rho',{},'p',{});
    for i = 1:size(pairs,1)
        x = pairs{i,4}; y = pairs{i,2};
        if all(isnan(x)) || all(isnan(y)), continue; end
        [rh, pv] = spearman_rp(y, x);
        fprintf('%-8s %-45s %10.3f %10.4f\n', pairs{i,1}, pairs{i,3}, rh, pv);
        j = numel(sp)+1; sp(j).y = pairs{i,1}; sp(j).x = pairs{i,3}; sp(j).rho = rh; sp(j).p = pv;
    end

    fprintf('\n---- the identity''s predicted FORM: L0 = s*beta_rms^2 with s = %.1f ----\n', s);
    kfit = sum(Br.^2 .* L0m) / sum(Br.^4);
    kfitL = sum(Br.^2 .* Lm) / sum(Br.^4);
    kfitI = sum(Bi.^2 .* L0m) / sum(Bi.^4);
    [me, ~] = loglog_fit(Br, L0m); me2 = loglog_fit(Br, Lm);
    fprintf('least squares through the origin on beta_rms:     k = %.2f (predicted %.1f)\n', kfit, s);
    fprintf('the same on L:                                    k = %.2f\n', kfitL);
    fprintf('on the draw-free beta_ident:                      k = %.2f\n', kfitI);
    fprintf('log-log slope of L0 on beta_rms: %.3f (predicted exponent 2); of L: %.3f\n', me, me2);
    fprintf('\n%-16s %10s %10s %12s %10s %11s %11s %11s\n', ...
        'cell','L0','s*beta^2','L0/s*beta^2','L_ident','L0/L_ident','L_identLin','L0/L_idLin');
    for k = 1:nC
        fprintf('%-16s %10.4f %10.4f %12.4f %10.4f %11.4f %11.4f %11.4f\n', lab{k}, L0m(k), ...
            s*Br(k)^2, L0m(k)/(s*Br(k)^2), Li(k), L0m(k)/max(Li(k),eps), Lil(k), L0m(k)/max(Lil(k),eps));
    end
    fprintf('\nbeta_ident against beta_lin_rms (the same closed form evaluated on the draws):\n');
    for k = 1:nC
        fprintf('%-16s  beta_ident = %8.4f   beta_lin_rms = %8.4f   ratio = %6.4f\n', ...
            lab{k}, Bi(k), Bl(k), Bl(k)/max(Bi(k),eps));
    end

    rel = struct('lab',{lab},'L',Lm,'L0',L0m,'Lpop',Lpm,'beta_rms',Br,'beta_mean',Bm, ...
                 'beta_ident',Bi,'L_ident',Li,'L_identLin',Lil,'beta_lin_rms',Bl,'w_par',Wp, ...
                 'L_ci',LmCI,'L_identLin_ci',LilCI,'betaX2',bSaved,'betaX2op',bOp, ...
                 'spearman',sp,'k_L0_on_beta2',kfit,'k_L_on_beta2',kfitL, ...
                 'k_L0_on_betaIdent2',kfitI,'loglog_slope_L0',me,'loglog_slope_L',me2,'s',s);
end


function v = mean_of(S, arms, aList, k, met)
    nF = numel(aList);
    ai = floor((k-1)/nF) + 1; f = mod(k-1, nF) + 1;
    v = mean(S.(arms{ai}).(met)(f,:), 'omitnan');
end


% =====================================================================
function report_confound(S, arms, aList)
    fprintf('\n\n============ LOAD AGAINST TASK DIFFICULTY ============\n');
    fprintf('%-8s %-6s %10s %10s %10s %10s %10s %12s\n', ...
        'arm','a','L','L0','L0 - L','w_par','err_rms','(L0-L)/L0');
    nC = numel(arms)*numel(aList); Lv = nan(nC,1); L0v = Lv; Wv = Lv; k = 0;
    for ai = 1:numel(arms)
        for f = 1:numel(aList)
            k = k+1; R = S.(arms{ai});
            Lv(k) = median(R.L(f,:)); L0v(k) = median(R.L0(f,:)); Wv(k) = median(R.w_par(f,:));
            fprintf('%-8s %-6.2f %10.4f %10.4f %10.4f %10.4f %10.4f %12.4f\n', arms{ai}, aList(f), ...
                Lv(k), L0v(k), L0v(k)-Lv(k), Wv(k), median(R.err_rms(f,:)), ...
                (L0v(k)-Lv(k))/max(abs(L0v(k)),eps));
        end
    end
    [r1,p1] = spearman_rp(Lv, Wv); [r2,p2] = spearman_rp(L0v, Wv);
    [r3,p3] = spearman_rp(L0v - Lv, Wv);
    fprintf('\nSpearman(L,  w_par)   = %+.3f (p = %.4f)\n', r1, p1);
    fprintf('Spearman(L0, w_par)   = %+.3f (p = %.4f)\n', r2, p2);
    fprintf('Spearman(L0-L, w_par) = %+.3f (p = %.4f)\n', r3, p3);
end


% =====================================================================
function s = tf(b)
    if b, s = 'yes'; else, s = 'NO'; end
end

function [rh, pv] = spearman_rp(y, x)
    ok = ~isnan(x) & ~isnan(y); x = x(ok); y = y(ok);
    try
        [rh, pv] = corr(x(:), y(:), 'type', 'Spearman');
    catch
        rh = NaN; pv = NaN;
    end
end

function [m, c] = loglog_fit(x, y)
    ok = x > 0 & y > 0 & ~isnan(x) & ~isnan(y);
    if sum(ok) < 3, m = NaN; c = NaN; return; end
    b = [log(x(ok)) ones(sum(ok),1)] \ log(y(ok));
    m = b(1); c = b(2);
end

function start_pool(n)
    try
        pl = gcp('nocreate');
        if ~isempty(pl) && pl.NumWorkers == n, return; end
        if ~isempty(pl), delete(pl); end
        parpool('Processes', n);
    catch ME
        fprintf('parpool failed (%s) -- running serial\n', ME.message);
    end
end

function diary_file(f, on)
    try
        if on
            d = fileparts(f); if ~exist(d,'dir'), mkdir(d); end
            if exist(f,'file'), delete(f); end
            diary(f); diary on;
        else
            diary off;
        end
    catch
    end
end
