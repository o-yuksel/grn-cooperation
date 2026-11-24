function out = tracking_responsiveness(cfg)
% TRACKING_RESPONSIVENESS  In the tracking task, on the best genotypes at generations 1 and 350 of
% each design, with the target moving and held still: responsiveness alpha, the share of it the
% cue channels (inputs 1-2) carry, ||x*||, and the identity's shift alpha ||x*||.
%
% Produces : Fig. 6 (through tracking_load_vs_identity, which reads the output); Table S13 (||x*|| of the generation-350
%            best genotype; the cue channels' share of alpha at generation 1)
% Inputs   : publication/tracking_runs/tracking_runs_N<nPop>.mat (the stored
%              genotypes; written by run_tracking_task)
%            at nPop = 80 only: publication/tracking_runs/tracking_runs_N80_alpha.mat (part
%              of the data archive), compared with the first replicates (`routing_check`)
%            cfg (optional): nPop (1000; 80 is the N = 80 run), betaDraws (4000), nBoot (10000),
%              rngSeed (20260917), nWorkers (5), repList (all), doPop (true), inFile, outFile,
%              txtFile, refFile
% Outputs  : publication/tracking_runs/tracking_runs_N<nPop>_alpha.mat (read by
%              tracking_load_vs_identity) and _alpha.txt (the console text)
% Run      : tracking_responsiveness()   (about 35 min on 5 workers)
% Seed     : draws on replicate r's final best genotype in design f (its index in aList):
%            rng(20260917 + 1000 f + r, 'twister'); on its generation-1 best:
%            rng(20260917 + 31 (1000 f + r), 'twister'); the bootstrap of metric m over
%            replicates: seed 20260917 + 17 m + f
%
% Every quantity is read at the probe input u = [0, 0.25, 0, 0], at a settled state x* solved
% until the squared residual of a step is below 1e-18. With Phi = (mu I - D_f W)^-1 D_f (D_f the
% gains of the genes that are on) and C selecting outputs 3 and 4: alpha = ||C Phi||_F over all
% four inputs, the cue share ||C Phi(:, 1:2)||_F / alpha, and beta_pred = alpha ||x*||. beta_meas
% is the measured mean shift per mutation under tracking.steady's stopping rule. Population
% medians are over all N final individuals.


    if nargin < 1 || isempty(cfg), cfg = struct(); end
    dd = fullfile(grn_data_file('publication'),'tracking_runs');
    if ~isfield(cfg,'nPop'),      cfg.nPop      = 1000; end
    if ~isfield(cfg,'inFile'),    cfg.inFile    = fullfile(dd, sprintf('tracking_runs_N%d.mat', cfg.nPop)); end
    if ~isfield(cfg,'outFile'),   cfg.outFile   = fullfile(dd, sprintf('tracking_runs_N%d_alpha.mat', cfg.nPop)); end
    if ~isfield(cfg,'txtFile'),   cfg.txtFile   = fullfile(dd, sprintf('tracking_runs_N%d_alpha.txt', cfg.nPop)); end
    if ~isfield(cfg,'refFile'),   cfg.refFile   = fullfile(grn_data_file('publication'),'tracking_runs','tracking_runs_N80_alpha.mat'); end
    if ~isfield(cfg,'betaDraws'), cfg.betaDraws = 4000; end   % draws per measured shift per mutation
    if ~isfield(cfg,'nBoot'),     cfg.nBoot     = 10000; end
    if ~isfield(cfg,'rngSeed'),   cfg.rngSeed   = 20260917; end
    if ~isfield(cfg,'nWorkers'),  cfg.nWorkers  = 5; end
    if ~isfield(cfg,'repList'),   cfg.repList   = []; end
    if ~isfield(cfg,'doPop'),     cfg.doPop     = true; end

    L = load(cfg.inFile); D = L.out;
    P = D.P; aList = D.aList; nA = numel(aList);
    repList = cfg.repList; if isempty(repList), repList = 1:D.nReps; end
    nR = numel(repList);
    % Stored arm names: 'moving' (target moving) and 'fixed' (target held still)
    arms = {'moving','fixed'};
    uProbe = [0; P.trajectory_radius; 0; 0];

    diary_file(cfg.txtFile, true);
    fprintf('\n=========== TRACKING TASK, ALPHA DECOMPOSITION, N = %d, %d replicates ===========\n', P.n_pop, nR);
    fprintf('source %s (run finished %s)\n', cfg.inFile, D.finishedAt);
    fprintf('N = %d, %d generations, mutation p = %.2f at sd = %.3f\n', ...
        P.n_pop, P.n_generations, P.mutation_prob, P.mutation_sd);
    fprintf('per-entry p * 16 entries * N = %.1f (nominal; W33 and W44 are not mutated)\n', 16*P.mutation_prob*P.n_pop);
    fprintf('beta_meas: the mutation distribution of tracking.fragility at %d draws, probe u = [0, %.2f, 0, 0].\n', ...
        cfg.betaDraws, P.trajectory_radius);
    fprintf('alpha definition: ||C Phi||_F over all four inputs (see the header).\n');

    start_pool(cfg.nWorkers);

    nJ = numel(arms)*nA*nR;
    jobs = cell(nJ,1); ja = zeros(nJ,1); jf = ja; jr = ja;
    q = 0;
    for ai = 1:numel(arms)
        for f = 1:nA
            for r = 1:nR, q = q+1; ja(q) = ai; jf(q) = f; jr(q) = r; end
        end
    end
    % slice the genotypes so the parfor does not ship the whole run to every worker
    ThB = cell(nJ,1); ThG = cell(nJ,1); PopT = cell(nJ,1); PopF = cell(nJ,1);
    for q = 1:nJ
        R = D.(arms{ja(q)}); rr = repList(jr(q));
        ThB{q} = R.ThBest{jf(q),rr}; ThG{q} = R.ThG1{jf(q),rr};
        PopT{q} = R.popTh{jf(q),rr}; PopF{q} = R.popFreq{jf(q),rr};
    end
    seedBase = cfg.rngSeed; nDraws = cfg.betaDraws; doPop = cfg.doPop;
    rl = repList; npop = P.n_pop;
    t0 = tic;
    parfor q = 1:nJ
        jobs{q} = one_cell(ThB{q}, ThG{q}, PopT{q}, PopF{q}, npop, P, uProbe, ...
                           seedBase, jf(q), rl(jr(q)), nDraws, doPop);
    end
    fprintf('\n%d cells measured in %.1f min\n', nJ, toc(t0)/60);

    mets = {'xnorm','alpha_all','alpha_cue','alpha_nocue','alpha_share', ...
            'beta_pred','beta_pred_rms','beta_meas','beta_meas_rms','beta_measFP', ...
            'beta_measFP_rms','beta_lin_rms','lam_min','fp_min','fp_out','nactive', ...
            'ratio_fs','ratio_q90','ratio_q99','p_tail3','q99_dp','hmin','resid','maxRe','fd_rel', ...
            'xnorm0','alpha_all0','alpha_cue0','beta_pred0', ...
            'xnorm_pop','alpha_all_pop','alpha_cue_pop','alpha_nocue_pop','beta_pred_pop'};
    S = struct();
    for ai = 1:numel(arms)
        for m = 1:numel(mets), S.(arms{ai}).(mets{m}) = nan(nA,nR); end
    end
    G1 = struct();
    g1m = {'xnorm','alpha_all','alpha_cue','alpha_nocue','alpha_share','beta_pred','beta_meas','fp_out','nactive'};
    for m = g1m, G1.(m{1}) = nan(nA,nR); end
    for q = 1:nJ
        z = jobs{q}; ar = arms{ja(q)}; f = jf(q); r = jr(q);
        for m = 1:numel(mets), S.(ar).(mets{m})(f,r) = z.(mets{m}); end
        if ja(q) == 1
            for m = g1m, G1.(m{1})(f,r) = z.(['g1_' m{1}]); end
        end
    end

    out = struct('cfg',cfg,'aList',aList,'nPop',P.n_pop,'nReps',nR,'repList',repList, ...
                 'arms',{arms},'S',S,'G1',G1,'sourceFinished',D.finishedAt);

    out.routing = routing_check(S, arms, aList, cfg, repList);
    report_cells(S, G1, arms, aList, nR, cfg);
    out.paired = report_paired(S, arms, aList, nR, cfg);
    report_decomposition(S, aList, nR);
    report_identity(S, arms, aList, nR);

    out.finishedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    save(cfg.outFile, 'out');
    fprintf('\nSaved %s\n', cfg.outFile);
    diary_file(cfg.txtFile, false);
end


% =====================================================================
function rt = routing_check(S, arms, aList, cfg, repList)
%ROUTING_CHECK  At nPop = 80 only (not in the paper's run): the values of the first replicates
%   against those stored in cfg.refFile, printed as the largest absolute difference.
    rt = struct('applicable', false);
    if cfg.nPop ~= 80 || ~exist(cfg.refFile,'file')
        fprintf('\n---- CHECK AGAINST THE STORED RUN: not applicable at N = %d ----\n', cfg.nPop);
        return
    end
    LA = load(cfg.refFile); A = LA.out;
    n0 = A.nReps;
    cols = find(repList <= n0);
    if isempty(cols), return; end
    fprintf('\n---- CHECK AGAINST THE STORED RUN: replicates 1..%d against %s ----\n', n0, cfg.refFile);
    fprintf('The geometry is deterministic and the draw seed is the same, so these must be 0.\n\n');
    mets = {'xnorm','alpha_all','alpha_cue','alpha_nocue','beta_pred','beta_pred_rms', ...
            'beta_meas','beta_measFP','hmin','lam_min','xnorm_pop','alpha_all_pop','alpha_cue_pop'};
    fprintf('%-8s %-16s %14s\n','arm','metric','max abs diff');
    for ai = 1:numel(arms)
        for m = 1:numel(mets)
            a = S.(arms{ai}).(mets{m})(:,cols);
            b = A.S.(arms{ai}).(mets{m})(:,repList(cols));
            d = max(abs(a(:) - b(:)));
            rt.(arms{ai}).(mets{m}) = d;
            fprintf('%-8s %-16s %14.3g\n', arms{ai}, mets{m}, d);
        end
    end
    rt.applicable = true;
end


% =====================================================================
function z = one_cell(Th, ThG1, popTh, popFreq, nPop, P, uProbe, seedBase, f, r, nDraws, doPop)
%ONE_CELL  One replicate of one design and arm: the geometry of the final best genotype at the
%   probe and at u = 0, its measured shift per mutation, the population medians, and the same
%   for the generation-1 best genotype. No evolution is run.
    g = geometry(Th, uProbe, P);
    z = struct();
    z.xnorm = g.xnorm; z.alpha_all = g.alpha_all; z.alpha_cue = g.alpha_cue;
    z.alpha_nocue = g.alpha_nocue; z.alpha_share = g.alpha_share;
    z.beta_pred = g.beta_pred; z.beta_pred_rms = g.beta_pred_rms;
    z.lam_min = g.lam_min; z.fp_min = min(g.fprime); z.fp_out = mean(g.fprime(3:4));
    z.nactive = sum(g.fprime > 0);
    z.hmin = g.hmin; z.resid = g.resid; z.maxRe = g.maxRe;
    z.fd_rel = fd_check(Th, uProbe, P, g);

    rng(seedBase + 1000*f + r, 'twister');
    fr = frag_and_lin(Th, uProbe, P, nDraws, g);
    z.beta_meas = fr.mean; z.beta_meas_rms = fr.rms; z.ratio_fs = fr.medRatio;
    z.beta_measFP = fr.meanFP; z.beta_measFP_rms = fr.rmsFP; z.beta_lin_rms = fr.linRms;
    z.ratio_q90 = fr.q90Ratio; z.ratio_q99 = fr.q99Ratio;
    z.p_tail3 = fr.pTail3; z.q99_dp = fr.q99dp;

    g0 = geometry(Th, zeros(P.total,1), P);
    z.xnorm0 = g0.xnorm; z.alpha_all0 = g0.alpha_all;
    z.alpha_cue0 = g0.alpha_cue; z.beta_pred0 = g0.beta_pred;

    % Population median over all nPop final individuals: each distinct wiring is measured once,
    % then repeated by its count (its frequency times nPop) before the median
    if doPop
        nd = size(popTh,1);
        cnt = round(popFreq(:) * nPop);
        xv = nan(1,nd); aa = xv; ac = xv; an = xv; bp = xv;
        for i = 1:nd
            gi = geometry(reshape(popTh(i,:), P.total, P.total), uProbe, P);
            xv(i) = gi.xnorm; aa(i) = gi.alpha_all; ac(i) = gi.alpha_cue;
            an(i) = gi.alpha_nocue; bp(i) = gi.beta_pred;
        end
        z.xnorm_pop      = wmedian(xv, cnt);
        z.alpha_all_pop  = wmedian(aa, cnt);
        z.alpha_cue_pop  = wmedian(ac, cnt);
        z.alpha_nocue_pop= wmedian(an, cnt);
        z.beta_pred_pop  = wmedian(bp, cnt);
    else
        z.xnorm_pop = NaN; z.alpha_all_pop = NaN; z.alpha_cue_pop = NaN;
        z.alpha_nocue_pop = NaN; z.beta_pred_pop = NaN;
    end

    % Generation 1: the best genotype of the founder population, the same in both arms (same
    % seed, same first position of the target)
    g1 = geometry(ThG1, uProbe, P);
    z.g1_xnorm = g1.xnorm; z.g1_alpha_all = g1.alpha_all; z.g1_alpha_cue = g1.alpha_cue;
    z.g1_alpha_nocue = g1.alpha_nocue; z.g1_alpha_share = g1.alpha_share;
    z.g1_beta_pred = g1.beta_pred; z.g1_fp_out = mean(g1.fprime(3:4));
    z.g1_nactive = sum(g1.fprime > 0);
    rng(seedBase + 31*(1000*f + r), 'twister');
    fr1 = frag_and_lin(ThG1, uProbe, P, nDraws, g1);
    z.g1_beta_meas = fr1.mean;
end


% =====================================================================
function m = wmedian(v, cnt)
%WMEDIAN  Median over the population, each distinct wiring repeated cnt times.
    ok = ~isnan(v) & cnt(:)' > 0;
    if ~any(ok), m = NaN; return; end
    m = median(repelem(v(ok), cnt(ok)));
end


% =====================================================================
function g = geometry(Th, u, P)
%GEOMETRY  The congruence identity's quantities at the fixed point x*: D_f, the gains f'(h) of
%   the genes that are on (h = W x* + u > 0); Phi = (mu I - D_f W)^-1 D_f; alpha = ||C Phi||_F,
%   split into the cue channels (inputs 1-2) and the others; beta_pred = alpha ||x*||. Also the
%   Jacobian's eigenvalues, the smallest |h_i| and the residual; beta_pred_rms is the RMS
%   first-order shift under the mutation distribution (W33 and W44 not mutated).
    n = P.total;
    x = steady_fp(Th, u, P);
    h = Th*x + u;
    act = h > 0;
    fp = zeros(n,1);
    fp(act) = P.K ./ (P.K + h(act)).^2;
    Df = diag(fp);
    A  = P.mu*eye(n) - Df*Th;
    Phi = A \ Df;
    C  = [0 0 1 0; 0 0 0 1];
    CP = C*Phi;

    g.x = x; g.xnorm = norm(x);
    g.fprime = fp; g.active = act;
    g.Phi = Phi; g.CP = CP;
    g.alpha_all   = norm(CP, 'fro');
    g.alpha_cue   = norm(CP(:,1:2), 'fro');
    g.alpha_nocue = norm(CP(:,3:4), 'fro');
    g.alpha_share = g.alpha_cue / max(g.alpha_all, eps);
    g.beta_pred = g.alpha_all * g.xnorm;

    J = Df*Th - P.mu*eye(n);
    ev = eig(J);
    g.lam_min = min(abs(ev));
    g.maxRe   = max(real(ev));
    g.hmin    = min(abs(h));
    g.resid   = norm(max(h,0)./(P.K+max(h,0)) - P.mu*x);

    ca = sum(CP.^2, 1)';
    xb = x.^2;
    mutable = true(n); mutable(3,3) = false; mutable(4,4) = false;
    acc = 0;
    for a = 1:n
        acc = acc + ca(a) * sum(xb(mutable(a,:)));
    end
    g.beta_pred_rms = sqrt(P.mutation_prob * P.mutation_sd^2 * acc);
end


% =====================================================================
function fr = frag_and_lin(Th, u, P, ndraws, g)
%FRAG_AND_LIN  Per draw of the mutation distribution, in tracking.fragility's order: the output
%   shift under tracking.steady (mean, rms), the same at the converged fixed point (meanFP,
%   rmsFP), and the identity's first-order shift C Phi dW x* (linRms); per draw, the ratio of
%   the converged shift to the first-order one (median, 90th and 99th percentiles, share above 3).
    xL = tracking.steady(Th, u, P); pL = xL([3 4]);
    x0 = g.x; p0 = x0([3 4]);
    accL = 0; accL2 = 0; accF = 0; accF2 = 0;
    rat = nan(1, ndraws); dls = nan(1, ndraws); des = nan(1, ndraws);
    CP = g.CP;
    for k = 1:ndraws
        Tm = Th; mask = rand(size(Tm)) < P.mutation_prob;
        Tm(mask) = Tm(mask) + randn(nnz(mask),1) * P.mutation_sd;
        Tm(3,3) = Th(3,3); Tm(4,4) = Th(4,4);
        dW = Tm - Th;
        xm = tracking.steady(Tm, u, P);
        dL = norm(xm([3 4]) - pL);
        accL = accL + dL; accL2 = accL2 + dL^2;

        xf = steady_fp(Tm, u, P);
        dF = norm(xf([3 4]) - p0);
        accF = accF + dF; accF2 = accF2 + dF^2;

        dl = norm(CP * (dW * x0));
        dls(k) = dl; des(k) = dF;
        if dl > 1e-12, rat(k) = dF / dl; end
    end
    fr.mean = accL / ndraws;  fr.rms  = sqrt(accL2 / ndraws);
    fr.meanFP = accF / ndraws; fr.rmsFP = sqrt(accF2 / ndraws);
    fr.linRms = sqrt(mean(dls.^2));
    rr = rat(~isnan(rat));
    fr.medRatio = median(rr);
    fr.q90Ratio = prctile(rr, 90);
    fr.q99Ratio = prctile(rr, 99);
    fr.pTail3   = mean(rr > 3);
    fr.q99dp    = prctile(des(~isnan(des)), 99);
end


% =====================================================================
function e = fd_check(Th, u, P, g)
%FD_CHECK  Central finite difference of the outputs' response to each input against the
%   analytic C*Phi, as a relative Frobenius error.
    d = 1e-6; FD = zeros(2, P.total);
    for i = 1:P.total
        up = u; up(i) = up(i) + d; xp = steady_fp(Th, up, P);
        um = u; um(i) = um(i) - d; xm = steady_fp(Th, um, P);
        FD(:,i) = (xp([3 4]) - xm([3 4])) / (2*d);
    end
    e = norm(g.CP - FD, 'fro') / max(norm(FD,'fro'), eps);
end


% =====================================================================
function report_cells(S, G1, arms, aList, nR, cfg)
    fprintf('\n\n=============== PER CELL, median [IQR] over %d replicates ===============\n', nR);
    show = {'xnorm','||x*||'; 'alpha_all','alpha_all'; 'alpha_cue','alpha_cue'; ...
            'alpha_nocue','alpha_nocue'; 'alpha_share','alpha_cue / alpha_all'; ...
            'beta_pred','beta_pred = alpha_all*||x*||'; 'beta_pred_rms','beta_pred_rms (mutation distribution)'; ...
            'beta_meas','beta_meas (tracking.steady)'; ...
            'beta_measFP','beta_meas (converged fixed point)'; ...
            'lam_min','min|eig J| (bowl depth)'; 'fp_out','f''(h) mean over outputs 3,4'; ...
            'nactive','n genes active'; 'hmin','min|h_i| (0 = on the rectifier kink)'; ...
            'ratio_fs','median per-draw exact/linear ratio'; ...
            'p_tail3','share of draws exceeding 3x the first-order shift'; ...
            'fd_rel','finite-difference cross-check on C*Phi (relative)'};
    for s = 1:size(show,1)
        fprintf('\n---- %s ----\n', show{s,2});
        fprintf('%-8s %-6s %12s %12s %26s\n','arm','a','mean','median','IQR');
        for ai = 1:numel(arms)
            for f = 1:numel(aList)
                v = S.(arms{ai}).(show{s,1})(f,:);
                qq = prctile(v,[25 75]);
                fprintf('%-8s %-6.2f %12.5g %12.5g    [%10.5g,%10.5g]\n', ...
                    arms{ai}, aList(f), mean(v), median(v), qq(1), qq(2));
            end
        end
    end

    fprintf('\n\n---- generation 1 (the SHARED start: identical population in both arms) ----\n');
    fprintf('%-6s %10s %10s %10s %11s %10s %10s %10s %10s %8s\n', ...
        'a','||x*||','alpha_all','alpha_cue','alpha_nocue','cue share','beta_pred','beta_meas','f''(h) out','n active');
    for f = 1:numel(aList)
        fprintf('%-6.2f %10.5g %10.5g %10.5g %11.5g %10.4f %10.5g %10.5g %10.5g %8.2f\n', aList(f), ...
            median(G1.xnorm(f,:)), median(G1.alpha_all(f,:)), median(G1.alpha_cue(f,:)), ...
            median(G1.alpha_nocue(f,:)), median(G1.alpha_share(f,:)), ...
            median(G1.beta_pred(f,:)), median(G1.beta_meas(f,:)), ...
            median(G1.fp_out(f,:)), median(G1.nactive(f,:)));
    end

    if cfg.doPop
        fprintf('\n---- population median over all %d final individuals ----\n', cfg.nPop);
        fprintf('%-8s %-6s %12s %12s %12s %12s\n','arm','a','||x*||','alpha_all','alpha_cue','beta_pred');
        for ai = 1:numel(arms)
            for f = 1:numel(aList)
                R = S.(arms{ai});
                fprintf('%-8s %-6.2f %12.5g %12.5g %12.5g %12.5g\n', arms{ai}, aList(f), ...
                    median(R.xnorm_pop(f,:)), median(R.alpha_all_pop(f,:)), ...
                    median(R.alpha_cue_pop(f,:)), median(R.beta_pred_pop(f,:)));
            end
        end
    end

    fprintf('\n---- alpha_cue exactly zero: how many replicates of %d ----\n', nR);
    fprintf('%-8s %-6s %14s %14s\n','arm','a','best genotype','population median');
    for ai = 1:numel(arms)
        for f = 1:numel(aList)
            R = S.(arms{ai});
            fprintf('%-8s %-6.2f %14d %14d\n', arms{ai}, aList(f), ...
                sum(R.alpha_cue(f,:) < 1e-12), sum(R.alpha_cue_pop(f,:) < 1e-12));
        end
    end
end


% =====================================================================
function pd = report_paired(S, arms, aList, nR, cfg)
    fprintf('\n\n============ PAIRED DIFFERENCE, TARGET HELD STILL MINUS MOVING (n = %d) ============\n', nR);
    fprintf('Replicate r shares its founder population across arms, so these are within-replicate.\n');
    fprintf('CI = percentile bootstrap of the mean paired difference, %d draws.\n', cfg.nBoot);
    fprintf('p(sign) exact two-sided sign test; p(srank) Wilcoxon signed rank.\n\n');
    mets = {'xnorm','alpha_all','alpha_cue','alpha_nocue','alpha_share','beta_pred', ...
            'beta_pred_rms','beta_meas','beta_measFP','lam_min','fp_out','hmin','p_tail3', ...
            'xnorm_pop','alpha_all_pop','alpha_cue_pop','beta_pred_pop'};
    pd = struct('metric',{},'a',{},'diffs',{},'meanDiff',{},'medDiff',{},'ci',{},'pSign',{},'pSRank',{},'ratio',{});
    fprintf('%-15s %-6s %12s %26s %9s %9s %9s\n', ...
        'metric','a','mean diff','bootstrap 95% CI','fix/mov','p(sign)','p(srank)');
    for m = 1:numel(mets)
        mm = mets{m};
        for f = 1:numel(aList)
            vMoving = S.moving.(mm)(f,:); vHeld = S.fixed.(mm)(f,:);
            dv = vHeld - vMoving;
            ci = boot_mean_ci(dv, cfg.nBoot, cfg.rngSeed + 17*m + f);
            ps = sign_test(dv); pr = srank_test(dv);
            rr = median(vHeld) / median(vMoving);
            k = numel(pd)+1;
            pd(k).metric = mm; pd(k).a = aList(f); pd(k).diffs = dv;
            pd(k).meanDiff = mean(dv); pd(k).medDiff = median(dv); pd(k).ci = ci;
            pd(k).pSign = ps; pd(k).pSRank = pr; pd(k).ratio = rr;
            fprintf('%-15s %-6.2f %12.5g    [%10.5g,%10.5g] %9.4f %9.4f %9.4f\n', ...
                mm, aList(f), mean(dv), ci(1), ci(2), rr, ps, pr);
        end
        fprintf('\n');
    end
end


% =====================================================================
function report_decomposition(S, aList, nR)
    M = S.moving; X = S.fixed;
    fprintf('\n\n============ WHERE DID alpha GO? THE QUADRATURE CHECK ============\n');
    fprintf('If the arm with the target held still lost ONLY the cue channels, its alpha_all would equal the moving arm''s\n');
    fprintf('alpha_nocue. Medians over %d replicates.\n\n', nR);
    fprintf('%-6s %11s %11s %11s %9s %11s %11s %11s %9s\n', ...
        'a','mov a_all','mov a_cue','mov a_nocue','cue sh.','fix a_all','fix a_cue','predicted','pred/obs');
    for f = 1:numel(aList)
        am = median(M.alpha_all(f,:)); cm = median(M.alpha_cue(f,:)); nm = median(M.alpha_nocue(f,:));
        af = median(X.alpha_all(f,:)); cf = median(X.alpha_cue(f,:));
        fprintf('%-6.2f %11.5g %11.5g %11.5g %9.4f %11.5g %11.5g %11.5g %9.4f\n', ...
            aList(f), am, cm, nm, cm/am, af, cf, nm, nm/max(af,eps));
    end

    fprintf('\n============ THE BETA RATIO, FACTORED ============\n');
    fprintf('%-6s %11s %11s %11s %13s %12s %12s %16s\n', ...
        'a','alpha_all','||x*||','beta_pred','beta_pred_rms','beta_meas','residual','tail share m|f');
    for f = 1:numel(aList)
        ra = median(X.alpha_all(f,:)) / median(M.alpha_all(f,:));
        rx = median(X.xnorm(f,:))     / median(M.xnorm(f,:));
        rp = median(X.beta_pred(f,:)) / median(M.beta_pred(f,:));
        rk = median(X.beta_pred_rms(f,:)) / median(M.beta_pred_rms(f,:));
        rb = median(X.beta_measFP(f,:))   / median(M.beta_measFP(f,:));
        fprintf('%-6.2f %11.4f %11.4f %11.4f %13.4f %12.4f %12.4f  %6.3f | %6.3f\n', ...
            aList(f), ra, rx, rp, rk, rb, rb/rk, median(M.p_tail3(f,:)), median(X.p_tail3(f,:)));
    end
end


% =====================================================================
function report_identity(S, arms, aList, nR)
    fprintf('\n============ THE IDENTITY: DOES beta = alpha*||x*|| HOLD IN BOTH ARMS? ============\n');
    fprintf('(a) beta_pred_rms against beta_lin_rms: algebra. (b) against beta_measFP_rms:\n');
    fprintf('same mutation distribution. (c) the median per-draw exact/linear ratio.\n\n');
    fprintf('%-8s %-6s %13s %13s %8s %15s %8s %12s %12s\n', ...
        'arm','a','beta_pred_rms','beta_lin_rms','(a)','beta_measFP_rms','(b)','(c) median','(c) min|max');
    for ai = 1:numel(arms)
        for f = 1:numel(aList)
            R = S.(arms{ai});
            pr = median(R.beta_pred_rms(f,:)); li = median(R.beta_lin_rms(f,:));
            me = median(R.beta_measFP_rms(f,:));
            fprintf('%-8s %-6.2f %13.5g %13.5g %8.4f %15.5g %8.4f %12.4f %6.3f|%.3f\n', ...
                arms{ai}, aList(f), pr, li, li/pr, me, me/pr, median(R.ratio_fs(f,:)), ...
                min(R.ratio_fs(f,:)), max(R.ratio_fs(f,:)));
        end
    end
end


% =====================================================================
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
