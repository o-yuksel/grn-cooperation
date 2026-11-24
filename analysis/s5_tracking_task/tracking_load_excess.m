function out = tracking_load_excess(cfg)
% TRACKING_LOAD_EXCESS  With the target held still at a = 0.70, the measured per-mutation load
% exceeds the identity's prediction: how much of the excess comes from the draws that turn a gene
% on or off? Two other cells (a = 0.70 moving, a = 0.28 held still) are measured beside it.
%
% Produces : SI S5 (the share of draws that carry the excess, and measured over predicted load on
%            the draws that leave every gene as it was)
% Inputs   : publication/tracking_runs/tracking_runs_N1000.mat (best genotypes and
%              scoring windows; written by run_tracking_task) and tracking_runs_N1000_load.mat
%              in the same folder (written by tracking_load_vs_identity; to check the draws match)
%            cfg (optional): nPop (1000), nDraws (6000), nTight (500),
%              tolFacs ([1e-2 1e-4 1e-8]), nBoot (10000), rngSeed (20260917), nWorkers (5),
%              cells ({2 1; 1 1; 2 2}: arm, 1 moving or 2 held still, and the design's index in
%              aList, a = 0.70, 0.28, 0.20, 0.15), inFile, loadFile, outFile, txtFile
% Outputs  : publication/tracking_runs/tracking_load_excess.mat and .txt (console)
% Run      : tracking_load_excess()   (about 5 min on 5 workers)
% Seed     : the draws on replicate r's best genotype in design f use rng(20260917 + 1000 f + r,
%            'twister'), the seeds of tracking_load_vs_identity, so the first 2000 draws are its
%            draws; the bootstraps use their own fixed streams
%
% A cell is one arm and one design. A draw turns a gene on or off when the gene's input
% h = W x* + u changes sign between parent and mutant (each at its own settled state) at any point
% of the scoring window; the identity, first order at the parent's active genes, cannot see such a
% draw. Per draw, L_k = 1 - w(mutant)/w(parent), w the fitness over the scoring window, and
% Lpred_k scores the identity's first-order shift the same way; their means over the first 2000
% draws are checked against L and L_identLin of the load file. Also printed: the load with nTight
% draws re-solved at tighter stopping tolerances, and how the excess spreads over replicates and
% draws.


    if nargin < 1 || isempty(cfg), cfg = struct(); end
    dd = fullfile(grn_data_file('publication'),'tracking_runs');
    if ~isfield(cfg,'nPop'),     cfg.nPop     = 1000; end
    if ~isfield(cfg,'inFile'),   cfg.inFile   = fullfile(dd, sprintf('tracking_runs_N%d.mat', cfg.nPop)); end
    if ~isfield(cfg,'loadFile'), cfg.loadFile = fullfile(dd, sprintf('tracking_runs_N%d_load.mat', cfg.nPop)); end
    if ~isfield(cfg,'outFile'),  cfg.outFile  = fullfile(dd, 'tracking_load_excess.mat'); end
    if ~isfield(cfg,'txtFile'),  cfg.txtFile  = fullfile(dd, 'tracking_load_excess.txt'); end
    if ~isfield(cfg,'nDraws'),   cfg.nDraws   = 6000; end     % draws per genotype
    if ~isfield(cfg,'nTight'),   cfg.nTight   = 500; end      % draws re-solved at tighter tolerances
    if ~isfield(cfg,'tolFacs'),  cfg.tolFacs  = [1e-2 1e-4 1e-8]; end  % factors on the tolerance
    if ~isfield(cfg,'nBoot'),    cfg.nBoot    = 10000; end
    if ~isfield(cfg,'rngSeed'),  cfg.rngSeed  = 20260917; end % base seed of the draws
    if ~isfield(cfg,'nWorkers'), cfg.nWorkers = 5; end
    % cells: {arm index (1 moving, 2 held still), design index in aList}
    if ~isfield(cfg,'cells'),    cfg.cells    = {2 1; 1 1; 2 2}; end

    LD = load(cfg.inFile); D = LD.out;
    P = D.P; aList = D.aList; arms = {'moving','fixed'};
    nR = D.nReps; nC = size(cfg.cells,1);

    diary_file(cfg.txtFile, true);
    fprintf('\n=========== TRACKING TASK, THE LOW-LOAD CELL, N = %d, %d replicates ===========\n', P.n_pop, nR);
    fprintf('source %s (run finished %s)\n', cfg.inFile, D.finishedAt);
    fprintf('%d draws per genotype at the seeds of tracking_load_vs_identity; tolerance probe on %d of them\n', ...
        cfg.nDraws, cfg.nTight);
    fprintf('mutation distribution: per-entry p = %.2f at sd = %.3f, self-activations (3,3), (4,4) not mutated; s = %.1f\n\n', ...
        P.mutation_prob, P.mutation_sd, P.selection_strength);

    start_pool(cfg.nWorkers);

    nJ = nC*nR; jc = zeros(nJ,1); jr = jc; q = 0;
    for k = 1:nC
        for r = 1:nR, q = q+1; jc(q) = k; jr(q) = r; end
    end
    Th = cell(nJ,1); Uu = Th; Ou = Th; Wt = Th;
    for q = 1:nJ
        ai = cfg.cells{jc(q),1}; f = cfg.cells{jc(q),2};
        R = D.(arms{ai});
        Th{q} = R.ThBest{f,jr(q)}; Uu{q} = R.winU{f}; Ou{q} = R.winO{f}; Wt{q} = R.winW{f};
    end
    seeds = zeros(nJ,1);
    for q = 1:nJ, seeds(q) = cfg.rngSeed + 1000*cfg.cells{jc(q),2} + jr(q); end
    nD = cfg.nDraws; nT = cfg.nTight; tf = cfg.tolFacs;
    res = cell(nJ,1);
    t0 = tic;
    parfor q = 1:nJ
        res{q} = one_cell(Th{q}, Uu{q}, Ou{q}, Wt{q}, P, seeds(q), nD, nT, tf);
    end
    fprintf('%d cell-replicates measured in %.1f min\n', nJ, toc(t0)/60);

    % ---- pack ----
    flds = fieldnames(res{1});
    C = struct();
    for k = 1:nC
        lab = sprintf('%s_a%02d', arms{cfg.cells{k,1}}, round(100*aList(cfg.cells{k,2})));
        C(k).label = lab;
        C(k).arm = arms{cfg.cells{k,1}}; C(k).a = aList(cfg.cells{k,2});
        for m = 1:numel(flds)
            v = nan(nR, numel(res{1}.(flds{m})));
            for r = 1:nR, v(r,:) = res{(k-1)*nR + r}.(flds{m}); end
            C(k).(flds{m}) = v;
        end
    end

    out = struct('cfg',cfg,'nPop',P.n_pop,'nReps',nR,'cells',{C}, ...
                 'sourceFinished',D.finishedAt);

    routing_check(C, cfg, aList, arms);
    report_sign_changes(C, cfg);
    report_ratios(C, cfg);
    report_skew(C, cfg);
    report_tolerance(C, cfg);
    report_margins(C);

    out.finishedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    save(cfg.outFile, 'out');
    fprintf('\nSaved %s\n', cfg.outFile);
    diary_file(cfg.txtFile, false);
end


% =====================================================================
function z = one_cell(Th, Uu, Ou, wt, P, seed, nD, nT, tolFacs)
%ONE_CELL  One stored genotype and its draws of the mutation distribution, split into the draws
%   that change the set of active genes and those that do not.
    n = P.total; s = P.selection_strength; nE = size(Uu,2);
    rng(seed, 'twister');
    TA = draw_mutants(Th, P, nD);

    % ---- parent at the runs' stopping tolerance; marg, the smallest |h| over genes and window
    % points, is how near a gene sits to turning on or off ----
    actP = false(n, nE); marg = inf; wpar = 0; xp = zeros(n, nE);
    for g = 1:nE
        x = tracking.steady(Th, Uu(:,g), P);
        h = Th*x + Uu(:,g);
        actP(:,g) = h > 0; marg = min(marg, min(abs(h))); xp(:,g) = x;
        wpar = wpar + wt(g) * exp(-s * sum((x([3 4]) - Ou(:,g)).^2));
    end

    % ---- mutants at the runs' stopping tolerance: the genes whose h changes sign, and fitness ----
    accW = zeros(nD,1); sw = false(nD,1); accD2 = zeros(nD,1);
    swGene = zeros(1,n); swOn = 0; swOff = 0; nSwEvt = zeros(nD,1);
    Xprod = cell(1,nE);
    for g = 1:nE
        X = steady_batch(TA, Uu(:,g), P, P.dev_tolerance, P.max_dev_steps);
        Xprod{g} = X([3 4],:);
        H = drive_field(TA, X, Uu(:,g));
        A = H > 0;
        dif = A ~= actP(:,g);
        sw = sw | any(dif,1)';
        nSwEvt = nSwEvt + sum(dif,1)';
        swGene = swGene + sum(dif,2)';
        swOn = swOn + nnz(dif & A); swOff = swOff + nnz(dif & ~A);
        d = X([3 4],:) - Ou(:,g);
        accW = accW + wt(g) * exp(-s * sum(d.^2,1))';
        dp = X([3 4],:) - xp([3 4],g);
        accD2 = accD2 + wt(g) * sum(dp.^2,1)';
    end
    Lk = 1 - accW / wpar;

    % ---- the identity's first-order prediction, per draw ----
    G = identity_beta(Th, Uu, wt, P);
    Pk = lin_load_draws(TA, Th, G, wt, P);

    ns = ~sw;
    z.nD = nD;
    z.nSwitch = sum(sw);
    z.switchFrac = mean(sw);
    z.meanSwEvents = mean(nSwEvt);
    z.swGene = swGene / nD;
    z.swOnFrac = swOn / max(swOn + swOff, 1);
    z.marginParent = marg;
    z.w_par = wpar;
    z.xnorm = norm(xp(:,1));

    % contributions: every sum is divided by nD, so the parts add to the load
    nRoute = min(2000, nD);         % the load file's draw count, for the check of the draws
    z.L_route = mean(Lk(1:nRoute)); z.Lpred_route = mean(Pk(1:nRoute));
    z.L      = sum(Lk) / nD;        z.Lpred    = sum(Pk) / nD;
    z.L_sw   = sum(Lk(sw)) / nD;    z.Lpred_sw = sum(Pk(sw)) / nD;
    z.L_ns   = sum(Lk(ns)) / nD;    z.Lpred_ns = sum(Pk(ns)) / nD;
    z.shareSw = z.L_sw / max(z.L, eps);
    z.Lmed_draw = median(Lk);
    z.top1share = sum(maxk(Lk, max(1, round(0.01*nD)))) / max(sum(Lk), eps);
    z.maxLk = max(Lk);
    z.meanLk_sw = mean(Lk(sw)); z.meanLk_ns = mean(Lk(ns));

    z.Lquad = s * mean(accD2);   % gamma times the mean squared shift: the small-shift form of L

    % ---- tighter tolerances: the parent and the first nTight draws re-solved at the runs'
    % tolerance times each factor, with up to 5 times the step limit ----
    m = min(nT, nD); ms = 5 * P.max_dev_steps; nL = numel(tolFacs);
    z.L_tol = nan(1,nL); z.wpar_tol = nan(1,nL); z.step_tol = nan(1,nL);
    dread = 0; dreadPar = 0;
    for L = 1:nL
        tol = P.dev_tolerance * tolFacs(L);
        wparT = 0; xparT = zeros(2,nE);
        for g = 1:nE
            xt = steady_tight(Th, Uu(:,g), P, tol, ms);
            xparT(:,g) = xt([3 4]);
            wparT = wparT + wt(g) * exp(-s * sum((xt([3 4]) - Ou(:,g)).^2));
        end
        accT = zeros(m,1); mxStep = 0; dr = 0;
        for g = 1:nE
            [Xt, st] = steady_batch(TA(:,:,1:m), Uu(:,g), P, tol, ms);
            mxStep = max(mxStep, st);
            d = Xt([3 4],:) - Ou(:,g);
            accT = accT + wt(g) * exp(-s * sum(d.^2,1))';
            dr = dr + wt(g) * mean(sqrt(sum((Xt([3 4],:) - Xprod{g}(:,1:m)).^2,1)));
        end
        z.L_tol(L) = mean(1 - accT / wparT);
        z.wpar_tol(L) = wparT;
        z.step_tol(L) = mxStep;
        if L == nL
            dread = dr;
            dreadPar = sum(wt(:)' .* sqrt(sum((xparT - xp([3 4],:)).^2,1)));
        end
    end
    z.dread = dread; z.dreadPar = dreadPar;
    z.L_loose_sub = mean(Lk(1:m));
    z.Lpred_sub = mean(Pk(1:m));
    z.nTight = m;
end


% =====================================================================
function H = drive_field(TA, X, u)
%DRIVE_FIELD  the rectifier argument (W x + u) of every draw at its own fixed
%   point, column by column. Same expression steady_batch integrates.
    [n, M] = size(X);
    H = repmat(u, 1, M);
    for j = 1:n
        H = H + reshape(TA(:,j,:), n, M) .* X(j,:);
    end
end


function x = steady_tight(Th, u, P, tol, maxSteps)
%STEADY_TIGHT  The same integrator on one wiring, at a given tolerance and step limit.
    x = P.x0; K = P.K; mu = P.mu; dt = P.dev_dt;
    for st = 1:maxSteps
        b = Th*x + u; b(b < 0) = 0;
        xn = x + dt*((b ./ (K + b)) - mu*x);
        d = sum((xn - x).^2); x = xn;
        if d < tol, break; end
    end
end


% =====================================================================
function g = identity_beta(Th, Uu, wt, P)
%IDENTITY_BETA  The identity's prediction, as in tracking_load_vs_identity: C Phi and
%   x* at each window point, and the RMS output shift under the mutation distribution.
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


function Pk = lin_load_draws(TA, Th, g, wt, P)
%LIN_LOAD_DRAWS  The predicted load of each draw: its first-order shift C Phi dW x* scored by the
%   task's fitness over the window. Its mean over draws is the L_identLin of
%   tracking_load_vs_identity.
    M = size(TA,3); s = P.selection_strength; nE = numel(g.CP);
    accW = zeros(M,1);
    for e = 1:nE
        CP = g.CP{e}; xs = g.x{e};
        d2 = zeros(M,1);
        for k = 1:M
            dp = CP * ((TA(:,:,k) - Th) * xs);
            d2(k) = sum(dp.^2);
        end
        accW = accW + wt(e) * exp(-s * d2);
    end
    Pk = 1 - accW;
end


% =====================================================================
function routing_check(C, cfg, aList, arms)
%ROUTING_CHECK  Over the first 2000 draws, the mean L_k and Lpred_k of each replicate against
%   the L and L_identLin stored by tracking_load_vs_identity (largest difference).
    fprintf('---- CHECK AGAINST THE STORED RUN: load file %s ----\n', cfg.loadFile);
    fprintf('Same genotype, same seed, same draw order: over the FIRST 2000 draws (the draw count of\n');
    fprintf('tracking_load_vs_identity) mean_k L_k IS L and mean_k Lpred_k IS L_identLin, so both columns must be 0.\n\n');
    LR = load(cfg.loadFile); o = LR.out;
    fprintf('%-18s %18s %18s\n', 'cell', 'max|dL|', 'max|dL_identLin|');
    for k = 1:numel(C)
        ai = cfg.cells{k,1}; f = cfg.cells{k,2};
        a = o.S.(arms{ai}).L(f,:)';           b = C(k).L_route;
        c = o.S.(arms{ai}).L_identLin(f,:)';  d = C(k).Lpred_route;
        fprintf('%-18s %18.3g %18.3g\n', C(k).label, max(abs(a-b)), max(abs(c-d)));
    end
    fprintf('(aList %s)\n\n', mat2str(aList));
end


% =====================================================================
function report_sign_changes(C, cfg)
    fprintf('\n============ ACTIVE-SET CHANGES ============\n');
    fprintf('A draw counts as a sign change when any gene''s (W x* + u) changes sign between parent and mutant,\n');
    fprintf('in any environment of the window, at the mutant''s own fixed point.\n\n');
    fprintf('%-18s %10s %12s %14s %14s %12s %12s\n', 'cell','draws','sign-change %', ...
        'sign-change % CI', 'load share %', 'mean L sw', 'mean L ns');
    for k = 1:numel(C)
        sf = C(k).switchFrac; sh = C(k).L_sw ./ max(C(k).L, eps);
        ci = boot_mean_ci(sf, cfg.nBoot, 31 + k);
        shb = boot_ratio_ci(C(k).L_sw, C(k).L, cfg.nBoot, 41 + k);
        fprintf('%-18s %10d %12.3f %14s %14s %12.5f %12.5f\n', C(k).label, C(k).nD(1), ...
            100*mean(sf), sprintf('[%.2f,%.2f]', 100*ci(1), 100*ci(2)), ...
            sprintf('%.1f [%.0f,%.0f]', 100*mean(sh), 100*shb(1), 100*shb(2)), ...
            mean(C(k).meanLk_sw,'omitnan'), mean(C(k).meanLk_ns,'omitnan'));
    end
    fprintf('\n%-18s %14s %12s %38s\n', 'cell','parent margin','on-switch %','sign-change events per 100 draws, genes 1..4');
    for k = 1:numel(C)
        fprintf('%-18s %14.4f %12.1f    %s\n', C(k).label, min(C(k).marginParent), ...
            100*mean(C(k).swOnFrac,'omitnan'), mat2str(round(100*mean(C(k).swGene,1),3)));
    end
    fprintf('\n"parent margin" is the smallest |W x* + u| over genes, environments and replicates:\n');
    fprintf('how far the nearest gene sits from its own rectifier corner. The per-gene column\n');
    fprintf('counts gene-by-environment events, so in a 36-environment window it can exceed 100.\n');
end


% =====================================================================
function report_ratios(C, cfg)
    fprintf('\n\n============ THE IDENTITY, WITH AND WITHOUT THE SIGN-CHANGE DRAWS ============\n');
    fprintf('Ratios are of sums over replicates (the ratio of means), paired bootstrap, %d draws.\n', cfg.nBoot);
    fprintf('Both directions are printed; tracking_load_vs_identity prints measured/predicted.\n\n');
    fprintf('%-18s %-14s %11s %11s %9s %20s %9s\n', 'cell','draws','measured','predicted', ...
        'meas/pred','meas/pred 95% CI','pred/meas');
    for k = 1:numel(C)
        sets = {'all', C(k).L, C(k).Lpred; ...
                'no sign change', C(k).L_ns, C(k).Lpred_ns; ...
                'sign change', C(k).L_sw, C(k).Lpred_sw};
        for i = 1:3
            y = sets{i,2}; x = sets{i,3};
            if sum(x) <= 0 || sum(y) <= 0
                fprintf('%-18s %-14s %11.5f %11.5f %9s  %20s %9s\n', ...
                    C(k).label, sets{i,1}, mean(y), mean(x), 'n/a', '(no such draw)', 'n/a');
                continue
            end
            r = sum(y)/sum(x);
            ci = boot_ratio_ci(y, x, cfg.nBoot, 7000 + 10*k + i);
            fprintf('%-18s %-14s %11.5f %11.5f %9.3f  [%7.3f, %7.3f] %9.3f\n', ...
                C(k).label, sets{i,1}, mean(y), mean(x), r, ci(1), ci(2), 1/max(r,eps));
        end
        exc = mean(C(k).L) - mean(C(k).Lpred);
        excSw = mean(C(k).L_sw) - mean(C(k).Lpred_sw);
        if exc > 1e-5
            fprintf('%-18s %-14s excess %+.5f, of which the sign-change draws carry %+.5f (%.0f%%)\n\n', ...
                C(k).label, '', exc, excSw, 100*excSw/exc);
        else
            fprintf('%-18s %-14s excess %+.5f (the identity does not under-predict here); the\n', ...
                C(k).label, '', exc);
            fprintf('%-18s %-14s sign-change draws'' own share of it is %+.5f\n\n', '', '', excSw);
        end
    end
end


% =====================================================================
function report_skew(C, cfg)
    fprintf('\n============ SKEW: ACROSS REPLICATES AND ACROSS DRAWS ============\n');
    fprintf('%-18s %10s %10s %10s %10s %12s %12s\n', 'cell','mean L','median L','max L', ...
        'min L','top-1% draws','max single draw');
    for k = 1:numel(C)
        fprintf('%-18s %10.5f %10.5f %10.5f %10.5f %12.1f%% %12.4f\n', C(k).label, ...
            mean(C(k).L), median(C(k).L), max(C(k).L), min(C(k).L), ...
            100*mean(C(k).top1share), mean(C(k).maxLk));
    end
    fprintf('\n"top-1%% draws" is the share of a replicate''s own measured load carried by its\n');
    fprintf('heaviest 1%% of draws, averaged over replicates.\n');

    fprintf('\n---- how many replicates carry the excess (L - Lpred, summed over replicates) ----\n');
    fprintf('%-18s %10s %10s %10s %10s %26s\n', 'cell','total','n>0','n for 50%','n for 80%','top 3 replicates (share)');
    for k = 1:numel(C)
        e = C(k).L - C(k).Lpred; tot = sum(e);
        [es, ord] = sort(e, 'descend'); cs = cumsum(es) / tot;
        n50 = find(cs >= 0.5, 1); n80 = find(cs >= 0.8, 1);
        if isempty(n50), n50 = NaN; end
        if isempty(n80), n80 = NaN; end
        fprintf('%-18s %10.4f %10d %10s %10s   r%d %.0f%%, r%d %.0f%%, r%d %.0f%%\n', ...
            C(k).label, tot, sum(e > 0), num2str(n50), num2str(n80), ...
            ord(1), 100*es(1)/tot, ord(2), 100*es(2)/tot, ord(3), 100*es(3)/tot);
    end

    fprintf('\n---- does a replicate''s sign-change rate track its excess? ----\n');
    fprintf('%-18s %12s %10s %12s %10s\n','cell','Spearman','p','Pearson','p');
    for k = 1:numel(C)
        e = C(k).L - C(k).Lpred; sf = C(k).switchFrac;
        [rs, ps] = corr_safe(sf, e, 'Spearman');
        [rp, pp] = corr_safe(sf, e, 'Pearson');
        fprintf('%-18s %12.3f %10.4f %12.3f %10.4f\n', C(k).label, rs, ps, rp, pp);
    end
    fprintf('\n(nBoot = %d everywhere above.)\n', cfg.nBoot);
end


% =====================================================================
function report_tolerance(C, cfg)
    fprintf('\n\n============ THE ALTERNATIVE: THE FIXED-POINT TOLERANCE ============\n');
    fprintf('The same draws re-solved with the stopping rule sum((x_{t+1}-x_t)^2) < tol tightened\n');
    fprintf('by the factors below, parent included. Production tol = 1e-6; a ladder, not one point,\n');
    fprintf('so a monotone approach to a limit is visible and a non-converged solve is not.\n\n');
    fac = cfg.tolFacs;
    fprintf('%-18s %7s %13s', 'cell', 'draws', 'L (prod)');
    for L = 1:numel(fac), fprintf(' %13s', sprintf('L (tol x%.0e)', fac(L))); end
    fprintf(' %13s %13s\n', 'predicted', 'tightest/pred');
    for k = 1:numel(C)
        fprintf('%-18s %7d %13.5f', C(k).label, C(k).nTight(1), mean(C(k).L_loose_sub));
        for L = 1:numel(fac), fprintf(' %13.5f', mean(C(k).L_tol(:,L))); end
        fprintf(' %13.5f %13.3f\n', mean(C(k).Lpred_sub), ...
            mean(C(k).L_tol(:,end)) / max(mean(C(k).Lpred_sub), eps));
    end
    fprintf('\n%-18s %14s', 'cell', 'w(parent)');
    for L = 1:numel(fac), fprintf(' %14s', sprintf('w (tol x%.0e)', fac(L))); end
    fprintf('\n');
    for k = 1:numel(C)
        fprintf('%-18s %14.5f', C(k).label, mean(C(k).w_par));
        for L = 1:numel(fac), fprintf(' %14.5f', mean(C(k).wpar_tol(:,L))); end
        fprintf('\n');
    end
    fprintf('\nRead that row FIRST. The tracking target was itself defined through the production\n');
    fprintf('integrator during the evolution run, so a tightened solve scores a different phenotype\n');
    fprintf('against an unchanged target. If w(parent) moves, the ladder is a sensitivity and not\n');
    fprintf('a correction.\n');

    fprintf('\n%-18s %22s %22s %16s\n', 'cell', 'readout move, mutants', 'readout move, parent', 'steps at tightest');
    for k = 1:numel(C)
        fprintf('%-18s %22.5f %22.5f %16.0f\n', C(k).label, mean(C(k).dread), ...
            mean(C(k).dreadPar), max(C(k).step_tol(:,end)));
    end
    fprintf('\n"readout move" is the mean distance the scored readout (genes 3 and 4) travels\n');
    fprintf('between the production solve and the tightest one. Compare it with the load''s own\n');
    fprintf('scale: the identity''s beta for these cells is 0.006 to 0.04.\n');

    fprintf('\n---- the quadratic form, for the saturation question ----\n');
    fprintf('%-18s %14s %14s %12s\n','cell','L (measured)','s*E[d^2]','ratio');
    for k = 1:numel(C)
        fprintf('%-18s %14.5f %14.5f %12.3f\n', C(k).label, mean(C(k).L), ...
            mean(C(k).Lquad), mean(C(k).L)/max(mean(C(k).Lquad),eps));
    end
end


% =====================================================================
function report_margins(C)
    fprintf('\n\n============ WHY THE SMOOTH PART IS SMALL HERE ============\n');
    fprintf('%-18s %12s %12s %14s\n','cell','w(parent)','||x*||','predicted load');
    for k = 1:numel(C)
        fprintf('%-18s %12.5f %12.4f %14.5f\n', C(k).label, mean(C(k).w_par), ...
            mean(C(k).xnorm), mean(C(k).Lpred));
    end
end


% =====================================================================
function [r, p] = corr_safe(x, y, type)
    ok = ~isnan(x) & ~isnan(y);
    try
        [r, p] = corr(x(ok), y(ok), 'type', type);
    catch
        r = NaN; p = NaN;
    end
end

function ci = boot_ratio_ci(y, x, nBoot, seed)
%BOOT_RATIO_CI  percentile bootstrap of sum(y)/sum(x), paired: one resample of
%   the replicate index serves both sums.
    ok = ~isnan(y) & ~isnan(x); y = y(ok); x = x(ok); n = numel(y);
    if n == 0, ci = [NaN NaN]; return; end
    st = RandStream('twister','Seed',seed);
    br = nan(nBoot,1);
    for b = 1:nBoot
        idx = randi(st, n, n, 1);
        br(b) = sum(y(idx)) / max(sum(x(idx)), eps);
    end
    ci = prctile(br, [2.5 97.5]);
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
