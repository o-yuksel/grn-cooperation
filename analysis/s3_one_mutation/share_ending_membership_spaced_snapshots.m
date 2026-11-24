function share_ending_membership_spaced_snapshots(cfg)
% SHARE_ENDING_MEMBERSHIP_SPACED_SNAPSHOTS  At saturating contrast, what share of single mutations
% to a member end membership when read only on snapshots stored at least ten generations before
% the next?
%
% Produces : SI S3.2, "Re-measured on snapshots stored at least ten generations before the next,
%            saturating contrast reads 0.340 [0.315, 0.366]" (q_now and its interval, printed)
% Inputs   : publication/saturating_contrast_100reps.mat   (written by run_signaling_architectures;
%              the newest matching file)
%            publication/member_entry_exit_flux.mat   (measured exit per lineage, for the comparison
%              printed last; part of the data archive)
%            cfg   optional overrides of the defaults below (modeName, nReps, nIntervals,
%              minIntervalGens, K, steadyStateFrom, lambda, nBoot, rngSeed, dataDir, ...)
% Outputs  : publication/member_exit_robustness.mat   (per replicate the weighted sums; qNow, qFwd,
%              sRob and their intervals ciQ; the predicted and measured exit; metadata); prints
%              the share of members, the member wirings read, the three shares and the comparison
% Run      : share_ending_membership_spaced_snapshots   (100 replicates, about 1 h)
% Seed     : rng(20260903), generator not named (MATLAB's default, Mersenne twister): the
%            mutation draws, then the two bootstraps, from one stream; membership_test seeds each
%            wiring itself and restores the stream
%
% Up to 5 stored generations per replicate from generation 100 on, chosen evenly by index among
% those stored at least 10 generations before the next. Each distinct member gets 30 draws of the
% mutation distribution, scored by the membership test; p(w) is its share ending membership.
% q_now (qNow) = sum_w n_w p(w) / sum_w n_w, n_w the individuals carrying w, pooled over snapshots
% and replicates; q_fwd (qFwd) weights w instead by its descendants at the next stored generation,
% so sRob = q_fwd / q_now asks whether members whose lineages persist lose membership to a
% mutation less often. Intervals: percentile bootstrap over replicates.

    thisDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(thisDir));
    addpath(thisDir);

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    if ~isfield(cfg,'dataDir') || isempty(cfg.dataDir)
        cfg.dataDir = fullfile(grn_data_file('publication'));
    end
    if ~isfield(cfg,'modeName'),        cfg.modeName = 'saturating_contrast'; end
    if ~isfield(cfg,'sensingMode'),     cfg.sensingMode = 'contrast'; end
    if ~isfield(cfg,'hillCoefficient'), cfg.hillCoefficient = 1; end
    if ~isfield(cfg,'discreteOutput'),  cfg.discreteOutput = false; end
    if ~isfield(cfg,'nReps'),           cfg.nReps = 100;      end
    if ~isfield(cfg,'nIntervals'),      cfg.nIntervals = 5;   end
    if ~isfield(cfg,'minIntervalGens'), cfg.minIntervalGens = 10; end
    if ~isfield(cfg,'K'),               cfg.K = 30;           end
    if ~isfield(cfg,'steadyStateFrom'), cfg.steadyStateFrom = 100; end
    if ~isfield(cfg,'lambda'),          cfg.lambda = 0.005;   end
    if ~isfield(cfg,'nBoot'),           cfg.nBoot = 10000;    end
    if ~isfield(cfg,'rngSeed'),         cfg.rngSeed = 20260903; end
    rng(cfg.rngSeed);

    params = struct('sensingMode', cfg.sensingMode, ...
        'hillCoefficient', cfg.hillCoefficient, ...
        'discreteOutput',  cfg.discreteOutput, ...
        'michaelisConstant',1.0,'degradationRate',0.1,'timeStep',1.0, ...
        'developmentSteps',600,'backend','simulate_pair');

    f = dir(fullfile(cfg.dataDir, sprintf('%s_100reps.mat', cfg.modeName)));
    if isempty(f), error('No file matching: %s_100reps.mat', cfg.modeName); end
    [~,ix] = max([f.datenum]);
    d = load(fullfile(cfg.dataDir, f(ix).name), 'results');
    useReps = min(cfg.nReps, numel(d.results));

    % per-replicate sums, so the bootstrap can resample replicates
    sNow = zeros(useReps,1); wNow = zeros(useReps,1);   % sum n_w p_exit, sum n_w
    sFwd = zeros(useReps,1); wFwd = zeros(useReps,1);   % sum d_w p_exit, sum d_w
    nMemW = zeros(useReps,1); nDeadW = zeros(useReps,1);
    shareRep = nan(useReps,1); dgAll = [];

    fprintf('=== %s: forward-weighted share of mutations ending membership (%d reps, %d intervals, K = %d) ===\n', ...
        cfg.modeName, useReps, cfg.nIntervals, cfg.K);
    t0 = tic;
    for r = 1:useReps
        rep = d.results{r};
        sg = double(rep.historyGenerations(:))';
        P  = rep.historyParents;
        popSize = size(rep.historyM,3);
        % dg, the interval length of the predicted exit, is the mean over every stored interval
        % starting at steadyStateFrom or later, as in the measured exit it is compared with. The
        % stored generations are not evenly spaced (some are one generation apart).
        allc = find(sg >= cfg.steadyStateFrom);
        allc = allc(allc < numel(sg));
        dgAll(end+1) = mean(sg(allc+1) - sg(allc)); %#ok<AGROW>
        % Only snapshots at least minIntervalGens before the next are read: over a one-generation
        % interval selection has no time to act on which lineages persist.
        cand = allc(sg(allc+1) - sg(allc) >= cfg.minIntervalGens);
        if isempty(cand), continue; end
        ks = cand(round(linspace(1, numel(cand), min(cfg.nIntervals, numel(cand)))));
        ks = unique(ks);
        shareAcc = [];
        for kk = ks
            ga = sg(kk); gb = sg(kk+1);
            M = double(rep.historyM(:,:,:,kk)); B = double(rep.historyB(:,:,:,kk));
            nG = size(M,1); N = size(M,3);
            Wall = reshape(M.*B, nG*nG, N)';
            [Wu, ~, ic] = unique(Wall, 'rows');
            nU = size(Wu,1);
            counts = accumarray(ic, 1);
            firstIdx = accumarray(ic, (1:N)', [], @min);

            isMem = false(nU,1);
            for u = 1:nU
                isMem(u) = membership_test(reshape(Wu(u,:), nG, nG), params).is_css;
            end
            shareAcc(end+1) = sum(counts(isMem)) / N; %#ok<AGROW>
            memIdx = find(isMem)';
            if isempty(memIdx), continue; end

            % descendant count per distinct wiring: each individual at the next stored generation
            % is traced back through the parent records to its ancestor at this one
            desc = zeros(nU,1);
            for j = 1:popSize
                p = j; bad = false;
                for gen = gb-1 : -1 : ga
                    p = double(P(p, gen));
                    if p < 1 || isnan(p), bad = true; break; end
                end
                if ~bad, desc(ic(p)) = desc(ic(p)) + 1; end
            end

            % K draws of the mutation distribution per distinct member (M and B of its first
            % carrier); p_exit is the share whose mutant fails the membership test
            for u = memIdx
                Mi = M(:,:,firstIdx(u)); Bi = B(:,:,firstIdx(u));
                outs = 0;
                for t = 1:cfg.K
                    Mm = Mi; Bm = Bi;
                    g1 = randi(nG); g2 = randi(nG);
                    if rand() < 0.2
                        Mm(g1,g2) = ~Mm(g1,g2);
                    else
                        Bm(g1,g2) = Bm(g1,g2) + randn()*0.2;
                    end
                    if ~membership_test(Mm.*Bm, params).is_css, outs = outs + 1; end
                end
                pe = outs / cfg.K;
                sNow(r) = sNow(r) + counts(u)*pe;  wNow(r) = wNow(r) + counts(u);
                sFwd(r) = sFwd(r) + desc(u)*pe;    wFwd(r) = wFwd(r) + desc(u);
                nMemW(r) = nMemW(r) + 1;
                if desc(u) == 0, nDeadW(r) = nDeadW(r) + 1; end
            end
        end
        if ~isempty(shareAcc), shareRep(r) = mean(shareAcc); end
        if mod(r,10) == 0
            fprintf('  rep %d/%d (%.0fs elapsed, ETA %.0fs)\n', r, useReps, ...
                toc(t0), toc(t0)*(useReps-r)/r);
        end
    end

    % The shares are ratios of sums pooled over replicates; the bootstrap resamples replicates
    dg = mean(dgAll);
    qNow = sum(sNow) / max(sum(wNow), eps);
    qFwd = sum(sFwd) / max(sum(wFwd), eps);
    sRob = qFwd / max(qNow, eps);

    bq = nan(cfg.nBoot,3);
    for b = 1:cfg.nBoot
        idx = randi(useReps, 1, useReps);
        a = sum(sNow(idx)) / max(sum(wNow(idx)), eps);
        c = sum(sFwd(idx)) / max(sum(wFwd(idx)), eps);
        bq(b,:) = [a, c, c/max(a,eps)];
    end
    ciQ = prctile(bq, [2.5 97.5], 1);

    % Comparison: each share composed into an exit probability per interval,
    % 1 - (1 - lambda q)^dg, against the measured exit of member_entry_exit_flux.mat (fields
    % exitN, exitD per architecture, perRepExitN, perRepExitD per replicate); no paper number
    % comes from it
    fx = load(fullfile(cfg.dataDir,'member_entry_exit_flux.mat'), 'exitN','exitD','perRepExitN','perRepExitD');
    modes = {'graded_separate','graded_contrast','saturating_separate','saturating_contrast'};
    mi = find(strcmp(modes, cfg.modeName), 1);
    Pmeas = fx.exitN(mi) / max(1, fx.exitD(mi));
    bm = nan(cfg.nBoot,1); nRepMeas = size(fx.perRepExitN,2);
    for b = 1:cfg.nBoot
        jdx = randi(nRepMeas, 1, nRepMeas);
        bm(b) = sum(fx.perRepExitN(mi,jdx)) / max(1, sum(fx.perRepExitD(mi,jdx)));
    end
    ciMeas = prctile(bm, [2.5 97.5]);

    Pnow  = 1 - (1 - cfg.lambda*qNow)^dg;
    Pfwd  = 1 - (1 - cfg.lambda*qFwd)^dg;
    inCInow = Pnow >= ciMeas(1) && Pnow <= ciMeas(2);
    inCIfwd = Pfwd >= ciMeas(1) && Pfwd <= ciMeas(2);
    if Pnow > Pmeas, gapClosed = (Pnow - Pfwd) / (Pnow - Pmeas); else, gapClosed = NaN; end

    fprintf('\n=== FORWARD-WEIGHTED SHARE OF MUTATIONS ENDING MEMBERSHIP, %s (dg = %.1f gens) ===\n', cfg.modeName, dg);
    fprintf('  share of members              %.3f\n', mean(shareRep,'omitnan'));
    fprintf('  distinct member wirings       %d  (of which no descendants: %d, %.1f%%)\n', ...
        sum(nMemW), sum(nDeadW), 100*sum(nDeadW)/max(1,sum(nMemW)));
    fprintf('  q_now (clone-weighted)        %.4f [%.4f, %.4f]\n', qNow, ciQ(1,1), ciQ(2,1));
    fprintf('  q_fwd (descendant-weighted)   %.4f [%.4f, %.4f]\n', qFwd, ciQ(1,2), ciQ(2,2));
    fprintf('  s_rob = q_fwd / q_now         %.4f [%.4f, %.4f]\n', sRob, ciQ(1,3), ciQ(2,3));
    fprintf('  P_pred from q_now             %.4f   (inside measured CI: %d)\n', Pnow, inCInow);
    fprintf('  P_pred from q_fwd             %.4f   (inside measured CI: %d)\n', Pfwd, inCIfwd);
    fprintf('  P_meas                        %.4f [%.4f, %.4f]\n', Pmeas, ciMeas(1), ciMeas(2));
    fprintf('  share of the gap closed       %.0f%%\n', 100*gapClosed);

    fprintf('GRAND TOTAL: %.1f min\n', toc(t0)/60);

    metadata = struct('cfg', cfg, 'developmentSteps', 600, 'backend', 'simulate_pair', ...
        'dg', dg);
    save(fullfile(cfg.dataDir,'member_exit_robustness.mat'), 'sNow','wNow','sFwd','wFwd', ...
        'nMemW','nDeadW','shareRep','qNow','qFwd','sRob','ciQ','Pnow','Pfwd','Pmeas', ...
        'ciMeas','gapClosed','metadata');
    fprintf('Saved: %s\n', fullfile(cfg.dataDir,'member_exit_robustness.mat'));
end
