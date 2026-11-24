function share_ending_membership(cfg)
% SHARE_ENDING_MEMBERSHIP  Per signaling architecture, what share of single mutations to a member
% end membership? Averaged within each snapshot, then over snapshots and replicates.
%
% Produces : SI S3.2, "Per mutation, the share of mutations that end membership falls from graded
%            separate to saturating contrast"; Table S3 (the row "Share of mutations ending
%            membership")
% Inputs   : publication/<architecture>_100reps.mat   (written by run_signaling_architectures)
%            publication/member_entry_exit_flux.mat   (measured exit per lineage, for the comparison
%              printed after the shares; part of the data archive)
% Outputs  : publication/member_exit_supply.mat   (qbar per architecture, replicate and snapshot;
%              `endsMembershipMean`, the share per architecture, and its interval ci; snapGen, the
%              generations read; the comparison)
%            prints per architecture the share with its interval, then the comparison
% Run      : share_ending_membership   (100 replicates; long, not timed here: 3 replicates take
%              about 3 min)
% Seed     : rng(20260903), generator not named (MATLAB's default, Mersenne twister): the
%            mutation draws, then the bootstrap, from one stream; membership_test seeds each
%            wiring itself and restores the stream
%
% Five stored generations per replicate from generation 100 on, spaced evenly by index (SI S3.2:
% three were stored one generation before the next, so the shares may run high). Every distinct
% member gets 30 draws of the mutation distribution (connection uniform; toggled with probability
% 0.2, else its weight gains N(0, 0.2^2)), each scored by the membership test. A snapshot's share
% is qbar = sum_w n_w p(w) / sum_w n_w over its members w (n_w individuals carry w; p(w) is its
% share of draws ending membership); the interval resamples replicates.

    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    if ~isfield(cfg,'dataDir') || isempty(cfg.dataDir)
        cfg.dataDir = fullfile(grn_data_file('publication'));
    end
    if ~isfield(cfg,'nReps'),           cfg.nReps = 100;      end
    if ~isfield(cfg,'nSnaps'),          cfg.nSnaps = 5;       end
    if ~isfield(cfg,'K'),               cfg.K = 30;           end
    if ~isfield(cfg,'steadyStateFrom'), cfg.steadyStateFrom = 100; end
    if ~isfield(cfg,'lambda'),          cfg.lambda = 0.005;   end
    if ~isfield(cfg,'pToggle'),         cfg.pToggle = 0.2;    end
    if ~isfield(cfg,'sdWeight'),        cfg.sdWeight = 0.2;   end
    if ~isfield(cfg,'nBoot'),           cfg.nBoot = 10000;    end
    if ~isfield(cfg,'rngSeed'),         cfg.rngSeed = 20260903; end
    rng(cfg.rngSeed);

    sensing_modes = { ...
        struct('name','graded_separate',     'sensingMode','separate',     'hillCoefficient',Inf,'discreteOutput',false), ...
        struct('name','graded_contrast',     'sensingMode','contrast', 'hillCoefficient',Inf,'discreteOutput',false), ...
        struct('name','saturating_separate', 'sensingMode','separate',     'hillCoefficient',1,  'discreteOutput',false), ...
        struct('name','saturating_contrast', 'sensingMode','contrast', 'hillCoefficient',1,  'discreteOutput',false)};
    nModes = numel(sensing_modes);

    qbar      = nan(nModes, cfg.nReps, cfg.nSnaps);   % weighted share of draws ending membership
    share     = nan(nModes, cfg.nReps, cfg.nSnaps);
    nDistinct = nan(nModes, cfg.nReps, cfg.nSnaps);
    nMemDist  = nan(nModes, cfg.nReps, cfg.nSnaps);   % distinct MEMBER wirings
    minPexit  = nan(nModes, cfg.nReps, cfg.nSnaps);
    maxPexit  = nan(nModes, cfg.nReps, cfg.nSnaps);
    snapGen   = nan(nModes, cfg.nReps, cfg.nSnaps);
    dgMode    = nan(nModes,1);

    totalStart = tic;
    for m = 1:nModes
        mn = sensing_modes{m}.name;
        f = dir(fullfile(cfg.dataDir, sprintf('%s_100reps.mat', mn)));
        if isempty(f), error('No file matching: %s_100reps.mat', mn); end
        [~,ix] = max([f.datenum]);
        d = load(fullfile(cfg.dataDir, f(ix).name), 'results');
        useReps = min(cfg.nReps, numel(d.results));
        params = struct('sensingMode', sensing_modes{m}.sensingMode, ...
            'hillCoefficient', sensing_modes{m}.hillCoefficient, ...
            'discreteOutput',  sensing_modes{m}.discreteOutput, ...
            'michaelisConstant',1.0,'degradationRate',0.1,'timeStep',1.0, ...
            'developmentSteps',600,'backend','simulate_pair');
        fprintf('\n=== %s (%d reps, %d snaps, K = %d) ===\n', mn, useReps, cfg.nSnaps, cfg.K);
        modeStart = tic;
        for r = 1:useReps
            rep = d.results{r};
            sg = double(rep.historyGenerations(:))';
            dgMode(m) = mean(diff(sg));
            ks = find(sg >= cfg.steadyStateFrom);
            ks = ks(round(linspace(1, numel(ks), cfg.nSnaps)));
            for s = 1:numel(ks)
                k = ks(s);
                M = double(rep.historyM(:,:,:,k)); B = double(rep.historyB(:,:,:,k));
                N = size(M,3); nG = size(M,1);
                Wall = reshape(M.*B, nG*nG, N)';
                [Wu, ~, ic] = unique(Wall, 'rows');
                counts = accumarray(ic, 1);
                nU = size(Wu,1);
                isMem = false(nU,1); pExit = nan(nU,1);
                for u = 1:nU
                    W = reshape(Wu(u,:), nG, nG);
                    isMem(u) = membership_test(W, params).is_css;
                end
                % Each distinct member (M and B of its first carrier) receives K draws of the
                % mutation distribution; p(w) is the share whose mutant fails the membership test
                firstIdx = accumarray(ic, (1:N)', [], @min);
                for u = find(isMem)'
                    Mi = M(:,:,firstIdx(u)); Bi = B(:,:,firstIdx(u));
                    outs = 0;
                    for t = 1:cfg.K
                        Mm = Mi; Bm = Bi;
                        g1 = randi(nG); g2 = randi(nG);
                        if rand() < cfg.pToggle
                            Mm(g1,g2) = ~Mm(g1,g2);
                        else
                            Bm(g1,g2) = Bm(g1,g2) + randn()*cfg.sdWeight;
                        end
                        if ~membership_test(Mm.*Bm, params).is_css, outs = outs + 1; end
                    end
                    pExit(u) = outs / cfg.K;
                end
                share(m,r,s)     = sum(counts(isMem)) / N;
                nDistinct(m,r,s) = nU;
                nMemDist(m,r,s)  = sum(isMem);
                snapGen(m,r,s)   = sg(k);
                if any(isMem)
                    qbar(m,r,s)     = sum(counts(isMem).*pExit(isMem)) / sum(counts(isMem));
                    minPexit(m,r,s) = min(pExit(isMem));
                    maxPexit(m,r,s) = max(pExit(isMem));
                end
            end
            if mod(r,5) == 0
                fprintf('  rep %d/%d (%.0fs elapsed, ETA %.0fs)\n', r, useReps, ...
                    toc(modeStart), toc(modeStart)*(useReps-r)/r);
            end
        end
        fprintf('%s: mean qbar = %.4g, mean share of members = %.3f, mean distinct members = %.1f (%.1f min)\n', ...
            mn, mean(qbar(m,:,:),'all','omitnan'), mean(share(m,:,:),'all','omitnan'), ...
            mean(nMemDist(m,:,:),'all','omitnan'), toc(modeStart)/60);
    end

    % The share per architecture: mean over snapshots within a replicate, then over replicates,
    % with a percentile bootstrap over replicates. The comparison that follows composes it into
    % an exit probability per interval, 1 - (1 - lambda qbar)^dg (dg the mean spacing of the
    % stored generations), against the measured exit of member_entry_exit_flux.mat; no paper number
    % comes from that comparison, but its resampling draws alternate with those of the share's
    % interval in one stream, so it stays: without it the interval would change.
    fx = load(fullfile(cfg.dataDir,'member_entry_exit_flux.mat'), 'exitN','exitD','perRepExitN','perRepExitD');
    Pmeas = fx.exitN(:) ./ max(1, fx.exitD(:));

    repEndsMembership = squeeze(mean(qbar, 3, 'omitnan'));           % nModes x nReps
    endsMembershipMean = mean(repEndsMembership, 2, 'omitnan');
    bootMean = nan(nModes, cfg.nBoot); bootR = nan(1, cfg.nBoot);
    bootMeas = nan(nModes, cfg.nBoot);
    nRepMeas = size(fx.perRepExitN, 2);
    for b = 1:cfg.nBoot
        idx = randi(useReps, 1, useReps);
        bm = mean(repEndsMembership(:, idx), 2, 'omitnan');
        bootMean(:,b) = bm;
        bootR(b) = bm(1) / max(bm(4), eps);
        jdx = randi(nRepMeas, 1, nRepMeas);
        bootMeas(:,b) = sum(fx.perRepExitN(:,jdx), 2) ./ max(1, sum(fx.perRepExitD(:,jdx), 2));
    end
    ci      = prctile(bootMean, [2.5 97.5], 2);
    ciMeas  = prctile(bootMeas, [2.5 97.5], 2);
    Ppred   = 1 - (1 - cfg.lambda .* endsMembershipMean).^dgMode;
    ratio   = Pmeas ./ Ppred;
    Rpred   = endsMembershipMean(1) / endsMembershipMean(4);
    RpredCI = prctile(bootR, [2.5 97.5]);
    Rmeas   = Pmeas(1) / Pmeas(4);
    inCI    = Ppred >= ciMeas(:,1) & Ppred <= ciMeas(:,2);

    fprintf('\n=== SHARE ENDING MEMBERSHIP vs MEASURED EXIT (per interval, dg ~ %.1f gens, lambda = %.3f) ===\n', mean(dgMode), cfg.lambda);
    fprintf('%-22s %9s %20s %8s %9s %9s %20s %7s %5s\n', ...
        'architecture','qbar','95% CI','share','P_pred','P_meas','P_meas 95% CI','ratio','inCI');
    for m = 1:nModes
        fprintf('%-22s %9.4g [%8.4g, %8.4g] %8.3f %9.4g %9.4g [%8.4g, %8.4g] %7.2f %5d\n', ...
            sensing_modes{m}.name, endsMembershipMean(m), ci(m,1), ci(m,2), ...
            mean(share(m,:,:),'all','omitnan'), Ppred(m), Pmeas(m), ciMeas(m,1), ciMeas(m,2), ...
            ratio(m), inCI(m));
    end
    fprintf('Range ratio (graded-separate / sat-contrast): predicted %.2f [%.2f, %.2f]   measured %.2f\n', ...
        Rpred, RpredCI(1), RpredCI(2), Rmeas);
    fprintf('Ratio spread across architectures: max/min = %.2f\n', max(ratio)/min(ratio));
    fprintf('GRAND TOTAL: %.1f min\n', toc(totalStart)/60);

    metadata = struct('cfg', cfg, 'developmentSteps', 600, 'backend', 'simulate_pair', ...
        'modes', {{sensing_modes{1}.name, sensing_modes{2}.name, sensing_modes{3}.name, sensing_modes{4}.name}});
    save(fullfile(cfg.dataDir,'member_exit_supply.mat'), 'qbar','share','nDistinct','nMemDist', ...
        'minPexit','maxPexit','snapGen','dgMode','repEndsMembership','endsMembershipMean','ci','ciMeas', ...
        'Ppred','Pmeas','ratio','Rpred','RpredCI','Rmeas','inCI','metadata');
    fprintf('Saved: %s\n', fullfile(cfg.dataDir,'member_exit_supply.mat'));
end
