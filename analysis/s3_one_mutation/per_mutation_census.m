function per_mutation_census(cfg)
% PER_MUTATION_CENSUS  What does one random mutation do to an evolved member? Gives every distinct
% member at five stored generations of each replicate 30 draws of the mutation distribution and
% records, per draw, whether the mutant ends membership and whether it ends cooperation: the
% per-mutation census.
%
% Produces : the census behind Fig. 4A (sample_census_members draws its members from it, through
%            run_fixed_opponent_census), Table S6 (tableS6_per_mutation_shares) and SI S3.3's comparison with
%            cooperators that do not retaliate
%            (cooperate_only_losses). Also read by retaliation_margin_vs_rho,
%            measure_retaliation_fixed_opponent, response_dial_networks and measure_joint_loss
% Inputs   : publication/<architecture>_100reps.mat   (written by run_signaling_architectures;
%              the newest match is read)
% Outputs  : publication/rho_by_architecture.mat   (`W_all`: per architecture, one row per member
%              wiring and stored generation: replicate, generation, count (individuals carrying
%              it), ok, readout r, ||x*||, alpha, beta, and per draw rho (the readout shift a
%              measured slope predicts, as a fraction of r), the weight change dW, the connection
%              (`slot`), toggle or not, and whether it ends membership (`ejected`) or cooperation
%              (`lostReadout`); summary tables; cfg); prints weighted medians per architecture
% Run      : per_mutation_census   (100 replicates; long, not timed here: 4 replicates per
%              architecture take about 3 min on 3 workers)
% Seed     : rng(20260904 + 1000 m + r) for architecture m and replicate r, generator not named
%            (in a pool, the worker's default, Threefry); membership_test seeds each wiring itself
%
% Kept members (ok): cooperation margin (readout r = x_C - x_D against the clone) at least 1e-6
% and ||x*|| in [0.5, 14]; with a margin near zero, rho would divide by a decayed residue.

    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    % Settings: dataDir the publication folder under the data root; nSnaps stored generations
    % from steadyStateFrom on; K draws per member; the mutation distribution (toggle with
    % probability pToggle, else a weight change N(0, sdWeight^2); Methods, "Game and
    % evolution"); the ok rule (minXnorm, maxXnorm, minR); deadRho, a cutoff for the printed
    % tables only, not a value the paper states
    if ~isfield(cfg,'dataDir') || isempty(cfg.dataDir)
        cfg.dataDir = fullfile(grn_data_file('publication'));
    end
    if ~isfield(cfg,'outDir') || isempty(cfg.outDir), cfg.outDir = cfg.dataDir; end
    if ~isfield(cfg,'nReps'),           cfg.nReps = 100;      end
    if ~isfield(cfg,'nSnaps'),          cfg.nSnaps = 5;       end
    if ~isfield(cfg,'K'),               cfg.K = 30;           end
    if ~isfield(cfg,'steadyStateFrom'), cfg.steadyStateFrom = 100; end
    if ~isfield(cfg,'pToggle'),         cfg.pToggle = 0.2;    end
    if ~isfield(cfg,'sdWeight'),        cfg.sdWeight = 0.2;   end
    if ~isfield(cfg,'rngSeed'),         cfg.rngSeed = 20260904; end
    if ~isfield(cfg,'useParallel'),     cfg.useParallel = true; end
    if ~isfield(cfg,'minXnorm'),        cfg.minXnorm = 0.5;   end
    if ~isfield(cfg,'maxXnorm'),        cfg.maxXnorm = 14;    end
    if ~isfield(cfg,'minR'),            cfg.minR = 1e-6;      end
    if ~isfield(cfg,'deadRho'),         cfg.deadRho = 1e-3;   end

    sensing_modes = { ...
        struct('name','graded_separate',     'sensingMode','separate',     'hillCoefficient',Inf,'discreteOutput',false), ...
        struct('name','graded_contrast',     'sensingMode','contrast', 'hillCoefficient',Inf,'discreteOutput',false), ...
        struct('name','saturating_separate', 'sensingMode','separate',     'hillCoefficient',1,  'discreteOutput',false), ...
        struct('name','saturating_contrast', 'sensingMode','contrast', 'hillCoefficient',1,  'discreteOutput',false)};
    nModes = numel(sensing_modes);
    modeNames = cell(nModes,1);
    for m = 1:nModes, modeNames{m} = sensing_modes{m}.name; end

    if cfg.useParallel && isempty(gcp('nocreate')), parpool('local'); end

    W_all = cell(nModes,1);   % per architecture: struct array, one row per member wiring
    totalStart = tic;

    for m = 1:nModes
        mn = sensing_modes{m}.name;
        f = dir(fullfile(cfg.dataDir, sprintf('%s_100reps.mat', mn)));
        if isempty(f), error('No file matching: %s_100reps.mat', mn); end
        [~,ix] = max([f.datenum]);
        fprintf('\n=== %s: extracting snapshots ===\n', mn);
        d = load(fullfile(cfg.dataDir, f(ix).name), 'results');
        useReps = min(cfg.nReps, numel(d.results));
        snaps = extract_snapshots(d.results, 1:useReps, cfg);
        clear d;

        params = struct('sensingMode', sensing_modes{m}.sensingMode, ...
            'hillCoefficient', sensing_modes{m}.hillCoefficient, ...
            'discreteOutput',  sensing_modes{m}.discreteOutput, ...
            'michaelisConstant',1.0,'degradationRate',0.1,'timeStep',1.0, ...
            'developmentSteps',600,'backend','simulate_pair');
        sim_params = struct('sensingMode', sensing_modes{m}.sensingMode, ...
            'hillCoefficient', sensing_modes{m}.hillCoefficient, ...
            'discreteOutput',  sensing_modes{m}.discreteOutput, ...
            'michaelisConstant',1.0,'degradationRate',0.1,'timeStep',1.0, ...
            'developmentSteps',600,'payoffDC',0,'sensorIntegrationRate',1.0, ...
            'externalSignal1',zeros(4,1),'externalSignal2',zeros(4,1));

        fprintf('=== %s: %d reps x %d snaps, K = %d ===\n', mn, useReps, cfg.nSnaps, cfg.K);
        modeStart = tic;
        rows = cell(useReps,1);
        % Replicate r of architecture m runs under rng(cfg.rngSeed + 1000 m + r), so the draws
        % do not depend on how parfor schedules the replicates
        seed0 = cfg.rngSeed + 1000*m;
        if cfg.useParallel
            parfor r = 1:useReps
                rows{r} = process_rep(snaps{r}, params, sim_params, cfg, seed0 + r);
            end
        else
            for r = 1:useReps
                rows{r} = process_rep(snaps{r}, params, sim_params, cfg, seed0 + r);
            end
        end
        W_all{m} = vertcat(rows{:});
        fprintf('%s: %d member wirings, %.1f min\n', mn, numel(W_all{m}), toc(modeStart)/60);
    end

    report_and_save(W_all, modeNames, cfg);
    fprintf('GRAND TOTAL: %.1f min\n', toc(totalStart)/60);
end


% =========================================================================
function rows = process_rep(snapRep, params, sim_params, cfg, seed)
% One replicate under its own seed: the membership test on each distinct wiring, then the
% measures and draws of wiring_metrics on each member
    rng(seed);
    rows = [];
    S = snapRep.S;
    for s = 1:numel(S)
        Wu = S(s).Wu; counts = S(s).counts; nU = size(Wu,1); nG = 4;
        isMem = false(nU,1);
        for u = 1:nU
            isMem(u) = membership_test(reshape(Wu(u,:),nG,nG), params).is_css;
        end
        for u = find(isMem)'
            Mi = S(s).Mfirst(:,:,u); Bi = S(s).Bfirst(:,:,u);
            rec = wiring_metrics(Mi, Bi, params, sim_params, cfg);
            if isempty(rec), continue; end
            rec.rep = snapRep.rep; rec.gen = S(s).gen;
            rec.count = counts(u);  rec.N = S(s).N;
            if isempty(rows), rows = rec; else, rows(end+1) = rec; end %#ok<AGROW>
        end
    end
    rows = rows(:);
end


% =========================================================================
function rec = wiring_metrics(Mi, Bi, params, sim_params, cfg)
% Everything measured on one member wiring, with the 600-step development from the fixed start
% of the membership test, plus its cfg.K drawn mutations.
    rec = [];
    nG = 4; Msq = ones(nG); Kv = [];
    W = Mi .* Bi;

    [~, fs, ~] = simulate_pair(Msq, W, Kv, Msq, W, Kv, sim_params);
    x = fs(1:nG);
    r = x(1) - x(2);
    if ~isfinite(r) || r <= 0, return; end     % a member cooperates with its clone (r > 0)
    xn = norm(x);

    % --- alpha_SI, beta_SI: responsiveness and mutation-sensitivity of the readout
    %     (alpha_and_beta; SI S1.3), the partner's signal held at what the clone delivers;
    %     beta_SI = alpha_SI ||x*|| by the congruence identity.
    sPart = partner_signal(x, sim_params);
    th = alpha_and_beta(W, x, sPart, sim_params.degradationRate, ...
                         sim_params.michaelisConstant, [1;-1;0;0]);
    if th.valid
        alpha_SI = th.alpha; beta_SI = th.beta; v = th.v;
    else
        alpha_SI = NaN; beta_SI = NaN; v = nan(nG,1);
    end

    % --- beta_state = ||dx*/dW||_F of the whole state, finite differences of 1e-4 with both
    %     networks of the clone pair changed; NaN if a step flips the readout's sign.
    eps_W = 1e-4; dsq = 0; flipped = false;
    for ii = 1:nG
        for jj = 1:nG
            Wp = W; Wp(ii,jj) = Wp(ii,jj) + eps_W;
            [~, fp, ~] = simulate_pair(Msq, Wp, Kv, Msq, Wp, Kv, sim_params);
            xp = fp(1:nG);
            if sign(xp(1)-xp(2)) ~= sign(r), flipped = true; break; end
            dsq = dsq + sum(((xp - x)/eps_W).^2);
        end
        if flipped, break; end
    end
    if flipped, beta_state = NaN; else, beta_state = sqrt(dsq); end

    % --- the directional slopes: the clone readout's change per unit weight over a step of
    %     +sigma_mut (gF) and -sigma_mut (gB) at each slot.
    sg = cfg.sdWeight;
    gF = zeros(nG*nG,1); gB = zeros(nG*nG,1);
    for q = 1:nG*nG
        Wup = W; Wup(q) = Wup(q) + sg;
        Wdn = W; Wdn(q) = Wdn(q) - sg;
        [~, fu, ~] = simulate_pair(Msq, Wup, Kv, Msq, Wup, Kv, sim_params);
        [~, fd, ~] = simulate_pair(Msq, Wdn, Kv, Msq, Wdn, Kv, sim_params);
        gF(q) = ((fu(1)-fu(2)) - r) / sg;
        gB(q) = (r - (fd(1)-fd(2))) / sg;
    end

    % --- cfg.K draws of the mutation distribution (Methods, "Game and evolution"): connection
    %     (g1, g2) uniform; toggled with probability pToggle, else its weight gains
    %     N(0, sdWeight^2). Per draw: rho_dir = -(slope on the side of dW) dW / r and rho_an =
    %     -(v_i x*_j dW) / r, the predicted shift as a fraction of the cooperation margin;
    %     `ejected` = the mutant fails the membership test (ends membership); `lostReadout` = the
    %     mutant's readout with its own clone is not positive (ends cooperation).
    K = cfg.K;
    rho_dir = nan(K,1); rho_an = nan(K,1); dWv = nan(K,1);
    slot = nan(K,1); isTog = false(K,1); ejected = false(K,1); lostR = false(K,1);
    for t = 1:K
        Mm = Mi; Bm = Bi;
        g1 = randi(nG); g2 = randi(nG);
        tog = rand() < cfg.pToggle;
        if tog
            Mm(g1,g2) = ~Mm(g1,g2);
        else
            Bm(g1,g2) = Bm(g1,g2) + randn()*cfg.sdWeight;
        end
        Wm = Mm .* Bm;
        q = sub2ind([nG nG], g1, g2);
        dl = Wm(g1,g2) - W(g1,g2);
        if dl >= 0, g = gF(q); else, g = gB(q); end
        rho_dir(t) = -(g * dl) / r;
        if isfinite(beta_SI), rho_an(t) = -(v(g1) * x(g2) * dl) / r; end
        dWv(t) = dl; slot(t) = q; isTog(t) = tog;
        ejected(t) = ~membership_test(Wm, params).is_css;
        [~, fm, ~] = simulate_pair(Msq, Wm, Kv, Msq, Wm, Kv, sim_params);
        lostR(t) = ~((fm(1) - fm(2)) > 0);
    end

    ok = (xn >= cfg.minXnorm) && (xn <= cfg.maxXnorm) && (r >= cfg.minR);
    rec = struct('rep',0,'gen',0,'count',0,'N',0, 'ok', ok, ...
        'r', r, 'xn', xn, 'r_over_xn', r/xn, ...
        'alpha_SI', alpha_SI, 'beta_SI', beta_SI, 'beta_state', beta_state, ...
        'theorem_valid', th.valid, ...
        'pExit', mean(ejected), 'pRhoGe1', mean(rho_dir >= 1), ...
        'pRhoAnGe1', mean(rho_an >= 1), 'pLostReadout', mean(lostR), ...
        'pEjectNotReadout', mean(ejected & ~lostR), ...
        'absRhoMed', median(abs(rho_dir)), ...
        'rho_dir', rho_dir, 'rho_an', rho_an, 'dW', dWv, 'slot', slot, ...
        'isToggle', isTog, 'ejected', ejected, 'lostReadout', lostR);
end


% =========================================================================
function s = partner_signal(x, sim_params)
% What the partner delivers into the sensor channels at the clone pair's fixed point, by the
% sensing rule of simulate_pair (Methods, "Model"). It is held fixed under mutation, the setting
% in which SI S1.3 states alpha and beta.
    nG = numel(x);
    intRate = 1.0;
    if isfield(sim_params,'sensorIntegrationRate') && ~isempty(sim_params.sensorIntegrationRate)
        intRate = sim_params.sensorIntegrationRate;
    end
    gC = 1; gD = 2; gS1 = 3; gS2 = 4;
    if isfield(sim_params,'discreteOutput') && sim_params.discreteOutput
        oC = double(x(gC) > x(gD)); oD = 1 - oC;
    elseif isfinite(sim_params.hillCoefficient)
        n = sim_params.hillCoefficient;
        cn = x(gC)^n; dn = x(gD)^n; den = cn + dn + 1e-12;
        oC = cn/den; oD = dn/den;
    else
        oC = x(gC); oD = x(gD);
    end
    s = zeros(nG,1);
    if strcmp(stored_value('sensingMode', sim_params.sensingMode), 'separate')
        s(gS1) = intRate*oC; s(gS2) = intRate*oD;
    else
        sig = intRate*(oC - oD);
        s(gS1) = sig; s(gS2) = -sig;
    end
end


% =========================================================================
function report_and_save(W_all, modeNames, cfg)
% Summary tables per architecture over the ok members, weighted by how many individuals carry
% each wiring; then the save.
    nModes = numel(modeNames);

    QN = {'beta_SI','beta_state','alpha_SI','r','xn','r_over_xn','beta_over_r','pRhoGe1','pExit', ...
          'pLostReadout','pEjectNotReadout','c2_only_share','pRhoAnGe1'};
    lbl = {'beta_SI    = ||dr/dW||_F  (alpha_and_beta)', ...
           'beta_state = ||dx*/dW||_F (whole state)', ...
           'alpha_SI   = ||(grad r) Phi||_2', ...
           'r          = margin x_C - x_D', ...
           '||x*||', ...
           'r / ||x*||', ...
           'beta_SI / r  (the composite that sets rho; = alpha / (r/||x*||))', ...
           'P(rho_dir >= 1) per wiring', ...
           'P(ends membership) per wiring', ...
           'P(ends cooperation) per wiring', ...
           'P(ends only retaliation) per wiring', ...
           'share of the mutations that end membership that end only retaliation, per wiring', ...
           'P(rho_an >= 1) per wiring'};
    stats = struct();
    for i = 1:numel(QN), stats.(QN{i}) = nan(nModes,6); end
    rhoStat = nan(nModes,6); absRhoStat = nan(nModes,6); rhoAnStat = nan(nModes,6);
    absRhoWeightOnly = nan(nModes,6); absRhoLive = nan(nModes,6); deadFrac = nan(nModes,1);
    nW = nan(nModes,1); nOK = nan(nModes,1); validFrac = nan(nModes,1);

    for m = 1:nModes
        Rw = W_all{m};
        nW(m) = numel(Rw);
        if nW(m) == 0, continue; end
        wAll = double([Rw.count])';

        keep = logical([Rw.ok])';
        Rk = Rw(keep); w = wAll(keep);
        nOK(m) = numel(Rk);
        if nOK(m) == 0, continue; end
        validFrac(m) = mean([Rk.theorem_valid]);
        for i = 1:numel(QN)
            if strcmp(QN{i}, 'beta_over_r')
                % beta / r: rho_an = -(v_i x*_j dW) / r, and the products v_i x*_j have
                % Frobenius norm beta, so beta / r scales rho; by the congruence identity
                % beta / r = alpha / (r / ||x*||).
                v = double([Rk.beta_SI])' ./ double([Rk.r])';
            elseif strcmp(QN{i}, 'c2_only_share')
                pe = double([Rk.pExit])';
                v = double([Rk.pEjectNotReadout])' ./ pe;
                v(pe == 0) = NaN;
            else
                v = double([Rk.(QN{i})])';
            end
            stats.(QN{i})(m,:) = wstats(v, w);
        end
        % Each member stores its K draws as a K-by-1 column; stacking them keeps a member's
        % draws together, in member order, so repelem(w, K, 1) gives each draw its member's weight
        % (a column also when one member is kept)
        rhoAll = vertcat(Rk.rho_dir);
        wRho = repelem(w, cfg.K, 1);
        rhoStat(m,:)    = wstats(rhoAll, wRho);
        absRhoStat(m,:) = wstats(abs(rhoAll), wRho);
        tog = vertcat(Rk.isToggle);
        absRhoWeightOnly(m,:) = wstats(abs(rhoAll(~tog)), wRho(~tog));
        % The rectifier splits the draws into those at a slot whose directional slope is 0
        % (|rho| near 0) and the rest; a median over both mixes the two groups, so the share
        % below deadRho and |rho| over the rest are reported separately.
        live = abs(rhoAll) >= cfg.deadRho;
        deadFrac(m) = sum(wRho(~live)) / sum(wRho); %#ok<AGROW>
        absRhoLive(m,:) = wstats(abs(rhoAll(live)), wRho(live)); %#ok<AGROW>
        rhoAnAll = vertcat(Rk.rho_an);
        rhoAnStat(m,:) = wstats(rhoAnAll, wRho);
    end

    % ---------------- printed report ----------------
    fprintf('\n================ rho BY ARCHITECTURE (nReps = %d, K = %d, seed %d) ================\n', ...
        cfg.nReps, cfg.K, cfg.rngSeed);
    fprintf('Clone-size-weighted median [IQR] (range) n, over DISTINCT member wirings.\n\n');
    for i = 1:numel(QN)
        fprintf('--- %s\n', lbl{i});
        for m = 1:nModes
            S = stats.(QN{i})(m,:);
            fprintf('   %-21s %10.4g [%9.4g, %9.4g]  (%9.4g .. %9.4g)  n=%d\n', ...
                modeNames{m}, S(1), S(2), S(3), S(4), S(5), S(6));
        end
    end
    blocks = {'|rho_dir| over ALL drawn mutations', absRhoStat; ...
              '|rho_dir| over LIVE-slot draws only (|rho| >= deadRho)', absRhoLive; ...
              '|rho_dir| over WEIGHT mutations only', absRhoWeightOnly; ...
              'rho_dir (signed) over ALL drawn mutations', rhoStat; ...
              'rho_an (signed) over ALL drawn mutations', rhoAnStat};
    for bi = 1:size(blocks,1)
        Sb = blocks{bi,2};
        fprintf('--- %s\n', blocks{bi,1});
        for m = 1:nModes
            fprintf('   %-21s %10.4g [%9.4g, %9.4g]  (%9.4g .. %9.4g)  n=%d\n', ...
                modeNames{m}, Sb(m,1), Sb(m,2), Sb(m,3), Sb(m,4), Sb(m,5), Sb(m,6));
        end
    end

    fprintf('Member wirings: %s total, %s kept after the ||x*|| and margin filters\n', ...
        mat2str(nW'), mat2str(nOK'));
    fprintf('alpha_and_beta valid fraction (kept members): %s\n', mat2str(validFrac',3));
    metadata = struct('cfg', cfg, 'modes', {modeNames}, ...
        'backends', ['beta_SI = alpha_and_beta readout Jacobian (SI S1.3 routine); ' ...
                     'beta_state = ||dx*/dW||_F of the whole state, finite differences of 1e-4 on the clone pair; ' ...
                     'rho_dir = -(directional one-sided slope) dW / r; rho_an = -(v_i x*_j dW) / r, the rank-one form.']);
    save(fullfile(cfg.outDir,'rho_by_architecture.mat'), 'W_all','modeNames','stats','QN', ...
        'rhoStat','absRhoStat','absRhoWeightOnly','rhoAnStat','nW','nOK','validFrac', ...
        'absRhoLive','deadFrac', ...
        'metadata','-v7.3');
    fprintf('Saved: %s\n', fullfile(cfg.outDir,'rho_by_architecture.mat'));
end


% =========================================================================
function S = wstats(v, w)
% [weighted median, weighted q1, weighted q3, min, max, n] over finite entries.
    ok = isfinite(v) & isfinite(w) & w > 0;
    v = v(ok); w = w(ok);
    if isempty(v), S = nan(1,6); return; end
    S = [wquantile(v,w,0.5), wquantile(v,w,0.25), wquantile(v,w,0.75), ...
         min(v), max(v), numel(v)];
end
