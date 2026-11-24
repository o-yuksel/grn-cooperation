function retaliation_margin_vs_rho(cfg)
% RETALIATION_MARGIN_VS_RHO  How does one mutation move a member's retaliation margin? Replays the
% per-mutation census and scores the retaliation margin of each member and of each of its 30
% mutants, against rho: the mutation's predicted shift of the readout as a fraction of the
% cooperation margin, taken from the census (rho >= 1 predicts that cooperation ends).
%
% Produces : members with no defecting mutant (SI S2.3, "Members with no defecting mutant"; printed
%            as "parent worst_retal finite: k of n"); the input of two_margin_plane (SI S2.3) and
%            of run_fixed_opponent_census (Fig. 4A). The tables by rho and the correlations print but
%            feed no paper item, nor does the printed split of draws into dead and live.
% Inputs   : publication/rho_by_architecture.mat   (per_mutation_census; per member and draw rho
%              and the scored outcomes, for the check and the join)
%            publication/<architecture>_100reps.mat   (run_signaling_architectures; the mutants
%              are rebuilt from these, as the census stores no mutated wirings)
%            cfg   optional struct of settings; the defaults at the top of the code give the
%                  paper's numbers; nReps takes the first replicates only
% Outputs  : publication/retaliation_margin_vs_rho.mat (R_all: per member and draw the margins) and
%            retaliation_margin_vs_rho_summary.mat (the tables); console
% Run      : retaliation_margin_vs_rho   (parfor on a local pool; 3 replicates take about 10 min;
%              the full run is not timed)
% Seed     : replicate r of architecture m: rng(20260904 + 1000 m + r), the census's seeds, generator
%            not named, so the same mutations are drawn only on pool workers, as the census ran;
%            rng(20260905) for the bootstrap and the shuffle; membership_test seeds each wiring
%
% The retaliation margin is minus worst_retal, the member's readout against its worst-case mutant
% (membership_test with gradedC2 = true); a mutant that fails the cooperation condition (C1) gets
% none. Members kept: ||x*|| (the size of the settled expression state) in [0.5, 14] and
% cooperation margin >= 1e-6, as in the census. The run checks that it reproduces the stored census
% draw for draw and prints FAILED otherwise.

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    if ~isfield(cfg,'dataDir') || isempty(cfg.dataDir)
        cfg.dataDir = fullfile(grn_data_file('publication'));
    end
    if ~isfield(cfg,'srcFile') || isempty(cfg.srcFile)
        cfg.srcFile = fullfile(cfg.dataDir,'rho_by_architecture.mat');
    end
    if ~isfield(cfg,'outMat') || isempty(cfg.outMat)
        cfg.outMat = fullfile(cfg.dataDir,'retaliation_margin_vs_rho.mat');
    end
    if ~isfield(cfg,'nReps'),           cfg.nReps = 100;      end
    if ~isfield(cfg,'nSnaps'),          cfg.nSnaps = 5;       end
    if ~isfield(cfg,'K'),               cfg.K = 30;           end   % draws per wiring (RNG replay)
    if ~isfield(cfg,'steadyStateFrom'), cfg.steadyStateFrom = 100; end
    if ~isfield(cfg,'pToggle'),         cfg.pToggle = 0.2;    end
    if ~isfield(cfg,'sdWeight'),        cfg.sdWeight = 0.2;   end
    if ~isfield(cfg,'rngSeed'),         cfg.rngSeed = 20260904; end
    if ~isfield(cfg,'rngSeedStat'),     cfg.rngSeedStat = 20260905; end
    if ~isfield(cfg,'useParallel'),     cfg.useParallel = true; end
    if ~isfield(cfg,'binEdges'),        cfg.binEdges = -0.5:0.125:2.5; end
    if ~isfield(cfg,'minBinN'),         cfg.minBinN = 40;     end
    if ~isfield(cfg,'nBoot'),           cfg.nBoot = 400;      end

    sensing_modes = { ...
        struct('name','graded_separate',     'sensingMode','separate',     'hillCoefficient',Inf,'discreteOutput',false), ...
        struct('name','graded_contrast',     'sensingMode','contrast', 'hillCoefficient',Inf,'discreteOutput',false), ...
        struct('name','saturating_separate', 'sensingMode','separate',     'hillCoefficient',1,  'discreteOutput',false), ...
        struct('name','saturating_contrast', 'sensingMode','contrast', 'hillCoefficient',1,  'discreteOutput',false)};
    nModes = numel(sensing_modes);
    modeNames = cell(nModes,1);
    for m = 1:nModes, modeNames{m} = sensing_modes{m}.name; end

    if cfg.useParallel && isempty(gcp('nocreate')), parpool('local'); end

    R_all = cell(nModes,1);
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

        fprintf('=== %s: %d reps x %d snaps, K = %d draws ===\n', ...
            mn, useReps, cfg.nSnaps, cfg.K);
        modeStart = tic;
        rows = cell(useReps,1);
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
        R_all{m} = vertcat(rows{:});
        fprintf('%s: %d member wirings, %.1f min\n', mn, numel(R_all{m}), toc(modeStart)/60);
    end
    fprintf('\nREPLAY TOTAL: %.1f min\n', toc(totalStart)/60);

    repro = check_reproduction(R_all, modeNames, cfg);
    metadata = struct('cfg', cfg, 'modes', {modeNames});
    save(cfg.outMat, 'R_all','modeNames','repro','metadata','-v7.3');
    fprintf('Saved: %s\n', cfg.outMat);

    analyse_and_report(R_all, modeNames, repro, cfg);
end


% =========================================================================
function rows = process_rep(snapRep, params, sim_params, cfg, seed)
% One replicate in the census's order: the membership test on each distinct wiring of each
% snapshot, then the margins of each member; the stream is seeded once per replicate.
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
            rec = wiring_margins(Mi, Bi, params, sim_params, cfg);
            if isempty(rec), continue; end
            rec.rep = snapRep.rep; rec.gen = S(s).gen;
            rec.count = counts(u);  rec.N = S(s).N;
            if isempty(rows), rows = rec; else, rows(end+1) = rec; end %#ok<AGROW>
        end
    end
    rows = rows(:);
end


% =========================================================================
function rec = wiring_margins(Mi, Bi, params, sim_params, cfg)
% One member: its own retaliation margin, then cfg.K mutations from the mutation distribution (a
% toggle of one connection with probability pToggle, else that weight plus N(0, sdWeight^2)), with
% each mutant's margin. Random numbers are drawn in the census's order (randi, randi, rand, and
% randn on the weight branch only), and the early return comes before any draw, as there.
    rec = [];
    nG = 4; Msq = ones(nG); Kv = [];
    W = Mi .* Bi;

    [~, fs, ~] = simulate_pair(Msq, W, Kv, Msq, W, Kv, sim_params);
    x = fs(1:nG);
    r = x(1) - x(2);
    if ~isfinite(r) || r <= 0, return; end     % r: the cooperation margin, positive in a member
    xn = norm(x);

    params_g = params; params_g.gradedC2 = true;
    p0 = membership_test(W, params_g);

    K = cfg.K;
    dWv = nan(K,1); slot = nan(K,1); isTog = false(K,1);
    ejected = false(K,1); lostR = false(K,1);
    wr_m = nan(K,1); wr_slot_m = nan(K,1); minrM_m = nan(K,1); dr_meas = nan(K,1);
    isCssM = false(K,1);
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
        dWv(t) = Wm(g1,g2) - W(g1,g2); slot(t) = q; isTog(t) = tog;

        pm = membership_test(Wm, params_g);
        isCssM(t) = pm.is_css;
        ejected(t) = ~pm.is_css;
        wr_m(t) = pm.worst_retal; wr_slot_m(t) = pm.worst_retal_slot;
        minrM_m(t) = pm.min_r_M;

        [~, fm, ~] = simulate_pair(Msq, Wm, Kv, Msq, Wm, Kv, sim_params);
        rm = fm(1) - fm(2);
        lostR(t) = ~(rm > 0);
        dr_meas(t) = rm - r;
    end

    % Members kept (ok): ||x*|| in [0.5, 14] and cooperation margin r >= 1e-6
    ok = (xn >= 0.5) && (xn <= 14) && (r >= 1e-6);
    rec = struct('rep',0,'gen',0,'count',0,'N',0, 'ok', ok, ...
        'r', r, 'xn', xn, ...
        'wr0', p0.worst_retal, 'wr0_slot', p0.worst_retal_slot, ...
        'minrM0', p0.min_r_M, 'is_css0', p0.is_css, ...
        'dW', dWv, 'slot', slot, 'isToggle', isTog, ...
        'ejected', ejected, 'lostReadout', lostR, 'is_css_m', isCssM, ...
        'wr_m', wr_m, 'wr_slot_m', wr_slot_m, 'minrM_m', minrM_m, ...
        'dr_meas', dr_meas);
end


% =========================================================================
function repro = check_reproduction(R_all, modeNames, cfg)
% Reproduction check: row order, per-member r and ||x*||, clone counts, ok, and the draws of `dW`,
% `slot`, `isToggle`, `ejected` and `lostReadout` must equal the stored census exactly. The
% census scored `ejected` (the mutation ends membership) with gradedC2 off, stopping at the first
% mutant that earns more than the clone;
% this run scores the whole grid, so equality also checks that the two give the same membership.
    L = load(cfg.srcFile, 'W_all', 'modeNames');
    nModes = numel(modeNames);
    L.W_all = subset_stored(L.W_all, cfg);
    repro = struct('mode',{},'nStored',{},'nReplay',{},'ok',{},'detail',{});
    fprintf('\n================ REPRODUCTION CHECK ================\n');
    allok = true;
    for m = 1:nModes
        A = L.W_all{m}; B = R_all{m};
        d = struct(); okm = true;
        d.nStored = numel(A); d.nReplay = numel(B);
        if d.nStored ~= d.nReplay
            okm = false;
            fprintf('%-22s ROW COUNT MISMATCH stored %d replay %d\n', modeNames{m}, d.nStored, d.nReplay);
        else
            Ku = cfg.K;
            cmp = {'r','xn'};
            for i = 1:numel(cmp)
                a = double([A.(cmp{i})])'; b = double([B.(cmp{i})])';
                d.(['max_' cmp{i}]) = max(abs(a-b));
                okm = okm && isequal(a,b);
            end
            % per-draw fields; members in rows, draws in columns
            fld = {'dW','slot','isToggle','ejected','lostReadout'};
            for i = 1:numel(fld)
                MA = reshape(vertcat(A.(fld{i})), cfg.K, [])'; MA = MA(:,1:Ku);
                MB = reshape(vertcat(B.(fld{i})), cfg.K, [])'; MB = MB(:,1:Ku);
                same = isequal(double(MA), double(MB));
                d.(['eq_' fld{i}]) = same;
                d.(['nDiff_' fld{i}]) = sum(double(MA(:)) ~= double(MB(:)));
                okm = okm && same;
            end
            % clone weights and the ok flag
            okm = okm && isequal(double([A.count])', double([B.count])');
            okm = okm && isequal(logical([A.ok])', logical([B.ok])');
            fprintf(['%-22s n = %d  r/xn identical %d  dW %d (%d diff)  slot %d  toggle %d  ' ...
                     'ends membership %d (%d diff)  lostReadout %d (%d diff)\n'], ...
                modeNames{m}, d.nStored, d.max_r == 0 && d.max_xn == 0, ...
                d.eq_dW, d.nDiff_dW, d.eq_slot, d.eq_isToggle, ...
                d.eq_ejected, d.nDiff_ejected, d.eq_lostReadout, d.nDiff_lostReadout);
        end
        allok = allok && okm;
        repro(m).mode = modeNames{m}; repro(m).nStored = d.nStored;
        repro(m).nReplay = d.nReplay; repro(m).ok = okm; repro(m).detail = d;
    end
    fprintf('REPRODUCTION: %s\n', ternary(allok, 'BIT FOR BIT', 'FAILED -- no downstream number is quotable'));
end


% =========================================================================
function analyse_and_report(R_all, modeNames, repro, cfg)
    nModes = numel(modeNames);
    L = load(cfg.srcFile, 'W_all');
    L.W_all = subset_stored(L.W_all, cfg);
    Ku = cfg.K;

    % Join the stored rho_dir onto the replayed margins, per architecture and pooled. Stored
    % fields: `ejected` = the mutation ends membership, `lostReadout` = it ends cooperation; C2 =
    % ends membership with cooperation kept. Members in rows, draws in columns.
    G = struct('name',{},'RHO',{},'WRM',{},'DWR',{},'DRM',{},'DRP',{}, ...
               'ENDS_MEM_C',{},'LOST',{},'C2',{},'POS',{},'POSGE',{},'TOG',{}, ...
               'w',{},'nAll',{},'nKept',{},'wr0',{});
    for m = 1:nModes
        A = L.W_all{m}; B = R_all{m};
        sel = logical([B.ok])';
        Ak = A(sel); Bk = B(sel);
        rho  = reshape(vertcat(Ak.rho_dir), cfg.K, [])'; rho = rho(:,1:Ku);
        wrm  = reshape(vertcat(Bk.wr_m),    cfg.K, [])'; wrm = wrm(:,1:Ku);
        endsMemC  = reshape(vertcat(Bk.ejected), cfg.K, [])'; endsMemC = logical(endsMemC(:,1:Ku));
        lst  = reshape(vertcat(Bk.lostReadout), cfg.K, [])'; lst = logical(lst(:,1:Ku));
        drm  = reshape(vertcat(Bk.dr_meas), cfg.K, [])'; drm = drm(:,1:Ku);
        tog  = reshape(vertcat(Bk.isToggle), cfg.K, [])'; tog = logical(tog(:,1:Ku));
        wr0  = double([Bk.wr0])'; rpar = double([Bk.r])';
        G(m).name = modeNames{m};
        G(m).RHO = rho; G(m).WRM = wrm;
        G(m).DWR = wrm - repmat(wr0, 1, Ku);
        G(m).DRM = drm;
        G(m).DRP = -rho .* repmat(rpar, 1, Ku);   % first-order numerator of rho
        G(m).ENDS_MEM_C = endsMemC; G(m).LOST = lst;
        G(m).C2  = endsMemC & ~lst;
        G(m).POS   = isfinite(wrm) & wrm >  0;    % (C2) fails: membership_test's rule
        G(m).POSGE = isfinite(wrm) & wrm >= 0;
        G(m).TOG = tog;
        G(m).w   = double([Bk.count])';
        G(m).wr0 = wr0;
        G(m).nAll = numel(B); G(m).nKept = numel(Bk);
    end
    P = struct('name','pooled', ...
        'RHO',vertcat(G.RHO),'WRM',vertcat(G.WRM),'DWR',vertcat(G.DWR), ...
        'DRM',vertcat(G.DRM),'DRP',vertcat(G.DRP),'ENDS_MEM_C',vertcat(G.ENDS_MEM_C), ...
        'LOST',vertcat(G.LOST),'C2',vertcat(G.C2),'POS',vertcat(G.POS), ...
        'POSGE',vertcat(G.POSGE),'TOG',vertcat(G.TOG),'w',vertcat(G.w), ...
        'nAll',sum([G.nAll]),'nKept',sum([G.nKept]),'wr0',vertcat(G.wr0));
    GG = [G, P]; nG = numel(GG);

    fprintf('\n================ SAMPLE ================\n');
    fprintf('member wirings: %s  (pooled %d); kept by the ok filter: %s (pooled %d)\n', ...
        mat2str([G.nAll]), P.nAll, mat2str([G.nKept]), P.nKept);
    fprintf('draws per wiring: %d; total draws %s (pooled %d)\n', ...
        Ku, mat2str([G.nKept]*Ku), P.nKept*Ku);
    fprintf('parent worst_retal finite: %s of %s wirings\n', ...
        mat2str(arrayfun(@(g) sum(isfinite(g.wr0)), G)), mat2str([G.nKept]));
    for g = 1:nG
        wr = draw_weights(GG(g).w, size(GG(g).RHO,2));
        f = isfinite(GG(g).WRM(:));
        fprintf('%-22s mutant margin finite (mutant passed C1) on %d of %d draws (%.4f raw, %.4f clone-weighted)\n', ...
            GG(g).name, sum(f), numel(f), mean(f), wmean(f, wr));
    end

    % Check: a mutant margin above zero is the same event as ending membership with cooperation
    % kept, so the two must agree on every draw
    fprintf('\n================ CHECK: mutant margin > 0 against the stored (C2) loss ================\n');
    fprintf(['A mutant margin above zero is the same event as a mutation that ends only retaliation. ' ...
             'Agreement must be 1.000.\n']);
    for g = 1:nG
        a = GG(g).C2(:); b = GG(g).POS(:); c = GG(g).POSGE(:);
        fprintf('%-22s agree(>0) %.6f (%d of %d disagree) | agree(>=0) %.6f (%d disagree) | k C2 %d, k pos %d\n', ...
            GG(g).name, mean(a==b), sum(a~=b), numel(a), mean(a==c), sum(a~=c), sum(a), sum(b));
    end
    lostNotFinite = sum(P.LOST(:) & isfinite(P.WRM(:)));
    fprintf('draws with readout lost but a finite margin: %d (must be 0: membership_test returns before scoring C2)\n', ...
        lostNotFinite);

    % Per rho bin, clone-weighted: pooled and per architecture
    fprintf('\n================ 1-2. PER rho BIN, POOLED (clone-weighted) ================\n');
    T = bin_table(P, cfg.binEdges, cfg, true);
    fprintf('\n--- bins: edges %g:%g:%g\n', ...
        cfg.binEdges(1), cfg.binEdges(2)-cfg.binEdges(1), cfg.binEdges(end));
    print_bin_table(T);

    Tm = cell(nModes,1);
    for m = 1:nModes
        fprintf('\n--- %s\n', modeNames{m});
        Tm{m} = bin_table(G(m), cfg.binEdges, cfg, true);
        print_bin_table(Tm{m});
    end

    % Spearman correlation and clone-weighted slope of the margin change on dr
    rng(cfg.rngSeedStat);
    fprintf('\n================ 3. d(worst_retal) AGAINST dr ================\n');
    fprintf('dr = measured r_mutant - r_parent on the clone pair at 600 steps.\n');
    Sp = relation_table(GG, 'DRM', cfg);
    fprintf('\nOther dr = the first-order numerator of rho, -rho_dir * r (robustness line).\n');
    Sp2 = relation_table(GG, 'DRP', cfg);

    reproOK = all([repro.ok]);
    if ~reproOK
        fprintf('\n*** REPRODUCTION CHECK FAILED: the replay does not match the stored census; do not use the numbers above. ***\n');
    end

    unused_split(GG, cfg);

    out = struct('T',T,'Tm',{Tm},'Sp',Sp,'Sp2',Sp2,'reproOK',reproOK, ...
        'modeNames',{modeNames},'cfg',cfg); %#ok<NASGU>
    save(fullfile(cfg.dataDir,'retaliation_margin_vs_rho_summary.mat'), 'out');
    fprintf('\nDone.\n');
end


% =========================================================================
function T = bin_table(S, edges, cfg, ~)
% Per rho bin, clone-weighted: draws, Kish n_eff, members, the share with a finite mutant margin,
% the margin and its change (median, IQR, range), and the shares ending membership and cooperation.
    Ku = size(S.RHO,2);
    w = draw_weights(S.w, Ku);
    rho = S.RHO(:); wrm = S.WRM(:); dwr = S.DWR(:);
    c2 = double(S.C2(:)); pos = double(S.POS(:)); endsMemC = double(S.ENDS_MEM_C(:)); lst = double(S.LOST(:));
    nb = numel(edges) - 1;
    idx = discretize(rho, edges);
    T = struct();
    T.edges = edges;
    T.bC = 0.5*(edges(1:end-1) + edges(2:end));
    z = nan(nb,1);
    T.n = z; T.nEff = z; T.nWir = z; T.nFinite = z; T.fracFinite = z;
    T.wrm_med = z; T.wrm_q1 = z; T.wrm_q3 = z; T.wrm_lo = z; T.wrm_hi = z;
    T.dwr_med = z; T.dwr_q1 = z; T.dwr_q3 = z; T.dwr_lo = z; T.dwr_hi = z;
    T.pC2 = z; T.pC2_lo = z; T.pC2_hi = z; T.pPos = z;
    T.pEndsMem = z; T.pLost = z; T.kC2 = z;
    idxW = discretize(S.RHO, edges);
    for b = 1:nb
        mm = (idx == b);
        if ~any(mm), continue; end
        T.n(b) = sum(mm);
        sw = sum(w(mm)); T.nEff(b) = sw^2 / max(sum(w(mm).^2), eps);
        T.nWir(b) = sum(any(idxW == b, 2));
        f = mm & isfinite(wrm);
        T.nFinite(b) = sum(f); T.fracFinite(b) = wmean(double(isfinite(wrm(mm))), w(mm));
        if any(f)
            T.wrm_med(b) = wq(wrm(f), w(f), 0.5);
            T.wrm_q1(b)  = wq(wrm(f), w(f), 0.25);
            T.wrm_q3(b)  = wq(wrm(f), w(f), 0.75);
            T.wrm_lo(b)  = min(wrm(f)); T.wrm_hi(b) = max(wrm(f));
            fd = mm & isfinite(dwr);
            T.dwr_med(b) = wq(dwr(fd), w(fd), 0.5);
            T.dwr_q1(b)  = wq(dwr(fd), w(fd), 0.25);
            T.dwr_q3(b)  = wq(dwr(fd), w(fd), 0.75);
            T.dwr_lo(b)  = min(dwr(fd)); T.dwr_hi(b) = max(dwr(fd));
        end
        T.pC2(b) = wmean(c2(mm), w(mm));
        T.pPos(b) = wmean(pos(mm), w(mm));
        T.pEndsMem(b) = wmean(endsMemC(mm), w(mm));
        T.pLost(b) = wmean(lst(mm), w(mm));
        T.kC2(b) = sum(c2(mm));
        [T.pC2_lo(b), T.pC2_hi(b)] = wilson(T.pC2(b), T.nEff(b));
    end
end


function print_bin_table(T)
    fprintf(['  %7s %8s %9s %8s %8s | %-28s | %-28s | %8s %8s %8s\n'], ...
        'rho','n draws','n_eff','n wir','fin frac', ...
        'mutant margin  med [IQR]', 'd margin  med [IQR]', 'P(C2)','P(m>0)','P(ends)');
    for b = 1:numel(T.bC)
        if ~isfinite(T.n(b)) || T.n(b) == 0, continue; end
        fprintf(['  %7.3f %8d %9.1f %8d %8.4f | %9.4f [%8.4f,%8.4f] | %9.4f [%8.4f,%8.4f] | ' ...
                 '%8.4f %8.4f %8.4f\n'], ...
            T.bC(b), T.n(b), T.nEff(b), T.nWir(b), T.fracFinite(b), ...
            T.wrm_med(b), T.wrm_q1(b), T.wrm_q3(b), ...
            T.dwr_med(b), T.dwr_q1(b), T.dwr_q3(b), ...
            T.pC2(b), T.pPos(b), T.pEndsMem(b));
    end
    fprintf('  ranges (min..max): mutant margin %.4g .. %.4g | d margin %.4g .. %.4g\n', ...
        min(T.wrm_lo), max(T.wrm_hi), min(T.dwr_lo), max(T.dwr_hi));
end


% =========================================================================
function Sp = relation_table(GG, drField, cfg)
% Spearman and the clustered slope of d(worst_retal) on dr, per group.
    nG = numel(GG);
    blank = struct('name','','n',0,'nWir',0,'rho_s',NaN,'rho_s_lo',NaN,'rho_s_hi',NaN, ...
                   'rho_shuf',NaN,'slope',NaN,'slope_lo',NaN,'slope_hi',NaN,'slope_unw',NaN, ...
                   'dr_med',NaN,'dr_q1',NaN,'dr_q3',NaN,'dr_lo',NaN,'dr_hi',NaN, ...
                   'dwr_med',NaN,'dwr_q1',NaN,'dwr_q3',NaN,'dwr_lo',NaN,'dwr_hi',NaN);
    Sp = repmat(blank, 1, nG);
    fprintf('%-22s %9s %8s %9s %20s %9s %11s %22s %11s\n', ...
        'group','n pairs','n wir','Spearman','bootstrap 95%','shuffled','slope (w)','bootstrap 95%','slope (unw)');
    for g = 1:nG
        DR = GG(g).(drField); DW = GG(g).DWR; Ku = size(DR,2);
        w = draw_weights(GG(g).w, Ku);
        ok = isfinite(DR(:)) & isfinite(DW(:));
        x = DR(ok); y = DW(ok); ww = w(ok);
        Sp(g).name = GG(g).name;
        Sp(g).n = double(sum(ok));
        Sp(g).nWir = double(sum(any(isfinite(DR) & isfinite(DW), 2)));
        if Sp(g).n < 10
            fprintf('%-22s %9d  (fewer than 10 usable pairs -- not estimable)\n', GG(g).name, Sp(g).n);
            continue;
        end
        Sp(g).rho_s = scal(corr(x, y, 'Type', 'Spearman'));
        Sp(g).slope = scal(wls_slope(x, y, ww));
        Sp(g).slope_unw = scal(wls_slope(x, y, ones(size(x))));
        Sp(g).dr_med = scal(wq(x, ww, 0.5)); Sp(g).dr_q1 = scal(wq(x, ww, 0.25)); Sp(g).dr_q3 = scal(wq(x, ww, 0.75));
        Sp(g).dr_lo = min(x); Sp(g).dr_hi = max(x);
        Sp(g).dwr_med = scal(wq(y, ww, 0.5)); Sp(g).dwr_q1 = scal(wq(y, ww, 0.25)); Sp(g).dwr_q3 = scal(wq(y, ww, 0.75));
        Sp(g).dwr_lo = min(y); Sp(g).dwr_hi = max(y);

        % wiring-clustered bootstrap: whole wirings, all K draws carried
        nW = size(DR,1); cb = nan(cfg.nBoot,1); sb = nan(cfg.nBoot,1);
        for t = 1:cfg.nBoot
            gi = randi(nW, nW, 1);
            Xr = DR(gi,:); Yr = DW(gi,:); Wr = repmat(GG(g).w(gi), Ku, 1);
            o = isfinite(Xr(:)) & isfinite(Yr(:));
            if sum(o) < 10, continue; end
            xr = Xr(:); yr = Yr(:);
            cb(t) = scal(corr(xr(o), yr(o), 'Type','Spearman'));
            sb(t) = scal(wls_slope(xr(o), yr(o), Wr(o)));
        end
        [Sp(g).rho_s_lo, Sp(g).rho_s_hi] = pctlCI(cb);
        [Sp(g).slope_lo, Sp(g).slope_hi] = pctlCI(sb);

        % control: shuffle dr WITHIN wiring (the K draws permuted among themselves)
        Xs = DR;
        for i = 1:nW, Xs(i,:) = DR(i, randperm(Ku)); end
        os = isfinite(Xs(:)) & isfinite(DW(:));
        if sum(os) >= 10
            xs = Xs(:); ys = DW(:);
            Sp(g).rho_shuf = scal(corr(xs(os), ys(os), 'Type','Spearman'));
        end

        fprintf('%-22s %9d %8d %9.4f [%8.4f, %8.4f] %9.4f %11.4g [%9.4g, %9.4g] %11.4g\n', ...
            Sp(g).name, Sp(g).n, Sp(g).nWir, Sp(g).rho_s, Sp(g).rho_s_lo, Sp(g).rho_s_hi, ...
            Sp(g).rho_shuf, Sp(g).slope, Sp(g).slope_lo, Sp(g).slope_hi, Sp(g).slope_unw);
    end
    fprintf('\n%-22s %-44s %-44s\n', 'group', 'dr  med [IQR] (range)', 'd worst_retal  med [IQR] (range)');
    for g = 1:nG
        if ~isfinite(Sp(g).dr_med), continue; end
        fprintf('%-22s %9.4g [%9.4g,%9.4g] (%8.3g..%8.3g) %9.4g [%9.4g,%9.4g] (%8.3g..%8.3g)\n', ...
            Sp(g).name, Sp(g).dr_med, Sp(g).dr_q1, Sp(g).dr_q3, Sp(g).dr_lo, Sp(g).dr_hi, ...
            Sp(g).dwr_med, Sp(g).dwr_q1, Sp(g).dwr_q3, Sp(g).dwr_lo, Sp(g).dwr_hi);
    end
end


function v = scal(x)
    if isempty(x), v = NaN; else, v = double(x(1)); end
end

function [lo, hi] = pctlCI(v)
    v = v(isfinite(v));
    if numel(v) < 20, lo = NaN; hi = NaN; return; end
    q = quantile(v(:), [0.025 0.975]);
    lo = q(1); hi = q(2);
end


function s = wls_slope(x, y, w)
    w = w(:) / sum(w);
    mx = sum(w.*x); my = sum(w.*y);
    vx = sum(w.*(x-mx).^2);
    if vx <= 0, s = NaN; return; end
    s = sum(w.*(x-mx).*(y-my)) / vx;
end


% =========================================================================
function unused_split(GG, cfg)
% Draws split at |dr| < cfg.deadDr: dead (the drawn weight has no path to the readout, so the
% mutation moves nothing) against live; per group and per pooled rho bin, the dead share and the
% live margin change. The split is on dr, never on the response.
    if ~isfield(cfg,'deadDr'), cfg.deadDr = 1e-6; end
    fprintf('\n================ the dead-slot split ================\n');
    fprintf(['Dead = |dr| < %g (predictor side; the rectifier leaves the drawn slot with no path\n' ...
             'to the readout). Never split on the response.\n'], cfg.deadDr);
    fprintf('%-22s %10s %12s %12s %34s %34s\n', 'group','n draws','dead frac','live n', ...
        'd worst_retal, LIVE  med [IQR]', 'Spearman LIVE [shuffled]');
    for g = 1:numel(GG)
        Ku = size(GG(g).DRM,2); w = draw_weights(GG(g).w, Ku);
        dr = GG(g).DRM(:); dw = GG(g).DWR(:);
        ok = isfinite(dr) & isfinite(dw);
        live = ok & abs(dr) >= cfg.deadDr;
        dfrac = 1 - wmean(double(live(ok)), w(ok));
        rl = NaN; rs = NaN;
        if sum(live) > 10
            rl = corr(dr(live), dw(live), 'Type','Spearman');
            DRl = GG(g).DRM; nW = size(DRl,1); Xs = DRl;
            for i = 1:nW, Xs(i,:) = DRl(i, randperm(Ku)); end
            xs = Xs(:); os = isfinite(xs) & isfinite(dw) & abs(xs) >= cfg.deadDr;
            if sum(os) > 10, rs = corr(xs(os), dw(os), 'Type','Spearman'); end
        end
        fprintf('%-22s %10d %12.4f %12d %10.4f [%9.4f,%9.4f] %14.4f [%9.4f]\n', ...
            GG(g).name, sum(ok), dfrac, sum(live), ...
            wq(dw(live), w(live), 0.5), wq(dw(live), w(live), 0.25), wq(dw(live), w(live), 0.75), ...
            rl, rs);
    end
    % Per pooled rho bin: dead fraction and the live median.
    P = GG(end); Ku = size(P.DRM,2); w = draw_weights(P.w, Ku);
    rho = P.RHO(:); dr = P.DRM(:); dw = P.DWR(:);
    idx = discretize(rho, cfg.binEdges);
    bC = 0.5*(cfg.binEdges(1:end-1) + cfg.binEdges(2:end));
    % A spike at zero plus a heavy tail: the weighted shares above 0.01 and below -0.01 and the
    % mean are printed beside the median.
    fprintf('\n%-9s %9s %10s %9s %32s %9s %9s %9s %13s\n', 'rho','n usable','dead frac','live n', ...
        'd worst_retal, LIVE  med [IQR]', 'mean', 'P(d>.01)', 'P(d<-.01)', 'live range');
    for b = 1:numel(bC)
        mm = (idx == b) & isfinite(dr) & isfinite(dw);
        if sum(mm) < cfg.minBinN, continue; end
        live = mm & abs(dr) >= cfg.deadDr;
        fprintf('%9.3f %9d %10.4f %9d %9.4f [%8.4f,%8.4f] %9.4f %9.4f %9.4f %6.3g..%-6.3g\n', ...
            bC(b), sum(mm), 1 - wmean(double(live(mm)), w(mm)), sum(live), ...
            wq(dw(live), w(live), 0.5), wq(dw(live), w(live), 0.25), wq(dw(live), w(live), 0.75), ...
            wmean(dw(mm), w(mm)), wmean(double(dw(mm) > 0.01), w(mm)), ...
            wmean(double(dw(mm) < -0.01), w(mm)), min(dw(live)), max(dw(live)));
    end
end


% =========================================================================
function wr = draw_weights(w, K)
% Clone weight per draw, in the order of M(:) for an nWiring x K matrix M (the wiring index runs
% fastest), so each draw carries its own wiring's weight.
    wr = repmat(w(:), K, 1);
end

function m = wmean(y, w)
    y = double(y(:)); w = double(w(:));
    ok = isfinite(y) & isfinite(w) & w > 0;
    if ~any(ok), m = NaN; return; end
    m = sum(w(ok).*y(ok)) / sum(w(ok));
end

function [lo, hi] = wilson(ph, n)
    if ~isfinite(ph) || ~isfinite(n) || n <= 0, lo = NaN; hi = NaN; return; end
    z = 1.959963984540054;
    den = 1 + z^2/n;
    c = (ph + z^2/(2*n)) / den;
    h = (z/den) * sqrt(ph*(1-ph)/n + z^2/(4*n^2));
    lo = max(0, c-h); hi = min(1, c+h);
end

function Wsub = subset_stored(W_all, cfg)
% A partial replay (cfg.nReps < 100) covers only the leading replicates; the stored rows are cut to
% those replicates so that the check and the join stay row-aligned. At 100 this changes nothing.
    Wsub = W_all;
    for m = 1:numel(W_all)
        rr = double([W_all{m}.rep])';
        Wsub{m} = W_all{m}(rr <= cfg.nReps);
    end
end
