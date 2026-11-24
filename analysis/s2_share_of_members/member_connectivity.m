function member_connectivity(cfg)
% MEMBER_CONNECTIVITY  Are members joined by chains of single-weight changes that never leave the
% conditional genotype network? For pairs of members of one architecture from different
% populations, searches for such a path and prints the share of pairs joined, a lower bound.
%
% Produces : the share of member pairs joined, per architecture and pooled, with Clopper-Pearson
%            95% intervals, printed under "adaptive-order search" (60 pairs per architecture;
%            SI S2.3; Table S3, "Member pairs joined")
% Inputs   : publication/<architecture>_100reps.mat   (run_signaling_architectures; members are
%              sampled as per_mutation_census samples them: 5 stored generations from 100 on)
%            cfg   optional struct of settings; the defaults at the top of the code give the
%                  paper's numbers; nPairsSearch (pairs per architecture) takes a smaller sample
% Outputs  : publication/member_connectivity.mat   (the search results and the member counts);
%            console summary
% Run      : member_connectivity   (parfor on a local pool; about 5 min on 3 workers)
% Seed     : rng(20260913) (MATLAB's default generator) before the pair draws; each pair's search
%            has its own RandStream('mt19937ar') seeded from that stream; membership_test seeds
%            each wiring itself
%
% Pairs are drawn among distinct members with ||x*|| (the size of the settled expression state)
% in [0.5, 14] and a cooperation margin of at least 1e-6. A step changes one weight by at most 0.6,
% the half-width of the membership test's grid; toggles are left out. A path counts when every
% wiring on it passes the membership test. A pair not joined may be a search failure, so every
% share is a lower bound.

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    % dataDir: the data root's publication folder
    if ~isfield(cfg,'dataDir') || isempty(cfg.dataDir)
        cfg.dataDir = fullfile(grn_data_file('publication'));
    end
    if ~isfield(cfg,'outDir') || isempty(cfg.outDir), cfg.outDir = cfg.dataDir; end
    if ~isfield(cfg,'nReps'),           cfg.nReps = 100;        end
    if ~isfield(cfg,'nSnaps'),          cfg.nSnaps = 5;         end
    if ~isfield(cfg,'steadyStateFrom'), cfg.steadyStateFrom = 100; end
    if ~isfield(cfg,'nPairsSearch'),    cfg.nPairsSearch = 60;  end
    if ~isfield(cfg,'maxEval'),         cfg.maxEval = 300;      end
    if ~isfield(cfg,'stepMax'),         cfg.stepMax = 0.6;      end   % 3*sigma_mut
    if ~isfield(cfg,'sigmaMut'),        cfg.sigmaMut = 0.2;     end
    if ~isfield(cfg,'developmentSteps'),cfg.developmentSteps = 600; end
    if ~isfield(cfg,'minXnorm'),        cfg.minXnorm = 0.5;     end
    if ~isfield(cfg,'maxXnorm'),        cfg.maxXnorm = 14;      end
    if ~isfield(cfg,'minR'),            cfg.minR = 1e-6;        end
    if ~isfield(cfg,'rngSeed'),         cfg.rngSeed = 20260913; end
    if ~isfield(cfg,'useParallel'),     cfg.useParallel = true; end

    modes = { ...
        struct('name','graded_separate',     'sensingMode','separate',     'hillCoefficient',Inf,'discreteOutput',false), ...
        struct('name','graded_contrast',     'sensingMode','contrast', 'hillCoefficient',Inf,'discreteOutput',false), ...
        struct('name','saturating_separate', 'sensingMode','separate',     'hillCoefficient',1,  'discreteOutput',false), ...
        struct('name','saturating_contrast', 'sensingMode','contrast', 'hillCoefficient',1,  'discreteOutput',false)};
    nModes = numel(modes);
    modeNames = cell(nModes,1);
    for m = 1:nModes, modeNames{m} = modes{m}.name; end
    P = cell(nModes,1);
    for m = 1:nModes, P{m} = membership_params(modes{m}, cfg); end

    if cfg.useParallel && isempty(gcp('nocreate')), parpool('local'); end
    totalStart = tic;

    MEM = extract_members(modes, P, cfg);
    SRCH = run_adaptive(MEM, P, modeNames, cfg);

    counts = struct();
    counts.distinctClassified = cellfun(@(z) z.nDistinctClassified, MEM);
    counts.members            = cellfun(@(z) size(z.W,1), MEM);
    counts.ok                 = cellfun(@(z) sum(z.meta(:,6)), MEM);
    counts.distinctW          = cellfun(@(z) size(unique(z.W,'rows'),1), MEM);
    counts.distinctWok        = cellfun(@(z) z.nDistinctOk, MEM);
    metadata = struct('cfg', cfg, 'modes', {modeNames}, 'generator', mfilename, ...
        'runDate', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
        'elapsedMin', toc(totalStart)/60);
    save(fullfile(cfg.outDir,'member_connectivity.mat'), ...
        'SRCH','counts','modeNames','metadata','-v7.3');
    fprintf('\nsaved %s\n', fullfile(cfg.outDir,'member_connectivity.mat'));

    report(SRCH, counts, modeNames);
    fprintf('\ntotal %.1f min\n', toc(totalStart)/60);
end

% ------------------------------------------------------------------ params ---
function p = membership_params(mode, cfg)
    p = struct('sensingMode',mode.sensingMode,'hillCoefficient',mode.hillCoefficient, ...
        'discreteOutput',mode.discreteOutput,'michaelisConstant',1.0,'degradationRate',0.1, ...
        'timeStep',1.0,'developmentSteps',cfg.developmentSteps,'backend','simulate_pair', ...
        'sigma_mut',cfg.sigmaMut);
end

% ----------------------------------------------------------- member wirings ---
function MEM = extract_members(modes, P, cfg)
% The census's sample, as per_mutation_census takes it: nSnaps snapshots per replicate at
% generation >= steadyStateFrom, every distinct wiring through the membership test; per member the
% replicate, generation, clone count, r (the cooperation margin), ||x*|| and ok (kept: ||x*|| in
% [minXnorm, maxXnorm], r >= minR). No random draws (membership_test restores the caller's stream).
    nG = 4; MEM = cell(numel(modes),1);
    for m = 1:numel(modes)
        mn = modes{m}.name;
        sp = P{m}; sp.payoffDC = 0; sp.sensorIntegrationRate = 1.0;
        sp.externalSignal1 = zeros(nG,1); sp.externalSignal2 = zeros(nG,1);
        L = load(fullfile(cfg.dataDir,sprintf('%s_100reps.mat',mn)),'results');
        R = L.results; clear L;
        nR = min(cfg.nReps, numel(R));
        snaps = cell(nR,1);
        for r = 1:nR
            rep = R{r}; sg = double(rep.historyGenerations(:))';
            ks = find(sg >= cfg.steadyStateFrom);
            ks = ks(round(linspace(1,numel(ks),cfg.nSnaps)));
            S = cell(numel(ks),1);
            for s = 1:numel(ks)
                k = ks(s);
                M = double(rep.historyM(:,:,:,k)); B = double(rep.historyB(:,:,:,k));
                N = size(M,3); Wall = reshape(M.*B, nG*nG, N)';
                [Wu,~,ic] = unique(Wall,'rows');
                S{s} = struct('gen',sg(k),'Wu',Wu,'counts',accumarray(ic,1),'N',N);
            end
            snaps{r} = S;
        end
        clear R;
        t0 = tic; rows = cell(nR,1); nDist = zeros(nR,1);
        params = P{m}; minXn = cfg.minXnorm; maxXn = cfg.maxXnorm; minR = cfg.minR;
        parfor r = 1:nR
            S = snaps{r}; W = []; meta = []; nd = 0;
            for s = 1:numel(S)
                Wu = S{s}.Wu; nd = nd + size(Wu,1);
                for u = 1:size(Wu,1)
                    Wm = reshape(Wu(u,:),nG,nG);
                    if ~membership_test(Wm, params).is_css, continue; end
                    [~,fs,~] = simulate_pair(ones(nG),Wm,[],ones(nG),Wm,[],sp);
                    x = fs(1:nG); rr = x(1)-x(2); xn = norm(x);
                    ok = (xn >= minXn) && (xn <= maxXn) && (rr >= minR);
                    W = [W; Wu(u,:)];
                    meta = [meta; r, S{s}.gen, S{s}.counts(u), rr, xn, double(ok)];
                end
            end
            rows{r} = struct('W',W,'meta',meta); nDist(r) = nd;
        end
        tW = cellfun(@(z) z.W, rows,'uni',0); Wc = vertcat(tW{:});
        tM = cellfun(@(z) z.meta, rows,'uni',0); Mc = vertcat(tM{:});
        nOk = size(unique(Wc(Mc(:,6)==1,:),'rows'),1);
        MEM{m} = struct('name',mn,'W',Wc,'meta',Mc,'nDistinctClassified',sum(nDist), ...
                        'nDistinctOk',nOk);
        fprintf('%-22s membership_tested=%6d  members=%6d  ok=%6d  distinctW=%6d  distinctW_ok=%6d  %.1f min\n', ...
            mn, sum(nDist), size(Wc,1), sum(Mc(:,6)), size(unique(Wc,'rows'),1), nOk, toc(t0)/60);
    end
end

function [Wu, ru] = ok_members(A)
% The distinct ok members and the replicate each was first seen in. The ordering
% is what the pair draws index into, so it must not change.
    ok = A.meta(:,6) == 1;
    [Wu,ia] = unique(A.W(ok,:),'rows');
    mt = A.meta(ok,:); ru = mt(ia,1);
end

% ---------------------------------------------------------- adaptive order ---
function SRCH = run_adaptive(MEM, P, modeNames, cfg)
    fprintf('\n=== adaptive-order search, across-replicate pairs, budget %d tests/pair ===\n', cfg.maxEval);
    rng(cfg.rngSeed);
    SRCH = cell(numel(MEM),1);
    for m = 1:numel(MEM)
        A = MEM{m}; [Wu,ru] = ok_members(A); n = size(Wu,1);
        pairs = zeros(0,2);
        while size(pairs,1) < cfg.nPairsSearch
            i = randi(n); j = randi(n);
            if i == j || ru(i) == ru(j), continue; end
            pairs(end+1,:) = [i j]; %#ok<AGROW>
        end
        np = cfg.nPairsSearch;
        seeds = randi(2^31-2, np, 1); res = zeros(np,4); t0 = tic;
        Pm = P{m}; step = cfg.stepMax; budget = cfg.maxEval;
        parfor q = 1:np
            res(q,:) = search_pair(Wu(pairs(q,1),:), Wu(pairs(q,2),:), Pm, step, budget, seeds(q));
        end
        SRCH{m} = struct('name',modeNames{m},'pairs',pairs,'res',res,'n',np,'secs',toc(t0));
        fprintf('%-22s found: %3d/%3d | moves med=%s | detours med=%s | tests med=%s | %.1f min\n', ...
            modeNames{m}, sum(res(:,1)), np, num2str(md(res(res(:,1)==1,2))), ...
            num2str(md(res(res(:,1)==1,3))), num2str(md(res(res(:,1)==1,4))), toc(t0)/60);
    end
end

function r = search_pair(A, B, P, step, maxEval, seed)
% r = [found nMoves nDetours nTests]. Tries the weights still to move in random order and takes the
% first move toward B (at most step) that stays inside the set; when none does, tries up to 12
% random single-weight steps inside the set (detours). Stops after more than maxDet = 60 detours or
% after maxEval membership tests.
    rs = RandStream('mt19937ar','Seed',seed);
    cur = A; open = find(abs(B-A) > 0); nEval = 0; nMove = 0; nDet = 0; maxDet = 60;
    while nEval < maxEval
        if isempty(open), r = [1 nMove nDet nEval]; return; end
        ord = open(randperm(rs,numel(open))); moved = false;
        for e = ord
            d = B(e)-cur(e); cand = cur;
            if abs(d) <= step, cand(e) = B(e); else, cand(e) = cur(e) + sign(d)*step; end
            nEval = nEval + 1;
            if nEval > maxEval, break; end
            if membership_test(reshape(cand,4,4),P).is_css
                cur = cand; nMove = nMove+1; moved = true;
                if cur(e) == B(e), open(open==e) = []; end
                break;
            end
        end
        if moved, continue; end
        ok = false;
        for t = 1:12
            e = randi(rs,16); dd = (rs.rand()*2-1)*step;
            cand = cur; cand(e) = cur(e)+dd; nEval = nEval+1;
            if nEval > maxEval, break; end
            if membership_test(reshape(cand,4,4),P).is_css
                cur = cand; nMove = nMove+1; nDet = nDet+1; ok = true;
                if cur(e) ~= B(e) && ~any(open==e), open(end+1) = e; end %#ok<AGROW>
                break;
            end
        end
        if ~ok || nDet > maxDet, break; end
    end
    r = [0 nMove nDet nEval];
end

% ----------------------------------------------------------------- report ----
function report(SRCH, counts, modeNames)
    nM = numel(modeNames);
    fprintf('\n=== member wirings (distinct, ok) ===\n');
    for m = 1:nM
        fprintf('%-22s membership_tested=%6d members=%6d ok=%6d distinctW=%6d distinctW_ok=%6d\n', ...
            modeNames{m}, counts.distinctClassified(m), counts.members(m), counts.ok(m), ...
            counts.distinctW(m), counts.distinctWok(m));
    end
    fprintf('%-22s membership_tested=%6d members=%6d ok=%6d distinctW=%6d distinctW_ok=%6d\n','TOTAL', ...
        sum(counts.distinctClassified), sum(counts.members), sum(counts.ok), ...
        sum(counts.distinctW), sum(counts.distinctWok));

    fprintf('\n=== P(a route found), entry order chosen adaptively, across-replicate pairs ===\n');
    tk = 0; tn = 0;
    for m = 1:nM
        R = SRCH{m}.res; n = SRCH{m}.n; k = sum(R(:,1)); [lo,hi] = cp(k,n); s = R(R(:,1)==1,:);
        fprintf('%-22s found: %3d/%3d = %.3f [%.3f %.3f] | moves med=%s (min 16) | detours med=%s | tests med=%s\n', ...
            SRCH{m}.name, k, n, k/n, lo, hi, num2str(md(s(:,2))), num2str(md(s(:,3))), num2str(md(s(:,4))));
        tk = tk+k; tn = tn+n;
    end
    [lo,hi] = cp(tk,tn);
    fprintf('%-22s found: %3d/%3d = %.3f [%.3f %.3f]\n','POOLED',tk,tn,tk/tn,lo,hi);

    fprintf('\nEvery fraction above is a LOWER BOUND: a failure is a search failure.\n');
end

function [lo,hi] = cp(k,n)
% Clopper-Pearson exact 95% interval for k of n
    if k == 0, lo = 0; else, lo = betaincinv(0.025,k,n-k+1); end
    if k == n, hi = 1; else, hi = betaincinv(0.975,k+1,n-k); end
end
function v = md(x), if isempty(x), v = NaN; else, v = median(x); end, end
