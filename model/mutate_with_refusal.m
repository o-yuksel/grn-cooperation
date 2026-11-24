function varargout = mutate_with_refusal(varargin)
% MUTATE_WITH_REFUSAL  The runs' mutation step, with a rule that refuses chosen mutations and
% draws them again; in the paper, a share q of the mutations that end cooperation is refused.
%
% Implements : refusing mutations that end cooperation with the clone (Methods, "Failure route";
%              SI S4.6, "All refused" and "Refusing a random share of these mutations"; Fig. 5).
%              As config.mutationOperator it replaces evolve_population's mutation step with the
%              same Poisson number of mutations per offspring and the same draws; a refused draw
%              is drawn again, up to 20 times, then that mutation is dropped.
% Inputs     : M, B, K, config   the offspring generation, as evolve_population passes it;
%                                config.refusalMode selects the refusal rule (below)
%              or an action name and its arguments ('reset', 'log', 'coop', 'coop_pair', 'gap',
%              'partners', 'mutation_class'), each listed beside the switch that runs it
% Outputs    : [M, B, K] the mutated generation; or the action's outputs
% Seed       : mutation draws use the caller's stream (evolve_population's, per replicate); every
%              cooperation value is computed under rng(20261022), the caller's stream saved and
%              restored around it; the sample of accepted mutations has its own stream (seed
%              20261024), and the coin of 'REFUSE_PARTIAL' its own, seeded by config.refusalSeed
%
% coop(W, P) is W's mean cooperation against P over 100 encounters under the runs' encounter
% rule; a mutation ends cooperation when coop(parent, parent) >= 0.5 and coop(mutant, mutant)
% < 0.5. Refusal rules (config.refusalMode), applied to every draw:
%   'REFUSE_PARTIAL'  the paper's share q refused: each draw that ends cooperation, first draws and
%               redraws alike, is refused when a coin u ~ U(0, 1) falls below q =
%               config.refusalShare (p in the code); q = 1 is all refused. With
%               config.refusalNeverRefuse true the coin is drawn but never refuses.
%   'REFUSE_ALL', 'REFUSE_RESPONSE', 'REFUSE_OTHER'   refuse the mutations that end
%               cooperation by class (mutation_class_end): all of them (all refused), those
%               through the response to the
%               partner, the rest
%   'REFUSE_LESS_CONDITIONAL', 'REFUSE_LESS_RETALIATING'   refuse a mutant less conditional
%               (Gap lower) or retaliating less (Ret higher) than its parent (gap below); 'off'
%               refuses nothing

    persistent LOG CACHE RS BC TS
    if isempty(CACHE), CACHE = containers.Map('KeyType', 'char', 'ValueType', 'double'); end
    if isempty(LOG), LOG = new_log(); end
    if isempty(RS), RS = RandStream('mt19937ar', 'Seed', 20261024); end
    if isempty(BC), BC = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
    if ischar(varargin{1})
        % Actions: 'reset' clears the log and the stored values; L = ('log'): counts of this
        % process's draws, refusals and redraws, and a sample of 50 accepted mutations;
        % c = ('coop', W, P, config): cooperation of W against P (coop below);
        % [c, b, fa, pa] = ('coop_pair', W, P, config): as 'coop', plus P's mean cooperation b and
        % both actions per encounter, fa (W) and pa (P), in the same encounters;
        % g = ('gap', W, S, config): Gap(W) against the partner set S; S = ('partners', M, B,
        % config): the partner set of a population; [cls, cP, cM, cPo, cMo] = ('mutation_class', Wp, Wm,
        % config): the class of the mutation Wp -> Wm (mutation_class_end below)
        switch varargin{1}
            case 'reset'
                LOG = new_log(); CACHE = containers.Map('KeyType', 'char', 'ValueType', 'double');
                RS = RandStream('mt19937ar', 'Seed', 20261024); varargout = {}; TS = [];
                BC = containers.Map('KeyType', 'char', 'ValueType', 'any');
            case 'log', varargout = {LOG};
            case 'coop'
                [W, P, config] = varargin{2:4}; [c, CACHE, LOG] = coop(W, P, config, CACHE, LOG); varargout = {c};
            case 'coop_pair'
                [W, P, config] = varargin{2:4}; [c, b, fa, pa] = coop_pair(W, P, config); varargout = {c, b, fa, pa};
            case 'gap'
                [W, S, config] = varargin{2:4}; [g, CACHE, LOG] = gap(W, S, config, CACHE, LOG); varargout = {g};
            case 'partners'
                [M, B, config] = varargin{2:4}; [S, CACHE, LOG] = partners(M, B, config, CACHE, LOG); varargout = {S};
            case 'mutation_class'
                [Wp, Wm, config] = varargin{2:4}; assert_split_config(config);
                [cP, CACHE, LOG] = coop(Wp, Wp, config, CACHE, LOG);
                [cM, CACHE, LOG] = coop(Wm, Wm, config, CACHE, LOG);
                [cls, v, CACHE, BC, LOG] = mutation_class_end(Wp, Wm, cP, cM, config, CACHE, BC, LOG);
                varargout = {cls, v(1), v(2), v(3), v(4)};
            otherwise
                error('mutate_with_refusal:action', 'unknown action %s', varargin{1});
        end
        return;
    end
    [M, B, K, config] = varargin{1:4};
    assert(~isfield(config, 'geneInputBias') || isempty(config.geneInputBias), 'mutate_with_refusal: run config carries geneInputBias');
    for fn = {'fixedMutualInhibition', 'frozenGenes', 'frozenConnections', 'evolvableKM_genes'}
        if isfield(config, fn{1}) && ~isempty(config.(fn{1})) && ~isequal(config.(fn{1}), false)
            error('mutate_with_refusal:unsupported', 'config.%s is set.', fn{1});
        end
    end
    if ~LOG.started, s = rng; LOG.firstState = double(s.State(1:8))'; LOG.started = true; end
    retOnly = strcmp(config.refusalMode, 'REFUSE_LESS_RETALIATING');
    split = find(strcmp(config.refusalMode, {'REFUSE_ALL', 'REFUSE_RESPONSE', 'REFUSE_OTHER'}));
    refusing_partial = strcmp(config.refusalMode, 'REFUSE_PARTIAL');
    refusing = strcmp(config.refusalMode, 'REFUSE_LESS_CONDITIONAL') || retOnly || ~isempty(split) || refusing_partial;
    assert(refusing || strcmp(config.refusalMode, 'off'), 'refusalMode must be off, REFUSE_LESS_CONDITIONAL, REFUSE_LESS_RETALIATING, REFUSE_ALL, REFUSE_RESPONSE, REFUSE_OTHER or REFUSE_PARTIAL');
    if ~isempty(split)
        assert_split_config(config); refSet = {[1 2 3], 1, [2 3]}; refSet = refSet{split};
    end
    if refusing_partial
        [TS, LOG] = refuse_stream(TS, LOG, config); S = {}; LOG.calls = LOG.calls + 1;
        [M, B, CACHE, LOG] = refuse_mutate(M, B, config, TS, CACHE, LOG, RS);
        varargout = {M, B, K}; return;
    end
    [S, CACHE, LOG] = partners(M, B, config, CACHE, LOG);
    LOG.calls = LOG.calls + 1; LOG.noPop = LOG.noPop + (numel(S) == 1);
    [nGenes, ~, N] = size(M);
    for i = 1:N
        nMutations = poissrnd(config.mutationRate);
        for m = 1:nMutations
            W0 = double(M(:, :, i)) .* double(B(:, :, i));
            [g0, CACHE, LOG, r0] = gap(W0, S, config, CACHE, LOG);
            % The draws of evolve_population's mutation step: a connection toggled with probability
            % mutationProbArchitecture, else a weight changed by N(0, mutationStdWeight^2); a
            % refused draw is drawn again, up to 20 times, then the mutation is dropped
            for attempt = 0:20
                if rand() < config.mutationProbArchitecture
                    g1 = randi(nGenes); g2 = randi(nGenes);
                    tog = true; z = 0;
                else
                    g1 = randi(nGenes); g2 = randi(nGenes);
                    tog = false; z = randn();
                end
                Mi = M(:, :, i); Bi = B(:, :, i);
                if tog, Mi(g1, g2) = ~Mi(g1, g2); else, Bi(g1, g2) = Bi(g1, g2) + z * config.mutationStdWeight; end
                Wm = double(Mi) .* double(Bi); dl = Wm(g1, g2) - W0(g1, g2);
                [gm, CACHE, LOG, rm] = gap(Wm, S, config, CACHE, LOG);
                % 'REFUSE_LESS_CONDITIONAL' refuses a lower Gap, 'REFUSE_LESS_RETALIATING' a higher
                % Ret; the tolerances 1e-9 are numerical, not model parameters
                refGap = gm < g0 - 1e-9; refRet = rm > r0 + 1e-9;
                if retOnly, refuse = refRet; else, refuse = refGap; end
                if ~isempty(split)    % both clone values are already stored by gap() above
                    [cP, CACHE, LOG] = coop(W0, W0, config, CACHE, LOG); [cM, CACHE, LOG] = coop(Wm, Wm, config, CACHE, LOG);
                    [cls, v, CACHE, BC, LOG] = mutation_class_end(W0, Wm, cP, cM, config, CACHE, BC, LOG);
                    refuse = any(cls == refSet);
                    if refuse, LOG.refCls(cls) = LOG.refCls(cls) + 1; end
                    if attempt == 0
                        LOG.endFirst = LOG.endFirst + (cls > 0); LOG.endResp = LOG.endResp + (cls == 1); LOG.endLoc = LOG.endLoc + (cls == 1 && isnan(v(3)));
                        LOG.endOrd = LOG.endOrd + (cls == 2); LOG.endUndef = LOG.endUndef + (cls == 3);
                    end
                end
                LOG.checks = LOG.checks + 1; LOG.refChecks = LOG.refChecks + refuse;
                if attempt == 0
                    LOG.mut = LOG.mut + 1; LOG.firstSq = LOG.firstSq + dl^2; LOG.refFirst = LOG.refFirst + refuse;
                    LOG.refFirstGap = LOG.refFirstGap + refGap; LOG.refFirstRet = LOG.refFirstRet + refRet;
                    LOG.refFirstGapOnly = LOG.refFirstGapOnly + (refGap && ~refRet);
                end
                if ~(refusing && refuse)
                    M(:, :, i) = Mi; B(:, :, i) = Bi;
                    LOG.accN = LOG.accN + 1; LOG.accSq = LOG.accSq + dl^2; LOG.accRef = LOG.accRef + refuse;
                    LOG = reservoir(LOG, RS, W0, Wm, S);
                    break;
                end
                LOG.turned = LOG.turned + 1;
                if attempt == 20, LOG.fail = LOG.fail + 1; break; end
                LOG.redraws = LOG.redraws + 1;
            end
        end
    end
    varargout = {M, B, K};
end


function L = new_log()
    L = struct('started', false, 'firstState', [], 'calls', 0, 'noPop', 0, 'mut', 0, 'refFirst', 0, ...
        'checks', 0, 'refChecks', 0, 'turned', 0, 'redraws', 0, 'fail', 0, 'accN', 0, 'accSq', 0, ...
        'accRef', 0, 'firstSq', 0, 'refFirstGap', 0, 'refFirstRet', 0, 'refFirstGapOnly', 0, 'nEval', 0, 'nSeen', 0, 'sampleP', {{}}, 'sampleM', {{}}, 'sampleS', {{}}, ...
        'endFirst', 0, 'endResp', 0, 'endLoc', 0, 'endOrd', 0, 'endUndef', 0, 'refCls', zeros(1, 3), ...
        'endFirstRef', 0, 'endAcc', 0, 'thinCoins', 0, 'thinSeed', NaN);
end


function [TS, LOG] = refuse_stream(TS, LOG, config)
% REFUSE_PARTIAL coin stream, created at the first call after 'reset' from config.refusalSeed.
    p = config.refusalShare;
    assert(isscalar(p) && p >= 0 && p <= 1, 'mutate_with_refusal: refusalShare must lie in [0, 1]');
    seed = config.refusalSeed;
    if isempty(TS)
        TS = RandStream('mt19937ar', 'Seed', seed); LOG.thinSeed = seed;
    else
        assert(LOG.thinSeed == seed, 'mutate_with_refusal: coin seed changed without reset');
    end
end


function [M, B, CACHE, LOG] = refuse_mutate(M, B, config, TS, CACHE, LOG, RS)
% REFUSE_PARTIAL: the main loop's draws and counts, with only the check that a mutation ends cooperation
% and the coin.
    [nGenes, ~, N] = size(M);
    for i = 1:N
        nMutations = poissrnd(config.mutationRate);
        for m = 1:nMutations
            W0 = double(M(:, :, i)) .* double(B(:, :, i));
            [cP, CACHE, LOG] = coop(W0, W0, config, CACHE, LOG);
            for attempt = 0:20
                if rand() < config.mutationProbArchitecture
                    g1 = randi(nGenes); g2 = randi(nGenes);
                    tog = true; z = 0;
                else
                    g1 = randi(nGenes); g2 = randi(nGenes);
                    tog = false; z = randn();
                end
                Mi = M(:, :, i); Bi = B(:, :, i);
                if tog, Mi(g1, g2) = ~Mi(g1, g2); else, Bi(g1, g2) = Bi(g1, g2) + z * config.mutationStdWeight; end
                Wm = double(Mi) .* double(Bi); dl = Wm(g1, g2) - W0(g1, g2);
                [refuse, ending, CACHE, LOG] = refuse_check(cP, Wm, config, TS, CACHE, LOG);
                LOG.checks = LOG.checks + 1; LOG.refChecks = LOG.refChecks + refuse;
                if attempt == 0
                    LOG.mut = LOG.mut + 1; LOG.firstSq = LOG.firstSq + dl^2; LOG.refFirst = LOG.refFirst + refuse;
                    LOG.endFirst = LOG.endFirst + ending; LOG.endFirstRef = LOG.endFirstRef + refuse;
                end
                if ~refuse
                    M(:, :, i) = Mi; B(:, :, i) = Bi;
                    LOG.accN = LOG.accN + 1; LOG.accSq = LOG.accSq + dl^2; LOG.endAcc = LOG.endAcc + ending;
                    LOG = reservoir(LOG, RS, W0, Wm, {});
                    break;
                end
                LOG.turned = LOG.turned + 1;
                if attempt == 20, LOG.fail = LOG.fail + 1; break; end
                LOG.redraws = LOG.redraws + 1;
            end
        end
    end
end


function [refuse, ending, CACHE, LOG] = refuse_check(cP, Wm, config, TS, CACHE, LOG)
% A draw ends cooperation when cP >= 0.5 and cM < 0.5; it is then refused when its coin u < p.
    [cM, CACHE, LOG] = coop(Wm, Wm, config, CACHE, LOG);
    ending = cP >= 0.5 && cM < 0.5; refuse = false;
    if ~ending, return; end
    u = rand(TS); LOG.thinCoins = LOG.thinCoins + 1;
    force = isfield(config, 'refusalNeverRefuse') && isequal(config.refusalNeverRefuse, true);
    refuse = u < config.refusalShare && ~force;
end


function L = reservoir(L, RS, W0, Wm, S)
% Reservoir sample of 50 accepted mutations (Algorithm R), on its own stream.
    L.nSeen = L.nSeen + 1;
    if numel(L.sampleP) < 50
        L.sampleP{end + 1} = W0; L.sampleM{end + 1} = Wm; L.sampleS{end + 1} = S;
    else
        j = randi(RS, L.nSeen);
        if j <= 50, L.sampleP{j} = W0; L.sampleM{j} = Wm; L.sampleS{j} = S; end
    end
end


function k = wkey(W)
    k = sprintf('%.17g,', W);
end


function assert_split_config(config)
% The class rules need C = gene 1 and D = gene 2 (response_off cuts W(1:2,3:4)) and no constant
% input in the run.
    assert(~isfield(config, 'geneInputBias') || isempty(config.geneInputBias), 'mutate_with_refusal: config carries geneInputBias');
    assert(config.geneCooperator == 1 && config.geneDefector == 2, 'mutate_with_refusal: split modes need C = gene 1, D = gene 2');
end


function [cls, v, CACHE, BC, LOG] = mutation_class_end(Wp, Wm, cP, cM, config, CACHE, BC, LOG)
% Class of a mutation: 0 if it does not end cooperation (coop(parent, parent) >= 0.5 and
% coop(mutant, mutant) < 0.5); 1, through the response to the partner; 2, ordinary; 3, undefined.
% A mutation that changes only entries of W(1:2,3:4) or W(3:4,:) is class 1. Otherwise both
% parent and mutant are switched off (response_off: W(1:2,3:4) = 0, with the constant input that
% holds the parent's resting state): class 3 if that input is undefined or the switched-off
% parent does not cooperate with its clone; class 1 if the switched-off mutant still does; class
% 2 otherwise. 'REFUSE_ALL' refuses classes 1-3, 'REFUSE_RESPONSE' class 1, 'REFUSE_OTHER' classes 2 and 3.
% v = [cP cM cPo cMo]. The constant
% input and the switched-off values are computed only for mutations that end cooperation; the
% input is stored per parent in BC.
    v = [cP cM NaN NaN];
    if cP < 0.5 || cM >= 0.5, cls = 0; return; end
    d = (Wm ~= Wp);    % a change only in W(1:2,3:4) or W(3:4,:) is class 1 by location
    if any(d(:)) && ~any(any(d(1:2, 1:2))), cls = 1; return; end
    key = wkey(Wp);
    if isKey(BC, key)
        row = BC(key);
    else
        sp = response_off('params');
        [b, ok] = response_off('bias', Wp, sp); row = [b(1) b(2) ok 0 0]; BC(key) = row;
    end
    if ~row(3), cls = 3; return; end
    cfgB = config; cfgB.geneInputBias = [row(1); row(2); 0; 0];
    Wpo = response_off('off', Wp); Wmo = response_off('off', Wm);
    [cPo, CACHE, LOG] = coop(Wpo, Wpo, cfgB, CACHE, LOG); [cMo, CACHE, LOG] = coop(Wmo, Wmo, cfgB, CACHE, LOG);
    v = [cP cM cPo cMo];
    if cPo < 0.5, cls = 3; elseif cMo >= 0.5, cls = 1; else, cls = 2; end
end


function Wd = default_defector(config)
% The zero wiring with D self-activation 1.0: a partner that defects by expressing D.
    Wd = zeros(config.numGenes); Wd(config.geneDefector, config.geneDefector) = 1.0;
end


function [c, CACHE, LOG] = coop(W, P, config, CACHE, LOG)
% coop(W, P): mean cooperation of W against P over NENC = 100 encounters under the runs' encounter
% rule (Methods, "Model"), played by evolve_population itself ('__test_evaluateFitness__') in a
% population of 100 copies of W paired with 100 copies of P (for P = W, copy i with copy i + 50),
% each pair developed for at most config.developmentSteps steps (300) from a random low start
% x_j ~ U(0, 0.01) drawn under rng(20261022), the same starts for every wiring. A network
% cooperates with its clone when coop(W, W) >= 0.5. Stored per wiring pair.
    key = [config.sensingMode num2str(config.hillCoefficient) '|' wkey(W) '|' wkey(P)];
    if isfield(config, 'geneInputBias') && ~isempty(config.geneInputBias)    % a constant input in config
        key = [key '|bias ' wkey(config.geneInputBias)];
    end
    if isKey(CACHE, key), c = CACHE(key); return; end
    NENC = 100; SEED = 20261022;
    s = rng; rng(SEED, 'twister');
    if isequal(W, P)
        n = NENC; Wp = repmat(W, 1, 1, n); opp = [n/2 + (1:n/2), 1:n/2]'; use = 1:n;
    else
        n = 2 * NENC; Wp = cat(3, repmat(W, 1, 1, NENC), repmat(P, 1, 1, NENC));
        opp = [NENC + (1:NENC), 1:NENC]'; use = 1:NENC;
    end
    cfg = config; cfg.precomputedOpponents = opp; cfg.populationSize = n;
    K = ones(config.numGenes, n, 'single') * config.michaelisConstant;
    R = evolve_population('__test_evaluateFitness__', true(config.numGenes, config.numGenes, n), single(Wp), K, cfg);
    rng(s);
    cv = R{3}; c = mean(cv(use));
    CACHE(key) = c; LOG.nEval = LOG.nEval + 1;
end


function [c, b, fa, pa] = coop_pair(W, P, config)
% coop(W, P) with the partner's side of the same encounters: the same population, pairing and
% starts as coop. evolve_population is called with an all-false external-partner mask, so no encounter
% is replaced and it also returns each pair's two actions (R{4}.focalCoop, R{4}.oppDefected);
% pair p has W's copy p as focal. Not stored.
    NENC = 100; SEED = 20261022;
    s = rng; rng(SEED, 'twister');
    if isequal(W, P)
        n = NENC; Wp = repmat(W, 1, 1, n); opp = [n/2 + (1:n/2), 1:n/2]'; use = 1:n;
    else
        n = 2 * NENC; Wp = cat(3, repmat(W, 1, 1, NENC), repmat(P, 1, 1, NENC));
        opp = [NENC + (1:NENC), 1:NENC]'; use = 1:NENC;
    end
    cfg = config; cfg.precomputedOpponents = opp; cfg.populationSize = n;
    K = ones(config.numGenes, n, 'single') * config.michaelisConstant;
    Kx = ones(config.numGenes, 1) * config.michaelisConstant;
    R = evolve_population('__test_evaluateFitness__', true(config.numGenes, config.numGenes, n), single(Wp), K, cfg, ...
        false(n, 1), zeros(config.numGenes), Kx);
    rng(s);
    cv = R{3}; c = mean(cv(use));
    fa = R{4}.focalCoop(use); pa = ~R{4}.oppDefected(use);
    assert(isequal(double(fa(:)), cv(use)), 'coop_pair: per-encounter focal actions differ from cooperationVec');
    b = mean(pa);
end


function [g, CACHE, LOG, cd] = gap(W, S, config, CACHE, LOG)
% Gap(W) = coop(W, W) - Ret(W); cd = Ret(W), the highest coop(W, P) over the partner set S, high
% when W retaliates little.
    [cc, CACHE, LOG] = coop(W, W, config, CACHE, LOG); cd = -inf;
    for q = 1:numel(S)
        [x, CACHE, LOG] = coop(W, S{q}, config, CACHE, LOG); cd = max(cd, x);
    end
    g = cc - cd;
end


function [S, CACHE, LOG] = partners(M, B, config, CACHE, LOG)
% Partner set: the default defector (zero wiring with D self-activation 1.0), and the commonest wiring of the
%              population that does not
% cooperate with its clone (if any).
    S = {default_defector(config)};
    N = size(M, 3); [Wu, ~, ic] = unique(reshape(double(M) .* double(B), 16, N)', 'rows');
    cnt = accumarray(ic, 1); [~, ord] = sort(cnt, 'descend');
    for u = ord'
        W = reshape(Wu(u, :), 4, 4); [cc, CACHE, LOG] = coop(W, W, config, CACHE, LOG);
        if cc < 0.5, S{2} = W; return; end
    end
end
