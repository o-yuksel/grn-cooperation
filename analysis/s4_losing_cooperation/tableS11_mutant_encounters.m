function tableS11_mutant_encounters(phase)
% TABLES11_MUTANT_ENCOUNTERS  How often does a member's defecting mutant earn the temptation T
% against the partners it meets in its population (its parent's wiring, other members,
% non-members), and does it earn less than its parent against the same partners?
%
% Produces : SI S4.4, "Members also defect back against other members' defecting mutants";
%            Table S11
% Inputs   : phase   'full' (the paper's numbers) or 'sample' (5 populations at t0, a timing run)
%            publication/mutation_off/branch.mat   (the 100 saturating-contrast populations at t0,
%              the generation their cooperation first exceeds 0.90; part of the data archive, a
%              copy of refusal_start_populations/branch.mat, which run_refusal_start_populations writes)
%            publication/<architecture>_100reps.mat   (run_signaling_architectures; the last
%              stored generation, 500)
% Outputs  : in publication/encounter_representativeness/: encounter_results.txt (the printed
%              report; log_sample.txt for 'sample') and, for 'full', encounter_rows.mat (one row
%              per encounter, and per-population summaries)
% Run      : tableS11_mutant_encounters('full')   (not timed; opens 3 workers if no pool is open
%            and closes the pool at the end); 'sample' takes about 1 min
% Seed     : population k under rng(20261109 + k, 'twister'), the same k at t0 and at
%            generation 500; membership_test and mutate_with_refusal('coop') use their own fixed
%            seeds and restore this stream
%
% A defecting mutant: one random mutation (the runs' mutation distribution) of a member parent
% drawn from the population, kept when the parent cooperates with its clone in at least half of
% 100 encounters under the runs' encounter rule and the mutant does not (mutate_with_refusal
% 'coop', the rule of SI S4.6). Each mutant and its parent meet the same 50 random partners.

    here = fullfile(grn_data_file('publication'), 'encounter_representativeness');
    pub = fullfile(grn_data_file('publication'));
    if strcmp(phase, 'sample'), logf = fullfile(here, 'log_sample.txt'); else, logf = fullfile(here, 'encounter_results.txt'); end
    if exist(logf, 'file'), delete(logf); end
    diary(logf); cl = onCleanup(@() diary('off'));
    fprintf('===== ENCOUNTERS OF DEFECTING MUTANTS %s  %s UTC =====\n', phase, char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss'));
    fprintf('Table S11: defecting mutants of member parents against partners drawn from their population\n');
    switch phase
        case 'sample', run_sample(pub);
        case 'full',   run_full(pub, here);
        otherwise, error('tableS11_mutant_encounters:phase', 'phase must be sample or full');
    end
    fprintf('\nDONE %s UTC\n', char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd HH:mm:ss'));
end


% =========================================================================
function run_sample(pub)
    [c, pops] = onset_pops(pub);
    pops = pops(1:5); t0 = tic;
    S = run_set(pops, c, false);
    fprintf('\nSAMPLE: 5 onset populations in %.1f s (serial)\n', toc(t0));
    for p = 1:numel(S)
        s = S(p);
        fprintf('  rep %3d: distinct %d, member share %.3f, draws %d (non-member parent draws %d, silent %d, parent clone<0.5 %d), defecting mutants %d, %.1f s\n', ...
            s.k, s.nDistinct, s.memShare, s.nDraw, s.nNonMemDraw, s.nSilent, s.nParentFail, s.nMut, s.secs);
        if s.nMut > 0
            o = s.rows(:, 7); fprintf('     mutant outcome mix CC %d CD %d DC %d DD %d; mean payoff diff %.3f\n', ...
                sum(o == 1), sum(o == 2), sum(o == 3), sum(o == 4), mean(s.dpay));
        end
    end
    report('onset (sample, 5)', S);
end


function run_full(pub, here)
% The sets of Table S11: saturating contrast at t0 and at generation 500 (successful replicates,
% final cooperation >= 0.95); the other three architectures at generation 500, successful and
% failed replicates separately (the table shows successful).
    T0 = tic;
    [c, pops] = onset_pops(pub);
    p = gcp('nocreate'); if isempty(p), parpool('Processes', 3); end
    sets = {}; tagRows = {};
    t = tic; S = run_set(pops, c, true); report(sprintf('SATURATING CONTRAST, onset (t0), %.1f min', toc(t) / 60), S);
    sets{end + 1} = {'satcon_onset', S};
    arch = {'saturating contrast', 'saturating_contrast_100reps.mat', true; ...
            'graded separate', 'graded_separate_100reps.mat', false; ...
            'graded contrast', 'graded_contrast_100reps.mat', false; ...
            'saturating separate', 'saturating_separate_100reps.mat', false};
    for a = 1:size(arch, 1)
        [ca, ps, succ] = gen500_pops(pub, arch{a, 2});
        for wantSucc = [true false]
            if arch{a, 3} && ~wantSucc, continue; end
            lab = 'failed'; if wantSucc, lab = 'successful'; end
            t = tic; S = run_set(ps(succ == wantSucc), ca, true);
            report(sprintf('%s, generation 500, %s replicates (%d), %.1f min', upper(arch{a, 1}), lab, sum(succ == wantSucc), toc(t) / 60), S);
            sets{end + 1} = {sprintf('%s_gen500_%s', strrep(arch{a, 2}, '_100reps.mat', ''), lab), S}; %#ok<AGROW>
        end
    end
    save_rows(here, sets);
    fprintf('\nTOTAL %.1f min\n', toc(T0) / 60);
    delete(gcp('nocreate'));
end


% =========================================================================
function [c, pops] = onset_pops(pub)
    d = load(fullfile(pub, 'saturating_contrast_100reps.mat'), 'config'); c = d.config; assert_config(c);
    L = load(fullfile(pub, 'mutation_off', 'branch.mat')); br = L.br;
    pops = struct('k', num2cell([br.k]), 'M', {br.M}, 'B', {br.B}, 'gen', num2cell([br.t0]));
end


function [c, pops, succ] = gen500_pops(pub, fname)
    d = load(fullfile(pub, fname)); c = d.config; assert_config(c);
    RR = d.results;
    n = numel(RR); succ = false(1, n); pops = struct('k', cell(1, n), 'M', [], 'B', [], 'gen', []);
    for r = 1:n
        s = size(RR{r}.historyM, 4); succ(r) = RR{r}.historyCooperation(end) >= 0.95;
        pops(r).k = r; pops(r).M = RR{r}.historyM(:, :, :, s); pops(r).B = RR{r}.historyB(:, :, :, s);
        pops(r).gen = RR{r}.historyGenerations(s);
    end
    if isfield(d, 'successRate'), fprintf('%s: successful %d, file successRate %g\n', fname, sum(succ), d.successRate); end
end


function assert_config(c)
% Stops if a config sets an evolve_population option that the reproduced mutation draw below omits.
    for fn = {'fixedMutualInhibition', 'frozenGenes', 'frozenConnections', 'evolvableKM_genes', 'geneInputBias'}
        if isfield(c, fn{1}) && ~isempty(c.(fn{1})) && ~isequal(c.(fn{1}), false)
            error('tableS11_mutant_encounters:config', 'config.%s is set; the mutation distribution reproduced here does not include it', fn{1});
        end
    end
    assert(c.geneCooperator == 1 && c.geneDefector == 2 && c.numGenes == 4, 'gene order');
end


function sp = pair_params(c)
% Encounters as in the runs: the file config's development settings (developmentSteps 300).
    sp = struct('sensingMode', c.sensingMode, 'hillCoefficient', c.hillCoefficient, 'discreteOutput', c.discreteOutput, ...
        'michaelisConstant', c.michaelisConstant, 'degradationRate', c.degradationRate, 'timeStep', c.timeStep, ...
        'developmentSteps', c.developmentSteps, 'payoffDC', c.payoffDC, 'sensorIntegrationRate', c.sensorIntegrationRate);
end


function cp = css_params(c)
% The membership test: 600 development steps from its fixed start (SI S2.3).
    cp = struct('sensingMode', c.sensingMode, 'hillCoefficient', c.hillCoefficient, 'discreteOutput', c.discreteOutput, ...
        'michaelisConstant', c.michaelisConstant, 'degradationRate', c.degradationRate, 'timeStep', c.timeStep, ...
        'developmentSteps', 600, 'backend', 'simulate_pair');
end


% =========================================================================
function S = run_set(pops, c, par)
    n = numel(pops); C = cell(1, n);
    if par
        parfor p = 1:n, C{p} = measure_pop(pops(p), c); end
    else
        for p = 1:n, C{p} = measure_pop(pops(p), c); end
    end
    S = [C{:}];
end


function s = measure_pop(pop, c)
% One population: membership of its distinct wirings (one row of M.*B per individual, unique
% rows; ic maps each individual to its wiring), then mutations drawn on parents drawn uniformly
% from the population until 30 defecting mutants are found or 3,000 draws are made on member
% parents; each defecting mutant and its parent play the same 50 partners.
    t = tic; rng(20261109 + pop.k, 'twister'); mutate_with_refusal('reset');
    M = pop.M; B = pop.B; N = size(M, 3);
    [Wu, ~, ic] = unique(reshape(double(M) .* double(B), 16, N)', 'rows'); nU = size(Wu, 1);
    cp = css_params(c); mem = false(nU, 1);
    for u = 1:nU, mem(u) = membership_test(reshape(Wu(u, :), 4, 4), cp).is_css; end
    s = struct('k', pop.k, 'gen', pop.gen, 'nDistinct', nU, 'memShare', mean(mem(ic)), 'nDraw', 0, 'nNonMemDraw', 0, ...
        'nSilent', 0, 'nParentFail', 0, 'nMut', 0, 'rows', zeros(0, 11), 'dpay', zeros(0, 1), 'Wm', zeros(0, 16), 'Wp', zeros(0, 16), 'secs', 0, 'skipped', ~any(mem));
    if s.skipped, s.secs = toc(t); return; end
    sp = pair_params(c); mid = 0;
    while s.nMut < 30 && s.nDraw < 3000
        i = randi(N);
        if ~mem(ic(i)), s.nNonMemDraw = s.nNonMemDraw + 1; continue; end
        [Mm, Bm, tog] = mutation_draw(M(:, :, i), B(:, :, i), c); s.nDraw = s.nDraw + 1;
        Wp = double(M(:, :, i)) .* double(B(:, :, i)); Wm = double(Mm) .* double(Bm);
        if isequal(Wm, Wp), s.nSilent = s.nSilent + 1; continue; end
        cP = mutate_with_refusal('coop', Wp, Wp, c);
        if cP < 0.5, s.nParentFail = s.nParentFail + 1; continue; end
        if mutate_with_refusal('coop', Wm, Wm, c) >= 0.5, continue; end
        mid = mid + 1; s.nMut = mid;
        [rows, dp] = play_partners(Wm, Wp, ic(i), randperm(N, 50), ic, mem, Wu, sp, c);
        s.rows = [s.rows; [repmat([pop.k pop.gen mid double(tog)], 50, 1) rows]];
        s.dpay(end + 1, 1) = dp; s.Wm(end + 1, :) = Wm(:)'; s.Wp(end + 1, :) = Wp(:)';
    end
    s.secs = toc(t);
end


function [Mi, Bi, tog] = mutation_draw(Mi, Bi, c)
% One mutation as evolve_population mutatePopulation draws it, in the same draw order: a
% toggle of one connection of M with probability mutationProbArchitecture, else a weight change
% of N(0, mutationStdWeight^2) on B. A weight change where M is 0 leaves W = M.*B unchanged.
    nG = size(Mi, 1);
    if rand() < c.mutationProbArchitecture
        g1 = randi(nG); g2 = randi(nG); Mi(g1, g2) = ~Mi(g1, g2); tog = true;
    else
        g1 = randi(nG); g2 = randi(nG); Bi(g1, g2) = Bi(g1, g2) + randn() * c.mutationStdWeight; tog = false;
    end
end


function [rows, dp] = play_partners(Wm, Wp, uParent, partners, ic, mem, Wu, sp, c)
% rows: [partnerWiringId partnerClass mutantOutcome parentOutcome mutantPay parentPay parentWiringId]
% outcome codes 1 CC, 2 CD, 3 DC, 4 DD (focal first; focal = mutant or parent); 3 is the mutant
% earning T. partnerClass 1 the parent's wiring, 2 another member wiring, 3 a non-member. The
% mutant and the parent meet each partner from the same two random starts U(0, 0.01).
    Msq = ones(4); n = numel(partners); rows = zeros(n, 7);
    for q = 1:n
        u = ic(partners(q));
        if u == uParent, cls = 1; elseif mem(u), cls = 2; else, cls = 3; end
        Wq = reshape(Wu(u, :), 4, 4);
        sp.initialState1 = rand(4, 1) * 0.01; sp.initialState2 = rand(4, 1) * 0.01;
        om = code(simulate_pair(Msq, Wm, [], Msq, Wq, [], sp));
        op = code(simulate_pair(Msq, Wp, [], Msq, Wq, [], sp));
        rows(q, :) = [u cls om op pay(om, c) pay(op, c) uParent];
    end
    dp = mean(rows(:, 5)) - mean(rows(:, 6));
end


function k = code(o)
    k = find(strcmp(o, {'CC', 'CD', 'DC', 'DD'}));
end


function v = pay(k, c)
    P = [c.payoffCC c.payoffCD c.payoffDC c.payoffDD]; v = P(k);
end


function report(label, S)
% One set: populations used (holding a member wiring and yielding a defecting mutant), partner
% class shares, exploitation (mutant D, partner C) pooled over encounters and by class, and per
% population the median over its mutants of mutant minus parent mean payoff; across populations
% the median, the count negative and a two-sided exact sign test (zeros dropped).
    fprintf('\n===== %s =====\n', label);
    used = S(~[S.skipped] & [S.nMut] > 0);
    fprintf('populations %d; skipped (no member wiring) %d; members but no defecting mutant in 3000 draws %d; used %d\n', ...
        numel(S), sum([S.skipped]), sum(~[S.skipped] & [S.nMut] == 0), numel(used));
    if isempty(used), return; end
    nD = sum([S.nDraw]); nM = sum([S.nMut]); R = vertcat(used.rows);
    fprintf('member share of individuals: median %.3f [%.3f, %.3f] over populations\n', median([S.memShare]), min([S.memShare]), max([S.memShare]));
    fprintf('mutation draws on member parents %d: silent %d, parent coop with clone < 0.5 %d, defecting mutants %d (%.4f per draw); non-member parent draws %d\n', ...
        nD, sum([S.nSilent]), sum([S.nParentFail]), nM, nM / max(nD, 1), sum([S.nNonMemDraw]));
    fprintf('defecting mutants per used population: median %g [%d, %d]; mutants by topology toggle %d of %d\n', ...
        median([used.nMut]), min([used.nMut]), max([used.nMut]), round(sum(R(:, 4)) / 50), nM);
    cls = R(:, 6); om = R(:, 7); op = R(:, 8);
    names = {'same wiring as parent', 'different member', 'non-member'};
    fprintf('encounters %d. Partner-class shares: same wiring as parent %.3f; different member %.3f; non-member %.3f\n', ...
        size(R, 1), mean(cls == 1), mean(cls == 2), mean(cls == 3));
    fprintf('%-24s %8s %9s %9s %9s %9s | %9s %9s | %8s %6s\n', 'partner class', 'n', 'mut CC', 'mut CD', 'mut DC(T)', 'mut DD', 'par DC', 'par CD', 'popMedT', 'nPops');
    print_class('all', 0, om, op, cls, used);
    for k = 1:3, print_class(names{k}, k, om, op, cls, used); end
    fprintf('  (mut DC(T) = exploitation share pooled over encounters; popMedT = median over populations of the per-population share)\n');
    md = arrayfun(@(s) median(s.dpay), used); nn = sum(md < 0); np = sum(md > 0); nz = sum(md == 0);
    pv = min(1, 2 * binocdf(min(nn, np), nn + np, 0.5));
    fprintf('payoff difference (mutant minus parent, mean over the same 50 partners), per-population median over mutants:\n');
    fprintf('  across %d populations: median %.4f, IQR [%.4f, %.4f]; negative %d (%.3f), zero %d, positive %d; sign test p = %.3g\n', ...
        numel(md), median(md), prctile(md, 25), prctile(md, 75), nn, nn / numel(md), nz, np, pv);
    fprintf('  pooled mean payoff: mutant %.3f, parent %.3f; per-mutant differences negative in %.3f of %d mutants\n', ...
        mean(R(:, 9)), mean(R(:, 10)), mean(vertcat(used.dpay) < 0), nM);
    for k = 1:3
        sel = cls == k; fprintf('  partner class %-22s mean payoff mutant %.3f, parent %.3f\n', names{k}, mean(R(sel, 9)), mean(R(sel, 10)));
    end
    fprintf('time: %.1f s summed over populations\n', sum([S.secs]));
end


function print_class(name, k, om, op, cls, used)
    if k == 0, sel = true(size(cls)); else, sel = cls == k; end
    sh = @(v, j) mean(v(sel) == j);
    pm = nan(1, numel(used));
    for p = 1:numel(used)
        r = used(p).rows; if k == 0, s2 = true(size(r, 1), 1); else, s2 = r(:, 6) == k; end
        if any(s2), pm(p) = mean(r(s2, 7) == 3); end
    end
    fprintf('%-24s %8d %9.4f %9.4f %9.4f %9.4f | %9.4f %9.4f | %8.4f %6d\n', name, sum(sel), sh(om, 1), sh(om, 2), sh(om, 3), sh(om, 4), ...
        sh(op, 3), sh(op, 2), median(pm, 'omitnan'), sum(~isnan(pm)));
end


function save_rows(here, sets)
% Per-encounter rows of every set, compact types, plus per-population summaries (without the rows).
    names = cellfun(@(x) x{1}, sets, 'UniformOutput', false); R = cell(numel(sets), 1); P = cell(numel(sets), 1);
    for j = 1:numel(sets)
        S = sets{j}{2}; r = vertcat(S.rows); R{j} = [repmat(j, size(r, 1), 1) r];
        P{j} = rmfield(S, {'rows'}); [P{j}.set] = deal(j);
    end
    R = vertcat(R{:});
    rows = struct('set', uint8(R(:, 1)), 'replicate', uint16(R(:, 2)), 'generation', uint16(R(:, 3)), 'mutant', uint8(R(:, 4)), ...
        'toggle', logical(R(:, 5)), 'partnerWiring', uint16(R(:, 6)), 'partnerClass', uint8(R(:, 7)), ...
        'mutantOutcome', uint8(R(:, 8)), 'parentOutcome', uint8(R(:, 9)), 'mutantPayoff', uint8(R(:, 10)), ...
        'parentPayoff', uint8(R(:, 11)), 'parentWiring', uint16(R(:, 12))); %#ok<NASGU>
    legend = ['outcome 1 CC 2 CD 3 DC 4 DD, focal (mutant or parent) first; partnerClass 1 same wiring as parent, 2 different member, ' ...
        '3 non-member; wiring ids index the distinct wirings of that population (unique rows of M.*B, sorted); parent is always a member; ' ...
        'pops: per-population summary, Wm / Wp = defecting mutant and parent wirings (16 = W(:) per mutant)']; %#ok<NASGU>
    pops = [P{:}]; %#ok<NASGU>
    save(fullfile(here, 'encounter_rows.mat'), 'rows', 'names', 'pops', 'legend', '-v7');
    fprintf('\nsaved encounter_rows.mat: %d encounter rows, %d sets\n', size(R, 1), numel(sets));
end
