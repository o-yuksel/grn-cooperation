function out = fig1_member_and_mutant(rebuild)
% FIG1_MEMBER_AND_MUTANT  Picks the member drawn in Fig. 1 (network 1) and its worst-case
% single-weight mutant, and traces the member's readout over development against its clone and
% against that mutant.
%
% Produces : the member of Fig. 1 (read by fig1_partner_data, for fig1_encounter_pair) and the
%            member and mutant of the worked example in the main text, "One summed input ties
%            responsiveness to mutation-sensitivity", and SI S3.3 (read by worked_example)
% Inputs   : rebuild   (optional) true recomputes; otherwise the cache is read when it exists
%            publication/saturating_contrast_100reps.mat   (written by run_signaling_architectures)
% Outputs  : output/figures/fig1_member_and_mutant.mat   (the cache: struct out, also returned);
%            prints the pick, the mutant and the readouts at 300 and 600 steps
% Run      : fig1_member_and_mutant(true)
% Seed     : none here (membership_test seeds each wiring itself; simulate_pair has no random draws)
%
% The member: in replicate 81 at generation 500, of the distinct wirings that are members (the
% membership test at 600 development steps, every mutant tried, `gradedC2`) and have a defecting
% single-weight mutant, the one with the median cooperation margin (readout x_C - x_D against its
% clone; the lower median for an even count). The mutant: at the weight where membership_test
% records the worst-case mutant, of the 10 nonzero changes of the test's grid (-0.6 to +0.6,
% i.e. -3 to +3 times the standard deviation of a weight mutation) whose mutant defects on the
% member, the one against which the member's readout is highest (printed beside membership_test's
% value). Each trace point is one simulate_pair call run for k steps; the runs develop for 300.

    if nargin < 1 || isempty(rebuild), rebuild = false; end

    here    = fileparts(mfilename('fullpath'));
    outDir  = fullfile(fileparts(here), 'output', 'figures');
    cacheF  = fullfile(outDir, 'fig1_member_and_mutant.mat');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    if ~rebuild && exist(cacheF, 'file')
        S = load(cacheF);
        if isfield(S.out, 'horizonTestMembership')
            out = S.out;
            fprintf('fig1_member_and_mutant: cache hit (%s)\n', cacheF);
            return;
        end
        fprintf('fig1_member_and_mutant: cache %s has the old field names; rebuilding\n', cacheF);
    end

    REP    = 81;
    HORIZ  = 600;          % trace length; Fig. 1B draws steps 1 to 300
    EVOLVE = 300;          % the horizon the population evolved under

    dataF = grn_data_file(fullfile('publication', ...
        'saturating_contrast_100reps.mat'));
    fprintf('Loading %s ...\n', dataF);
    d = load(dataF, 'results', 'config', 'coopValues');
    cfgRun = d.config;
    rep    = d.results{REP};

    classP = membership_test_params(cfgRun);

    % --- distinct wirings in the gen-500 population -------------------------
    N = size(rep.finalNetworkM, 3);
    Wall = zeros(16, N);
    for k = 1:N
        Wk = rep.finalNetworkM(:,:,k) .* rep.finalNetworkB(:,:,k);
        Wall(:,k) = Wk(:);
    end
    [~, firstIdx] = unique(round(Wall', 12), 'rows', 'stable');
    fprintf('Replicate %d: %d individuals, %d distinct wirings.\n', ...
        REP, N, numel(firstIdx));

    cand = struct('idx', {}, 'W', {}, 'M', {}, 'B', {}, 'coop_margin', {}, ...
        'retal_margin', {}, 'testResult', {});
    for k = firstIdx(:)'
        M = rep.finalNetworkM(:,:,k);
        B = rep.finalNetworkB(:,:,k);
        res = membership_test(M .* B, classP);
        if res.is_css && isfinite(res.worst_retal_slot)
            cand(end+1) = struct('idx', k, 'W', M .* B, 'M', M, 'B', B, ...
                'coop_margin', res.x_star(1) - res.x_star(2), ...
                'retal_margin', -res.worst_retal, 'testResult', res);  %#ok<AGROW>
        end
    end
    if isempty(cand)
        error('fig1_member_and_mutant:no_member', ...
            'No member of the conditional genotype network with a defecting worst-case mutant among the distinct wirings of replicate %d.', REP);
    end
    [~, ord] = sort([cand.coop_margin]);
    pick = cand(ord(max(1, ceil(numel(ord)/2))));
    fprintf('%d distinct wirings are members with a defecting worst-case mutant; drawing the one with the median cooperation margin.\n', ...
        numel(cand));
    fprintf('  picked individual %d: cooperation margin %.4f, retaliation margin %.4f\n', ...
        pick.idx, pick.coop_margin, pick.retal_margin);

    % --- the worst-case single-weight mutant at membership_test's recorded slot --
    [W_worst_mutant, worst_mutant] = locate_worst_mutant(pick, classP);
    fprintf('  worst-case mutant: slot (%d,%d) [%s], weight change %+.4f; retaliation readout %.4f (membership test %.4f)\n', ...
        worst_mutant.ii, worst_mutant.jj, gene_pair_str(worst_mutant.ii, worst_mutant.jj), worst_mutant.delta, ...
        worst_mutant.reconstructed, pick.testResult.worst_retal);

    % --- readout traces, one simulate_pair call per horizon -----------------
    sp = sim_params(classP);
    Mones = ones(4);
    [tr_clone_self, tr_clone_partner] = readout_trace(pick.M, pick.B, pick.M, pick.B, sp, HORIZ);
    [tr_esc_self,   tr_esc_partner]   = readout_trace(pick.M, pick.B, Mones,  W_worst_mutant, sp, HORIZ);

    out = struct();
    out.rep             = REP;
    out.individual      = pick.idx;
    out.architecture    = 'saturating contrast';
    out.W               = pick.W;
    out.M               = pick.M;
    out.B               = pick.B;
    out.W_worst_mutant  = W_worst_mutant;
    out.worst_mutant    = worst_mutant;
    out.coop_margin     = pick.coop_margin;      % cooperation margin: readout against the clone (600 steps)
    out.retal_margin    = pick.retal_margin;     % retaliation margin, -worst_retal (600 steps)
    out.testResult      = pick.testResult;
    out.t               = (1:HORIZ)';
    out.r_clone         = tr_clone_self;         % member's readout, meeting its clone
    out.r_clone_partner = tr_clone_partner;
    out.r_vs_mut        = tr_esc_self;           % member's readout, meeting the mutant
    out.r_mut           = tr_esc_partner;        % the mutant's own readout
    out.horizonTrace    = HORIZ;
    out.horizonEvolve   = EVOLVE;
    out.horizonTestMembership = classP.developmentSteps;
    out.nCandidates     = numel(cand);
    out.nDistinct       = numel(firstIdx);
    out.coopValue       = d.coopValues(REP);
    out.config          = cfgRun;

    fprintf('\nReadout at the two horizons (member''s own readout):\n');
    for h = [EVOLVE HORIZ]
        fprintf('  %3d steps: vs clone %+.4f (%s), vs mutant %+.4f (%s); mutant''s readout %+.4f (%s)\n', ...
            h, out.r_clone(h), act(out.r_clone(h)), out.r_vs_mut(h), act(out.r_vs_mut(h)), ...
            out.r_mut(h), act(out.r_mut(h)));
    end

    save(cacheF, 'out');
    fprintf('Cached: %s\n', cacheF);
end

% =========================================================================
function p = membership_test_params(cfgRun)
% The parameters of the membership test (membership_test): 600 development steps (the runs develop
% for 300) and gradedC2, so worst_retal and worst_retal_slot are recorded; the mutation grid is
% 11 values from -3 to +3 sigma_mut (sigma_mut, the standard deviation of a weight mutation).
    p = struct( ...
        'sensingMode',           cfgRun.sensingMode, ...
        'hillCoefficient',       cfgRun.hillCoefficient, ...
        'discreteOutput',        cfgRun.discreteOutput, ...
        'michaelisConstant',     cfgRun.michaelisConstant, ...
        'degradationRate',       cfgRun.degradationRate, ...
        'timeStep',              cfgRun.timeStep, ...
        'developmentSteps',      600, ...
        'sigma_mut',             cfgRun.mutationStdWeight, ...
        'sensorIntegrationRate', cfgRun.sensorIntegrationRate, ...
        'backend',               'simulate_pair', ...
        'gradedC2',              true, ...
        'payoffs', [cfgRun.payoffDC, cfgRun.payoffCC, cfgRun.payoffDD, cfgRun.payoffCD]);
    p.cssMagGrid = linspace(-3*p.sigma_mut, 3*p.sigma_mut, 11);
end

function sp = sim_params(classP)
    sp = struct('sensingMode', classP.sensingMode, ...
        'hillCoefficient', classP.hillCoefficient, ...
        'discreteOutput', classP.discreteOutput, ...
        'michaelisConstant', classP.michaelisConstant, ...
        'degradationRate', classP.degradationRate, ...
        'timeStep', classP.timeStep, ...
        'developmentSteps', classP.developmentSteps, ...
        'sensorIntegrationRate', classP.sensorIntegrationRate, ...
        'externalSignal1', zeros(4,1), 'externalSignal2', zeros(4,1));
end

function [W_worst_mutant, info] = locate_worst_mutant(pick, classP)
% At membership_test's recorded slot, try every nonzero weight change of the grid; among the mutants
% that defect against the member, keep the one against which the member's readout is highest.
% The stored field was_connected says whether the slot is wired in the member's M.
    slot = pick.testResult.worst_retal_slot;
    jj = ceil(slot / 4); ii = slot - 4*(jj-1);
    sp = sim_params(classP);
    Mones = ones(4);
    bestR = -inf; bestDelta = NaN; W_worst_mutant = pick.W;
    for delta = classP.cssMagGrid
        if delta == 0, continue; end
        Wm = pick.W; Wm(ii,jj) = Wm(ii,jj) + delta;
        [outc, fin, ~] = simulate_pair(Mones, pick.W, [], Mones, Wm, [], sp);
        if outc(2) == 'D'
            rR = fin(1) - fin(2);
            if rR > bestR, bestR = rR; bestDelta = delta; W_worst_mutant = Wm; end
        end
    end
    info = struct('ii', ii, 'jj', jj, 'slot', slot, 'delta', bestDelta, ...
        'reconstructed', bestR, 'membership_test_value', pick.testResult.worst_retal, ...
        'was_connected', pick.M(ii,jj) ~= 0);
end

function [r1, r2] = readout_trace(M1, B1, M2, B2, sp, H)
% Readout x_C - x_D of both networks after k development steps, k = 1..H, each from one
% simulate_pair call at that horizon.
    r1 = zeros(H,1); r2 = zeros(H,1);
    spk = sp;
    for k = 1:H
        spk.developmentSteps = k;
        [~, fs, ~] = simulate_pair(M1, B1, [], M2, B2, [], spk);
        r1(k) = fs(1) - fs(2);
        r2(k) = fs(5) - fs(6);
    end
end

function s = gene_pair_str(ii, jj)
    g = {'C', 'D', 'S_C', 'S_D'};
    s = sprintf('%s <- %s', g{ii}, g{jj});
end

function c = act(r)
    if r > 0, c = 'C'; else, c = 'D'; end
end
