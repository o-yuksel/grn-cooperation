function out = fig1_partner_data(memberFile, rebuild)
% FIG1_PARTNER_DATA  Network 1's readout over development when network 2 is a copy of it
% (conditional) and when network 2 is an unconditional defector: the two curves of Fig. 1B and the
% two meetings marked in Fig. 1C.
%
% Implements : the data of Fig. 1B and 1C (read by fig1_encounter_pair)
% Inputs     : memberFile   (optional) a saved cache of fig1_member_and_mutant (struct out), read
%                           only; default: fig1_member_and_mutant(), which reads or builds its cache
%                           output/figures/fig1_member_and_mutant.mat
%              rebuild      (optional) true rebuilds; otherwise the cache is read when it exists
% Outputs    : output/figures/fig1_partner_data.mat   (the cache: struct out, also returned);
%              prints the partner chosen and the readouts at 300 steps
% Seed       : rng(12345 + 81, 'twister'), the runs' seed for replicate 81, before the partner
%              search (the stream is not restored); simulate_pair has no random draws
%
% Network 1 is the member that fig1_member_and_mutant picks (replicate 81 of the saturating-contrast
% runs, the median cooperation margin). No evolution here: pairs are replayed with simulate_pair,
% one call per horizon. The copy's trace must equal the cached r_clone exactly (error otherwise).
% The unconditional defector is chosen by a fixed rule: wirings are drawn one at a time
% as the runs draw their starting wirings (evolve_population: M = rand(4) < initConnectivity,
% B = randn(4) * initSigma, single), and the first that (i) defects clearly, readout x_C - x_D
% below -0.1, when its sensors receive a constant partner signal for full cooperation (s = +1,
% S_C +1 and S_D -1 under contrast sensing) and again for full defection (s = -1), each against a
% silent partner (M = B = 0) for 300 steps, and (ii) defects when paired with network 1 for 300
% steps is the partner. Network 1 must defect back against it (error otherwise, printing against
% how many of the first 50 such defectors network 1 defects back). The curve is evidence for that
% one pair only, and "unconditional" is tested at the two extreme inputs only.

    here = fileparts(mfilename('fullpath'));
    outDir = fullfile(fileparts(here), 'output', 'figures');
    cacheF = fullfile(outDir, 'fig1_partner_data.mat');
    if nargin < 1, memberFile = ''; end
    if nargin < 2 || isempty(rebuild), rebuild = false; end
    if ~rebuild && exist(cacheF, 'file')
        S = load(cacheF); out = S.out;
        fprintf('fig1_partner_data: cache hit (%s)\n', cacheF);
        return;
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    if isempty(memberFile)
        L = fig1_member_and_mutant();
    else
        S = load(memberFile, 'out'); L = S.out;
    end
    cfg = L.config;
    H = L.horizonEvolve;                 % 300
    sp = sim_params(cfg, H);
    MA = L.M; BA = L.B;

    % --- the conditional partner: a copy of network 1 --------------------------
    [rCond, rCondPartner] = readout_trace(MA, BA, MA, BA, sp, L.horizonTrace);
    dCopy = max(abs(rCond - L.r_clone));
    fprintf('Network 1 vs a copy: max |trace - cached r_clone| = %.3g over %d steps\n', dCopy, L.horizonTrace);
    if dCopy ~= 0
        error('fig1_partner_data:copy', 'the copy trace does not reproduce the cached r_clone (max diff %.3g)', dCopy);
    end

    % --- the unconditional partner ---------------------------------------------
    REP = 81; seed = 12345 + REP;
    rng(seed, 'twister');
    sigC = [0; 0; +1; -1]; sigD = [0; 0; -1; +1];
    MARGIN = 0.1;                        % test (i): readout below -MARGIN to both inputs
    nWant = 50; found = struct('draw', {}, 'M', {}, 'B', {}, 'rToC', {}, 'rToD', {}, ...
        'rPair', {}, 'rA', {});
    draw = 0;
    while numel(found) < nWant
        draw = draw + 1;
        M = rand(4) < cfg.initConnectivity;
        B = single(randn(4) * cfg.initSigma);
        rToC = solo_readout(M, B, sigC, sp);
        rToD = solo_readout(M, B, sigD, sp);
        if ~(rToC < -MARGIN && rToD < -MARGIN), continue; end
        [~, fs] = simulate_pair(MA, BA, [], M, B, [], sp);
        rA = fs(1) - fs(2); rP = fs(5) - fs(6);
        if ~(rP < 0), continue; end
        found(end+1) = struct('draw', draw, 'M', M, 'B', B, 'rToC', rToC, 'rToD', rToD, ...
            'rPair', rP, 'rA', rA); %#ok<AGROW>
        if numel(found) == 1
            fprintf('First unconditional defector: draw %d; readout %+.4f to s = +1, %+.4f to s = -1; with network 1 %+.4f; network 1 %+.4f\n', ...
                draw, rToC, rToD, rP, rA);
            if rA < 0, break; end
        end
    end
    P1 = found(1);
    if P1.rA >= 0
        nBack = sum([found.rA] < 0);
        fprintf('Network 1 defects back against %d of the first %d unconditional defectors (draws 1-%d).\n', ...
            nBack, numel(found), draw);
        error('fig1_partner_data:no_defect_back', ...
            'Network 1 does not defect back against the first unconditional defector (draw %d, readout %+.4f)', P1.draw, P1.rA);
    end

    [rUncond, rUncondPartner] = readout_trace(MA, BA, P1.M, P1.B, sp, L.horizonTrace);
    fprintf('At %d steps: network 1 vs conditional %+.4f, vs unconditional %+.4f; unconditional partner %+.4f\n', ...
        H, rCond(H), rUncond(H), rUncondPartner(H));
    if abs(rUncond(H) - P1.rA) > 0 || abs(rUncondPartner(H) - P1.rPair) > 0
        error('fig1_partner_data:trace', 'trace end differs from the selection run');
    end

    out = struct();
    out.rep             = L.rep;
    out.individual      = L.individual;
    out.M               = MA;
    out.B               = BA;
    out.t               = (1:L.horizonTrace)';
    out.r_cond          = rCond;            % network 1, partner a copy of itself
    out.r_cond_partner  = rCondPartner;
    out.r_uncond        = rUncond;          % network 1, partner the unconditional defector
    out.r_uncond_partner = rUncondPartner;
    out.partnerM        = P1.M;
    out.partnerB        = P1.B;
    out.partnerSeed     = seed;
    out.partnerDraw     = P1.draw;
    out.partnerSolo     = [P1.rToC P1.rToD];
    out.sensorInputs    = [sigC sigD];
    out.horizonTrace    = L.horizonTrace;
    out.horizonEvolve   = H;
    out.config          = cfg;
    save(cacheF, 'out');
    fprintf('Cached: %s\n', cacheF);
end

% =========================================================================
function sp = sim_params(cfg, H)
% The simulate_pair parameters of the run, H development steps.
    sp = struct('sensingMode', cfg.sensingMode, ...
        'hillCoefficient', cfg.hillCoefficient, ...
        'discreteOutput', cfg.discreteOutput, ...
        'michaelisConstant', cfg.michaelisConstant, ...
        'degradationRate', cfg.degradationRate, ...
        'timeStep', cfg.timeStep, ...
        'developmentSteps', H, ...
        'sensorIntegrationRate', cfg.sensorIntegrationRate, ...
        'externalSignal1', zeros(4,1), 'externalSignal2', zeros(4,1));
end

function r = solo_readout(M, B, sig, sp)
% Readout of a wiring whose sensors receive the constant input sig, as network 1 against a silent
% partner (which adds no signal).
    spk = sp; spk.externalSignal1 = sig;
    [~, fs] = simulate_pair(M, B, [], zeros(4), zeros(4), [], spk);
    r = fs(1) - fs(2);
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
