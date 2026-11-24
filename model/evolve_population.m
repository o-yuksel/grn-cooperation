function results = evolve_population(varargin)
% EVOLVE_POPULATION  Evolves populations of four-gene networks that play the Prisoner's Dilemma.
%
% Implements : Methods, "Model" and "Game and evolution": each generation every individual meets
%              one partner drawn at random, the two networks develop as one coupled system
%              (Eq. 3 and Eq. 4) and each one's action is read from the settled expression of
%              C against D (the runs' encounter rule); fitness is the payoff of the game played
%              as focal; parents are drawn by the selection weight of "Game and evolution" and offspring mutate.
%              Also the ring
%              (local replacement) of SI Extended Materials and Methods, "Spatial
%              reproduction", and group selection by fission and its random-split control
%              ('random_split') of "Group selection".
% Inputs     : config   struct, as built by publication_config (fields set and explained in
%                       base_config, setting_config and publication_config); optional fields
%                       are described where read; model variants that no paper result uses stop
%                       the run (rejectUnusedOptions)
%              or '__test_evaluateFitness__', M, B, K, config [, extMask, W_ext, K_ext]:
%                       one round of encounters for a given population, returning
%                       {fitness, cooperation rate, cooperation per individual [, encounters]};
%                       the analyses use it to read cooperation under the runs' encounter rule
% Outputs    : results  cell array, one struct per replicate: finalNetworkM, finalNetworkB,
%                       finalNetworkK (the last generation), historyCooperation (population
%                       cooperation per generation), historyParents (per generation, the parent
%                       of each slot, a position in the population), historyM, historyB and
%                       historyGenerations (the stored generations), snapshotCoopVectors and
%                       snapshotCoopGenerations
% Seed       : replicate r runs under rng(config.rngSeedBase + r, 'twister'), base 12345 when
%              config.rngSeedBase is absent; replicates run in a parfor loop
%
% The encounter rule in full (checks every 40 steps, at most 300 steps) is above
% simulateInteractionBatch; starting wirings are explained where read.

    % '__test_evaluateFitness__': one round of encounters for a given population (header, Inputs)
    if nargin >= 1 && ischar(varargin{1}) && strcmp(varargin{1}, '__test_evaluateFitness__')
        M_test = varargin{2};
        B_test = varargin{3};
        K_test = varargin{4};
        cfg_test = stored_spellings(varargin{5});
        rejectUnusedOptions(cfg_test);
        if length(varargin) >= 8
            % with an external partner W_ext (K_ext) for the individuals where extMask is true
            [f, cr, cv, enc] = evaluateFitness(M_test, B_test, K_test, cfg_test, ...
                varargin{6}, varargin{7}, varargin{8});
            results = {f, cr, cv, enc};
        else
            [f, cr, cv] = evaluateFitness(M_test, B_test, K_test, cfg_test);
            results = {f, cr, cv};
        end
        return;
    end

    config = stored_spellings(varargin{1});
    rejectUnusedOptions(config);
    validateConfig(config);

    % Combinations of population structures that are not defined stop the run
    if isfield(config, 'groupSize') && config.groupSize > 0 && ...
       isfield(config, 'spatialStructure') && strcmp(config.spatialStructure, 'ring')
        error('evolve_population:IncompatibleConfig', ...
            'groupSize > 0 is incompatible with spatialStructure=ring; run them separately');
    end

    % Groups are run only with fission (SI Extended Materials and Methods, "Group selection"),
    % which needs play within groups, an even group size and at least two groups
    if isfield(config, 'groupSize') && config.groupSize > 0 && ...
       ~(isfield(config, 'groupSelectionMode') && ...
         any(strcmp(config.groupSelectionMode, {'fission', 'random_split'})))
        error('evolve_population:IncompatibleConfig', ...
            'groupSize > 0 requires groupSelectionMode = fission or random_split');
    end
    if isfield(config, 'groupSelectionMode') && ...
       any(strcmp(config.groupSelectionMode, {'fission', 'random_split'})) && ...
       ~(isfield(config, 'groupSize') && config.groupSize > 0 && mod(config.groupSize, 2) == 0 && ...
         config.populationSize / config.groupSize >= 2 && ...
         isfield(config, 'groupPairingMode') && strcmp(config.groupPairingMode, 'within'))
        error('evolve_population:IncompatibleConfig', ...
            'groupSelectionMode=fission/random_split requires an even groupSize > 0, at least 2 groups and groupPairingMode=within');
    end

    nReps = config.numReplicates;
    results = cell(nReps, 1);

    sampleGens = unique([1:config.historySampleInterval:config.numGenerations, ...
                         config.alwaysSaveGenerations(config.alwaysSaveGenerations <= config.numGenerations)]);
    nSamples = length(sampleGens);

    fprintf('\n=== GRN EVOLUTION ===\n');
    fprintf('Genes: %d | Sensing: %s | Payoff: %s | Stability: %s | Selection: %s\n', ...
        config.numGenes, config.sensingMode, config.payoffMode, ...
        mat2str(config.requireStability), config.selectionAlgorithm);
    if isfield(config, 'spatialStructure')
        fprintf('Spatial: %s', config.spatialStructure);
        if strcmp(config.spatialStructure, 'ring')
            fprintf(' (k=%d)', config.neighborhoodSize);
        end
        fprintf('\n');
    end
    if config.developmentSteps <= 1
        fprintf('SETTING: single development step (instantaneous imposed)\n');
    end
    if isfield(config, 'initPopulationM')
        nSeed = size(config.initPopulationM, 3);
        if nSeed < config.populationSize
            fprintf('SEEDED START: %d seed individuals TILED to N = %d\n', nSeed, config.populationSize);
        else
            fprintf('SEEDED START: %d seed individuals (N = %d, no tiling)\n', nSeed, config.populationSize);
        end
    end
    fprintf('DevSteps: %d | Reps: %d | Gens: %d | Pop: %d\n', ...
        config.developmentSteps, nReps, config.numGenerations, config.populationSize);
    fprintf('VECTORIZED: pagemtimes batch ODE');
    fprintf('\n============================\n\n');

    tic;
    parfor repIdx = 1:nReps
        results{repIdx} = runSingleReplicate(repIdx, config, sampleGens, nSamples);
    end
    fprintf('\nSimulation complete (%.1f minutes)\n', toc/60);
end


%% =========================================================================
%% SINGLE REPLICATE
%% =========================================================================

function result = runSingleReplicate(repIdx, config, sampleGens, nSamples)
    seedBase = 12345;
    if isfield(config, 'rngSeedBase'), seedBase = config.rngSeedBase; end
    rng(seedBase + repIdx, 'twister');

    N = config.populationSize;
    nGenes = config.numGenes;
    nGens = config.numGenerations;

    % Starting population: fully connected at initConnectivity = 1, weights N(0, initSigma^2)
    % (Methods, "Game and evolution")
    connectivity = 1.0;
    if isfield(config, 'initConnectivity')
        connectivity = config.initConnectivity;
    end

    M = false(nGenes, nGenes, N);
    B = zeros(nGenes, nGenes, N, 'single');
    for i = 1:N
        M(:,:,i) = rand(nGenes) < connectivity;
        B(:,:,i) = randn(nGenes) * config.initSigma;
    end

    % Given starting wirings (config.initPopulationM, config.initPopulationB, nGenes x nGenes x n,
    % n <= N) replace the population after the loop above has drawn it, so the random stream
    % that follows is the same as without them; with n < N individual k copies given wiring
    % mod(k - 1, n) + 1
    if isfield(config, 'initPopulationM')
        nSeed = size(config.initPopulationM, 3);
        assert(nSeed <= N, 'initPopulationM holds %d individuals > N = %d', nSeed, N);
        seedIdx = mod(0:N-1, nSeed) + 1;       % identity when nSeed == N
        M = logical(config.initPopulationM(:,:,seedIdx));
        B = single(config.initPopulationB(:,:,seedIdx));
    end

    % K_M per gene and individual (nGenes x N): config.michaelisConstants per gene where given,
    % else michaelisConstant on every gene
    if isfield(config, 'michaelisConstants') && length(config.michaelisConstants) == nGenes
        K = repmat(single(config.michaelisConstants(:)), 1, N);
    else
        K = ones(nGenes, N, 'single') * config.michaelisConstant;
    end

    % Stored history
    historyCooperation = zeros(nGens, 1, 'single');
    historyParents = zeros(N, nGens, 'uint16');
    historyM = false(nGenes, nGenes, N, nSamples);
    historyB = zeros(nGenes, nGenes, N, nSamples, 'single');
    historyGenerations = sampleGens;
    sampleIdx = 1;

    logInterval = max(50, floor(nGens / 10));

    % Optional per-individual cooperation at the generations in config.snapshotCoopGenerations
    if isfield(config, 'snapshotCoopGenerations') && ~isempty(config.snapshotCoopGenerations)
        snapshotCoopGens = config.snapshotCoopGenerations;
        nCoopSnaps = numel(snapshotCoopGens);
        snapshotCoopVectors = zeros(N, nCoopSnaps);
        coopSnapIdx = 1;
    else
        snapshotCoopGens = [];
        nCoopSnaps = 0;
        snapshotCoopVectors = [];
        coopSnapIdx = 1;
    end

    % Evolution loop
    for gen = 1:nGens
        % Store the population before it plays
        if sampleIdx <= nSamples && gen == sampleGens(sampleIdx)
            historyM(:,:,:,sampleIdx) = M;
            historyB(:,:,:,sampleIdx) = B;
            sampleIdx = sampleIdx + 1;
        end

        % Encounters and fitness
        [fitness, cooperationRate, cooperationVec] = evaluateFitness(M, B, K, config);

        % Per-individual cooperation at the requested generations
        if nCoopSnaps > 0 && coopSnapIdx <= nCoopSnaps && gen == snapshotCoopGens(coopSnapIdx)
            snapshotCoopVectors(:, coopSnapIdx) = cooperationVec;
            coopSnapIdx = coopSnapIdx + 1;
        end

        historyCooperation(gen) = cooperationRate;

        % Parents (selection weight of Methods "Game and evolution"; within groups and fission in selectParents)
        parents = selectParents(fitness, config);
        historyParents(:, gen) = parents;

        % Clonal reproduction, then mutation (config.mutationOperator, if set, replaces
        % mutatePopulation)
        M_new = M(:,:,parents);
        B_new = B(:,:,parents);
        K_new = K(:, parents);
        if isfield(config, 'mutationOperator') && ~isempty(config.mutationOperator)
            [M, B, K] = config.mutationOperator(M_new, B_new, K_new, config);
        else
            [M, B, K] = mutatePopulation(M_new, B_new, K_new, config);
        end

        % Progress
        if mod(gen, logInterval) == 0 || gen == 1 || gen == nGens
            fprintf('  Rep %d | Gen %d/%d | Coop: %.1f%%\n', ...
                repIdx, gen, nGens, cooperationRate * 100);
        end
    end

    % Results of the replicate
    result.finalNetworkM = M;
    result.finalNetworkB = B;
    result.finalNetworkK = K;
    result.historyCooperation = historyCooperation;
    result.historyParents = historyParents;
    result.historyM = historyM;
    result.historyB = historyB;
    result.historyGenerations = historyGenerations;
    result.snapshotCoopVectors = snapshotCoopVectors;
    result.snapshotCoopGenerations = snapshotCoopGens;

    fprintf('  Rep %d COMPLETE | Final Coop: %.1f%%\n', repIdx, historyCooperation(end)*100);
end


%% =========================================================================
%% ENCOUNTERS AND FITNESS
%% =========================================================================

function [fitness, cooperationRate, cooperationVec, enc] = evaluateFitness(M, B, K, config, extMask, W_ext, K_ext)
    % extMask, W_ext, K_ext (optional): focal i meets the fixed external wiring W_ext where
    % extMask(i) is true
    hasExt = nargin >= 5 && ~isempty(extMask);
    N = size(M, 3);
    nOpp = config.numOpponents;

    enc = [];

    % Regulatory matrices W = M .* B (Methods, "Model") and the encounters
    W = double(M) .* double(B);  % nGenes x nGenes x N

    % Partners given by the caller (config.precomputedOpponents), one per individual
    if isfield(config, 'precomputedOpponents') && ~isempty(config.precomputedOpponents)
        persistOpp = config.precomputedOpponents;
        focalIdx = (1:N)';
        oppIdx = persistOpp(:);
        pairCount = N;
    else
        % numOpponents partners per individual, drawn at random (Methods, "Game and evolution")
        maxPairs = N * nOpp;
        focalIdx = zeros(maxPairs, 1);
        oppIdx = zeros(maxPairs, 1);
        pairCount = 0;

        for i = 1:N
            opponents = selectOpponents(i, N, nOpp, config);
            for k = 1:length(opponents)
                j = opponents(k);
                if j == i
                    continue;
                end
                pairCount = pairCount + 1;
                focalIdx(pairCount) = i;
                oppIdx(pairCount) = j;
            end
        end
        focalIdx = focalIdx(1:pairCount);
        oppIdx = oppIdx(1:pairCount);
    end
    nPairs = pairCount;

    % Development of every pair in one batch
    W_focal = W(:,:,focalIdx);  % nGenes x nGenes x nPairs
    W_opp = W(:,:,oppIdx);      % nGenes x nGenes x nPairs

    % K_M of each pair's two networks (nGenes x nPairs)
    km_focal = double(K(:, focalIdx));
    km_opp = double(K(:, oppIdx));

    % External partner: focal p meets W_ext instead of individual oppIdx(p), with the same dynamics
    % and the same random start as every other pair
    if hasExt
        assert(nPairs == N && isequal(focalIdx(:), (1:N)') && nOpp == 1, ...
            'externalDefector requires focal payoff, one opponent, one pair per focal');
        extPairs = find(extMask(:))';
        if ~isempty(extPairs)
            W_opp(:,:,extPairs) = repmat(W_ext, 1, 1, numel(extPairs));
            km_opp(:,extPairs) = repmat(double(K_ext(:)), 1, numel(extPairs));
        end
    end

    % One round of the game
    [coopFocal, coopOpp] = simulateInteractionBatch(W_focal, W_opp, config, km_focal, km_opp);
    coopFocal = double(coopFocal(:));
    coopOpp = double(coopOpp(:));

    % Fitness is the payoff of the game played as focal; games played as another's partner count
    % only for that other individual (Methods, "Game and evolution")
    fitness = zeros(N, 1);
    cooperation = zeros(N, 1);
    counts = zeros(N, 1);

    for p = 1:nPairs
        [payI, ~] = computePayoffs(coopFocal(p), coopOpp(p), config);
        i = focalIdx(p);
        fitness(i) = fitness(i) + payI;
        cooperation(i) = cooperation(i) + coopFocal(p);
        counts(i) = counts(i) + 1;
    end

    validMask = counts > 0;
    fitness(validMask) = fitness(validMask) ./ nOpp;
    cooperation(validMask) = cooperation(validMask) ./ nOpp;
    cooperationRate = mean(cooperation);
    cooperationVec = cooperation;

    if hasExt
        % Population cooperation counts encounters within the population only: defecting back at the
        % external defector is retaliation, not a failure to cooperate
        isExt = logical(extMask(:));
        if any(isExt)
            cooperationVec(isExt) = NaN;
            cooperationRate = mean(cooperation(~isExt));
        end
        enc.external = isExt;
        enc.oppDefected = ~logical(coopOpp(:));
        enc.focalCoop = logical(coopFocal(:));
    end
end


%% =========================================================================
%% DEVELOPMENT OF THE ENCOUNTERS
%% =========================================================================
% Every pair of a generation develops in one batch (pagemtimes). An encounter: both networks
% start from a random low state x_j ~ U(0, 0.01) and develop for at most config.developmentSteps
% Euler steps (300 in the paper). Every stabilityCheckInterval (40) steps each network is tested
% against ||x(t) - x(t-40)||_2 / max_i x_i(t) < stabilityThreshold (0.01). A network that passes
% while its partner does not is recorded and keeps developing, and is recorded again at every
% later check it passes; when both pass the same check, both are recorded there and the pair
% stops. A network's action is read from its last record (at the step limit if it never passed):
% cooperate when x_C > x_D, else defect.

function [coopFocal, coopOpp] = simulateInteractionBatch(W_focal, W_opp, config, km_focal, km_opp)

    nPairs = size(W_focal, 3);
    nGenes = config.numGenes;
    maxSteps = config.developmentSteps;
    dt = config.timeStep;
    intRate = config.sensorIntegrationRate;
    threshold = config.stabilityThreshold;
    interval = config.stabilityCheckInterval;
    gC = config.geneCooperator;
    gD = config.geneDefector;
    gS1 = config.geneSensor1;
    gS2 = config.geneSensor2;

    muDt = ones(nGenes, 1) * config.degradationRate;   % mu per gene (Delta t applied in the Euler step)
    thresholdSq = threshold^2;

    % The batched products below need at least two genes
    assert(nGenes >= 2, 'simulateInteractionBatch requires nGenes >= 2');

    useSeparate = strcmp(config.sensingMode, 'separate');
    % Output passed to the partner (Methods, "Model"): hillCoefficient 1 is saturating output
    % x_C / (x_C + x_D); Inf (or absent) graded output, x_C and x_D as they are
    hillCoeff = Inf;
    if isfield(config, 'hillCoefficient'), hillCoeff = config.hillCoefficient; end
    useHill = isfinite(hillCoeff);
    % Noise in gene expression (Table S5): standard deviation noiseSigma, 0 in the standard runs
    noiseSigma = 0;
    if isfield(config, 'noiseSigma'), noiseSigma = config.noiseSigma; end

    % Random low start x_j ~ U(0, 0.01) (Methods, "Model"), nGenes x nPairs
    xF = rand(nGenes, nPairs) * 0.01;   % focal networks
    xO = rand(nGenes, nPairs) * 0.01;   % partner networks

    % Convergence records (above simulateInteractionBatch)
    xF_snap = xF;  % state each action is read from
    xO_snap = xO;
    stableF = false(1, nPairs);
    stableO = false(1, nPairs);
    lastCheck = 0;

    % Optional constant input to every gene of every network (config.geneInputBias)
    useInputBias = isfield(config, 'geneInputBias') && ~isempty(config.geneInputBias);
    if useInputBias
        inputBias = double(config.geneInputBias(:));
        assert(numel(inputBias) == nGenes, 'geneInputBias must be numGenes x 1');
    end

    % Pairs still developing; a pair leaves the batch once both networks converge
    aIdx = 1:nPairs;               % indices of the pairs still developing
    nActive = nPairs;
    % Their wirings, states and K_M
    Wa_f = W_focal; Wa_o = W_opp;
    xa_f = xF;      xa_o = xO;
    xa_f_old = xF;  xa_o_old = xO;
    kma_f = km_focal;  kma_o = km_opp;  % nGenes x pairs still developing

    for t = 1:maxSteps
        % Regulatory input sum_j W_ij x_j for the pairs still developing
        inputF = squeeze(pagemtimes(Wa_f, reshape(xa_f, nGenes, 1, nActive)));
        inputO = squeeze(pagemtimes(Wa_o, reshape(xa_o, nGenes, 1, nActive)));

        % With one pair left, squeeze returns a row; keep a column
        if nActive == 1
            inputF = inputF(:);
            inputO = inputO(:);
        end

        % What each network's sensors receive of its partner's current outputs
        if useHill
            % Saturating output x_C^n / (x_C^n + x_D^n) and its complement,
            % n = hillCoefficient (1 in the paper)
            oC_n = xa_o(gC,:).^hillCoeff;
            oD_n = xa_o(gD,:).^hillCoeff;
            denom_o = oC_n + oD_n + 1e-12;  % defined when x_C = x_D = 0
            opp_CF = oC_n ./ denom_o;
            opp_DF = oD_n ./ denom_o;
            fC_n = xa_f(gC,:).^hillCoeff;
            fD_n = xa_f(gD,:).^hillCoeff;
            denom_f = fC_n + fD_n + 1e-12;
            opp_CO = fC_n ./ denom_f;
            opp_DO = fD_n ./ denom_f;
        else
            % Graded output: x_C and x_D as they are
            opp_CF = xa_o(gC,:);
            opp_DF = xa_o(gD,:);
            opp_CO = xa_f(gC,:);
            opp_DO = xa_f(gD,:);
        end

        % Sensor input (Eq. 4): separate sensing S_C = g x_C, S_D = g x_D; contrast sensing
        % S_C = -S_D = g (x_C - x_D)
        if useSeparate
            inputF(gS1,:) = inputF(gS1,:) + intRate * opp_CF;
            inputF(gS2,:) = inputF(gS2,:) + intRate * opp_DF;
            inputO(gS1,:) = inputO(gS1,:) + intRate * opp_CO;
            inputO(gS2,:) = inputO(gS2,:) + intRate * opp_DO;
        else
            sigF = intRate * (opp_CF - opp_DF);
            sigO = intRate * (opp_CO - opp_DO);
            inputF(gS1,:) = inputF(gS1,:) + sigF;
            inputF(gS2,:) = inputF(gS2,:) - sigF;
            inputO(gS1,:) = inputO(gS1,:) + sigO;
            inputO(gS2,:) = inputO(gS2,:) - sigO;
        end

        % Constant input (config.geneInputBias), added after sensing and before f, as in
        % simulate_pair
        if useInputBias
            inputF = inputF + inputBias;
            inputO = inputO + inputBias;
        end

        % Production f(u) = max(u, 0) / (K_M + max(u, 0)) (Eq. 3)
        inputF = max(inputF, 0);
        inputO = max(inputO, 0);
        prodF = inputF ./ (kma_f + inputF);
        prodO = inputO ./ (kma_o + inputO);

        % Euler step of dx_i/dt = f(u_i) - mu x_i; with noiseSigma > 0 an Euler-Maruyama step
        % with noise of standard deviation noiseSigma
        xa_f = xa_f + dt * (prodF - muDt .* xa_f);
        xa_o = xa_o + dt * (prodO - muDt .* xa_o);
        if noiseSigma > 0
            xa_f = xa_f + noiseSigma * sqrt(dt) * randn(size(xa_f));
            xa_o = xa_o + noiseSigma * sqrt(dt) * randn(size(xa_o));
        end

        % Expression kept non-negative
        xa_f = max(xa_f, 0);
        xa_o = max(xa_o, 0);

        % Convergence check every stabilityCheckInterval steps (above simulateInteractionBatch); 1e-9
        % keeps the test defined when every x_i = 0
        if t - lastCheck >= interval
            maxF = max(max(xa_f, [], 1), 1e-9);
            maxO = max(max(xa_o, [], 1), 1e-9);
            diffF = sum((xa_f - xa_f_old).^2, 1) ./ maxF.^2;
            diffO = sum((xa_o - xa_o_old).^2, 1) ./ maxO.^2;
            convF = diffF < thresholdSq;
            convO = diffO < thresholdSq;
            bothConv = convF & convO;

            % Both networks converged: record both
            if any(bothConv)
                snapIdx = aIdx(bothConv);
                xF_snap(:, snapIdx) = xa_f(:, bothConv);
                xO_snap(:, snapIdx) = xa_o(:, bothConv);
                stableF(snapIdx) = true;
                stableO(snapIdx) = true;
            end
            % One network converged: record it, it keeps developing
            newF = convF & ~bothConv;
            if any(newF)
                fIdx = aIdx(newF);
                xF_snap(:, fIdx) = xa_f(:, newF);
                stableF(fIdx) = true;
            end
            newO = convO & ~bothConv;
            if any(newO)
                oIdx = aIdx(newO);
                xO_snap(:, oIdx) = xa_o(:, newO);
                stableO(oIdx) = true;
            end

            % Pairs where both converged leave the batch
            stillActive = ~bothConv;
            if any(bothConv)
                aIdx = aIdx(stillActive);
                nActive = length(aIdx);
                if nActive == 0, break; end
                Wa_f = Wa_f(:,:,stillActive);
                Wa_o = Wa_o(:,:,stillActive);
                xa_f = xa_f(:,stillActive);
                xa_o = xa_o(:,stillActive);
                xa_f_old = xa_f;
                xa_o_old = xa_o;
                kma_f = kma_f(:,stillActive);
                kma_o = kma_o(:,stillActive);
            else
                xa_f_old = xa_f;
                xa_o_old = xa_o;
            end
            lastCheck = t;
        end
    end

    % State at the step limit of the pairs still developing
    for ii = 1:length(aIdx)
        p = aIdx(ii);
        xF(:, p) = xa_f(:, ii);
        xO(:, p) = xa_o(:, ii);
    end

    % A network that never converged is read at the step limit
    xF_snap(:, ~stableF) = xF(:, ~stableF);
    xO_snap(:, ~stableO) = xO(:, ~stableO);

    % Action: cooperate when x_C > x_D (a tie is defection)
    coopFocal = xF_snap(gC,:) > xF_snap(gD,:);   % logical 1 x nPairs
    coopOpp = xO_snap(gC,:) > xO_snap(gD,:);     % logical 1 x nPairs
end


%% =========================================================================
%% SELECTION
%% =========================================================================

function parents = selectParents(fitness, config)
    N = length(fitness);
    fitness = double(fitness);

    % Groups (SI Extended Materials and Methods, "Group selection"): each slot's parent is drawn
    % by the selection weight from its own group; groups are contiguous blocks of N_g slots; then one fission
    isGroup = isfield(config, 'groupSize') && config.groupSize > 0;
    if isGroup
        N_g = config.groupSize;
        if mod(N, N_g) ~= 0
            error('evolve_population:BadGroupSize', ...
                'populationSize (%d) must be divisible by groupSize (%d)', N, N_g);
        end
        G = N / N_g;
        parents = zeros(N, 1, 'uint16');
        for g = 1:G
            slots = ((g-1)*N_g + 1) : (g*N_g);
            local_fit = max(fitness(slots), 0) .^ config.selectionPower;
            local_fit = max(local_fit, 1e-10);
            probs = local_fit / sum(local_fit);
            cumProbs = cumsum(probs);
            cumProbs(end) = 1.0;
            for k = 1:N_g
                chosen = find(cumProbs >= rand(), 1, 'first');
                parents(slots(k)) = slots(chosen);
            end
        end
        parents = applyFission(parents, fitness, config, N_g, G, ...
            strcmp(config.groupSelectionMode, 'random_split'));
        return;
    end

    % Ring: local replacement; well-mixed (or no spatialStructure): the whole population
    isRing = isfield(config, 'spatialStructure') && ...
             strcmp(config.spatialStructure, 'ring');

    if isRing
        % Local replacement (SI Extended Materials and Methods, "Spatial reproduction"): slot i's
        % parent is drawn by the selection weight from its 2k ring neighbors and itself, {i-k, ..., i+k}
        nh = config.neighborhoodSize;
        offsets = -nh:nh;                      % includes self
        parents = zeros(N, 1, 'uint16');
        for i = 1:N
            neighborhood = mod(i - 1 + offsets, N) + 1;
            local_fit = max(fitness(neighborhood), 0) .^ config.selectionPower;
            local_fit = max(local_fit, 1e-10);
            probs = local_fit / sum(local_fit);
            cumProbs = cumsum(probs);
            cumProbs(end) = 1.0;
            chosen = find(cumProbs >= rand(), 1, 'first');
            parents(i) = neighborhood(chosen);
        end
        return;
    end

    % Well-mixed: parents drawn by the selection weight from the whole population, weights
    % max(fitness, 0)^zeta floored at 1e-10
    fitnessPower = max(fitness, 0) .^ config.selectionPower;
    fitnessPower = max(fitnessPower, 1e-10);
    probs = fitnessPower / sum(fitnessPower);
    cumProbs = cumsum(probs);
    cumProbs(end) = 1.0;
    parents = zeros(N, 1, 'uint16');
    for i = 1:N
        parents(i) = find(cumProbs >= rand(), 1, 'first');
    end
end


function parents = applyFission(parents, fitness, config, N_g, G, sham)
% One group fission per generation (groupSelectionMode = 'fission'), the
% generational form of Traulsen & Nowak (2006) with fixed group sizes. Acts on
% the next generation's parent vector after within-group reproduction:
%   1. splitting group h drawn with probability proportional to the mean
%      selection weight of the individuals now in it (max(fitness, 0)^selectionPower,
%      floored at 1e-10); sham = true (the random-split control, 'random_split'):
%      h drawn uniformly;
%   2. target group g drawn uniformly among the other G - 1 groups;
%   3. h's N_g next-generation entries randomly permuted and halved; daughter A
%      (first half, each entry twice) fills h's block, daughter B (second half,
%      each entry twice) fills g's block, replacing g.
% Group selection and its random-split control take the same draws: one rand() for h, one for g,
% one randperm.
    w = max(fitness, 0) .^ config.selectionPower;
    w = max(w, 1e-10);
    if sham
        groupCum = cumsum(ones(1, G) / G);
    else
        gw = mean(reshape(w, N_g, G), 1);
        groupCum = cumsum(gw / sum(gw));
    end
    groupCum(end) = 1.0;
    h = find(groupCum >= rand(), 1, 'first');
    others = [1:h-1, h+1:G];
    g = others(min(G - 1, floor(rand() * (G - 1)) + 1));
    hSlots = (h-1)*N_g + (1:N_g);
    gSlots = (g-1)*N_g + (1:N_g);
    entries = parents(hSlots);
    entries = entries(randperm(N_g));
    half = N_g / 2;
    parents(hSlots) = repelem(entries(1:half), 2);
    parents(gSlots) = repelem(entries(half+1:end), 2);
end


%% =========================================================================
%% MUTATION
%% =========================================================================

function [M, B, K] = mutatePopulation(M, B, K, config)
% Each offspring receives Poisson(mutationRate) mutations, each a toggle of one connection of M
% (probability mutationProbArchitecture) or a weight change of N(0, mutationStdWeight^2)
% (Methods, "Game and evolution"); K is returned unchanged.
    [nGenes, ~, N] = size(M);

    for i = 1:N
        nMutations = poissrnd(config.mutationRate);

        for m = 1:nMutations
            if rand() < config.mutationProbArchitecture
                g1 = randi(nGenes); g2 = randi(nGenes);
                M(g1, g2, i) = ~M(g1, g2, i);
            else
                g1 = randi(nGenes); g2 = randi(nGenes);
                B(g1, g2, i) = B(g1, g2, i) + randn() * config.mutationStdWeight;
            end
        end
    end
end


%% =========================================================================
%% HELPER FUNCTIONS
%% =========================================================================

function opponents = selectOpponents(i, N, nInteractions, config)
    if isfield(config, 'spatialStructure') && strcmp(config.spatialStructure, 'ring')
        nh = config.neighborhoodSize;
        offsets = [-nh:-1, 1:nh];
        neighbors = mod(i - 1 + offsets, N) + 1;
        neighbors(neighbors == i) = [];
        neighbors = unique(neighbors);
        opponents = neighbors(randperm(length(neighbors), min(nInteractions, length(neighbors))));
    elseif isfield(config, 'groupPairingMode') && strcmp(config.groupPairingMode, 'within') && ...
           isfield(config, 'groupSize') && config.groupSize > 0
        % Within-group play: the partner is drawn from the focal's own group
        N_g = config.groupSize;
        g = ceil(double(i) / N_g);
        groupSlots = ((g-1)*N_g + 1) : (g*N_g);
        groupSlots(groupSlots == i) = [];
        opponents = groupSlots(randperm(length(groupSlots), min(nInteractions, length(groupSlots))));
    else
        available = setdiff(1:N, i);
        opponents = available(randperm(length(available), min(nInteractions, length(available))));
    end
end


function [payI, payJ] = computePayoffs(coopI, coopJ, config)
    if coopI && coopJ
        payI = config.payoffCC; payJ = config.payoffCC;
    elseif coopI && ~coopJ
        payI = config.payoffCD; payJ = config.payoffDC;
    elseif ~coopI && coopJ
        payI = config.payoffDC; payJ = config.payoffCD;
    else
        payI = config.payoffDD; payJ = config.payoffDD;
    end
end


%% =========================================================================
%% VALIDATION
%% =========================================================================

function validateConfig(config)
    required = {'populationSize', 'numGenerations', 'numReplicates', ...
                'numGenes', 'geneCooperator', 'geneDefector', ...
                'selectionPower', 'mutationRate', 'developmentSteps', ...
                'sensingMode', 'payoffMode', 'selectionAlgorithm'};

    for i = 1:length(required)
        if ~isfield(config, required{i})
            error('evolve_population:MissingConfig', 'Required config field missing: %s', required{i});
        end
    end

    validSensing = {'separate', 'contrast'};
    if ~ismember(config.sensingMode, validSensing)
        error('evolve_population:InvalidConfig', 'sensingMode must be ''separate'' or ''contrast'', got ''%s''', config.sensingMode);
    end

    if ~strcmp(config.selectionAlgorithm, 'cumulative')
        error('evolve_population:InvalidConfig', 'selectionAlgorithm must be ''cumulative'', got ''%s''', config.selectionAlgorithm);
    end
end


function rejectUnusedOptions(config)
% Settings of model variants that no result in the paper uses stop the run instead of being
% ignored: each must be absent, empty or off (false, 0 or NaN), or hold the value named below.
    off = {'useGPU', 'staticSensorSignal', 'outputBlend', 'lateralAlpha', 'sensingDelay', ...
        'competitiveBinding', 'sensorCompetitive', 'multiSitePromoter', 'metabolicCost', ...
        'connectionCost', 'preWireSD_C', 'preWiredConnections', 'frozenConnections', ...
        'fixedMutualInhibition', 'frozenGenes', 'disableSensing', 'discreteOutput', ...
        'evolvableKM_genes', 'externalDefectorRate', 'recordEncounters', 'requireStability', ...
        'persistentPairingProb'};
    for k = 1:numel(off)
        if isfield(config, off{k}) && ~is_off(config.(off{k}))
            error('evolve_population:unusedOption', ...
                'config.%s is set, but no result in the paper uses it and it is not supported', off{k});
        end
    end
    named = {'numRounds', 1; 'sensorChannels', 'both'; 'payoffMode', 'focal'};
    for k = 1:size(named, 1)
        if isfield(config, named{k, 1}) && ~isequal(config.(named{k, 1}), named{k, 2})
            error('evolve_population:unusedOption', ...
                'config.%s must be %s; no result in the paper uses another value', ...
                named{k, 1}, mat2str(named{k, 2}));
        end
    end
    if isfield(config, 'degradationRates') && any(config.degradationRates(:) ~= config.degradationRate)
        error('evolve_population:unusedOption', ...
            'config.degradationRates must equal degradationRate on every gene');
    end
    if isfield(config, 'groupSelectionMode') && ...
       ~any(strcmp(config.groupSelectionMode, {'fission', 'random_split'}))
        error('evolve_population:unusedOption', ...
            'config.groupSelectionMode must be fission or random_split');
    end
end


function config = stored_spellings(config)
% A config read from a data file holds the older spellings of sensingMode and groupSelectionMode
% (stored_value); the run uses the current ones.
    for f = {'sensingMode', 'groupSelectionMode'}
        if isfield(config, f{1}), config.(f{1}) = stored_value(f{1}, config.(f{1})); end
    end
end


function tf = is_off(v)
% True for an empty value, or a single false, 0 or NaN.
    tf = isempty(v) || (isscalar(v) && (islogical(v) || isnumeric(v)) && (v == 0 || isnan(v)));
end
