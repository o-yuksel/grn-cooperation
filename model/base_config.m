function config = base_config()
% BASE_CONFIG  The model, game and mutation parameters every setting of the game starts from.
%
% Implements : the parameters of Table S1 that every signaling architecture shares (Methods,
%              "Model" and "Game and evolution"); setting_config adds the gene layout and
%              publication_config sets the values the paper runs
% Inputs     : none
% Outputs    : config   struct of model, game and mutation parameters
% Seed       : none (no random numbers)
%
% Genes, in the order of Methods, "Model": C (cooperation output, 1), D (defection output, 2),
% S_C (cooperation sensor, 3), S_D (defection sensor, 4). M(i,j) = 1 when gene j regulates gene i.

    %% Population and evolution
    config.populationSize = 1000;     % N (Methods, "Game and evolution")
    config.numGenerations = 500;      % generations per run
    config.numReplicates = 100;       % replicate runs per setting

    %% Mutation (Methods, "Game and evolution")
    config.mutationRate = 0.005;              % lambda: Poisson mean of mutations per offspring
    config.mutationProbArchitecture = 0.2;    % share of mutations that toggle one connection of M
    config.mutationStdWeight = 0.2;           % sigma_mut: standard deviation of a weight mutation
                                              % (B = B + randn * 0.2), so the variance is 0.04,
                                              % the N(0, 0.04) of the Methods
    config.selectionPower = 1.0;              % selection power zeta (selection weight, Methods "Game and evolution")

    %% Development (Eq. 3)
    config.timeStep = 1.0;                    % Euler step Delta t
    config.degradationRate = 0.10;            % degradation rate mu
    config.michaelisConstant = 1.0;           % Michaelis constant K_M of f
    config.sensorIntegrationRate = 1.0;       % sensor gain g (Eq. 4)
    config.stabilityThreshold = 1e-2;         % convergence tolerance (SI Extended Materials and
                                              % Methods, "Developmental integration and convergence")
    config.stabilityCheckInterval = 40;       % steps between convergence checks (same passage)

    %% Prisoner's Dilemma payoffs (Methods, "Game and evolution")
    config.payoffCC = 3;    % R, reward for mutual cooperation
    config.payoffCD = 0;    % S, sucker's payoff
    config.payoffDC = 5.0;  % T, temptation
    config.payoffDD = 1;    % P, punishment for mutual defection

    % The payoffs must form a Prisoner's Dilemma: T > R > P > S and 2R > T + S
    assert(config.payoffDC > config.payoffCC, 'PD violated: T must be > R');
    assert(config.payoffCC > config.payoffDD, 'PD violated: R must be > P');
    assert(config.payoffDD > config.payoffCD, 'PD violated: P must be > S');
    assert(2*config.payoffCC > config.payoffDC + config.payoffCD, 'PD violated: 2R must be > T+S');

    %% Encounters
    config.numOpponents = 1;          % partners drawn per individual per generation

    %% Stored generations
    config.historySampleInterval = 10;                    % store the population every 10 generations
    config.alwaysSaveGenerations = [1, 50, 100, 200, 500]; % and at these generations
end
