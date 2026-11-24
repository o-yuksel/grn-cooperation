function config = setting_config(conditionName)
% SETTING_CONFIG  Gene layout and run settings, with the partner's output as the network passes it
% to its sensors named by the setting: 'continuous_output' (continuous expression, as in the model)
% or 'binary_output' (a binary sensed signal imposed). Sensing is on and the dynamics run in full in
% both.
%
% Implements : the four-gene network of Methods, "Model" (C, D, S_C, S_D), focal payoff and
%              the selection weight of Methods, "Game and evolution"; publication_config, the only caller in the
%              paper, replaces several values set here with the ones the paper runs
% Inputs     : conditionName   'continuous_output' (discreteOutput false) or 'binary_output'
%                              (discreteOutput true); developmentSteps 150, disableSensing false in
%                              both. The paper's scripts use 'continuous_output' only, through
%                              publication_config, which then sets discreteOutput false and
%                              developmentSteps 300. (Data files store older names of this
%                              setting; see stored_value.)
% Outputs    : config          struct for evolve_population
% Seed       : none (no random numbers)

    config = base_config();

    % Four genes in the order of Methods, "Model"; focal payoff (each individual's fitness is the
    % payoff of the game it plays as focal; Methods, "Game and evolution")
    config.numGenes = 4;
    config.geneNames = {'C', 'D', 'S_C', 'S_D'};
    config.geneCooperator = 1;
    config.geneDefector = 2;
    config.geneSensor1 = 3;
    config.geneSensor2 = 4;
    config.frozenGenes = [];
    config.sensingMode = 'contrast';
    config.payoffMode = 'focal';
    config.selectionAlgorithm = 'cumulative';
    config.requireStability = false;
    config.spatialStructure = 'wellmixed';

    % Run settings; publication_config replaces populationSize, selectionPower and numOpponents
    config.numReplicates = 100;
    config.populationSize = 200;
    config.numGenerations = 500;
    config.mutationRate = 0.005;
    config.mutationProbArchitecture = 0.2;
    config.mutationStdWeight = 0.2;
    config.selectionPower = 2.0;
    config.timeStep = 1.0;
    config.sensorIntegrationRate = 1.0;
    config.stabilityThreshold = 1e-2;
    config.stabilityCheckInterval = 40;
    config.numOpponents = 3;
    config.initConnectivity = 1.0;    % networks start fully connected (Methods, "Game and evolution")
    config.initSigma = 0.1;           % sigma_init of the starting weights (same passage)
    config.fixedMutualInhibition = false;
    config.historySampleInterval = 10;
    config.alwaysSaveGenerations = [1, 50, 100, 200, 500];

    % Setting name: the partner's output as passed to the sensors, continuous (D = 1) or binary
    % (D = 0); sensing on (I = 1) and the full dynamics (T = 1) in both
    switch conditionName
        case 'continuous_output', flags = '111';
        case 'binary_output',     flags = '011';
        otherwise
            error('setting_config:InvalidName', ...
                'Expected continuous_output or binary_output, got: %s', conditionName);
    end
    config.discreteOutput = flags(1) == '0';
    config.disableSensing = false;
    config.developmentSteps = 150;

    config.conditionName = conditionName;

    % Label of the three flags, stored with the run
    labels = {'discrete', 'independent', 'instantaneous'};
    states = {'IMPOSED', 'RELAXED'};
    desc = '';
    for k = 1:3
        stateIdx = str2double(flags(k)) + 1;
        desc = [desc, labels{k}, '=', states{stateIdx}];
        if k < 3, desc = [desc, ', ']; end
    end
    config.ablationDescription = desc;
end
