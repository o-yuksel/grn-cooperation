function config = publication_config(conditionName, varargin)
% PUBLICATION_CONFIG  The paper's run settings, with optional overrides of single fields.
%
% Implements : the parameters of Table S1 and Methods, "Game and evolution" (N = 1,000, one
%              partner per generation, zeta = 1.0, T = 5.0, lambda = 0.005, mu = 0.10,
%              K_M = 1.0, up to 300 development steps of Delta t = 1.0), over setting_config
% Inputs     : conditionName   'continuous_output' in the paper's scripts (setting_config; data
%                              files store the older names of this setting, see stored_value)
%              varargin       name-value pairs, applied last, e.g. 'payoffDC', 4.0
% Outputs    : config          struct for evolve_population; with overrides, each pair is
%                              appended to config.conditionName as _<name>=<value>
% Seed       : none (no random numbers)
%
% Signaling architecture: sensingMode 'contrast' (contrast sensing) and hillCoefficient 1
% (saturating output) are set here, the standard setting; a caller sets sensingMode 'separate'
% (separate sensing) or hillCoefficient Inf (graded output) for the other three architectures.

    config = setting_config(conditionName);

    % Population, game and mutation (Methods, "Game and evolution"; Table S1)
    config.populationSize  = 1000;
    config.numOpponents    = 1;
    config.selectionPower  = 1.0;
    config.payoffDC        = 5.0;
    config.mutationRate    = 0.005;
    config.degradationRate = 0.10;
    config.numGenerations  = 500;
    config.michaelisConstant = 1.0;

    % Saturating contrast, continuous expression (Methods, "Model")
    config.sensingMode     = 'contrast';
    config.discreteOutput  = false;
    config.hillCoefficient = 1;

    % Per-gene Michaelis constants, where the base config defines them, all K_M = 1.0
    if isfield(config, 'michaelisConstants')
        config.michaelisConstants = ones(size(config.michaelisConstants));
    end

    % Up to 300 development steps (setting_config sets 150)
    if config.developmentSteps > 1
        config.developmentSteps = 300;
    end

    % Convergence checked every 40 steps (SI Extended Materials and Methods, "Developmental
    % integration and convergence")
    config.stabilityCheckInterval = 40;

    % Stored generations: every 25th, and these
    config.historySampleInterval  = 25;
    config.alwaysSaveGenerations  = [1, 50, 100, 200, 300, 500];

    % Overrides, applied last
    if mod(numel(varargin), 2) ~= 0
        error('publication_config:InvalidOverrides', ...
            'Overrides must be name-value pairs. Got %d arguments.', numel(varargin));
    end

    for i = 1:2:numel(varargin)
        config.(varargin{i}) = varargin{i+1};
    end

    % The stored name records the overrides
    if ~isempty(varargin)
        nPairs = numel(varargin) / 2;
        overrideParts = cell(1, nPairs);
        for i = 1:2:numel(varargin)
            val = varargin{i+1};
            if ischar(val) || isstring(val)
                valStr = char(val);
            elseif islogical(val)
                valStr = mat2str(val);
            else
                valStr = num2str(val);
            end
            overrideParts{(i+1)/2} = ['_' varargin{i} '=' valStr];
        end
        config.conditionName = [conditionName, overrideParts{:}];
    end
end
