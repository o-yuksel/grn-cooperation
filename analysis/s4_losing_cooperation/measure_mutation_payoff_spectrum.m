function measure_mutation_payoff_spectrum(nMembers)
% MEASURE_MUTATION_PAYOFF_SPECTRUM  What does each single-weight mutant of a member earn against
% that member? Every mutant on the membership test's grid, for final saturating-contrast members.
%
% Produces : SI S4.4, "Members cooperate with nearly all of their own mutants that cooperate":
%            the shares of the 9,600 mutants of 60 members earning T, R, S and P; prints them
% Inputs   : nMembers (default 60, the paper's run)
%            publication/saturating_contrast_100reps.mat   (written by run_signaling_architectures)
% Outputs  : publication/mutation_payoff_spectrum.mat   (piAll, each mutant's payoff; nMembers_scored, the
%              number of members; Rp = R = 3); the console: the shares earning more than,
%              exactly, and less than R, and at each payoff value
% Run      : measure_mutation_payoff_spectrum   (3 members took about 1 min)
% Seed     : rng(7), the session's generator type: the order of the replicates and the
%            individuals tried in each
%
% Members: at most one per replicate, in a random order of replicates, the first of up to 6
% random final individuals that passes the membership test. Each mutant shifts one of the 16
% weights by one of the test's 10 nonzero steps in [-0.6, 0.6] and plays its member once from
% the test's start (0.005, 600 steps). A mutant earning more than the clone's R would be favored.

    %% Parameters
    if nargin < 1 || isempty(nMembers), nMembers = 60; end

    pub  = fullfile(grn_data_file('publication'));
    rng(7);

    %% Inputs
    files = dir(fullfile(pub, 'saturating_contrast_100reps.mat'));
    assert(~isempty(files), 'standard setting file not found');
    [~, mx] = max([files.datenum]);
    d = load(fullfile(pub, files(mx).name));
    cfg = d.config;

    params = struct('sensingMode', cfg.sensingMode, ...
        'hillCoefficient', cfg.hillCoefficient, ...
        'discreteOutput', cfg.discreteOutput, ...
        'michaelisConstant', cfg.michaelisConstant, ...
        'degradationRate', cfg.degradationRate, ...
        'timeStep', cfg.timeStep, 'developmentSteps', 600, ...
        'backend', 'simulate_pair');

    %% Compute
    % Members: at most one per replicate (see header)
    resid = {};
    nrep = numel(d.results);
    order = randperm(nrep);
    for r = order
        if numel(resid) >= nMembers, break; end
        B = d.results{r}.finalNetworkB; M = d.results{r}.finalNetworkM;
        na = size(B, 3);
        for a = randperm(na, min(6, na))
            W = double(M(:,:,a)) .* double(B(:,:,a));
            res = membership_test(W, params);
            if res.is_css
                resid{end+1} = W; %#ok<AGROW>
                break;
            end
        end
    end
    fprintf('members collected: %d\n', numel(resid));

    % Every single-weight mutant of each member on the membership test's grid. payoff_I_table is
    % the mutant's payoff, indexed by (member cooperates, mutant cooperates).
    P = struct('T',5,'R',3,'P',1,'S',0);
    payoff_I_table = [P.P, P.S; P.T, P.R];
    grid = css_mag_grid();
    Msq = ones(4); Kv = [];
    sp = struct('sensingMode', params.sensingMode, ...
        'hillCoefficient', params.hillCoefficient, ...
        'discreteOutput', params.discreteOutput, ...
        'michaelisConstant', params.michaelisConstant, ...
        'degradationRate', params.degradationRate, ...
        'timeStep', params.timeStep, 'developmentSteps', params.developmentSteps, ...
        'payoffDC', 0, 'sensorIntegrationRate', 1.0, ...
        'externalSignal1', zeros(4,1), 'externalSignal2', zeros(4,1));

    piAll = [];
    for k = 1:numel(resid)
        W = resid{k};
        for ii = 1:4
            for jj = 1:4
                for dlt = grid
                    if dlt == 0, continue; end
                    Wm = W; Wm(ii,jj) = Wm(ii,jj) + dlt;
                    o = simulate_pair(Msq, W, Kv, Msq, Wm, Kv, sp);
                    aR = double(o(1) == 'C'); aM = double(o(2) == 'C');
                    piAll(end+1) = payoff_I_table(aR+1, aM+1); %#ok<AGROW>
                end
            end
        end
        if mod(k, 10) == 0, fprintf('  %d/%d members\n', k, numel(resid)); end
    end

    %% Report
    Rp = P.R;
    fDefector = mean(piAll >  Rp);
    fNeu = mean(piAll == Rp);
    fWor = mean(piAll <  Rp);

    fprintf('\n=== single-weight mutants of members (n = %d mutants) ===\n', numel(piAll));
    fprintf('  invade   (pi_M >  3): %6.2f%%\n', 100*fDefector);
    fprintf('  neutral  (pi_M == 3): %6.2f%%\n', 100*fNeu);
    fprintf('  worse    (pi_M <  3): %6.2f%%\n', 100*fWor);
    u = unique(piAll);
    fprintf('  payoff values realized: %s\n', mat2str(u));
    for v = u
        fprintf('    pi_M = %g : %6.2f%%\n', v, 100*mean(piAll == v));
    end

    %% Save
    nMembers_scored = numel(resid); %#ok<NASGU>
    save(fullfile(pub, 'mutation_payoff_spectrum.mat'), 'piAll', 'nMembers_scored', 'Rp');
    fprintf('Saved -> %s\n', fullfile(pub, 'mutation_payoff_spectrum.mat'));
end

% The membership test's grid of weight steps: 11 points over +/- 3 sigma (sigma = 0.2).
function g = css_mag_grid()
    s = 0.2;
    g = linspace(-3*s, 3*s, 11);
end
