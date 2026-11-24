function membership_exits_split(nWire, nMut, condIdx)
% MEMBERSHIP_EXITS_SPLIT  Of the single-weight mutations that end membership, what share end
% cooperation (C1), and what share keep it and end only retaliation (C2)? Per architecture.
%
% Produces : SI S4.1 and Table S9: the share of mutations ending membership that leave
%            cooperation intact (printed as "C2 lost"; the second factor of u, which
%            measure_defector_encounter_frequency holds as fixed values); prints it itself
% Inputs   : nWire (members per architecture, default 25; SI S4.1 reports 20, 25, 25, 25),
%            nMut (mutations per member, 150), condIdx (architectures, 1:4)
%            publication/<architecture>_100reps.mat   (final networks and settings; written by
%              run_signaling_architectures)
% Outputs  : publication/member_exit_c1_c2_weight.mat   (counts per architecture)
% Run      : membership_exits_split()   (about 10 to 15 min)
% Seed     : rng(2600 + a) for architecture a (the session's generator type): the replicate
%            order, the members and the mutations; membership_test seeds each wiring itself
%            and restores this stream
%
% Members: one per replicate, in random replicate order, the first of its distinct wirings
% (rounded to 6 decimals, tried in random order) that passes the membership test. A mutation
% ends membership when the mutant fails the test; it ends cooperation when the mutant no longer
% cooperates with its clone, and otherwise ends only retaliation (some grid mutant of the
% mutant earns more than R = 3). The mutations are weight changes only (0.8 of the mutation
% distribution; SI S4.1, "without toggles"): one picks an entry at random and draws a change
% N(0, sigma^2); as in the runs (mutatePopulation in model/evolve_population.m changes B, and W = M .* B), the
% change applies
% only where the connection is present (M = 1), and on an absent connection W stays unchanged.

    if nargin < 1 || isempty(nWire),   nWire   = 25;  end
    if nargin < 2 || isempty(nMut),    nMut    = 150; end
    if nargin < 3 || isempty(condIdx), condIdx = 1:4; end

    pub  = fullfile(grn_data_file('publication'));

    SETTINGS = {
        'graded-separate',      'graded_separate_100reps.mat'
        'graded-contrast',      'graded_contrast_100reps.mat'
        'saturating-separate',  'saturating_separate_100reps.mat'
        'saturating-contrast',  'saturating_contrast_100reps.mat'
    };
    % A gene counts as off below OFF_TOL (a numerical tolerance the paper does not state)
    OFF_TOL = 1e-4;
    Msq = ones(4); Kv = [];
    OUT = struct();

    for a = condIdx(:)'
        d = load(fullfile(pub, SETTINGS{a,2})); c = d.config;
        sigma = c.mutationStdWeight;
        cp = struct('sensingMode', c.sensingMode, 'hillCoefficient', c.hillCoefficient, ...
            'discreteOutput', c.discreteOutput, 'michaelisConstant', c.michaelisConstant, ...
            'degradationRate', c.degradationRate, 'timeStep', c.timeStep, ...
            'developmentSteps', 600, 'backend', 'simulate_pair');
        sp = struct('sensingMode', c.sensingMode, 'hillCoefficient', c.hillCoefficient, ...
            'discreteOutput', c.discreteOutput, 'michaelisConstant', c.michaelisConstant, ...
            'degradationRate', c.degradationRate, 'timeStep', c.timeStep, ...
            'developmentSteps', 600, 'payoffDC', 0, 'sensorIntegrationRate', 1.0, ...
            'externalSignal1', zeros(4,1), 'externalSignal2', zeros(4,1));
        grid = linspace(-3*sigma, 3*sigma, 11);

        nTot=0; nEx=0; nC1=0; nC2=0; nB=0; nA=0; nS=0;
        marginC2 = [];    % readout of the most strongly defecting grid mutant met with cooperation, per event
        rng(2600 + a); nW = 0; t0 = tic;

        for k = randperm(numel(d.results))
            if nW >= nWire, break; end
            rep = d.results{k};
            M = rep.finalNetworkM; B = rep.finalNetworkB; N = size(M,3);
            flat = zeros(N,16);
            for j = 1:N
                flat(j,:) = reshape(double(M(:,:,j)) .* double(B(:,:,j)), 1, 16);
            end
            [~, ia] = unique(round(flat,6), 'rows', 'stable');
            Mall = double(M);
            took = false;
            for j = ia(randperm(numel(ia)))'
                if took, break; end
                W = reshape(flat(j,:), 4, 4);
                r0 = membership_test(W, cp);
                if ~r0.is_css, continue; end          % the membership test
                took = true; nW = nW + 1;

                % The member's clone-pair state x* and its genes that are off, the start and
                % the reference for judging whether a mutation ends cooperation
                [~, fCl, ~] = simulate_pair(Msq, W, Kv, Msq, W, Kv, sp);
                xr = fCl(1:4); off = xr(:) < OFF_TOL;
                spCont = sp; spCont.initialState1 = xr; spCont.initialState2 = fCl(5:8);

                for q = 1:nMut
                    i = randi(4); jj = randi(4);
                    % As in the runs (mutatePopulation): the change is to B, and W = M .* B, so a
                    % draw on an absent connection (M = 0) leaves W unchanged
                    dW = randn()*sigma;
                    Wm = W;
                    if Mall(i,jj,j) > 0.5
                        Wm(i,jj) = Wm(i,jj) + dW;
                    end
                    nTot = nTot + 1;
                    rm = membership_test(Wm, cp);
                    if rm.is_css, continue; end
                    nEx = nEx + 1;

                    % Ends cooperation (C1): basin if the mutant's clone pair still cooperates
                    % when started from the member's x*; active set if the set of genes off
                    % changes; smooth otherwise. Else ends retaliation only (C2).
                    if ~rm.is_weak
                        nC1 = nC1 + 1;
                        [~, fk, ~] = simulate_pair(Msq, Wm, Kv, Msq, Wm, Kv, spCont);
                        rCon = fk(1) - fk(2);
                        o2 = fk(1:4) < OFF_TOL;
                        if isfinite(rCon) && rCon > 0
                            nB = nB + 1;
                        elseif any(o2(:) ~= off)
                            nA = nA + 1;
                        else
                            nS = nS + 1;
                        end
                    else
                        nC2 = nC2 + 1;
                        m = worst_mutant(Wm, grid, sp, Msq, Kv);
                        if isfinite(m), marginC2(end+1) = m; end   %#ok<AGROW>
                    end
                end
            end
        end

        fprintf('\n=== %s, weight mutations: %d members, %d mutations, %d end membership (%.1f%%)  [%.1f min] ===\n', ...
            SETTINGS{a,1}, nW, nTot, nEx, 100*nEx/max(1,nTot), toc(t0)/60);
        if nEx > 0
            fprintf('  %-34s %7s\n', 'mutations ending membership', 'share');
            fprintf('  %-34s %6.1f%%\n', 'C1 lost (clone stops cooperating)', 100*nC1/nEx);
            fprintf('      %-30s %6.1f%%\n', 'of which basin',      100*nB/nEx);
            fprintf('      %-30s %6.1f%%\n', 'of which active-set', 100*nA/nEx);
            fprintf('      %-30s %6.1f%%\n', 'of which smooth',     100*nS/nEx);
            fprintf('  %-34s %6.1f%%   <- retaliation, intrinsic\n', ...
                'C2 lost (a mutant of it invades)', 100*nC2/nEx);
            if ~isempty(marginC2)
                fprintf('      invading mutant''s readout: median %.3f, p10 %.3f, p90 %.3f\n', ...
                    median(marginC2), quantile(marginC2,0.1), quantile(marginC2,0.9));
            end
        end
        OUT.(matlab.lang.makeValidName(SETTINGS{a,1})) = struct('nW',nW,'nTot',nTot,'nEx',nEx, ...
            'nC1',nC1,'nC2',nC2,'nBasin',nB,'nAct',nA,'nSmooth',nS,'marginC2',marginC2);
    end

    save(fullfile(pub, 'member_exit_c1_c2_weight.mat'), '-struct', 'OUT');
    fprintf('\nSaved: %s\n', fullfile(pub, 'member_exit_c1_c2_weight.mat'));
end


function m = worst_mutant(W, grid, sp, Msq, Kv)
% Readout of the grid mutant that defects most strongly while met with cooperation: the most negative mutant readout
% among the single-weight grid mutants against which the wiring W cooperates (readout > 0)
% and the mutant defects (readout <= 0); NaN if none.
    m = NaN; best = inf;
    for i = 1:4
        for j = 1:4
            for dg = grid
                if dg == 0, continue; end
                Wn = W; Wn(i,j) = Wn(i,j) + dg;
                [~, fs, ~] = simulate_pair(Msq, W, Kv, Msq, Wn, Kv, sp);
                rR = fs(1) - fs(2); rM = fs(5) - fs(6);
                if isfinite(rR) && isfinite(rM) && rR > 0 && rM <= 0 && rM < best
                    best = rM; m = rM;
                end
            end
        end
    end
end
