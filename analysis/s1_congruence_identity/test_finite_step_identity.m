function test_finite_step_identity(nWire, condIdx, seedBase)
% TEST_FINITE_STEP_IDENTITY  Does the congruence identity hold for finite steps? On members,
% compares the readout change after a weight step dw of W_ij with that after input x*_j dw at gene i.
%
% Produces : the median ratio of the two changes, printed per architecture: SI S1.3, "At finite
%            steps the weight and signal responses also agree, in median"
% Inputs   : nWire      members per architecture (default 20, one per sampled replicate)
%            condIdx    architectures to run (default 1:4)
%            seedBase   base seed (default 8800)
%            publication/<architecture>_100reps.mat   (results, config; written by
%              run_signaling_architectures)
% Outputs  : printed per architecture and stratum (no file): draws, median ratio of the input-
%            driven to the weight-driven change with its IQR, the share within 5% of 1, the
%            Spearman correlation of the two, and the median ratio of the first-order prediction
%            v_i x*_j dw (v_i: the readout's response to an input at gene i) to the weight change
% Run      : test_finite_step_identity   and   test_finite_step_identity([], [], 7700)
%              (the two independent member samples of SI S1.3; not timed)
% Seed     : rng(seedBase + a), default generator, per architecture a; it orders the replicates
%            and the distinct wirings drawn; membership_test seeds each wiring itself
%
% x* is the member's settled state with its clone (600 steps). The mutant meets the unmutated
% wiring (the partner held fixed): a clone pair would move both sides under a weight mutation but
% only one under the added input. Steps are -2, -1, -0.5, 0.5, 1 and 2 times the run's mutation
% size sigma_mut. Strata: whether gene i is off at x* (input not positive) and whether the step
% changes which genes are on.

    if nargin < 1 || isempty(nWire),   nWire   = 20;    end
    if nargin < 2 || isempty(condIdx), condIdx = 1:4;   end
    if nargin < 3 || isempty(seedBase), seedBase = 8800; end

    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    pub  = fullfile(grn_data_file('publication'));

    SETTINGS = {
        'graded-separate',      'graded_separate_100reps.mat'
        'graded-contrast',      'graded_contrast_100reps.mat'
        'saturating-separate',  'saturating_separate_100reps.mat'
        'saturating-contrast',  'saturating_contrast_100reps.mat'
    };
    EPS = 1e-6; Msq = ones(4); Kv = [];
    % Cut below which a gene counts as off (a numerical tolerance the paper does not state): a
    % gene whose input is not positive decays as (1 - mu*dt)^600 to about 1.8e-30, never to
    % exactly 0, so a test for > 0 would call it on.
    OFFCUT = 1e-12;
    % Step sizes, in units of the mutation size sigma_mut
    STEPMULT = [-2 -1 -0.5 0.5 1 2];

    for a = condIdx(:)'
        d = load(fullfile(pub, SETTINGS{a,2})); c = d.config;
        mu = c.degradationRate; KM = c.michaelisConstant; sigma = c.mutationStdWeight;
        cp = struct('sensingMode', c.sensingMode, 'hillCoefficient', c.hillCoefficient, ...
            'discreteOutput', c.discreteOutput, 'michaelisConstant', KM, ...
            'degradationRate', mu, 'timeStep', c.timeStep, ...
            'developmentSteps', 600, 'backend', 'simulate_pair');
        sp = struct('sensingMode', c.sensingMode, 'hillCoefficient', c.hillCoefficient, ...
            'discreteOutput', c.discreteOutput, 'michaelisConstant', KM, ...
            'degradationRate', mu, 'timeStep', c.timeStep, ...
            'developmentSteps', 600, 'payoffDC', 0, 'sensorIntegrationRate', 1.0, ...
            'externalSignal1', zeros(4,1), 'externalSignal2', zeros(4,1));
        % Every measured pair develops the full 600 steps: the stop on convergence is turned
        % off, since the identity holds at the settled state
        spD = sp; spD.stabilityThreshold = 0; spD.stabilityCheckInterval = 1e9;

        dWp = []; dSp = []; dAn = []; offG = []; asChg = [];
        rng(seedBase + a); nW = 0; t0 = tic;

        for k = randperm(numel(d.results))
            if nW >= nWire, break; end
            rep = d.results{k};
            M = rep.finalNetworkM; B = rep.finalNetworkB; N = size(M,3);
            flat = zeros(N,16);
            for j = 1:N
                flat(j,:) = reshape(double(M(:,:,j)) .* double(B(:,:,j)), 1, 16);
            end
            [~, ia] = unique(round(flat,6), 'rows', 'stable');
            took = false;
            for j = ia(randperm(numel(ia)))'
                if took, break; end
                W = reshape(flat(j,:), 4, 4);
                % A member (membership test) with a positive readout against its clone at x*
                r0s = membership_test(W, cp); if ~r0s.is_css, continue; end
                [~, fB, ~] = simulate_pair(Msq, W, Kv, Msq, W, Kv, spD);
                xb = fB(1:4); xb = xb(:); r0 = fB(1) - fB(2);
                if ~(r0 > 0), continue; end
                z0 = mu*xb; u0 = KM * z0 ./ (1 - z0);
                isOn = (xb > OFFCUT) & (z0 < 1) & (u0 > 0);
                onPat = xb > OFFCUT;

                % v_i: first-order response of the readout to an input on gene i, by a central
                % difference at the settled state
                v = zeros(4,1);
                for kk = 1:4
                    e = zeros(4,1); e(kk) = EPS;
                    sP = spD; sP.externalSignal1 =  e;
                    sM = spD; sM.externalSignal1 = -e;
                    [~, f1, ~] = simulate_pair(Msq, W, Kv, Msq, W, Kv, sP);
                    [~, f2, ~] = simulate_pair(Msq, W, Kv, Msq, W, Kv, sM);
                    v(kk) = ((f1(1)-f1(2)) - (f2(1)-f2(2))) / (2*EPS);
                end
                took = true; nW = nW + 1;

                for i = 1:4
                    for jj = 1:4
                        if xb(jj) <= OFFCUT, continue; end   % from an off gene j, dw x*_j is ~0
                        for sm = STEPMULT
                            dw = sm * sigma;
                            Wm = W; Wm(i,jj) = Wm(i,jj) + dw;
                            [~, fW, ~] = simulate_pair(Msq, Wm, Kv, Msq, W, Kv, spD);
                            sS = spD; sS.externalSignal1 = zeros(4,1);
                            sS.externalSignal1(i) = xb(jj) * dw;
                            [~, fS, ~] = simulate_pair(Msq, W, Kv, Msq, W, Kv, sS);

                            xW = fW(1:4); xW = xW(:);
                            dWp(end+1)  = (fW(1)-fW(2)) - r0;               %#ok<AGROW>
                            dSp(end+1)  = (fS(1)-fS(2)) - r0;               %#ok<AGROW>
                            dAn(end+1)  = v(i) * xb(jj) * dw;               %#ok<AGROW>
                            offG(end+1) = ~isOn(i);                         %#ok<AGROW>
                            asChg(end+1)= any((xW > OFFCUT) ~= onPat);      %#ok<AGROW>
                        end
                    end
                end
            end
        end

        dWp = dWp(:); dSp = dSp(:); dAn = dAn(:);
        offG = logical(offG(:)); asChg = logical(asChg(:));
        fprintf('\n=== %s: %d members, %d draws (%.1f min) ===\n', ...
            SETTINGS{a,1}, nW, numel(dWp), toc(t0)/60);
        fprintf('  %-34s %6s %9s %17s %7s %9s %9s\n', 'stratum', 'n', ...
            'ratio_S', 'IQR_S', 'w/in5%', 'rho_S', 'ratio_A');
        report('all draws',                     dWp, dSp, dAn, true(size(dWp)));
        report('gene i ON, active set held',    dWp, dSp, dAn, ~offG & ~asChg);
        report('gene i ON, active set changed', dWp, dSp, dAn, ~offG &  asChg);
        report('gene i off (input not positive)', dWp, dSp, dAn, offG);
        report('  ... of those, |dW| > 1e-9',   dWp, dSp, dAn, offG & abs(dWp) > 1e-9);

        nz = offG & abs(dWp) > 1e-9;
        fprintf('  rectifier-off draws with a real W effect: %d / %d (%.1f%%)\n', ...
            nnz(nz), nnz(offG), 100*nnz(nz)/max(1,nnz(offG)));
        fprintf('  of those, signal probe nonzero: %d (%.1f%%) | analytic nonzero: %d (%.1f%%)\n', ...
            nnz(nz & abs(dSp) > 1e-9), 100*nnz(nz & abs(dSp) > 1e-9)/max(1,nnz(nz)), ...
            nnz(nz & abs(dAn) > 1e-9), 100*nnz(nz & abs(dAn) > 1e-9)/max(1,nnz(nz)));
    end
end


% One stratum: over draws whose weight-driven change exceeds 1e-9 in size, the median ratio of
% the input-driven to the weight-driven change, its IQR and the share within 5% of 1, and the
% median ratio of the first-order prediction; the Spearman correlation over all its draws. A
% stratum with fewer than 8 draws prints dashes.
function report(name, dW, dS, dA, sel)
    n = nnz(sel);
    if n < 8
        fprintf('  %-34s %6d %9s %17s %7s %9s %9s\n', name, n, '-', '-', '-', '-', '-');
        return;
    end
    w = dW(sel); s = dS(sel); A = dA(sel);
    ok = abs(w) > 1e-9;
    rS = NaN; rA = NaN; q1 = NaN; q3 = NaN; w5 = NaN;
    if nnz(ok) >= 8
        rr = s(ok) ./ w(ok);
        rS = median(rr); rA = median(A(ok) ./ w(ok));
        q1 = prctile(rr, 25); q3 = prctile(rr, 75);
        w5 = 100 * mean(abs(rr - 1) < 0.05);
    end
    pS = corr(s, w, 'type', 'Spearman');
    fprintf('  %-34s %6d %9.3f  [%6.3f,%6.3f] %6.1f%% %9.3f %9.3f\n', ...
        name, n, rS, q1, q3, w5, pS, rA);
end
