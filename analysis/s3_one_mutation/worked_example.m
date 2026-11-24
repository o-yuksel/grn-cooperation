function worked_example()
% WORKED_EXAMPLE  The worked example's numbers for the member of Fig. 1B: the input its drawn
% mutation adds to D, D's settled input, its tipping input for cooperation at D, the readouts.
%
% Produces : main text, "One summed input ties responsiveness to mutation-sensitivity" (0.24,
%            0.73) and "Conditional cooperation stays within a mutation's reach of the threshold"
%            (0.73 against about 0.01); SI S3.3, "The worked example"
% Inputs   : the Fig. 1 member and its worst-case mutant, from fig1_member_and_mutant (its cache
%              output/figures/fig1_member_and_mutant.mat, built from
%              publication/saturating_contrast_100reps.mat, run_signaling_architectures)
%            publication/tipping_input/tipping_input_7000.mat   (tipping_input_flag, stage
%              'run'; read for one check)
% Outputs  : publication/example_numbers.txt   (the numbers, then the checks passed); printed
% Run      : worked_example   (about 30 s)
% Seed     : none (no random draws)
%
% The member: saturating contrast, replicate 81, generation 500; its mutant raises W(2, 3), the
% weight by which S_C regulates D, by delta (gene order C, D, S_C, S_D). The tipping input is
% tipping_input_flag's (shared/flip_threshold): the smallest constant input at one gene, added to
% both networks of the clone pair at every development step, that ends cooperation with the clone.
    out = fig1_member_and_mutant();
    sp = sim_params(4);                       % saturating contrast: contrast sensing, Hill 1, 600 steps
    cfg = out.config;
    assert(strcmp(stored_value('sensingMode', cfg.sensingMode), sp.sensingMode) && cfg.hillCoefficient == sp.hillCoefficient, ...
        'sensing config of the cache differs from sim_params(4)');
    assert(out.horizonTestMembership == sp.developmentSteps, 'horizon differs');
    W = out.W; Wm = out.W_worst_mutant; E = ones(4);
    gi = struct('C', 1, 'D', 2, 'S_C', 3, 'S_D', 4);
    nb = out.worst_mutant;
    assert(nb.ii == gi.D && nb.jj == gi.S_C, 'recorded slot is not D <- S_C');
    delta = nb.delta;
    assert(abs(Wm(gi.D, gi.S_C) - W(gi.D, gi.S_C) - delta) < 1e-6 && nnz(Wm ~= W) == 1, 'mutant is not a single-weight change');

    % The member's settled state x* with its clone; the input the mutation adds to D,
    % delta x x*_SC; D's settled input W(2, :) x*; the readouts x_C - x_D of the member with its
    % clone, of the mutant with its clone, and of the member against the mutant
    [~, fc] = simulate_pair(E, W, [], E, W, [], sp); xs = fc(1:4);
    x_SC = xs(gi.S_C);
    added_input = delta * x_SC;
    D_rest_input = W(gi.D, :) * xs(:);
    readout_member_clone = fc(1) - fc(2);
    [~, fmm] = simulate_pair(E, Wm, [], E, Wm, [], sp);
    readout_mutant_clone = fmm(1) - fmm(2);
    [~, fmem] = simulate_pair(E, W, [], E, Wm, [], sp);
    readout_member_vs_mutant = fmem(1) - fmem(2);
    tip = flip_threshold(@(d) flip_coop(W, sp, gi.D, d));
    % passes sets the input one network's mutation adds against a tipping input defined for an
    % input added to both networks; the flip itself is read from readout_mutant_clone
    passes = added_input >= tip;

    % Checks: the state, inputs and readouts against the cache and the stated values, and the
    % tipping inputs of member 1 of tipping_input_7000.mat against the stored ones; a failure is
    % reported in the control line and as a warning
    ctrl = {};
    ctrl{end+1} = chk('x_star vs cache', max(abs(xs(:) - out.testResult.x_star(:))), 1e-3);
    ctrl{end+1} = chk('x_star vs builder [1.668 1.2e-6 3.046 1.4e-6]', max(abs(xs(:) - [1.668; 1.2e-6; 3.046; 1.4e-6])), 1e-3);
    ctrl{end+1} = chk('added_input vs 0.731', abs(added_input - 0.731), 1e-3);
    ctrl{end+1} = chk('D_rest_input vs -0.286', abs(D_rest_input + 0.286), 1e-3);
    ctrl{end+1} = chk('readout_member_clone vs 1.668', abs(readout_member_clone - 1.668), 1e-3);
    ctrl{end+1} = chk('readout_mutant_clone vs -1.447', abs(readout_mutant_clone + 1.447), 1e-2);
    d = load(fullfile(grn_data_file('publication'), 'tipping_input', 'tipping_input_7000.mat'), 'S', 'R');
    s1 = d.S(1); sp1 = sim_params(s1.arch); thr1 = nan(4, 2);
    for i = 1:4
        for sgn = 1:2
            sv = 3 - 2 * sgn;
            thr1(i, sgn) = flip_threshold(@(dd) flip_coop(s1.W, sp1, i, sv * dd));
        end
    end
    stored = d.R{1}.thrC; fin = isfinite(stored);
    assert(isequal(fin, isfinite(thr1)), 'finite pattern of thrC differs');
    ctrl{end+1} = chk('thrC member 1 of tipping_input_7000 to 4 decimals', max(abs(thr1(fin) - stored(fin))), 5e-5);

    f = fopen(fullfile(grn_data_file('publication'), 'example_numbers.txt'), 'w');
    fprintf(f, 'x_SC = %.4f\nadded_input = %.4f\nD_rest_input = %.4f\n', x_SC, added_input, D_rest_input);
    fprintf(f, 'readout_member_clone = %.4f\nreadout_mutant_clone = %.4f\n', readout_member_clone, readout_mutant_clone);
    fprintf(f, 'readout_member_vs_mutant = %.4f\n', readout_member_vs_mutant);
    fprintf(f, 'tipping_input_D_plus = %.6f\npasses = %s\nhorizon = %d\n', tip, lower(string(passes)), sp.developmentSteps);
    fprintf(f, 'x_star = [%.4g %.4g %.4g %.4g]\nweight_change_D_from_SC = %.4f\n', xs, delta);
    fprintf(f, 'control = %s\n', strjoin(ctrl, '; '));
    fclose(f);
    type(fullfile(grn_data_file('publication'), 'example_numbers.txt'));
end

function s = chk(name, err, tol)
    if err <= tol, s = sprintf('%s ok (err %.2g)', name, err); return; end
    s = sprintf('%s FAILED (err %.3g > %.3g)', name, err, tol);
    warning('worked_example:control', 'CONTROL FAILED: %s (err %.3g > %.3g)', name, err, tol);
end

% sim_params: the same code as in tipping_input_flag
function sp = sim_params(m)
    T = {'separate', Inf; 'contrast', Inf; 'separate', 1; 'contrast', 1};
    sp = struct('sensingMode', T{m,1}, 'hillCoefficient', T{m,2}, 'discreteOutput', false, ...
        'michaelisConstant', 1.0, 'degradationRate', 0.1, 'timeStep', 1.0, 'developmentSteps', 600, ...
        'payoffDC', 0, 'sensorIntegrationRate', 1.0, 'externalSignal1', zeros(4, 1), ...
        'externalSignal2', zeros(4, 1));
end
