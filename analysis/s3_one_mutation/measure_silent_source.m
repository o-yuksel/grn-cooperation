function measure_silent_source()
% MEASURE_SILENT_SOURCE  Among mutations that add under a thousandth of the tipping input, how often
% is the source gene off in the member and on in the mutant, with and without the loss?
%
% Produces : SI S3.3, "Mutations that add almost no input mostly act on connections from genes
%            that are off" (the source gene off before and on after)
% Inputs   : publication/tipping_input/tipping_input_7000.mat   (tipping_input_flag, stage 'run')
% Outputs  : publication/tipping_input/silent_source.txt   (the lines this run prints)
% Run      : measure_silent_source   (about 1 min)
% Seed     : none (no random draws)
%
% Per draw, q = |x_j dW_ij| over the member's tipping input at gene i of the same sign (q = 0 when
% the mutation adds no input or the gene has no tipping input of that sign); x_j is the source
% gene's expression in the member, with its clone (cooperation) or against its worst-case mutant
% H (retaliation). "Silent then on": below 0.01 in the member and at least max(10 times that,
% 0.01) in the mutant, in the same encounter. Printed per side: among losses with q < 1, and among
% draws with 0 < q < 1e-3 split by whether they cause the loss. Shares are weighted by the
% individuals carrying the wiring; n counts draws.
    outFile = fullfile(grn_data_file('publication'), 'tipping_input', 'silent_source.txt');
    L = load(grn_data_file('publication/tipping_input/tipping_input_7000.mat'));

    %% Report
    report = @() silent_report(L); %#ok<NASGU> called inside evalc
    rep = evalc('report();');
    fprintf('%s', rep);

    %% Save
    fid = fopen(outFile, 'w'); fprintf(fid, '%s', rep); fclose(fid);
end


function silent_report(L)
% The shares of each side; prints the report's lines.
    S = L.S; R = L.R; q = {[], []}; b = {[], []}; x0 = {[], []}; x1 = {[], []}; dw = []; w = [];
    % Per draw, side c = 1 cooperation (loss: the mutant's clone readout <= 0), c = 2 retaliation
    % (loss: the mutant's readout against H > 0): q, the loss flag, the source gene's expression in
    % the member (x0) and in the mutant (x1), dW, and the member's clone count as weight
    for k = 1:numel(S)
        s = S(k); o = R{k}; [i, ~] = ind2sub([4 4], s.slot(:)); K = numel(s.dW);
        e = {s.dW(:) .* o.xjC(:), s.dW(:) .* o.xjR(:)}; th = {o.thrC, o.thrR};
        out = {s.lost(:), s.rH(:) > 0}; xb = {o.xjC(:), o.xjR(:)}; xa = {o.xjCm(:), o.xjRm(:)};
        for c = 1:2
            qq = zeros(K, 1);
            for t = 1:K
                if e{c}(t) == 0 || ~isfinite(e{c}(t)), continue; end
                qq(t) = abs(e{c}(t)) / th{c}(i(t), 1 + (e{c}(t) < 0));
            end
            q{c} = [q{c}; qq]; b{c} = [b{c}; out{c}]; x0{c} = [x0{c}; xb{c}]; x1{c} = [x1{c}; xa{c}];
        end
        dw = [dw; s.dW(:)]; w = [w; repmat(s.count, K, 1)];
    end
    nm = {'cooperation', 'retaliation'};
    % Per side: among losses with q < 1, the share silent then on; then, among draws with
    % 0 < q < 1e-3, split by whether they cause the loss, the share silent then on and its parts
    for c = 1:2
        g = q{c} > 0 & q{c} < 1e-3;
        sil = x0{c} < 0.01 & x1{c} >= max(10 * x0{c}, 0.01);
        ms = q{c} < 1 & b{c}; Wm = w(ms);
        fprintf('%s MISSED losses (q<1 & loss): share source off then on %.4f | share q==0 %.4f | share 0<q<1e-2 %.4f\n', nm{c}, ...
            sum(Wm .* sil(ms)) / sum(Wm), sum(Wm .* (q{c}(ms) == 0)) / sum(Wm), sum(Wm .* (q{c}(ms) > 0 & q{c}(ms) < 1e-2)) / sum(Wm));
        for br = [true false]
            m = g & (b{c} == br); W = w(m);
            fprintf('%s, adds under 1e-3 of the tipping input, causes loss=%d: n=%d | source off then on %.3f | before<0.01 %.3f | after>=10x %.3f | |dW|<0.01 %.3f | median before %.4g after %.4g\n', ...
                nm{c}, br, sum(m), sum(W .* sil(m)) / sum(W), sum(W .* (x0{c}(m) < 0.01)) / sum(W), ...
                sum(W .* (x1{c}(m) >= 10 * x0{c}(m))) / sum(W), sum(W .* (abs(dw(m)) < 0.01)) / sum(W), median(x0{c}(m)), median(x1{c}(m)));
        end
    end
end
