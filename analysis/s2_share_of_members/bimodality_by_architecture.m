function bimodality_by_architecture()
% BIMODALITY_BY_ARCHITECTURE  How bimodal is cooperation at generation 500 in each of the four
% signaling architectures? The bimodality coefficient behind the 0.95 success cutoff.
%
% Produces : the bimodality coefficients, the shares of populations at the two ends and between,
%            and the 20th percentile, behind the success cutoff of SI S2.4
% Inputs   : publication/<architecture>_100reps.mat   (100 populations per architecture; written
%            by run_signaling_architectures)
% Outputs  : publication/bimodality_by_architecture.txt and the same lines on the console, one per
%            architecture: n, BC, the shares of populations at or above 0.95, at or below 0.05
%            and in between, and the 20th percentile
% Run      : bimodality_by_architecture   (not timed)
% Seed     : none (no random draws)
%
% BC = (g^2 + 1) / (k + 3(n-1)^2 / ((n-2)(n-3))), with g the bias-corrected skewness and k the
% bias-corrected excess kurtosis of the n populations' mean cooperation in the last generation
% (Pfister et al. 2013); above 0.555, the uniform distribution's value, suggests bimodality.

    pub  = fullfile(grn_data_file('publication'));
    arch = {'graded_separate', 'graded_contrast', 'saturating_separate', 'saturating_contrast'};
    fid = fopen(grn_data_file('publication/bimodality_by_architecture.txt', 'new'), 'w');
    cleaner = onCleanup(@() fclose(fid));
    for k = 1:numel(arch)
        d = load(fullfile(pub, [arch{k} '_100reps.mat']));
        res = d.results;
        n = numel(res);
        % Each population's mean cooperation in its last generation (500)
        f = zeros(n, 1);
        for r = 1:n
            f(r) = res{r}.historyCooperation(end);
        end
        g = skewness(f, 0);
        ex = kurtosis(f, 0) - 3;
        bc = (g^2 + 1) / (ex + 3*(n-1)^2 / ((n-2)*(n-3)));
        txt = sprintf('%-20s n=%d BC=%.3f  >=0.95 %3.0f%%  <=0.05 %3.0f%%  between %3.0f%%  p20=%.3f\n', ...
            arch{k}, n, bc, 100*mean(f >= 0.95), 100*mean(f <= 0.05), ...
            100*mean(f > 0.05 & f < 0.95), prctile(f, 20));
        fprintf('%s', txt); fprintf(fid, '%s', txt);
    end
end
