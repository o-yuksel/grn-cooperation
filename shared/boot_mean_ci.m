function ci = boot_mean_ci(x, nBoot, seed)
%BOOT_MEAN_CI  Percentile bootstrap 95% interval of a mean.
%
% Implements : nBoot resamples with replacement of the non-NaN values; used by the tracking analyses and
%              fig6_load_cells
% Inputs     : x (values), nBoot, seed
% Outputs    : ci ([2.5 97.5] percentiles)
% Seed       : a RandStream(twister) built from the caller's seed argument; the global stream is not used

    x = x(~isnan(x)); n = numel(x);
    if n == 0, ci = [NaN NaN]; return; end
    st = RandStream('twister','Seed',seed);
    bm = nan(nBoot,1);
    for b = 1:nBoot, bm(b) = mean(x(randi(st, n, n, 1))); end
    ci = prctile(bm, [2.5 97.5]);
end
