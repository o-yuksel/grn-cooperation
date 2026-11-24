function share_interval_effective_n()
% SHARE_INTERVAL_EFFECTIVE_N  How much too narrow is an interval on the generation-500 share of
% members that treats every sampled wiring as independent? Compares it, per architecture, with an
% interval that resamples replicate populations.
%
% Produces : the comparison behind SI S2.3, "Generation-500 wirings are not independent draws"
%            (the paper's intervals resample replicates); printed only
% Inputs   : publication/share_of_members_cloud_sensing.mat   (gen500_share_of_members_per_replicate; per
%              architecture the share of members of each replicate, and N_SUB, the wirings
%              sampled per replicate)
% Outputs  : console table: per architecture the share, the number of replicates, both 95%
%              intervals and how many times wider the replicate interval is
% Run      : share_interval_effective_n   (about 1 min)
% Seed     : rng(42, 'twister') for the bootstrap (10,000 resamples)
%
% Populations are nearly clonal, so the replicate, not the wiring, is the independent unit. The
% two intervals: a Wilson interval with n = replicates x N_SUB wirings, and a percentile bootstrap
% over the per-replicate shares.

    pubDir = fullfile(grn_data_file('publication'));
    R = load(fullfile(pubDir, 'share_of_members_cloud_sensing.mat'));

    nBoot = 10000;
    rng(42, 'twister');   % bootstrap seed

    if isfield(R, 'N_SUB'), nSub = R.N_SUB; else, nSub = NaN; end
    fprintf('wirings sampled per replicate (N_SUB): %s\n\n', num2str(nSub));

    fprintf('%-22s %8s %8s   %-22s   %-22s  %7s\n', ...
        'architecture', 'reach', 'reps', 'pooled Wilson 95%', ...
        'replicate bootstrap 95%', 'x wider');
    for m = 1:numel(R.conds)
        rr = R.conds(m).rep_reach(:);
        rr = rr(~isnan(rr));
        nRep = numel(rr);
        phat = mean(rr);

        % (a) pooled Wilson, treating every sampled wiring as independent
        nPooled = nRep * nSub;
        [loW, hiW] = wilson(phat, nPooled);

        % (b) percentile bootstrap over replicates (the independent unit)
        bs = nan(nBoot, 1);
        for b = 1:nBoot
            bs(b) = mean(rr(randi(nRep, nRep, 1)));
        end
        loB = prctile(bs, 2.5); hiB = prctile(bs, 97.5);

        widthW = hiW - loW; widthB = hiB - loB;
        fprintf('%-22s %8.4f %8d   [%6.4f, %6.4f]       [%6.4f, %6.4f]      %6.1f\n', ...
            R.conds(m).name, phat, nRep, loW, hiW, loB, hiB, widthB / widthW);
    end

    fprintf(['\nThe right-hand column is how many times wider the replicate interval is.\n' ...
             'A value near 1 would mean pooling was harmless.\n']);
end

function [lo, hi] = wilson(p, n)
% Wilson score interval at 95% for a share p of n independent draws.
    z = 1.959963984540054;
    d = 1 + z^2/n;
    c = p + z^2/(2*n);
    s = z * sqrt(p*(1-p)/n + z^2/(4*n^2));
    lo = max(0, (c - s)/d);
    hi = min(1, (c + s)/d);
end
