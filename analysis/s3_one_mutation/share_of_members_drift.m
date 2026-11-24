function out = share_of_members_drift(cfg)
% SHARE_OF_MEMBERS_DRIFT  Does the share of members still drift after generation 300? Per-replicate
% slopes, their mean with a cluster bootstrap, and each replicate's autocorrelation time.
%
% Produces : SI S3.2, "After generation 300, the share of members shows no resolved drift in three
%            architectures"; Table S3 (the slope rows)
% Inputs   : publication/member_share_by_generation_architectures.mat, _standard100.mat and
%              _satsep100.mat   (the share of members at every stored generation of each
%              replicate; share_of_members_by_generation)
%            publication/<the run each of those names>   (cooperation history, for success at
%              generation 300; run_signaling_architectures)
% Outputs  : publication/share_of_members_drift.mat   (out, one row per setting; metadata);
%            prints median and mean slopes with intervals, overall and split by success (at
%              generation 500 for the medians, 300 for the means), the window contrast (mean share
%              over generations 400 to 500 minus 300 to 400), and the lag-1 autocorrelation and
%              its time
% Run      : share_of_members_drift   (about 15 s)
% Seed     : rng(20260917), generator not named (MATLAB's default, Mersenne twister), for the
%            median bootstraps; RandStream('mt19937ar', 'Seed', 20260922) for the mean bootstraps
%
% Slope: least squares of one replicate's share of members on generation, over the stored
% generations from 300 on, per 100 generations. The mean slope equals the slope of the
% replicate-averaged share (every replicate has the same stored generations); its interval
% resamples whole replicates (10,000 times). Autocorrelation time tau = -dg / log(rho1): rho1 is
% the lag-1 autocorrelation of each replicate's deviation from its own mean over the stored
% generations from 100 on (the later of two stored one generation apart dropped); dg, their mean
% spacing.

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    if ~isfield(cfg,'dataDir'), cfg.dataDir = fullfile(grn_data_file('publication')); end
    if ~isfield(cfg,'from'),    cfg.from = 300;    end   % last 200 generations
    if ~isfield(cfg,'nBoot'),   cfg.nBoot = 10000; end
    if ~isfield(cfg,'rngSeed'), cfg.rngSeed = 20260917; end
    if ~isfield(cfg,'rngSeedMean'), cfg.rngSeedMean = 20260922; end
    if ~isfield(cfg,'fateGen'), cfg.fateGen = 300; end
    rng(cfg.rngSeed);
    % The mean bootstraps draw from a stream of their own, so the median intervals drawn from the
    % global stream do not depend on them
    meanStream = RandStream('mt19937ar', 'Seed', cfg.rngSeedMean);

    % Read in this order: the bootstraps consume the streams in file and setting order, so the
    % order fixes every interval; a missing file is skipped with a warning
    sources = {'member_share_by_generation_architectures.mat'; 'member_share_by_generation_standard100.mat'; ...
               'member_share_by_generation_satsep100.mat'};

    out = struct('name',{},'source',{},'sourceDate',{},'N',{},'nReps',{},'success',{}, ...
                 'slope100',{},'slopeMed',{},'slopeCI',{},'fracPos',{}, ...
                 'slopeMedSucc',{},'slopeCISucc',{},'slopeMedFail',{},'slopeCIFail',{}, ...
                 'winDiff',{},'winMed',{},'winCI',{},'rho1',{},'tau',{},'dg',{},'occEnd',{}, ...
                 'slopeMean',{},'slopeMeanCI',{},'slopeIQR',{},'succFate',{},'nSuccFate',{}, ...
                 'slopeMeanSuccFate',{},'slopeMeanCISuccFate',{}, ...
                 'slopeMeanFailFate',{},'slopeMeanCIFailFate',{},'transitions',{});

    for s = 1:numel(sources)
        f = fullfile(cfg.dataDir, sources{s});
        if ~exist(f,'file'), warning('missing %s', sources{s}); continue; end
        fi = dir(f);
        S = load(f, 'out'); res = S.out;
        for c = 1:numel(res)
            r = res(c);
            sg = double(r.gens(:))';
            share = r.occ;
            nReps = size(share,1);
            succ = r.finalCoop(:) >= 0.95;

            % ---- (a) per-replicate slope of the share of members, generations >= cfg.from ----
            w = sg >= cfg.from;
            gw = sg(w);
            slope = nan(nReps,1);
            for i = 1:nReps
                y = share(i,w)';
                ok = ~isnan(y);
                if nnz(ok) < 3, continue; end
                p = polyfit(gw(ok)', y(ok), 1);
                slope(i) = p(1) * 100;      % change in the share per 100 generations
            end

            % ---- (b) paired window contrast: mean share over [400, 500] minus [300, 400) ----
            wA = sg >= 300 & sg < 400;  wB = sg >= 400;
            winDiff = mean(share(:,wB),2,'omitnan') - mean(share(:,wA),2,'omitnan');

            % ---- (c) lag-1 autocorrelation of within-replicate deviation ----
            rep1 = [false, diff(sg) <= 1];
            keep = (sg >= 100) & ~rep1;
            gk = sg(keep);
            Y = share(:, keep);
            Y = Y - mean(Y, 2, 'omitnan');         % within-replicate deviation
            a = Y(:,1:end-1); b = Y(:,2:end);
            ok = ~isnan(a) & ~isnan(b);
            rho1 = sum(a(ok).*b(ok)) / sqrt(sum(a(ok).^2) * sum(b(ok).^2));
            dg = mean(diff(gk));
            tau = -dg / log(max(rho1, eps));

            k = numel(out) + 1;
            out(k).name   = r.name;
            out(k).source = sources{s};
            out(k).sourceDate = datestr(fi.datenum, 'yyyy-mm-dd HH:MM');
            out(k).N = r.N; out(k).nReps = nReps; out(k).success = mean(succ);
            out(k).slope100 = slope;
            [out(k).slopeMed, out(k).slopeCI] = boot_median(slope, cfg.nBoot);
            out(k).fracPos = mean(slope(~isnan(slope)) > 0);
            [out(k).slopeMedSucc, out(k).slopeCISucc] = boot_median(slope(succ), cfg.nBoot);
            [out(k).slopeMedFail, out(k).slopeCIFail] = boot_median(slope(~succ), cfg.nBoot);
            out(k).winDiff = winDiff;
            [out(k).winMed, out(k).winCI] = boot_median(winDiff, cfg.nBoot);
            out(k).rho1 = rho1; out(k).tau = tau; out(k).dg = dg;
            out(k).occEnd = mean(share(:,end),'omitnan');

            % ---- (a') mean slope, cluster bootstrap; split by fate at cfg.fateGen
            % (cooperation >= 0.95 then), which does not condition on the end point ----
            [out(k).slopeMean, out(k).slopeMeanCI] = boot_mean(slope, cfg.nBoot, meanStream);
            out(k).slopeIQR = prctile(slope(~isnan(slope)), [25 75]);
            sF = fate_at_gen(cfg.dataDir, r.file, r.finalCoop(:), cfg.fateGen);
            out(k).succFate = sF;
            out(k).nSuccFate = nnz(sF);
            [out(k).slopeMeanSuccFate, out(k).slopeMeanCISuccFate] = boot_mean(slope(sF), cfg.nBoot, meanStream);
            [out(k).slopeMeanFailFate, out(k).slopeMeanCIFailFate] = boot_mean(slope(~sF), cfg.nBoot, meanStream);
            out(k).transitions = [nnz(sF & succ), nnz(sF & ~succ), nnz(~sF & succ), nnz(~sF & ~succ)];
        end
    end

    % ------------------------------- report -------------------------------
    fprintf('\nDRIFT IN THE SHARE OF MEMBERS over the last 200 generations (gen >= %d)\n', cfg.from);
    fprintf('slope in SHARE-OF-MEMBERS UNITS PER 100 GENERATIONS; CI = percentile bootstrap of the median over replicates (%d draws)\n\n', cfg.nBoot);
    fprintf('%-24s %5s %5s %6s %9s %26s %7s %26s\n','setting','reps','N','succ','share g500','median slope /100 gen','frac+','window [400,500]-[300,400)');
    for k = 1:numel(out)
        o = out(k);
        fprintf('%-24s %5d %5d %5.0f%% %9.4f   %+8.5f [%+.5f,%+.5f] %6.2f   %+8.5f [%+.5f,%+.5f]\n', ...
            o.name, o.nReps, o.N, 100*o.success, o.occEnd, o.slopeMed, o.slopeCI, o.fracPos, o.winMed, o.winCI);
    end
    fprintf('\nSPLIT BY FATE (the pooled median is a mixture)\n');
    fprintf('%-24s %30s %30s\n','setting','successful reps','failed reps');
    for k = 1:numel(out)
        o = out(k);
        fprintf('%-24s   %+8.5f [%+.5f,%+.5f]   %+8.5f [%+.5f,%+.5f]\n', o.name, ...
            o.slopeMedSucc, o.slopeCISucc, o.slopeMedFail, o.slopeCIFail);
    end
    fprintf('\nMEAN SLOPE = slope of the replicate-averaged share of members; CI = cluster bootstrap over replicates (%d draws, seed %d)\n', cfg.nBoot, cfg.rngSeedMean);
    fprintf('%-24s %5s %29s %22s\n','setting','reps','mean slope /100 gen','IQR of rep slopes');
    for k = 1:numel(out)
        o = out(k);
        fprintf('%-24s %5d   %+8.5f [%+.5f,%+.5f]   [%+.4f,%+.4f]\n', o.name, o.nReps, o.slopeMean, o.slopeMeanCI, o.slopeIQR);
    end
    fprintf('\nMEAN SLOPE SPLIT BY FATE AT GENERATION %d (cooperation >= 0.95 then), and %d -> 500 fate transitions\n', cfg.fateGen, cfg.fateGen);
    fprintf('%-24s %36s %36s   %s\n','setting','successful at that generation','failed at that generation','S->S S->F F->S F->F');
    for k = 1:numel(out)
        o = out(k);
        fprintf('%-24s n=%3d %+8.5f [%+.5f,%+.5f]  n=%3d %+8.5f [%+.5f,%+.5f]  %4d %4d %4d %4d\n', o.name, ...
            o.nSuccFate, o.slopeMeanSuccFate, o.slopeMeanCISuccFate, o.nReps - o.nSuccFate, ...
            o.slopeMeanFailFate, o.slopeMeanCIFailFate, o.transitions);
    end
    fprintf('\nAUTOCORRELATION of within-replicate deviation of the share of members, steady window (gen >= 100)\n');
    fprintf('%-24s %5s %8s %11s %12s\n','setting','reps','spacing','rho(lag 1)','tau (gens)');
    for k = 1:numel(out)
        o = out(k);
        fprintf('%-24s %5d %8.1f %11.3f %12.1f\n', o.name, o.nReps, o.dg, o.rho1, o.tau);
    end
    fprintf('\nsources: %s\n', strjoin(unique({out.source}), ', '));

    outFile = fullfile(cfg.dataDir, 'share_of_members_drift.mat');
    metadata = cfg; save(outFile, 'out', 'metadata');
    fprintf('Saved: %s\n', outFile);
end


function [m, ci] = boot_median(x, nBoot)
% Median over replicates with a percentile bootstrap over replicates (global stream).
    x = x(~isnan(x));
    if isempty(x), m = NaN; ci = [NaN NaN]; return; end
    m = median(x);
    n = numel(x);
    bm = nan(nBoot,1);
    for b = 1:nBoot, bm(b) = median(x(randi(n, n, 1))); end
    ci = prctile(bm, [2.5 97.5]);
end


function [m, ci] = boot_mean(x, nBoot, stream)
%BOOT_MEAN  Mean over replicates with a cluster (whole-replicate) percentile bootstrap.
    x = x(~isnan(x));
    if isempty(x), m = NaN; ci = [NaN NaN]; return; end
    m = mean(x);
    n = numel(x);
    bm = nan(nBoot,1);
    for b = 1:nBoot, bm(b) = mean(x(randi(stream, n, n, 1))); end
    ci = prctile(bm, [2.5 97.5]);
end


function succ = fate_at_gen(dataDir, file, finalCoop, gen)
%FATE_AT_GEN  Success (cooperation >= 0.95) at generation GEN, per replicate.
%   The share-of-members files carry only the final cooperation, so the fate at
%   an earlier generation is read from the source run's cooperation history.
%   They cover the first numel(finalCoop) replicates of that run; the run stops
%   unless the final cooperation of both agrees to 1e-9.
    S = load(fullfile(dataDir, file), 'results');
    n = numel(finalCoop);
    C = cell2mat(cellfun(@(x) x.historyCooperation(:)', reshape(S.results(1:n), [], 1), ...
        'UniformOutput', false));
    assert(max(abs(C(:,end) - finalCoop)) < 1e-9, ...
        'fate_at_gen: %s does not match the traced replicates', file);
    succ = C(:, gen) >= 0.95;
end
