function out = share_of_members_by_generation(conds, cfg)
% SHARE_OF_MEMBERS_BY_GENERATION  The share of members at every stored generation of every
% replicate, testing each individual with the membership test (once per distinct wiring).
%
% Produces : SI S3.2 and Table S3, the drift of the share of members (input: the traces that
%            share_of_members_drift reads; run_saturating_separate_to_gen1000 and
%            measure_defector_encounter_frequency, for Table S9, also read them)
% Inputs   : conds   struct array, one element per setting: name, file (a run file or pattern
%                    under the publication folder, the newest match read; written by
%                    run_signaling_architectures), params (sensingMode, hillCoefficient,
%                    discreteOutput)
%            cfg     dataDir (the publication folder), nReps (Inf: all), outFile ('' saves
%                    nothing)
% Outputs  : cfg.outFile (publication/member_share_by_generation_architectures.mat, _standard100.mat or
%              _satsep100.mat in the paper's runs): out, per setting the file read, the stored
%              generations, the share of members (replicate x stored generation), the members'
%              and non-members' mean offspring numbers, final cooperation, success and the mean
%              share of members from generation 100 on; prints one row per setting
% Run      : three calls: the four architectures with nReps 50 (_architectures), saturating
%              contrast with 100 (_standard100), saturating separate with 100 (_satsep100); e.g.
%              share_of_members_by_generation(struct('name', 'sat_sep', 'file',
%              'saturating_separate_100reps.mat', 'params', struct('sensingMode', 'separate',
%              'hillCoefficient', 1, 'discreteOutput', false)), struct('nReps', 100, 'outFile',
%              grn_data_file('publication/member_share_by_generation_satsep100.mat', 'new')))
%              (long, not timed here: 3 replicates of one setting take about 25 s)
% Seed     : none (no random draws; membership_test seeds each wiring itself)
%
% Membership uses the test's own settings (600 steps from the fixed start).

    if nargin < 2 || isempty(cfg), cfg = struct(); end
    if ~isfield(cfg,'dataDir') || isempty(cfg.dataDir), cfg.dataDir = fullfile(grn_data_file('publication')); end
    if ~isfield(cfg,'nReps'),   cfg.nReps = Inf;   end
    if ~isfield(cfg,'outFile'), cfg.outFile = '';  end
    base = struct('michaelisConstant',1.0,'degradationRate',0.1,'timeStep',1.0,'developmentSteps',600,'backend','simulate_pair');

    nC = numel(conds); out = struct('name',{conds.name});
    for c = 1:nC
        f = dir(fullfile(cfg.dataDir, conds(c).file)); if isempty(f), error('No file: %s', conds(c).file); end
        [~,ix] = max([f.datenum]); d = load(fullfile(f(ix).folder, f(ix).name), 'results');
        useReps = min(cfg.nReps, numel(d.results));
        params = base; pf = fieldnames(conds(c).params); for i = 1:numel(pf), params.(pf{i}) = conds(c).params.(pf{i}); end
        rep1 = d.results{1}; nK = numel(rep1.historyGenerations); N = size(rep1.historyM,3);
        share = nan(useReps, nK); wMem = nan(useReps, nK); wNon = nan(useReps, nK); nMem = zeros(useReps, nK);
        fc = nan(useReps,1); sg = double(rep1.historyGenerations(:))';
        t0 = tic; fprintf('\n=== %s: %s (%d reps, N = %d, %d snapshots) ===\n', conds(c).name, f(ix).name, useReps, N, nK);
        % Each stored generation of each replicate: the membership test once per distinct
        % wiring, the share of members, and the offspring numbers of members and non-members.
        % Stored generation k holds the population at the start of generation g = sg(k); its
        % offspring are counted from historyParents(:, g), the parent map from g to g + 1.
        for r = 1:useReps
            rep = d.results{r}; fc(r) = rep.historyCooperation(end); P = rep.historyParents; nG = size(P,2);
            for k = 1:nK
                g = sg(k);
                M = double(rep.historyM(:,:,:,k)); B = double(rep.historyB(:,:,:,k));
                Wall = reshape(M.*B, 16, N)'; [Wu,~,ic] = unique(Wall,'rows');
                mem = false(size(Wu,1),1);
                for u = 1:size(Wu,1), mem(u) = membership_test(reshape(Wu(u,:),4,4), params).is_css; end
                isMem = mem(ic);
                share(r,k) = mean(isMem); nMem(r,k) = sum(isMem);
                if g >= 1 && g <= nG
                    off = accumarray(double(P(:,g)), 1, [N 1]);   % offspring of each slot at generation g
                    if any(isMem),  wMem(r,k) = mean(off(isMem));  end
                    if any(~isMem), wNon(r,k) = mean(off(~isMem)); end
                end
            end
            if mod(r,10) == 0, fprintf('  rep %d/%d (%.0fs, ETA %.0fs)\n', r, useReps, toc(t0), toc(t0)*(useReps-r)/r); end
        end
        steady = sg >= 100;
        out(c).file = f(ix).name; out(c).N = N; out(c).nReps = useReps; out(c).gens = sg;
        out(c).occ = share; out(c).wMem = wMem; out(c).wNon = wNon; out(c).nMem = nMem; out(c).finalCoop = fc;
        out(c).success = mean(fc >= 0.95);
        out(c).meanShare100_500 = mean(share(:,steady),'all');
        fprintf('%s: success %.0f%% | mean share of members (100-500) %.3f (%.1f min)\n', ...
            conds(c).name, 100*out(c).success, out(c).meanShare100_500, toc(t0)/60);
    end
    fprintf('\n%-14s %5s %8s %8s\n','setting','N','success','share');
    for c = 1:nC
        fprintf('%-14s %5d %7.0f%% %8.3f\n', out(c).name, out(c).N, 100*out(c).success, out(c).meanShare100_500);
    end
    if ~isempty(cfg.outFile), metadata = cfg; save(cfg.outFile, 'out', 'metadata'); fprintf('Saved: %s\n', cfg.outFile); end
end
