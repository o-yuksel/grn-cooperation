function tableS6_per_mutation_shares(pubDir)
% TABLES6_PER_MUTATION_SHARES  From the per-mutation census: the share of mutations to a member
% that end membership, end cooperation, or end retaliation only, by architecture and pooled.
%
% Produces : Table S6 (printed); main text, "Conditional cooperation stays within a mutation's
%            reach of the threshold" ("more than a third end membership"); SI S3.2, "Most
%            mutations that end membership keep cooperation" (the pooled shares; 65.3% is their
%            ratio)
% Inputs   : publication/rho_by_architecture.mat   (written by per_mutation_census); pubDir
%              defaults to the publication folder under the data root
% Outputs  : publication/census_overall_shares/census_overall_shares.txt   (the lines this run
%              prints: per architecture and pooled, the shares with 95% intervals; the
%              unweighted draw counts and shares; a check of the pooled row against fixed values)
% Run      : tableS6_per_mutation_shares   (about 20 s)
% Seed     : rng(20261029), generator not named (MATLAB's default, Mersenne twister), for the
%            bootstrap
%
% Members kept: those the census marks ok (cooperation margin at least 1e-6, ||x*|| in [0.5, 14]).
% Shares are weighted by the individuals carrying each wiring, with intervals from 2,000
% resamples of members within architecture; the draw counts printed after them are unweighted.
    if nargin < 1 || isempty(pubDir)
        cfg = figure_config();
        pubDir = fullfile(cfg.dataDir, 'publication');
    end
    seed = 20261029; nBoot = 2000;
    matFile = fullfile(pubDir, 'rho_by_architecture.mat');
    d = dir(matFile);
    outDir = fullfile(pubDir, 'census_overall_shares');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    outFile = fullfile(outDir, 'census_overall_shares.txt');

    %% Report
    report = @() shares_report(matFile, d, seed, nBoot); %#ok<NASGU> called inside evalc
    rep = evalc('report();');
    fprintf('%s', rep);

    %% Save
    fid = fopen(outFile, 'w'); fprintf(fid, '%s', rep); fclose(fid);
end


function shares_report(matFile, d, seed, nBoot)
% The table's rows, the raw counts and the check; prints the report's lines.
    fprintf('# tableS6_per_mutation_shares.m | input %s (file timestamp %s) | seed %d | %s\n', ...
        matFile, datestr(d.datenum, 'yyyy-mm-dd HH:MM:SS'), seed, datestr(now, 'yyyy-mm-dd'));
    L = load(matFile, 'W_all', 'modeNames', 'metadata');
    K = L.metadata.cfg.K;
    nM = numel(L.modeNames);
    % Per kept member: its weight (individuals carrying the wiring) and how many of its K draws
    % end membership (stored field 'ejected'), end cooperation ('lostReadout'), or end
    % membership with cooperation kept (retaliation only). Each member's K draws are stored as
    % one K-by-1 column; stacking them, reshaping to K rows and transposing puts members in rows
    % and their draws in columns.
    S = cell(1, nM);   % per member: weight, and draw counts of each event
    for m = 1:nM
        R = L.W_all{m};
        Rk = R(logical([R.ok])');
        S{m}.w   = double([Rk.count])';
        S{m}.endsMem  = sum(reshape(vertcat(Rk.ejected),     K, [])', 2);
        S{m}.lo  = sum(reshape(vertcat(Rk.lostReadout), K, [])', 2);
        S{m}.c2  = sum(reshape(vertcat(Rk.ejected) & ~vertcat(Rk.lostReadout), K, [])', 2);
        S{m}.n   = numel(S{m}.w);
    end
    names = L.modeNames; clear L;
    fprintf('K = %d draws per member; weights = individuals carrying each wiring; members kept: census ok\n', K);

    % Weighted totals of the three events and of draws over the members idx; the bootstrap
    % resamples members (all K draws and the weight together) within each architecture, and the
    % pooled row sums the architectures' resampled totals, so it is stratified by architecture
    rng(seed);
    est = @(s, idx) [sum(s.w(idx).*s.endsMem(idx)), sum(s.w(idx).*s.lo(idx)), ...
                     sum(s.w(idx).*s.c2(idx)), sum(s.w(idx))*K];
    fprintf('%-14s %7s %8s | %-24s | %-24s | %-24s\n', 'architecture', 'members', 'draws', ...
        'P(ends membership) [95% CI]', 'P(ends cooperation) [95% CI]', 'P(ends only retaliation) [95% CI]');
    tot = zeros(1, 4); cnt = zeros(1, 3); nMem = 0;
    boot = zeros(nBoot, 4, nM);
    for m = 1:nM
        s = S{m};
        t = est(s, 1:s.n);
        for b = 1:nBoot
            boot(b, :, m) = est(s, randi(s.n, s.n, 1));
        end
        tot = tot + t; nMem = nMem + s.n;
        cnt = cnt + [sum(s.endsMem), sum(s.lo), sum(s.c2)];
        print_row(names{m}, s.n, s.n*K, t, boot(:, :, m));
    end
    print_row('pooled', nMem, nMem*K, tot, sum(boot, 3));
    fprintf('\nraw counts pooled (unweighted draws): ends membership k = %d, ends cooperation k = %d, ends only retaliation k = %d, of %d draws\n', ...
        cnt(1), cnt(2), cnt(3), nMem*K);
    fprintf('raw (unweighted) shares: %.4f / %.4f / %.4f\n', cnt/(nMem*K));

    % Check of the pooled shares and counts against fixed values
    p = tot(1:3)/tot(4);
    ok = all(abs(round(p, 4) - [0.3597 0.1248 0.2349]) < 1e-12) && nMem*K == 457320 && ...
         nMem == 15244 && cnt(1) == 154371 && cnt(3) == 99861;
    fprintf('check (0.3597 / 0.1248 / 0.2349, 457320 draws, 15244 members, 154371, 99861): %s\n', ...
        ternary(ok, 'PASS', 'FAIL'));
end

function print_row(name, nMem, nDraw, t, bt)
% One row: the weighted shares and their 95% percentile intervals over the bootstrap totals bt
    p = t(1:3)/t(4);
    q = bt(:, 1:3) ./ bt(:, 4);
    ci = prctile(q, [2.5 97.5]);
    fprintf('%-14s %7d %8d | %.4f [%.4f, %.4f] | %.4f [%.4f, %.4f] | %.4f [%.4f, %.4f]\n', ...
        name, nMem, nDraw, p(1), ci(1,1), ci(2,1), p(2), ci(1,2), ci(2,2), p(3), ci(1,3), ci(2,3));
end
