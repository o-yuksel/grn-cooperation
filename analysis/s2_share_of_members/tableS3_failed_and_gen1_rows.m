function tableS3_failed_and_gen1_rows()
% TABLES3_FAILED_AND_GEN1_ROWS  Two parts of Table S3 for the four signaling architectures: the
% failed populations and the generation-1 share of members.
%
% Produces : Table S3, "Failed populations" (number, exceeded 0.90, median peak and final) and
%            "Share of members, Generation 1"
% Inputs   : publication/<architecture>_100reps.mat   (run_signaling_architectures; cooperation by
%              generation)
%            publication/scope_checks/scope_checks_report.txt   (scope_checks_summary; the
%              generation-1 member counts)
% Outputs  : publication/table_s4_failed_and_gen1.txt   (despite its name, it holds these Table S3
%              rows: the lines this run prints)
% Run      : tableS3_failed_and_gen1_rows   (about 30 s)
% Seed     : none (no random draws)
%
% Failed: final mean cooperation below 0.95; "exceeded 0.90" counts those whose cooperation rose
% above 0.90 at some generation, the SI's rule ("after exceeding 0.90"). Generation 1: members
% among 10,000 random wirings, with Clopper-Pearson 95% intervals; each count is the first integer
% after the architecture's name in the report's 'GENERATION-1 SHARES' block.

outfile = grn_data_file('publication/table_s4_failed_and_gen1.txt', 'new');
fid = fopen(outfile, 'w');
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

arch = {'graded separate', 'graded contrast', 'saturating separate', 'saturating contrast'};
files = {'graded_separate', 'graded_contrast', 'saturating_separate', 'saturating_contrast'};
say(fid, 'Table S3 (tab:byarch): failed populations (final cooperation < 0.95)\n');
for a = 1:4
    S = load(grn_data_file(['publication/' files{a} '_100reps.mat']), 'results');
    n = numel(S.results); fin = zeros(n, 1); pk = zeros(n, 1);
    for r = 1:n
        c = double(S.results{r}.historyCooperation);
        fin(r) = c(end); pk(r) = max(c);
    end
    f = fin < 0.95;
    nExc = sum(pk(f) > 0.90);
    mPk = median(pk(f)); mFin = median(fin(f));
    say(fid, '  %-20s failed %3d of %d | exceeded 0.90 %2d | median peak %.3f | median final %.3f\n', ...
        arch{a}, sum(f), n, nExc, mPk, mFin);
end

say(fid, '\nGeneration-1 share of members, 10,000 random wirings per architecture (Table S3), Clopper-Pearson 95%%\n');
% Member counts: the first integer after the architecture's name in the GENERATION-1 SHARES block
txt = fileread(grn_data_file('publication/scope_checks/scope_checks_report.txt'));
i1 = strfind(txt, 'GENERATION-1 SHARES'); txt = txt(i1(1):end);
for a = 1:4
    tok = regexp(txt, [arch{a} '\s+(\d+)\s+\d+\s+'], 'tokens', 'once');
    k = str2double(tok{1});
    [~, ci] = binofit(k, 10000, 0.05);
    s = sprintf('%.2f [%.2f, %.2f]', 100*k/10000, 100*ci(1), 100*ci(2));
    say(fid, '  %-20s %5d of 10000 -> %s\n', arch{a}, k, s);
end
end

function say(fid, varargin)
fprintf(1, varargin{:}); fprintf(fid, varargin{:});
end
