function tracking_load_rank_correlation()
% TRACKING_LOAD_RANK_CORRELATION  With the target held still, does the cut in the mutation load
% from lowering responsiveness alpha follow the shift one mutation causes? The Spearman correlation of the two
% over the eleven design-and-size points of Table S15.
%
% Produces : SI S5.1, "The cut follows the shift one mutation causes, not the design" (the
%            Spearman correlations over the eleven points, on the median shift and on the mean)
% Inputs   : under publication/tracking_load_regime/ (or GRN_RESULTS if set):
%            shift_by_size.txt     median and mean shift of the first seven points (written by
%                              run_tracking_cut_vs_shift, mode 'x')
%            cut_vs_shift_summary.txt   their cuts, "CUT VS SHIFT FINAL" table (run_tracking_cut_vs_shift, 'final')
%            shift_by_size_design1.txt     median and mean shift of four more design-1 sizes
%                              (run_tracking_cut_vs_shift, 'x3')
%            cut_vs_size_summary.txt   their cuts, "CUT VS SIZE FINAL" (run_tracking_cut_vs_mutation_size, 'final')
% Outputs  : publication/tracking_load_regime/load_rank_correlation.txt (the printed lines; the run
%              stops if the file exists)
% Run      : tracking_load_rank_correlation()   (seconds; corr needs the Statistics and Machine
%              Learning Toolbox)
% Seed     : none (no random draws)
%
% The first seven points are design 1 at sd 0.020, 0.060 and 0.200 and the four design-2 rows
% (sd 0.020, 0.0063, 0.0018 and 0.0003); the four more are design 1 at sd 0.006, 0.010, 0.040 and
% 0.100. The shift is the median over the 30 generation-350 best genotypes and, second, their mean;
% the cut is the alpha-lowered (`ALPHA`) cut as a percentage of the control's median. The same
% correlation over the first seven points is printed to compare with the Spearman lines of
% cut_vs_shift_summary.txt.

root = getenv('GRN_RESULTS');
if isempty(root), root = fileparts(grn_data_file('publication')); end
base = fullfile(root, 'publication', 'tracking_load_regime');
outfile = fullfile(base, 'load_rank_correlation.txt');
if exist(outfile, 'file'), error('tracking_load_rank_correlation:exists', 'Refusing to overwrite %s', outfile); end

x2 = fileread(fullfile(base, 'shift_by_size.txt'));
s2 = fileread(fullfile(base, 'cut_vs_shift_summary.txt'));
x3 = fileread(fullfile(base, 'shift_by_size_design1.txt'));
s3 = fileread(fullfile(base, 'cut_vs_size_summary.txt'));
s2 = s2(strfind(s2, 'CUT VS SHIFT FINAL'):end);
s3 = s3(strfind(s3, 'CUT VS SIZE FINAL'):end);

% Cut-vs-shift contexts: design 1 at sd 0.020, 0.060, 0.200; the four design-2 rows.
lab2 = {'d1 sd0.020', 'd1 sd0.060', 'd1 sd0.200', 'd2 sd0.020', 'd2 sd0.0063', 'd2 sd0.001821', 'd2 A-mean'};
% Cut-vs-size runs: design 1 at sd 0.006, 0.010, 0.040, 0.100.
lab3 = {'0.3x', '0.5x', '2x', '5x'};
run3 = {'E', 'F', 'G', 'H'};

n = numel(lab2) + numel(lab3);
P = struct('label', cell(n, 1), 'sd', [], 'med', [], 'mean', [], 'cut', []);
num = '([-+]?[\d.]+)';
% Cut-vs-shift rows: in shift_by_size.txt the label is followed by the design's arch, then sd, p, median and
% mean shift; in the cut-vs-shift table by the same, then the ALPHA cut. The run stops if the two files
% disagree on sd, median or mean.
for i = 1:numel(lab2)
    L = regexptranslate('escape', lab2{i});
    t = regexp(x2, ['(?m)^' L '\s+\d\s+' num '\s+' num '\s+' num '\s+' num], 'tokens', 'once');
    u = regexp(s2, ['(?m)^' L '\s+\d\s+' num '\s+' num '\s+' num '\s+' num '\s+' num '\s+\['], 'tokens', 'once');
    assert(~isempty(t) && ~isempty(u), 'tracking_load_rank_correlation:parse', 'cut-vs-shift row not found: %s', lab2{i});
    t = str2double(t); u = str2double(u);
    assert(t(1) == u(1) && t(3) == u(3) && t(4) == u(4), 'tracking_load_rank_correlation:mismatch', ...
        'shift_by_size.txt and cut_vs_shift_summary.txt disagree on %s', lab2{i});
    P(i) = struct('label', lab2{i}, 'sd', t(1), 'med', t(3), 'mean', t(4), 'cut', u(5));
end
% Cut-vs-size rows: in shift_by_size_design1.txt the size label is followed by sd, median and mean shift; in the
% cut-vs-size table the label and run letter by the size multiple, then the ALPHA cut
for j = 1:numel(lab3)
    L = regexptranslate('escape', lab3{j});
    t = regexp(x3, ['(?m)^' L '\s+' num '\s+' num '\s+' num], 'tokens', 'once');
    u = regexp(s3, ['(?m)^' L ' \(' run3{j} '\)\s+' num '\s+' num '\s+\['], 'tokens', 'once');
    assert(~isempty(t) && ~isempty(u), 'tracking_load_rank_correlation:parse', 'cut-vs-size row not found: %s', lab3{j});
    t = str2double(t); u = str2double(u);
    P(numel(lab2) + j) = struct('label', ['d1 ' lab3{j} ' (run ' run3{j} ')'], 'sd', t(1), ...
        'med', t(2), 'mean', t(3), 'cut', u(2));
end

fid = fopen(outfile, 'w');
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
say(fid, 'SI S5.1: Spearman correlation of the ALPHA-DOWN cut against the shift per mutation\n');
say(fid, 'inputs: %s/{shift_by_size,cut_vs_shift_summary,shift_by_size_design1,cut_vs_size_summary}.txt\n\n', base);
say(fid, '%-24s %10s %10s %10s %8s\n', 'point', 'sd', 'median x', 'mean x', 'cut %');
for k = 1:n
    say(fid, '%-24s %10.6f %10.5f %10.5f %8.1f\n', P(k).label, P(k).sd, P(k).med, P(k).mean, P(k).cut);
end

y = [P.cut]'; xm = [P.med]'; xa = [P.mean]';
rM = corr(xm, y, 'type', 'Spearman');
rA = corr(xa, y, 'type', 'Spearman');
say(fid, '\nEleven points:\n');
say(fid, '  Spearman(median shift, cut) = %+.4f  (n = %d)\n', rM, n);
say(fid, '  Spearman(mean shift, cut)   = %+.4f  (n = %d)\n', rA, n);

i7 = 1:numel(lab2);
say(fid, '\nCheck, the seven cut-vs-shift points (as cut_vs_shift_summary.txt prints them):\n');
say(fid, '  Spearman(median shift, cut) = %+.4f\n', corr(xm(i7), y(i7), 'type', 'Spearman'));
say(fid, '  Spearman(mean shift, cut)   = %+.4f\n', corr(xa(i7), y(i7), 'type', 'Spearman'));
say(fid, '\nTies: median shift %d, mean shift %d, cut %d (a tie would make the rank correlation depend on rounding)\n', ...
    n - numel(unique(xm)), n - numel(unique(xa)), n - numel(unique(y)));
end


function say(fid, varargin)
fprintf(varargin{:});
fprintf(fid, varargin{:});
end
