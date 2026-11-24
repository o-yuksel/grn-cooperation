function tracking_load_exponents()
% TRACKING_LOAD_EXPONENTS  With the target held still and alpha lowered to 0.70, the median
% per-replicate ratio of the mutation load to the control's in designs 1 and 2, and the exponents
% of the mutation load and of the per-generation mutation cost on the factor 0.7.
%
% Produces : SI S5.1, "Design 1 sits nearer the small-effect regime" (printed beside the SI values)
% Inputs   : publication/tracking_cut_after_150/design1/ and design2/main.mat (written by
%              tracking_cut_after_150): out.analysis.Y, replicates x arms x outcomes; arm 1 the
%              control (`'CONTROL'`), arm 3 alpha lowered (`'ALPHA-DOWN'`); the outcomes mutation
%              load (`mutLoad`) and per-generation mutation cost (`mutCost`), found by name
% Outputs  : publication/tracking_cut_after_150/load_exponents.txt (the printed
%              lines; the run stops if the file exists); read and written under
%              GRN_RESULTS when that variable is set, else under the data root
% Run      : tracking_load_exponents()   (seconds)
% Seed     : none (no random draws)
%
% The exponent of a quantity on the factor 0.7 is ln(median ratio) / ln(0.7), the ratio taken per
% replicate, alpha lowered over control.

root = getenv('GRN_RESULTS');
if isempty(root), root = fileparts(grn_data_file('publication')); end
base = fullfile(root, 'publication', 'tracking_cut_after_150');
outfile = fullfile(base, 'load_exponents.txt');
if exist(outfile, 'file'), error('tracking_load_exponents:exists', 'Refusing to overwrite %s', outfile); end
fid = fopen(outfile, 'w');
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

say(fid, 'Target held still: alpha lowered to 0.7 against control, median per-replicate ratio\n');
r = zeros(1, 2); eL = r; eC = r;
for d = 1:2
    L = load(fullfile(base, sprintf('design%d', d), 'main.mat'), 'out'); A = L.out.analysis;
    iL = find(strcmp(A.names, 'mutLoad')); iC = find(strcmp(A.names, 'mutCost'));
    assert(strcmp(L.out.cfg.arms{1}, 'CONTROL') && strcmp(L.out.cfg.arms{3}, 'ALPHA-DOWN'), 'arm order');
    rl = A.Y(:, 3, iL) ./ A.Y(:, 1, iL); rc = A.Y(:, 3, iC) ./ A.Y(:, 1, iC);
    r(d) = median(rl, 'omitnan'); eL(d) = log(r(d)) / log(0.7);
    eC(d) = log(median(rc, 'omitnan')) / log(0.7);
    say(fid, '  design %d (%d replicates): load ratio %.4f, load exponent %.3f, cost exponent %.3f\n', ...
        d, size(A.Y, 1), r(d), eL(d), eC(d));
end
say(fid, '\n                       design 1 / design 2\n');
say(fid, '  median ratio          %.3f / %.3f\n', r);
say(fid, '  load exponents        %.2f / %.2f\n', eL);
say(fid, '  cost exponents        %.2f / %.2f\n', eC);
end

function say(fid, varargin)
fprintf(1, varargin{:}); fprintf(fid, varargin{:});
end

