function repair_test_without_defector()
% REPAIR_TEST_WITHOUT_DEFECTOR  Does the repair's effect need the added defector? Without it,
% prints sham minus repair, observed minus sham, observed minus repair, and the part specific to
% restoring retaliation, (sham minus repair) / (observed minus repair); then the same four with the
% defector, as a check.
%
% Produces : SI S4.3, "The repair's effect does not need the added defector"; Table S10
% Inputs   : publication/repair_matched/repair_matched_build_targets.mat and repair_matched_main.mat
%              (run_repair_test, modes 'build' and 'main')
%            publication/repair_test_observed/repair_test_observed.mat   (run_repair_test_observed)
%            (under the data root, or under GRN_RESULTS when that variable is set)
% Outputs  : publication/repair_matched/repair_no_later_defector.txt   (the printed lines; the run stops
%              if the file exists)
% Run      : repair_test_without_defector   (about 1 min)
% Seed     : rng(20261013, 'twister'), then randi(65, 65, 10000): the same resamples as
%            repair_test_summary('main') and figS1_draw
%
% An arm's collapse share is the mean over its runs; 95% percentile bootstrap intervals over the
% 65 populations (checked).

root = getenv('GRN_RESULTS');
if isempty(root), root = fileparts(grn_data_file('publication')); end
pub = fullfile(root, 'publication');
dir0 = fullfile(pub, 'repair_matched');
outfile = fullfile(dir0, 'repair_no_later_defector.txt');
if exist(outfile, 'file'), error('repair_test_without_defector:exists', 'Refusing to overwrite %s', outfile); end
fid = fopen(outfile, 'w');
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

Bd = load(fullfile(dir0, 'repair_matched_build_targets.mat'), 'out'); Bd = Bd.out;
R = load(fullfile(dir0, 'repair_matched_main.mat'), 'out'); R = R.out;
I = load(fullfile(pub, 'repair_test_observed', 'repair_test_observed.mat'), 'out');
T = Bd.targets;
T = T([T.inRun] & ~[T.excluded]);
assert(isequal(sort({T.id}), sort(R.events)), 'run events differ from the build''s included events');
n = numel(T);
assert(n == 65, 'expected 65 included events, found %d', n);

X.pA = share(T, I.out.runs, 'A');   X.pA0 = share(T, I.out.runs, 'A0');
X.pSp = share(T, R.runs, 'Sp');     X.pS0p = share(T, R.runs, 'S0p');
X.pBp = share(T, R.runs, 'Bp');     X.pB0p = share(T, R.runs, 'B0p');

rng(20261013, 'twister');
bi = randi(n, n, 10000);

say(fid, 'Repair test, matched start: %d populations, 10 runs per arm, 10,000 bootstrap resamples over populations\n', n);
say(fid, 'means: defector added: observed %.3f, sham %.3f, repair %.3f | no defector added: observed %.3f, sham %.3f, repair %.3f\n\n', ...
    mean(X.pA), mean(X.pSp), mean(X.pBp), mean(X.pA0), mean(X.pS0p), mean(X.pB0p));
diffline(fid, 'no defector added: sham - repair    ', X.pS0p - X.pB0p, bi);
diffline(fid, 'no defector added: observed - sham  ', X.pA0 - X.pS0p, bi);
diffline(fid, 'no defector added: observed - repair', X.pA0 - X.pB0p, bi);
ratioline(fid, 'no defector added: specific share   ', X.pS0p - X.pB0p, X.pA0 - X.pB0p, bi);
say(fid, '\nCheck, with the defector added\n');
diffline(fid, 'defector added: sham - repair       ', X.pSp - X.pBp, bi);
diffline(fid, 'defector added: observed - sham     ', X.pA - X.pSp, bi);
diffline(fid, 'defector added: observed - repair   ', X.pA - X.pBp, bi);
ratioline(fid, 'defector added: specific share      ', X.pSp - X.pBp, X.pA - X.pBp, bi);
end

function p = share(T, pool, arm)
% Collapse share per event: mean of the collapse flags of the one job for that event and arm.
p = nan(numel(T), 1);
for i = 1:numel(T)
    r = pool([pool.event] == T(i).num & strcmp({pool.arm}, arm));
    assert(numel(r) == 1, 'event %s arm %s: expected one job', T(i).id, arm);
    p(i) = mean(r.r.collapse);
end
end

function diffline(fid, name, x, bi)
ci = quantile(mean(x(bi), 1), [0.025 0.975]);
say(fid, '  %s %+.3f [%+.3f, %+.3f]\n', name, mean(x), ci);
end

function ratioline(fid, name, num, den, bi)
ci = quantile(mean(num(bi), 1) ./ mean(den(bi), 1), [0.025 0.975]);
say(fid, '  %s %.3f [%.3f, %.3f]\n', name, mean(num) / mean(den), ci);
end

function say(fid, varargin)
fprintf(1, varargin{:}); fprintf(fid, varargin{:});
end
