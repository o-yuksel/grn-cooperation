function tableS5_parameter_sweeps()
% TABLES5_PARAMETER_SWEEPS  How does success change when one parameter moves away from the standard
% setting? Recomputes Table S5: success, its interval and mean cooperation under five one-parameter
% sweeps (temptation T, expression noise, population size N, selection power zeta, mutation size
% sigma_mut), and the range of graded-separate success over the same sigma_mut values.
%
% Produces : Table S5 and its caption's graded-separate range (SI S2.6)
% Inputs   : publication/sweep_T_<T>_100reps_*.mat   (T = 2.5, 3.0: run_low_temptation; 3.5 to
%              5.8: run_signaling_architectures; 5.8 is the run at T = 5.75, named to one decimal)
%            publication/sweep_noise_<sigma>_100reps_*.mat and
%              sweep_selection_s<zeta>_100reps_*.mat   (run_selection_noise_sweeps)
%            publication/sweep_N_<N>_50reps_*.mat   (run_signaling_architectures)
%            publication/sweep_mutation_size_saturating_contrast_<sigma_mut>_100reps_*.mat and
%              sweep_mutation_size_graded_separate_<sigma_mut>_100reps_*.mat   (run_mutation_size_sweep)
%            publication/saturating_contrast_100reps.mat   (run_signaling_architectures; the
%              standard setting)
% Outputs  : publication/table_s3_sweeps.txt   (despite its name, this file holds Table S5) and the
%              same lines on the console
% Run      : tableS5_parameter_sweeps   (not timed)
% Seed     : none (no random draws)
%
% Success: the share of replicates whose last-generation cooperation is at least 0.95. Where
% several files match one pattern, their success counts must agree or the run stops.

pub = grn_data_file('publication');
outfile = grn_data_file('publication/table_s3_sweeps.txt', 'new');
fid = fopen(outfile, 'w');
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

sat = 'saturating_contrast_100reps.mat';
% {block, value, file pattern}
R = {
 'T','2.5','sweep_T_2.5_100reps_*.mat'
 'T','3.0','sweep_T_3.0_100reps_*.mat'
 'T','3.5','sweep_T_3.5_100reps_*.mat'
 'T','4.0','sweep_T_4.0_100reps_*.mat'
 'T','4.5','sweep_T_4.5_100reps_*.mat'
 'T','5.0',sat
 'T','5.5','sweep_T_5.5_100reps_*.mat'
 'T','5.8','sweep_T_5.8_100reps_*.mat'
 'noise','0.00','sweep_noise_0.00_100reps_*.mat'
 'noise','0.01','sweep_noise_0.01_100reps_*.mat'
 'noise','0.05','sweep_noise_0.05_100reps_*.mat'
 'noise','0.10','sweep_noise_0.10_100reps_*.mat'
 'N','50','sweep_N_50_50reps_*.mat'
 'N','100','sweep_N_100_50reps_*.mat'
 'N','200','sweep_N_200_50reps_*.mat'
 'N','500','sweep_N_500_50reps_*.mat'
 'N','1000',sat
 'zeta','0.5','sweep_selection_s0.5_100reps_*.mat'
 'zeta','1.0','sweep_selection_s1.0_100reps_*.mat'
 'zeta','2.0','sweep_selection_s2.0_100reps_*.mat'
 'zeta','4.0','sweep_selection_s4.0_100reps_*.mat'
 'zeta','8.0','sweep_selection_s8.0_100reps_*.mat'
 'sigma_mut','0.01','sweep_mutation_size_saturating_contrast_0.010_100reps_*.mat'
 'sigma_mut','0.02','sweep_mutation_size_saturating_contrast_0.020_100reps_*.mat'
 'sigma_mut','0.05','sweep_mutation_size_saturating_contrast_0.050_100reps_*.mat'
 'sigma_mut','0.10','sweep_mutation_size_saturating_contrast_0.100_100reps_*.mat'
 'sigma_mut','0.20','sweep_mutation_size_saturating_contrast_0.200_100reps_*.mat'
 'sigma_mut','0.50','sweep_mutation_size_saturating_contrast_0.500_100reps_*.mat'
 'sigma_mut','1.00','sweep_mutation_size_saturating_contrast_1.000_100reps_*.mat'
 };

say(fid, 'Table S5 (tab:tsensitivity). Per row: success %%, 95%% CI, mean coop %%\n');
blk = '';
for k = 1:size(R,1)
    if ~strcmp(blk, R{k,1}), blk = R{k,1}; say(fid, '\n[%s]\n', blk); end
    [cnt, n, m, copies, srStored] = read_row(pub, R{k,3});
    [~, pci] = binofit(cnt, n, 0.05); pci = pci * 100;
    s = sprintf('%d', round(100*cnt/n));
    c = sprintf('[%.1f, %.1f]', pci(1), pci(2));
    mm = sprintf('%.1f', m);
    say(fid, '  %-5s %3d/%-3d  %s %s %s  (files=%d, stored successRate=%s)\n', ...
        R{k,2}, cnt, n, s, c, mm, copies, srStored);
end

% The caption's range of graded-separate success over sigma_mut = 0.01 to 1.00
say(fid, '\nGraded-separate success over sigma_mut 0.01 to 1.00\n');
g = {'0.010','0.020','0.050','0.100','0.200','0.500','1.000'}; pct = zeros(1, 7);
for k = 1:7
    [cnt, n] = read_row(pub, sprintf('sweep_mutation_size_graded_separate_%s_100reps_*.mat', g{k}));
    pct(k) = 100*cnt/n;
end
say(fid, '  graded separate success %%: %s  -> range %d-%d%%\n', mat2str(pct), min(pct), max(pct));
end

% Success count, replicate count and mean cooperation (percent) of the files matching pattern,
% with the stored successRate where a file has one (printed beside the row as a check)
function [cnt, n, m, copies, srStored] = read_row(pub, pattern)
d = dir(fullfile(pub, pattern));
if isempty(d), error('No file for %s', pattern); end
copies = numel(d); cnts = zeros(1, copies); ms = zeros(1, copies);
srStored = 'none';
for i = 1:copies
    S = load(fullfile(pub, d(i).name));
    if isfield(S, 'coopValues')
        v = double(S.coopValues(:));
    else
        v = cellfun(@(r) double(r.historyCooperation(end)), S.results(:));
    end
    cnts(i) = sum(v >= 0.95); ms(i) = mean(v) * 100; n = numel(v);
    if isfield(S, 'successRate'), srStored = sprintf('%g', S.successRate); end
end
if any(cnts ~= cnts(1)), error('Copies of %s disagree: %s', pattern, mat2str(cnts)); end
cnt = cnts(1); m = ms(1);
end

function say(fid, varargin)
fprintf(1, varargin{:}); fprintf(fid, varargin{:});
end
