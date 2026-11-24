function D = fig2_data(cfg)
% FIG2_DATA  The data of Fig. 2: per signaling architecture, mean cooperation and success over
% generations 1 to 500, and each population's cooperation at generation 500.
%
% Implements : the data of Fig. 2, read by fig2_draw (called by fig2_signaling_architecture)
% Inputs     : cfg   figure_config struct (default figure_config()); cfg.dataDir, or its
%                    publication folder, holds publication/<architecture>_100reps.mat (written by
%                    run_signaling_architectures)
% Outputs    : D   struct, one row or cell per architecture: gens; meanCoopOverTime (percent, mean
%                  over replicates) with coopCiLo, coopCiHi (mean -/+ 1.96 standard errors, cut to
%                  [0, 100]); successOverTime (percent of replicates at or above 0.95);
%                  finalCoopByMode (generation-500 cooperation of each replicate); mode_labels,
%                  short_labels, mode_colors (figure_utils.architecture_palette); nModes
% Seed       : none (no random draws)

    if nargin < 1, cfg = figure_config(); end

    dataDir = cfg.dataDir;
    if isfolder(fullfile(dataDir, 'publication'))   % the data root, or already its publication folder
        dataDir = fullfile(dataDir, 'publication');
    end

    modes = {'graded_separate', 'graded_contrast', ...
             'saturating_separate', 'saturating_contrast'};
    palette = figure_utils.architecture_palette();

    nModes = length(modes);
    nGensB = 500;
    meanCoopOverTime = zeros(nModes, nGensB);
    meanCoopSEM      = zeros(nModes, nGensB);
    successOverTime  = zeros(nModes, nGensB);
    finalCoopByMode  = cell(nModes, 1);
    for m = 1:nModes
        % The pattern holds no wildcard, so the newest-file pick below finds the one file
        pattern = sprintf('%s_100reps.mat', modes{m});
        files = dir(fullfile(dataDir, pattern));
        if isempty(files), error('No file matching: %s', pattern); end
        [~, idx] = max([files.datenum]);
        d = load(fullfile(dataDir, files(idx).name), 'results');
        nReps = length(d.results);
        coopMat = zeros(nReps, nGensB);
        for r = 1:nReps
            coopMat(r, :) = d.results{r}.historyCooperation(1:nGensB);
        end
        meanCoopOverTime(m, :) = mean(coopMat, 1) * 100;
        meanCoopSEM(m, :)      = (std(coopMat, 0, 1) ./ sqrt(nReps)) * 100;
        successOverTime(m, :)  = mean(coopMat >= 0.95, 1) * 100;
        finalCoopByMode{m}     = coopMat(:, end);
    end

    D.gens             = 1:nGensB;
    D.coopCiLo         = max(meanCoopOverTime - 1.96 * meanCoopSEM, 0);
    D.coopCiHi         = min(meanCoopOverTime + 1.96 * meanCoopSEM, 100);
    D.meanCoopOverTime = meanCoopOverTime;
    D.finalCoopByMode  = finalCoopByMode;
    D.successOverTime  = successOverTime;
    D.mode_labels      = palette.labels;
    D.short_labels     = palette.short_labels;
    D.mode_colors      = palette.colors;
    D.nModes           = nModes;
end
