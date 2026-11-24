function params = publication_params()
%PUBLICATION_PARAMS  The paper's run settings for its standard setting, as one struct.
%
%   PARAMS = PUBLICATION_PARAMS() returns the settings of saturating contrast (saturating
%   output, contrast sensing): N = 1000, 500 generations, 100 replicates, mutation rate 0.005,
%   selection power 1.0, T = 5, degradation rate 0.10, K_M = 1.0, up to 300 development steps.
%   The values are not set here: they come from publication_config('continuous_output'), the
%   file the runs read, so the two cannot disagree.

    params = publication_config('continuous_output');
end
