function [S, SG, SL, nd] = parts_of(par, z, idx)
%PARTS_OF  Splits the selection on cooperation in one generation into its inherited part.
%
% Implements : S = Cov(n, z), S_G = Cov(n, g), S_L = Cov(n, z - g) (all 1/N), with n_i the offspring count of
%              individual i, z_i its cooperation and g_i the mean z over the individuals that share its wiring;
%              used by selection_inherited, selection_around_collapse and fig2_draw
% Inputs     : par (parent index of each individual of the next generation), z (cooperation), idx (wiring
%              index), all of one generation
% Outputs    : S, SG, SL, nd (the number of distinct wirings); a caller that asks for two outputs gets S and SG
% Seed       : none

    N = numel(z); z = double(z(:));
    n = accumarray(double(par(:)), 1, [N 1]);
    [u, ~, ic] = unique(double(idx(:)));
    gm = accumarray(ic, z) ./ accumarray(ic, 1);
    g = gm(ic);
    cv = @(a, b) mean((a - mean(a)) .* (b - mean(b)));
    S = cv(n, z); SG = cv(n, g); SL = cv(n, z - g); nd = numel(u);
end
