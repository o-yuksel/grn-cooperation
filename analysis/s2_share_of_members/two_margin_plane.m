function G = two_margin_plane(matFile, cfg)
% TWO_MARGIN_PLANE  Where do members sit between losing cooperation and losing retaliation?
% Computes a clone-weighted density over the two margins, their medians, and how far one mutation
% moves a member in that plane.
%
% Produces : SI S2.3, "In median, members sit nearer to losing retaliation than to losing
%            cooperation" (the medians, G.med, printed); one console line on how one mutation
%            moves a member
% Inputs   : matFile   publication/retaliation_margin_vs_rho.mat   (retaliation_margin_vs_rho; per
%                      member the two margins and the clone count, per mutation the shifts)
%            cfg       optional struct (defaults in the code): pHi (0.99, the clone-weighted
%                      percentile of each margin where the box ends), nDens (72 cells per axis),
%                      bwDiv (14: smoothing width = box size / 14), encl ([0.25 0.50 0.75 0.90], the
%                      weight shares the contours enclose), deadDr (1e-6, the shift of the
%                      cooperation margin below which a mutation counts as moving nothing)
% Outputs  : G, a struct: per member the two margins, clone count and architecture; the density and
%              contour levels; the medians, pooled (med) and per architecture (archMed); the slope
%              and size of one mutation's move
% Run      : G = two_margin_plane(fullfile(grn_data_file('publication'), ...
%              'retaliation_margin_vs_rho.mat'), struct())   (about 1 min)
% Seed     : none (no random draws)
%
% The cooperation margin r is the readout x_C - x_D against the wiring's own clone (r > 0 is the
% cooperation condition (C1)); the retaliation margin y = -worst_retal (the field of that name in
% retaliation_margin_vs_rho) is minus the readout against
% the worst-case mutant (y > 0 is the retaliation condition (C2)), so members fill the positive
% quadrant. Members kept are the stored ones with ||x*|| in [0.5, 14] and r >= 1e-6; a member with
% no defecting mutant has no y and is left out of the density and the medians.

    if ~exist(matFile, 'file')
        error('two_margin_plane:noData', ['%s not found. Pass the stored retaliation_margin_vs_rho.mat ' ...
            '(publication folder of the data archive) as matFile.'], matFile);
    end
    if ~isfield(cfg, 'pHi'),   cfg.pHi = 0.99; end
    if ~isfield(cfg, 'nDens'), cfg.nDens = 72; end
    if ~isfield(cfg, 'bwDiv'), cfg.bwDiv = 14; end
    if ~isfield(cfg, 'encl'),  cfg.encl = [0.25 0.50 0.75 0.90]; end
    if ~isfield(cfg, 'deadDr'), cfg.deadDr = 1e-6; end

    G = load_members(matFile);
    G = build_density(G, cfg);
    G = build_coupling(G, cfg);
    fprintf('SI S2.3 medians (clone-weighted): cooperation margin %.4f, retaliation margin %.4f\n', G.med);
end


% =========================================================================
function G = load_members(matFile)
% Per kept member: the two margins, the clone weight and the architecture; per draw: dr and the
% change in the retaliation margin, for the coupling.
    fprintf('Loading %s\n', matFile);
    L = load(matFile, 'R_all', 'modeNames');
    G.modeNames = L.modeNames;
    G.nM = numel(L.R_all);
    r = []; y = []; w = []; a = []; dx = []; dy = []; dw = []; df = [];
    G.nKept = zeros(G.nM, 1); G.nAll = zeros(G.nM, 1); G.nNoDef = zeros(G.nM, 1);
    for m = 1:G.nM
        B = L.R_all{m};
        Bk = B(logical([B.ok])');
        rm = double([Bk.r])';
        ym = -double([Bk.wr0])';       % the sign flip, applied once
        wm = double([Bk.count])';
        G.nAll(m) = numel(B); G.nKept(m) = numel(Bk); G.nNoDef(m) = sum(~isfinite(ym));
        r = [r; rm]; y = [y; ym]; w = [w; wm]; a = [a; m*ones(numel(rm),1)]; %#ok<AGROW>

        K = numel(Bk(1).dr_meas);
        % K draws per member: members in rows, draws in columns
        drm = reshape(vertcat(Bk.dr_meas), K, [])';
        wrm = reshape(vertcat(Bk.wr_m),    K, [])';
        dcu = -(wrm + repmat(ym, 1, K));   % change in the retaliation margin, -(wr_m - wr0)
        wd  = repmat(wm, 1, K);            % clone weight per draw, same layout as drm
        fd  = repmat(isfinite(ym), 1, K);  % draws of a member that has a retaliation margin
        dx = [dx; drm(:)]; dy = [dy; dcu(:)]; dw = [dw; wd(:)]; df = [df; fd(:)]; %#ok<AGROW>
    end
    G.r = r; G.y = y; G.w = w; G.a = a;
    G.dx = dx; G.dy = dy; G.dw = dw; G.dfin = df;
end


% =========================================================================
function G = build_density(G, cfg)
    s = isfinite(G.y);
    r = G.r(s); y = G.y(s); w = G.w(s);
    xhi = wq(r, w, cfg.pHi);
    yhi = wq(y, w, cfg.pHi);
    G.xhi = xhi; G.yhi = yhi;
    G.xlim = [-0.055*xhi, xhi];
    G.ylim = [-0.085*yhi, yhi];

    n = cfg.nDens;
    G.xe = linspace(0, xhi, n+1);
    G.ye = linspace(0, yhi, n+1);
    % Members outside the box are left out of the lattice rather than placed in the edge cells,
    % which with the edge correction below would inflate the top and right edges.
    ix = discretize(r, G.xe);
    iy = discretize(y, G.ye);
    ok = isfinite(ix) & isfinite(iy);
    G.nDrop = sum(~ok);
    G.dens = accumarray([iy(ok), ix(ok)], w(ok), [n n], @sum, 0);
    G.xc = 0.5*(G.xe(1:end-1) + G.xe(2:end));
    G.yc = 0.5*(G.ye(1:end-1) + G.ye(2:end));

    % Gaussian smoothing, bandwidth one bwDiv-th of each axis. Edge correction: the smoothed weight
    % is divided by the smoothed indicator of the box, so a cell at a margin's zero edge is not
    % pulled down by the zero padding. Contour levels: the lowest density whose superlevel set
    % holds each share encl of the clone weight.
    G.hx = xhi / cfg.bwDiv;
    G.hy = yhi / cfg.bwDiv;
    kx = gauss_kernel(G.hx / (G.xe(2)-G.xe(1)));
    ky = gauss_kernel(G.hy / (G.ye(2)-G.ye(1)));
    num = conv2(ky, kx, G.dens, 'same');
    den = conv2(ky, kx, ones(n), 'same');
    G.kde = num ./ max(den, eps);
    G.kdeN = G.kde / max(G.kde(:));
    G.encl = cfg.encl;
    G.levels = weight_levels(G.kdeN, G.dens, G.encl);
    G.bwDiv = cfg.bwDiv;

    G.med = [wq(r, w, 0.5), wq(y, w, 0.5)];
    G.archMed = nan(G.nM, 2);
    for m = 1:G.nM
        sm = (G.a == m) & isfinite(G.y);
        G.archMed(m,:) = [wq(G.r(sm), G.w(sm), 0.5), wq(G.y(sm), G.w(sm), 0.5)];
    end
end


% =========================================================================
function G = build_coupling(G, cfg)
% How one mutation moves a member: one readout serves both margins, so a shift dr of the
% cooperation margin comes with a change in the retaliation margin.
%   slope:   clone-weighted least-squares slope of the change in the retaliation margin on dr,
%            over the draws where both are defined (the mutant keeps (C1), so the membership test
%            scored its retaliation margin).
%   length:  clone-weighted median, quartiles and mean of |dr| over live draws (|dr| >= deadDr) of
%            members with a retaliation margin; below deadDr the drawn weight has no path to the
%            readout and the mutation moves nothing. The split is on dr, not on the response, and
%            the length is not restricted to draws with a defined retaliation change.
    ok = isfinite(G.dx) & isfinite(G.dy) & G.dw > 0;
    x = G.dx(ok); y = G.dy(ok); w = G.dw(ok) / sum(G.dw(ok));
    mx = sum(w.*x); my = sum(w.*y);
    G.slope = sum(w.*(x-mx).*(y-my)) / sum(w.*(x-mx).^2);
    G.nPair = sum(ok);

    live = isfinite(G.dx) & G.dfin & G.dw > 0 & abs(G.dx) >= cfg.deadDr;
    G.deadDr = cfg.deadDr;
    G.nLive = sum(live);
    v = abs(G.dx(live)); vw = G.dw(live);
    G.drMed  = wq(v, vw, 0.5);
    G.drQ1   = wq(v, vw, 0.25);
    G.drQ3   = wq(v, vw, 0.75);
    G.drMean = sum(vw.*v) / sum(vw);
    fprintf(['coupling: slope %.4f over %d draws; one mutation |dr| median %.4f ' ...
             '[%.4f, %.4f] MEAN %.4f over %d live draws (|dr| >= %g)\n'], ...
        G.slope, G.nPair, G.drMed, G.drQ1, G.drQ3, G.drMean, G.nLive, cfg.deadDr);
end


% =========================================================================
function lv = weight_levels(kdeN, dens, encl)
% Lowest density threshold whose superlevel set holds `encl` of the total
% weight, one per target, returned IN THE ORDER OF `encl` so a caller can pair
% a level with the target that produced it. A contour plot needs the levels ascending, so a
% caller that draws them sorts a copy: sorting here would mispair levels and targets silently,
% because a bigger enclosed share needs a LOWER threshold.
    [v, ord] = sort(kdeN(:), 'descend');
    c = cumsum(dens(ord)); c = c / c(end);
    lv = zeros(1, numel(encl));
    for i = 1:numel(encl)
        j = find(c >= encl(i), 1);
        if isempty(j), j = numel(v); end
        lv(i) = v(j);
    end
end


function k = gauss_kernel(sigmaCells)
% Gaussian, truncated at 3 sigma, normalized. Row vector; conv2(ky, kx, ...)
% applies it separably down the columns and across the rows.
    hw = max(1, ceil(3*sigmaCells));
    t = -hw:hw;
    k = exp(-0.5*(t/sigmaCells).^2);
    k = k / sum(k);
end
