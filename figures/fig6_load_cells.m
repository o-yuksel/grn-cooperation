function varargout = fig6_load_cells(what, varargin)
% FIG6_LOAD_CELLS  The per-mutation load cells of Fig. 6A.
%
% Implements : mode 'data', the load cells of Fig. 6A (called by fig6_figS2_data, for
%              fig6_tracking_load)
% Inputs     : what   'data'; cfg (optional): dataDir (default publication/tracking_runs
%                     under the data root), nBoot (10,000)
%              <dataDir>/tracking_runs_N1000_load.mat   (written by tracking_load_vs_identity)
% Outputs    : TD with L (the load cells); prints each cell
% Seed       : each bootstrap draws from its own RandStream('twister') with a fixed seed set by
%              the quantity, arm and design (e.g. 44003 + 10 arm + design)
%
% A load cell is one arm (target moving or held still) and one design (output self-activation a):
% x is the mean predicted load (L_identLin, from alpha and ||x*|| alone) and y the mean measured
% load over the replicates, each with a 95% percentile bootstrap interval.

    switch what
        case 'data', [varargout{1:nargout}] = tracking_data(varargin{:});
        otherwise, error('fig6_load_cells:request', 'unknown request: %s', what);
    end
end


function TD = tracking_data(cfg)
%TRACKING_DATA  the load cells of the N = 1000 scale-up.
    if nargin < 1 || isempty(cfg), cfg = struct(); end
    if ~isfield(cfg,'dataDir'), cfg.dataDir = fullfile(grn_data_file('publication'), ...
            'tracking_runs'); end
    if ~isfield(cfg,'nBoot'),   cfg.nBoot   = 10000; end
    TD.L = load_panel_data(cfg);
end


% ===================== data =========================================
function D = load_panel_data(cfg)
%LOAD_PANEL_DATA  the eight N = 1000 load cells, measured against predicted.
    D.arms = {'moving','fixed'};
    D.x = []; D.xlo = []; D.xhi = []; D.L = []; D.Llo = []; D.Lhi = [];
    D.r = []; D.rlo = []; D.rhi = []; D.f = []; D.ai = [];
    L = load(fullfile(cfg.dataDir,'tracking_runs_N1000_load.mat'));
    o = L.out; D.aList = o.aList(:)'; D.nReps = o.nReps; D.nPop = 1000;
    fprintf('\nLoad panel, eight cells at N = %d:\n', D.nPop);
    fprintf('%-7s %-6s %9s %18s %11s %18s %8s %18s\n', ...
        'arm','a','L','L 95% CI','predicted','predicted 95% CI','ratio','ratio 95% CI');
    for ai = 1:2
        for f = 1:numel(D.aList)
            yv = o.S.(D.arms{ai}).L(f,:);
            xv = o.S.(D.arms{ai}).L_identLin(f,:);
            mx = mean(xv,'omitnan'); my = mean(yv,'omitnan');
            cx = boot_mean_ci(xv, cfg.nBoot, 44003 + 10*ai + f);
            cy = boot_mean_ci(yv, cfg.nBoot, 45003 + 10*ai + f);
            ci = boot_ratio_ci(yv, xv, cfg.nBoot, 66003 + 10*ai + f);
            D.x(end+1) = mx;  D.xlo(end+1) = cx(1); D.xhi(end+1) = cx(2);
            D.L(end+1) = my;  D.Llo(end+1) = cy(1); D.Lhi(end+1) = cy(2);
            D.r(end+1) = my/mx; D.rlo(end+1) = ci(1); D.rhi(end+1) = ci(2);
            D.f(end+1) = f; D.ai(end+1) = ai;
            fprintf('%-7s %-6.2f %9.4f  [%7.4f,%7.4f] %11.4f  [%7.4f,%7.4f] %8.3f   [%6.3f, %6.3f]\n', ...
                D.arms{ai}, D.aList(f), my, cy(1), cy(2), mx, cx(1), cx(2), D.r(end), ci(1), ci(2));
        end
    end
    fprintf('  ratio over the eight cells: %.3f to %.3f\n', min(D.r), max(D.r));
    fprintf('  measured load over the eight cells: %.4f to %.4f (a factor of %.0f)\n', ...
        min(D.L), max(D.L), max(D.L)/min(D.L));
end


function ci = boot_ratio_ci(y, x, nBoot, seed)
%BOOT_RATIO_CI  paired: one resample of the replicate index serves both means.
    ok = ~isnan(y) & ~isnan(x); y = y(ok); x = x(ok); n = numel(y);
    if n == 0, ci = [NaN NaN]; return; end
    st = RandStream('twister','Seed',seed);
    br = nan(nBoot,1);
    for b = 1:nBoot
        idx = randi(st, n, n, 1);
        br(b) = mean(y(idx)) / mean(x(idx));
    end
    ci = prctile(br, [2.5 97.5]);
end
