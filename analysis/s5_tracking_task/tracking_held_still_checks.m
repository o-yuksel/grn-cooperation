function R = tracking_held_still_checks()
% TRACKING_HELD_STILL_CHECKS  Three checks, without simulation, on the stored population tests with
% the target held still in designs 1 and 2: has the load leveled off, why did failed solves fail,
% and how near is each design to a load set by the mutation rate?
%
% Produces : SI S5.1 (the slopes of load and mean fitness over the load window; fitness variance
%            not recorded)
% Inputs   : the runs with 30 generations without mutation: publication/tracking_held_still_30gen/design1/
%              main.mat (design 1, a = 0.70; written by run_tracking_held_still) and
%              publication/tracking_held_still_30gen/design2/main.mat (design 2, a = 0.20; written by
%              run_tracking_held_still_design2)
% Outputs  : publication/tracking_held_still_30gen/design2/nosim_checks.txt (printed lines; replaced)
% Run      : tracking_held_still_checks()   (about 30 s; on Windows it lowers the process priority through
%            PowerShell, on other systems it leaves the priority alone)
% Seed     : none (no random draws)
%
% Per design: 1. the slopes of load and mean fitness over the load window, per 100 generations
% relative to the window mean; 2. for each replicate whose solve for (b, delta) failed (the output
% bias and the change to W33 and W44 that lower alpha or ||x*||), and one solved replicate, the
% solve repeated from 1 to 0.70 in steps of 0.01, and the first level that fails and why; 3. U, the
% mutations per offspring, the ceiling 1 - exp(-U) of a load set by the mutation rate, and the
% exponents of the per-generation mutation cost and the loads on the root-mean-square shift per
% mutation (`beta_rms`) and on the identity's factor 0.7.

    pub = fullfile(grn_data_file('publication'));
    nThreads = maxNumCompThreads(1); restoreThreads = onCleanup(@() maxNumCompThreads(nThreads));
    if ispc
        system(sprintf(['powershell -NoProfile -Command "(Get-Process -Id %d).PriorityClass' ...
            ' = ''BelowNormal''"'], feature('getpid')));      % this process only, by PID
    end
    files = {fullfile(pub, 'tracking_held_still_30gen/design1', 'main.mat'), ...
             fullfile(pub, 'tracking_held_still_30gen/design2', 'main.mat')};
    names = {'DESIGN 1 (a = 0.70)', 'DESIGN 2 (a = 0.20)'};
    txt = fullfile(pub, 'tracking_held_still_30gen/design2', 'nosim_checks.txt');
    if exist(txt, 'file'), delete(txt); end
    diary(txt); diary on;
    fprintf('=========== NO-SIMULATION CHECKS (descriptive: load plateau, failed solves, mutation-rate ceiling) ===========\n');
    for d = 1:2
        o = getfield(load(files{d}), 'out');
        fprintf('\n################ %s: %s ################\n', names{d}, files{d});
        R(d).plateau = check_plateau(o);               %#ok<AGROW>
        R(d).dropped = check_dropped(o);               %#ok<AGROW>
        R(d).regime = check_regime(o);                 %#ok<AGROW>
        clear o
    end
    diary off;
end


% Check 1: plateau of the load and of mean fitness over the load window
function S = check_plateau(o)
    cfg = o.cfg; J = [o.jobs{:}]; rs = cfg.repList; keep = o.verif.rep; G = cfg.nGen;
    fprintf('\n---- 1. PLATEAU OF LOAD AND MEAN FITNESS over the load window (slope per 100 generations / window mean; median [IQR] over the arm''s contrast replicates) ----\n');
    fprintf('population fitness variance: NOT RECORDED by the generator (job fields: %s); not recoverable without re-simulation.\n', ...
        strjoin(fieldnames(J(1))', ', '));
    lab = {'load per window point (12 pts)', 'load per generation 1 - meanFit/wMax (60 gens)', 'mean fitness per generation (60 gens)'};
    S = struct();
    for a = 1:numel(cfg.arms)
        sl = nan(numel(rs), 3);
        for q = 1:numel(rs)
            z = J([J.r] == rs(q) & [J.arm] == a);
            if isempty(z) || z.void || (a >= 3 && ~keep(q, a-2)), continue; end
            T = z.track([z.track.k] > G - cfg.loadWin & isfinite([z.track.wBest]));
            k = [T.k]; y = 1 - [T.wBar] ./ [T.wBest];
            sl(q, 1) = rel_slope(k, y);
            kg = G - cfg.loadWin + 1 : G;
            sl(q, 2) = rel_slope(kg, 1 - z.meanFit(kg) ./ z.wMax(kg));
            sl(q, 3) = rel_slope(kg, z.meanFit(kg));
        end
        S.(matlab.lang.makeValidName(cfg.arms{a})) = sl;
        fprintf('%-10s', cfg.arms{a});
        for v = 1:3
            fprintf(' | %s %+.4f [%+.4f, %+.4f]', lab{v}, median(sl(:, v), 'omitnan'), ...
                prctile(sl(:, v), 25), prctile(sl(:, v), 75));
        end
        fprintf(' | n %d | within +-5%%: %s\n', sum(isfinite(sl(:, 1))), ...
            yesno(all(abs(median(sl, 1, 'omitnan')) <= 0.05)));
    end
end


function s = rel_slope(k, y)
% REL_SLOPE  Least-squares slope of y on k, per 100 generations, relative to the mean of y.
    ok = isfinite(y);
    if sum(ok) < 3, s = NaN; return; end
    c = polyfit(k(ok), y(ok), 1); s = 100 * c(1) / mean(y(ok));
end


% Check 2: replicates left out because the (b, delta) solve on the best genotype failed, each
% with one solved replicate for comparison
function D = check_dropped(o)
    sp = [o.splits{:}]; P = o.P; E = o.E; cfg = o.cfg;
    fprintf('\n---- 2. DROPPED REPLICATES (knob solve on ThBest failed at m = %.1f) ----\n', cfg.m);
    nm = {'ALPHA-DOWN', 'XNORM-DOWN'}; arms = 'AX'; D = struct('r', {}, 'arm', {}, 'res', {});
    for k = 1:2
        bad = find(~[sp.(['ok' arms(k)])]);
        fprintf('%s: solved %d/%d; unsolved replicates %s\n', nm{k}, numel(sp) - numel(bad), numel(sp), ...
            mat2str([sp(bad).r]));
        cmp = find([sp.(['ok' arms(k)])], 1);
        for q = [bad, cmp]
            res = march(sp(q), arms(k), E, P);
            D(end+1) = struct('r', sp(q).r, 'arm', arms(k), 'res', res); %#ok<AGROW>
            fprintf('  replicate %d (%s): alpha0 %.4f, ||x*||0 %.4f | %s\n', sp(q).r, ...
                ternary(ismember(q, bad), 'UNSOLVED', 'solved, comparator'), sp(q).alpha0, sp(q).xn0, res.msg);
        end
    end
end


function res = march(s, arm, E, P)
% MARCH  Continuation from m = 1 to 0.70 in steps of 0.01, with run_tracking_held_still's Newton
% solve, each level checked from x0; reports the first level that fails and why.
    g = struct('Th', s.dom, 'U', E.U, 'wt', E.wt);
    geo0 = geom(g.Th, g.U, g.wt, 0, 1, P, []);
    z = [0; 0]; res.lastOK = 1; res.failM = NaN; res.msg = '';
    for m = 0.99:-0.01:0.695
        if arm == 'A', tA = m*geo0.alpha; tX = geo0.xn; else, tA = geo0.alpha; tX = m*geo0.xn; end
        [z2, ok, why] = newton_bd(g, P, tA, tX, z);
        if ok
            e = geom(with_knob(g.Th, z2), g.U, g.wt, z2(1), 1, P, []);
            if ~(e.valid && abs(e.alpha/tA - 1) <= 1e-6 && abs(e.xn/tX - 1) <= 1e-6)
                ok = false; why = ['v1 from x0 fails: ' valid_reason(e)];
            end
        end
        if ~ok
            res.failM = m; res.why = why;
            eL = geom(with_knob(g.Th, z), g.U, g.wt, z(1), 1, P, []);
            res.msg = sprintf(['first failing level m = %.2f (last solved %.2f at b = %.4f, delta = %.4f; ' ...
                'there x* = %s, min|h| %.3g) -- %s'], m, res.lastOK, z(1), z(2), ...
                mat2str(eL.X(:, 1)', 3), eL.hmin, why);
            return
        end
        z = z2; res.lastOK = m;
    end
    e = geom(with_knob(g.Th, z), g.U, g.wt, z(1), 1, P, []);
    res.msg = sprintf('solved to m = %.2f (b = %.4f, delta = %.4f; x* = %s)', res.lastOK, z(1), z(2), mat2str(e.X(:, 1)', 3));
end


function [z, ok, why] = newton_bd(g, P, tA, tX, z)
% NEWTON_BD  run_tracking_held_still's damped Newton solve on (b, delta), returning why it stops.
    ok = false; why = '';
    [r, geo] = bd_resid(g, P, z, tA, tX, []);
    if ~geo.valid, why = ['start invalid: ' valid_reason(geo)]; return; end
    for it = 1:40
        if max(abs(r)) < 1e-9, ok = true; return; end
        J = zeros(2); hs = 1e-6*max(1, abs(z));
        for k = 1:2
            zp = z; zp(k) = zp(k) + hs(k);
            J(:, k) = (bd_resid(g, P, zp, tA, tX, geo.X) - r) / hs(k);
        end
        s = -J \ r; lam = 1; moved = false; lastWhy = '';
        for h = 1:30
            zn = z + lam*s;
            [rn, gn] = bd_resid(g, P, zn, tA, tX, geo.X);
            if gn.valid && norm(rn) < norm(r), moved = true; break; end
            if ~gn.valid, lastWhy = valid_reason(gn); else, lastWhy = 'no residual decrease'; end
            lam = lam/2;
        end
        if ~moved
            why = sprintf(['Newton stalled at iteration %d, residual (log alpha, log ||x*||) = %s, ' ...
                'cond(J) %.3g; full step (%+.3g, %+.3g); smallest trial step rejected: %s'], ...
                it, mat2str(r', 3), cond(J), s(1), s(2), lastWhy);
            return
        end
        z = zn; r = rn; geo = gn;
    end
    ok = max(abs(r)) < 1e-9;
    if ~ok, why = sprintf('Newton not converged in 40 iterations, residual %s', mat2str(r', 3)); end
end


function w = valid_reason(geo)
    X = geo.X; w = {};
    if ~geo.conv, w{end+1} = 'fixed point not converged'; end
    if ~(geo.maxRe < 0), w{end+1} = sprintf('unstable (max Re eig %.3g)', geo.maxRe); end
    if ~(min(X(3, :)) > 0.5 && min(X(4, :)) > 0.5)
        w{end+1} = sprintf('output gene below 0.5 (x3 %.3g, x4 %.3g)', min(X(3, :)), min(X(4, :)));
    end
    if ~all(isfinite(X(:))), w{end+1} = 'non-finite'; end
    if isempty(w), w = 'valid'; else, w = strjoin(w, '; '); end
end


function [r, geo] = bd_resid(g, P, z, tA, tX, Xinit)
    geo = geom(with_knob(g.Th, z), g.U, g.wt, z(1), 1, P, Xinit);
    r = [log(max(geo.alpha, realmin)/tA); log(max(geo.xn, realmin)/tX)];
    if ~geo.valid, r = [Inf; Inf]; end
end


function geo = geom(Th, U, wt, b, gain, P, Xinit)
% GEOM  run_tracking_held_still's window-weighted alpha and ||x*|| (copied), also returning conv
% and hmin, the smallest |input| of any gene.
    n = P.total; C = [0 0 1 0; 0 0 0 1]; nE = size(U, 2);
    Ue = gain*U + [0; 0; b; b];
    if isempty(Xinit), X = integrate_fp(Th, Ue, P); else, X = Xinit; end
    [X, conv] = polish_fp(Th, Ue, P, X);
    a = 0; xn = 0; mr = -Inf; hmin = Inf;
    for e = 1:nE
        x = X(:, e); h = Th*x + Ue(:, e); act = h > 0; hmin = min(hmin, min(abs(h)));
        fp = zeros(n, 1); fp(act) = P.K ./ (P.K + h(act)).^2; Df = diag(fp);
        CP = C * ((P.mu*eye(n) - Df*Th) \ Df);
        a = a + wt(e)*norm(CP, 'fro'); xn = xn + wt(e)*norm(x);
        mr = max(mr, max(real(eig(Df*Th - P.mu*eye(n)))));
    end
    on = min(X(3, :)) > 0.5 && min(X(4, :)) > 0.5;
    geo = struct('alpha', a, 'xn', xn, 'maxRe', mr, 'X', X, 'conv', conv, 'hmin', hmin, ...
        'valid', conv && mr < 0 && on && all(isfinite(X(:))));
end


% Check 3: mutation rate, the ceiling of a load set by the rate, and the exponents on the shift
function Q = check_regime(o)
    cfg = o.cfg; P = o.P; J = [o.jobs{:}]; rs = cfg.repList; keep = o.verif.rep; G = cfg.nGen;
    nMut = P.total^2 - 2; p = P.mutation_prob; U = nMut * p;
    fprintf('\n---- 3. MUTATION RATE AND REGIME ----\n');
    fprintf('per offspring: %d mutable entries x p = %.2f -> U = %.3f mutations (P(>= 1) = %.3f); Haldane-Muller load 1 - exp(-U) = %.3f (upper bound: every mutation deleterious)\n', ...
        nMut, p, U, 1 - (1 - p)^nMut, 1 - exp(-U));
    nm = {'mutCost', 'mutLoad', 'load'}; lab = {'per-generation mutation cost', 'mutation load', 'load carried'};
    V = nan(numel(rs), 4, 4);                   % replicate x arm x [cost, mutation load, load carried, rms shift per mutation]
    for q = 1:numel(rs)
        for a = 1:4
            z = J([J.r] == rs(q) & [J.arm] == a);
            if isempty(z) || z.void || (a >= 3 && ~keep(q, a-2)), continue; end
            s = wstats(z, cfg);
            T = z.track([z.track.k] > G - cfg.loadWin & isfinite([z.track.beta_rms]));
            V(q, a, :) = [s.mutCost, s.mutLoad, s.load, mean([T.beta_rms])];
        end
    end
    fprintf('CONTROL (median over replicates): load carried %.4f, mutation load %.4f, per-generation cost %.4f, beta_rms %.4f | load carried / Haldane-Muller = %.3f\n', ...
        median(V(:,1,3), 'omitnan'), median(V(:,1,2), 'omitnan'), median(V(:,1,1), 'omitnan'), ...
        median(V(:,1,4), 'omitnan'), median(V(:,1,3), 'omitnan') / (1 - exp(-U)));
    Q.V = V;
    for a = 2:4
        rb = V(:, a, 4) ./ V(:, 1, 4);
        fprintf('%-10s beta_rms ratio vs CONTROL median %.3f [IQR %.3f, %.3f] (identity alpha*||x*|| ratio %s)', ...
            cfg.arms{a}, median(rb, 'omitnan'), prctile(rb, 25), prctile(rb, 75), ternary(a == 2, '1', '0.7'));
        for v = 1:3
            ro = V(:, a, v) ./ V(:, 1, v);
            eb = log(ro) ./ log(rb); ei = log(ro) ./ log(0.7);
            if a == 2
                fprintf(' | %s ratio %.3f', lab{v}, median(ro, 'omitnan'));
            else
                fprintf('\n           %-29s ratio %.3f [IQR %.3f, %.3f] | exponent on beta_rms %.2f [IQR %.2f, %.2f] | on the identity ratio 0.7 %.2f', ...
                    lab{v}, median(ro, 'omitnan'), prctile(ro, 25), prctile(ro, 75), ...
                    median(eb, 'omitnan'), prctile(eb, 25), prctile(eb, 75), median(ei, 'omitnan'));
            end
            Q.(sprintf('%s_%d', nm{v}, a)) = ro;
        end
        fprintf('\n');
    end
    fprintf('reading: Gaussian small-effect (Lande-type) predicts exponent 2 for the per-generation cost and 1 for the load; Haldane-Muller predicts 0 for the load.\n');
end


function s = wstats(z, cfg)
% WSTATS  run_tracking_held_still's window_stats (copied).
    G = numel(z.meanFit); idx = G - cfg.loadWin + 1 : G;
    tk = z.track([z.track.k] > G - cfg.loadWin & isfinite([z.track.wBest]));
    s.load = 1 - sum([tk.wBar]) / sum([tk.wBest]);
    n = numel(z.off.wBar); io = n - cfg.settleWin + 1 : n;
    s.offLoad = 1 - sum(z.off.wBar(io)) / sum(z.off.wBest(io));
    s.mutLoad = s.load - s.offLoad;
    s.mutCost = 1 - sum(z.meanFit(idx)) / sum(z.wPar(idx));
end


function s = yesno(c)
    s = ternary(c, 'yes', 'no');
end
