function ring_assortment()
% RING_ASSORTMENT  On rings, how much more often members meet members than their share
% predicts, and how much that adds to their payoff advantage over non-members against pairing
% the same population at random.
%
% Produces : assortment on rings, printed: SI Extended Materials and Methods, "Spatial reproduction"
% Inputs   : publication/sweep_spatial_N50_<setting>_100reps_<time>.mat   (ring_k1, ring_k2,
%              ring_k5, ring_k10, wellmixed; written by run_spatial_n50)
% Outputs  : publication/assortment_route.mat   (struct array out, one row per file, saved after
%            each file); one printed RESULT line per file
% Run      : ring_assortment   (the file names, time stamps included, are fixed in the code)
% Seed     : rng(20260921), default generator, set once; the bootstrap draws in file order
%            (membership_test seeds each wiring itself; simulate_pair draws nothing)
%
% Nothing evolves here: each stored snapshot at generations 100, 200, 300, 401 and 500 is scored
% twice, with each individual's ring neighbors (radius k) and with all other individuals. A
% replicate is the mean over its snapshots that hold both members and non-members; replicates
% with none are left out, and intervals are percentile bootstrap over replicates (10,000
% resamples). The well-mixed file is scored as if on a ring with k = 2, as a control.
    d0 = fullfile(grn_data_file('publication'));
    S = { 'ring k=1 (in the paper)',  'sweep_spatial_N50_ring_k1_100reps_20260411_104446.mat', 1; ...
          'ring k=2 (in the paper)',  'sweep_spatial_N50_ring_k2_100reps_20260411_104544.mat', 2; ...
          'ring k=5',                 'sweep_spatial_N50_ring_k5_100reps_20260411_104638.mat', 5; ...
          'ring k=10',                'sweep_spatial_N50_ring_k10_100reps_20260411_104737.mat', 10; ...
          'CONTROL well-mixed, k=2 pseudo', 'sweep_spatial_N50_wellmixed_100reps_20260411_104345.mat', 2};
    wantG = [100 200 300 401 500]; Msq = ones(4); nBoot = 10000; rng(20260921);
    out = struct('name',{},'nRep',{},'a',{},'kin',{},'dLoc',{},'dGlob',{},'gain',{},'success',{});
    for f = 1:size(S,1)
        d = load(fullfile(d0,S{f,2}),'results','config'); c = d.config; k = S{f,3}; N = c.populationSize;
        R = c.payoffCC; Sp = c.payoffCD; T = c.payoffDC; P = c.payoffDD;
        sp = struct('sensingMode',c.sensingMode,'hillCoefficient',c.hillCoefficient,'discreteOutput',c.discreteOutput, ...
            'michaelisConstant',c.michaelisConstant,'degradationRate',c.degradationRate,'timeStep',c.timeStep,'developmentSteps',600, ...
            'payoffDC',T,'sensorIntegrationRate',1.0,'externalSignal1',zeros(4,1),'externalSignal2',zeros(4,1));
        cssp = struct('sensingMode',c.sensingMode,'hillCoefficient',c.hillCoefficient,'discreteOutput',c.discreteOutput, ...
            'michaelisConstant',c.michaelisConstant,'degradationRate',c.degradationRate,'timeStep',c.timeStep,'developmentSteps',600,'backend','simulate_pair');
        off = [-k:-1, 1:k]; nR = numel(d.results);
        A = nan(nR,1); K = A; DL = A; DG = A; t0 = tic;
        for r = 1:nR
            rep = d.results{r}; gs = double(rep.historyGenerations(:))';
            acc = zeros(0,4);
            for g = wantG
                s = find(gs == g, 1); if isempty(s), continue; end
                M = double(rep.historyM(:,:,:,s)); B = double(rep.historyB(:,:,:,s));
                [Wu,~,ic] = unique(reshape(M.*B,16,N)','rows'); nU = size(Wu,1);
                % Members: distinct wirings passing the membership test (600 development steps)
                isMu = false(nU,1); for u = 1:nU, isMu(u) = membership_test(reshape(Wu(u,:),4,4),cssp).is_css; end
                isM = isMu(ic); nm = sum(isM); if nm == 0 || nm == N, continue; end
                % Payoff of one encounter for each pair of distinct wirings, both developed for
                % 600 steps from the fixed start (simulate_pair), at the file's payoffs
                Pi = nan(nU);                                   % Pi(u,v) = payoff of type u against type v
                for u = 1:nU
                    for v = u:nU
                        o = simulate_pair(Msq,reshape(Wu(u,:),4,4),[],Msq,reshape(Wu(v,:),4,4),[],sp);
                        switch o, case 'CC', pu=R; pv=R; case 'CD', pu=Sp; pv=T; case 'DC', pu=T; pv=Sp; otherwise, pu=P; pv=P; end
                        Pi(u,v) = pu; Pi(v,u) = pv;
                    end
                end
                % Per snapshot: a = P(neighbor is a member | focal is a member) - (n_m - 1)/(N - 1)
                % (assortment); kin = P(neighbor has the focal's wiring) - that wiring's share of
                % the others (clone clustering); dLoc = mean payoff of members minus non-members
                % against their neighbors; dGlob = the same against all others; gain = dLoc - dGlob
                fr = accumarray(ic,1,[nU 1]);
                piLoc = nan(N,1); piGlob = nan(N,1); nbM = nan(N,1); nbK = nan(N,1);
                for i = 1:N
                    nb = unique(mod(i-1+off, N)+1); nb(nb==i) = [];
                    piLoc(i) = mean(Pi(ic(i), ic(nb)));
                    w = fr; w(ic(i)) = w(ic(i)) - 1;            % all others
                    piGlob(i) = (Pi(ic(i),:) * w) / (N-1);
                    nbM(i) = mean(isM(nb)); nbK(i) = mean(ic(nb) == ic(i)) - (fr(ic(i))-1)/(N-1);
                end
                acc(end+1,:) = [mean(nbM(isM)) - (nm-1)/(N-1), mean(nbK), ...
                                mean(piLoc(isM)) - mean(piLoc(~isM)), mean(piGlob(isM)) - mean(piGlob(~isM))]; %#ok<AGROW>
            end
            if ~isempty(acc), m = mean(acc,1); A(r)=m(1); K(r)=m(2); DL(r)=m(3); DG(r)=m(4); end
        end
        ok = ~isnan(A); G = DL - DG;
        q = numel(out)+1; out(q).name = S{f,1}; out(q).nRep = sum(ok);
        out(q).a = bootci1(A(ok),nBoot); out(q).kin = bootci1(K(ok),nBoot);
        out(q).dLoc = bootci1(DL(ok),nBoot); out(q).dGlob = bootci1(DG(ok),nBoot); out(q).gain = bootci1(G(ok),nBoot);
        out(q).success = mean(cellfun(@(x) x.historyCooperation(end), d.results) >= 0.95);
        fprintf('RESULT %-44s reps %3d succ %.2f | a %+.3f [%+.3f,%+.3f] | kin %+.3f [%+.3f,%+.3f] | dLoc %+.3f [%+.3f,%+.3f] | dGlob %+.3f [%+.3f,%+.3f] | gain %+.3f [%+.3f,%+.3f]  (%.1f min)\n', ...
            S{f,1}, sum(ok), out(q).success, out(q).a, out(q).kin, out(q).dLoc, out(q).dGlob, out(q).gain, toc(t0)/60);
        save(fullfile(d0,'assortment_route.mat'),'out');
    end
    fprintf('DONE\n');
end

% Mean of x and its 95% percentile bootstrap interval over its entries (replicates)
function v = bootci1(x, nBoot)
    n = numel(x); if n == 0, v = [NaN NaN NaN]; return; end
    b = mean(x(randi(n, nBoot, n)), 2); v = [mean(x), prctile(b,[2.5 97.5])];
end
