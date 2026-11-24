function out = tracking_cut_after_150(cfg)
% TRACKING_CUT_AFTER_150  With the target held still, how much does lowering alpha or ||x*|| to
% 0.70 cut the mutation load, read after 150 generations without mutation? Designs 1 and 2.
%
% Produces : the cuts in the mutation load (Table S14; SI S5.1; Fig. 6; Table S15; main text,
%            "Responsiveness carries a mutation load in a non-social trait"; Methods, "Non-social
%            tracking task"; Abstract and Significance); the two outputs are also read by
%            run_tracking_cut_vs_shift, run_tracking_cut_vs_mutation_size, tracking_load_exponents
%            and fig6_figS2_data (for fig6_tracking_load)
% Inputs   : publication/tracking_held_still_150gen/design1/main.mat (design 1, a = 0.70) and
%              design2/main.mat (design 2, a = 0.20): the runs with 150 generations without
%              mutation (written by run_tracking_generations_without_mutation, mode 'run')
%            through run_tracking_held_still: publication/tracking_runs/
%              tracking_runs_N1000.mat and publication/tracking_single_genotypes/lowering_settings.mat
%            cfg (optional): T, the generations without mutation kept (150); outRoot
%              (publication/tracking_cut_after_150); designs ([1 2])
% Outputs  : <outRoot>/design1/ and design2/: main.mat (the run with its analysis) and main.txt
%              (two header lines, then the analysis run_tracking_held_still prints)
% Run      : tracking_cut_after_150()   (about 15 s)
% Seed     : none here; run_tracking_held_still's paired bootstrap uses seed 20261002
%
% Each stored run is cut to its first T generations without mutation and handed unchanged to
% run_tracking_held_still's reanalysis path (`cfg.reanalyse`), which reads the load left over the
% last 10 of them (the mutation load is the load carried minus the load left) and computes every
% cut and its paired bootstrap (2,000 resamples). No evolution is run.

    if nargin < 1, cfg = struct(); end
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    pub = fullfile(grn_data_file('publication'));
    def = struct('T', 150, 'outRoot', fullfile(pub, 'tracking_cut_after_150'), 'designs', [1 2]);
    f = fieldnames(def);
    for i = 1:numel(f), if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end, end
    % Design 1 is `designIndex` 1 (a = 0.70), design 2 `designIndex` 3 (a = 0.20); tol is the tolerance each
    % design's run used to check that its control reproduces the stored fitness of generations
    % 291-350 (run_tracking_held_still)
    D = struct('designIndex', {1, 3}, 'tol', {0.0106, 0.0110}, 'src', {'design1', 'design2'}, 'name', {'design1', 'design2'});
    out = cell(1, 2);
    for k = cfg.designs
        out{k} = one_design(cfg, D(k), pub);
    end
end


function R = one_design(cfg, d, pub)
    S = load(fullfile(pub, 'tracking_held_still_150gen', d.src, 'main.mat')); o = S.out;
    nOffSource = o.cfg.nOff;
    o.cfg.nOff = cfg.T;
    % Keep each job's first T generations without mutation
    for j = 1:numel(o.jobs)
        z = o.jobs{j};
        if isempty(z.off), continue; end
        for fn = {'meanFit', 'wBar', 'wBest', 'varFit'}
            z.off.(fn{1}) = z.off.(fn{1})(1:cfg.T);
        end
        o.jobs{j} = z;
    end
    od = fullfile(cfg.outRoot, d.name); if ~exist(od, 'dir'), mkdir(od); end
    tmp = [tempname '.mat']; S2 = struct('out', o); save(tmp, '-struct', 'S2');
    c = struct('reanalyse', tmp, 'tolerance', d.tol, 'designIndex', d.designIndex, 'outDir', od);
    txt = evalc('R = run_tracking_held_still(c);');
    delete(tmp);
    hdr = sprintf(['=========== TRACKING, TARGET HELD STILL: REANALYSIS WITH A %d-GENERATION TAIL WITHOUT MUTATION ===========\n' ...
        'source %s/main.mat (nOff = %d, cut to the first %d generations without mutation) | runner reanalyse path, unchanged | a = %.2f | %s\n'], ...
        cfg.T, d.src, nOffSource, cfg.T, o.E.a, datestr(now, 'yyyy-mm-dd'));
    fid = fopen(fullfile(od, 'main.txt'), 'w'); fprintf(fid, '%s%s', hdr, txt); fclose(fid);
    out = R; save(fullfile(od, 'main.mat'), 'out'); %#ok<NASGU>
end
