function A = group_selection_arms(files, tag)
% GROUP_SELECTION_ARMS  Does group selection let more populations hold cooperation? Success and the
% share of members in populations of 1,000 in 20 groups of 50 (group selection; the random-split
% control) or undivided, with the tests between them.
%
% Produces : the group-selection results of SI S2.7, Table S2 and Fig. 3; Methods, "Game and
%            evolution"; main text, "Mutual sensing lets gene networks escape the dilemma"
% Inputs   : files   optional struct of full paths, fields `fission` (group selection), `fission_sham`
%                    (the random-split control), `control` (undivided); default publication/
%                    group_two_level/group_fission_<arm>_N1000_185reps_500gen.mat
%                    (written by run_group_selection)
%            tag     names the output (default 'N1000')
% Outputs  : publication/group_two_level/group_fission_analysis_<tag>.mat   (struct A); prints
%            group_selection_tests's report
% Run      : group_selection_arms   (not timed; run_group_selection_arms ends by calling it)
% Seed     : as group_selection_tests: rng(42) for the sampled wirings of each arm, then
%            rng(20261007) (the default generator) for the bootstrap
%
% It checks that each file is an N = 1,000 run of its arm, then calls group_selection_tests. The
% comparison of success with the line's prediction interval is reported, not tested.

    dDir = fullfile(grn_data_file('publication'), 'group_two_level');
    if nargin < 1 || isempty(files)
        f = @(a) fullfile(dDir, sprintf('group_fission_%s_N1000_185reps_500gen.mat', a));
        files = struct('fission', f('fission'), 'fission_sham', f('random_split'), 'control', f('control'));
    end
    if nargin < 2 || isempty(tag), tag = 'N1000'; end
    outFile = fullfile(dDir, sprintf('group_fission_analysis_%s.mat', tag));
    check_design(files);

    A = group_selection_tests(files, tag);
    A.S2_role = 'reported only, not a prediction';
    save(outFile, 'A');

    fprintf('\n----- N = 1000 arms -----\n');
    fprintf('The comparison of fission success with the line''s prediction interval, printed above, is reported only.\n');
end


function check_design(files)
% Every file must be an N = 1,000 run of the arm it is filed under, all at the same generation
% count; the two group arms in groups of 50, the undivided arm without groups.
    arms = fieldnames(files);
    g = zeros(numel(arms), 1);
    for a = 1:numel(arms)
        d = load(files.(arms{a}), 'cfg');
        assert(d.cfg.populationSize == 1000, '%s: populationSize %d, not 1000', arms{a}, d.cfg.populationSize);
        isGroup = isfield(d.cfg, 'groupSize') && d.cfg.groupSize > 0;
        if strcmp(arms{a}, 'control')
            assert(~isGroup, 'control file has group structure');
        else
            assert(isGroup && d.cfg.groupSize == 50 && ...
                strcmp(stored_value('groupSelectionMode', d.cfg.groupSelectionMode), stored_value('groupSelectionMode', arms{a})), ...
                '%s: not a groupSize-50 %s run', arms{a}, arms{a});
        end
        g(a) = d.cfg.numGenerations;
    end
    assert(all(g == g(1)), 'arms differ in numGenerations');
end
