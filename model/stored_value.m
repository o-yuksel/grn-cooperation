function v = stored_value(field, v, direction)
% STORED_VALUE  The current spelling of an option value read from a data file, or the stored
% spelling of a current value.
%
% Implements : no model step. The data files of the archive were written before the option values
%              took their present names; the configs and labels stored in them hold the old
%              spellings listed below. Code that compares a value it read from a data file passes
%              it through this function first; values the code sets itself use the current names.
% Inputs     : field      'sensingMode', 'conditionName', 'groupSelectionMode', 'refusalMode',
%                         'snapMode' or 'arm' (the label of a refusal arm or a run class, as stored
%                         in a run file's lab)
%              v          a char value, or a cell array of them (each mapped)
%              direction  'current' (default): a stored spelling -> the current one; a value that
%                         is already current is returned unchanged.
%                         'stored': a current value -> the spelling the data files use (for
%                         building the name of a stored file); a value with no older spelling is
%                         returned unchanged.
% Outputs    : v          the mapped value; for 'sensingMode' the result must be 'contrast' or
%                         'separate', and any other value stops with an error
% Seed       : none (no random numbers)

    if nargin < 3 || isempty(direction), direction = 'current'; end
    if iscell(v)
        for k = 1:numel(v), v{k} = stored_value(field, v{k}, direction); end
        return;
    end
    if isstring(v), v = char(v); end

    % The spellings stored in the data files (first column) and the current ones (second column)
    switch field
        case 'sensingMode'
            T = {'differential', 'contrast'; 'absolute', 'separate'};
        case 'groupSelectionMode'
            T = {'fission_sham', 'random_split'};
        case 'refusalMode'
            T = {'COOPTHIN', 'REFUSE_PARTIAL'; 'COOPEND', 'REFUSE_ALL'; 'RESPONSE', 'REFUSE_RESPONSE';
                 'ORDINARY', 'REFUSE_OTHER'; 'RATCHET', 'REFUSE_LESS_CONDITIONAL';
                 'RETALIATION', 'REFUSE_LESS_RETALIATING'};
        case 'snapMode'
            T = {'engine', 'encounter_rule'};
        case 'arm'
            T = {'thin050', 'refuse050'; 'thin082', 'refuse082'; 'thin090', 'refuse090';
                 'thin095', 'refuse095'; 'thin098', 'refuse098'; 'refuseall', 'all_refused';
                 'ratchet', 'refuse_less_conditional'; 'fission_sham', 'random_split';
                 'c1only', 'cooperate_only'};
        case 'conditionName'
            v = condition_name(v, direction);
            return;
        otherwise
            error('stored_value:field', 'no stored spellings for field %s', field);
    end
    if strcmp(direction, 'stored'), T = T(:, [2 1]); end
    k = find(strcmp(v, T(:, 1)), 1);
    if ~isempty(k), v = T{k, 2}; end
    if strcmp(field, 'sensingMode') && strcmp(direction, 'current') && ~any(strcmp(v, {'contrast', 'separate'}))
        error('stored_value:sensingMode', 'sensingMode must be ''contrast'' or ''separate'', got ''%s''', v);
    end
end


function v = condition_name(v, direction)
% conditionName records the setting and the overrides of a run (publication_config). Stored forms:
% 'ablation_011' and 'ablation_111' (both are the paper's setting, continuous expression, once
% publication_config has set its values), each optionally followed by _<field>=<value> overrides
% (sensingMode values among them); 'spatial_N50v2_<structure>'; 'group_fission_<arm>_confirm_N1000'.
    if strcmp(direction, 'stored')
        error('stored_value:direction', 'conditionName has no single stored spelling');
    end
    v = regexprep(v, '^ablation_(011|111)(?=_|$)', 'continuous_output');
    v = strrep(v, 'sensingMode=differential', 'sensingMode=contrast');
    v = strrep(v, 'sensingMode=absolute', 'sensingMode=separate');
    v = regexprep(v, '^spatial_N50v2_', 'spatial_N50_');
    v = regexprep(v, '^group_fission_fission_sham_', 'group_fission_random_split_');
    v = regexprep(v, '^(group_fission_\w+?)_confirm_N1000$', '$1_main_N1000');
end
