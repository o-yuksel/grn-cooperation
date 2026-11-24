function drop_final_block(f, marker)
%DROP_FINAL_BLOCK  Removes the block a rerun of 'final' would repeat, so the table is not appended twice.
%
% Implements : cuts the file from the first occurrence of marker to its end; used by run_tracking_cut_vs_shift
%              and run_tracking_cut_vs_mutation_size
% Inputs     : f (report file), marker
% Outputs    : none (rewrites f)
% Seed       : none

    if ~exist(f, 'file'), return; end
    txt = fileread(f); i = strfind(txt, marker);
    if isempty(i), return; end
    fid = fopen(f, 'w'); fprintf(fid, '%s', regexprep(txt(1:i(1)-1), '\n+$', '')); fclose(fid);
end
