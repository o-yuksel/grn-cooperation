function p = grn_data_file(relpath, mode)
%GRN_DATA_FILE  Full path of a result file in the data archive.
%
%   P = GRN_DATA_FILE(RELPATH) returns the full path of RELPATH (a file or
%   folder, relative to the data root) and stops with an error naming the
%   file if it is not there. GRN_DATA_FILE('') returns the data root.
%   P = GRN_DATA_FILE(RELPATH, 'new') returns the path without the check,
%   for a file or folder a script is about to write.
%
%   The data root is the folder named by the environment variable GRN_DATA,
%   or <release root>/data when GRN_DATA is not set. It holds the data
%   archive (see README) with its folder layout unchanged.

    if nargin < 1 || isempty(relpath), relpath = ''; end
    if nargin < 2, mode = ''; end

    root = getenv('GRN_DATA');
    if isempty(root)
        root = fullfile(fileparts(mfilename('fullpath')), 'data');
    end
    p = char(fullfile(root, relpath));

    if ~strcmp(mode, 'new') && ~isfile(p) && ~isfolder(p)
        error('grn_data_file:missing', ...
            'Missing data file %s: download the data archive (see README) into %s', ...
            relpath, root);
    end
end
