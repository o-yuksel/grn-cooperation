function nW = open_workers(nW)
%OPEN_WORKERS  Opens a pool of worker processes with one thread each.
%
% Implements : an open pool of another size is closed first; used by response_dial and response_dial_networks
% Inputs     : nW (number of workers; default 3)
% Outputs    : nW
% Seed       : none

    if nargin < 1 || isempty(nW), nW = 3; end
    p = gcp('nocreate');
    if ~isempty(p) && p.NumWorkers ~= nW, delete(p); p = []; end
    if isempty(p), p = parpool('Processes', nW); end
    wait(parfevalOnAll(p, @maxNumCompThreads, 0, 1));
    fprintf('pool: %d workers, one thread each\n', nW);
end
