function setup_process()
%SETUP_PROCESS  Sets this MATLAB process to one computational thread and, on Windows, BelowNormal priority.
%
% Implements : used by the three tracking runners of the factor-lowering test (alpha or ||x*|| lowered at the split)
% Inputs     : none
% Outputs    : none (prints the process id)
% Seed       : none

    maxNumCompThreads(1);
    pid = feature('getpid');
    if ispc
        system(sprintf(['powershell -NoProfile -Command "(Get-Process -Id %d).PriorityClass' ...
            ' = ''BelowNormal''"'], pid));
        fprintf('client PID %d set to BelowNormal, 1 computational thread\n', pid);
    else
        fprintf('client PID %d, 1 computational thread (priority not changed on this system)\n', pid);
    end
end
