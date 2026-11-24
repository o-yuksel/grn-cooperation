function tf = is_window_point(k, cfg)
%IS_WINDOW_POINT  Whether generation k is a recorded point of the tracking runs.
%
% Implements : every trajEvery generations before the load window, every loadEvery in it; used by
%              run_tracking_held_still and run_tracking_moving_target
% Inputs     : k (generation), cfg (nGen, loadWin, trajEvery, loadEvery)
% Outputs    : tf
% Seed       : none

    inWin = k > cfg.nGen - cfg.loadWin;
    tf = (~inWin && mod(k, cfg.trajEvery) == 0) || (inWin && mod(k, cfg.loadEvery) == 0);
end
