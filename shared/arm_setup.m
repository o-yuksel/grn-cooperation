function [knob, shift, seed] = arm_setup(s, arm, cfg)
%ARM_SETUP  The dial values, optimum shift and seed of one arm of the factor-lowering test (alpha or ||x*||
%lowered at the split).
%
% Implements : arm 1 CONTROL, 2 the second control (seed2), 3 ALPHA, 4 XNORM; used by run_tracking_held_still
%              and run_tracking_moving_target
% Inputs     : s (one split), arm, cfg
% Outputs    : knob, shift, seed
% Seed       : none

    switch arm
        case 1, knob = [0 0]; shift = [0; 0]; seed = cfg.seed + s.r;
        case 2, knob = [0 0]; shift = [0; 0]; seed = cfg.seed2 + s.r;
        case 3, knob = s.knobA; shift = s.shiftA; seed = cfg.seed + s.r;
        case 4, knob = s.knobX; shift = s.shiftX; seed = cfg.seed + s.r;
    end
end
