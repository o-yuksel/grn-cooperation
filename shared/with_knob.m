function Th2 = with_knob(Th, z)
%WITH_KNOB  A tracking wiring with the output genes' self-weights shifted.
%
% Implements : W33 and W44 (the self-weights of the two output genes) each raised by z(2), the second of the two
%              dial values; used
%              by run_tracking_held_still, run_tracking_moving_target and tracking_held_still_checks
% Inputs     : Th (wiring), z (the two dial values)
% Outputs    : Th2
% Seed       : none

    Th2 = Th; Th2(3,3) = Th2(3,3) + z(2); Th2(4,4) = Th2(4,4) + z(2);
end
