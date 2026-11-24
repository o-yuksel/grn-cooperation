function sh = phenotype_shift(Th, knob, E, P)
%PHENOTYPE_SHIFT  The window-mean output shift of a genotype caused by the two dials.
%
% Implements : on the evolution's own developmental map (steady_batch: Euler from x0 to the dev_tolerance stop),
%              which is what fitness is scored on; used by run_tracking_held_still and
%              run_tracking_moving_target
% Inputs     : Th (wiring), knob (the two dial values), E (environments), P (tracking parameters)
% Outputs    : sh (2 x 1 shift; NaN if knob is NaN)
% Seed       : none

    if any(isnan(knob)), sh = [NaN; NaN]; return; end
    p0 = window_pheno(Th, [0 0], E, P); p1 = window_pheno(with_knob(Th, knob), knob, E, P);
    sh = (p1 - p0) * E.wt(:);
end
