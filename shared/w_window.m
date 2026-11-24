function w = w_window(Th, knob, shift, E, P)
%W_WINDOW  Window-weighted tracking fitness of one genotype.
%
% Implements : with the optimum shifted by shift, scored as the evolution scores it; used by
%              run_tracking_held_still and run_tracking_moving_target
% Inputs     : Th (wiring), knob (the two dial values), shift (optimum shift), E (environments), P (tracking
%              parameters)
% Outputs    : w (NaN if knob or shift is NaN)
% Seed       : none

    if any(isnan(knob)) || any(isnan(shift)), w = NaN; return; end
    p = window_pheno(with_knob(Th, knob), knob, E, P);
    d2 = sum((p - (E.O + shift)).^2, 1);
    w = exp(-P.selection_strength * d2) * E.wt(:);
end
