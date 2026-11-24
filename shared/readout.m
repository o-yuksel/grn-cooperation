function r = readout(Wa, Wb, ea, eb, sp)
%READOUT  The settled expression states of a network and its partner under given external inputs.
%
% Implements : one simulate_pair encounter with externalSignal1 and externalSignal2 set; used by response_dial,
%              response_dial_from_rest and fig4b_retaliation_from_rest (model/response_off keeps its own copy,
%              so the model folder needs
%              no other folder)
% Inputs     : Wa, Wb (the two wirings), ea, eb (external inputs), sp (simulate_pair settings)
% Outputs    : r (the settled state of the first network and its partner, as simulate_pair returns it; the
%              readout is the C entry minus the D entry of the first network's part)
% Seed       : none

    sp.externalSignal1 = ea; sp.externalSignal2 = eb; E = ones(4);
    [~, fs] = simulate_pair(E, Wa, [], E, Wb, [], sp); r = fs(:);
end
