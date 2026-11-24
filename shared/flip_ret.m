function y = flip_ret(W, H, sp, i, d)
%FLIP_RET  Whether an input of size d to one gene makes a network cooperate with a defector.
%
% Implements : one encounter of the network with the defector H, the input added to gene i of the network; used
%              by safety_search_losses and tipping_input_flag
% Inputs     : W (wiring), H (the defector's wiring), sp (simulate_pair settings), i (gene), d (input size,
%              signed)
% Outputs    : y (true if the network cooperates, so it no longer retaliates)
% Seed       : none

    sp.externalSignal1 = zeros(4, 1); sp.externalSignal1(i) = d;
    E = ones(4); [~, f] = simulate_pair(E, W, [], E, H, [], sp); y = (f(1) - f(2)) > 0;
end
