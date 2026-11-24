function y = flip_coop(W, sp, i, d)
%FLIP_COOP  Whether an input of size d to one gene ends a network's cooperation with its clone.
%
% Implements : one encounter of the network with its clone, the input added to gene i of both; used by
%              safety_search_losses, tipping_input_flag and worked_example
% Inputs     : W (wiring), sp (simulate_pair settings), i (gene), d (input size, signed)
% Outputs    : y (true if the network no longer cooperates)
% Seed       : none

    sp.externalSignal1 = zeros(4, 1); sp.externalSignal1(i) = d; sp.externalSignal2 = sp.externalSignal1;
    E = ones(4); [~, f] = simulate_pair(E, W, [], E, W, [], sp); y = ~((f(1) - f(2)) > 0);
end
