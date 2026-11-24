function [Ws, bias, Bs] = dial(M, B, W, xstar, s)
%DIAL  Sets the response dial of SI S3.5 on a network.
%
% Implements : the four sensor-to-output weights scaled by s, with the constant input (1 - s) * W(1:2, 3:4) *
%              x*(3:4) to genes C and D; used by response_dial, response_dial_as_input, response_dial_from_rest,
%              response_dial_recount and fig4b_retaliation_from_rest
% Inputs     : M, B, W (the network), xstar (its resting state), s (the dial setting)
% Outputs    : Ws (dialled wiring), bias (constant input), Bs (dialled weights)
% Seed       : none

    Bs = B; Bs(1:2, 3:4) = s * B(1:2, 3:4); Ws = M .* Bs;
    bias = [(1 - s) * W(1:2, 3:4) * xstar(3:4); 0; 0];
end
