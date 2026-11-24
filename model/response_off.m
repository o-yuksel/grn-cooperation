function varargout = response_off(action, varargin)
% RESPONSE_OFF  A network with its response to the partner cut (scaled by s = 0), held at its
% resting state.
%
% Implements : the scaled network of Methods, "Mutation tests in the game" (Fig. 4C), at s = 0:
%              the four connections from the sensor genes to C and D, W(1:2,3:4), are set to 0,
%              and a constant input to C and D replaces what they carried at the network's
%              resting state with its clone, so that this state is unchanged
% Inputs     : action    'off'     Woff = response_off('off', W): W with W(1:2,3:4) = 0
%                        'bias'    [b, ok, xstar] = response_off('bias', W, sp): the constant
%                                  input b = [W(1:2,3:4) * xstar(3:4); 0; 0] at the clone pair's
%                                  resting state xstar
%                        'params'  sp = response_off('params'): settings of the resting-state
%                                  search (saturating contrast, 600 development steps)
% Outputs    : as listed per action. ok: the clone pair settled (no gene changes by 1e-10 or more
%              in one step); if not, 'bias' returns b = NaN(4, 1), and callers class the network
%              as one whose resting state is undefined.
% Seed       : 'bias' saves the caller's random stream and restores it on return
    switch action
        case 'off'
            varargout{1} = cut_links(varargin{1});
        case 'bias'
            [varargout{1:3}] = rest_bias(varargin{1}, varargin{2});
        case 'params'
            varargout{1} = census_sp();
        otherwise
            error('response_off:action', 'unknown action %s', action);
    end
end


function Wo = cut_links(W)
% The four sensor-to-output connections W(1:2,3:4) set to 0
    Wo = W; Wo(1:2, 3:4) = 0;
end


function [b, ok, xstar] = rest_bias(W, sp)
% The clone pair's resting state xstar: develop the pair with simulate_pair, then settle it with
% the Euler map; b is the input the cut connections carried there.
    s = rng;
    f0 = readout(W, W, zeros(4, 1), zeros(4, 1), sp); x0 = f0(1:4);
    [xstar, ~, ok] = refine(W, zeros(4, 1), x0, sp);
    rng(s);
    if ok, b = [W(1:2, 3:4) * xstar(3:4); 0; 0]; else, b = NaN(4, 1); end
end


function sp = census_sp()
% Saturating contrast, K_M = 1.0, mu = 0.10, Delta t = 1.0, 600 development steps, g = 1.0
% (Table S1); no constant input.
    sp = struct('sensingMode', 'contrast', 'hillCoefficient', 1, 'discreteOutput', false, ...
        'michaelisConstant', 1.0, 'degradationRate', 0.1, 'timeStep', 1.0, 'developmentSteps', 600, ...
        'payoffDC', 0, 'sensorIntegrationRate', 1.0, 'externalSignal1', zeros(4, 1), 'externalSignal2', zeros(4, 1));
end


function [x, it, ok] = refine(W, bias, x, sp)
% Euler map of Eq. 3 for a network paired with its own clone (the partner's state equals x),
% plus a constant input, until the largest change of any gene in one step is below 1e-10
% (ok = true), or at most 1e6 steps.
    km = sp.michaelisConstant; mu = sp.degradationRate; dt = sp.timeStep; ok = false; it = 0;
    cap = 1e6;
    for it = 1:cap
        u = max(W * x + partner_signal(x, sp) + bias, 0);
        xn = max(x + dt * (u ./ (km + u) - mu * x), 0);
        ch = max(abs(xn - x)); x = xn;
        if ch < 1e-10, ok = true; return; end
    end
end


function s = partner_signal(x, sp)
% Sensor input of Eq. 4 under saturating contrast when the partner's state is x:
% S_C = -S_D = g (x_C^n - x_D^n) / (x_C^n + x_D^n), n = sp.hillCoefficient; 1e-12 keeps it
% defined when x_C = x_D = 0.
    n = sp.hillCoefficient; cn = x(1)^n; dn = x(2)^n; den = cn + dn + 1e-12;
    g = sp.sensorIntegrationRate * (cn / den - dn / den); s = [0; 0; g; -g];
end


function r = readout(Wa, Wb, ea, eb, sp)
% Settled state of the pair (Wa, Wb), with constant inputs ea and eb, from simulate_pair
    sp.externalSignal1 = ea; sp.externalSignal2 = eb; E = ones(4);
    [~, fs] = simulate_pair(E, Wa, [], E, Wb, [], sp); r = fs(:);
end
