function P = tracking_config()
% TRACKING_CONFIG  Base parameters of the non-social tracking task.
%
% Implements : the tracking network of Methods, "Non-social tracking task": the dynamics of Eq. 3
%              (dx_i/dt = f(sum_j W_ij x_j + u_i) - mu x_i, f(v) = max(v, 0) / (K + max(v, 0)))
%              on four genes, inputs 1-2 receiving the environmental signal u and the phenotype
%              read from outputs 3-4, which must track a circling optimum. The caller sets its
%              own population size, generations, horizon and mutation size.
% Inputs     : none
% Outputs    : P   struct of task parameters
% Seed       : none (no random numbers)

    % Network (Eq. 3)
    P.total              = 4;            % genes: inputs 1-2, outputs 3-4
    P.mu                 = 0.10;         % degradation rate mu
    P.K                  = 1.0;          % Michaelis constant K
    P.input_selfinhib    = -0.65;        % founder input-gene self-inhibition (SI Extended
                                         % Materials and Methods, "The non-social tracking task")

    % Task (Methods, "Non-social tracking task")
    P.trajectory_radius  = 0.25;         % radius of the optimum's circle
    P.selection_strength = 50.0;         % gamma: fitness = exp(-gamma d^2)

    % Development: Euler steps to the settled state (numerical settings, not model parameters)
    P.x0                 = ones(4,1)*0.1;
    P.max_dev_steps      = 5000;
    P.dev_dt             = 0.1;
    P.dev_tolerance      = 1e-6;

    % Evolution
    P.n_pop              = 80;
    P.n_generations      = 350;
    P.mutation_prob      = 0.08;         % per-entry mutation probability
    P.mutation_sd        = 0.012;        % standard deviation of a mutation
    P.horizon            = 60;           % last generations over which the load is read

    % Angular step omega of the optimum per generation; omega = 0 holds it still
    P.omega_static       = 0.0;
    P.omega_veryslow     = (2*pi)/600;
    P.omega_slow         = (2*pi)/300;
    P.omega_mid          = (2*pi)/100;
    P.omega_fast         = (2*pi)/36;
end
