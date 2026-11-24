% TRACKING  The network of the non-social tracking task: founder wiring, settled state x*,
% responsiveness alpha and the mean shift per mutation.
%
% Implements : the tracking network of Methods, "Non-social tracking task" and SI Extended
%              Materials and Methods, "The non-social tracking task": the dynamics of Eq. 3,
%              inputs 1-2 receiving the environmental signal u, the phenotype read from outputs 3-4
% Inputs     : Th   4-by-4 wiring, Th(i,j) the influence of gene j on gene i
%              u    environmental input to each gene (4 values)
%              P    parameters from tracking_config
% Outputs    : per method, below (alpha and the mean shift per mutation as in SI S5)
% Seed       : the mean shift per mutation (`fragility`) draws its mutations from the caller's stream
classdef tracking

    methods (Static)

        function Th = build_theta(a, w, P)
            % Founder wiring: input self-inhibition P.input_selfinhib, output self-activation a
            % (the design's fixed parameter, which sets ||x*||), input-to-output gain w
            Th = zeros(P.total);
            Th(1,1) = P.input_selfinhib; Th(2,2) = P.input_selfinhib;
            Th(3,3) = a; Th(4,4) = a;
            Th(3,1) = w; Th(4,2) = w;
        end

        function x = steady(Th, u, P)
            % Settled state x* under constant input u: Euler steps of Eq. 3 from P.x0 until the
            % squared change in one step is below P.dev_tolerance, or P.max_dev_steps steps
            x = P.x0; K = P.K; mu = P.mu; dt = P.dev_dt;
            for s = 1:P.max_dev_steps
                xo = x; b = Th*x + u; b(b < 0) = 0;
                dx = (b ./ (K + b)) - mu*x; x = x + dt*dx;
                if sum((x - xo).^2) < P.dev_tolerance; break; end
            end
        end

        function [u_path, pos, opt] = circle(omega, cc, n_gen, P)
            % Moving optimum: the signal circles at radius P.trajectory_radius with angular step
            % omega per generation and enters inputs 1-2 (u = pos); the optimum is the founder's
            % unstimulated output cc offset by the same displacement
            ang = (0:n_gen-1) * omega;
            pos = [P.trajectory_radius*cos(ang); P.trajectory_radius*sin(ang)];
            u_path = zeros(P.total, n_gen); u_path(1:2,:) = pos;
            opt = cc + pos;
        end

        function alpha = responsiveness(Th, u, P)
            % Responsiveness alpha: the phenotype's response to an input at gene 1 and at gene 2,
            % by finite differences of size 1e-3, averaged over the two
            p0 = tracking.steady(Th, u, P); p0 = p0([3 4]);
            d = 1e-3; acc = 0;
            for i = 1:2
                up = u; up(i) = up(i) + d;
                xp = tracking.steady(Th, up, P);
                acc = acc + norm(xp([3 4]) - p0) / d;
            end
            alpha = acc / 2;
        end

        function beta = fragility(Th, u, P, freeze_diag, ndraws)
            % Mean shift of the phenotype under one application of the task's mutation
            % distribution (each entry mutated with probability P.mutation_prob by
            % N(0, P.mutation_sd^2)), over ndraws draws. freeze_diag = true keeps the output
            % self-activations W33 = W44 fixed, as in the designs.
            p0 = tracking.steady(Th, u, P); p0 = p0([3 4]);
            acc = 0;
            for k = 1:ndraws
                Tm = Th; mask = rand(size(Tm)) < P.mutation_prob;
                Tm(mask) = Tm(mask) + randn(nnz(mask), 1) * P.mutation_sd;
                if freeze_diag; Tm(3,3) = Th(3,3); Tm(4,4) = Th(4,4); end
                xm = tracking.steady(Tm, u, P);
                acc = acc + norm(xm([3 4]) - p0);
            end
            beta = acc / ndraws;
        end

        function w0 = pick_w0(a, P)
            % Founder input gain w (grid 0.1 to 9) that maximizes generation-1 fitness
            % exp(-gamma d^2) over 36 positions of the optimum, so that every design starts
            % from its best founder
            npd = 36; ang = (0:npd-1) * (2*pi)/36;
            pos = [P.trajectory_radius*cos(ang); P.trajectory_radius*sin(ang)];
            u1  = zeros(P.total, npd); u1(1:2,:) = pos;
            wg = 0.1:0.1:9; best = -inf; w0 = 1;
            for w = wg
                Th = tracking.build_theta(a, w, P);
                cc = tracking.steady(Th, zeros(P.total,1), P); cc = cc([3 4]);
                d2 = zeros(1, npd);
                for g = 1:npd
                    xe = tracking.steady(Th, u1(:,g), P);
                    d2(g) = sum((xe([3 4]) - (cc + pos(:,g))).^2);
                end
                f = mean(exp(-P.selection_strength * d2));
                if f > best; best = f; w0 = w; end
            end
        end

    end
end
