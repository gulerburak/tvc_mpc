function [E, U, X, T_end, infeasible_k, U_seqs, E_hat] = run_mpc_sim( ...
        e0, x_ref_0, Ad, Bd, Q, R, P_inf, beta, N, ...
        e_lb, e_ub, u_lb, u_ub, ...
        T_sim, vz_nom, m, g, J, l_tvc, T0, Ts, opts)
    % This function runs closed-loop MPC simulation.
    % Parameters:
    % e0, x_ref_0        initial error state and world reference state
    % Ad, Bd             discrete TVC model
    % Q, R, P_inf        stage and terminal cost matrices
    % beta               terminal cost scaling (>= 1) applied to P_inf
    % N                  prediction horizon
    % e_lb/ub, u_lb/ub   state and input constraints
    % T_sim              simulation steps
    % vz_nom             nominal vz
    % m,g,J,l_tvc,T0,Ts  simulation parameters
    %
    % Optional parameters:
    % ConvTol        convergence tolerance
    % D              state disturbance
    % Lkf            Kalman gain
    % Cd             measurement matrix for original state
    % Aaug           augmented A matrix
    % Baug           augmented B matrix
    % Caug           augmented output matrix
    % Rn             measurement noise covariance
    % LiveAnimation  show animation during simulation
    % Rp             rocket visual parameters struct
    % SaveVideo      write the animation to a Motion JPEG AVI file
    % VideoFile      output filename
    % Verbose        print the applied input at every step
    %
    % Returns:
    % E            true error state trajectory
    % U            input trajectory
    % X            world-frame trajectory
    % T_end        index of last completed step
    % infeasible_k step of first infeasibility   (0 if never)
    % U_seqs       predicted input sequences
    % E_hat        KF estimated error state (equals E when no observer is used)

    arguments
        e0
        x_ref_0
        Ad
        Bd
        Q
        R
        P_inf
        beta
        N
        e_lb
        e_ub
        u_lb
        u_ub
        T_sim
        vz_nom
        m
        g
        J
        l_tvc
        T0
        Ts
        opts.ConvTol (1, 1) double = 0.1
        opts.D = []
        opts.Lkf = []
        opts.Cd = []
        opts.Aaug = []
        opts.Baug = []
        opts.Caug = []
        opts.Rn = []
        opts.LiveAnimation (1, 1) logical = false
        opts.Rp = []
        opts.SaveVideo (1, 1) logical = false
        opts.VideoFile = 'tvc_animation.avi'
        opts.Verbose (1, 1) logical = false
    end

    use_observer = ~isempty(opts.Lkf);

    nx = size(Ad, 1);
    nu = size(Bd, 2);

    E = zeros(nx, T_sim + 1); E(:, 1) = e0;
    E_hat = zeros(nx, T_sim + 1); E_hat(:, 1) = e0;
    U = zeros(nu, T_sim);
    X = zeros(nx, T_sim + 1); X(:, 1) = x_ref_0 + e0;
    z_ref_0 = x_ref_0(2); % base altitude of reference (0 for climb from origin, z_target for hover)
    U_seqs = zeros(nu, N, T_sim);
    T_end = T_sim;
    infeasible_k = 0;

    e_hat = e0;
    use_aug = use_observer && ~isempty(opts.Aaug);

    if use_aug
        xi_hat = [e0; zeros(size(opts.Aaug, 1) - nx, 1)];
    end

    if opts.LiveAnimation
        fig = figure('Color', 'w', 'Position', [40 60 1400 560], 'Name', 'TVC MPC');
    end

    vw = [];

    if opts.SaveVideo

        if ~opts.LiveAnimation
            fig = figure('Color', 'w', 'Position', [40 60 1400 560], 'Name', 'TVC MPC', 'Visible', 'off');
        end

        [fdir, fname] = fileparts(opts.VideoFile);
        avi_file = fullfile(fdir, [fname '.avi']);
        vw = VideoWriter(avi_file, 'Motion JPEG AVI');
        vw.FrameRate = round(1 / Ts);
        open(vw);
    end

    for k = 1:T_sim
        [u_k, U_seq, ok] = solve_mpc_tvc(e_hat, Ad, Bd, Q, R, P_inf, beta, N, ...
            e_lb, e_ub, u_lb, u_ub);

        if ~ok
            fprintf('--- INFEASIBLE at k=%d (t=%.2f s) ---\n', k, k * Ts);
            fprintf('  e_hat : ey=%.3f  ez=%.3f  vy=%.3f  vz=%.3f  th=%.3f deg  q=%.3f deg/s\n', ...
                e_hat(1), e_hat(2), e_hat(3), e_hat(4), rad2deg(e_hat(5)), rad2deg(e_hat(6)));
            fprintf('  E_true: ey=%.3f  ez=%.3f  vy=%.3f  vz=%.3f  th=%.3f deg  q=%.3f deg/s\n', ...
                E(1, k), E(2, k), E(3, k), E(4, k), rad2deg(E(5, k)), rad2deg(E(6, k)));
            e_next0 = Ad * e_hat;
            fprintf('  Ad*e_hat (u=0): ey=%.3f  [lb=%.3f ub=%.3f]\n', e_next0(1), e_lb(1), e_ub(1));
            warning('MPC infeasible at k=%d (t=%.2f s)', k, k * Ts);
            infeasible_k = k;
            T_end = k - 1;
            break;
        end

        U(:, k) = u_k;
        U_seqs(:, :, k) = U_seq;

        if opts.Verbose
            fprintf('k=%3d, t=%5.2f s, u=[%7.2f, %7.2f]\n', k, k * Ts, u_k(1), u_k(2));
        end

        x_next = tvc_nonlinear_step(X(:, k), u_k, m, g, J, l_tvc, T0, Ts);

        if ~isempty(opts.D)
            x_next = x_next + opts.D(:, k);
        end

        X(:, k + 1) = x_next;
        x_ref_kp1 = [0; z_ref_0 + vz_nom * k * Ts; 0; vz_nom; 0; 0];
        E(:, k + 1) = x_next - x_ref_kp1;

        if use_observer
            ny_obs = size(opts.Cd, 1);
            noise = sqrt(diag(opts.Rn)) .* randn(ny_obs, 1);
            y_meas = opts.Cd * E(:, k + 1) + noise;

            if use_aug
                xi_hat_pred = opts.Aaug * xi_hat + opts.Baug * u_k;
                xi_hat = xi_hat_pred + opts.Lkf * (y_meas - opts.Caug * xi_hat_pred);
                e_hat = xi_hat(1:nx);
            else
                e_hat_pred = Ad * e_hat + Bd * u_k;
                e_hat = e_hat_pred + opts.Lkf * (y_meas - opts.Cd * e_hat_pred);
            end

        else
            e_hat = E(:, k + 1);
        end

        E_hat(:, k + 1) = e_hat;

        if opts.LiveAnimation || opts.SaveVideo
            pred_world = zeros(2, N + 1);
            pred_world(:, 1) = X(1:2, k);
            e_tmp = E(:, k);

            for j = 1:N
                e_tmp = Ad * e_tmp + Bd * U_seq(:, j);
                pred_world(1, j + 1) = e_tmp(1);
                pred_world(2, j + 1) = z_ref_0 + vz_nom * (k - 1 + j) * Ts + e_tmp(2);
            end

            clf(fig);
            ax_rocket = subplot(1, 2, 1);
            ax_states = subplot(1, 2, 2);
            draw_tvc_frame(ax_rocket, ax_states, X(:, 1:k + 1), U(:, 1:k), ...
                pred_world, k * Ts, N, T_sim, Ts, opts.Rp);
            drawnow;

            if ~isempty(vw)
                writeVideo(vw, getframe(fig));
            end

            if opts.LiveAnimation
                pause(0.01);
            end

        end

        if norm(E(:, k + 1)) < opts.ConvTol
            T_end = k;
            break;
        end

    end

    if ~isempty(vw)
        close(vw);
        [fdir, fname] = fileparts(opts.VideoFile);
        fprintf('Animation saved to %s\n', fullfile(fdir, [fname '.avi']));

        if ~opts.LiveAnimation
            close(fig);
        end

    end

    E = E(:, 1:T_end + 1);
    E_hat = E_hat(:, 1:T_end + 1);
    U = U(:, 1:T_end);
    X = X(:, 1:T_end + 1);
    U_seqs = U_seqs(:, :, 1:T_end);
end
