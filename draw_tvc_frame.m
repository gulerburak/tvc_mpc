function draw_tvc_frame(ax_rocket, ax_states, X_world, U_hist, ...
        pred_world, t_now, N, T_sim, Ts, rp)
    % DRAW_TVC_FRAME  Render one animation frame for the TVC MPC simulation
    %
    % ax_rocket  - 2-D (y, z) rocket-trajectory axes (left panel)
    % ax_states  - pitch / gimbal time-history axes  (right panel)
    % X_world    - world-frame state history  (6 x k)  [y,z,vy,vz,theta,q]
    % U_hist     - input history  (2 x k-1) or []     [delta, dT]
    % pred_world - MPC predicted (y,z) positions (2 x N+1) or []
    % t_now      - current time [s]
    % N          - MPC horizon
    % T_sim      - total simulation steps
    % Ts         - sample period [s]
    % rp         - struct: rp.L, rp.W, rp.cg_frac, rp.l_tvc (rocket visuals)

    k = size(X_world, 2);
    t_vec = (0:k - 1) * Ts;

    c_past = [0.20 0.45 0.85]; % blue   - history
    c_horiz = [0.92 0.45 0.10]; % orange - MPC horizon
    c_ref = [0.10 0.72 0.22]; % green  - reference

    %% ---- Left panel: 2-D rocket trajectory ----
    cla(ax_rocket);
    hold(ax_rocket, 'on');
    grid(ax_rocket, 'on');

    z_min = min(X_world(2, :)) - rp.L * 3;
    z_max = max(X_world(2, :)) + rp.L * 4;
    y_pad = max(max(abs(X_world(1, :))) + rp.L * 2, rp.L * 3);

    % Reference vertical axis (y = 0)
    plot(ax_rocket, [0 0], [z_min z_max], '--', 'Color', c_ref, ...
        'LineWidth', 1.0, 'DisplayName', 'Ref: y = 0');

    % Closed-loop flight path
    plot(ax_rocket, X_world(1, :), X_world(2, :), '-', ...
        'Color', c_past, 'LineWidth', 2.0, 'DisplayName', 'Flight path');

    % MPC predicted trajectory
    if ~isempty(pred_world)
        plot(ax_rocket, pred_world(1, :), pred_world(2, :), '--o', ...
            'Color', c_horiz, 'MarkerSize', 3, 'MarkerFaceColor', c_horiz, ...
            'LineWidth', 1.2, 'DisplayName', sprintf('MPC horizon (N=%d)', N));
    end

    % Draw rocket body at current position
    theta_now = X_world(5, k);
    delta_now = 0;

    if ~isempty(U_hist) && size(U_hist, 2) >= k - 1 && k > 1
        delta_now = U_hist(1, k - 1);
    end

    draw_rocket(ax_rocket, X_world(1, k), X_world(2, k), theta_now, delta_now, rp);

    % Axes
    xlabel(ax_rocket, 'y  [m]  (lateral)');
    ylabel(ax_rocket, 'z  [m]  (altitude)');
    title(ax_rocket, sprintf('TVC Rocket  |  t = %.2f s  |  \\theta = %.2f°  |  \\delta = %.2f°', ...
        t_now, rad2deg(theta_now), rad2deg(delta_now)), 'FontSize', 10);
    xlim(ax_rocket, [-y_pad, y_pad]);
    ylim(ax_rocket, [z_min, z_max]);
    legend(ax_rocket, 'Location', 'northeast', 'FontSize', 8);
    set(ax_rocket, 'FontSize', 10);

    %% ---- Right panel: pitch angle, lateral error, gimbal vs time ----
    cla(ax_states);
    hold(ax_states, 'on');
    grid(ax_states, 'on');

    % Pitch angle [deg]
    plot(ax_states, t_vec, rad2deg(X_world(5, :)), '-', ...
        'Color', [0.85 0.15 0.15], 'LineWidth', 1.8, 'DisplayName', '\theta [deg]');

    % Lateral error y [m]
    plot(ax_states, t_vec, X_world(1, :), '-', ...
        'Color', c_past, 'LineWidth', 1.8, 'DisplayName', 'e_y [m]');

    % Gimbal angle [deg]
    if ~isempty(U_hist) && size(U_hist, 2) > 0
        t_u = (0:size(U_hist, 2) - 1) * Ts;
        plot(ax_states, t_u, rad2deg(U_hist(1, :)), ':', ...
            'Color', c_horiz, 'LineWidth', 1.8, 'DisplayName', '\delta [deg]');
    end

    % Pitch constraint lines (±15 deg)
    plot(ax_states, [0 T_sim * Ts], [15 15], 'r--', 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
    plot(ax_states, [0 T_sim * Ts], [-15 -15], 'r--', 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
    % Gimbal constraint lines (±5 deg)
    plot(ax_states, [0 T_sim * Ts], [5 5], ':', 'Color', c_horiz, 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
    plot(ax_states, [0 T_sim * Ts], [-5 -5], ':', 'Color', c_horiz, 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
    plot(ax_states, [0 T_sim * Ts], [0 0], 'k:', 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');

    xlabel(ax_states, 't  [s]');
    ylabel(ax_states, '[deg]  /  [m]');
    title(ax_states, '\theta [deg],  e_y [m],  \delta [deg]', 'FontSize', 10);
    xlim(ax_states, [0, T_sim * Ts]);
    legend(ax_states, 'Location', 'northeast', 'FontSize', 8);
    set(ax_states, 'FontSize', 10);
end

function draw_rocket(ax, y_cg, z_cg, theta, delta, rp)
    % Draw rocket body at world position (y_cg, z_cg) with body pitch theta
    % and nozzle gimbal angle delta.
    %
    %  rp.L        total rocket length [m]
    %  rp.W        body width [m]
    %  rp.cg_frac  fraction of L below CG  (CG at cg_frac*L from base)
    %  rp.l_tvc    distance from CG to nozzle gimbal [m]

    L = rp.L;
    W = rp.W;
    cf = rp.cg_frac;

    h_above = (1 - cf) * L; % body length above CG
    h_below = cf * L; % body length below CG

    % Body-frame to world-frame transform:
    %   y_w = y_cg + xb*cos(theta) + zb*sin(theta)
    %   z_w = z_cg - xb*sin(theta) + zb*cos(theta)
    b2w = @(xb, zb) [y_cg + xb .* cos(theta) + zb .* sin(theta); ...
                         z_cg - xb .* sin(theta) + zb .* cos(theta)];

    % ---- Fuselage ----
    fu = [b2w(-W / 2, -h_below), b2w(W / 2, -h_below), ...
              b2w(W / 2, h_above), b2w(-W / 2, h_above)];
    fill(ax, fu(1, :), fu(2, :), [0.72 0.72 0.78], 'EdgeColor', 'k', 'LineWidth', 1.0, 'HandleVisibility', 'off');

    % ---- Nose cone ----
    nose = [b2w(-W / 2, h_above), b2w(W / 2, h_above), b2w(0, h_above + 0.18 * L)];
    fill(ax, nose(1, :), nose(2, :), [0.88 0.22 0.22], 'EdgeColor', 'k', 'LineWidth', 1.0, 'HandleVisibility', 'off');

    % ---- Fins ----
    fin_h = 0.38 * h_below;
    fin_w = 0.85 * W;
    fin_L = [b2w(-W / 2, -h_below), ...
                 b2w(-W / 2 - fin_w, -h_below), ...
                 b2w(-W / 2, -h_below + fin_h)];
    fin_R = [b2w(W / 2, -h_below), ...
                 b2w(W / 2 + fin_w, -h_below), ...
                 b2w(W / 2, -h_below + fin_h)];
    fill(ax, fin_L(1, :), fin_L(2, :), [0.48 0.48 0.54], 'EdgeColor', 'k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    fill(ax, fin_R(1, :), fin_R(2, :), [0.48 0.48 0.54], 'EdgeColor', 'k', 'LineWidth', 0.8, 'HandleVisibility', 'off');

    % ---- Exhaust plume ----
    t_dir = [sin(theta + delta); cos(theta + delta)]; % thrust direction (world)
    p_dir = [cos(theta + delta); -sin(theta + delta)]; % perpendicular

    noz = b2w(0, -h_below); % gimbal pivot (tail base)

    L_pl = 1.6 * W;
    W_pl = 0.50 * W;

    pl_tip = noz - L_pl * t_dir;
    pl_left = noz + W_pl / 2 * p_dir;
    pl_right = noz - W_pl / 2 * p_dir;

    fill(ax, [pl_left(1) pl_right(1) pl_tip(1)], ...
        [pl_left(2) pl_right(2) pl_tip(2)], ...
        [1.00 0.55 0.10], 'EdgeColor', 'none', 'FaceAlpha', 0.82, 'HandleVisibility', 'off');

    pl_tip2 = noz - L_pl * 0.58 * t_dir;
    pl_left2 = noz + W_pl * 0.30 * p_dir;
    pl_right2 = noz - W_pl * 0.30 * p_dir;
    fill(ax, [pl_left2(1) pl_right2(1) pl_tip2(1)], ...
        [pl_left2(2) pl_right2(2) pl_tip2(2)], ...
        [1.00 0.98 0.72], 'EdgeColor', 'none', 'FaceAlpha', 0.95, 'HandleVisibility', 'off');

    % ---- CG marker ----
    plot(ax, y_cg, z_cg, '+', 'Color', 'k', 'MarkerSize', 6, 'LineWidth', 1.5, 'HandleVisibility', 'off');
end
