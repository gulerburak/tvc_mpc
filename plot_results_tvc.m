function plot_results_tvc(t_mpc, E_mpc, U_mpc, t_lqr, E_lqr, U_lqr, ...
        e_lb, e_ub, u_lb, u_ub, vz_nom, z_target)
    % PLOT_RESULTS_TVC  Static summary figures after TVC MPC simulation
    %
    % E_mpc, E_lqr are error-state trajectories (6 x k); world-frame positions
    % are reconstructed below by adding the reference back on.

    clr_safe = [0.85 1.00 0.85];
    clr_mpc = [0.20 0.45 0.85];
    clr_lqr = [0.85 0.20 0.20];

    tmax = max(t_mpc(end), t_lqr(end));

    %% Figure 1 — 2-D world trajectories
    figure('Name', '2-D World Trajectories', 'NumberTitle', 'off', 'Color', 'w');
    hold on; grid on; box on;

    % Reconstruct world-frame positions
    y_mpc = E_mpc(1, :);
    z_mpc = E_mpc(2, :) + z_target + vz_nom .* t_mpc;
    y_lqr = E_lqr(1, :);
    z_lqr = E_lqr(2, :) + z_target + vz_nom .* t_lqr;

    % Vertical reference line (y = 0)
    z_all = [z_mpc, z_lqr];
    plot([0 0], [min(z_all) - 2, max(z_all) + 2], '--', ...
        'Color', [0.10 0.72 0.22], 'LineWidth', 1.2, 'DisplayName', 'y = 0  (target)');

    plot(y_mpc, z_mpc, '-', 'Color', clr_mpc, 'LineWidth', 2.0, 'DisplayName', 'MPC');
    plot(y_lqr, z_lqr, '--', 'Color', clr_lqr, 'LineWidth', 2.0, 'DisplayName', 'LQR');

    % Start and end markers
    scatter(y_mpc(1), z_mpc(1), 100, 'k', 'filled', 'DisplayName', 'Start');
    scatter(y_mpc(end), z_mpc(end), 100, clr_mpc, 'filled', 'HandleVisibility', 'off');
    scatter(y_lqr(end), z_lqr(end), 100, clr_lqr, 'filled', 'HandleVisibility', 'off');

    xlabel('y  [m]  (lateral)');
    ylabel('z  [m]  (altitude)');
    title('World-frame 2-D trajectories: MPC vs LQR');
    legend('Location', 'northwest', 'FontSize', 9);
    axis equal;
    set(gca, 'FontSize', 10);

    save_figure(gcf, 'fig_2d_world_trajectories');

    %% Figure 2 — Lateral position and pitch angle
    figure('Name', 'Lateral & Pitch', 'NumberTitle', 'off', 'Color', 'w');
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    labels_sp = {'\ite_y\rm  [m]  (lateral error)', ...
                 '\it\theta\rm  [deg]  (pitch angle)'};
    idx_sp = [1, 5];
    lb_sp = [e_lb(1), rad2deg(e_lb(5))];
    ub_sp = [e_ub(1), rad2deg(e_ub(5))];
    scale_sp = [1, 180 / pi];

    for i = 1:2
        si = idx_sp(i);
        nexttile; hold on; grid on;
        patch([0 tmax tmax 0], [lb_sp(i) lb_sp(i) ub_sp(i) ub_sp(i)], ...
            clr_safe, 'EdgeColor', 'none', 'DisplayName', 'Feasible region');
        yline(ub_sp(i), 'k-', 'LineWidth', 1.2, 'DisplayName', 'Constraint');
        yline(lb_sp(i), 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        plot(t_mpc, scale_sp(i) * E_mpc(si, :), '-', 'Color', clr_mpc, 'LineWidth', 1.8, 'DisplayName', 'MPC');
        plot(t_lqr, scale_sp(i) * E_lqr(si, :), '--', 'Color', clr_lqr, 'LineWidth', 1.5, 'DisplayName', 'LQR');
        xlabel('t  [s]'); ylabel(labels_sp{i});
        title(labels_sp{i}, 'FontSize', 10);
        legend('Location', 'northeast', 'FontSize', 7);
    end

    save_figure(gcf, 'lateral_pos_and_pitch_angle');

    %% Figure 3 — Velocities
    figure('Name', 'Velocities', 'NumberTitle', 'off', 'Color', 'w');
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    vel_labels = {'\ite_{vy}\rm  [m/s]  (lateral velocity error)', ...
                  '\ite_{vz}\rm  [m/s]  (vertical velocity error)'};
    vel_idx = [3, 4];

    for i = 1:2
        si = vel_idx(i);
        nexttile; hold on; grid on;
        patch([0 tmax tmax 0], [e_lb(si) e_lb(si) e_ub(si) e_ub(si)], ...
            clr_safe, 'EdgeColor', 'none', 'DisplayName', 'Feasible');
        yline(e_ub(si), 'k-', 'LineWidth', 1.2, 'DisplayName', 'Constraint');
        yline(e_lb(si), 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        plot(t_mpc, E_mpc(si, :), '-', 'Color', clr_mpc, 'LineWidth', 1.8, 'DisplayName', 'MPC');
        plot(t_lqr, E_lqr(si, :), '--', 'Color', clr_lqr, 'LineWidth', 1.5, 'DisplayName', 'LQR');
        xlabel('t [s]'); ylabel(vel_labels{i}); title(vel_labels{i}, 'FontSize', 10);
        legend('Location', 'northeast', 'FontSize', 7);
    end

    save_figure(gcf, 'velocities');

    %% Figure 4 — Control inputs
    figure('Name', 'Control Inputs', 'NumberTitle', 'off', 'Color', 'w');
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    u_labels = {'\delta  [deg]  (gimbal angle)', '\DeltaT  [N]  (thrust deviation)'};
    u_scales = [180 / pi, 1];

    for i = 1:2
        nexttile; hold on; grid on;
        lb_u = u_lb(i) * u_scales(i);
        ub_u = u_ub(i) * u_scales(i);
        patch([0 tmax tmax 0], [lb_u lb_u ub_u ub_u], clr_safe, ...
            'EdgeColor', 'none', 'DisplayName', 'Feasible');
        yline(ub_u, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Constraint');
        yline(lb_u, 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

        if size(U_mpc, 2) > 0
            stairs(t_mpc(1:end - 1), u_scales(i) * U_mpc(i, :), '-', ...
                'Color', clr_mpc, 'LineWidth', 1.8, 'DisplayName', 'MPC');
        end

        if size(U_lqr, 2) > 0
            U_lqr_sat = min(u_ub(i), max(u_lb(i), U_lqr(i, :)));
            % saturated input (what plant receives)
            stairs(t_lqr(1:end - 1), u_scales(i) * U_lqr_sat, '-', ...
                'Color', clr_lqr, 'LineWidth', 2.0, 'DisplayName', 'LQR (saturated)');
            % out-of-bounds portion only (unsaturated command beyond limits)
            U_lqr_raw = u_scales(i) * U_lqr(i, :);
            U_oob = U_lqr_raw;
            U_oob(U_lqr_raw >= lb_u & U_lqr_raw <= ub_u) = NaN;
            stairs(t_lqr(1:end - 1), U_oob, '-', ...
                'Color', [0.55 0.00 0.85], 'LineWidth', 2.0, 'DisplayName', 'LQR (unsaturated)');
        end

        xlabel('t [s]'); ylabel(u_labels{i}); title(u_labels{i}, 'FontSize', 10);
        legend('Location', 'northeast', 'FontSize', 7);
    end

    save_figure(gcf, 'control_inputs');

    %% Figure 5 — Error-norm convergence
    figure('Name', 'Convergence', 'NumberTitle', 'off', 'Color', 'w');
    hold on; grid on;
    semilogy(t_mpc, vecnorm(E_mpc), '-', 'Color', clr_mpc, 'LineWidth', 1.8, 'DisplayName', 'MPC');
    semilogy(t_lqr, vecnorm(E_lqr), '--', 'Color', clr_lqr, 'LineWidth', 1.5, 'DisplayName', 'LQR');
    xlabel('t  [s]'); ylabel('log  ||e||_2');
    title('Error-norm convergence  (error-state coordinates)');
    legend; set(gca, 'FontSize', 10);

    save_figure(gcf, 'error_norm_convergence');

    %% Figure 6 — Phase portrait: lateral position vs pitch angle
    figure('Name', 'Phase Portrait', 'NumberTitle', 'off', 'Color', 'w');
    hold on; grid on;
    plot(rad2deg(E_mpc(5, :)), E_mpc(1, :), '-', 'Color', clr_mpc, 'LineWidth', 1.8, 'DisplayName', 'MPC');
    plot(rad2deg(E_lqr(5, :)), E_lqr(1, :), '--', 'Color', clr_lqr, 'LineWidth', 1.5, 'DisplayName', 'LQR');
    scatter(rad2deg(E_mpc(5, 1)), E_mpc(1, 1), 80, 'k', 'filled', 'DisplayName', 't=0');
    scatter(0, 0, 100, [0.1 0.7 0.2], 'p', 'filled', 'DisplayName', 'Setpoint');
    xline(rad2deg(e_lb(5)), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xline(rad2deg(e_ub(5)), 'k--', 'LineWidth', 0.8, 'DisplayName', 'Constraints');
    xlabel('\theta  [deg]'); ylabel('e_y  [m]');
    title('Phase portrait: lateral error vs pitch angle');
    legend('Location', 'northwest'); set(gca, 'FontSize', 10);

    save_figure(gcf, 'phase_portrait');
end
