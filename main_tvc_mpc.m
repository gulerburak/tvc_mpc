clear; clc; close all;
live_animation = true; % set to false to skip live animation during simulation
%%
% Figures are written by save_figure.m: vector PDF into report/figures/ and a
% 150 dpi PNG into assets/ for the README.

%% Physical Parameters
m = 50; % rocket mass [kg]
g = 9.81; % gravitational acceleration    [m/s^2]
T0 = m * g; % nominal thrust (= weight)     [N]   490.5 N
J = 80; % moment of inertia about CG    [kg.m^2]
l_tvc = 1.3; % CG to nozzle gimbal distance  [m]
vz_nom = 0; % nominal speed: 0 = hover [m/s]
Ts = 0.05; % sampling period  (20 Hz)      [s]
N = 40; % MPC prediction horizon  (2 s)

% Rocket visual parameters for animation
rp.L = 4.0; % total length   [m]
rp.W = 0.4; % body width     [m]
rp.cg_frac = 0.375; % fraction of L below CG  (CG at 1.5 m from base)
rp.l_tvc = l_tvc;

%% Linearization
%  Equilibrium:  theta=0, delta=0, vy=0, vz=vz_nom, dT=0
%  Error state:  e = x - x_ref(t)   where x_ref(t) = [0;vz_nom*t;0;vz_nom;0;0]
%
%  Linearised equations (small-angle):
%    ey_dot   = evy
%    ez_dot   = evz
%    evy_dot  = g*(etheta + delta)   [lateral force from pitch+gimbal]
%    evz_dot  = dT/m                 [axial force from thrust deviation]
%    etheta_dot = eq
%    eq_dot   = (T0*l_tvc/J)*delta  [pitch torque from gimbal]

nx = 6;
nu = 2;

Ac = zeros(6);
Ac(1, 3) = 1; % ey_dot = evy
Ac(2, 4) = 1; % ez_dot = evz
Ac(3, 5) = g; % evy_dot = g*etheta
Ac(5, 6) = 1; % etheta_dot = eq

Bc = zeros(6, 2);
Bc(3, 1) = g; % evy_dot += g*delta
Bc(4, 2) = 1 / m; % evz_dot += dT/m
Bc(6, 1) = T0 * l_tvc / J; % eq_dot += (T0*l_tvc/J)*delta

eigs_c = eig(Ac);
fprintf('Continuous-time eigenvalues:\n');
disp(real(eigs_c))

%% Discretisation
sysc = ss(Ac, Bc, eye(nx), zeros(nx, nu));
sysd = c2d(sysc, Ts, 'zoh');
Ad = sysd.A;
Bd = sysd.B;

eigs_d = abs(eig(Ad));
fprintf('Discrete-time eigenvalue magnitudes: min=%.6f  max=%.6f\n', ...
    min(eigs_d), max(eigs_d));

%% Controllability check
Cm = ctrb(Ad, Bd);
rk = rank(Cm);
fprintf('Controllability matrix rank: %d  (need %d)\n', rk, nx);

%% Cost matrices Q and R
%         ey    ez    evy  evz  etheta  eq
Q = diag([50, 500, 10, 50, 200, 20]);
%         delta  dT
R = diag([50, 0.01]);

%% Terminal ingredients
[P_inf, ~, ~] = dare(Ad, Bd, Q, R);
%             ey ez evy evz etheta eq
Q_lqr = diag([5, 50, 2, 0.1, 200, 20]);
%              delta  dT
R_lqr = diag([500, 0.1]);

[~, ~, Kd] = dare(Ad, Bd, Q_lqr, R_lqr);
K_lqr = -Kd; % sign convention

% verification of LQR closed-loop stability
A_cl = Ad + Bd * K_lqr;
rho_cl = max(abs(eig(A_cl)));
fprintf('\n--- Terminal ingredients ---\n');
fprintf('LQR closed-loop spectral radius: %.6f  (must be < 1)\n', rho_cl);
assert(rho_cl < 1, 'LQR is not stabilising');

%% Constraints
%  State constraints (error coordinates):
%    ey in [-1, 4] m  (lateral position error — see note below)
%    |ez| <= 100 m    (altitude deviation)
%    |evy| <= 30 m/s  (lateral velocity)
%    |evz| <= 30 m/s  (vertical velocity deviation)
%    |theta| <= 15°   (linearisation validity)
%    |q| <= 20°/s     (pitch rate)
%
%  Input constraints (hardware limits):
%    |delta| <= 5°    (gimbal mechanical limit)
%    |dT| <= 150 N    (thrust deviation limit)
%
%  Note on the lateral bound: the rocket starts 3 m off-axis and the corridor is
%  deliberately tight and asymmetric — only 1 m of overshoot past the target is
%  allowed. This is what separates the two controllers. MPC sees the wall coming
%  over its horizon and shapes the approach to respect it, whereas LQR has no
%  notion of the constraint and only stops when its command saturates against the
%  5° gimbal limit.

e_lb = [-1; -100; -30; -30; -15 * pi / 180; -20 * pi / 180];
e_ub = [4; 100; 30; 30; 15 * pi / 180; 20 * pi / 180];

u_lb = [-5 * pi / 180; -150];
u_ub = [5 * pi / 180; 150];

% Initial state w.r.t reference
% rocket starts at z=0, which is z_target below the hover reference
z_target = 20; % hover altitude [m]
e0 = [3; -z_target; 0; 0; 5 * pi / 180; 0];

%% MPC Simulation
T_sim = 250;
beta = 3;

x_ref_0 = [0; z_target; 0; vz_nom; 0; 0];

fprintf('\n Starting MPC Simulation (N=%d, T_sim=%d, beta=%g) ---\n', N, T_sim, beta);
tic;
[E_mpc, U_mpc, X_world, T_end_mpc, infeasible_mpc] = run_mpc_sim( ...
    e0, x_ref_0, Ad, Bd, Q, R, P_inf, beta, N, ...
    e_lb, e_ub, u_lb, u_ub, T_sim, vz_nom, m, g, J, l_tvc, T0, Ts, ...
    'LiveAnimation', live_animation, 'Rp', rp, 'SaveVideo', true, 'VideoFile', 'tvc_animation.avi');
t_mpc_sim = toc;

fprintf('MPC simulation done in %.2f \n', t_mpc_sim);
fprintf('Final error norm: %.6f\n', norm(E_mpc(:, end)));

%% LQR Simulation
T_end_lqr = T_sim;
E_lqr = zeros(nx, T_sim + 1);
U_lqr = zeros(nu, T_sim);
X_lqr_w = zeros(nx, T_sim + 1);

X_lqr_w(:, 1) = x_ref_0 + e0;
E_lqr(:, 1) = e0;

for k = 1:T_sim
    ek = E_lqr(:, k);
    u_k = K_lqr * ek;
    U_lqr(:, k) = u_k;

    % saturate to hardware limits to mimic reality for
    % for simulation but LQR still computed without constraints
    u_k = max(u_lb, min(u_ub, u_k));

    x_next = tvc_nonlinear_step(X_lqr_w(:, k), u_k, m, g, J, l_tvc, T0, Ts);
    X_lqr_w(:, k + 1) = x_next;
    x_ref_kp1 = [0; z_target + vz_nom * k * Ts; 0; vz_nom; 0; 0];
    E_lqr(:, k + 1) = x_next - x_ref_kp1;

    if norm(E_lqr(:, k + 1)) < 0.1
        T_end_lqr = k;
        break;
    end

end

E_lqr = E_lqr(:, 1:T_end_lqr + 1);
U_lqr = U_lqr(:, 1:T_end_lqr);
X_lqr_w = X_lqr_w(:, 1:T_end_lqr + 1);
fprintf('\nLQR simulation done.\n');

%% Horizon study

% horizons to study
N_study = [10, 20, 30, 40, 50];
conv_times = nan(1, numel(N_study));
peak_theta = nan(1, numel(N_study));
solve_t = zeros(1, numel(N_study));
infeas_step = zeros(1, numel(N_study)); % 0 means no infeasibility
results_h = cell(numel(N_study), 1); % for plotting

for ni = 1:numel(N_study)
    Ni = N_study(ni);
    t0h = tic;

    [E_h, ~, ~, T_end_h, inf_k] = run_mpc_sim(e0, x_ref_0, Ad, Bd, Q, R, P_inf, beta, Ni, ...
        e_lb, e_ub, u_lb, u_ub, T_sim, vz_nom, m, g, J, l_tvc, T0, Ts);
    solve_t(ni) = toc(t0h);
    results_h{ni} = E_h;
    infeas_step(ni) = inf_k;

    peak_theta(ni) = max(abs(rad2deg(E_h(5, :))));

    if inf_k > 0
        conv_times(ni) = inf_k * Ts;
        fprintf('  N=%2d: INFEASIBLE at k=%d (t=%.2f s)  peak|theta|=%.2f deg  wall=%.2f s\n', ...
            Ni, inf_k, conv_times(ni), peak_theta(ni), solve_t(ni));
    else
        conv_times(ni) = T_end_h * Ts;
        fprintf('  N=%2d: converged t=%.2f s  peak|theta|=%.2f deg  wall=%.2f s\n', ...
            Ni, conv_times(ni), peak_theta(ni), solve_t(ni));
    end

end

%%

% Horizon study plots
clrs_h = lines(numel(N_study));
figure('Name', 'Horizon Study', 'NumberTitle', 'off', 'Color', 'w');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
bar(N_study, conv_times, 0.5, 'FaceColor', [0.2 0.45 0.85]);
xlabel('Horizon N'); ylabel('Time [s]');
title('Settling / Infeasibility time vs N'); grid on;

nexttile;
bar(N_study, peak_theta, 0.5, 'FaceColor', [0.85 0.2 0.2]);
xlabel('Horizon N'); ylabel('Peak |\theta|  [deg]');
title('Peak pitch excursion vs N'); grid on;

nexttile;
bar(N_study, solve_t, 0.5, 'FaceColor', [0.2 0.72 0.2]);
xlabel('Horizon N'); ylabel('Total solve time [s]');
title('Computation time vs N'); grid on;

% State trajectory plots
ax_hey = nexttile;
hold on;
grid on;
ax_hth = nexttile;
hold on;
grid on;
ax_hnm = nexttile;
hold on;
grid on;

for ni = 1:numel(N_study)
    t_h = (0:size(results_h{ni}, 2) - 1) * Ts;
    lbl = sprintf('N=%d', N_study(ni));

    if infeas_step(ni) > 0
        ls = '--';
    else
        ls = '-';
    end

    plot(ax_hey, t_h, results_h{ni}(1, :), ls, 'Color', clrs_h(ni, :), 'LineWidth', 1.5, 'DisplayName', lbl);
    plot(ax_hth, t_h, rad2deg(results_h{ni}(5, :)), ls, 'Color', clrs_h(ni, :), 'LineWidth', 1.5, 'DisplayName', lbl);
    semilogy(ax_hnm, t_h, vecnorm(results_h{ni}), ls, 'Color', clrs_h(ni, :), 'LineWidth', 1.5, 'DisplayName', lbl);
end

yline(ax_hey, e_ub(1), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax_hey, e_lb(1), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax_hth, rad2deg(e_ub(5)), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax_hth, rad2deg(e_lb(5)), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
xlabel(ax_hey, 't [s]'); ylabel(ax_hey, 'e_y [m]'); title(ax_hey, 'Lateral error vs horizon'); legend(ax_hey);
xlabel(ax_hth, 't [s]'); ylabel(ax_hth, '\theta [deg]'); title(ax_hth, 'Pitch angle vs horizon'); legend(ax_hth);
xlabel(ax_hnm, 't [s]'); ylabel(ax_hnm, '||e||_2'); title(ax_hnm, 'Error norm vs horizon'); legend(ax_hnm);

save_figure(gcf, 'horizon_study');

%% Weight tuning study

% Tuning scenarios
W_cases = {
           diag([50, 500, 10, 50, 200, 20]), diag([50, 0.01]), '1';
           diag([50, 2000, 10, 50, 200, 20]), diag([50, 0.01]), '2';
           diag([50, 500, 10, 50, 200, 20]), diag([5, 0.01]), '3';
           diag([50, 500, 10, 50, 200, 20]), diag([500, 0.01]), '4';
           diag([50, 500, 10, 50, 20, 20]), diag([50, 0.01]), '5';
           };

% plot settings
clrs_w = num2cell(lines(size(W_cases, 1)), 2);

fig_wt = figure('Name', 'Weight Tuning', 'NumberTitle', 'off', 'Color', 'w');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax_wy = nexttile;
hold(ax_wy, 'on');
grid(ax_wy, 'on');
title(ax_wy, 'Lateral position  e_y');
ylabel(ax_wy, 'e_y [m]');
xlabel(ax_wy, 't [s]');

ax_wth = nexttile;
hold(ax_wth, 'on');
grid(ax_wth, 'on');
title(ax_wth, 'Pitch angle  \theta');
ylabel(ax_wth, '\theta [deg]');
xlabel(ax_wth, 't [s]');

ax_wz = nexttile; hold(ax_wz, 'on');
grid(ax_wz, 'on');
title(ax_wz, 'Altitude error  e_z');
ylabel(ax_wz, 'e_z [m]');
xlabel(ax_wz, 't [s]');

ax_wd = nexttile; hold(ax_wd, 'on');
grid(ax_wd, 'on');
title(ax_wd, 'Gimbal angle  \delta');
ylabel(ax_wd, '\delta [deg]');
xlabel(ax_wd, 't [s]');

ax_wn = nexttile; hold(ax_wn, 'on');
grid(ax_wn, 'on');
title(ax_wn, 'Error norm  ||e||_2');
ylabel(ax_wn, '||e||_2  (log)');
xlabel(ax_wn, 't [s]');

nW = size(W_cases, 1);
wt_settling = zeros(nW, 1);
wt_peak_pitch = zeros(nW, 1);
wt_peak_ey = zeros(nW, 1);
wt_gimbal_sat = zeros(nW, 1);
wt_thrust_sat = zeros(nW, 1);
wt_infeas = zeros(nW, 1); % step of first infeasibility, 0 if the case converged

% runs MPC simulation for each case and plot results
for wi = 1:nW
    Qi = W_cases{wi, 1}; Ri = W_cases{wi, 2}; lbl = W_cases{wi, 3};
    [Pi, ~, Kdi] = dare(Ad, Bd, Qi, Ri);
    Ki = -Kdi;
    rho_i = max(abs(eig(Ad + Bd * Ki)));
    fprintf('  Case %d (%s): rho=%.4f\n', wi, lbl, rho_i);

    [E_w, U_w, ~, T_end_w, inf_kw] = run_mpc_sim(e0, x_ref_0, Ad, Bd, Qi, Ri, Pi, beta, N, ...
        e_lb, e_ub, u_lb, u_ub, T_sim, vz_nom, m, g, J, l_tvc, T0, Ts);

    if inf_kw > 0
        fprintf('  Case %d (%s): INFEASIBLE at k=%d (t=%.2f s)\n', wi, lbl, inf_kw, inf_kw * Ts);
    end

    t_w = (0:T_end_w) * Ts;

    % metrics. When a case goes infeasible the run stops early, so T_end_w*Ts is
    % the time it FAILED at, not a settling time — record the distinction rather
    % than reporting a failure as if it were the fastest case.
    wt_infeas(wi) = inf_kw;
    wt_settling(wi) = T_end_w * Ts;
    wt_peak_pitch(wi) = max(abs(rad2deg(E_w(5, :))));
    wt_peak_ey(wi) = max(abs(E_w(1, :)));

    if ~isempty(U_w)
        wt_gimbal_sat(wi) = sum(abs(U_w(1, :)) >= u_ub(1) * 0.99);
        wt_thrust_sat(wi) = sum(abs(U_w(2, :)) >= u_ub(2) * 0.99);
    end

    plot(ax_wy, t_w, E_w(1, :), 'Color', clrs_w{wi}, 'LineWidth', 1.8, 'DisplayName', lbl);
    plot(ax_wth, t_w, rad2deg(E_w(5, :)), 'Color', clrs_w{wi}, 'LineWidth', 1.8, 'DisplayName', lbl);
    plot(ax_wz, t_w, E_w(2, :), 'Color', clrs_w{wi}, 'LineWidth', 1.8, 'DisplayName', lbl);

    if ~isempty(U_w)
        stairs(ax_wd, (0:T_end_w - 1) * Ts, rad2deg(U_w(1, :)), 'Color', clrs_w{wi}, 'LineWidth', 1.8, 'DisplayName', lbl);
    end

    semilogy(ax_wn, t_w, vecnorm(E_w), 'Color', clrs_w{wi}, 'LineWidth', 1.8, 'DisplayName', lbl);
end

fprintf('\nWeight tuning metrics:\n');
fprintf('%-6s %-16s %-16s %-12s %-20s %-18s\n', 'Case', 'Outcome [s]', 'Peak |th| [deg]', 'Peak ey [m]', 'Gimbal sat steps', 'Thrust sat steps');

for wi = 1:nW

    if wt_infeas(wi) > 0
        outcome = sprintf('INFEAS @%.2f', wt_infeas(wi) * Ts);
    else
        outcome = sprintf('settled %.2f', wt_settling(wi));
    end

    fprintf('%-6s %-16s %-16.2f %-12.2f %-20d %-18d\n', W_cases{wi, 3}, outcome, wt_peak_pitch(wi), wt_peak_ey(wi), wt_gimbal_sat(wi), wt_thrust_sat(wi));
end

%%
% constraint lines
yline(ax_wy, e_ub(1), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax_wy, e_lb(1), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax_wth, rad2deg(e_ub(5)), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax_wth, rad2deg(e_lb(5)), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax_wz, e_ub(2), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax_wz, e_lb(2), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax_wd, rad2deg(u_ub(1)), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax_wd, rad2deg(u_lb(1)), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');

lg1 = legend(ax_wy, 'Location', 'northeast', 'FontSize', 8); lg1.ItemTokenSize = [15 8];
lg2 = legend(ax_wth, 'Location', 'northeast', 'FontSize', 8); lg2.ItemTokenSize = [15 8];
lg3 = legend(ax_wz, 'Location', 'northeast', 'FontSize', 8); lg3.ItemTokenSize = [15 8];
lg4 = legend(ax_wd, 'Location', 'northeast', 'FontSize', 8); lg4.ItemTokenSize = [15 8];
lg5 = legend(ax_wn, 'Location', 'northeast', 'FontSize', 8); lg5.ItemTokenSize = [15 8];

save_figure(fig_wt, 'weight_tuning_study');

%% Disturbance rejection

% define a lateral gust disturbance
gust_k = round(8.0 / Ts); % step index, applied after t=8s
gust_force = 200; % force
gust_dur = round(0.25 / Ts); % duration

% Gust acts on vy: delta_vy = F/m * Ts
D_gust = zeros(nx, T_sim);
D_gust(3, gust_k:gust_k + gust_dur - 1) = gust_force / m * Ts;

[E_dist, U_dist, ~, T_end_dist, infeas_dist] = run_mpc_sim( ...
    e0, x_ref_0, Ad, Bd, Q, R, P_inf, beta, N, ...
    e_lb, e_ub, u_lb, u_ub, T_sim, vz_nom, m, g, J, l_tvc, T0, Ts, ...
    'ConvTol', 5e-3, 'D', D_gust);

figure('Name', 'Disturbance Rejection', 'NumberTitle', 'off', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; grid on;
plot((0:T_end_dist) * Ts, E_dist(1, :), '-', 'Color', [0.85 0.5 0.2], 'LineWidth', 1.8, 'DisplayName', 'MPC + wind gust');
xline(gust_k * Ts, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Gust onset');
yline(e_ub(1), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(e_lb(1), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
ylabel('e_y [m]'); xlabel('t [s]');
title('Lateral error — disturbance rejection (200 N lateral gust for 0.25 s)');
legend('Location', 'northeast');

nexttile; hold on; grid on;
plot((0:T_end_dist) * Ts, rad2deg(E_dist(5, :)), '-', 'Color', [0.85 0.5 0.2], 'LineWidth', 1.8, 'DisplayName', 'MPC + wind gust');
xline(gust_k * Ts, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(rad2deg(e_ub(5)), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(rad2deg(e_lb(5)), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
ylabel('\theta [deg]'); xlabel('t [s]');
title('Pitch angle during gust');

save_figure(gcf, 'disturbance_rejection');

%% Observer: Measurement noise on output feedback
%

fprintf('\nAugmented Kalman filter \n');

Cd = [1 0 0 0 0 0; % ey
      0 1 0 0 0 0; % ez
      0 0 0 0 1 0; % etheta
      0 0 0 0 0 1]; % eq

ny = size(Cd, 1);

% Observability check
Ob = obsv(Ad, Cd);
fprintf('Observability rank: %d  (need %d)\n', rank(Ob), nx);

% Disturbance model
nd = 1;
Bd_dist = zeros(nx, nd);
Cd_dist = [0; 0; 1; 0];

% Augmented system
A_aug = [Ad, Bd_dist;
         zeros(nd, nx), eye(nd)];
B_aug = [Bd; zeros(nd, nu)];
C_aug = [Cd, Cd_dist];

% Observability check
Ob_aug = obsv(A_aug, C_aug);
fprintf('Augmented observability rank: %d  (need %d)\n', rank(Ob_aug), nx + nd);

% Process noise covariance
Qn_aug = blkdiag(diag([0.01; 0.01; 0.01; 0.01; (0.5 * pi / 180) ^ 2; (1 * pi / 180) ^ 2]), ...
    (1e-6) * eye(nd));

% Measurement noise covariance
sigma_y = 0.10;
sigma_th = 0.1 * pi / 180;
sigma_q = 0.5 * pi / 180;
Rn = diag([sigma_y ^ 2; sigma_y ^ 2; sigma_th ^ 2; sigma_q ^ 2]);

% Steady-state Kalman gain
[L_aug, ~] = dlqe(A_aug, eye(nx + nd), C_aug, Qn_aug, Rn);
L1 = L_aug(1:nx, :);
L2 = L_aug(nx + 1:end, :);

% Observer stability
eigs_obs = eig(A_aug - L_aug * C_aug);
fprintf('Kalman gain max singular value:    %.4f\n', max(svd(L_aug)));
fprintf('Augmented observer max|eig|:       %.6f  (must be < 1)\n', max(abs(eigs_obs)));

% constraints
e_lb_obs = e_lb; e_lb_obs(1) = -3; e_lb_obs(5) = -20 * pi / 180; e_lb_obs(6) = -25 * pi / 180;
e_ub_obs = e_ub; e_ub_obs(5) = 20 * pi / 180; e_ub_obs(6) = 25 * pi / 180;

% run MPC simulation with observer and measurement noise
rng(42);
[E_obs, U_obs, ~, T_end_obs, infeas_obs, ~, E_hat] = run_mpc_sim( ...
    e0, x_ref_0, Ad, Bd, Q, R, P_inf, beta, N, ...
    e_lb_obs, e_ub_obs, u_lb, u_ub, T_sim, vz_nom, m, g, J, l_tvc, T0, Ts, ...
    'Lkf', L_aug, 'Cd', Cd, 'Aaug', A_aug, 'Baug', B_aug, 'Caug', C_aug, ...
    'Rn', Rn, 'LiveAnimation', false, 'Rp', rp);
fprintf('Observer-based MPC converged at t = %.2f s\n', T_end_obs * Ts);

figure('Name', 'Observer MPC with Measurement Noise', 'NumberTitle', 'off', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Constraint lines must show the bounds this run actually enforced (the relaxed
% observer set), not the nominal ones — otherwise the true state appears to
% violate a limit that was never imposed.
nexttile; hold on; grid on;
plot((0:T_end_obs) * Ts, E_obs(1, :), '-', 'Color', [0.10 0.72 0.22], 'LineWidth', 1.8, 'DisplayName', 'True state');
plot((0:T_end_obs) * Ts, E_hat(1, :), '--', 'Color', [0.92 0.45 0.10], 'LineWidth', 1.2, 'DisplayName', 'KF estimate');
yline(e_ub_obs(1), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(e_lb_obs(1), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
ylabel('e_y [m]'); xlabel('t [s]');
title('Lateral error — observer-based MPC with measurement noise');
legend('Location', 'northeast');

nexttile; hold on; grid on;
plot((0:T_end_obs) * Ts, rad2deg(E_obs(5, :)), '-', 'Color', [0.10 0.72 0.22], 'LineWidth', 1.8, 'DisplayName', 'True state');
plot((0:T_end_obs) * Ts, rad2deg(E_hat(5, :)), '--', 'Color', [0.92 0.45 0.10], 'LineWidth', 1.2, 'DisplayName', 'KF estimate');
yline(rad2deg(e_ub_obs(5)), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(rad2deg(e_lb_obs(5)), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
ylabel('\theta [deg]'); xlabel('t [s]');
title('Pitch angle — observer tracking');
legend('Location', 'northeast');

save_figure(gcf, 'observer_mpc_noise');

%%  Final plots
fprintf('\n--- Generating final summary plots ---\n');

t_mpc_v = (0:T_end_mpc) * Ts;
t_lqr_v = (0:T_end_lqr) * Ts;

plot_results_tvc(t_mpc_v, E_mpc, U_mpc, t_lqr_v, E_lqr, U_lqr, ...
    e_lb, e_ub, u_lb, u_ub, vz_nom, z_target);
