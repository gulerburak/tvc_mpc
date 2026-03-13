%% main_tvc_mpc.m — Thrust-Vector Control MPC (SC421 Project)
%%
%% System: 2-D rocket (lateral y, altitude z) controlled by gimbal angle
%%         (delta) and thrust deviation (dT) during vertical ascent.
%%
%% State   e = [ey; ez; evy; evz; etheta; eq]  — error from reference trajectory
%% Input   u = [delta; dT]                     — gimbal angle [rad], thrust deviation [N]
%%
%% Reference trajectory: x_ref(t) = [0; vz_nom*t; 0; vz_nom; 0; 0]
%% Error dynamics:       e(k+1) = Ad*e(k) + Bd*u(k)   (exact for linearised plant)
%%
%% Run:  addpath(genpath('../mpt3/tbxmanager')); tbxmanager restorepath
%%       cd tvc_mpc; main_tvc_mpc
%%
%% Requires: MATLAB + YALMIP + quadprog + MPT3

clear; clc; close all;

%% =========================================================================
%%  1. PHYSICAL PARAMETERS
%% =========================================================================
m      = 50;        % rocket mass                   [kg]
g      = 9.81;      % gravitational acceleration    [m/s^2]
T0     = m * g;     % nominal thrust (= weight)     [N]   490.5 N
J      = 80;        % moment of inertia about CG    [kg.m^2]
l_tvc  = 1.3;       % CG to nozzle gimbal distance  [m]
vz_nom = 15;        % nominal climb speed            [m/s]
Ts     = 0.05;      % sampling period  (20 Hz)      [s]
N      = 20;        % MPC prediction horizon  (1 s)

% Rocket visual parameters for animation
rp.L       = 4.0;    % total length   [m]
rp.W       = 0.4;    % body width     [m]
rp.cg_frac = 0.375;  % fraction of L below CG  (CG at 1.5 m from base)
rp.l_tvc   = l_tvc;

fprintf('=== TVC Rocket MPC ===\n');
fprintf('  m=%.0f kg  T0=%.1f N  J=%.0f kg.m^2  l_tvc=%.2f m\n', ...
        m, T0, J, l_tvc);
fprintf('  vz_nom=%.0f m/s  Ts=%.3f s  N=%d\n', vz_nom, Ts, N);

%% =========================================================================
%%  2. LINEARISED CONTINUOUS-TIME MODEL  (around vertical-flight equilibrium)
%%
%%  Equilibrium:  theta=0, delta=0, vy=0, vz=vz_nom, dT=0
%%  Error state:  e = x - x_ref(t)   where x_ref(t) = [0;vz_nom*t;0;vz_nom;0;0]
%%
%%  Linearised equations (small-angle):
%%    ey_dot   = evy
%%    ez_dot   = evz
%%    evy_dot  = g*(etheta + delta)   [lateral force from pitch+gimbal]
%%    evz_dot  = dT/m                 [axial force from thrust deviation]
%%    etheta_dot = eq
%%    eq_dot   = (T0*l_tvc/J)*delta  [pitch torque from gimbal]
%% =========================================================================
nx = 6;   nu = 2;

Ac         = zeros(6);
Ac(1,3)    = 1;          % ey_dot = evy
Ac(2,4)    = 1;          % ez_dot = evz
Ac(3,5)    = g;          % evy_dot = g*etheta
Ac(5,6)    = 1;          % etheta_dot = eq

Bc         = zeros(6,2);
Bc(3,1)    = g;           % evy_dot += g*delta
Bc(4,2)    = 1/m;         % evz_dot += dT/m
Bc(6,1)    = T0*l_tvc/J;  % eq_dot += (T0*l_tvc/J)*delta

fprintf('\n--- Linearised model ---\n');
eigs_c = eig(Ac);
fprintf('Continuous-time eigenvalues:');
fprintf('  %.4f', real(eigs_c)');
fprintf('\n  (All zero: chain of integrators — marginally stable)\n');

%% =========================================================================
%%  3. DISCRETISE  (zero-order hold)
%% =========================================================================
sysc = ss(Ac, Bc, eye(nx), zeros(nx,nu));
sysd = c2d(sysc, Ts, 'zoh');
Ad   = sysd.A;
Bd   = sysd.B;

eigs_d = abs(eig(Ad));
fprintf('Discrete-time eigenvalue magnitudes: min=%.6f  max=%.6f\n', ...
        min(eigs_d), max(eigs_d));

%% =========================================================================
%%  4. CONTROLLABILITY CHECK
%% =========================================================================
Cm   = ctrb(Ad, Bd);
rk   = rank(Cm);
fprintf('Controllability matrix rank: %d  (need %d)  — ', rk, nx);
if rk == nx
    fprintf('CONTROLLABLE\n');
else
    error('System is NOT controllable — check model parameters.');
end

%% =========================================================================
%%  5. COST MATRICES  (Q, R)
%%
%%  Design rationale:
%%    - Pitch angle (theta) is the most safety-critical state:  Q(5,5) = 200
%%    - Lateral position (y) is the guidance objective:          Q(1,1) = 5
%%    - Gimbal angle (delta) is physically limited to ±5°:       R(1,1) = 500
%%      Large R_delta discourages saturation and reduces wear.
%%    - Thrust deviation (dT) is cheap and effective:            R(2,2) = 0.1
%% =========================================================================
Q = diag([5,   0.5,  2,   0.1,  200,  20]);
%         ey   ez   evy  evz  etheta  eq
R = diag([500,  0.1]);
%         delta  dT

%% =========================================================================
%%  6. TERMINAL INGREDIENTS  (DARE  →  LQR terminal cost + invariant set)
%% =========================================================================
[P_inf, ~, Kd] = dare(Ad, Bd, Q, R);
K_lqr          = -Kd;    % u = K_lqr * e  (optimal LQR gain, u = -Kd*e)

% Verify LQR closed-loop stability
A_cl   = Ad + Bd*K_lqr;
rho_cl = max(abs(eig(A_cl)));
fprintf('\n--- Terminal ingredients ---\n');
fprintf('LQR closed-loop spectral radius: %.6f  (must be < 1)\n', rho_cl);
assert(rho_cl < 1, 'LQR is not stabilising — adjust Q, R.');

% Verify P_inf satisfies discrete Lyapunov / DARE condition
% A_cl'*P*A_cl - P + Q + K_lqr'*R*K_lqr = 0  (DARE residual)
DARE_residual = norm(A_cl'*P_inf*A_cl - P_inf + Q + K_lqr'*R*K_lqr, 'fro');
fprintf('DARE residual (Frobenius norm): %.2e  (should be ~0)\n', DARE_residual);

fprintf('\nComputing maximal LQR-invariant terminal set (MPT3) ... ');

%% =========================================================================
%%  7. CONSTRAINTS  (physical actuator and linearisation-validity limits)
%%
%%  State constraints (error coordinates):
%%    |ey|  <= 100 m   (lateral position error)
%%    |ez|  <= 100 m   (altitude deviation — decoupled, loose bound)
%%    |evy| <= 30 m/s  (lateral velocity)
%%    |evz| <= 30 m/s  (vertical velocity deviation)
%%    |theta| <= 15°   (linearisation validity; ≈0.262 rad)
%%    |q|   <= 20°/s   (pitch rate;  ≈0.349 rad/s)
%%
%%  Input constraints (hardware limits):
%%    |delta| <= 5°    (gimbal mechanical travel limit)
%%    |dT|    <= 150 N (~30% of T0 — fuel-flow range)
%% =========================================================================
e_lb = [-100; -100; -30; -30; -15*pi/180; -20*pi/180];
e_ub = [ 100;  100;  30;  30;  15*pi/180;  20*pi/180];

u_lb = [-5*pi/180;  -150];    % [delta_min;  dT_min]
u_ub = [ 5*pi/180;   150];    % [delta_max;  dT_max]

% Check K_lqr*e satisfies input constraints for e in the constraint set
% (needed for terminal set positive invariance)
u_at_corner  = abs(K_lqr * e_ub);   % rough check at one constraint-set corner (2x1)
fprintf('\n--- Constraint summary ---\n');
fprintf('|theta| <= %.1f deg,  |delta| <= %.1f deg\n', ...
        rad2deg(e_ub(5)), rad2deg(u_ub(1)));
fprintf('|K_lqr*e_ub| at corner: delta=%.4f rad  dT=%.1f N\n', ...
        u_at_corner(1), u_at_corner(2));
fprintf('  (terminal set X_f will be smaller than constraint set)\n');

%% =========================================================================
%%  7b. TERMINAL SET via MPT3
%% =========================================================================
tvc_sys              = LTISystem('A', Ad, 'B', Bd, 'Ts', Ts);
tvc_sys.x.min        = e_lb;
tvc_sys.x.max        = e_ub;
tvc_sys.u.min        = u_lb;
tvc_sys.u.max        = u_ub;
tvc_sys.x.penalty    = QuadFunction(Q);
tvc_sys.u.penalty    = QuadFunction(R);

tvc_sys.x.with('terminalPenalty');
tvc_sys.x.terminalPenalty = QuadFunction(P_inf);

tvc_sys.x.with('terminalSet');
X_f = tvc_sys.LQRSet;

fprintf('done.\n');
fprintf('Terminal set: Chebyshev radius = %.4f  (in error-state space)\n', ...
        X_f.chebyCenter.r);

% Extract H-rep for YALMIP constraints
Xf_A = X_f.A;
Xf_b = X_f.b;

%% =========================================================================
%%  8. INITIAL CONDITION
%%
%%  The rocket starts with:
%%    y  = 3 m lateral offset  (e.g. wind during ignition / pad misalignment)
%%    theta = 5°  (≈0.087 rad) initial tilt
%%    All velocities at nominal  →  vy=0, vz=vz_nom
%%
%%  In error coordinates (reference x_ref(0) = [0;0;0;vz_nom;0;0]):
%%    e0 = [3; 0; 0; 0; 0.087; 0]
%% =========================================================================
e0 = [3; 0; 0; 0; 5*pi/180; 0];

in_Xf = all(Xf_A * e0 <= Xf_b);
fprintf('\nInitial error e0 in X_f: %s  (N=%d should be sufficient)\n', ...
        mat2str(in_Xf), N);
assert(all(e0 >= e_lb) && all(e0 <= e_ub), 'e0 violates state constraints');

%% =========================================================================
%%  9. MPC CLOSED-LOOP SIMULATION  with live animation
%% =========================================================================
T_sim      = 150;                         % simulation steps  (7.5 s)
E_mpc      = zeros(nx, T_sim+1);
U_mpc      = zeros(nu, T_sim);
X_world    = zeros(nx, T_sim+1);          % world-frame states

% Initial world state: x = x_ref(0) + e0
x_ref_0      = [0; 0; 0; vz_nom; 0; 0];
X_world(:,1) = x_ref_0 + e0;
E_mpc(:,1)   = e0;

fprintf('\n--- MPC closed-loop simulation (N=%d, T_sim=%d) ---\n', N, T_sim);

% Figure layout: wide window
fig = figure('Color','w', 'Position',[40 60 1400 560], ...
             'Name','TVC MPC — Thrust Vector Control Rocket');

ax_rocket = subplot(1, 2, 1);
ax_states = subplot(1, 2, 2);

T_end_mpc = T_sim;
tic;

for k = 1:T_sim
    ek = E_mpc(:,k);

    % ---- Solve MPC QP ----
    % Terminal set constraint omitted: nonlinear model-plant mismatch can push
    % the state outside the linearised feasibility domain; DARE terminal cost
    % alone provides convergence (Rawlings & Mayne 2017, Remark 2.21).
    [u_k, U_seq, ok] = solve_mpc_tvc(ek, Ad, Bd, Q, R, P_inf, N, ...
                                      e_lb, e_ub, u_lb, u_ub, [], []);

    if ~ok
        warning('MPC infeasible at k=%d — switching to LQR fallback', k);
        u_k = max(u_lb, min(u_ub, K_lqr * ek));
        U_seq = repmat(u_k, 1, N);
    end

    U_mpc(:,k) = u_k;

    % ---- Propagate NONLINEAR plant one step ----
    x_real_k   = X_world(:,k);
    x_real_kp1 = tvc_nonlinear_step(x_real_k, u_k, m, g, J, l_tvc, T0, Ts);
    X_world(:,k+1) = x_real_kp1;

    % ---- Update error state (reference advances with time) ----
    x_ref_kp1     = [0; vz_nom*(k)*Ts; 0; vz_nom; 0; 0];
    E_mpc(:,k+1)  = x_real_kp1 - x_ref_kp1;

    % ---- Build predicted world-frame trajectory for animation ----
    pred_world = zeros(2, N+1);
    pred_world(:,1) = x_real_k(1:2);
    e_tmp = ek;
    for j = 1:N
        e_tmp = Ad*e_tmp + Bd*U_seq(:,j);
        z_ref_j = vz_nom * (k-1+j) * Ts;
        pred_world(1, j+1) = e_tmp(1);          % y_world = e_y  (y_ref = 0)
        pred_world(2, j+1) = e_tmp(2) + z_ref_j;% z_world = e_z + z_ref
    end

    % ---- Animation ----
    clf(fig);
    ax_rocket = subplot(1, 2, 1);
    ax_states = subplot(1, 2, 2);
    draw_tvc_frame(ax_rocket, ax_states, X_world(:,1:k+1), U_mpc(:,1:k), ...
                   pred_world, k*Ts, N, T_sim, Ts, rp);
    drawnow limitrate;
    pause(0.01);

    % ---- Convergence check ----
    if norm(E_mpc(:,k+1)) < 1e-3
        T_end_mpc = k;
        fprintf('  Converged at k=%d (t=%.2f s)\n', k, k*Ts);
        break;
    end
end

t_mpc_sim = toc;
fprintf('MPC simulation done in %.2f s  (wall clock)\n', t_mpc_sim);

E_mpc   = E_mpc(:,   1:T_end_mpc+1);
U_mpc   = U_mpc(:,   1:T_end_mpc);
X_world = X_world(:, 1:T_end_mpc+1);

fprintf('Final error norm: %.6f\n', norm(E_mpc(:,end)));
fprintf('Final pitch: %.4f deg\n', rad2deg(X_world(5,end)));

%% =========================================================================
%%  10. LQR BASELINE SIMULATION
%% =========================================================================
fprintf('\n--- LQR baseline ---\n');
T_end_lqr = T_sim;
E_lqr     = zeros(nx, T_sim+1);
U_lqr     = zeros(nu, T_sim);
X_lqr_w   = zeros(nx, T_sim+1);

X_lqr_w(:,1) = x_ref_0 + e0;
E_lqr(:,1)   = e0;

for k = 1:T_sim
    ek  = E_lqr(:,k);
    u_k = max(u_lb, min(u_ub, K_lqr * ek));   % LQR + hard saturation
    U_lqr(:,k) = u_k;

    x_next            = tvc_nonlinear_step(X_lqr_w(:,k), u_k, m, g, J, l_tvc, T0, Ts);
    X_lqr_w(:,k+1)   = x_next;
    x_ref_kp1         = [0; vz_nom*k*Ts; 0; vz_nom; 0; 0];
    E_lqr(:,k+1)      = x_next - x_ref_kp1;

    if norm(E_lqr(:,k+1)) < 1e-3
        T_end_lqr = k;
        break;
    end
end

E_lqr   = E_lqr(:, 1:T_end_lqr+1);
U_lqr   = U_lqr(:, 1:T_end_lqr);
X_lqr_w = X_lqr_w(:, 1:T_end_lqr+1);
fprintf('LQR converged at k=%d (t=%.2f s)\n', T_end_lqr, T_end_lqr*Ts);

%% =========================================================================
%%  11. HORIZON STUDY  (N = 5, 10, 20, 30)  — no animation, quick
%% =========================================================================
fprintf('\n--- Horizon study ---\n');
N_study    = [5, 10, 20, 30];
conv_times = zeros(1, numel(N_study));
peak_theta = zeros(1, numel(N_study));
solve_t    = zeros(1, numel(N_study));

for ni = 1:numel(N_study)
    Ni = N_study(ni);
    E_h = zeros(nx, T_sim+1);  E_h(:,1) = e0;
    X_h = zeros(nx, T_sim+1);  X_h(:,1) = x_ref_0 + e0;
    T_end_h = T_sim;

    t0h = tic;
    for k = 1:T_sim
        ek = E_h(:,k);
        % No terminal set constraint here — small N may not steer into X_f
        [u_k, ~, ok] = solve_mpc_tvc(ek, Ad, Bd, Q, R, P_inf, Ni, ...
                                      e_lb, e_ub, u_lb, u_ub, [], []);
        if ~ok; u_k = max(u_lb, min(u_ub, K_lqr*ek)); end
        x_next        = tvc_nonlinear_step(X_h(:,k), u_k, m, g, J, l_tvc, T0, Ts);
        X_h(:,k+1)    = x_next;
        x_ref_kp1     = [0; vz_nom*k*Ts; 0; vz_nom; 0; 0];
        E_h(:,k+1)    = x_next - x_ref_kp1;
        if norm(E_h(:,k+1)) < 1e-3; T_end_h = k; break; end
    end
    solve_t(ni)    = toc(t0h);
    conv_times(ni) = T_end_h * Ts;
    peak_theta(ni) = max(abs(rad2deg(E_h(5,1:T_end_h+1))));
    fprintf('  N=%2d: converged t=%.2f s  peak|theta|=%.2f deg  wall=%.2f s\n', ...
            Ni, conv_times(ni), peak_theta(ni), solve_t(ni));
end

% Horizon study plot
figure('Name','Horizon Study','NumberTitle','off','Color','w');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile; bar(N_study, conv_times, 0.5, 'FaceColor',[0.2 0.45 0.85]);
xlabel('Horizon N'); ylabel('Convergence time [s]');
title('Settling time vs N'); grid on;

nexttile; bar(N_study, peak_theta, 0.5, 'FaceColor',[0.85 0.2 0.2]);
xlabel('Horizon N'); ylabel('Peak |\theta|  [deg]');
title('Peak pitch excursion vs N'); grid on;

nexttile; bar(N_study, solve_t, 0.5, 'FaceColor',[0.2 0.72 0.2]);
xlabel('Horizon N'); ylabel('Total solve time [s]');
title('Computation time vs N'); grid on;

%% =========================================================================
%%  12. WEIGHT TUNING STUDY
%% =========================================================================
fprintf('\n--- Weight tuning study ---\n');

% Three tuning scenarios
W_cases = {
    diag([5,   0.5,  2,   0.1, 200, 20]),  diag([500, 0.1]),  'Nominal  (Q_\theta=200)';
    diag([5,   0.5,  2,   0.1,  50, 10]),  diag([500, 0.1]),  'Low \theta penalty  (Q_\theta=50)';
    diag([5,   0.5,  2,   0.1, 500, 50]),  diag([500, 0.1]),  'High \theta penalty  (Q_\theta=500)';
};

clrs_w = {[0.20 0.45 0.85]; [0.85 0.20 0.20]; [0.10 0.72 0.22]};

fig_wt = figure('Name','Weight Tuning','NumberTitle','off','Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ax_wy = nexttile; hold(ax_wy,'on'); grid(ax_wy,'on');
title(ax_wy, 'Lateral position  e_y');
ylabel(ax_wy, 'e_y [m]'); xlabel(ax_wy, 't [s]');

ax_wth = nexttile; hold(ax_wth,'on'); grid(ax_wth,'on');
title(ax_wth, 'Pitch angle  \theta');
ylabel(ax_wth, '\theta [deg]'); xlabel(ax_wth, 't [s]');

ax_wd = nexttile; hold(ax_wd,'on'); grid(ax_wd,'on');
title(ax_wd, 'Gimbal angle  \delta');
ylabel(ax_wd, '\delta [deg]'); xlabel(ax_wd, 't [s]');

ax_wn = nexttile; hold(ax_wn,'on'); grid(ax_wn,'on');
title(ax_wn, 'Error norm  \|e\|_2');
ylabel(ax_wn, '\|e\|_2  (log)'); xlabel(ax_wn, 't [s]');

for wi = 1:size(W_cases,1)
    Qi = W_cases{wi,1};  Ri = W_cases{wi,2};  lbl = W_cases{wi,3};
    [Pi, ~, Kdi] = dare(Ad, Bd, Qi, Ri);
    Ki = -Kdi;
    rho_i = max(abs(eig(Ad + Bd*Ki)));
    fprintf('  Case %d (%s): rho=%.4f\n', wi, lbl, rho_i);

    E_w = zeros(nx, T_sim+1);  E_w(:,1) = e0;
    U_w = zeros(nu, T_sim);
    X_w = zeros(nx, T_sim+1);  X_w(:,1) = x_ref_0 + e0;
    T_end_w = T_sim;
    for k = 1:T_sim
        ek = E_w(:,k);
        [u_k, ~, ok] = solve_mpc_tvc(ek, Ad, Bd, Qi, Ri, Pi, N, ...
                                      e_lb, e_ub, u_lb, u_ub, [], []);
        if ~ok; u_k = max(u_lb, min(u_ub, Ki*ek)); end
        U_w(:,k) = u_k;
        x_next        = tvc_nonlinear_step(X_w(:,k), u_k, m, g, J, l_tvc, T0, Ts);
        X_w(:,k+1)    = x_next;
        x_ref_kp1     = [0; vz_nom*k*Ts; 0; vz_nom; 0; 0];
        E_w(:,k+1)    = x_next - x_ref_kp1;
        if norm(E_w(:,k+1)) < 1e-3; T_end_w = k; break; end
    end
    E_w = E_w(:, 1:T_end_w+1);
    U_w = U_w(:, 1:T_end_w);
    t_w = (0:T_end_w)*Ts;

    plot(ax_wy,  t_w, E_w(1,:),                      'Color',clrs_w{wi},'LineWidth',1.8,'DisplayName',lbl);
    plot(ax_wth, t_w, rad2deg(E_w(5,:)),              'Color',clrs_w{wi},'LineWidth',1.8,'DisplayName',lbl);
    if ~isempty(U_w)
        stairs(ax_wd, (0:T_end_w-1)*Ts, rad2deg(U_w(1,:)), 'Color',clrs_w{wi},'LineWidth',1.8,'DisplayName',lbl);
    end
    semilogy(ax_wn, t_w, vecnorm(E_w),               'Color',clrs_w{wi},'LineWidth',1.8,'DisplayName',lbl);
end

% Add constraint lines
yline(ax_wy,  e_ub(1),'k--','LineWidth',0.8);   yline(ax_wy,  e_lb(1),'k--','LineWidth',0.8);
yline(ax_wth, rad2deg(e_ub(5)),'k--','LineWidth',0.8); yline(ax_wth, rad2deg(e_lb(5)),'k--','LineWidth',0.8);
yline(ax_wd,  rad2deg(u_ub(1)),'k--','LineWidth',0.8); yline(ax_wd,  rad2deg(u_lb(1)),'k--','LineWidth',0.8);

legend(ax_wy,'Location','northeast','FontSize',7);
legend(ax_wth,'Location','northeast','FontSize',7);
legend(ax_wd,'Location','northeast','FontSize',7);
legend(ax_wn,'Location','northeast','FontSize',7);

%% =========================================================================
%%  13. DISTURBANCE REJECTION TEST  (lateral wind gust)
%% =========================================================================
fprintf('\n--- Disturbance rejection (wind gust at t=1 s) ---\n');
gust_k      = round(1.0 / Ts);   % step index for gust
gust_force  = 200;                % lateral gust force [N]  (= 0.4*T0)
gust_dur    = round(0.25 / Ts);  % gust duration [steps]

E_dist = zeros(nx, T_sim+1);  E_dist(:,1) = e0;
X_dist = zeros(nx, T_sim+1);  X_dist(:,1) = x_ref_0 + e0;
U_dist = zeros(nu, T_sim);
T_end_dist = T_sim;

for k = 1:T_sim
    ek = E_dist(:,k);
    [u_k, ~, ok] = solve_mpc_tvc(ek, Ad, Bd, Q, R, P_inf, N, ...
                                  e_lb, e_ub, u_lb, u_ub, [], []);
    if ~ok; u_k = max(u_lb, min(u_ub, K_lqr*ek)); end
    U_dist(:,k) = u_k;

    % Nonlinear plant step with optional wind gust disturbance
    % Gust modelled as additional lateral force: add delta_f to vy_dot
    x_k = X_dist(:,k);
    x_next = tvc_nonlinear_step(x_k, u_k, m, g, J, l_tvc, T0, Ts);
    if k >= gust_k && k < gust_k + gust_dur
        % Impulsive lateral force: approximate as delta_vy = F/m * Ts
        x_next(3) = x_next(3) + gust_force/m * Ts;
    end
    X_dist(:,k+1)  = x_next;
    x_ref_kp1      = [0; vz_nom*k*Ts; 0; vz_nom; 0; 0];
    E_dist(:,k+1)  = x_next - x_ref_kp1;
    if norm(E_dist(:,k+1)) < 5e-3; T_end_dist = k; break; end
end

E_dist = E_dist(:, 1:T_end_dist+1);
U_dist = U_dist(:, 1:T_end_dist);

figure('Name','Disturbance Rejection','NumberTitle','off','Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile; hold on; grid on;
plot((0:T_end_mpc)*Ts,  E_mpc(1,:),   '-',  'Color',[0.2 0.45 0.85],'LineWidth',1.8,'DisplayName','MPC (no disturbance)');
plot((0:T_end_dist)*Ts, E_dist(1,:),  '--', 'Color',[0.85 0.5 0.2],'LineWidth',1.8,'DisplayName','MPC + wind gust');
xline(gust_k*Ts, 'k:', 'LineWidth',1.5, 'DisplayName','Gust onset');
yline(e_ub(1),'k--','LineWidth',0.8,'HandleVisibility','off');
ylabel('e_y [m]'); xlabel('t [s]');
title('Lateral error — disturbance rejection (200 N lateral gust for 0.25 s)');
legend('Location','northeast');

nexttile; hold on; grid on;
plot((0:T_end_mpc)*Ts,  rad2deg(E_mpc(5,:)),   '-',  'Color',[0.2 0.45 0.85],'LineWidth',1.8);
plot((0:T_end_dist)*Ts, rad2deg(E_dist(5,:)),  '--', 'Color',[0.85 0.5 0.2],'LineWidth',1.8);
xline(gust_k*Ts, 'k:', 'LineWidth',1.5);
yline(rad2deg(e_ub(5)),'k--','LineWidth',0.8);
yline(rad2deg(e_lb(5)),'k--','LineWidth',0.8);
ylabel('\theta [deg]'); xlabel('t [s]');
title('Pitch angle during gust');

%% =========================================================================
%%  14. OBSERVER DESIGN + MEASUREMENT NOISE SIMULATION
%%
%%  Kalman filter (discrete-time LQE, dual of LQR) estimates the error
%%  state from noisy sensor readings. The MPC then closes the loop on the
%%  estimate, not the true state.
%%
%%  Sensor noise model (1-sigma):
%%    y, z   : ±0.50 m     (GPS horizontal/vertical)
%%    vy, vz : ±0.32 m/s   (inertial + Doppler fusion)
%%    theta  : ±1°          (IMU pitch angle)
%%    q      : ±2°/s        (rate gyroscope)
%% =========================================================================
fprintf('\n--- Observer design (Kalman filter / discrete LQE) ---\n');

Cd = eye(nx);   % full-state measurement (all 6 states observed with noise)

% Process noise covariance Qn — model uncertainty (unmodelled dynamics, gusts)
Qn = diag([0.01; 0.01; 0.05; 0.05; (0.5*pi/180)^2; (1*pi/180)^2]);

% Measurement noise covariance Rn — sensor accuracy (1-sigma values above)
sigma_y    = 0.50;          % GPS position noise [m]
sigma_v    = 0.32;          % velocity noise [m/s]
sigma_th   = 1*pi/180;      % IMU pitch noise [rad]
sigma_q    = 2*pi/180;      % gyro rate noise [rad/s]
Rn = diag([sigma_y^2; sigma_y^2; sigma_v^2; sigma_v^2; sigma_th^2; sigma_q^2]);

% Steady-state Kalman gain via dlqe (dual of LQR / DARE on transposed system)
[L_kf, ~] = dlqe(Ad, eye(nx), Cd, Qn, Rn);

% Observer error dynamics: eigenvalues of (I - L_kf*Cd)*Ad (corrector form)
eigs_obs = eig((eye(nx) - L_kf*Cd)*Ad);
fprintf('Kalman gain max singular value:    %.4f\n', max(svd(L_kf)));
fprintf('Observer error max|eig|:           %.6f  (must be < 1)\n', max(abs(eigs_obs)));

% --- Observer-based MPC simulation with measurement noise ---
rng(42);   % reproducible random seed
T_end_obs = T_sim;

E_obs  = zeros(nx, T_sim+1);  E_obs(:,1)  = e0;   % TRUE error state
E_hat  = zeros(nx, T_sim+1);  E_hat(:,1)  = e0;   % ESTIMATED error state
X_obs  = zeros(nx, T_sim+1);  X_obs(:,1)  = x_ref_0 + e0;
U_obs  = zeros(nu, T_sim);
e_hat  = e0;   % initial estimate = true state (known at ignition)

for k = 1:T_sim
    % ---- MPC uses ESTIMATED state ----
    [u_k, ~, ok] = solve_mpc_tvc(e_hat, Ad, Bd, Q, R, P_inf, N, ...
                                  e_lb, e_ub, u_lb, u_ub, [], []);
    if ~ok; u_k = max(u_lb, min(u_ub, K_lqr*e_hat)); end
    U_obs(:,k) = u_k;

    % ---- Propagate TRUE nonlinear plant ----
    x_next = tvc_nonlinear_step(X_obs(:,k), u_k, m, g, J, l_tvc, T0, Ts);
    X_obs(:,k+1) = x_next;
    x_ref_kp1   = [0; vz_nom*k*Ts; 0; vz_nom; 0; 0];
    E_obs(:,k+1) = x_next - x_ref_kp1;

    % ---- Noisy measurement of new error state ----
    noise_v = sqrt(diag(Rn)) .* randn(nx,1);
    y_meas  = E_obs(:,k+1) + noise_v;

    % ---- Kalman filter: time update (predict) then measurement update ----
    e_hat_pred = Ad*e_hat + Bd*u_k;                          % prediction
    e_hat      = e_hat_pred + L_kf*(y_meas - Cd*e_hat_pred); % correction
    E_hat(:,k+1) = e_hat;

    if norm(E_obs(:,k+1)) < 1e-3; T_end_obs = k; break; end
end

E_obs = E_obs(:, 1:T_end_obs+1);
E_hat = E_hat(:, 1:T_end_obs+1);
U_obs = U_obs(:, 1:T_end_obs);
fprintf('Observer-based MPC converged at t = %.2f s\n', T_end_obs*Ts);

% Quick plot
figure('Name','Observer MPC with Measurement Noise','NumberTitle','off','Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile; hold on; grid on;
plot((0:T_end_mpc)*Ts,  E_mpc(1,:),   '-',  'Color',[0.20 0.45 0.85],'LineWidth',1.8,'DisplayName','MPC (no noise)');
plot((0:T_end_obs)*Ts,  E_obs(1,:),   '-',  'Color',[0.10 0.72 0.22],'LineWidth',1.8,'DisplayName','True state (noisy)');
plot((0:T_end_obs)*Ts,  E_hat(1,:),   '--', 'Color',[0.92 0.45 0.10],'LineWidth',1.2,'DisplayName','KF estimate');
yline(e_ub(1),'k--','LineWidth',0.8,'HandleVisibility','off');
yline(e_lb(1),'k--','LineWidth',0.8,'HandleVisibility','off');
ylabel('e_y [m]'); xlabel('t [s]');
title('Lateral error — observer-based MPC with measurement noise');
legend('Location','northeast');

nexttile; hold on; grid on;
plot((0:T_end_mpc)*Ts,  rad2deg(E_mpc(5,:)),  '-',  'Color',[0.20 0.45 0.85],'LineWidth',1.8,'DisplayName','MPC (no noise)');
plot((0:T_end_obs)*Ts,  rad2deg(E_obs(5,:)),  '-',  'Color',[0.10 0.72 0.22],'LineWidth',1.8,'DisplayName','True state (noisy)');
plot((0:T_end_obs)*Ts,  rad2deg(E_hat(5,:)),  '--', 'Color',[0.92 0.45 0.10],'LineWidth',1.2,'DisplayName','KF estimate');
yline(rad2deg(e_ub(5)),'k--','LineWidth',0.8,'HandleVisibility','off');
yline(rad2deg(e_lb(5)),'k--','LineWidth',0.8,'HandleVisibility','off');
ylabel('\theta [deg]'); xlabel('t [s]');
title('Pitch angle — observer tracking');

%% =========================================================================
%%  16. FINAL SUMMARY PLOTS
%% =========================================================================
fprintf('\n--- Generating final summary plots ---\n');

t_mpc_v = (0:T_end_mpc)  * Ts;
t_lqr_v = (0:T_end_lqr)  * Ts;

plot_results_tvc(t_mpc_v, E_mpc, U_mpc, t_lqr_v, E_lqr, U_lqr, ...
                 e_lb, e_ub, u_lb, u_ub, vz_nom);

%% =========================================================================
%%  17. FINAL REPORT
%% =========================================================================
fprintf('\n========== FINAL REPORT ==========\n');
fprintf('System:  m=%g kg  T0=%.1f N  J=%g kg.m^2  l_tvc=%g m\n', m, T0, J, l_tvc);
fprintf('LQR spectral radius:         %.6f\n', rho_cl);
fprintf('DARE residual:               %.2e\n', DARE_residual);
fprintf('Terminal set Chebyshev r:    %.4f\n', X_f.chebyCenter.r);
fprintf('Initial error in X_f:        %s\n',   mat2str(in_Xf));
fprintf('\nMPC  convergence: t = %.2f s   ||e_final|| = %.6f\n', ...
        T_end_mpc*Ts, norm(E_mpc(:,end)));
fprintf('LQR  convergence: t = %.2f s   ||e_final|| = %.6f\n', ...
        T_end_lqr*Ts, norm(E_lqr(:,end)));
fprintf('\nHorizon study:  N = %s\n', mat2str(N_study));
fprintf('  Settling:     %s s\n', mat2str(conv_times, 3));
fprintf('  Peak theta:   %s deg\n', mat2str(peak_theta, 3));
fprintf('  Solve time:   %s s\n', mat2str(solve_t, 3));
fprintf('=====================================\n');
