%% save_figures_tvc.m — Generate and save all report figures
%%
%% Run AFTER main_tvc_mpc.m so workspace variables are available.
%% Saves publication-quality PDF figures to  report/figures/.
%%
%% Figures produced:
%%   fig_mpc_vs_lqr.pdf   — MPC vs LQR comparison (3-panel)
%%   fig_animation.pdf    — Animation snapshot at t = 1 s (2-panel)
%%   fig_horizon.pdf      — Horizon study bar charts (3-panel)
%%   fig_weights.pdf      — Weight-tuning comparison (4-panel)
%%   fig_disturbance.pdf  — Disturbance rejection (2-panel)

%% Check workspace
required = {'E_mpc','U_mpc','E_lqr','U_lqr','X_world','Ts','e_lb','e_ub', ...
            'u_lb','u_ub','N_study','conv_times','peak_theta','solve_t', ...
            'Ad','Bd','Q','R','P_inf','K_lqr','m','g','J','l_tvc','T0', ...
            'vz_nom','N','rp','Xf_A','Xf_b', ...
            'E_obs','E_hat','L_kf','Rn','T_end_mpc','T_end_obs'};
for i = 1:numel(required)
    if ~exist(required{i}, 'var')
        error('Variable ''%s'' not found. Run main_tvc_mpc.m first.', required{i});
    end
end

FIG_DIR = fullfile(fileparts(mfilename('fullpath')), 'report', 'figures');
if ~exist(FIG_DIR, 'dir'); mkdir(FIG_DIR); end
fprintf('Saving figures to:  %s\n', FIG_DIR);

%% Shared style helpers
clr_mpc  = [0.20 0.45 0.85];
clr_lqr  = [0.85 0.20 0.20];
clr_safe = [0.85 1.00 0.85];
clr_horiz= [0.92 0.45 0.10];

LW   = 2.0;   % main line width
LWc  = 1.0;   % constraint line width
FS   = 14;    % axis font size
FStitle = 15;

% Use vector PDF export (requires R2020a+; fallback: print -dpdf)
savepdf = @(fig, name) exportgraphics(fig, fullfile(FIG_DIR, name), ...
    'ContentType','vector','BackgroundColor','white');

%% Helper: consistent figure + layout
function h = make_fig(w_cm, h_cm)
    h = figure('Color','w','Units','centimeters');
    h.Position = [2 2 w_cm h_cm];
end

% =========================================================================
%% FIGURE 1 — MPC vs LQR comparison
% =========================================================================
T_mpc = (0:size(E_mpc,2)-1)*Ts;
T_lqr = (0:size(E_lqr,2)-1)*Ts;
tmax  = max(T_mpc(end), T_lqr(end));

fig1 = make_fig(16, 18);
tiledlayout(fig1, 3, 1, 'TileSpacing','compact','Padding','compact');

% Panel 1: lateral error
ax = nexttile; hold on; grid on;
patch([0 tmax tmax 0],[e_lb(1) e_lb(1) e_ub(1) e_ub(1)], clr_safe, ...
      'EdgeColor','none','DisplayName','Feasible region');
plot(ax, [0 tmax],[e_ub(1) e_ub(1)],'k--','LineWidth',LWc,'HandleVisibility','off');
plot(ax, [0 tmax],[e_lb(1) e_lb(1)],'k--','LineWidth',LWc,'HandleVisibility','off');
plot(ax, T_mpc, E_mpc(1,:),  '-',  'Color',clr_mpc,'LineWidth',LW,'DisplayName','MPC');
plot(ax, T_lqr, E_lqr(1,:),  '--', 'Color',clr_lqr,'LineWidth',LW,'DisplayName','LQR');
yline(0,'k:','LineWidth',0.8,'HandleVisibility','off');
xlabel('$t$ [s]','Interpreter','latex');
ylabel('$e_y$ [m]','Interpreter','latex');
title('Lateral position error','FontSize',FStitle);
legend('Location','northeast','FontSize',FS-1);
xlim([0 tmax]); set(ax,'FontSize',FS);

% Panel 2: pitch angle
ax = nexttile; hold on; grid on;
patch([0 tmax tmax 0],[rad2deg(e_lb(5)) rad2deg(e_lb(5)) ...
      rad2deg(e_ub(5)) rad2deg(e_ub(5))], clr_safe, 'EdgeColor','none');
plot(ax, [0 tmax],rad2deg(e_ub(5))*[1 1],'k--','LineWidth',LWc,'DisplayName','Constraint $\pm15^{\circ}$');
plot(ax, [0 tmax],rad2deg(e_lb(5))*[1 1],'k--','LineWidth',LWc,'HandleVisibility','off');
plot(ax, T_mpc, rad2deg(E_mpc(5,:)),  '-',  'Color',clr_mpc,'LineWidth',LW,'DisplayName','MPC');
plot(ax, T_lqr, rad2deg(E_lqr(5,:)),  '--', 'Color',clr_lqr,'LineWidth',LW,'DisplayName','LQR');
yline(0,'k:','LineWidth',0.8,'HandleVisibility','off');
xlabel('$t$ [s]','Interpreter','latex');
ylabel('$\theta$ [deg]','Interpreter','latex');
title('Pitch angle $e_\theta$','Interpreter','latex','FontSize',FStitle);
legend('Location','northeast','FontSize',FS-1);
xlim([0 tmax]); set(ax,'FontSize',FS);

% Panel 3: gimbal angle
ax = nexttile; hold on; grid on;
patch([0 tmax tmax 0],[rad2deg(u_lb(1)) rad2deg(u_lb(1)) ...
      rad2deg(u_ub(1)) rad2deg(u_ub(1))], clr_safe, 'EdgeColor','none');
plot(ax, [0 tmax],rad2deg(u_ub(1))*[1 1],'k--','LineWidth',LWc,'DisplayName','Constraint $\pm5^{\circ}$');
plot(ax, [0 tmax],rad2deg(u_lb(1))*[1 1],'k--','LineWidth',LWc,'HandleVisibility','off');
if size(U_mpc,2)>0
    stairs(ax, T_mpc(1:end-1), rad2deg(U_mpc(1,:)), '-',  'Color',clr_mpc,'LineWidth',LW,'DisplayName','MPC');
end
if size(U_lqr,2)>0
    stairs(ax, T_lqr(1:end-1), rad2deg(U_lqr(1,:)), '--', 'Color',clr_lqr,'LineWidth',LW,'DisplayName','LQR');
end
yline(0,'k:','LineWidth',0.8,'HandleVisibility','off');
xlabel('$t$ [s]','Interpreter','latex');
ylabel('$\delta$ [deg]','Interpreter','latex');
title('Gimbal angle $\delta$','Interpreter','latex','FontSize',FStitle);
legend('Location','northeast','FontSize',FS-1);
xlim([0 tmax]); set(ax,'FontSize',FS);

savepdf(fig1, 'fig_mpc_vs_lqr.pdf');
fprintf('  Saved fig_mpc_vs_lqr.pdf\n');

% =========================================================================
%% FIGURE 2 — Animation snapshot at t = 1 s (k = 20)
% =========================================================================
k_snap = min(round(1.0/Ts) + 1, size(X_world,2));  % column index at t ≈ 1 s

fig2 = make_fig(16, 10);
ax_r = subplot(1,2,1);
ax_s = subplot(1,2,2);

% Reconstruct predicted trajectory at k_snap using linearised model
ek_snap = E_mpc(:, k_snap);
U_snap  = zeros(2, N);   % approximate: use zeros (actual U_seq not stored)
% Better: replay the solver to get the horizon; fall back to straight prediction
pred_w = zeros(2, N+1);
pred_w(:,1) = X_world(1:2, k_snap);
e_tmp = ek_snap;
for j = 1:N
    u_j = max(u_lb, min(u_ub, K_lqr * e_tmp));   % LQR prediction for horizon visualisation
    e_tmp = Ad*e_tmp + Bd*u_j;
    z_ref_j = vz_nom * (k_snap-1+j) * Ts;
    pred_w(1,j+1) = e_tmp(1);
    pred_w(2,j+1) = e_tmp(2) + z_ref_j;
end

draw_tvc_frame(ax_r, ax_s, X_world(:, 1:k_snap), U_mpc(:, 1:k_snap-1), ...
               pred_w, (k_snap-1)*Ts, N, size(X_world,2)-1, Ts, rp);

sgtitle(fig2, sprintf('Animation snapshot  |  t = %.1f s', (k_snap-1)*Ts), ...
        'FontSize', FStitle+2, 'FontWeight','bold');

savepdf(fig2, 'fig_animation.pdf');
fprintf('  Saved fig_animation.pdf\n');

% =========================================================================
%% FIGURE 3 — Horizon study
% =========================================================================
fig3 = make_fig(16, 9);
tiledlayout(fig3, 1, 3, 'TileSpacing','compact','Padding','compact');

nexttile;
bar(N_study, conv_times, 0.55, 'FaceColor',clr_mpc, 'EdgeColor','k');
xlabel('Horizon $N$','Interpreter','latex');
ylabel('Settling time [s]');
title('Settling time vs $N$','Interpreter','latex','FontSize',FStitle);
grid on; set(gca,'FontSize',FS);

nexttile;
bar(N_study, peak_theta, 0.55, 'FaceColor',clr_lqr, 'EdgeColor','k');
xlabel('Horizon $N$','Interpreter','latex');
ylabel('Peak $|\theta|$ [deg]','Interpreter','latex');
title('Peak pitch vs $N$','Interpreter','latex','FontSize',FStitle);
grid on; set(gca,'FontSize',FS);

nexttile;
bar(N_study, solve_t, 0.55, 'FaceColor',[0.20 0.72 0.20], 'EdgeColor','k');
xlabel('Horizon $N$','Interpreter','latex');
ylabel('Computation time [s]');
title('Solve time vs $N$','Interpreter','latex','FontSize',FStitle);
grid on; set(gca,'FontSize',FS);

savepdf(fig3, 'fig_horizon.pdf');
fprintf('  Saved fig_horizon.pdf\n');

% =========================================================================
%% FIGURE 4 — Weight tuning  (re-run 3 quick simulations)
% =========================================================================
W_cases = {
    diag([5, 0.5, 2, 0.1,  50, 10]), diag([500, 0.1]), 'Low  $Q_\theta=50$';
    diag([5, 0.5, 2, 0.1, 200, 20]), diag([500, 0.1]), 'Nominal  $Q_\theta=200$';
    diag([5, 0.5, 2, 0.1, 500, 50]), diag([500, 0.1]), 'High  $Q_\theta=500$';
};
clrs_w = {[0.85 0.20 0.20]; [0.20 0.45 0.85]; [0.10 0.72 0.22]};
lstyles = {'--', '-', ':'};

nx = 6; nu = 2;
e0 = E_mpc(:,1);   % same initial condition as main simulation
x_ref_0 = [0; 0; 0; vz_nom; 0; 0];
T_sim_w = round(min(size(E_mpc,2), 200)) - 1;

fig4 = make_fig(16, 18);
tiledlayout(fig4, 2, 2, 'TileSpacing','compact','Padding','compact');
ax4 = gobjects(4,1);
for i = 1:4; ax4(i) = nexttile; hold(ax4(i),'on'); grid(ax4(i),'on'); end

titles4 = {'Lateral error  $e_y$ [m]', 'Pitch angle  $\theta$ [deg]', ...
           'Gimbal angle  $\delta$ [deg]', 'Error norm  $\|e\|_2$'};
for i = 1:4; title(ax4(i), titles4{i}, 'Interpreter','latex','FontSize',FStitle); end

% Constraint lines on each panel
for i = 1:4
    if i == 4; set(ax4(i),'YScale','log'); end
end
plot(ax4(1), [0 T_sim_w*Ts], [1 1]*e_ub(1),     'k--','LineWidth',LWc,'HandleVisibility','off');
plot(ax4(1), [0 T_sim_w*Ts], [1 1]*e_lb(1),     'k--','LineWidth',LWc,'HandleVisibility','off');
plot(ax4(2), [0 T_sim_w*Ts], [1 1]*rad2deg(e_ub(5)),'k--','LineWidth',LWc,'HandleVisibility','off');
plot(ax4(2), [0 T_sim_w*Ts], [1 1]*rad2deg(e_lb(5)),'k--','LineWidth',LWc,'HandleVisibility','off');
plot(ax4(3), [0 T_sim_w*Ts], [1 1]*rad2deg(u_ub(1)),'k--','LineWidth',LWc,'DisplayName','Constraint','HandleVisibility','off');
plot(ax4(3), [0 T_sim_w*Ts], [1 1]*rad2deg(u_lb(1)),'k--','LineWidth',LWc,'HandleVisibility','off');

for wi = 1:3
    Qi = W_cases{wi,1};  Ri = W_cases{wi,2};  lbl = W_cases{wi,3};
    [Pi, ~, Kdi] = dare(Ad, Bd, Qi, Ri);
    Ki = -Kdi;

    E_w = zeros(nx, T_sim_w+1);  E_w(:,1) = e0;
    U_w = zeros(nu, T_sim_w);
    X_w = zeros(nx, T_sim_w+1);  X_w(:,1) = x_ref_0 + e0;
    T_end_w = T_sim_w;
    for k = 1:T_sim_w
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
    E_w = E_w(:,1:T_end_w+1);
    U_w = U_w(:,1:T_end_w);
    t_w = (0:T_end_w)*Ts;

    ls = lstyles{wi};  clr = clrs_w{wi};
    plot(ax4(1), t_w, E_w(1,:),             ls,'Color',clr,'LineWidth',LW,'DisplayName',lbl);
    plot(ax4(2), t_w, rad2deg(E_w(5,:)),    ls,'Color',clr,'LineWidth',LW,'DisplayName',lbl);
    if ~isempty(U_w)
        stairs(ax4(3),(0:T_end_w-1)*Ts, rad2deg(U_w(1,:)), ls,'Color',clr,'LineWidth',LW,'DisplayName',lbl);
    end
    semilogy(ax4(4), t_w, max(vecnorm(E_w),1e-8), ls,'Color',clr,'LineWidth',LW,'DisplayName',lbl);
    fprintf('  Weight case %d (Q_th=%g): settled at t=%.2f s\n', wi, Qi(5,5), T_end_w*Ts);
end

for i = 1:4
    xlabel(ax4(i),'$t$ [s]','Interpreter','latex');
    xlim(ax4(i),[0 T_sim_w*Ts]);
    set(ax4(i),'FontSize',FS);
    legend(ax4(i),'Location','northeast','FontSize',FS-1,'Interpreter','latex');
end
ylabel(ax4(1),'$e_y$ [m]','Interpreter','latex');
ylabel(ax4(2),'$\theta$ [deg]');
ylabel(ax4(3),'$\delta$ [deg]');
ylabel(ax4(4),'$\|e\|_2$','Interpreter','latex');

savepdf(fig4, 'fig_weights.pdf');
fprintf('  Saved fig_weights.pdf\n');

% =========================================================================
%% FIGURE 5 — Disturbance rejection
% =========================================================================
% Re-run disturbance simulation with same parameters as main script
gust_k    = round(1.0 / Ts);
gust_F    = 200;
gust_dur  = round(0.25 / Ts);
T_sim_d   = T_sim_w;

E_dist = zeros(nx, T_sim_d+1);  E_dist(:,1) = e0;
X_dist = zeros(nx, T_sim_d+1);  X_dist(:,1) = x_ref_0 + e0;
T_end_d = T_sim_d;

for k = 1:T_sim_d
    ek = E_dist(:,k);
    [u_k, ~, ok] = solve_mpc_tvc(ek, Ad, Bd, Q, R, P_inf, N, ...
                                  e_lb, e_ub, u_lb, u_ub, [], []);
    if ~ok; u_k = max(u_lb, min(u_ub, K_lqr*ek)); end
    x_next = tvc_nonlinear_step(X_dist(:,k), u_k, m, g, J, l_tvc, T0, Ts);
    if k >= gust_k && k < gust_k + gust_dur
        x_next(3) = x_next(3) + gust_F/m * Ts;
    end
    X_dist(:,k+1)  = x_next;
    x_ref_kp1      = [0; vz_nom*k*Ts; 0; vz_nom; 0; 0];
    E_dist(:,k+1)  = x_next - x_ref_kp1;
    if norm(E_dist(:,k+1)) < 5e-3; T_end_d = k; break; end
end
E_dist = E_dist(:, 1:T_end_d+1);
T_mpc_d = (0:size(E_mpc,2)-1)*Ts;
T_dist  = (0:T_end_d)*Ts;

fig5 = make_fig(16, 13);
tiledlayout(fig5, 2, 1, 'TileSpacing','compact','Padding','compact');

ax = nexttile; hold on; grid on;
patch([0 max(T_dist(end),T_mpc(end)) max(T_dist(end),T_mpc(end)) 0], ...
      [e_lb(1) e_lb(1) e_ub(1) e_ub(1)], clr_safe, 'EdgeColor','none');
plot([0 max(T_dist(end),T_mpc(end))],[e_ub(1) e_ub(1)],'k--','LineWidth',LWc,'HandleVisibility','off');
plot(T_mpc_d, E_mpc(1,:), '-',  'Color',clr_mpc,'LineWidth',LW,'DisplayName','No disturbance');
plot(T_dist,  E_dist(1,:),'--', 'Color',clr_horiz,'LineWidth',LW,'DisplayName','Wind gust  (200 N, 0.25 s)');
xline(gust_k*Ts,'k:','LineWidth',1.2,'DisplayName','Gust onset');
yline(0,'k:','LineWidth',0.8,'HandleVisibility','off');
xlabel('$t$ [s]','Interpreter','latex');
ylabel('$e_y$ [m]','Interpreter','latex');
title('Lateral error — disturbance rejection','FontSize',FStitle);
legend('Location','northeast','FontSize',FS-1);
xlim([0 max(T_dist(end), T_mpc(end))]); set(ax,'FontSize',FS);

ax = nexttile; hold on; grid on;
plot([0 max(T_dist(end),T_mpc(end))],rad2deg(e_ub(5))*[1 1],'k--','LineWidth',LWc,'HandleVisibility','off');
plot([0 max(T_dist(end),T_mpc(end))],rad2deg(e_lb(5))*[1 1],'k--','LineWidth',LWc,'HandleVisibility','off');
plot(T_mpc_d, rad2deg(E_mpc(5,:)), '-',  'Color',clr_mpc,'LineWidth',LW,'DisplayName','No disturbance');
plot(T_dist,  rad2deg(E_dist(5,:)),'--', 'Color',clr_horiz,'LineWidth',LW,'DisplayName','Wind gust');
xline(gust_k*Ts,'k:','LineWidth',1.2,'HandleVisibility','off');
yline(0,'k:','LineWidth',0.8,'HandleVisibility','off');
xlabel('$t$ [s]','Interpreter','latex');
ylabel('$\theta$ [deg]');
title('Pitch angle $e_\theta$','Interpreter','latex','FontSize',FStitle);
legend('Location','northeast','FontSize',FS-1);
xlim([0 max(T_dist(end), T_mpc(end))]); set(ax,'FontSize',FS);

savepdf(fig5, 'fig_disturbance.pdf');
fprintf('  Saved fig_disturbance.pdf\n');

% =========================================================================
%% FIGURE 6 — Observer (Kalman filter) with measurement noise
% =========================================================================
% 1-sigma noise levels (from main_tvc_mpc.m Rn diagonal)
sigma_str = sprintf('GPS $\\pm%.2f$ m, IMU $\\pm%.0f^{\\circ}$', ...
                    sqrt(Rn(1,1)), rad2deg(sqrt(Rn(5,5))));

T_kf = (0:T_end_obs)*Ts;
T_cl = (0:T_end_mpc)*Ts;

fig6 = make_fig(16, 15);
tiledlayout(fig6, 2, 1, 'TileSpacing','compact','Padding','compact');

% Panel 1: lateral error
ax = nexttile; hold on; grid on;
patch([0 max(T_kf(end),T_cl(end)) max(T_kf(end),T_cl(end)) 0], ...
      [e_lb(1) e_lb(1) e_ub(1) e_ub(1)], clr_safe, 'EdgeColor','none','HandleVisibility','off');
plot(ax,[0 max(T_kf(end),T_cl(end))],[e_ub(1) e_ub(1)],'k--','LineWidth',LWc,'HandleVisibility','off');
plot(ax,[0 max(T_kf(end),T_cl(end))],[e_lb(1) e_lb(1)],'k--','LineWidth',LWc,'HandleVisibility','off');
plot(ax, T_cl, E_mpc(1,:),  '-',  'Color',clr_mpc, 'LineWidth',LW, 'DisplayName','MPC (no noise)');
plot(ax, T_kf, E_obs(1,:),  '-',  'Color',[0.10 0.68 0.22], 'LineWidth',LW, 'DisplayName','True state (noisy plant)');
plot(ax, T_kf, E_hat(1,:),  '--', 'Color',clr_horiz, 'LineWidth',LW-0.4, 'DisplayName','KF estimate');
yline(0,'k:','LineWidth',0.8,'HandleVisibility','off');
xlabel('$t$ [s]','Interpreter','latex');
ylabel('$e_y$ [m]','Interpreter','latex');
title(['Lateral error — observer-based MPC  (' sigma_str ')'], ...
      'Interpreter','latex','FontSize',FStitle);
legend('Location','northeast','FontSize',FS-1,'Interpreter','latex');
xlim([0 max(T_kf(end),T_cl(end))]); set(ax,'FontSize',FS);

% Panel 2: pitch angle
ax = nexttile; hold on; grid on;
plot(ax,[0 max(T_kf(end),T_cl(end))],rad2deg(e_ub(5))*[1 1],'k--','LineWidth',LWc,'HandleVisibility','off');
plot(ax,[0 max(T_kf(end),T_cl(end))],rad2deg(e_lb(5))*[1 1],'k--','LineWidth',LWc,'HandleVisibility','off');
plot(ax, T_cl, rad2deg(E_mpc(5,:)),  '-',  'Color',clr_mpc,         'LineWidth',LW,     'DisplayName','MPC (no noise)');
plot(ax, T_kf, rad2deg(E_obs(5,:)),  '-',  'Color',[0.10 0.68 0.22],'LineWidth',LW,     'DisplayName','True state (noisy plant)');
plot(ax, T_kf, rad2deg(E_hat(5,:)),  '--', 'Color',clr_horiz,       'LineWidth',LW-0.4, 'DisplayName','KF estimate');
yline(0,'k:','LineWidth',0.8,'HandleVisibility','off');
xlabel('$t$ [s]','Interpreter','latex');
ylabel('$\theta$ [deg]','Interpreter','latex');
title('Pitch angle — KF tracking','FontSize',FStitle);
legend('Location','northeast','FontSize',FS-1,'Interpreter','latex');
xlim([0 max(T_kf(end),T_cl(end))]); set(ax,'FontSize',FS);

savepdf(fig6, 'fig_observer.pdf');
fprintf('  Saved fig_observer.pdf\n');

% =========================================================================
fprintf('\nAll figures saved to:  %s\n', FIG_DIR);
