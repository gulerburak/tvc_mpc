function x_next = tvc_nonlinear_step(x, u, m, g, J, l_tvc, T0, Ts)
% TVC_NONLINEAR_STEP  RK4 integration of tvc_ode over one sample period Ts
%
% Integrates the nonlinear 2-D rocket equations one step forward using the
% classic 4th-order Runge-Kutta method.  Used as the plant model during
% closed-loop simulation (replaces the linearised prediction model).
%
% Inputs / outputs match tvc_ode conventions.

k1 = tvc_ode(x,            u, m, g, J, l_tvc, T0);
k2 = tvc_ode(x + Ts/2*k1,  u, m, g, J, l_tvc, T0);
k3 = tvc_ode(x + Ts/2*k2,  u, m, g, J, l_tvc, T0);
k4 = tvc_ode(x + Ts*k3,    u, m, g, J, l_tvc, T0);

x_next = x + (Ts/6) * (k1 + 2*k2 + 2*k3 + k4);
end
