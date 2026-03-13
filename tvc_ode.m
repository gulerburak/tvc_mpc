function xdot = tvc_ode(x, u, m, g, J, l_tvc, T0)
% TVC_ODE  2-D nonlinear equations of motion for a TVC rocket
%
% Physical setup:
%   The rocket flies in the vertical (y-z) plane.
%   The nozzle can be gimballed by angle delta relative to the body axis
%   to redirect thrust, creating both a lateral force and a pitch torque.
%
% State  x = [y; z; vy; vz; theta; q]
%   y, z     = lateral (horizontal) and axial (vertical) position [m]
%   vy, vz   = corresponding velocities [m/s]
%   theta    = pitch angle from vertical, positive rightward (+y) [rad]
%   q        = pitch rate  dtheta/dt  [rad/s]
%
% Input  u = [delta; dT]
%   delta    = gimbal deflection angle relative to body axis [rad]
%              positive delta -> thrust deflected rightward -> torque increases theta
%   dT       = thrust deviation above the nominal T0 [N]
%
% Torque about CG from TVC:
%   tau = (T0+dT) * l_tvc * sin(delta)
%   (derived from moment arm: nozzle at -l_tvc along body axis,
%    force component perpendicular to moment arm = T*sin(delta))

vy    = x(3);
vz    = x(4);
theta = x(5);
q     = x(6);

delta = u(1);
dT    = u(2);

T = T0 + dT;                         % total thrust [N]
phi = theta + delta;                  % total thrust angle from vertical

xdot = [ vy;
         vz;
         T * sin(phi) / m;           % y acceleration
         T * cos(phi) / m - g;       % z acceleration (gravity opposes)
         q;
         T * l_tvc * sin(delta) / J];% pitch angular acceleration
end
