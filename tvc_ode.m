function xdot = tvc_ode(x, u, m, g, J, l_tvc, T0)
    % TVC_ODE  2-D nonlinear equations of motion for a TVC rocket

    vy = x(3);
    vz = x(4);
    theta = x(5);
    q = x(6);

    delta = u(1);
    dT = u(2);

    T = T0 + dT;
    phi = theta + delta;

    xdot = [vy;
            vz;
            T * sin(phi) / m;
            T * cos(phi) / m - g;
            q;
            T * l_tvc * sin(delta) / J];
end
