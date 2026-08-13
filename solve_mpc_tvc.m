function [u_opt, U_seq, feasible] = solve_mpc_tvc(e0, Ad, Bd, Q, R, P_inf, beta, N, ...
        e_lb, e_ub, u_lb, u_ub, solver)
    % SOLVE_MPC_TVC
    % Arguments
    %   e0          initial error state (nx x 1)
    %   Ad, Bd      discrete-time LTI matrices
    %   Q, R        stage weight matrices
    %   P_inf       terminal cost matrix (DARE solution)
    %   beta        terminal cost weight (>= 1)
    %   N           prediction horizon (steps)
    %   e_lb, e_ub  state constraint bounds (nx x 1)
    %   u_lb, u_ub  input constraint bounds (nu x 1)
    %   solver      optional YALMIP solver name; '' lets YALMIP choose
    %
    % Returns
    %   u_opt    first optimal input  (nu x 1)
    %   U_seq    full optimal input sequence  (nu x N)
    %   feasible true if QP solved to optimality

    if nargin < 13 || isempty(solver)
        % Prefer quadprog when the Optimization Toolbox is installed, otherwise
        % let YALMIP pick whatever QP solver is on the path.
        if isempty(which('quadprog'))
            solver = '';
        else
            solver = 'quadprog';
        end
    end

    nx = size(Ad, 1);
    nu = size(Bd, 2);

    E = sdpvar(nx, N + 1, 'full');
    U = sdpvar(nu, N, 'full');

    % The Code Analyzer warnings suppressed below are false positives on YALMIP's
    % overloaded operators: growing a constraint list by concatenation and the
    % double-sided inequality lb <= x <= ub are both standard YALMIP syntax.
    cost = beta * (E(:, N + 1)' * P_inf * E(:, N + 1));
    constraints = [E(:, 1) == e0]; %#ok<NBRAK2>

    for k = 1:N
        cost = cost + E(:, k)' * Q * E(:, k) + U(:, k)' * R * U(:, k);
        constraints = [constraints, ...
                           E(:, k + 1) == Ad * E(:, k) + Bd * U(:, k), ...
                           e_lb <= E(:, k + 1) <= e_ub, ... %#ok<CHAIN>
                           u_lb <= U(:, k) <= u_ub]; %#ok<AGROW,CHAIN>
    end

    opts = sdpsettings('solver', solver, 'verbose', 0);

    if strcmp(solver, 'quadprog')
        opts = sdpsettings(opts, 'quadprog.Display', 'off');
    end

    diagnostics = optimize(constraints, cost, opts);

    if diagnostics.problem == 0
        u_opt = value(U(:, 1));
        U_seq = value(U);
        feasible = true;
    else
        u_opt = zeros(nu, 1);
        U_seq = zeros(nu, N);
        feasible = false;
    end

end
