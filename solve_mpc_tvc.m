function [u_opt, U_seq, feasible] = solve_mpc_tvc(e0, Ad, Bd, Q, R, P, N, ...
                                                    e_lb, e_ub, u_lb, u_ub, ...
                                                    Xf_A, Xf_b)
% SOLVE_MPC_TVC  YALMIP QP for the TVC MPC (error-state coordinates)
%
% Solves:
%   min   e_N' P e_N  +  sum_{k=0}^{N-1} [ e_k'Q e_k + u_k'R u_k ]
%   s.t.  e_{k+1} = Ad*e_k + Bd*u_k
%         e_lb <= e_k <= e_ub          (hard state constraints)
%         u_lb <= u_k <= u_ub          (hard input constraints)
%         Xf_A * e_N  <= Xf_b          (terminal set, if provided)
%         e_0 = e0
%
% Arguments
%   e0          initial error state (nx x 1)
%   Ad, Bd      discrete-time LTI matrices
%   Q, R, P     stage / terminal weight matrices
%   N           prediction horizon (steps)
%   e_lb, e_ub  state constraint bounds (nx x 1)
%   u_lb, u_ub  input constraint bounds (nu x 1)
%   Xf_A, Xf_b  H-rep of terminal set: Xf_A*x <= Xf_b  (pass [] to omit)
%
% Returns
%   u_opt    first optimal input  (nu x 1)
%   U_seq    full optimal input sequence  (nu x N)
%   feasible true if QP solved to optimality

nx = size(Ad, 1);
nu = size(Bd, 2);

E = sdpvar(nx, N+1, 'full');
U = sdpvar(nu, N,   'full');

cost        = E(:,N+1)' * P * E(:,N+1);   % terminal cost
constraints = [ E(:,1) == e0 ];

for k = 1:N
    cost = cost + E(:,k)'*Q*E(:,k) + U(:,k)'*R*U(:,k);
    constraints = [ constraints, ...
        E(:,k+1) == Ad*E(:,k) + Bd*U(:,k), ...
        e_lb <= E(:,k+1) <= e_ub, ...   % constrain predicted states only (not E(:,1)=e0)
        u_lb <= U(:,k) <= u_ub ];
end

% Optional terminal set
if ~isempty(Xf_A)
    constraints = [ constraints, Xf_A * E(:,N+1) <= Xf_b ];
end

opts        = sdpsettings('solver','quadprog','verbose',0, ...
                          'quadprog.Display','off');
diagnostics = optimize(constraints, cost, opts);

if diagnostics.problem == 0
    u_opt    = value(U(:,1));
    U_seq    = value(U);
    feasible = true;
else
    u_opt    = zeros(nu,1);
    U_seq    = zeros(nu,N);
    feasible = false;
end
end
