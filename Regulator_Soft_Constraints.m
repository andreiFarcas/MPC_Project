function u0 = Regulator_Soft_Constraints(A, B, Bv, C, D, N, xinit, uprev, zref, dist, Qz, Hs, u_min, u_max, delta_u_min, delta_u_max, z_min, z_max, rho2, rho1)
    % Initialization of dimensions
    n = length(A); % Number of state variables
    m = length(B(1,:)); % Number of control variables
    o = length(C(:,1)); % Number of output variables
    q = length(Bv(1,:)); % Number of disturbance variables

    % Matrices phi, gamma, gamma_d
    phi = cell(N, 1); 
    for i = 1:N
        phi{i} = C * (A^i);
    end
    phi = cell2mat(phi); % Shape (oN x n)

    markov = MarkovParam(A, B, C, D, N, []); 
    markov_d = MarkovParam(A, Bv, C, D, N, []);
    gamma = cell(N, N); 
    gamma_d = cell(N, N);
    for i = 1:N
        for j = 1:N
            if (i - j + 1) > 0
                gamma{i, j} = markov{i - j + 1}; 
                gamma_d{i, j} = markov_d{i - j + 1}; 
            else
                gamma{i, j} = zeros(size(markov{1}));
                gamma_d{i, j} = zeros(size(markov_d{1}));
            end
        end
    end
    gamma = cell2mat(gamma); % Shape (oN x mN)
    gamma_d = cell2mat(gamma_d); % Shape (oN x qN)

    % Matrices Mx0, Mr, Md
    Mx0 = gamma' * Qz * phi; % Shape (mN x n)
    Mr = -gamma' * Qz; % Shape (mN x oN)
    Md = gamma' * Qz * gamma_d; % Shape (mN x qN)

    % Update the cost function to account for soft constraints
    H = gamma' * Qz * gamma + Hs;  % Original cost (mN x mN)
    grad = Mx0 * xinit + Mr * zref + Md * dist;  % Gradient (mN x 1)

    % Add penalty terms for the soft constraints
    H_soft = blkdiag(H, rho2 * eye(o * N)); % Add weights for eta_k
    grad_soft = [grad; rho1 * ones(o * N, 1)]; % Updated gradient with L1 term

    % Constraints on U
    Umin = repmat(u_min, N, 1); % Shape (mN x 1)
    Umax = repmat(u_max, N, 1);

    % Constraints on DeltaU
    Delta_u_min = repmat(delta_u_min, N, 1);
    Delta_u_max = repmat(delta_u_max, N, 1);
    Delta_u_min(1:2) = Delta_u_min(1:2) + uprev;
    Delta_u_max(1:2) = Delta_u_max(1:2) + uprev;

    % Matrix Lambda for DeltaU
    Lambda = zeros(m * N, m * N);
    for i = 1:N
        if i == 1
            Lambda((i-1)*m+1:i*m, (i-1)*m+1:i*m) = eye(m); % First step
        else
            Lambda((i-1)*m+1:i*m, (i-2)*m+1:(i-1)*m) = -eye(m); % -u[k-1]
            Lambda((i-1)*m+1:i*m, (i-1)*m+1:i*m) = eye(m); % u[k]
        end
    end    

    % Constraints on z (with eta_k)
    Etamin =  zeros(o * N, 1); % eta_k > 0
    Etamax = inf(o * N, 1);

    Zmin = repmat(z_min, N, 1); % Lower bounds on outputs
    Zmax = repmat(z_max, N, 1); % Upper bounds on outputs
    z_bar_min = Zmin - phi*xinit - gamma_d*dist; % Account for initial state and disturbance
    z_bar_max = Zmax - phi*xinit - gamma_d*dist;    

    G_u = [Lambda, zeros(m * N, o * N)];
    bl_u = Delta_u_min;
    bu_u = Delta_u_max;

    G_eta = [gamma, -eye(o * N); gamma, eye(o * N)];
    bl_eta = [1e-6 * ones(o * N, 1); z_bar_min];
    bu_eta = [z_bar_max; inf(o * N, 1)];

    % Combine all 
    Umin_tot = [Umin;Etamin];
    Umax_tot = [Umax;Etamax];
    G_total = [G_u; G_eta]; % Add eta_k constraints
    bl_total = [bl_u; bl_eta];
    bu_total = [bu_u; bu_eta];

    % Solve the QP with soft constraints
    [U_eta, info] = qpsolver(H_soft, grad_soft, G_total, [], Umin_tot, Umax_tot, bl_total, bu_total, 1);

    if isempty(U_eta) || numel(U_eta) < m
        error('QP solver did not return a valid solution.');
    end

    % Extract the first control input
    u0 = U_eta(1:m, 1);
end
