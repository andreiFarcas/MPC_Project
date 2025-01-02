function u0 = Regulator_InputConstrained(A, B, Bv, C, D, N, xinit, uprev, zref, dist, Qz, Hs, u_min, u_max)
    n = length(A); % Number of state variables
    m = length(B(1,:)); % Number of control variables
    o = length(C(:,1)); % Number of output variables
    q = length(Bv(1,:)); % Number of perturbation variables

    % Compute psi, gamma, gamma_d
    phi = cell(N, 1); 
    for i = 1:N
        phi{i} = C * (A^i);
    end
    phi = cell2mat(phi); % shape (pN x n)
    
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
    gamma = cell2mat(gamma); % Shape (mN x pN)
    gamma_d = cell2mat(gamma_d); % Shape (mN x qN)
    
    % Build the matrice Muprev
    Muprev = zeros(m*N, m); % shape (mN x m)
    Muprev(1:m, :) = eye(m);
    
    % Compute Mxo, Mr, Md
    Mx0 = gamma * Qz * phi; % shape (mN x n)
    Mr = -gamma' * Qz; % shape (mN x pN)
    Md = gamma' * Qz * gamma_d; % shape (mN x qN)
    
    % Compute H and g
    H = gamma' * Qz * gamma + Hs;  % shape (mN x mN)
    grad = Mx0 * xinit + Mr * zref + Muprev * uprev + Md * dist;  % shape (mN x 1)
    
    % Constraints on U
    % Lower and upper bounds for all time steps
    Umin = repmat(u_min, N, 1); 
    Umax = repmat(u_max, N, 1); 
    
    % Define G_u and h_u for box constraints
    %G_u = [eye(m * N); -eye(m * N)]; % shape (2mN x mN)
    G_u = eye(m * N);
    h_u = [Umax; -Umin]; % shape (2mN x 1)
   
    % Solve the QP
    [U, info] = qpsolver(H, grad, G_u, [], Umin, Umax, Umin, Umax, 1);
    
    %options = optimoptions('quadprog', 'Display', 'none'); % Optional settings
    %U = quadprog(H, grad, G_u, h_u, [], [], [], [], [], options);

    % Extract the first control variable
    u0 = U(1:m, 1);
end