function u0 = OldRegulator_IO_Constraints(A, B, Bv, C, D, N, xinit, uprev, zref, dist, Qz, Hs, u_min, u_max,delta_u_min,delta_u_max,z_min,z_max)
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
    Mx0 = gamma' * Qz * phi; % shape (mN x n)
    Mr = -gamma' * Qz; % shape (mN x pN)
    Md = gamma' * Qz * gamma_d; % shape (mN x qN)
    
    % Compute H and g
    H = gamma' * Qz * gamma + Hs;  % shape (mN x mN)
    grad = Mx0 * xinit + Mr * zref + Muprev * uprev + Md * dist;  % shape (mN x 1)
    
    % Constraints on U
    % Repeat the constraints over the prediction horizon (N steps):
    Umin = repmat(u_min, N, 1);  % Shape: (m * N) x 1
    Umax = repmat(u_max, N, 1);  
       
    % Constraints on DeltaU
    % G_u = [eye(m * N); -eye(m * N)]; % shape (2mN x mN)
    % G_u = eye(m * N);
    Delta_u_min = repmat(delta_u_min, N, 1);
    Delta_u_max = repmat(delta_u_max, N, 1);

    Lambda = zeros(m * N, m * N); % Shape: (mN x mN)
    for i = 1:N
        if i == 1
            Lambda((i-1)*m+1:i*m, (i-1)*m+1:i*m) = eye(m); % First step
        else
            Lambda((i-1)*m+1:i*m, (i-2)*m+1:(i-1)*m) = -eye(m); % -u[k-1]
            Lambda((i-1)*m+1:i*m, (i-1)*m+1:i*m) = eye(m); % u[k]
        end
    end    

    % Constraints on z
    Zmin = repmat(z_min, N, 1);
    Zmax = repmat(z_max, N, 1);
    z_bar_min = Zmin - phi*xinit - gamma_d*dist;
    z_bar_max = Zmax - phi*xinit - gamma_d*dist;


    % Combine all constraints
    G_u = [Lambda; gamma]; % Combine Δ and z constraints
    bl = [Delta_u_min; z_bar_min];
    bu = [Delta_u_max; z_bar_max];

    % Solve the QP
    [U, info] = qpsolver(H, grad, G_u, [], Umin, Umax, bl, bu, 1);

    if isempty(U) || numel(U) < m
    error('QP solver did not return a valid solution.');
    end

   
    
    %options = optimoptions('quadprog', 'Display', 'none'); % Optional settings
    %U = quadprog(H, grad, G_u, h_u, [], [], [], [], [], options);

    % Extract the first control variable
    u0 = U(1:m, 1);
end