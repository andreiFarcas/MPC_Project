function u0 = Regulator(A,B,Bv,C,D,N,xinit,uprev,zref,dist)
    n = length(A); % nbr of state variables
    m = length(B(1,:)); % nbr of control variables
    o = length(C(:,1)); % nbr of output variables
    q = length(Bv(1,:)); % nbr of perturbation variables

    % Compute psi, gamma, gamma_d
    phi = cell(N, 1); 
    for i = 1:N
        phi{i} = C * (A^i); % Calculer C * A^i et stocker dans phi
    end
    phi = cell2mat(phi); % shape (pN x n)
    
    markov = MarkovParam(A,B,C,D,N,[]);
    markov_d = MarkovParam(A,Bv,C,D,N,[]);
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
    
    % Build the matrices Qz, Hs, Muprev
    Qz = blkdiag(eye(m), zeros(m*(N-1)));
    
    central_diag = 2 * ones(1, m*N); 
    central_diag(N) = 1;
    inf_diag = -1 * ones(1, m*N-1);       
    sup_diag = -1 * ones(1, m*N-1);        
    Hs = diag(central_diag) + diag(inf_diag, -1) + diag(sup_diag, 1);
    Hs(2*N,2*N)=1; % shape (mN x mN) 
    
    Muprev = zeros(m*N,m); % shape (mN x m)
    Muprev(1,:) = 1;
    
    % Compute Mxo, Mr, Md
    Mx0 = gamma*Qz*phi; % shape (mN x n)
    Mr = -gamma'*Qz; % shape (mN x pN)
    Md = gamma'*Qz*gamma_d;  % shape (mN x qN)
    
    % Compute H and g
    H = gamma'*Qz*gamma + Hs;  % shape (mN x mN)
    grad = Mx0*xinit + Mr*zref + Muprev*uprev + Md*dist;  % shape (mN x 1)
    
    % Solve the QP
    [U,~]= qpsolver(H, grad, A, xinit,[],[],[],[],0);
    
    % Extract the first control variable
    u0 = U(1:2,1);
end