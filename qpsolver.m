function [x, info] = qpsolver(hessian, gradient, A, xinit, lx, ux, bl, bu,constraints)
    %  Solves a convex quadratic programming problem using quadprog
    n = length(gradient);
    
    % quadprog expects inequality constraints in the form A*x <= b
    % we have bl <= A*x <= bu
    % so we make two inequatlities
    %   A*x <= bu
    %   -A*x <= -bl
    
    options = optimoptions('quadprog', 'Display', 'off');

    if constraints == 1 % constraints
        % A*x <= bu
        % -A*x <= -bl
        A_ineq = [A; -A];
        b_ineq = [bu; -bl];
    else % no constraints
        A_ineq = [];
        b_ineq = [];
        A = [];
        bl = [];
        bu = [];
    end


    [x, fval, ~, ~] = quadprog(hessian, gradient, A_ineq, b_ineq, [], [], lx, ux, xinit, options);

    
    info.fval = fval; % objective function value

end
