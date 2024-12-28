function [xs,ys,zs,A,B,Bv,C,D] = StateSpaceModeling(p,xs0,us,d,v)
    % xs0: initial guess on x0
    % us: flow rates
    % d: model disturbances
    % v: measurement noise

    steadyStateFunction = @(x) ModifiedFourTankSystem(0, x, us, d, p); % Define the function handle for steady-state
    xs = fsolve(steadyStateFunction, xs0); % Solve using fsolve
    ys = FourTankSystemSensor(xs,p,v); % [cm]
    zs = FourTankSystemOutput(xs(1:2),p); % [cm]
    h1s = ys(1); h2s = ys(2); h3s = ys(3); h4s = ys(4);
    
    % Linearization
    a1 = p(1); a2 = p(2); a3 = p(3); a4 = p(4);
    A1 = p(5); A2 = p(6); A3 = p(7); A4 = p(8);
    gamma1 = p(9); gamma2 = p(10);
    g=p(11); rho=p(12);
    
    A = zeros(4,4);
    A(1,1) = -(a1/A1) * sqrt(g*A1*rho / (2 * xs(1)));
    A(1,3) = (a3/A3) * sqrt(g*A3*rho / (2 * xs(3)));
    A(2,2) = -(a2/A2) * sqrt(g*A2*rho / (2 * xs(2)));
    A(2,4) = (a4/A4) * sqrt(g*A4*rho / (2 * xs(4)));
    A(3,3) = -(a3/A3) * sqrt(g*A3*rho / (2 * xs(3)));
    A(4,4) = -(a4/A4) * sqrt(g*A4*rho / (2 * xs(4)));
    
    B=[rho*gamma1 0;0 rho*gamma2; 0 rho*(1-gamma2); rho*(1-gamma1) 0];
    
    Bv = [0 0; 0 0; rho 0; 0 rho];
    C = [1/(rho*A1) 0 0 0; 0 1/(rho*A2) 0 0];
    
    D = 0;
end
