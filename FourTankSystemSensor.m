function y = FourTankSystemSensor(x,p,v)
    A = p(5:8,1);
    rho = p(12);
    y = x./(rho*A) + v; % Liquid level in each tank [cm]
    % qout = a.*sqrt(2*g*h);  % Outflow from each tank [cm3/s]
end