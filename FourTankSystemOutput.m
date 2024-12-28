function z = FourTankSystemOutput(x,p)
    A = p(5:6,1);
    rho = p(12);
    z = x./(rho*A); % Liquid level in each tank [cm]
end