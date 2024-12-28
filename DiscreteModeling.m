function [taus, gains,Ts,Ff,Gg,Gv,Cc,Dd] = DiscreteModeling(A,B,Bv,C,D,F1,F2,modelname)
    sys = ss(A,[B, Bv],C,D); 
    G = tf(sys);

    if ~isempty(modelname)
        figure;
        step(0.1*F1*G(1,1), 0.1*F2*G(1,2), 0.1*F1*G(2,1), 0.1*F2*G(2,2));
        legend('G11', 'G12', 'G21', 'G22');
        titleStr = [modelname,': Step response of 0.1*F1 and 0.1*F2']; 
        title(titleStr); 
    end

    gains = dcgain(G); % Steady-state gain (DC gain)
    poles = pole(G);
    taus = -1 ./ real(poles); % Time constants

    % Discrete-time state space models
    tau_min = min(taus);
    Ts = tau_min/10; %[s]
    [Ff,Gg] = c2d(A,B,Ts);
    [~,Gv] = c2d(A,Bv,Ts);
    Cc = C;
    Dd=D; 

end