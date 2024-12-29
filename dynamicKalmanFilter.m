function [x_hat, P] = dynamicKalmanFilter(Ff, Gg, Gv, Cc, Q, R, y, u, x_hat_prev, P_prev)
    % Prediction Step
    x_pred = Ff * x_hat_prev+Gg * u; % Predicted state
    P_pred = Ff * P_prev * Ff'+ Gv * Q * Gv'; % Predicted covariance

    % Update step
    K = P_pred * Cc' / (Cc * P_pred * Cc'+R); % Kalman Gain
    x_hat = x_pred + K* (y - Cc * x_pred); % Updated state estimate
    P = (eye(size(P_prev)) - K* Cc) * P_pred; % Updated covariance
end