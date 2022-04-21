% function error_rate: return the error rate vs. sigma plot of a pulse
% shape
% parameters:
% p: pulse shape
% dt: time spacing between pulse shape sampling
% Tp: time range of pulse shape (time vector = -Tp : dt : Tp)
% Ts: symbol time
% xn: sent message
% fc: carrier frequency in hertz
% 
% WARNING: COMMENT OUT ALL THE PLOT() FUNCTION IN PAM_SIMULATION() OR YOUR
% COMPUTER WILL EXPLOTE

function [] = error_rate(p, dt, Tp, Ts, xn, fc)
    
    % range of sigma
    s = 0:.01:1;
    error = zeros(size(s));

    for i = 1:length(s)
        xn_est = pam_simulation(p, dt, Tp, Ts, xn, fc, s(i));
        error(i) = sum(xn_est ~= xn)/length(xn);
    end

    figure;
    plot(s,error);
    title('sigma vs. error rate'); xlabel('sigma'); ylabel('error rate');

end