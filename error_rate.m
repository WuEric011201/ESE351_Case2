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

function [] = error_rate(p, dt, Tp, Ts, xn, fc, count)
    
    % range of sigma
    s = 0:.01:1;
    error = zeros(size(s));

    for num = 1: count
        for i = 1:length(s)
            [~, ~, ~, ~, ~, ~, xn_est]= pam(p,xn, dt, Tp, Ts, fc, s(i), count);
            error(i) = sum(xn_est(num, : ) ~= xn(num, : ))/size(xn, 2);
        end
    
        figure;
        plot(s,error);
%         title(['sigma vs. error rate of ', num2str(num), 'th signal']); xlabel('sigma'); ylabel('error rate');
        title('sigma vs. error rate');
    end

end