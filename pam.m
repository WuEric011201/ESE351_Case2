% function pam_simulation
% parameters:
% p: pulse shape
% N: length of the signal
% dt: time spacing between pulse shape sampling
% Tp: time range of pulse shape (time vector = -Tp : dt : Tp)
% Ts: symbol time
% fc: 1 by 3, carrier frequency in hertz
% sigma: related to noise
% xn: three row input signal 
% 
% return:
% xn_est: 1 by 3, decoded xn

function [tImp, r, y, y_total, y_up, y_rec, xn_est] = pam(p, xn,  dt, Tp, Ts, fc, sigma, count)

    % assuming length of xn1, xn2, xn3 are the same

    % initialization
    fs = 1/dt; % sampling frequency
    N = size(xn, 2);  % Length of xn
    tImp = 0 : dt : Ts*(N-1); % time vector of the impulse train
    y_up = zeros(3, length(tImp)); 
    y = zeros(3, length(tImp)); 
    xn_est= zeros(3, N); 
    y_rec = zeros(3, length(tImp)); 
    
    % generate the impulse train
    impTrain = zeros(3, size(tImp, 2)); % initialization

%     ----------------- UP CONVERT -------------------
    for num= 1:count
        impCnt = 1; % pointer: the index of the next pulse
        xnCnt = 1; % pointer: the next message
    
        for i = 1 : size(tImp, 2)
            if (i == impCnt)
                impTrain(num, i) = xn(num, xnCnt);
                impCnt = impCnt + Ts/dt;
                xnCnt = xnCnt + 1;
            end
        end

        % convolve pulse shape with the impulse train
        result_conv = conv(impTrain(num, : ), p);
        y(num, : ) = result_conv(Tp/dt+1 : length(result_conv)-Tp/dt); % Make the y same length as the impulse train

        % up convert
        y_up(num, :) = y(num, : ) .* cos(2*pi*fc(1, num)*tImp);
    end

    y_total = sum(y_up, 1); 


% -------- Transmit the signal -------------------------


    % add noise
    r = y_total + sigma .*randn(1,length(tImp));


% ------------- Down convert ----------------
    for num = 1: count

        % down convert 
        y_down = r .* cos(2*pi*fc(1, num)*tImp) * 2;
        y_rec(num, : ) = lowpass(y_down, 5, fs); % Passing through a low pass filter of frequency 5Hz
    
        % matched filter receiver
        z = conv(y_rec(num, : ),flip(p));
        z = z(Tp/dt+1 : length(z)-Tp/dt);
        
        % calculate xn_est
        
        for i = 1:N
            if (z(int64(1+(i-1)*Ts/dt)) > 0)
                xn_est(num, i) = 1;
            else
                xn_est(num, i) = -1;
            end
        end

    end

end