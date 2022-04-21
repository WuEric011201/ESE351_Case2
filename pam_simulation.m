% function pam_simulation
% parameters:
% p: pulse shape
% dt: time spacing between pulse shape sampling
% Tp: time range of pulse shape (time vector = -Tp : dt : Tp)
% Ts: symbol time
% xn: sent message
% fc: carrier frequency in hertz
% sigma: related to noise
% 
% return:
% xn_est: decoded xn

function [xn_est] = pam_simulation(p, dt, Tp, Ts, xn, fc, sigma)
    
    N = length(xn); % length of the sent message
    fs = 1/dt; % sampling frequency

    tImp = 0 : dt : Ts*(N-1); % time vector of the impulse train
    
    % generate the impulse train
    impTrain = zeros(size(tImp)); % initialization
    impCnt = 1; % pointer: the index of the next pulse
    xnCnt = 1; % pointer: the next message
    
    for i = 1:length(tImp)
        if (i == impCnt)
            impTrain(i) = xn(xnCnt);
            impCnt = impCnt + Ts/dt;
            xnCnt = xnCnt + 1;
        end
    end

    % convolve pulse shape with the impulse train
    y = conv(impTrain, p);
    y = y(Tp/dt+1 : length(y)-Tp/dt);
    
    figure;
    plot(tImp,y);
    title('Noise-free PAM signal'); xlabel('time'); ylabel('y(t)');
    
    % up convert
    y_up = y .* cos(2*pi*fc*tImp);

    figure; 
    plot(tImp,y_up);
    title('Up-converted signal'); xlabel('time'); ylabel('y(t)cos(wc*t)');
    
    % add noise
    r = y_up + sigma*randn(1,length(y));
    
    figure;
    plot(tImp,r);
    title('Noisy received signal'); xlabel('time'); ylabel('r(t)');
    
    % down convert
    y_down = r .* cos(2*pi*fc*tImp) * 2;
    y_rec = lowpass(y_down, 5, fs);
    
    figure();
    plot(tImp,y_rec);
    title('Down-converted signal'); xlabel('time'); ylabel('yrec(t)');
    
    % matched filter receiver
    z = conv(y_rec,flip(p));
    z = z(Tp/dt+1 : length(z)-Tp/dt);
    
    % calculate xn_est
    xn_est = zeros(1, N); % initialization
    
    for i = 1:N
        if (z(int64(1+(i-1)*Ts/dt)) > 0)
            xn_est(i) = 1;
        else
            xn_est(i) = -1;
        end
    end

    figure;
    subplot(2,1,1); stem(xn); title('sent message');
    subplot(2,1,2); stem(xn_est); title('decoded message');
    
end