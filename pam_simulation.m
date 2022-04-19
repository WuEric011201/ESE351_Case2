function [xn_est] = pam_simulation(p, dt, Tp, Ts, xn)

    fb = 1/Ts;
    tp = -Tp : dt : Tp;
    
    N = length(xn);
    
    tImp = 0 : dt : Ts*(N-1);
    
    impTrain = zeros(size(tImp));
    
    impCnt = 1;
    xnCnt = 1;
    
    for i = 1:length(tImp)
        if (i == impCnt)
            impTrain(i) = xn(xnCnt);
            impCnt = impCnt + Ts/dt;
            xnCnt = xnCnt + 1;
        end
    end
    
    y = conv(impTrain, p);
    y = y(Tp/dt+1 : length(y)-Tp/dt);
        
    figure;
    plot(f,fft(y,1000));
    
    y = y .* cos(2*pi*20*tImp);
    
    r = y + randn(1,length(y));
    
    H = double(f <= 10 | f(length(f)) - f <= 10);
    h = ifft(H) * length(f);
    
    figure;
    plot(h);
    
    y_rec = r .* cos(2*pi*20*tImp) * 2;
        
    z = conv(r,flip(p));
    z = z(Tp/dt+1 : length(z)-Tp/dt);
    
    xn_est = zeros(1, N);
    
    for i = 1:N
        if (z(int64(1+(i-1)*Ts/dt)) > 0)
            xn_est(i) = 1;
        else
            xn_est(i) = -1;
        end
    end
    
end