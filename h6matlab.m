%% Homework 6 Matlab
% * Author:                   Jasmine Cheng, Tong Wu, Zhongli Tong
% * Class:                    ESE 351
% * Date:                     4/7, 2022

%% part i: pulse shape p(t)

Tp = 0.1; % fixed 
dt = Tp/50; % fixed 
tp = -Tp : dt : Tp; % time interval of p
p = (1-abs(tp)./Tp);
P = fftshift(fft(p)); % fourier transform of p

figure;
subplot(3,1,1);
plot(tp,p); title('pulse shape p(t)');
xlabel('time(s)'); ylabel('amplitude');

subplot(3,1,2);
plot(linspace(-pi,pi,length(tp)),abs(P)); title('magnitude of P(jw)');
xlabel('frequency'); ylabel('magnitude');

subplot(3,1,3);
plot(linspace(-pi,pi,length(tp)),angle(P)); title('phase of P(jw)');
xlabel('frequency'); ylabel('phase');

%% part ii, iii, iv : y(t), r(t), xn

% two bit rates
fb = [0.5/Tp 1/Tp];

N = 11; % choose a large N for accuracy
xn = 2*((rand(1,N)>0.5)-0.5); % binary message

for f = 1:length(fb)

    Ts = 1/fb(f);

    % generate the impulse train which convolve with p(t)
    tImp = -Ts*(N/2)+Tp : dt : Ts*(N/2)-Tp;  % time interval of impulse train
    impTrain = zeros(size(tImp));
    tCnt = 1; % a counter saving the index of next impulse
    xnCnt = 1; % a counter saving the next value of xn
    for i = 1:length(tImp)
        if (i == tCnt)
            impTrain(i) = xn(xnCnt);
            tCnt = tCnt + Ts/dt;
            xnCnt = xnCnt + 1;
        end
    end
    
    % generate y(t)
    ty = -Ts*(N/2) : dt : Ts*(N/2); % time interval of y
    y = conv(p,impTrain); % the noise-free transmitted signal y(t)
    r = y + randn(1,length(y)); % the noisy received signal r(t) 

    figure;
    subplot(2,1,1);
    plot(y); title(sprintf('the noise-free transmitted signal y(t) with fb = %d', fb(f)));
    xlabel('time(s)'); ylabel('amplitude');
    
    subplot(2,1,2);
    plot(r); title(sprintf('the noisy received signal r(t) with fb = %d', fb(f)));
    xlabel('time(s)'); ylabel('amplitude');

    % initialize the estimated xn
    xnEst_r = zeros(1,N);
    xnEst_z = zeros(1,N);
    
    % generate z(t)
    z = conv(r,flip(p));
    
    % xnEst using r(t)
    % calculate xn based on certain value of r
    for n = 1:N
    
        % estimate the nth value using r(t)
        if (r(Tp/dt+(n-1)*(1/fb(f)/dt)) > 0)
            xnEst_r(n) = 1;
        else
            xnEst_r(n) = -1;
        end
        
        % estimate the nth value using z(t)
        if (z(2*Tp/dt+(n-1)*(1/fb(f)/dt)) > 0)
            xnEst_z(n) = 1;
        else
            xnEst_z(n) = -1;
        end
    
    end
    
    figure;
    subplot(3,1,1);
    stem(xn); title('xn'); xlabel('n');

    subplot(3,1,2);
    stem(xnEst_r); title(sprintf('estimated xn using a sign-based receiver with fb = %d', fb(f))); xlabel('n');

    subplot(3,1,3);
    stem(xnEst_z); title(sprintf('estimated xn using a matched filter with fb = %d', fb(f))); xlabel('n');

end

%% performance analysis

% two bit rates
fb = [0.5/Tp 1/Tp];

N = 10001; % choose a large N for accuracy
xn = 2*((rand(1,N)>0.5)-0.5); % binary message

% range of sigma
sigma = 0:.01:1;

% initialize the error rate
error_r = zeros(length(fb),length(sigma)); % using r(t)
error_z = zeros(length(fb),length(sigma)); % using z(t)

for s = 1:length(sigma)

    for f = 1:length(fb)
    
        Ts = 1/fb(f);
    
        % generate the impulse train which convolve with p(t)
        tImp = -Ts*(N/2)+Tp : dt : Ts*(N/2)-Tp;  % time interval of impulse train
        impTrain = zeros(size(tImp));
        tCnt = 1; % a counter saving the index of next impulse
        xnCnt = 1; % a counter saving the next value of xn
        for i = 1:length(tImp)
            if (i == tCnt)
                impTrain(i) = xn(xnCnt);
                tCnt = tCnt + Ts/dt;
                xnCnt = xnCnt + 1;
            end
        end
        
        % generate y(t)
        ty = -Ts*(N/2) : dt : Ts*(N/2); % time interval of y
        y = conv(p,impTrain); % the noise-free transmitted signal y(t)
        r = y + sigma(s)*randn(1,length(y)); % the noisy received signal r(t) 

        % initialize the estimated xn
        xnEst_r = zeros(1,N);
        xnEst_z = zeros(1,N);
        
        % generate z(t)
        z = conv(r,flip(p));
        
        % xnEst using r(t)
        % calculate xn based on certain value of r
        for n = 1:N
        
            % estimate the nth value using r(t)
            if (r(Tp/dt+(n-1)*(1/fb(f)/dt)) > 0)
                xnEst_r(n) = 1;
            else
                xnEst_r(n) = -1;
            end
            
            % estimate the nth value using z(t)
            if (z(2*Tp/dt+(n-1)*(1/fb(f)/dt)) > 0)
                xnEst_z(n) = 1;
            else
                xnEst_z(n) = -1;
            end
        
        end
        
        % calculate the error rate
        error_r(f,s) = sum(xnEst_r ~= xn)/N;
        error_z(f,s) = sum(xnEst_z ~= xn)/N;
    
    end

end

% plot the error rate of different bit rate and receiver method
for i = 1:length(fb)

    figure;
    plot(sigma,error_r(i,:)); xlabel('sigma'); ylabel('error rate');
    title(sprintf('error rate of sign-based receiver at fb = %d', fb(i)));
    
    figure;
    plot(sigma,error_z(i,:)); xlabel('sigma'); ylabel('error rate');
    title(sprintf('error rate of matched filter at fb = %d', fb(i)));

end

% In this matlab homework, we use triangle pulse. In this part, we compare
% two different bit rates - 1/Tp and 1/2Tp - and two receiver methods -
% sign-based receiver and matched filter - over a range of sigma from 0 to
% 1. We use N=10001 to ensure the accuracy of the result.
% From the graphs, we observe matched filter performs far better than
% sign-based receiver under the condition of a large N: there is no error
% when fb is 1/2Tp and a little error at high sigma when fb is 1/Tp using
% matched receiver while there is an obvious increase in error rate when
% sigma increases using sign-based receiver under both bit rates. Simply
% processing the noisy signal by convolving with a fliped pulse vastly
% increase the accuracy of decoding the signal. 
% Bit rates also influence the error rate of decoding. 1/2Tp, where Tp is
% half of the bandwidth of triangular pulse in time domain, has a better
% performance over bit rate 1/Tp. This is because convolving the triangular
% pulse with impulse train with frequency 1/2Tp remains the intact shape of
% the pulse. On the other hand, using a impulse train of 1/Tp, parts of the
% pulse overlapped which makes the signal harder to decode. Using the
% sign-based receiver, error rate is always above 0.2 when sigma is 1 when
% fb = 1/Tp while it's usually about 0.16 when fb = 1/2Tp. Using the
% matched filter, error rates remain 0 when fb = 1/2Tp while there is some
% error when fb = 1/Tp. 
% An increase in sigma, sqaure root of the variance of the Gaussian
% additive noise, leads to an increase in error rate as well. There is an
% obviouse positive correlation between sigma and error rate when using the
% sign-based receiver. Though this relationship is not apparent when using
% the matched filter, we can still observe some error at high sigma at the
% bit rate of 1/Tp while there usually is no error when sigma is low, below
% 0.7.
