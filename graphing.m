% This function graphs all the signals being processed in the
% pam_simulation file. 

function []  = graphing(tImp, fs, y, y_up, y_total, r, y_rec, xn, xn_est, num)

     f = -fs/2 : fs/1000 : fs/2-fs/1000; % frequency range % fs/L is one step in frequency

    for i = 1: num
        % Graph the signals in the process in the time domain
        figure;
        subplot(2, 2, 1); 
        plot(tImp,y(num, : ));
        title(['Noise-free PAM of signal x', num2str(num)] ); xlabel('time'); ylabel('y(t)');
        subplot(2, 2, 2); 
        plot(tImp,y_up(num, : ));
        title(['Up-converted signal x', num2str(num)] ); xlabel('time'); ylabel('y(t)cos(wc*t)');
        subplot(2, 2, 3); 
        plot(tImp,r);
        title(['Noisy received signal x',  num2str(num)] ); xlabel('time'); ylabel('r(t)');
        subplot(2, 2, 4); 
        plot(tImp,y_rec(num, : ));
        title(['Down-converted signal x', num2str(num)] ); xlabel('time'); ylabel('yrec(t)');
    
        % Graph the signals in the process in the frequency domain
        figure;
        subplot(2, 2, 1); 
        plot(f,fftshift(fft(y(num, : ),1000))); 
        title(['Fourier Transform of noise-free PAM of signal x', num2str(num)] ); xlabel('Frequency'); ylabel('Y(jw)');
        subplot(2, 2, 2); 
        plot(f,fftshift(fft(y_up(num, : ),1000))); 
        title(['Fourier Transform of up-converted of signal x' , num2str(num)] ); xlabel('Frequency'); ylabel('Y_up(jw)');
        subplot(2, 2, 3); 
        plot(f,fftshift(fft(r,1000))); 
        title(['Fourier Transform of noisy received of signal x' , num2str(num)] ); xlabel('Frequency'); ylabel('R(jw)');
        subplot(2, 2, 4); 
        plot(f,fftshift(fft(y_rec(num, : ),1000))); 
        title(['Fourier Transform of down-converted of signal x' , num2str(num)] ); xlabel('Frequency'); ylabel('Yrec(jw)');

        % Graph the final comparison
        figure; 
        subplot(2,1,1); plot(xn(num, : )); title('sent message');
        subplot(2,1,2); plot(xn_est(num, : )); title('decoded message');
    end
    
    % Graph the total signal in the transmission
    figure; 
    subplot(2, 1, 1); 
    plot(tImp,y_total);
    title('Total signal in the transmission'); xlabel('time'); ylabel('y_total(t)');
    subplot(2, 1, 2); 
    plot(f,fftshift(fft(y_total,1000))); 
    title('Fourier Transform of total signal in the transmission'); xlabel('Frequency'); ylabel('Y_total(jw)');

end