#! /usr/bin/octave

# read audio files

[y1, fs1] = audioread('audio1.wav');
[y2, fs2] = audioread('audio2.wav');
[y3, fs3] = audioread('audio3.wav');

[y1, y2, y3] = deal(y1(:,1), y2(:,1), y3(:,1));
disp(length(y1));
disp(length(y2));
disp(length(y3));

max_length = max([length(y1), length(y2), length(y3)]);
append_zeros = @(signal, max_length) [signal; zeros(max_length - length(signal), 1)];
[y1, y2, y3] = deal(append_zeros(y1, max_length), append_zeros(y2, max_length), append_zeros(y3, max_length));

fs = min([fs1, fs2, fs3]);

function plot_signal (signal, fs)
    n = length(signal);
    t = linspace(0, n/fs, n);
    figure;
    subplot(3, 1, 1);
    plot(t, signal);
    title('signal in time domain');
    xlabel('time');
    ylabel('x(t)');

    # plot spectrum
    freq_x = fftshift(fft(signal, n));
    f = linspace(-fs/2, fs/2, n);
    subplot(3, 1, 2);
    plot(f, abs(freq_x));
    title('signal in frequency domain');
    xlabel('frequency');
    ylabel('|X(W)|');

    subplot(3, 1, 3);
    plot(f, unwrap(angle(freq_x)));
    title('signal in frequency domain');
    xlabel('frequency');
    ylabel('Phase');
endfunction

function modulated_signal = modulate(signal)
    fc = 25000;
    carrier = cos(2*pi*fc*t);
    modulated_signal = x.*carrier';
endfunction


plot_signal(y1, fs1);
plot_signal(y2, fs2);
plot_signal(y3, fs3);





while(waitforbuttonpress()==0) pause(1) end