#! /usr/bin/octave
clear all;
clc;
close all;

pkg load signal

# function to plot signal in time domain and its spectrum
function plot_signal(signal, fs)

    n = length(signal);
    t = linspace(0, n/fs, n);
    figure;
    subplot(3, 1, 1);
    plot(t, signal);
    title('signal in time domain');
    xlabel('time');
    ylabel('x(t)');

    # plot spectrum
    # obtain fourier transform 
    freq_x = fftshift(fft(signal));
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

# function to plot two signals on the same graph for comparison
function compare_signals(signal1, signal2, fs)

    n = length(signal1);
    t = linspace(0, n/fs, n);
    
    figure;
    subplot(3, 1, 1);
    plot(t, signal1, 'r', t, signal2, 'g');
    title('signal in time domain');
    xlabel('time');
    ylabel('x(t)');

    # plot spectrum
    # obtain fourier transform 
    freq_x1 = fftshift(fft(signal1));
    freq_x2 = fftshift(fft(signal2));
    f = linspace(-fs/2, fs/2, n);
    subplot(3, 1, 2);
    plot(f, abs(freq_x1), 'r', f, abs(freq_x2), 'g');
    title('signal in frequency domain');
    xlabel('frequency');
    ylabel('|X(W)|');

    subplot(3, 1, 3);
    plot(f, unwrap(angle(freq_x1)), 'r', f, unwrap(angle(freq_x2)), 'g');
    title('signal in frequency domain');
    xlabel('frequency');
    ylabel('Phase');

endfunction


# function to modulate the signal by multiplying it by cos with
# frequency fc
function modulated_signal = modulate(signal, fs, fc)
    t = linspace(0, length(signal)/fs, length(signal));
    carrier = cos(2*pi*fc*t);
    # s(t) = m(t) * cos(Wc)
    modulated_signal = signal .* carrier';
endfunction
# function to modulate the signal by multiplying it by the same
# carrier in the modulation step
# with ability to perform phase shift by deltaPhi and/or frequency
# shift by deltaFc
function demodulated_signal = demodulate(signal, fs, fc, deltaPhi = 0, deltaFc = 0)
    t = linspace(0, length(signal)/fs, length(signal));
    carrier = cos(2*pi*(fc + deltaFc) * t + deltaPhi);
    # s(t) = m(t) * cos(2pi (fc + delta)t + deltaPhi)
    demodulated_signal = signal .* carrier';
endfunction

# function to filter the passed signal by bandbass filter centered at fc
# and has a bandwidth BW
function y_filtered = bandPassFilter(y_modulated, order, BW, fc, fs)
    [b, a] = butter(order, [fc - BW/2, fc + BW/2]/(fs/2));
    y_filtered = filter(b, a, y_modulated) * 2;
endfunction 

#function to filter the passed signal by low pass filter with bancwidth BW
function baseband_demodulated = lpf(y_demodulated, order, BW, fs)
fc = BW/2;
[b, a] = butter(order, fc/(fs/2));
baseband_demodulated = filter(b, a, y_demodulated) *2 ;
endfunction

# read audio files
[y1, fs1] = audioread('audio1.wav');
[y2, fs2] = audioread('audio2.wav');
[y3, fs3] = audioread('audio3.wav');

# pick first channel of each signal if any has more than one channel
[y1, y2, y3] = deal(y1(:,1), y2(:,1), y3(:,1));

# resample signals with rate=4 to increase the frequency period to avoid overlapping
[y1, y2, y3] = deal(resample(y1, 4, 1), resample(y2, 4, 1), resample(y3, 4, 1));

# make signals have the same length = maximum length of the three
max_length = max([length(y1), length(y2), length(y3)]);
append_zeros = @(signal, max_length) [signal; zeros(max_length - length(signal), 1)];
[y1, y2, y3] = deal(append_zeros(y1, max_length), append_zeros(y2, max_length), append_zeros(y3, max_length));

# pick sampling frequency = fs * 4, all fsx are the same = 441000
# but just in case they are different, pick the minimum as reference
fs = 4 * min([fs1, fs2, fs3]);
plot_signal(y1, fs);
plot_signal(y2, fs);
plot_signal(y3, fs);

# all signals don't go beyong 4700
bandwidth = 4700;
# carrier frequency >> signal frequency to avoid overlapping
fc = 15000;

# modulation
y1_mod = modulate(y1, fs, fc);
y2_mod = modulate(y2, fs, fc + 2 * bandwidth);
y3_mod = modulate(y3, fs, fc + 4 * bandwidth);

# FDM system
y = y1_mod + y2_mod + y3_mod;

f = linspace(-fs/2, fs/2, length(y));
figure;
plot_signal(y, fs);

# filter modulated signals by bandpass filters
# first signal is centered at fc = fc
y1_filtered = bandPassFilter(y, 3, bandwidth, fc, fs);
# first signal is centered at fc = fc + 2* bandwidth
y2_filtered = bandPassFilter(y, 3, bandwidth, fc + 2 * bandwidth, fs);
# first signal is centered at fc = fc + 4 * bandwidth
y3_filtered = bandPassFilter(y, 3, bandwidth, fc + 4 * bandwidth, fs);

# plot results
plot_signal(y1_filtered, fs);
plot_signal(y2_filtered, fs);
plot_signal(y3_filtered, fs);

# demodulated signals
# first signal carrier had frequency  = fc
y1_demodulated = demodulate(y1_filtered, fs, fc);
# first signal carrier had frequency  = fc + 2* bandwidth
y2_demodulated = demodulate(y2_filtered, fs, fc + 2 * bandwidth);
# first signal carrier had frequency  = fc + 4 * bandwidth
y3_demodulated = demodulate(y3_filtered, fs, fc + 4 * bandwidth);

# plot results
plot_signal(y1_demodulated, fs);
plot_signal(y2_demodulated, fs);
plot_signal(y3_demodulated, fs);

# extract demodulated signals with low pass filter 
y1_lpf = lpf(y1_demodulated, 3, bandwidth, fs);
y2_lpf = lpf(y2_demodulated, 3, bandwidth, fs);
y3_lpf = lpf(y3_demodulated, 3, bandwidth, fs);

# compare original signals and demodulated
compare_signals(y1_lpf, y1, fs);
compare_signals(y2_lpf, y2, fs);
compare_signals(y3_lpf, y3, fs);

# phase shift by 10 degrees to signal 1
y1_demodulated_10 = demodulate(y1_filtered, fs, fc, deltaPhi = 10*pi/180);
# phase shift by 30 degrees to signal 1
y1_demodulated_30 = demodulate(y1_filtered, fs, fc, deltaPhi = 30*pi/180);
# phase shift by 90 degrees to signal 1
y1_demodulated_90 = demodulate(y1_filtered, fs, fc, deltaPhi = 90*pi/180);

# filter signals
y1_10_lpf = lpf(y1_demodulated_10, 3, bandwidth, fs);
y1_30_lpf = lpf(y1_demodulated_30, 3, bandwidth, fs);
y1_90_lpf = lpf(y1_demodulated_90, 3, bandwidth, fs);

# plot results
plot_signal(y1_10_lpf, fs);
plot_signal(y1_30_lpf, fs);
plot_signal(y1_90_lpf, fs);


# frequency shift by 2Hz
y1_demodulated_F_2 = demodulate(y1_filtered, fs, fc, deltaFc = 2);
# frequency shift by 10Hz
y1_demodulated_F_10 = demodulate(y1_filtered, fs, fc, deltaFc = 10);

# filter signals
y1_2_F_lpf = lpf(y1_demodulated_F_2, 3, bandwidth, fs);
y1_10_F_lpf = lpf(y1_demodulated_F_2, 3, bandwidth, fs);

# plot results
plot_signal(y1_2_F_lpf, fs);
plot_signal(y1_10_F_lpf, fs);
