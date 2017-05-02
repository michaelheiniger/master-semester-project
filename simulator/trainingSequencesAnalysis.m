clear all;
close all;
clc;


NFFT = 64;
Fs = 0.5e6;
Ts = 1/Fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Short Training Sequence (STS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('stsFreq.mat');

fourStsTime = ifft(fftshift(stsFreq), NFFT);
stsTime = fourStsTime(1:16);

tenStsTime = repmat(stsTime, 10, 1);

% STS in FREQUENCY domain
f = linspace(-Fs/2, Fs/2-Fs/NFFT, NFFT)/1000; % kHz

figure;
subplot(1,3,1),plot(f, abs(stsFreq));
xlabel('f [kHz]');
ylabel('|STS_F|(f)');
subplot(1,3,2),plot(f, real(stsFreq));
xlabel('f [kHz]');
ylabel('Re(STS_F)(f)');
subplot(1,3,3),plot(f, imag(stsFreq));
xlabel('f [kHz]');
ylabel('Im(STS_F)(f)');

% STS in TIME domain
t = (0:length(stsTime)-1)*Ts*1e6;

figure;
subplot(1,3,1),plot(t, abs(stsTime));
xlabel('t [\mus]');
ylabel('|STS|(t)');
subplot(1,3,2),plot(t, real(stsTime));
xlabel('t [\mus]');
ylabel('Re(STS)(t)');
subplot(1,3,3),plot(t, imag(stsTime));
xlabel('t [\mus]');
ylabel('Im(STS)(t)');

% STS auto-correlation in time domain
figure;
plot(abs(xcorr(stsTime, stsTime)));
ylabel('|xcorr(STS, STS)|')
title('Auto-correlation of STS in time domain');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Long Training Sequence (LTS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('ltsFreq.mat');

ltsTime = ifft(fftshift(ltsFreq), NFFT);

% LTS in FREQUENCY domain
f = linspace(-Fs/2, Fs/2-Fs/NFFT, NFFT)/1000; % kHz

figure;
subplot(1,3,1),plot(f, abs(ltsFreq));
xlabel('f [kHz]');
ylabel('|LTS_F|(f)');
subplot(1,3,2),plot(f, real(ltsFreq));
xlabel('f [kHz]');
ylabel('Re(LTS_F)(f)');
subplot(1,3,3),plot(f, imag(ltsFreq));
xlabel('f [kHz]');
ylabel('Im(LTS_F)(f)');

% LTS in TIME domain
t = (0:length(ltsTime)-1)*Ts*1e6;

figure;
subplot(1,3,1),plot(t, abs(ltsTime));
xlabel('t [\mus]');
ylabel('|LTS|(t)');
subplot(1,3,2),plot(t, real(ltsTime));
xlabel('t [\mus]');
ylabel('Re(LTS)(t)');
subplot(1,3,3),plot(t, imag(ltsTime));
xlabel('t [\mus]');
ylabel('Im(LTS)(t)');

% STS auto-correlation in time domain
figure;
plot(abs(xcorr(ltsTime, ltsTime)));
ylabel('|xcorr(LTS, LTS)|')
title('Auto-correlation of LTS in time domain');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full preamble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

preambleTime = [tenStsTime; ltsTime(end-31:end); ltsTime; ltsTime];

t = (0:length(preambleTime)-1)*Ts*1e6;
figure;
subplot(1,3,1),plot(t, abs(preambleTime));
xlabel('t [\mus]');
ylabel('|Preamble|(t)');
subplot(1,3,2),plot(t, real(preambleTime));
xlabel('t [\mus]');
ylabel('Re(Preamble)(t)');
subplot(1,3,3),plot(t, imag(preambleTime));
xlabel('t [\mus]');
ylabel('Im(Preamble)(t)');


