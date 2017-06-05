clear all;
close all;
clc;

load('training-sequences/ltsFreq.mat');

ltsTime = ifft(fftshift(ltsFreq),64);

cfo = 456;

Ts = 1/(0.5e6);
% Add CFO
t = (0:2*length(ltsTime)-1)*Ts;
signalCFO = [ltsTime; ltsTime] .* exp(1j*2*pi*cfo*t).';

ltsTime1CFO = signalCFO(1:64);
ltsTime2CFO = signalCFO(65:end);

cfoEst = angle(doInnerProduct(ltsTime2CFO, ltsTime1CFO))/(2*pi*64*Ts)

ts = (0:length(ltsTime)-1)*Ts;
ltsTime1Corr = ltsTime1CFO  .* exp(-1j*2*pi*cfoEst*ts).';
ltsTime2Corr = ltsTime2CFO .* exp(-1j*2*pi*cfoEst*ts).';

% ltsTime1Corr = signalCorr(1:64);
% ltsTime2Corr = signalCorr(65:end);

angle(doInnerProduct(ltsTime2Corr, ltsTime1Corr))/(2*pi*64*Ts)





