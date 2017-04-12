function [intialPhaseEst, totalfrequencyOffsetEst] = frequencyOffsetEstimation(frame, shortTrainingSeqIfft, longTrainingSeqIfft)
%FREQUENCYOFFSETESTIMATION Summary of this function goes here
%   Detailed explanation goes here

global decoder;

Fs = 0.5e6;
Ts = 1/Fs;

preamble2CPLength = 32;

lengthSTSIfft = length(shortTrainingSeqIfft);
lengthLTSIfft = length(longTrainingSeqIfft);

firstPreambleRx = transpose(reshape(frame(1:10*lengthSTSIfft), lengthSTSIfft, []));

% Coarse frequency offset estimation
coarseFrequencyOffset = 1/(2*pi*lengthSTSIfft*Ts)*angle(...
                        doInnerProduct(firstPreambleRx(2,:), firstPreambleRx(1,:))...
                        + doInnerProduct(firstPreambleRx(3,:), firstPreambleRx(2,:))...
                        + doInnerProduct(firstPreambleRx(4,:), firstPreambleRx(3,:))...
                        + doInnerProduct(firstPreambleRx(5,:), firstPreambleRx(4,:))...
                        + doInnerProduct(firstPreambleRx(6,:), firstPreambleRx(5,:))...
                        + doInnerProduct(firstPreambleRx(7,:), firstPreambleRx(6,:))...
                        + doInnerProduct(firstPreambleRx(8,:), firstPreambleRx(7,:))...
                        + doInnerProduct(firstPreambleRx(9,:), firstPreambleRx(8,:)))

% Extract second preamble (prefix excluded)
startPreamble2 = 1 + numel(firstPreambleRx) + preamble2CPLength;
endPreamble2 = startPreamble2+2*lengthLTSIfft-1;
preamble2Rx = frame(startPreamble2:endPreamble2);

% Correct for Doppler using coarse estimate
t = (0:length(preamble2Rx)-1)*Ts;
preamble2RxCorrected = preamble2Rx .* exp(-1j*2*pi*coarseFrequencyOffset*t).';

% Fine frequency offset estimation
startT0 = 1;
endT0 = startT0+lengthLTSIfft-1;
startT1 = endT0+1;
endT1 = startT1+lengthLTSIfft-1;
fineFrequencyOffset = 1/(2*pi*lengthLTSIfft*Ts)*angle(doInnerProduct(preamble2RxCorrected(startT1:endT1), preamble2RxCorrected(startT0:endT0)))

% Correct for Doppler using fine estimate
t = (0:length(preamble2RxCorrected)-1)*Ts;
preamble2RxCorrected = preamble2RxCorrected .* exp(-1j*2*pi*fineFrequencyOffset*t).';

% Initial phase estimate (inner product with clean LTS)
secondPreambleClean = [longTrainingSeqIfft; longTrainingSeqIfft];
intialPhaseEst = angle(doInnerProduct(preamble2RxCorrected, secondPreambleClean))

% Total frequency offset correction
totalfrequencyOffsetEst = coarseFrequencyOffset + fineFrequencyOffset

end