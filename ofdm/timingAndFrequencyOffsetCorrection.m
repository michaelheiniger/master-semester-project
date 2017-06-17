function [frameCorrected, timingEst] = timingAndFrequencyOffsetCorrection(systemConfig, receiverConfig, coarseFrame, stsTime, ltsTime)
%TIMINGANDFREQUENCYOFFSET Summary of this function goes here
%   Detailed explanation goes here

sc = systemConfig;
rc = receiverConfig;

Ts = 1/sc.Fs;

% Timing and Frequency Offset Estimation
% - timingEst is the estimation of the index of the first sample of the frame
% - frequencyOffsetEst is the frequency offset the time-dependent complex exponential introducing
%   phase-shift drifting over time due to carrier frequency offset between
%   transmitter and receiver.
% - intialPhaseEst is the constant phase-shift part of the frequency offset
%   The point of reference for the inital phase is the start of the frame
% Here are several methods to compare, see the MATLAB functions header for
% more details
switch rc.timingAndFrequencyOffsetMethod
    case 'stsLtsOfdmDemod'
        [timingEst, frequencyOffsetEst, initialPhaseEst] = timingAndFrequencyOffsetEstimation1(sc, rc, coarseFrame, stsTime, ltsTime);
    case 'ideal'
        timingEst = rc.manualTiming;
        frequencyOffsetEst = rc.manualCFO;
        initialPhaseEst = angle(exp(1j*2*pi*frequencyOffsetEst*(timingEst-1)*Ts));
    otherwise
        error('Timing and frequency offset correction method unknown.');
end

disp(['Timing: ', num2str(timingEst)]);
disp(['Total CFO: ', num2str(frequencyOffsetEst)]);
disp(['Initial phase: ', num2str(initialPhaseEst)]);

% Allows to introduce manually a timing offset (see Receiver configuration)
beginningFrame = timingEst + rc.timingOffset;
disp(['Frame beginning (effectively used) :', num2str(beginningFrame)]);

% Frame extraction
frame = coarseFrame(beginningFrame:beginningFrame+sc.ofdmFrameLength-1);

frameCorrected = frame;

% Carrier frequency offset correction
if rc.cfoCorrection
    % Remove frequency offset and intial phase of the whole frame
    t = (0:length(frameCorrected)-1)*Ts;
    frameCorrected = frameCorrected .* exp(-1j*2*pi*frequencyOffsetEst*t).';
    frameCorrected = frameCorrected * exp(-1j*initialPhaseEst).';
end
    
end