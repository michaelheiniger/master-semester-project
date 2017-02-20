% Modified version of code from course "Software-defined radio: A hands-on course" from EPFL.

% SUFFICIENTSTATISTICS Process the output of the channel to generate
%  sufficient statistics about the transmitted symbols.
% X = SUFFICIENTSTATISTICS(R, PULSE, USF) produces sufficient statistics
% about the transmitted symbols, given the signal received in vector R,
% the impulse response PULSE of the basic pulse (transmitting filter), and
% the integer USF, which is the ratio between the sampling rate and the
% symbol rate (upsampling factor)

function x = sufficientStatistics(r, pulse, USF)

% here we assume h is a row vector

pulse_matched = conj(fliplr(pulse));
% in our case h_matched = h, since h is real and symmetric

% the matched filter output before downsampling is
y = conv(r, pulse_matched);

% Picking the desired matched filter outputs:
% Comment: the matched filter output (before downsampling) is the correlation
% between the received signal and h. The first useful sample is the 
% length(h) th correlation result. At the end there are length(h)-1 unused
% correlation results. 
   
x = y(length(pulse):USF:end-length(pulse)+1);


