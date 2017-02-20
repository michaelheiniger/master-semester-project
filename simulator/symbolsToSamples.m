% Modified version of code from course "Software-defined radio: A hands-on course" from EPFL.

% SYMBOLSTOSAMPLES Produces the samples of the modulated signal
% Z = SYMBOLSTOSAMPLES(SYMBOLS, PULSE, USF) produces the samples of a
% modulated pulse train. The sampled pulse is given in vector H,
% and the symbols modulating the pulse are contained in vector Y.
% USF is the upsampling factor, i.e., the relationship between the
% sampling frequency Fs and the symbol (data) frequency Fd, or equivalently,
% the number of samples per symbol: USF=Fs/Fd.

function z = symbolsToSamples(symbols, pulse, USF)

symbols = symbols(:); % make sure it is a column vector

% Inserting USF-1 zeros between every two consecutive samples of y
symbols_up = upsample(symbols, USF); 

% Convolve upsampled symbol vector with the shaping pulse h
 z =  conv(symbols_up, pulse);
 z = transpose(z);
