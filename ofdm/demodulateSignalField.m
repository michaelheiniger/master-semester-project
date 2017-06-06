function [bits] = demodulateSignalField(signalField, systemConfig)
%DEMODULATESIGNALFIELD Summary of this function goes here
%   Detailed explanation goes here

sc = systemConfig;

signalFieldSymbols = fftshift(fft(signalField, sc.signalFieldLength));

% Remove guard bands
signalFieldSymbols = signalFieldSymbols(sc.numZerosTop+1:end-sc.numZerosBottom);

% Remove zero-frequency subcarrier if needed
if sc.zeroFreqSubcarrier
    signalFieldSymbols = [signalFieldSymbols(1:sc.numTotalCarriers/2-sc.numZerosTop,:);...
        signalFieldSymbols(2+sc.numTotalCarriers/2-sc.numZerosTop:end,:)];
end

signalFieldSymbols(real(signalFieldSymbols) >= 0) = 0;
signalFieldSymbols(real(signalFieldSymbols) < 0) = 1;
bits = signalFieldSymbols;

end

