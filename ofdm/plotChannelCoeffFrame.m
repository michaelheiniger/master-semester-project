function [] = plotChannelCoeffFrame(systemConfig, ofdmSymbolsRxCorrected, dataFrame, description)
%PLOTCHANNELCOEFFFRAME Summary of this function goes here
%   Detailed explanation goes here

if numel(dataFrame) == 0
    return
end

sc = systemConfig;


% Plot of the channel coefficients for the whole frame (using division and
% NOT MMSE)
if sc.zeroFreqSubcarrier
    usedSubcarriersIndices = [-sc.numTotalCarriers/2+sc.numZerosTop:-1, 1:sc.numTotalCarriers/2-1-sc.numZerosBottom];
    
    % Compute channel coefficients, DC subcarrier excluded (it is only zeros)
    lambdaAllSymbols = ofdmSymbolsRxCorrected./[dataFrame(1+sc.numZerosTop:sc.numTotalCarriers/2,:);dataFrame(2+sc.numTotalCarriers/2:end-sc.numZerosBottom,:)];
else    
    usedSubcarriersIndices = -sc.numTotalCarriers/2+sc.numZerosTop:sc.numTotalCarriers/2-1-sc.numZerosBottom;
    
    % Compute channel coefficients, DC subcarrier included
    lambdaAllSymbols = ofdmSymbolsRxCorrected./dataFrame(1+sc.numZerosTop:end-sc.numZerosBottom,:);
%     lambdaAllSymbols = ofdmSymbolsRxCorrected(:,end)./dataFrame(1+c.numZerosTop:end-c.numZerosBottom,end);
end

figure;
xAxis = usedSubcarriersIndices;
subplot(2,1,1),plot(xAxis, angle(lambdaAllSymbols),'.-');
ylabel('Phase [rad]');
xlabel('Subcarriers index');
title(description);
grid on;
axis([-sc.numTotalCarriers/2-5,sc.numTotalCarriers/2+5,-3.5,3.5]);
subplot(2,1,2),plot(xAxis, abs(lambdaAllSymbols),'.-');
ylabel('Magnitude');
xlabel('Subcarriers index');
% title(['Magnitude evolution over subcarriers for the whole frame ']);
grid on;
% axis([-35,35,0,1.5]);

end

