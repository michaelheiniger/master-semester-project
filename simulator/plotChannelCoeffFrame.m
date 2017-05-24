function [] = plotChannelCoeffFrame(config, ofdmSymbolsRxCorrected, dataFrame, description)
%PLOTCHANNELCOEFFFRAME Summary of this function goes here
%   Detailed explanation goes here

if numel(dataFrame) == 0
    return
end

c = config;


% Plot of the channel coefficients for the whole frame (using division and
% NOT MMSE)
if c.dcSubcarrier
    usedSubcarriersIndices = [-c.numTotalCarriers/2+c.numZerosTop:-1, 1:c.numTotalCarriers/2-1-c.numZerosBottom];
    
    % Compute channel coefficients, DC subcarrier excluded (it is only zeros)
    lambdaAllSymbols = ofdmSymbolsRxCorrected./[dataFrame(1+c.numZerosTop:c.numTotalCarriers/2,:);dataFrame(2+c.numTotalCarriers/2:end-c.numZerosBottom,:)];
else    
    usedSubcarriersIndices = -c.numTotalCarriers/2+c.numZerosTop:c.numTotalCarriers/2-1-c.numZerosBottom;
    
    % Compute channel coefficients, DC subcarrier included
    lambdaAllSymbols = ofdmSymbolsRxCorrected./dataFrame(1+c.numZerosTop:end-c.numZerosBottom,:);
%     lambdaAllSymbols = ofdmSymbolsRxCorrected(:,end)./dataFrame(1+c.numZerosTop:end-c.numZerosBottom,end);
end

figure;
xAxis = usedSubcarriersIndices;
subplot(2,1,1),plot(xAxis, angle(lambdaAllSymbols),'.-');
ylabel('Phase [rad]');
xlabel('Subcarriers index');
title(description);
grid on;
axis([-c.numTotalCarriers/2-5,c.numTotalCarriers/2+5,-3.5,3.5]);
subplot(2,1,2),plot(xAxis, abs(lambdaAllSymbols),'.-');
ylabel('Magnitude');
xlabel('Subcarriers index');
% title(['Magnitude evolution over subcarriers for the whole frame ']);
grid on;
% axis([-35,35,0.8,1.1]);

end

