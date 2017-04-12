function [fineTau] = timingEstimation(signal, shortTrainingSeqIfft)
%TIMINGESTIMATION provide timing estimate of the begining of the next
%received frame

global decoder;

% Coarse timing estimation
% Cross-correlation between clean short STS with the signal
[R, ~] = xcorr(signal, shortTrainingSeqIfft);
R = R(length(signal):end);
[~, coarseTau] = max(abs(R))

figure;
plot(abs(R)); % Normalization needed ??
title('Cross-correlation with short training sequence ');

% Fine timing estimation
% Search before and after the coarse timing estimate
coarseWindowLength = 145; 
if decoder.upsample
    coarseWindowLength = coarseWindowLength*decoder.USF;
end
fineWindow = (coarseTau-coarseWindowLength+1:coarseTau+coarseWindowLength-1);
innerProducts = zeros(1,length(fineWindow));
firstPreambleClean = repmat(shortTrainingSeqIfft, 10,1);
for d = 1:length(fineWindow)
    firstPreambleRx = signal(fineWindow(d):fineWindow(d)+length(firstPreambleClean)-1);
    innerProducts(d) = doInnerProduct(firstPreambleRx, firstPreambleClean);
end

[~, pos] = max(abs(innerProducts));
fineTau = fineWindow(pos)



figure;
plot(abs(innerProducts),'.')
title('Inner-products with long training sequence');

end

