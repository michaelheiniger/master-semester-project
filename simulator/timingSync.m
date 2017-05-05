function [tau] = timingSync(signal, stsTime, ltsTime, twoLtsCpLength)
%TIMINGESTIMATION provide timing estimate of the begining of the next
%received frame

global decoder;

% Coarse timing estimation
% Cross-correlation between 10 clean shorts STS with the signal
% [R, ~] = xcorr(signal, stsTime);
[R, ~] = xcorr(signal, repmat(stsTime,10,1));
R = R(length(signal):end);
[~, coarseTau] = max(abs(R)); % Normalization needed !!

figure;
plot(abs(R));
xlabel('Offset [samples]');
ylabel('Magnitude');
title('Cross-correlation (signal, 10 STS)');

% % Fine timing estimation
% coarseWindowLength = 145; 
% if decoder.upsample
%     coarseWindowLength = coarseWindowLength*decoder.USF;
% end
% frameOffsetStartTwoLts = 10*length(stsTime) + twoLtsCpLength;
% fineWindow = (coarseTau-coarseWindowLength+1:16:coarseTau+coarseWindowLength-1);
% innerProducts = zeros(1,length(fineWindow));
% twoLtsTime = repmat(ltsTime, 2,1);
% for d = 1:length(fineWindow)
%     twoLtsTimeRx = signal(fineWindow(d)+frameOffsetStartTwoLts:fineWindow(d)+frameOffsetStartTwoLts+length(twoLtsTime)-1);
%     innerProducts(d) = doInnerProduct(twoLtsTimeRx, twoLtsTime);
% end
% 
% [~, pos] = max(abs(innerProducts));
% fineTau = fineWindow(pos)
% 
% figure;
% plot(abs(innerProducts),'.')
% title('Inner-products with long training sequence');

tau = coarseTau + decoder.timingOffset
% tau = fineTau;

end

