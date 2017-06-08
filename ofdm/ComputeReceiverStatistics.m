close all;
clear all;
clc;

% Concatenate all receiver1 and receiver2 results from all batches:
% a = ones(400,6,23);
% b = ones(400,6,216);
% c = cat(3,a,b);
load('rxPerfData/receiver1Results_batch3.mat');
receiver1ConcatResults = receiver1Results;
clear receiver1Results;

load('rxPerfData/receiver2Results_batch3.mat');
receiver2ConcatResults = receiver2Results;
clear receiver2Results;
for i = 4:20 % Note: batches 1 and 2 don't have the same format, so not used
    load(['rxPerfData/receiver1Results_batch' num2str(i) '.mat']);
    receiver1ConcatResults = cat(1,receiver1ConcatResults, receiver1Results);
    clear receiver1Results;

    load(['rxPerfData/receiver2Results_batch' num2str(i) '.mat']);
    receiver2ConcatResults = cat(1,receiver2ConcatResults, receiver2Results);
    clear receiver2Results;
end

% Concatenate all receiver1 timings
load('rxPerfData/timingResultFileReceiver1_batch3.mat');
timingsReceiver1Concat = receiver1Timings;
clear receiver1Timings;
for i = 4:20 % Note: batches 1 and 2 don't exist
    load(['rxPerfData/timingResultFileReceiver1_batch' num2str(i) '.mat']);
    timingsReceiver1Concat = cat(1,timingsReceiver1Concat, receiver1Timings);
    clear receiver1Timings;
end

% Compute the timing offset error
% The correct timing is 1001
timingsReceiver1Concat = timingsReceiver1Concat - 1001;

%% Compute mean over main runs for SER and timings
meanOverData1 = mean(receiver1ConcatResults,1);
meanOverData2 = mean(receiver2ConcatResults,1);
meanTimingOverData1 = mean(timingsReceiver1Concat,1);

% Compute variance over SNR runs (i.e. variance of noise and multipath)
varResultsReceiver1 = var(meanOverData1,0,3);
varResultsReceiver2 = var(meanOverData2,0,3);
varTimingsReceiver1 = var(meanTimingOverData1,0,3);

% Compute mean over SNR runs
meanOverDataThenSnr1 = mean(meanOverData1, 3);
meanOverDataThenSnr2 = mean(meanOverData2, 3);
meanTimingOverDataThenSnr1 = mean(meanTimingOverData1, 3);

% Compute min/max timings per SNR 
maxOffsetPerSNr = max(max(timingsReceiver1Concat, [],3),[],1);
minOffsetPerSNr = min(min(timingsReceiver1Concat, [],3),[],1);

%%
clc
size(timingsReceiver1Concat)
% allTimingsPerSnr = reshape(transpose(reshape(timingsReceiver1Concat,25,[])),36,[]); % FAUX !!???!!!
% allTimingsPerSnr = transpose(reshape(transpose(reshape(timingsReceiver1Concat,18,[])), 25*18,[]));
allTimingsPerSnr = transpose(reshape(transpose(reshape(timingsReceiver1Concat, 18*36, 25)),18*25,[]));
% size(allTimingsPerSnr)
a = allTimingsPerSnr(1,:);
b = a(a >=-4)
c = b(b < 1);
numel(c)/numel(a)

%%

snrValues = 20:55;

figure;
errorbar(snrValues,meanOverDataThenSnr1,varResultsReceiver1)
hold on;
errorbar(snrValues,meanOverDataThenSnr2,varResultsReceiver2)
% axis([15,55]);
xlabel('SNR [dB]');
ylabel('Mean SER');
legend('Actual Rx','Ideal Rx');
title('Mean SER vs SNR');

figure;
hax=axes;
errorbar(snrValues,meanTimingOverDataThenSnr1,varTimingsReceiver1)
hold on;
plot(snrValues, maxOffsetPerSNr,'.');
plot(snrValues, minOffsetPerSNr,'.');
xlabel('SNR [dB]');
ylabel('Timing offset error');
HL1 = 0;
line(get(hax,'XLim'),[HL1 HL1], 'Color', [0 1 0]);
legend('Mean timing', 'Max timing', 'Min timing');
title('Timing vs SNR');

% figure;
% plot(snrValues(1), allTimingsPerSnr(1,:),'.')
% hold on;
% for i = 2:25
%     plot(snrValues(i), allTimingsPerSnr(i,:),'.')
% end
% axis([15,50,-20,20]);
%% 
close all;
figure;
hax=axes;
histogram(allTimingsPerSnr(36,:));
% histogram(timingsReceiver1Concat(4,25,:));
hold on;
HL1 = -4.5;
HL2 = 0.5;
line([HL1 HL1], get(hax,'YLim'), 'Color', [0 1 0]);
line([HL2 HL2], get(hax,'YLim'), 'Color', [0 1 0]);

% histogram(timingsReceiver1Concat(4,1,:));

%%
figure;
plot(allTimingsPerSnr(1,:), snrValues(1),'.')


