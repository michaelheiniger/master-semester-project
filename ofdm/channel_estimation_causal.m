function [channel_est] = channel_estimation_causal(s_comp, training)
% Credits to Can Kanbak, EPFL

% Estimates the channel from the received and original training symbols by 
% solving a linear equation. The output channel is causal and the main
% tap is at the first tap.
%
% INPUT
% s_comp        = received training symbols
% training      = original training symbols
%
% OUTPUT
% channel_est   = estimated channel impulse response

L = floor((length(training)+1)/2);
training = training(:).';

X = toeplitz(training(L:2*L-1),training(L:-1:1));
y = s_comp(L:2*L-1);
% size(X)
% size(y)

channel_est = X\y;
end