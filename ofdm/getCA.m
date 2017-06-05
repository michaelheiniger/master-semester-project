function [ca] = getCA()
%GETCA Summary of this function goes here
%   Detailed explanation goes here

% C/A codes (from GPS)
load('training-sequences/CAcodes.mat');

% C/A code of satelite 1 in time domain
ca = satCAcodes(1,:);

end

