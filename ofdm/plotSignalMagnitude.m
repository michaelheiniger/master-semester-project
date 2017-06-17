function [] = plotSignalMagnitude(signal, xLabel, plotTitle, varargin)
%PLOTSIGNALMAGNITUDE Summary of this function goes here
%   Detailed explanation goes here

if length(varargin) == 3
    startVertBar = varargin{1};
    endVertBar = varargin{2};
    barColor = varargin{3};
end

figure;
hax=axes;
plot(abs(signal));
xlabel(xLabel);
ylabel('Magnitude');
title(plotTitle);

if exist ('startVertBar', 'var') && exist ('endVertBar', 'var')
    switch barColor
        case 'red'
            color = [1 0 0];
        case 'green'
            color = [0 1 0];
        case 'blue'
            color = [0 0 1];
    end
    VL1 = startVertBar;
    line([VL1 VL1],get(hax,'YLim'), 'Color', color);
    VL2 = endVertBar;
    line([VL2 VL2],get(hax,'YLim'), 'Color', color);
end

end

