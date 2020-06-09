function filtData = channelFilt(originalData,order,normFreq,plotOn,name)
% Apply non-causal low passfilter to logged channel
% Plot results if desired
% Optional
if ~exist('plotOn','var') || isempty(plotOn)
    plotOn = 0;
end
if ~exist('name','var') || isempty(name)
    name = '';
end
[b,a] = butter(order,normFreq,'low'); % IIR filter design
filtData = filtfilt(b,a,originalData); % zero-phase filtering

if plotOn
    figure
    plot(originalData)
    hold on
    plot(filtData)
    legend('Original', 'Filtered')
    xlabel('Sample Num')
    ylabel('State')
    title(string(name))
end