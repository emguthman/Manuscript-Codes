%%function to load raw data from NPM FPv3002
function [bonsai, ts, leds] = extractNPM(fName, ROI)

%bonsai csv fiber imaging data
for ii = 1:length(ROI.labels)
    bonsai.signal470(ii,:) = readmatrix(fName, 'Range', [ROI.gcamp_col{ii} ':' ROI.gcamp_col{ii}], 'NumHeaderLines',1);
    bonsai.signal560(ii,:) = readmatrix(fName, 'Range', [ROI.rcamp_col{ii} ':' ROI.rcamp_col{ii}], 'NumHeaderLines',1);
end

%time stamps + led states
ts = readmatrix(fName, 'Range', 'B:B', 'NumHeaderLines',1);
leds = readmatrix(fName, 'Range', 'C:C', 'NumHeaderLines',1);

%check for second header row and remove
if isnan(bonsai.signal470(1, 1)) || isnan(bonsai.signal560(1, 1))
    bonsai.signal470(:, 1) = [];
    bonsai.signal560(:, 1) = [];
    ts(1) = [];
    leds(1) = [];
end
end