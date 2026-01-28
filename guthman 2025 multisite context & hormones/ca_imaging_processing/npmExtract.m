%extract data from NPM version FP3002
%created 08-02-22, eartha mae guthman phd

%% start script
close all
clearvars

%% inits
%ROI
% roi.labels = {'TTd','BNSTp','MPOA','AHA','LSv','VMHvl','lPAG','MeAp','ArcN','PMv','PA','X'}; %pre-230406
% roi.labels = {'TTd','BNSTp','MPOA','AHA','LSv','VMHvl','lPAG','MeAp','ArcN','PA','PMv','X'}; %post-230406/pre230601
% roi.labels = {'AHA','VMHvl','PA','TTd','PMv','LSv','ArcN','lPAG','MeAp','BNSTp','MPOA','X'};  %post-230601/pre240101
% roi.labels = {'AHA','MeAp','TTd','lPAG','MPOA','PA','ArcN','LSv','BNSTp','X','PMv','VMHvl'}; %post-240101/pre260121
roi.labels = {'MPOA', 'lPAG', 'TTd', 'MeAp', 'AHA', 'BNSTp', 'LSv', 'ArcN', 'PA', 'VMHvl', 'PMv', 'X'}; %post-260121
% roi.gcamp_col = {'J' 'L' 'N' 'P' 'R' 'T' 'V' 'X' 'Z' 'AB' 'AD' 'AF'}; %pre-260121
% roi.rcamp_col = {'I' 'K' 'M' 'O' 'Q' 'S' 'U' 'W' 'Y' 'AA' 'AC' 'AE'}; %pre-260121
roi.gcamp_col = {'F' 'H' 'J' 'L' 'N' 'P' 'R' 'T' 'V' 'X' 'Z' 'AB'}; %post-260121
roi.rcamp_col = {'E' 'G' 'I' 'K' 'M' 'O' 'Q' 'S' 'U' 'W' 'Y' 'AA'}; %post-260121

%colors
cs.rc = [160 113 120] ./ 255;
cs.gc = [127 183 190] ./ 255;
cs.isos = [200 182 255] ./ 255;

%save toggle
saveOn = 1; %save output

%directory and filenames
file_prefix = '260122_L89_home_balbc';
loadingDir = '/Volumes/falkner/Mae/multisite_hormone_experiments/data/sham_gdx/L89/260122/';
saveDir = '/Volumes/falkner/Mae/multisite_hormone_experiments/data/sham_gdx/L89/260122/';
daqDir = dir([loadingDir file_prefix]);
daqFile = daqDir(3).name;
caImgFull = [loadingDir file_prefix '_ca_imaging0.csv'];
daqFull = [loadingDir file_prefix '/' daqFile];

%other
daqFs = 350; %Hz
npmFs = 40;
lsq_opt = optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',1e4,'MaxIterations',1e4,'Display','off');

%% load data
disp(strjoin(['~ ' string(datetime) '~ ////loading data////'], ''))

%FP3002
[bonsai, timestamps.npm.sec, led_state.all] = extractNPM(caImgFull, roi);

%DAQ
daqData = readmatrix(daqFull);

%% partition vectors
[daq, timestamps, led_state, bonsai, lVideo, noFrames] = partitionVectorsFP3002(daqData, daqFs, timestamps, led_state, bonsai);

%% de-interleave
led_state = deinterleave_leds(led_state);
bonsai = deinterleave_signals(bonsai, led_state);
timestamps = deinterleave_timestamps(timestamps, led_state);

%equalize signal lengths
[bonsai, timestamps] = equalize_vector_lengths(bonsai, timestamps);

%% match npm and behavior video timestamps
indexing = behavior_npm_indexing_alignment(timestamps);

%% plot figures to check for frame assignment
frameAlignmentQC(timestamps, indexing, cs)

%% plots to compare de-interleaved traces
deinterleaveQC(bonsai, timestamps, roi.labels, cs)

%% process signal for display
% remove first ten seconds of data to delete light artifacts at start of experiment
dur = 10; %duration to cut in seconds
[dispSignal, dispTimestamps] = cutInitialDur(bonsai, timestamps, dur, npmFs);

%detrend (linear decay)
[detrendedRCaMP, detrendedGCaMP] = detrend_signals(dispSignal, roi.labels);

%dff -- percent change in fluorescence --> 1+(signal - signal_mean)/(signal_mean)
[dffRCaMP, dffGCaMP] = compute_dff(detrendedRCaMP, detrendedGCaMP);

%isosbestic correct fiber (regression using red/green isosbestic that best correlates with signal
[dffRCaMP, dffGCaMP] = isosCorrect(dffRCaMP, dffGCaMP);

%filter (10 Hz lowpass - ~1.67X GCaMP6f decay [100 ms vs 140 ms]) filter
cutoff = 10; %lowpass Freq
o = 4; %order
Rs = o * 20; %decibels of stopband attenuation
Ws = cutoff / (npmFs / 2); %normalized cutoff requency
[lowpassRCaMP, lowpassGCaMP] = filter_signal(dffRCaMP.corrected, dffGCaMP.corrected, o, Rs, Ws);

%zscore
zRCaMP = zeros(size(lowpassRCaMP)); zGCaMP = zeros(size(lowpassGCaMP));
for ii = 1:height(lowpassRCaMP)
    zRCaMP(ii, :) = zscore(lowpassRCaMP(ii, :));
    zGCaMP(ii, :) = zscore(lowpassGCaMP(ii, :));
end

%% plot output
%inits
data_stream = "zscore"; %dff or zscore
%plot
if strcmp(data_stream, "dff")
    plotInitiallyProcessedSignal(lowpassRCaMP, lowpassGCaMP, dispTimestamps, roi.labels, data_stream, cs)
elseif strcmp(data_stream, "zscore")
    plotInitiallyProcessedSignal(zRCaMP, zGCaMP, dispTimestamps, roi.labels, data_stream, cs)
end

%% save
%pause here to view figures
close all

%save
if saveOn
    try
        save([saveDir file_prefix '_ca_imaging.mat'])
    catch
        mkdir(loadingDir,'fp_data');
        save([saveDir file_prefix '_ca_imaging.mat'])
    end
    disp(strjoin(['~ ' string(datetime) '~ ////Ca2+ imaging data saved////'], ''))
end

%% subroutines
function leds = deinterleave_leds(leds)
leds.off = find(leds.all == 0);
leds.isos = find(leds.all == 1);
leds.blue = find(leds.all == 2);
leds.green = find(leds.all == 4);
end

function signal = deinterleave_signals(signal, leds)
signal.gcamp.raw = signal.signal470(:,leds.blue);
signal.rcamp.raw = signal.signal560(:,leds.green);
signal.gcamp.isos = signal.signal470(:,leds.isos);
signal.rcamp.isos = signal.signal560(:,leds.isos);
signal.gcamp.rcamp = signal.signal470(:,leds.green);
signal.rcamp.gcamp = signal.signal560(:,leds.blue);
end

function ts = deinterleave_timestamps(ts, leds)
ts.gcamp.sec = ts.npm.sec(leds.blue);
ts.rcamp.sec = ts.npm.sec(leds.green);
ts.isos.sec = ts.npm.sec(leds.isos);
ts.gcamp.ind = 1:length(ts.gcamp.sec);
ts.rcamp.ind = 1:length(ts.rcamp.sec);
ts.isos.ind = 1:length(ts.isos.sec);
end

function frameAlignmentQC(ts, idx, mapping)
%initialize figure
frameassignmentFig = figure;
frameassignmentFig.Position = [600 600 600 600];

%RCaMP to Behavior
subplot(2,2,1)
x_start = randperm(round(ts.rcamp.sec(end)) - 10, 1);
scatter(ts.rcamp.sec, idx.rcamp_to_behv(2, :), 60, mapping.rc, 'filled')
ylabel('behavior cam (s)')
xlabel('rcamp (s)')
xlim([x_start x_start+4]); ylim([x_start x_start+4]);
setAx(gca);

%GCaMP to Behavior
subplot(2,2,2)
x_start = randperm(round(ts.gcamp.sec(end)) - 10, 1);
scatter(ts.gcamp.sec,idx.gcamp_to_behv(2, :), 60, mapping.gc, 'filled')
xlim([x_start x_start+4]); ylim([x_start x_start+4]);
ylabel('behavior cam (s)')
xlabel('gcamp (s)')
setAx(gca);

%Isosbestic to Behavior
subplot(2,2,3)
x_start = randperm(round(ts.isos.sec(end)) - 10, 1);
scatter(ts.isos.sec,idx.isos_to_behv(2, :), 60, mapping.isos, 'filled')
ylabel('behavior cam (s)')
xlabel('isos (s)')
xlim([x_start x_start+4]); ylim([x_start x_start+4]);
setAx(gca);

%Behavior to Neural Activity
subplot(2,2,4)
x_start = randperm(round(ts.behavior_cam.sec(end)) - 10, 1);
scatter(ts.behavior_cam.sec, idx.behv_to_npm(2, :), 60, 'k', 'filled')
ylabel('rcamp (s)')
xlabel('behavior cam (s)')
xlim([x_start x_start+4]); ylim([x_start x_start+4]);
setAx(gca);

end

function deinterleaveQC(signal, ts, roi, mapping)
%initialize figure
diFig = figure;
diFig.Position = [280 280 3200 1800];
diFig.Color = [1 1 1];

for ii = 1:length(roi) %loop thru fibers
    %raw RCaMP
    subplot(length(roi), 4, 1 + (ii-1)*4) %red channel, 560nm
    plot(ts.rcamp.sec, signal.rcamp.raw(ii, :), 'color', mapping.rc)
    if ii == 1
        title('jRCaMP1b (raw au)', 'FontSize', 18)
    elseif ii == length(roi)
        xlabel('time (s)', 'FontSize', 14)
    end
    xlim([0 ts.rcamp.sec(end)])
    ylabel(roi{ii}, 'fontsize', 16, 'fontweight', 'bold')
    axA = gca;
    axA.LineWidth = 2; axA.XColor = 'k'; axA.YColor = 'k'; axA.Box = 'off'; axA.TickDir = "out";

    %isosbestic
    subplot(length(roi), 4, 2 + (ii-1)*4) %red channel, 415nm
    plot(ts.isos.sec, signal.rcamp.isos(ii, :),'Color','k')
    if ii == 1
        title('isosbestic signal (415 nm)', 'FontSize', 18)
    elseif ii == length(roi)
        xlabel('time (s)', 'FontSize', 14)
    end
    xlim([0 ts.isos.sec(end)])
    axC = gca;
    axC.LineWidth = 2; axC.XColor = 'k'; axC.YColor = 'k'; axC.Box = 'off'; axC.TickDir = "out";

    %raw GCaMP
    subplot(length(roi), 4, 3 + (ii-1)*4) %green channel, 470nm
    plot(ts.gcamp.sec, signal.gcamp.raw(ii, :), 'color', mapping.gc)
    if ii == 1
        title('GCaMP6f (raw au)', 'FontSize', 18)
    elseif ii == length(roi)
        xlabel('time (s)', 'FontSize', 14)
    end
    xlim([0 ts.gcamp.sec(end)])
    axB = gca;
    axB.LineWidth = 2; axB.XColor = 'k'; axB.YColor = 'k'; axB.Box = 'off'; axB.TickDir = "out";

    %isosbestic
    subplot(length(roi), 4, 4 + (ii-1)*4) %green channel, 415nm
    plot(ts.isos.sec, signal.gcamp.isos(ii,:), 'Color', 'k')
    if ii == 1
        title('isosbestic signal (415 nm)', 'FontSize', 18)
    elseif ii == length(roi)
        xlabel('time (s)', 'FontSize', 14)
    end
    xlim([0 ts.isos.sec(end)])
    axC = gca;
    axC.LineWidth = 2; axC.XColor = 'k'; axC.YColor = 'k'; axC.Box = 'off'; axC.TickDir = "out";
end
end

function [tmp_signal, tmp_ts] = cutInitialDur(signal, ts, dur, fs)
%signal
tmp_signal.rcamp.raw = signal.rcamp.raw;
tmp_signal.rcamp.raw(:, 1:(dur*fs)) = [];
tmp_signal.rcamp.isos = signal.rcamp.isos;
tmp_signal.rcamp.isos(:, 1:(dur*fs)) = [];
tmp_signal.rcamp.gcamp = signal.rcamp.gcamp;
tmp_signal.rcamp.gcamp(:, 1:(dur*fs)) = [];
tmp_signal.gcamp.raw = signal.gcamp.raw;
tmp_signal.gcamp.raw(:, 1:(dur*fs)) = [];
tmp_signal.gcamp.isos = signal.gcamp.isos;
tmp_signal.gcamp.isos(:, 1:(dur*fs)) = [];
tmp_signal.gcamp.rcamp = signal.gcamp.rcamp;
tmp_signal.gcamp.rcamp(:, 1:(dur*fs)) = [];

%timestamps
tmp_ts.rcamp.sec = ts.rcamp.sec;
tmp_ts.rcamp.sec(1:(dur*fs)) = [];
tmp_ts.rcamp.ind = ts.rcamp.ind;
tmp_ts.rcamp.ind((end - (dur*fs - 1)):end) = [];
tmp_ts.gcamp.sec = ts.gcamp.sec;
tmp_ts.gcamp.sec(1:(dur*fs)) = [];
tmp_ts.gcamp.ind = ts.gcamp.ind;
tmp_ts.gcamp.ind((end - (dur*fs - 1)):end) = [];
tmp_ts.isos.sec = ts.isos.sec;
tmp_ts.isos.sec(1:(dur*fs)) = [];
tmp_ts.isos.ind = ts.isos.ind;
tmp_ts.isos.ind((end - (dur*fs - 1)):end) = [];
end

function [detrendedRCaMP, detrendedGCaMP] = detrend_signals(bonsai, roi)
for ii = 1:length(roi)
    detrendedRCaMP.signal(ii, :) = detrend(bonsai.rcamp.raw(ii, :)) + mean(bonsai.rcamp.raw(ii, :));
    detrendedRCaMP.isos(ii, :) = detrend(bonsai.rcamp.isos(ii, :)) + mean(bonsai.rcamp.isos(ii, :));
    detrendedGCaMP.signal(ii, :) = detrend(bonsai.gcamp.raw(ii, :)) + mean(bonsai.gcamp.raw(ii, :));
    detrendedGCaMP.isos(ii, :) = detrend(bonsai.gcamp.isos(ii,: )) + mean(bonsai.gcamp.isos(ii, :));
end
end

function [ddf0, dff1] = compute_dff(signal0, signal1)
for ii = 1:height(signal0.signal)
    ddf0.signal(ii, :) = 1 + (signal0.signal(ii, :) - mean(signal0.signal(ii, :))) ./ mean(signal0.signal(ii, :));
    ddf0.isos(ii, :) = 1 + (signal0.isos(ii, :) - mean(signal0.isos(ii, :))) ./ mean(signal0.isos(ii, :));
    dff1.signal(ii, :) = 1 + (signal1.signal(ii, :) - mean(signal1.signal(ii, :))) ./ mean(signal1.signal(ii, :));
    dff1.isos(ii, :) = 1 + (signal1.isos(ii, :) - mean(signal1.isos(ii, :))) ./ mean(signal1.isos(ii, :));
end
end

function [dff0, dff1] = isosCorrect(dff0, dff1)
for ii = 1:height(dff0.signal)
    %ID isosbestic to use (the one with most correlation to the signal
    s0_mat = [dff0.signal(ii, :)' dff0.isos(ii, :)' dff1.isos(ii, :)'];
    s1_mat = [dff1.signal(ii, :)' dff0.isos(ii, :)' dff1.isos(ii, :)'];
    coeff0 = corrcoef(s0_mat); coeff1 = corrcoef(s1_mat);
    if coeff0(1, 2) >= coeff0(1, 3)
        isos0 = dff0.isos(ii, :);
    else
        isos0 = dff1.isos(ii, :);
    end
    if coeff1(1, 2) > coeff0(1, 3)
        isos1 = dff0.isos(ii, :);
    else
        isos1 = dff1.isos(ii, :);
    end

    %isosbestic correct signal
    regrCoeff0 = robustfit(isos0', dff0.signal(ii, :)'); 
    dff0.corrected(ii, :) = dff0.signal(ii, :) - (regrCoeff0(1) + isos0 .* regrCoeff0(2));
    regrCoeff1 = robustfit(isos1', dff1.signal(ii, :)'); 
    dff1.corrected(ii, :) = dff1.signal(ii, :) - (regrCoeff1(1) + isos1 .* regrCoeff1(2));
end
end

function [lowpass0, lowpass1] = filter_signal(dff0, dff1, o, Rs, Ws)
%pick filter type
[b, a] = cheby2(o, Rs, Ws, 'low'); %chebyshev, ripple cutoff, lowpass

%filter
lowpass0 = zeros(size(dff0)); lowpass1 = zeros(size(dff1));
for ii = 1:height(dff0)
    lowpass0(ii, :) = filtfilt(b, a, dff0(ii, :));
    lowpass1(ii, :) = filtfilt(b, a, dff1(ii, :));
end
end

function plotInitiallyProcessedSignal(rcamp, gcamp, ts, roi, ds, mapping)
%init figure
figure('Position', [300 300 3200 900], 'Color', [1 1 1]);

%plot
for ii = 1:height(rcamp)
    subplot(3, 4, ii)
    hold on

    %plot signal
    plot(ts.gcamp.sec, gcamp(ii, :), 'color',  mapping.gc, 'linewidth',0.67)
    plot(ts.rcamp.sec, rcamp(ii, :), 'color', mapping.rc,  'linewidth',0.67)
    
    %aesthetics
    title([roi{ii} ' signal'],'FontSize',24)
    xlim([0 ts.rcamp.sec(end)])
    xlabel('time (s)', 'FontSize', 18)
    if strcmp(ds, "dff")
        ylabel('dF/F', 'fontsize', 18)
    elseif strcmp(ds, "zscore")
        ylabel('Z(dF/F)', 'fontsize', 18)
    end
    axA = gca;
    setAx(axA);
end
end