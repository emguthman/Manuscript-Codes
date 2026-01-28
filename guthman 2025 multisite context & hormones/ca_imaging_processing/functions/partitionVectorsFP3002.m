%%function to partition vectors for daq
function [daq, ts, leds, bonsai, vidlength, nFrames, fs] = partitionVectorsFP3002(data, Fs, ts, leds, bonsai)
%daq vectors
daq.t = data(:, 1) ./ Fs; %time in second (dp divided by sampling rate of DAQ//time does not import)
daq.ttl = data(:, 3); %ttl pulse (V)
daq.npmVolt = data(:, 4); %voltage from npm box (V)
daq.flirVolt = data(:, 6); %voltage from FLIR cam

%TTL trigger time
[~, daq.ttlOn] =  max(diff(daq.ttl));
daq.ttlOn = daq.ttlOn + 1; %time ttl on in data point
daq.ttlOnSec = daq.t(daq.ttlOn); %time ttl on in sec
[~, daq.ttlOff] =  min(diff(daq.ttl));
daq.ttlOff = daq.ttlOff + 1; %time ttl on in data point
daq.ttlOffSec = daq.t(daq.ttlOff); %time ttl on in sec

%NPM vectors
[~, daq.npmInd] =  findpeaks(diff(daq.npmVolt));
daq.npmOn = daq.npmInd(1) + 1;
daq.npmOnSec = daq.t(daq.npmOn); %time npm on in sec

%ID NPM frames and adjust NPM vectors if needed
firstFrameNpm = find(diff(ts.npm.sec) > .01) + 1; %first npm frame (if nec., cuts initial 2-3 frames captured by NPM prior to TTL by finding the pause in imaging)
bonsai.signal470(:,1:firstFrameNpm-1) = []; %remove data before first frame, if any
bonsai.signal560(:,1:firstFrameNpm-1) = [];
ts.npm.sec(1:firstFrameNpm-1) = [];
leds.all(1:firstFrameNpm-1) = [];
ts.npm.sec = ts.npm.sec - ts.npm.sec(1) + daq.npmOnSec - daq.ttlOnSec; %readjusts NPM ts to daq
ts.npm.ind = 1:length(ts.npm.sec);
vidlength.npm.sec = ts.npm.sec(end) - ts.npm.sec(1); %get npm video stats
nFrames.npm = length(ts.npm.sec);
fs.npm = nFrames.npm/vidlength.npm.sec;

%ID behavior camera frames and build vectors
[daq.flir.peaks, daq.flir.ind] = findpeaks(daq.flirVolt);
daq.flir.ind(daq.flir.peaks < 2) = [];
ts.behavior_cam.sec = daq.t(daq.flir.ind) - daq.ttlOnSec;
ts.behavior_cam.ind = 1:length(daq.flir.ind);
vidlength.behavior_cam.sec = ts.behavior_cam.sec(end) - ts.behavior_cam.sec(1);
nFrames.behavior_cam = length(daq.flir.ind);
fs.behavior_cam = nFrames.behavior_cam / vidlength.behavior_cam.sec;

%plot frame rate over time to check for dropped frames
plotFR(ts)

%output stats
out_stats(fs, nFrames, vidlength)

end

%% subroutines
function plotFR(ts)
%initialize
figure('Position', [-2475 650 275 325], 'Color', [1 1 1]);
hold on;
plot(ts.npm.sec(2:end), diff(ts.npm.sec).*1000) %npm
plot(ts.behavior_cam.sec(2:end), diff(ts.behavior_cam.sec).*1000) %behavior cam
xlim([0 ts.npm.sec(end)])

%aesthetics
title('Dropped Frame Check')
ylabel('inter-frame-interval (ms)')
xlabel('time (s)')
legend('npm','behv','location','best')
setAx(gca);
end

function out_stats(fs, n, l)
disp(strjoin(['~' string(datetime) '~ ////video stats////'], ''))
disp('                                             ')
disp([' *** NPM framerate: ' num2str(fs.npm) ' Hz'])
disp([' *** NPM frames: ' num2str(n.npm)])
disp([' *** NPM length: ' num2str(l.npm.sec/60) ' min'])
disp([' *** Flir behavior framerate: ' num2str(fs.behavior_cam) ' Hz'])
disp([' *** Flir behavior cam frames: ' num2str(n.behavior_cam)])
disp([' *** Flir behavior cam length: ' num2str(l.behavior_cam.sec/60) ' min'])
disp('                                             ')
end