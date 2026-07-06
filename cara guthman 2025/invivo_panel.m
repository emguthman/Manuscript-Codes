% cara, guthman et al 2025, in vivo figure
close all
clearvars

%% inits, exogenous
e2_mice = {"f3965" "f3972" "f3965"};
e2_dates = {"20230319" "20230327" "20230524"};
oil_mice = {"f3965" "f3972"};
oil_dates = {"20230319" "20230326"};
exogenous_timestamps = [0 2 4 6 16];

%% load exogeous data
cd("/Volumes/AlQahwa/SeeER_paper/exogenous_data");
[e2_raw_green, e2_raw_red, T_e2, T_inj_e2, ZG_e2] = load_exogenous_data(e2_mice, e2_dates, "e2");
[oil_raw_green, oil_raw_red, T_oil, T_inj_oil, ZG_oil] = load_exogenous_data(oil_mice, oil_dates, "oil");
[one_microgram_green, one_microgram_red, T_one, T_inj_one, ZG_one] = load_exogenous_data("m006", "20241214", "e2low");

%% exogenous data panels
%raw plot: f3965, experiment1
plotRawTrace(e2_raw_green{1, 2}, e2_raw_red{1, 2}, "raw_injection_panel_vmh0.pdf", [80 240], T_inj_e2{1})
plotRawTrace(e2_raw_green{2, 2}, e2_raw_red{2, 2}, "raw_injection_panel_vmh1.pdf", [80 240], T_inj_e2{2})
plotRawTrace(e2_raw_green{3, 2}, e2_raw_red{3, 2}, "raw_injection_panel_vmh2.pdf", [80 240], T_inj_e2{3})
%%
plotRawTrace(one_microgram_green{1, 2}, one_microgram_red{1, 2}, "1ug_e2_injection_panel_vmh0.pdf", [75 360], T_inj_one{1})
%%
plotRawTrace(oil_raw_green{1, 2}, oil_raw_red{1, 2}, "oil_injection_panel_vmh0.pdf", [80 240], T_inj_oil{1})
plotRawTrace(oil_raw_green{2, 2}, oil_raw_red{2, 2}, "oil_injection_panel_vmh1.pdf", [80 240], T_inj_oil{2})
plotRawTrace(e2_raw_green{1, 1}, e2_raw_red{1, 1}, "raw_injection_panel_mpo0.pdf", [90 160], T_inj_e2{1})
plotRawTrace(e2_raw_green{2, 1}, e2_raw_red{2, 1}, "raw_injection_panel_mpo1.pdf", [90 160], T_inj_e2{1})
plotRawTrace(e2_raw_green{3, 1}, e2_raw_red{3, 1}, "raw_injection_panel_mpo2.pdf", [90 160], T_inj_e2{1})
%%
plotRawTrace(one_microgram_green{1, 1}, one_microgram_red{1, 1}, "1ug_e2_injection_panel_mpo0.pdf", [90 160], T_inj_one{1})
%%
plotRawTrace(oil_raw_green{1, 1}, oil_raw_red{1, 1}, "oil_injection_panel_mpo0.pdf", [90 160], T_inj_oil{1})
plotRawTrace(oil_raw_green{2, 1}, oil_raw_red{2, 1}, "oil_injection_panel_mpo1.pdf", [90 160], T_inj_oil{2})

%plot population onset response
onset_fig = figure("Position", [400 400 800 300], "Color", [1 1 1]);
plotExogenousOnset(e2_raw_green, T_inj_e2{1}, 1/12);
exportgraphics(onset_fig, "exogenous_onset.pdf", "ContentType", "vector")
close all

%scatter plot
[exogenous_lme, exogenous_alpha] = plotExogenousScatter(e2_raw_green, e2_raw_red, oil_raw_green, oil_raw_red, exogenous_timestamps, T_inj_e2, 1);
close all

%plot mean response to e2 in mpoa
mpoa_fig = figure("Position", [400 400 300 300], "Color", [1 1 1]);
pop_idx = 1;
plotMeanExogenous(e2_raw_green, pop_idx, T_inj_e2{1}, 1/12);
exportgraphics(mpoa_fig, "exogenous_mpoa.pdf", "ContentType", "vector")
close all

%plot mean response to e2 in vmhvl
vmhvl_fig = figure("Position", [400 400 300 300], "Color", [1 1 1]);
pop_idx = 2;
plotMeanExogenous(e2_raw_green, pop_idx, T_inj_e2{1}, 1/12);
exportgraphics(vmhvl_fig, "exogenous_vmhvl.pdf", "ContentType", "vector")
close all

%plot mean response to oil in mpoa, vmhvl
oil_fig = figure("Position", [400 400 300 300], "Color", [1 1 1]);
pop_idx = [1 2];
plotMeanExogenous(oil_raw_green, pop_idx, T_inj_oil{1}, 1/12);
exportgraphics(oil_fig, "exogenous_oil.pdf", "ContentType", "vector")
close all

%% inits, endogenous
%basic headers
estrous_mice = {"f3965", "m003", "m003", "m004", "m005", "m006"}; %, "m007"};
estrous_dates = {"20230129", "20240404", "20240709", "20240807", "20241024", "20241025"}; %, , "20250114"};
M = length(estrous_mice);
estrous_alignment = [141, 441, -10, 567, 2842, 3457]; %, 490]; %490 or 778
down_sample_flag = [0 1 0 0 0 0 0]; %0, first m003 recording taken at 2x speed (1/2x frame rate, 0.00175 Hz instead of 0.0035 Hz)

%% load endogeous data
%inits
cd("/Volumes/AlQahwa/SeeER_paper/endogenous_data");
T_estrous = (1/12:1/12:24*4);
redblue_map = [linspace(120, 255, 128)', linspace(160, 255, 128)', linspace(193, 255, 128)'; ...
    255, 255, 255; linspace(255, 224, 128)', linspace(255, 97, 128)', linspace(255, 97, 128)'] ./ 255;

%load
[b2b_estrous_raw_green, b2b_estrous_raw_red] = ...
    load_endogenous_data("m005", "20241024", 2842, 0, "estrous");

%% plot back to back
plotBackToBack(b2b_estrous_raw_green, b2b_estrous_raw_red)

%load
[estrous_raw_green, estrous_raw_red, estrous_tmac, estrous_motion] = ...
    load_endogenous_data(estrous_mice, estrous_dates, estrous_alignment, down_sample_flag, "estrous");
[peak_raw_green, peak_raw_red, peak_tmac, peak_motion] = ...
    load_endogenous_data(estrous_mice, estrous_dates, estrous_alignment, down_sample_flag, "peak");

%get dff based on metestrus/first 12hr baseline
[estrous_dff, estrous_z] = compute_estrous_dff(estrous_raw_green);
[estrous_red_dff, estrous_red_z] = compute_estrous_dff(estrous_raw_red);
[peak_dff, peak_z] = compute_estrous_dff(peak_raw_green);

%fill nans
% estrous_tmac{1}(m, isnan(estrous_tmac{1}(m, :))) = eps;
% estrous_tmac{2}(m, isnan(estrous_tmac{2}(m, :))) = eps;

%get max of findpeaks results
estrous_idx{1} = zeros(height(estrous_dff{1}), 1); estrous_idx{2} = zeros(height(estrous_dff{2}), 1);
peak_idx{1} = zeros(height(peak_dff{1}), 1); peak_idx{2} = zeros(height(peak_dff{2}), 1);
peak_val{1} = zeros(height(peak_dff{1}), 1); peak_val{2} = zeros(height(peak_dff{2}), 1);
for m = 1:M
    [~, estrous_idx{1}(m)] = max(estrous_z{1}(m, :));
    [~, estrous_idx{2}(m)] = max(estrous_z{2}(m, :));
    peak_val{1}(m) = max(peak_z{1}(m, :));
    peak_val{2}(m) = max(peak_z{2}(m, :));
end
[~, peak_idx{1}] = sort(peak_val{1}, "descend");
[~, peak_idx{2}] = sort(peak_val{2}, "descend");

%% compute peak value within each window and save in table
[estrousTbl, estrous_rm_anovas{1}, estrous_rm_posthoc{1}, estrous_rm_anovas{2}, estrous_rm_posthoc{2}]  = getEstrousSignal(estrous_z, "estrous_SeeER.pdf");
[peakTbl, peak_rm_anovas{1}, peak_rm_posthoc{1}, peak_rm_anovas{2}, peak_rm_posthoc{2}]  = getEstrousSignal(peak_z, "peak_SeeER.pdf");

%% endogenous panels
%raw plot: f3965, experiment1
plotRawEstrous(estrous_dff{1}(1, :), estrous_red_dff{1}(1, :), "raw_estrous_panel_f3965_mpo.pdf", [50 170])
plotRawEstrous(estrous_dff{2}(1, :), estrous_red_dff{2}(1, :), "raw_estrous_panel_f3965_vmh.pdf", [50 170])

% peak time correlations
[pktime_rho, pktime_corr_pval] = compute_peakTimeCorrelation(estrous_idx);
[~, p_estrous_latency, latency_stat, latency_test] = compute_peakTimeDifference(estrous_idx);

%% vmh estrous aligned heatmap
% get data to plot
mpoa_plot = estrous_z{1};
mpoa_plot(isnan(mpoa_plot)) = eps;
vmh_plot = estrous_z{2};
vmh_plot(isnan(vmh_plot)) = eps;

% plot heatmap
estrousFig = figure("Position", [400, 1200, 660, 380], "Color", [1 1 1]);
colormap(redblue_map)
[~, sortIdx] = sort(estrous_idx{2});

%mpoa
subplot(2, 1, 1) %mpoa
imagesc(mpoa_plot(sortIdx, :));
clim([-2 2])
mpoAx = gca;
setAx(mpoAx);
set(mpoAx, "XTick", 12*12:12*12:width(estrous_z{1}))
set(mpoAx, "XTickLabel", 12:12:width(estrous_z{1})/12)
colorbar
title("MPO + MPN")
xlabel("time (hr)")
ylabel("mice")

%vmhvl
subplot(2, 1, 2) %vmh
imagesc(vmh_plot(sortIdx, :));
clim([-2 2])
vmhAx = gca;
setAx(vmhAx);
set(vmhAx, "XTick", 12*12:12*12:width(estrous_z{2}))
set(vmhAx, "XTickLabel", 12:12:width(estrous_z{2})/12)
colorbar
title("VMHvl")
xlabel("time (hr)")
ylabel("mice")

%save and close
exportgraphics(estrousFig, "estrous_heatmap.pdf", "ContentType", "vector")

%% vmh peak aligned heatmap
% get data to plot
mpoa_plot = peak_z{1};
mpoa_plot(isnan(mpoa_plot)) = eps;
vmh_plot = peak_z{2};
vmh_plot(isnan(vmh_plot)) = eps;

% plot heatmap
peakFig = figure("Position", [400, 1200, 660, 380], "Color", [1 1 1]);
colormap(redblue_map)

%mpoa
subplot(2, 1, 1) %mpoa
imagesc(mpoa_plot(sortIdx, :));
clim([-2 2])
mpoAx = gca;
setAx(mpoAx);
set(mpoAx, "XTick", 12*12:12*12:width(peak_z{1}))
set(mpoAx, "XTickLabel", 12:12:width(peak_z{1})/12)
colorbar
title("MPO + MPN")
xlabel("time (hr)")
ylabel("mice")

%vmhvl
subplot(2, 1, 2) %vmh
imagesc(vmh_plot(sortIdx, :));
clim([-2 2])
vmhAx = gca;
setAx(vmhAx);
set(vmhAx, "XTick", 12*12:12*12:width(peak_z{2}))
set(vmhAx, "XTickLabel", 12:12:width(peak_z{2})/12)
colorbar
title("VMHvl")
xlabel("time (hr)")
ylabel("mice")

%save and close
exportgraphics(peakFig, "peak_heatmap.pdf", "ContentType", "vector")
close all

%% correlate experimentally aligned vs mismatched MPO+MPN and VMHvl signals
[estrous_ks_stat, estrous_ks_p, estrous_correlations] = experimental_correlations(estrous_z{1}, estrous_z{2}, "estrous_correlations.pdf");
[peak_ks_stat, peak_ks_p, peak_correlations] = experimental_correlations(peak_z{1}, peak_z{2}, "peak_correlations.pdf");

%% subroutines
function [h, p, derStat, testID] = compute_peakTimeDifference(peak_times)
%inits
fs = 1/12;

%get arrays
mpo_latency = peak_times{1} * fs;
vmh_latency = peak_times{2} * fs;

%normality check
mpo_h = swtest(mpo_latency);
vmh_h = swtest(vmh_latency);

%compute
if ~mpo_h && ~vmh_h %if both normal (h = 0)
    [h, p, ~, derStat] = ttest(mpo_latency, vmh_latency);
    testID = "ttest";
else %non-parametric
    [h, p, derStat] = signrank(mpo_latency, vmh_latency);
    testID = "wilcoxon";
end


%plot
%inits
latencyFig = figure("Position",[400 400 160 320],"Color",[1 1 1]);
hold on

%plot
plot([mpo_latency vmh_latency]', "Color", [0 0 0])
% shadedErrorBar([1 2], mean([mpo_latency vmh_latency]), sem([mpo_latency vmh_latency], 1))

%aesthetics
latAx = gca;
setAx(latAx);
xlim([0.67 2.33])
ylim([24 84])
set(latAx, "YTick", (0:12:24*4))
set(latAx, "XTick", (1:2), "XTickLabel", {'MPO', 'VMH'})
xlabel("Population")
ylabel("Peak Latency (hr)")
title("Peaktime Latency")

%save and close
exportgraphics(latencyFig, "estrous_peak_latency.pdf", "ContentType", "vector")
close all
end

function [rho, p] = compute_peakTimeCorrelation(peak_times)
%inits
fs = 1/12;

%get arrays
mpo_latency = peak_times{1};
vmh_latency = peak_times{2};
xvals = 0:4*42/fs;

%compute regression
[beta, ~, ~, ~, ~] = regress(vmh_latency, [ones(size(mpo_latency)) mpo_latency]);
regression_line = polyval([beta(2) beta(1)], xvals);
[rho, p] = corr([mpo_latency vmh_latency]);

%plot
%inits
corrFig = figure("Position",[400 400 320 320],"Color",[1 1 1]);
hold on

%plot
plot(xvals,xvals, "LineStyle", "--", "LineWidth", 0.5, "Color", [0 0 0])
scatter(mpo_latency, vmh_latency, 75, [0 0 0], 'filled')
plot(xvals, regression_line, "LineWidth", 1, "Color", [0 0 0])

%aesthetics
corrAx = gca;
setAx(gca);
xlim([0 (4*24)/fs])
ylim([0 (4*24)/fs])
set(corrAx, "XTick", (0:24:24*4)./fs, "YTick", (0:24:24*4)./fs)
set(corrAx, "XTickLabels", 0:24:24*4, "YTickLabels", 0:24:24*4)
xlabel("MPO/MPN Peak Latency (hr)")
ylabel("VMHvl Peak Latency (hr)")
title("Peaktime Latency")

%save and close
exportgraphics(corrFig, "peaktime_correlation.pdf", "ContentType", "vector")
close all
end

function plotMotionCorrection(groyn, royz, tmac, motion)
%inits
motionFig = figure("Position", [400 400 640 640], "Color", [1 1 1]);
green_clr = [62 162 44] ./ 255;

%convert to dff, e makes me think tmac is dff
groyn_dff = 100 .* (groyn ./ mean(groyn));
royz_dff = 100 .* (royz ./ mean(royz));
tmac_dff = 100 .* tmac;
motion_dff = 100 .* motion;

%raw green
subplot(2, 2, 1)
plot(groyn_dff, "Color", green_clr)
groynAx = gca;
setAx(groynAx);
xlim([0 288])
set(groynAx, "XTick", 0:72:length(groyn_dff))

%raw red
subplot(2, 2, 3)
plot(royz_dff, "Color", [0 0 0])
royzAx = gca;
setAx(royzAx);
ylim([40 140])
xlim([0 288])
set(royzAx, "XTick", 0:72:length(groyn_dff));

%tmac
subplot(2, 2, 2)
plot(tmac_dff, "Color", green_clr)
tmacAx = gca;
setAx(tmacAx);
ylim([-10 90])
xlim([0 288])
set(tmacAx, "XTick", 0:72:length(groyn_dff));

%motion
subplot(2, 2, 4)
plot(motion_dff, "Color", [0 0 0])
motionAx = gca;
setAx(motionAx);
ylim([-40 60])
xlim([0 288])
set(motionAx, "XTick", 0:72:length(groyn_dff));

%save
exportgraphics(motionFig, "motion_correction_panel.pdf", "ContentType", "vector")
end

function plotRawTrace(groyn, royz, filename, y_limits, t_inj)
%inits
fs = 1/12;
t_inj = floor(fs^(-1) * t_inj);
rawDataFig = figure("Position", [400 400 640 640], "Color", [1 1 1]);
green_clr = [62 162 44] ./ 255;

%convert to dff
groyn_z = 100 + 100 .* ((groyn - mean(groyn(t_inj-fs^(-1)+1:t_inj), "omitmissing")) ./ mean(groyn(t_inj-fs^(-1)+1:t_inj), "omitmissing"));
royz_z = 100 + 100 .* ((royz - mean(royz(t_inj-fs^(-1)+1:t_inj), "omitmissing")) ./ mean(royz(t_inj-fs^(-1)+1:t_inj), "omitmissing"));

%raw green
subplot(2, 1, 1)
plot(groyn_z, "Color", green_clr)
groynAx = gca;
setAx(groynAx);
ylim(y_limits)
xlim([13.5*12 (24+13.5)*12]) %1ug experiment
%xlim([0 24*12])
% set(groynAx, "XTick", 0:48:lengt(groyn_z))
set(groynAx, "XTick", t_inj:10*12:length(groyn_z)) %1ug experiment
set(groynAx, "XTickLabels", 0:10:length(groyn_z)/12)

%raw red
subplot(2, 1, 2)
plot(royz_z, "Color", [0 0 0])
royzAx = gca;
setAx(royzAx);
ylim(y_limits)
xlim([13.5*12 (24+13.5)*12]) %1ug experiment
%xlim([0 24*12])
% set(royzAx, "XTick", 0:48:length(groyn_z));
set(royzAx, "XTick", t_inj:10*12:length(groyn_z)) %1ug experiment
set(royzAx, "XTickLabels", 0:10:length(groyn_z)/12)

%save
exportgraphics(rawDataFig, filename, "ContentType", "vector")
end

function plotRawEstrous(groyn, royz, filename, y_limits)
%inits
fs = 1/12;
rawDataFig = figure("Position", [400 400 640 640], "Color", [1 1 1]);
green_clr = [62 162 44] ./ 255;

%raw green
subplot(2, 1, 1)
plot(groyn, "Color", green_clr)
groynAx = gca;
setAx(groynAx);
ylim(y_limits)
set(groynAx, "XTick", 0:(1/fs)^(2):length(groyn));
set(groynAx, "XTickLabels", 0:fs^-1:length(groyn)/(fs^-1))

%raw red
subplot(2, 1, 2)
plot(royz, "Color", [0 0 0])
royzAx = gca;
setAx(royzAx);
ylim(y_limits)
set(royzAx, "XTick", 0:12^2:length(royz));
set(royzAx, "XTickLabels", 0:12:length(royz)/12)

%save
exportgraphics(rawDataFig, filename, "ContentType", "vector")
end

function [lme, alpha] = plotExogenousScatter(e2_green, e2_red, oil_green, oil_red, ts, t_e2, idx)
%inits
scatterPlot = figure("Position", [400 600 1200 420], "Color", [1 1 1]);
green_clr = [62 162 44] ./ 255;
ts = ts + t_e2{idx};
time_windows = {"0h", "2h", "4h", "6h", "16h"};
fs = 1/12;
xvals = 0:300;
anovaMat = cell(2, 1);
anova_p = cell(2, 1);
anova_tbl = cell(2, 1);
anova_stats = cell(2, 1);
posthoc = cell(2, 1);
e2_green_mat = cell(2, length(ts));
e2_green_baseline = cell(2, 3);
oil_green_mat = cell(2, length(ts));
oil_green_baseline = cell(2, 2);
e2_red_mat = cell(2, length(ts));
e2_red_baseline = cell(2, 3);
oil_red_mat = cell(2, length(ts));
oil_red_baseline = cell(2, 2);
pops = {"MPO" "VMH"};
lme = cell(size(pops));
alpha = cell(size(pops));
exps = [1 1 2]; mice = [1 2 1];
headers = {'Mouse', 'Experiment', 'Treatment', 'TimeWindow', 'Population', 'DFF'};
vartypes = {'categorical', 'categorical', 'categorical', 'categorical', 'categorical', 'double'};
dataTbl = table('Size', [0 numel(headers)], 'VariableTypes', vartypes, 'VariableNames', headers);

%collapse all sessions into one matrix that's normalized to the hour before injection
for pop = [1 2] %mpoa, vmh
    for loop = 1:length(ts)
        %inits
        ts_stop = ts(loop) * fs^(-1);

        %build matrices
        for versye = 1:size(e2_green, 1)
            e2_green_mat{pop, loop} = cat(1, e2_green_mat{pop, loop}, e2_green{versye, pop}(ts_stop-11:ts_stop)');
            if loop == 1
                e2_green_baseline{pop, versye} = mean(e2_green{versye, pop}(ts_stop-11:ts_stop), "omitmissing");
            end
            e2_green_mat{pop, loop}(versye, :) = 100 + 100 .* ((e2_green_mat{pop, loop}(versye, :) - e2_green_baseline{pop, versye}) ...
                ./ e2_green_baseline{pop, versye}); %normalize
            e2_red_mat{pop, loop} = cat(1, e2_red_mat{pop, loop}, e2_red{versye, pop}(ts_stop-11:ts_stop)');
            if loop == 1
                e2_red_baseline{pop, versye} = mean(e2_red{versye, pop}(ts_stop-11:ts_stop), "omitmissing");
            end
            e2_red_mat{pop, loop}(versye, :) = 100 + 100 .* ((e2_red_mat{pop, loop}(versye, :) - e2_red_baseline{pop, versye}) ...
                ./ e2_red_baseline{pop, versye});
        end
        for versye = 1:size(oil_green, 1)
            oil_green_mat{pop, loop} = cat(1, oil_green_mat{pop, loop}, oil_green{versye, pop}(ts_stop-11:ts_stop)');
            if loop == 1
                oil_green_baseline{pop, versye} = mean(oil_green{versye, pop}(ts_stop-11:ts_stop), "omitmissing");
            end
            oil_green_mat{pop, loop}(versye, :) = 100 + 100 .* ((oil_green_mat{pop, loop}(versye, :) - oil_green_baseline{pop, versye}) ...
                ./ oil_green_baseline{pop, versye}); %normalize
            oil_red_mat{pop, loop} = cat(1, oil_red_mat{pop, loop}, oil_red{versye, pop}(ts_stop-11:ts_stop)');
            if loop == 1
                oil_red_baseline{pop, versye} = mean(oil_red{versye, pop}(ts_stop-11:ts_stop), "omitmissing");
            end
            oil_red_mat{pop, loop}(versye, :) = 100 + 100 .* ((oil_red_mat{pop, loop}(versye, :) - oil_red_baseline{pop, versye}) ...
                ./ oil_red_baseline{pop, versye});
        end
    end
end

for pop = [1 2] %mpoa, vmh
    %inits
    anovaMat{pop} = zeros(24, length(ts));

    %get oil line
    beta_oil = regress(oil_green_mat{pop, 1}(:), [ones(size(oil_green_mat{pop, 1}(:))) oil_red_mat{pop, 1}(:)]);
    oil_line = beta_oil(1) + beta_oil(2) .* xvals;

    %get e2 line
    beta_e2 = regress(e2_green_mat{pop, 1}(:), [ones(size(e2_green_mat{pop, 1}(:))) e2_red_mat{pop, 1}(:)]);
    e2_line = beta_e2(1) + beta_e2(2) .* xvals;

    for loop = 1:length(ts)
        %inits
        subplot(2, length(ts), length(ts)*(pop-1) + loop)

        %plot
        hold on
        plot(xvals, oil_line, "LineStyle", "--", "Color", [0 0 0])
        plot(xvals, e2_line, "LineStyle", "--", "Color", green_clr)
        for o = 1:size(oil_red_mat{pop, loop}, 1)
            scatter(oil_red_mat{pop, loop}(o, :), oil_green_mat{pop, loop}(o, :), 75, [0.5 0.5 0.5] ./ o, "filled")
        end
        for g = 1:size(e2_red_mat{pop, loop}, 1)
            scatter(e2_red_mat{pop, loop}(g, :), e2_green_mat{pop, loop}(g, :), 75, green_clr ./ g, "filled")
        end

        %aesthetics
        scatterAx = gca;
        setAx(scatterAx);
        xlim([90 115])
        if pop == 1
            ylim([90 160])
        elseif pop == 2
            ylim([80 250])
        end
        if loop == 1
            ylabel("green fluor. (a.u.)")
            if pop == 1
                legend({'Oil baseline', 'E2 baseline', 'Oil', 'E2'})
            end
        end
        if pop == 2
            xlabel("red fluor. (a.u.)")
        end

        %build table
        for versye = 1:size(e2_green_mat{pop, loop}, 1)
            for data_point = 1:length(e2_green_mat{pop, loop}(versye, :))
                this_row = {categorical(mice(versye)), categorical(exps(versye)), "E2", time_windows{loop}, ...
                    pops{pop}, e2_green_mat{pop, loop}(versye, data_point)};
                dataTbl = cat(1, dataTbl, this_row);
            end
        end
        for versye = 1:size(oil_green_mat{pop, loop}, 1)
            for data_point = 1:length(oil_green_mat{pop, loop}(versye, :))
                this_row = {categorical(mice(versye)), categorical(exps(versye)), "Oil", time_windows{loop}, ...
                    pops{pop}, oil_green_mat{pop, loop}(versye, data_point)};
                dataTbl = cat(1, dataTbl, this_row);
            end
        end
    end
end

%run stats, linear mixed effects model separately by population
formula = 'DFF ~ Treatment * TimeWindow + (Experiment | Mouse)';

%MPO
mpoTbl = dataTbl(dataTbl.Population == "MPO", :);
[lme{1}, alpha{1}] = runLME(mpoTbl, formula);

%VMH
vmhTbl = dataTbl(dataTbl.Population == "VMH", :);
[lme{2}, alpha{2}] = runLME(vmhTbl, formula);

%save figure
exportgraphics(scatterPlot, "exogenous_results.pdf", "ContentType", "vector")
end

function [lme, lme_alpha] = runLME(tabelle, formula)
% treatments = string(unique(tabelle.Treatment));
windows = string(unique(tabelle.TimeWindow));
lme = cell(size(windows));
lme_pvals = [];

for bb = 1:length(windows)
    if bb == 1
        window_order = windows([1 2 3 4 5]);
    elseif bb == 2
        window_order = windows([2 3 4 5 1]);
    elseif bb == 3
        window_order = windows([3 4 5 1 2]);
    elseif bb == 4
        window_order = windows([4 5 1 2 3]);
    elseif bb == 5
        window_order = windows([5 1 2 3 4]);
    end

    %reorder categories based on desired order
    tabelle.TimeWindow = reordercats(tabelle.TimeWindow, window_order);

    %sort tabel based on column
    tabelle = sortrows(tabelle, 'TimeWindow');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% fit linear mixed effects model %%%
    lme{bb} = fitlme(tabelle, formula, 'FitMethod', 'REML');

    %collate  LME p values
    lme_pvals = cat(1, lme_pvals, lme{bb}.Coefficients.pValue(2));
end
lme_alpha = drsFDRpval(lme_pvals);
end

function plotMeanExogenous(datn, pop, t_inj, fs)
%inits
t_inj = fs^(-1) * t_inj;
if size(datn{1, 1}, 1) >= 576
    end_idx = 48 * fs^(-1); %24hr in
else
    end_idx = 24 * fs^(-1); %24hr in
end
M = size(datn, 1); P = length(pop);
datnMat = [];

%get matrix of data
for p = 1:P
    this_pop = pop(p);
    for m = 1:M
        %get data
        try
            dizDatn = datn{m, this_pop}(1:end_idx)';
        catch
            dizDatn = datn{m, this_pop}';
            dizDatn(length(dizDatn)+1:end_idx) = NaN;
        end
        u_baseline = mean(dizDatn(t_inj-fs^(-1):t_inj), "omitmissing");
        dffDatn = 100 + 100 .* ((dizDatn - u_baseline)/u_baseline);
        datnMat = cat(1, datnMat, dffDatn);
    end
end

%plot shaded error bar
hold on
line([-4 48], [100 100], 'linestyle', '--', 'linewidth', 1, 'color', [0 0 0])
shadedErrorBar((-4+fs):fs:(size(datnMat, 2)*fs - 4), mean(datnMat, 1, 'omitnan'), sem(datnMat, 1));

%aesthetics
setAx(gca);
ylim([90 210])
if size(datn{1, 1}, 1) >= 576
    xlim([-4 44])
else
    xlim([-4 20])
end
ylabel("SeeER (\DeltaF/F)")
xlabel("time (hr, rel to injection)")
end

function plotExogenousOnset(datn, t_inj, fs)
%inits
t_inj = fs^(-1) * t_inj;
if size(datn{1, 1}, 1) >= 576
    end_idx = 48 * fs^(-1); %24hr in
else
    end_idx = 24 * fs^(-1); %24hr in
end
[M, P] = size(datn);
mapping{1} = [255 184 28] ./ 255;
mapping{2} = [62 162 44] ./ 255;
datnMat = [];
for m = 1:M+1
    subplot(1, M+1, m)
    hold on
    line([-4 48], [100 100], 'linestyle', '--', 'linewidth', 1, 'color', [0 0 0])
end

%get matrix of data
for m = 1:M
    for p = 1:P
        %get data
        try
            dizDatn = datn{m, p}(1:end_idx)';
        catch
            dizDatn = datn{m, p}';
            dizDatn(length(dizDatn)+1:end_idx) = NaN;
        end
        %dizDatn = 100 .* dizDatn + 100 - mean(100 .* dizDatn, "omitmissing"); %convert to dff centered around 100
        dffDatn = 100 + 100 .* ((dizDatn - mean(dizDatn(t_inj-fs^(-1)+1:t_inj), "omitmissing")) ./ mean(dizDatn(t_inj-fs^(-1)+1:t_inj), "omitmissing")); %convert to dff centered around 100
        % u_baseline = mean(dizDatn(1:t_inj), "omitmissing");
        % baseline_shift = 100 - u_baseline;

        %dff relative to baseline
        % dffDatn = dizDatn + baseline_shift;

        %plot data, all sessions
        subplot(1, M+1, 1)
        plot((-4+fs):fs:(length(dffDatn)*fs - 4), dffDatn, 'linewidth', 1, 'color', mapping{p})

        %plot data, this recording only
        subplot(1, M+1, 1+m)
        plot((-4+fs):fs:(length(dffDatn)*fs - 4), dffDatn, 'linewidth', 1, 'color', mapping{p})

    end
end

%aesthetics
for m = 1:M+1
    %inits
    subplot(1, M+1, m)
    setAx(gca);

    %data range
    ylim([90 130])
    xlim([-2 4])

    %labels
    ylabel("SeeER (\DeltaF/F)")
    xlabel("time (hr, rel to injection)")
    if m == 1
        title("All Sessions")
    else
        title(strjoin(["Session #" num2str(m-1)], ''))
    end
end
end