function [raw_green, raw_red, tmac, motion] = load_endogenous_data(mahz, dates, alignment, flag, peak_align)
%inits
M = length(mahz);
raw_green{1} = NaN(M, 24 * 12 * 4);
raw_green{2} = NaN(M, 24 * 12 * 4);
raw_red{1} = NaN(M, 24 * 12 * 4);
raw_red{2} = NaN(M, 24 * 12 * 4);
tmac{1} = NaN(M, 24 * 12 * 4);
tmac{2} = NaN(M, 24 * 12 * 4);
motion{1} = NaN(M, 24 * 12 * 4);
motion{2} = NaN(M, 24 * 12 * 4);
warning('off')

for m = 1:M
    %inits
    mouse = mahz{m};
    tog = dates{m};
    fileheader = strjoin(["estrous_" mouse "_" tog "_"], '');

    %load
    raw_mpoa_green = readtable(strjoin([fileheader "mpoa_raw_green.csv"], ''));
    raw_vmh_green = readtable(strjoin([fileheader "vmh_raw_green.csv"], ''));
    raw_mpoa_red = readtable(strjoin([fileheader "mpoa_raw_red.csv"], ''));
    raw_vmh_red = readtable(strjoin([fileheader "vmh_raw_red.csv"], ''));
    mpoa_tmac = readtable(strjoin([fileheader "mpoa_tmac.csv"], ''));
    vmh_tmac = readtable(strjoin([fileheader "vmh_tmac.csv"], ''));
    mpoa_motion = readtable(strjoin([fileheader "mpoa_motion.csv"], ''));
    vmh_motion = readtable(strjoin([fileheader "vmh_motion.csv"], ''));

    %put in outputs
    tmp_mpoa_green = raw_mpoa_green.RawGreenData; %mpoa
    tmp_vmh_green = raw_vmh_green.RawGreenData; %vmh
    tmp_mpoa_red = raw_mpoa_red.RawRedData;
    tmp_vmh_red = raw_vmh_red.RawRedData;
    tmp_mpoa_tmac= mpoa_tmac.TMACData; %mpoa
    tmp_vmh_tmac = vmh_tmac.TMACData; %vmh
    tmp_mpoa_motion = mpoa_motion.MotionData;
    tmp_vmh_motion = vmh_motion.MotionData;
    
    if flag(m)
        %downsample
        tmp_mpoa_green(2:2:length(tmp_mpoa_green)) = [];
        tmp_vmh_green(2:2:length(tmp_vmh_green)) = [];
        tmp_mpoa_red(2:2:length(tmp_mpoa_red)) = [];
        tmp_vmh_red(2:2:length(tmp_vmh_red)) = [];
        tmp_mpoa_tmac(2:2:length(tmp_mpoa_green)) = [];
        tmp_vmh_tmac(2:2:length(tmp_vmh_green)) = [];
        tmp_mpoa_motion(2:2:length(tmp_mpoa_red)) = [];
        tmp_vmh_motion(2:2:length(tmp_vmh_red)) = [];
    end

    %get subsection, inits
    if strcmp(peak_align, "estrous")
        % t_window(1) = alignment(m);
        t_window(1) = alignment(m) - 12*24*4; %back to back
    elseif strcmp(peak_align, "peak")
        try
            [~, t_window(1)] = max(tmp_vmh_green(alignment(m):end));
            t_window(1) = t_window(1) - 12*24*2 + alignment(m) - 1;
        catch
            [~, t_window(1)] = max(tmp_vmh_green);
            t_window(1) = t_window(1) - 12*24*2 - 1;
        end
    end
    t_shift = 0;
    t_window(2) = t_window(1) + 12*24*4*2 - 1;
    if t_window(1) <= 0 
        tmp_alignment = t_window(1);
        t_window(1) = 1;
        t_window(2) = 12*24*4 + tmp_alignment;
        t_shift =  -1*tmp_alignment;
    end
    t_stop = 12 * 24 * 4 * 2;
    if t_window(2) > length(tmp_mpoa_green)
        t_stop = t_stop - (t_window(2) - length(tmp_mpoa_green));
        t_window(2) = length(tmp_mpoa_green);
    end

    %fill arrays
    % raw_green{1}(m, (1 + t_shift):(t_stop)) = tmp_mpoa_green(t_window(1):t_window(2)); %mpoa
    % raw_green{2}(m, (1 + t_shift):(t_stop)) = tmp_vmh_green(t_window(1):t_window(2)); %vmh
    % raw_red{1}(m, (1 + t_shift):(t_stop)) = tmp_mpoa_red(t_window(1):t_window(2)); %mpoa
    % raw_red{2}(m, (1 + t_shift):(t_stop)) = tmp_vmh_red(t_window(1):t_window(2)); %vmh
    % tmac{1}(m, (1 + t_shift):(t_stop)) = tmp_mpoa_tmac(t_window(1):t_window(2)); %mpoa
    % tmac{2}(m, (1 + t_shift):(t_stop)) = tmp_vmh_tmac(t_window(1):t_window(2)); %vmh
    % motion{1}(m, (1 + t_shift):(t_stop)) = tmp_mpoa_motion(t_window(1):t_window(2)); %mpoa
    % motion{2}(m, (1 + t_shift):(t_stop)) = tmp_vmh_motion(t_window(1):t_window(2)); %vmh

    %back to back
    raw_green{1}(m, (1 + t_shift):(t_stop)) = tmp_mpoa_green(t_window(1):t_window(2)); %mpoa
    raw_green{2}(m, (1 + t_shift):(t_stop)) = tmp_vmh_green(t_window(1):t_window(2)); %vmh
    raw_red{1}(m, (1 + t_shift):(t_stop)) = tmp_mpoa_red(t_window(1):t_window(2)); %mpoa
    raw_red{2}(m, (1 + t_shift):(t_stop)) = tmp_vmh_red(t_window(1):t_window(2)); %vmh
    tmac{1}(m, (1 + t_shift):(t_stop)) = tmp_mpoa_tmac(t_window(1):t_window(2)); %mpoa
    tmac{2}(m, (1 + t_shift):(t_stop)) = tmp_vmh_tmac(t_window(1):t_window(2)); %vmh
    motion{1}(m, (1 + t_shift):(t_stop)) = tmp_mpoa_motion(t_window(1):t_window(2)); %mpoa
    motion{2}(m, (1 + t_shift):(t_stop)) = tmp_vmh_motion(t_window(1):t_window(2)); %vmh
end
end