function [raw_green, raw_red, T, T_inj, ZG] = load_exogenous_data(mahz, dates, inj)
%inits
M = length(mahz);
raw_green = cell(M, 2);
raw_red = cell(M, 2);
tmac = cell(M, 2);
motion = cell(M, 2);
T = cell(M, 1);
T_inj = cell(M, 1);
ZG = cell(M, 1);
warning('off')

for m = 1:M
    %inits
    mouse = mahz{m};
    tog = dates{m};
    fileheader = strjoin([inj "_" mouse "_" tog "_"], '');

    %load
    raw_mpoa_green = readtable(strjoin([fileheader "mpoa_raw_green.csv"], ''));
    raw_vmh_green = readtable(strjoin([fileheader "vmh_raw_green.csv"], ''));
    raw_mpoa_red = readtable(strjoin([fileheader "mpoa_raw_red.csv"], ''));
    raw_vmh_red = readtable(strjoin([fileheader "vmh_raw_red.csv"], ''));
    timestamps = readtable(strjoin([fileheader "timestamps.csv"], ''));
    inj_time = readtable(strjoin([fileheader "injection_time.csv"], ''));
    zg_time = readtable(strjoin([fileheader "zeitgeber_time.csv"], ''));

    %put in outputs
    raw_green{m, 1} = raw_mpoa_green.RawGreenData; %mpoa
    raw_green{m, 2} = raw_vmh_green.RawGreenData; %vmh
    raw_red{m, 1} = raw_mpoa_red.RawRedData;
    raw_red{m, 2} = raw_vmh_red.RawRedData;
    T{m} = timestamps.Time_hr_;
    try
        T_inj{m} = inj_time.zeitgeber_time;
    catch
        T_inj{m} = inj_time.timestamp;
    end
    try
        ZG{m} = zg_time.zeitgeber_time;
    catch
        ZG{m} = zg_time.timestamp;
    end
end
end