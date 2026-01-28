%% select random videos for behavior scoring to train supervised behavior detection
%created by eartha mae guthman, phd 2026

%% load data
%load header table
analysis_headers = readtable('/Volumes/Kavve/analysis_headers_fully_processed.csv');

%inits
this_path = "/Volumes/falkner/Mae/multisite_hormone_experiments/data/intact";
rng(247);
teg = 6:8; %columns in header for days
shtaygn = {"home", "partner"};
partners = {"balbc", "cdm", "cdf"};
no_mice = size(analysis_header,1);
mahz = randperm(no_mice);
disp(['~* N mice (total): ' num2str(no_mice) ' *~'])

%% build arrays
fldrPath = cell(no_mice, 1); territoryLbl = cell(no_mice, 1); partnerLbl = cell(no_mice, 1);
for ii = 1:no_mice
    %inits
    moyz = string(analysis_headers{mahz(ii), 1});
    tog = num2str(analysis_headers{mahz(ii), teg(randperm(length(teg), 1))});
    folder = string(analysis_headers{mahz(ii), 5});
    shtayg = shtaygn(randperm(length(shtaygn), 1));   
    partner = partners(randperm(length(partners), 1));

    %strings for csv
    fldrPath{ii} = categorical(strjoin([this_path "/" folder "/" tog], ''));
    territoryLbl{ii} = categorical(shtayg{1});
    partnerLbl{ii} = categorical(partner{1});
end

%build table to export
outputTbl = table(fldrPath, territoryLbl, partnerLbl, 'VariableNames', {'Folder', 'Territory', 'Partner'});

%% export to csv
cd(this_path)
writetable(outputTbl, "videos_to_train_supervised_behavior_detection_model.csv")