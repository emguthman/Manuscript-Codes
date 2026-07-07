function [tbl, mpo_ra, mpo_posthoc, vmh_ra, vmh_posthoc] = getEstrousSignal(signal, filename)
%inits
t_metestrus = 1:24*12; t_diestrus = 24*12+1:24*12*2; t_proestrus = 24*12*2+1:12*12*5; t_estrus = 12*12*5+1:24*12*4;
headers = {'population', 'metestrus', 'diestrus', 'proestrus', 'estrus'};
vartypes = {'categorical', 'double', 'double', 'double', 'double'};
tbl = table('Size', [0 numel(headers)], 'VariableTypes', vartypes, 'VariableNames', headers);

for m = 1:height(signal{1})
    %get data
    this_metestrus{1} = signal{1}(m, t_metestrus); 
    this_metestrus{2} = signal{2}(m, t_metestrus); 
    this_diestrus{1} = signal{1}(m, t_diestrus); 
    this_diestrus{2} = signal{2}(m, t_diestrus); 
    this_proestrus{1} = signal{1}(m, t_proestrus); 
    this_proestrus{2} = signal{2}(m, t_proestrus); 
    this_estrus{1} = signal{1}(m, t_estrus); 
    this_estrus{2} = signal{2}(m, t_estrus); 

    % %pull indices
    % [~, mpo_metestrus_idx] = max(abs(this_metestrus{1}), [], "omitmissing");
    % [~, vmh_metestrus_idx] = max(abs(this_metestrus{2}), [], "omitmissing");
    % [~, mpo_diestrus_idx] = max(abs(this_diestrus{1}), [], "omitmissing");
    % [~, vmh_diestrus_idx] = max(abs(this_diestrus{2}), [], "omitmissing");
    % [~, mpo_proestrus_idx] = max(abs(this_proestrus{1}), [], "omitmissing");
    % [~, vmh_proestrus_idx] = max(abs(this_proestrus{2}), [], "omitmissing");
    % [~, mpo_estrus_idx] = max(abs(this_estrus{1}), [], "omitmissing");
    % [~, vmh_estrus_idx] = max(abs(this_estrus{2}), [], "omitmissing");
    % 
    % %get window and correct for edges
    % mpo_metestrus_window(1) = mpo_metestrus_idx - 5;
    % mpo_metestrus_window(2) = mpo_metestrus_idx + 6;
    % if mpo_metestrus_window(1) <= 0
    %     mpo_metestrus_window = [1 12];
    % elseif mpo_metestrus_window(2) > length(this_metestrus{1})
    %     mpo_metestrus_window = [length(this_metestrus{1})-11 length(this_metestrus{1})];
    % end
    % vmh_metestrus_window(1) = vmh_metestrus_idx - 5;
    % vmh_metestrus_window(2) = vmh_metestrus_idx + 6;
    % if vmh_metestrus_window(1) <= 0
    %     vmh_metestrus_window = [1 12];
    % elseif vmh_metestrus_window(2) > length(this_metestrus{2})
    %     vmh_metestrus_window = [length(this_metestrus{2})-11 length(this_metestrus{2})];
    % end
    % mpo_diestrus_window(1) = mpo_diestrus_idx - 5;
    % mpo_diestrus_window(2) = mpo_diestrus_idx + 6;
    % if mpo_diestrus_window(1) <= 0
    %     mpo_diestrus_window = [1 12];
    % elseif mpo_diestrus_window(2) > length(this_diestrus{1})
    %     mpo_diestrus_window = [length(this_diestrus{1})-11 length(this_diestrus{1})];
    % end
    % vmh_diestrus_window(1) = vmh_diestrus_idx - 5;
    % vmh_diestrus_window(2) = vmh_diestrus_idx + 6;
    % if vmh_diestrus_window(1) <= 0
    %     vmh_diestrus_window = [1 12];
    % elseif vmh_diestrus_window(2) > length(this_diestrus{2})
    %     vmh_diestrus_window = [length(this_diestrus{2})-11 length(this_diestrus{2})];
    % end
    % mpo_proestrus_window(1) = mpo_proestrus_idx - 5;
    % mpo_proestrus_window(2) = mpo_proestrus_idx + 6;
    % if mpo_proestrus_window(1) <= 0
    %     mpo_proestrus_window = [1 12];
    % elseif mpo_proestrus_window(2) > length(this_proestrus{1})
    %     mpo_proestrus_window = [length(this_proestrus{1})-11 length(this_proestrus{1})];
    % end
    % vmh_proestrus_window(1) = vmh_proestrus_idx - 5;
    % vmh_proestrus_window(2) = vmh_proestrus_idx + 6;
    % if vmh_proestrus_window(1) <= 0
    %     vmh_proestrus_window = [1 12];
    % elseif vmh_proestrus_window(2) > length(this_proestrus{2})
    %     vmh_proestrus_window = [length(this_proestrus{2})-11 length(this_proestrus{2})];
    % end
    % mpo_estrus_window(1) = mpo_estrus_idx - 5;
    % mpo_estrus_window(2) = mpo_estrus_idx + 6;
    % if mpo_estrus_window(1) <= 0
    %     mpo_estrus_window = [1 12];
    % elseif mpo_estrus_window(2) > length(this_estrus{1})
    %     mpo_estrus_window = [length(this_estrus{1})-11 length(this_estrus{1})];
    % end
    % vmh_estrus_window(1) = vmh_estrus_idx - 5;
    % vmh_estrus_window(2) = vmh_estrus_idx + 6;
    % if vmh_estrus_window(1) <= 0
    %     vmh_estrus_window = [1 12];
    % elseif vmh_estrus_window(2) > length(this_estrus{2})
    %     vmh_estrus_window = [length(this_estrus{2})-11 length(this_estrus{2})];
    % end

    % %get mean signal during hour with max deflection 
    % mpo_metestrus = mean(this_metestrus{1}(mpo_metestrus_window), "omitmissing");
    % vmh_metestrus = mean(this_metestrus{2}(vmh_metestrus_window), "omitmissing");
    % mpo_diestrus = mean(this_diestrus{1}(mpo_diestrus_window), "omitmissing");
    % vmh_diestrus = mean(this_diestrus{2}(vmh_diestrus_window), "omitmissing");
    % mpo_proestrus = mean(this_proestrus{1}(mpo_proestrus_window), "omitmissing");
    % vmh_proestrus = mean(this_proestrus{2}(vmh_proestrus_window), "omitmissing");
    % mpo_estrus = mean(this_estrus{1}(mpo_estrus_window), "omitmissing");
    % vmh_estrus = mean(this_estrus{2}(vmh_estrus_window), "omitmissing");

    %get mean signal during stage
    mpo_metestrus = mean(this_metestrus{1}, "omitmissing");
    vmh_metestrus = mean(this_metestrus{2}, "omitmissing");
    mpo_diestrus = mean(this_diestrus{1}, "omitmissing");
    vmh_diestrus = mean(this_diestrus{2}, "omitmissing");
    mpo_proestrus = mean(this_proestrus{1}, "omitmissing");
    vmh_proestrus = mean(this_proestrus{2}, "omitmissing");
    mpo_estrus = mean(this_estrus{1}, "omitmissing");
    vmh_estrus = mean(this_estrus{2}, "omitmissing");

    %rows
    mpo_row = {"MPOA", mpo_metestrus, mpo_diestrus, mpo_proestrus, mpo_estrus};
    vmh_row = {"VMHvl", vmh_metestrus, vmh_diestrus, vmh_proestrus, vmh_estrus};

    %input
    tbl = cat(1, tbl, mpo_row, vmh_row);
end

%run stats
withinDesign = table([1 2 3 4]', 'VariableNames', {'estrous_stage'});
mpoTbl = tbl(tbl.population == "MPOA", :);
mpo_rm = fitrm(mpoTbl, 'metestrus-estrus~1', 'WithinDesign', withinDesign);
mpo_ra = ranova(mpo_rm);
mpo_posthoc = multcompare(mpo_rm, "estrous_stage");
vmhTbl = tbl(tbl.population == "VMHvl", :);
vmh_rm = fitrm(vmhTbl, 'metestrus-estrus~1', 'WithinDesign', withinDesign);
vmh_ra = ranova(vmh_rm);
vmh_posthoc = multcompare(vmh_rm, "estrous_stage");

%plot data
stage_figure = figure("Position", [400, 1200, 520, 260], "Color", [1 1 1]);

%put mpo and vmh in separate plots
subplot(1, 2, 1);
hold on
plotEstrousData(mpoTbl)
title("MPO/MPN")
subplot(1, 2, 2);
hold on
plotEstrousData(vmhTbl)
title("VMHvl")

%save and close
exportgraphics(stage_figure, filename, "ContentType", "vector")
close all
end

%% subroutines
function plotEstrousData(tbl)
%get data
estrous_datn = zeros(height(tbl), 4);
estrous_datn(:, 1) = tbl.metestrus;
estrous_datn(:, 2) = tbl.diestrus;
estrous_datn(:, 3) = tbl.proestrus;
estrous_datn(:, 4) = tbl.estrus;

%shaded error bar on top of individual data
line([0 5], [0 0], "linestyle", "--", "linewidth", 1, "color", [0 0 0])
plot(estrous_datn', "Color", [0 0 0], "LineWidth", 0.5)
shadedErrorBar(1:4, mean(estrous_datn, 1), sem(estrous_datn, 1));

%aesthetics
estrousAx = gca;
ylim([-1.5 2])
xlim([0.25 4.75])
setAx(estrousAx);
set(estrousAx, "XTick", 1:1:4, "XTickLabels", {"met", "di", "pro", "est"})
xlabel("estrous stage")
ylabel("SeeER (zscored \DeltaF/F)")
end