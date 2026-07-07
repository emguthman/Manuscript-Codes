function [dff, z] = compute_estrous_dff(signal_obj)
%inits
fs = 1/12;
n_fibers  = 2;
M = height(signal_obj{1});
dff{1} = NaN(size(signal_obj{1}));
dff{2} = NaN(size(signal_obj{2}));
z{1} = NaN(size(signal_obj{1}));
z{2} = NaN(size(signal_obj{2}));

for moyz = 1:M
    for versye = 1:n_fibers
        %get baseline data
        tmp_baseline = signal_obj{versye}(moyz, 1:(fs^(-1)*24));

        %get dff
        % dff{versye}(moyz, :) = 100 + 100 .* (signal_obj{versye}(moyz, :) - mean(tmp_baseline, "omitmissing")) ./ mean(tmp_baseline, "omitmissing");
        dff{versye}(moyz, :) = 100 + 100 .* (signal_obj{versye}(moyz, :) - mean(signal_obj{versye}(moyz, :), "omitmissing")) ./ mean(signal_obj{versye}(moyz, :), "omitmissing");
        z{versye}(moyz, :) = (signal_obj{versye}(moyz, :) - mean(signal_obj{versye}(moyz, :), "omitmissing")) ./ std(signal_obj{versye}(moyz, :), [], "omitmissing");
    end
end
end