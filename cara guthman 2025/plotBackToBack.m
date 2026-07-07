function plotBackToBack(signal, red)
%inits
F = figure("Color", [1 1 1]);
hold on

for s = 1:length(signal)
    % %detrend signals
    % detrended_green_signal{s} = detrend(signal{s}) + mean(signal{s});
    % detrended_red_signal{s} = detrend(red{s}) + mean(red{s});
    % 
    % %regress signals
    % beta = robustfit(detrended_red_signal{s}, detrended_green_signal{s});
    % regressed_signal{s} = detrended_green_signal{s} - beta(2)*detrended_red_signal{s} + beta(1);
    % regressed_signal{s} = regressed_signal{s} + abs(min(regressed_signal{s})) + mean(detrended_green_signal{s});
    % 
    % %whole trace dff
    % dff_signal{s} = 100 + 100 .* (regressed_signal{s} - mean(regressed_signal{s}))/mean(regressed_signal{s});
    % 
    % %low pass filter data
    % filtered_signal{s} = lowpass_filter_signal(dff_signal{1}, 2, 12, 12 * 12);

    %plot
    plot(signal{s})
    % plot(100 + 100 .* (signal{s} - mean(signal{s}))/mean(signal{s}))
    % plot(dff_signal{s})
    % plot(filtered_signal{s})
end

%aesthetics
b2bAx = gca;
b2bAx.XTick = 0:12*12:8*12*24;

%save
exportgraphics(F, "back_to_back_estrous.pdf", "ContentType", "vector");

end

function trans_signal = lowpass_filter_signal(signal, cutoff, fs, N)
%pad signals
padded_signal = [flip(signal(1:N)), signal, flip(signal((end - N + 1):end))]'; %reflective padding

%filter
[b, a] = butter(8, (cutoff)/(fs / 2), 'low');
butter_signal = filtfilt(b, a, padded_signal);

%remove padding
trans_signal = butter_signal(N+1:N+length(signal));
end
