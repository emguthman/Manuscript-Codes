function [signal, ts] = equalize_vector_lengths(signal, ts)
if (length(ts.gcamp.sec) ~= length(ts.rcamp.sec)) || (length(ts.gcamp.sec) ~= length(ts.isos.sec)) || (length(ts.rcamp.sec) ~= length(ts.isos.sec))
    if length(ts.gcamp.sec) > length(ts.rcamp.sec) %if gcamp has more frames than rcamp
        if length(ts.gcamp.sec) > length(ts.isos.sec) %if gcamp has more frames than isos
            signal.gcamp.raw(:,1) = []; signal.rcamp.gcamp(:,1) = []; 
            ts.gcamp.sec(1) = []; ts.gcamp.ind(end) = [];
        else %if equal # of gcamp and isos frames, but more than rcamp
            signal.gcamp.raw(:,1) = []; signal.rcamp.gcamp(:,1) = []; signal.gcamp.isos(:,1) = []; signal.rcamp.isos(:,1) = []; 
            ts.gcamp.sec(1) = []; ts.gcamp.ind(end) = []; ts.isos.sec(1) = []; ts.isos.ind(end) = [];
        end
    elseif length(ts.rcamp.sec) > length(ts.isos.sec) %if rcamp has more frames than isos
        if length(ts.rcamp.sec) > length(ts.gcamp.sec) %if rcamp has more frames than gcamp
            signal.rcamp.raw(:,1) = []; signal.gcamp.rcamp(:,1) = [];
            ts.rcamp.sec(1) = []; ts.rcamp.ind(end) = [];
        else %if equal # of gcamp and rcamp frames, but more than isos
            signal.rcamp.raw(:,1) = []; signal.gcamp.rcamp(:,1) = []; signal.gcamp.raw(:,1) = []; signal.rcamp.gcamp(:,1) = [];
            ts.rcamp.sec(1) = []; ts.rcamp.ind(end) = []; ts.gcamp.sec(1) = []; ts.gcamp.ind(end) = [];
        end
    elseif length(ts.isos.sec) > length(ts.gcamp.sec) %if isos has more frames than gcamp
        if length(ts.isos.sec) > length(ts.rcamp.sec) %if isos has more frames than rcamp
            signal.gcamp.isos(:,1) = []; signal.rcamp.isos(:,1) = [];
            ts.isos.sec(1) = []; ts.isos.ind(end) = [];
        else %if equal # of isos and rcamp frames, but more than gcamp
            signal.gcamp.isos(:,1) = []; signal.rcamp.isos(:,1) = []; signal.rcamp.raw(:,1) = []; signal.gcamp.rcamp(:,1) = [];
            ts.isos.sec(1) = []; ts.isos.ind(end) = []; ts.rcamp.sec(1) = []; ts.rcamp.ind(end) = [];
        end
    end
end
end