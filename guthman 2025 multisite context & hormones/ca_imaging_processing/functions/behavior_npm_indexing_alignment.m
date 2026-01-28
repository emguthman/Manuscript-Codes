function indexing = behavior_npm_indexing_alignment(ts)
% matching behavior frame to npm frame (find nearest beahavior frame)
for ii = 1:length(ts.rcamp.sec)
    [~, idx] = min(abs(ts.rcamp.sec(ii) - ts.behavior_cam.sec));
    indexing.rcamp_to_behv(:, ii) = [ts.behavior_cam.ind(idx) ts.behavior_cam.sec(idx)];
    [~, idx] = min(abs(ts.gcamp.sec(ii) - ts.behavior_cam.sec));
    indexing.gcamp_to_behv(:, ii) = [ts.behavior_cam.ind(idx) ts.behavior_cam.sec(idx)];
    [~, idx] = min(abs(ts.isos.sec(ii) - ts.behavior_cam.sec));
    indexing.isos_to_behv(:, ii) = [ts.behavior_cam.ind(idx) ts.behavior_cam.sec(idx)];
end

% matching npm frame to behavior frame (find nearest npm frame)
for ii = 1:length(ts.behavior_cam.sec)
    [~, idx] = min(abs(ts.behavior_cam.sec(ii) - ts.rcamp.sec));
    indexing.behv_to_npm(:, ii) = [ts.rcamp.ind(idx) ts.rcamp.sec(idx)];
end
end