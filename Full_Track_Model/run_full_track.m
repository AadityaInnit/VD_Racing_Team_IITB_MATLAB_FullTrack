%% Independent four-wheel RWD full-track demonstration.
folder = fileparts(mfilename('fullpath')); addpath(folder);
c = full_track_config();
disp(c.assumptions);
% Replace these points with an ordered XY driving line in metres.
theta = linspace(0,2*pi,81).';
xy = [35*sin(theta),22*(1-cos(theta))];
track = make_full_track(xy,true,600);
cache = fullfile(folder,'data','four_wheel_ggv.mat');
% A mismatched cache is rejected; set rebuildGGV=true after changing physics.
rebuildGGV = ~isfile(cache);
if rebuildGGV
    ggv = full_track_model(c,[],BuildGGVOnly=true);
    save(cache,'ggv');
else
    data = load(cache,'ggv'); ggv = data.ggv;
end
out = full_track_model(c,track,GGV=ggv);
save(fullfile(folder,'full_track_results.mat'),'out');
if c.makePlots
    exportgraphics(gcf,fullfile(folder,'full_track_overview.png'),'Resolution',150);
end
fprintf('%s | simulated %.3f s | plan %.3f s | max path error %.3f m\n', ...
    out.status,out.lapTime,out.plan.lapTime,max(abs(out.e_y)));
fprintf('Tyre-fit range exceeded in %.1f%% of recorded samples.\n',100*mean(out.tireFitRangeExceeded));
