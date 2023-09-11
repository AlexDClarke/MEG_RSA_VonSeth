
% Load options
load('/rds/project/rds-6yHdsDfiMLk/MEG_objects/OscRSA_2016/MEG/Subjects/meg_0247/megdata/sensortest/MEGPLANARsource_dabMefffbasic_0247.mat','option');
addpath('/rds/project/rds-6yHdsDfiMLk/MEG_objects/OscRSA_2016/scripts/scripts_jvs');

% Set delay time
delay_ms = 30;
dt = round(delay_ms/option.srate);

% Load models
optionz = option; % for parfor
load([option.var_dir option.models '.mat']);  % Models
mod_names = fieldnames(Models);
nmods = size(mod_names,1);   

% Loop over models
for mod = 1:nmods

% Loop over subjects
for sub = 1:length(optionz.subs)
    
    % Load Ant, Pos data RDMs
    posterior = load([option.datadir option.sub_beg option.subs{sub} '/megdata/sensor_post/' option.masknic{1} 'source_spatiotemporal_RDMs_60ms_sTW.mat']);
    frontal = load([option.datadir option.sub_beg option.subs{sub} '/megdata/sensor_ant/' option.masknic{1} 'source_spatiotemporal_RDMs_60ms_sTW.mat']);
    
    % Convert RDMS to vectors
    meg_data_p = [];
    for time = 1:length(posterior.ROI_RDMs(:,1,1,1))
        meg_data_tmp = squeeze(posterior.ROI_RDMs(time,:,:));  % get MEG RDM for this timepoint
        meg_data_p(:,time) = single(vectorizeRDM(meg_data_tmp)'); clear meg_data_tmp
    end
    meg_data_a = [];
    for time = 1:length(frontal.ROI_RDMs(:,1,1,1))
        meg_data_tmp = squeeze(frontal.ROI_RDMs(time,:,:));  % get MEG RDM for this timepoint
        meg_data_a(:,time) = single(vectorizeRDM(meg_data_tmp)'); clear meg_data_tmp
    end
    models = vectorizeRDM(Models.(mod_names{mod}))';
    
    % Reject incorrect trials
    x=meg_data_p;
    meg_data_p(isnan(x(:,200)),:) = [];
    meg_data_a(isnan(x(:,200)),:) = [];
    models(isnan(x(:,200)),:) = [];
    clear x
    
    % Loop through time, creating FF and FB timecourses
    for t = (dt+1):length(meg_data_p(1,:))
        % FF
        out = squeeze(mean(meg_data_p(:,t-dt:t-1),2));
        in = squeeze(meg_data_a(:,t));
        r(t,:) = corr(models,in,'type',option.dist,'rows','pairwise');
        pc = partialcorr(models,in,out,'type', option.dist, 'rows', 'pairwise')';
        rp(t,:) = pc;
        
        % FB
        out = squeeze(mean(meg_data_a(:,t-dt:t-1),2));
        in = squeeze(meg_data_p(:,t));
        r2(t,:) = corr(models,in,'type',option.dist,'rows','pairwise');
        pc = partialcorr(models,in,out,'type', option.dist, 'rows', 'pairwise')';
        rp2(t,:) = pc;        
    end
    
    ff = r - rp;
    fb = r2 - rp2;
    
    % Save FF and FB responses for each subject
    
    outfile_ff = ([option.datadir option.sub_beg option.subs{sub} '/megdata/' 'connect/' option.masknic{1} mod_names{mod} 'ff.mat']); 
    outfile_fb = ([option.datadir option.sub_beg option.subs{sub} '/megdata/' 'connect/' option.masknic{1} mod_names{mod} 'fb.mat']);
    
    save(outfile_ff,'ff');
    save(outfile_fb,'fb');
    
    % Structure for all participants
    ff_all(:,sub) = ff;
    fb_all(:,sub) = fb;
    
end 
    
% Average over subjects    
mean_ff = mean(ff_all,2);
mean_fb = mean(fb_all,2); 

% Save ff/fb for all subjects
outfile_ff_all = ([option.datadir 'all/' option.masknic{1} mod_names{mod} 'FF_FB.fig']);
save(outfile_ff,'ff_all');

outfile_fb_all = ([option.datadir 'all/' option.masknic{1} mod_names{mod} 'FF_FB.fig']);
save(outfile_fb,'fb_all');

% Plot
figure;
plot(ff, 'LineWidth', 2);
hold on; 
plot(fb, 'LineWidth', 2);
set(gca,'XTick',[1:option.jump/option.srate:option.epoch_length/option.srate]);
set(gca,'XTickLabel',[round(option.start):option.jump:round(option.stop)]);
set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')
set(get(gca,'YLabel'),'String','Information Flow', 'fontweight','b')
xlim(option.range)
legend('feedforward','feedback') 
outfile_fig = ([option.datadir 'figures/' option.masknic{1} mod_names{mod} 'FF_FB.fig']);
savefig(outfile_fig)

end