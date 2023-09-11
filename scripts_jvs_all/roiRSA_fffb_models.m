
% Load options
load('/rds/project/rds-6yHdsDfiMLk/MEG_objects/OscRSA_2016/MEG/Subjects/meg_0247/megdata/sensortest/MEGPLANARsource_dabMefffbasic_0247.mat','option')

% Load models
optionz = option; % for parfor
load([option.var_dir option.models '.mat']);  % Models
mod_names = fieldnames(Models);
nmods = size(mod_names,1);

% Loop over subjects
for sub = 2%:length(optionz.subs)
    
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
    models = vectorizeRDM(Models.C1feats)';
    
    % Reject incorrect trials
    x=meg_data_p;
    meg_data_p(isnan(x(:,200)),:) = [];
    meg_data_a(isnan(x(:,200)),:) = [];
    models(isnan(x(:,200)),:) = [];
    clear x
    
    % Loop through time, creating FF and FB timecourses
    for t = 16:length(meg_data_p(1,:))
        % FF
        out = squeeze(mean(meg_data_p(:,t-15:t-1),2));
        in = squeeze(meg_data_a(:,t));
        r(t,:) = corr(models,in,'type',option.dist,'rows','pairwise');
        pc = partialcorr(models,in,out,'type', option.dist, 'rows', 'pairwise')';
        rp(t,:) = pc;
        
        % FB
        out = squeeze(mean(meg_data_a(:,t-15:t-1),2));
        in = squeeze(meg_data_p(:,t));
        r2(t,:) = corr(models,in,'type',option.dist,'rows','pairwise');
        pc = partialcorr(models,in,out,'type', option.dist, 'rows', 'pairwise')';
        rp2(t,:) = pc;        
    end
    
    ff = r - rp;
    fb = r2 - rp2;
    
    % Save FF and FB responses
    
    outfile_ff = ([option.datadir option.sub_beg option.subs{sub} '/megdata/' option.masknic{1} 'ff.mat']);
    outfile_fb = ([option.datadir option.sub_beg option.subs{sub} '/megdata/' option.masknic{1} 'fb.mat']);
    
    save(outfile_ff,'ff');
    save(outfile_fb,'fb');
    
end 
    
    % Collate over subjects
