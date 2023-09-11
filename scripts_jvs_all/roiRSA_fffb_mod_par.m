
% Load options
load('/rds/project/rds-6yHdsDfiMLk/MEG_objects/OscRSA_2016/MEG/Subjects/meg_0247/megdata/sens_post_4conn_gap/MEGPLANARsource_pdabMefffbasic_0247.mat','option');
addpath('/rds/project/rds-6yHdsDfiMLk/MEG_objects/OscRSA_2016/scripts/scripts_jvs_all');

% Set delay time
delay_ms = 30;
dt = round(delay_ms/option.srate);

% Load models
optionz = option; % for parfor
load([option.var_dir 'RSA_sim_modelRDMs_cornet_semfeats' '.mat']);  % Models
mod_names = fieldnames(Models);
nmods = size(mod_names,1); 

% Loop over models
for mod = 1:nmods 

% Loop over subjects
parfor sub = 1:length(optionz.subs)
    
    % Load Ant, Pos data RDMs
    posterior = load([option.datadir option.sub_beg option.subs{sub} '/megdata/sens_post_4conn_gap/' option.masknic{1} 'source_spatiotemporal_RDMs_2ms_sTW.mat']);
    frontal = load([option.datadir option.sub_beg option.subs{sub} '/megdata/sens_ant_4conn_gap/' option.masknic{1} 'source_spatiotemporal_RDMs_2ms_sTW.mat']);
    
    % Convert RDMS to vectors
    dimRDM = length(posterior.ROI_RDMs(1,1,:,1));
    dimTime = length(posterior.ROI_RDMs(:,1,1,1));
    meg_data_p = zeros(((dimRDM*dimRDM)/2)-(dimRDM/2),dimTime);
    for time = 1:length(posterior.ROI_RDMs(:,1,1,1))
        meg_data_tmp = squeeze(posterior.ROI_RDMs(time,:,:));  % get MEG RDM for this timepoint
        meg_data_p(:,time) = single(vectorizeRDM(meg_data_tmp)'); % clear meg_data_tmp
    end
    meg_data_a = zeros(((dimRDM*dimRDM)/2)-(dimRDM/2),dimTime);
    for time = 1:length(frontal.ROI_RDMs(:,1,1,1))
        meg_data_tmp = squeeze(frontal.ROI_RDMs(time,:,:));  % get MEG RDM for this timepoint
        meg_data_a(:,time) = single(vectorizeRDM(meg_data_tmp)'); % clear meg_data_tmp
    end
    models = vectorizeRDM(Models.(mod_names{mod}))';
    
    % for parcorr
    % Vect
    for modll = 1:nmods    
        model_tmp = Models.(mod_names{modll});
        model(:,modll) = vectorizeRDM(model_tmp)';
    end
    
    model(:,mod) = [];
    %out_mod = out_mod';
    
    % Reject incorrect trials
    x=meg_data_p;
    meg_data_p(isnan(x(:,200)),:) = [];
    meg_data_a(isnan(x(:,200)),:) = [];
    models(isnan(x(:,200)),:) = [];
    
    exc = [];
    for i = 1:min(size(model))
        y = model(:,i);
        y(isnan(x(:,200)),:) = [];
        exc(:,i)= y;
    end
    x=[]; y =[]; model = [];
%   clear x
%   clear y
%   clear model
    
    % Loop through time, creating FF and FB timecourses
    r = zeros(length(meg_data_p(1,:)),1);
    r2 = zeros(length(meg_data_p(1,:)),1);
    rp = zeros(length(meg_data_p(1,:)),1);
    rp2 = zeros(length(meg_data_p(1,:)),1);

    for t = (dt+1):length(meg_data_p(1,:))
        % FF
        out = squeeze(mean(meg_data_p(:,t-dt:t-1),2));
        in = squeeze(meg_data_a(:,t));
        exc_out = [exc out];
        rpc = partialcorr(models,in,exc,'type',option.dist,'rows','pairwise');
        r(t,:) = rpc;        
        pc = partialcorr(models,in,exc_out,'type', option.dist, 'rows', 'pairwise')';
        rp(t,:) = pc;
        
        % FB
        out = squeeze(mean(meg_data_a(:,t-dt:t-1),2));
        in = squeeze(meg_data_p(:,t));
        exc_out = [exc out];
        r2p = partialcorr(models,in,exc,'type',option.dist,'rows','pairwise');
        r2(t,:) = r2p;
        pc = partialcorr(models,in,exc_out,'type', option.dist, 'rows', 'pairwise')';
        rp2(t,:) = pc;        
    end
    
    ff = r - rp;
    fb = r2 - rp2;
    
    % Save FF and FB responses for each subject    
    outfile_ff = ([option.datadir option.sub_beg option.subs{sub} '/megdata/connect_sensor_4conn_gap/' option.masknic{1} mod_names{mod} 'ff_par.mat']); 
    outfile_fb = ([option.datadir option.sub_beg option.subs{sub} '/megdata/connect_sensor_4conn_gap/' option.masknic{1} mod_names{mod} 'fb_par.mat']);
    
    mkdir([option.datadir option.sub_beg option.subs{sub} '/megdata/connect_sensor_4conn_gap/'])

%   save(outfile_ff,'ff');
%   save(outfile_fb,'fb');
    parsave2(outfile_ff,ff);
    parsave3(outfile_fb,fb);
    
%    % Structure for all participants
%    ff_all(:,sub) = ff;
%    fb_all(:,sub) = fb;
    
%   clear exc
    
end
end

option.rfxdir = '/rds/project/rds-6yHdsDfiMLk/MEG_objects/OscRSA_2016/MEG/RSA_ROI/connect_sensor_4conn_gap/';

% Compile subject data
for sub = 1:length(option.subs)
    for mod = 1:size(mod_names,1);   

    load([option.datadir option.sub_beg option.subs{sub} '/megdata/connect_sensor_4conn_gap/' option.masknic{1} mod_names{mod} 'ff_par.mat']);     
    load([option.datadir option.sub_beg option.subs{sub} '/megdata/connect_sensor_4conn_gap/' option.masknic{1} mod_names{mod} 'fb_par.mat']);
    subff(1,sub,mod,:) = ff;
    subfb(1,sub,mod,:) = fb;    

    end
end

for m = 1:size(mod_names,1);
    rsa_out.(mod_names{m}) = squeeze(subff(1,:,m,:));  % store RSA timecourses    
end
option.rsafront = 'RS_connectivity_FF_';
outname = [option.rfxdir option.rsafront option.masknic{1} option.midname num2str(option.tw*option.srate) 'ms_sTW.mat'];
save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out

for m = 1:size(mod_names,1);
    rsa_out.(mod_names{m}) = squeeze(subfb(1,:,m,:));  % store RSA timecourses    
end
option.rsafront = 'RS_connectivity_FB_';
outname = [option.rfxdir option.rsafront option.masknic{1} option.midname num2str(option.tw*option.srate) 'ms_sTW.mat'];
save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out


% % Average over subjects    
% mean_ff = mean(ff_all,2);
% mean_fb = mean(fb_all,2); 

% % Save ff/fb for all subjects
% outfile_ff_all = ([option.datadir 'all/' option.masknic{1} mod_names{mod} 'FF_par.mat']);
% save(outfile_ff_all,'ff_all');
% 
% outfile_fb_all = ([option.datadir 'all/' option.masknic{1} mod_names{mod} 'FB_par.mat']);
% save(outfile_fb_all,'fb_all');

% % Plot
% figure;
% plot(ff, 'LineWidth', 2);
% hold on; 
% plot(fb, 'LineWidth', 2);
% set(gca,'XTick',[1:option.jump/option.srate:option.epoch_length/option.srate]);
% set(gca,'XTickLabel',[round(option.start):option.jump:round(option.stop)]);
% set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')
% set(get(gca,'YLabel'),'String','Information Flow', 'fontweight','b')
% xlim(option.range)
% legend('feedforward','feedback') 
% outfile_fig = ([option.datadir 'figures/' option.masknic{1} mod_names{mod} 'FF_FB_partial.fig']);
% savefig(outfile_fig)
