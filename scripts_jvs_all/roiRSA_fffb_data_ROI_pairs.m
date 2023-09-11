% Partial correlations connectivity between ROIs, combined hemispheres,

% Load options
%load('/rds/project/rds-6yHdsDfiMLk/MEG_objects/OscRSA_2016/MEG/Subjects/meg_0247/megdata/source/R_V13source_wmabMefffbasic_0247.mat','option');
%addpath('/rds/project/rds-6yHdsDfiMLk/MEG_objects/OscRSA_2016/scripts/scripts_ROI_all');
load('/rds/project/rds-6yHdsDfiMLk/MEG_objects/OscRSA_2016/MEG/Subjects/meg_0247/megdata/sens_post_4conn_gap/MEGPLANARsource_pdabMefffbasic_0247.mat','option');
addpath('/rds/project/rds-6yHdsDfiMLk/MEG_objects/OscRSA_2016/scripts/scripts_jvs_all');

% Set delay time
delay_ms = 30;
dt = round(delay_ms/option.srate);

% Load models
optionz = option; % for parfor

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
    
    % Reject incorrect trials
    x=meg_data_p;
    meg_data_p(isnan(x(:,200)),:) = [];
    meg_data_a(isnan(x(:,200)),:) = [];
        
    x=[]; ff =[]; fb = [];
    
    % Loop through time, creating FF and FB timecourses
    r = zeros(length(meg_data_p(1,:)),1);
    r2 = zeros(length(meg_data_p(1,:)),1);
    rp = zeros(length(meg_data_p(1,:)),1);
    rp2 = zeros(length(meg_data_p(1,:)),1);
    
    % Loop through time
    for t = (dt+1):length(meg_data_p(1,:))
        % FF (P->A)
        out = squeeze(mean(meg_data_p(:,t-dt:t-1),2));
        exc = squeeze(mean(meg_data_a(:,t-dt:t-1),2));
        in = squeeze(meg_data_a(:,t));       
        rpc = partialcorr(in,out,exc,'type',option.dist,'rows','pairwise');
        ff(t,:) = rpc;

        % FB (A->P)
        out = squeeze(mean(meg_data_a(:,t-dt:t-1),2));
        exc = squeeze(mean(meg_data_p(:,t-dt:t-1),2));
        in = squeeze(meg_data_p(:,t));       
        rpc = partialcorr(in,out,exc,'type',option.dist,'rows','pairwise');
        fb(t,:) = rpc;
    end
       
    % Save FF and FB responses for each subject    
    outfile_ff = ([option.datadir option.sub_beg option.subs{sub} '/megdata/connect_sensor_4conn_gap/' option.masknic{1} '_data_ff_par.mat']); 
    outfile_fb = ([option.datadir option.sub_beg option.subs{sub} '/megdata/connect_sensor_4conn_gap/' option.masknic{1} '_data_fb_par.mat']);

    mkdir([option.datadir option.sub_beg option.subs{sub} '/megdata/connect_sensor_4conn_gap/'])

    parsave2(outfile_ff,ff);
    parsave3(outfile_fb,fb);
    
end 

% Compile subject data    
option.rfxdir = '/rds/project/rds-6yHdsDfiMLk/MEG_objects/OscRSA_2016/MEG/RSA_ROI/connect_sensor_4conn_gap/';
if ~exist(option.rfxdir)
    mkdir(option.rfxdir)
end

for sub = 1:length(option.subs)
    load([option.datadir option.sub_beg option.subs{sub} '/megdata/connect_sensor_4conn_gap/' option.masknic{1} '_data_ff_par.mat']);     
    load([option.datadir option.sub_beg option.subs{sub} '/megdata/connect_sensor_4conn_gap/' option.masknic{1} '_data_fb_par.mat']);
    subff(1,sub,:) = ff;
    subfb(1,sub,:) = fb;    
end

rsa_out.data = squeeze(subff(1,:,:));  % store RSA timecourses    
option.rsafront = 'RS_data_connectivity_FF__';
outname = [option.rfxdir option.rsafront option.masknic{1} option.midname num2str(option.tw*option.srate) 'ms_sTW.mat'];
save(outname, 'rsa_out', 'option'); clear rsa_out

rsa_out.data = squeeze(subfb(1,:,:));  % store RSA timecourses    
option.rsafront = 'RS_data_connectivity_FB_';
outname = [option.rfxdir option.rsafront option.masknic{1} option.midname num2str(option.tw*option.srate) 'ms_sTW.mat'];
save(outname, 'rsa_out', 'option'); clear rsa_out
