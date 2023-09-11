function roiRSA_timecourses(option,s)

% Correlates each subjects data RDMs with model RDMs and saves RSA timecourses
%
% AC (11/2012)(10/2014) Feb 2015
% This may be out of date - use with caution


% Begin
for mask = 1:length(option.masknic)
    
    masknam = [option.masknic{mask}];
    
    % Load model to get model number
    load([option.var_dir option.models '.mat']);  % Models
    mod_names = fieldnames(Models);
    nmods = size(mod_names,1);
    
%   tmp_out = zeros(length(option.subs),nmods,option.epoch_length/option.srate,'single');     
    if option.doTFbands
        %   start = (1000/min(option.fsc))*(round(option.cycles/2));
        tmp_out = zeros(length(option.subs),nmods,size(option.fs,2),(option.epoch_length)/option.tfstep+1);
    elseif option.doTF
        %   start = (1000/min(option.fsc))*(round(option.cycles/2));
        tmp_out = zeros(length(option.subs),nmods,option.nfs,(option.epoch_length)/option.tfstep+1);
    else
        tmp_out = zeros(length(option.subs),nmods,option.epoch_length/option.srate+1);
    end
    
    for sub = 1:length(option.subs)
        
        cd([option.datadir option.sub_beg option.subs{sub} option.subdir]);
        
        options = optionsfile(s,1);
        option = setfield(option,'models',options.models);
        clear options
        
        % Load subject-specific model RDMs
        load([option.var_dir option.models '.mat']);  % Models
        mod_names = fieldnames(Models);
        nmods = size(mod_names,1);
        for m = 1:nmods
            model_tmp = Models.(mod_names{m});
            model(:,m) = vectorizeRDM(model_tmp)'; clear model_tmp
        end
        
        sprintf('......Subject %s, Region %s......', option.subs{sub},num2str(mask))
        
        % Load data RDMs
        if option.doTF && option.doTFbands==0
            if sub == 1 && mask == 1
            option.tw = round((option.tw*option.srate)/option.tfstep);
            end
            if option.tw > 1
                if mod(option.tw,2)  % must be even number
                    option.tw = option.tw+1;
                end
            end
            if option.doavg
                if option.doPhase
                    rdms = matfile([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                else
                    rdms = matfile([option.tf_pre option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                end
            else
                if option.doPhase
                    rdms = matfile([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                else
                    rdms = matfile([option.tf_pre option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                end
            end
        elseif option.doTFbands
            if sub == 1 && mask == 1
            if option.doTWs
                option.tw = 1;
            else
                option.tw = round((option.tw*option.srate)/option.tfstep);
            end
            end
            if option.tw > 1
                if mod(option.tw,2)  % must be even number
                    option.tw = option.tw+1;
                end
            end
            if option.doavg
                if option.doPhase
                    rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                else
                    rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                end
            else
                if option.doPhase
                    if option.doTWs
                        rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'TW_RDMs.mat']);
                    else
                        rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                    end
                else
                    if option.doTWs
                        rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'TW_RDMs.mat']);
                    else
                        rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                    end
                end
            end
        else
            if option.tw > 1
                if mod(option.tw,2)  % must be even number
                    option.tw = option.tw+1;
                end
            end
            if option.doavg
                rdms = matfile([option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat']);
            else
                rdms = matfile([option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat']);
            end
        end
        
        % Begin model-brain RSA
        if or(option.doTF,option.doTFbands)
            if option.parmodels  % run all models together
                vs = length(rdms.ROI_RDMs(1,1,:,1));
                meg_data = zeros((((vs*vs)-vs)/2),length(rdms.ROI_RDMs(1,:,1,1)),length(rdms.ROI_RDMs(:,1,1,1)),'single');
                for f = 1:length(rdms.ROI_RDMs(1,:,1,1))
                    for time = 1:length(rdms.ROI_RDMs(:,1,1,1))
                        meg_data_tmp = squeeze(rdms.ROI_RDMs(time,f,1:vs,1:vs));  % get MEG RDM for this timepoint
                        meg_data(:,f,time) = single(vectorizeRDM(meg_data_tmp)'); clear meg_data_tmp
                    end
                end % f
                meg_data = reshape(meg_data,size(meg_data,1),[]);
                         
                % Remove effect of time
                if option.remtime
                    warning off
                    sprintf('......Removing effect of time......')
                    times= option.timevector;
                    drift_reg = [];
                    for n = 1:option.ndrifts
                        % Drift models
                        drift_reg = [drift_reg, (pdist(times').^n)'];
                    end
                    % Remove NaNs
                    x = [];
                    x = [model, meg_data];  clear meg_data
                    drift_reg(any(isnan(x(:,[1:nmods round(size(x,2)/2)])),2),:) = [];
                    x(any(isnan(x(:,[1:nmods round(size(x,2)/2)])),2),:) = [];
                    for i = nmods+1:size(x,2)
                        [~,~,r] = regress(x(:,i),[ones(size(x,1),1) drift_reg]);
                        x(:,i) = r;
                    end % i
                    clear drift_reg r i
                    warning on
                else                    
                    % Remove NaNs
                    x = [];
                    x = [model, meg_data]; clear meg_data % combine models and data in same matrix
                    x(any(isnan(x(:,[1:nmods round(size(x,2)/2)])),2),:) = [];
                end
                
                % correlate
                sprintf('......model correlations......')
                y = x(:,1:nmods); y = double(y);  % seperate for memory saving
                x(:,1:nmods) = []; x = double(x);
                r = corr(y,x,'type',option.dist,'rows','pairwise');
                r=single(r);
                tmp_out(sub,:,:,:) = single(reshape(r,nmods,length(rdms.ROI_RDMs(1,:,1,1)),length(rdms.ROI_RDMs(:,1,1,1))));
            
            else % run each model seperately
                for m = 1:nmods
                    vs = length(rdms.ROI_RDMs(1,1,:,1));
                    meg_data = zeros((((vs*vs)-vs)/2),length(rdms.ROI_RDMs(1,:,1,1)),length(rdms.ROI_RDMs(:,1,1,1)),'single');
                    for f = 1:length(rdms.ROI_RDMs(1,:,1,1))
                        for time = 1:length(rdms.ROI_RDMs(:,1,1,1))
                            meg_data_tmp = squeeze(rdms.ROI_RDMs(time,f,1:vs,1:vs));  % get MEG RDM for this timepoint
                            meg_data(:,f,time) = single(vectorizeRDM(meg_data_tmp)'); clear meg_data_tmp
                        end
                    end % f                    
                    meg_data = reshape(meg_data,size(meg_data,1),[]);
                    
                    % Remove effect of time
                    if option.remtime
                        warning off
                        sprintf('......Removing effect of time......')
                        times=option.timevector;
                        drift_reg = [];
                        for n = 1:option.ndrifts
                            % Drift models
                            drift_reg = [drift_reg, (pdist(times').^n)'];
                        end
                        % Remove NaNs
                        x = [];
                        x = [model(:,m), meg_data];  clear meg_data
                        drift_reg(any(isnan(x(:,[1 round(size(x,2)/2)])),2),:) = [];
                        x(any(isnan(x(:,[1 round(size(x,2)/2)])),2),:) = [];
                        for i = 2:size(x,2)
                            [~,~,r] = regress(x(:,i),[ones(size(x,1),1) drift_reg]);
                            x(:,i) = r;
                        end % i
                        clear drift_reg r i
                        warning on
                    else
                        % Remove NaNs
                        x = [];
                        x = [model(:,m), meg_data]; clear meg_data % combine models and data in same matrix
                        x(any(isnan(x(:,[1 round(size(x,2)/2)])),2),:) = []; % only check model and first time for NaNs                        
                    end
                    
                    % correlate
                    sprintf('......model correlations......')
                    y=double(x(:,1)); % seperate model and data for memory saving
                    x(:,1) = [];
                    x=double(x);
                    tmpr = corr(y,x,'type',option.dist,'rows','pairwise');
                    tmp_out(sub,m,:,:) = single(reshape(tmpr,length(rdms.ROI_RDMs(1,:,1,1)),length(rdms.ROI_RDMs(:,1,1,1)))); clear tmpr x y
                end                
            end
            
            % Update options            
            options = rdms.option;
            option.epoch_length = options.epoch_length;
            option.baseline_length = options.baseline_length;
            option.tfepoch = options.tfepoch;
            option.baseline = options.baseline;
            option.fsc = options.fsc;
            clear options
            
        else
            vs = length(rdms.ROI_RDMs(1,1,:,1));
            meg_data = zeros((((vs*vs)-vs)/2),length(rdms.ROI_RDMs(:,1,1,1)),'single');
            for time = 1:length(rdms.ROI_RDMs(:,1,1))
                meg_data_tmp = squeeze(rdms.ROI_RDMs(time,1:vs,1:vs));  % get MEG RDM for this timepoint
                meg_data(:,time) = single(vectorizeRDM(meg_data_tmp)'); clear meg_data_tmp
            end
            
            if option.parmodels  % run all models together
                                
                % Remove effect of time
                if option.remtime
                    warning off
                    sprintf('......Removing effect of time......')
                    times=option.timevector;
                    drift_reg = [];
                    for n = 1:option.ndrifts
                        % Drift models
                        drift_reg = [drift_reg, (pdist(times').^n)'];
                    end
                    % Remove NaNs
                    x = [];
                    x = [model, meg_data];  clear meg_data
                    drift_reg(any(isnan(x(:,[1:nmods round(size(x,2)/2)])),2),:) = [];
                    x(any(isnan(x(:,[1:nmods round(size(x,2)/2)])),2),:) = [];
                    for i = nmods+1:size(x,2)
                        [~,~,r] = regress(x(:,i),[ones(size(x,1),1) drift_reg]);
                        x(:,i) = r;
                    end % i
                    clear drift_reg r i
                    warning on
                else                    
                    % Remove NaNs
                    x = [];
                    x = [model, meg_data]; clear meg_data % combine models and data in same matrix
                    x(any(isnan(x(:,[1:nmods round(size(x,2)/2)])),2),:) = [];
                end                                              
                
                % correlate
                sprintf('......model correlations......')
                y = x(:,1:nmods); y = double(y);  % seperate for memory saving
                x(:,1:nmods) = []; x = double(x);
                
                if option.partial
                    % partial correlation controlling for effect of all models on one model
                    tot_mods = 1:nmods;
                    for i = tot_mods
                        % Get model sets
                        exc = tot_mods; inc = tot_mods(i); exc(i) = [];
                        pc = partialcorr(x,y(:,inc),y(:,exc),'type', option.dist, 'rows', 'pairwise')';
                        r(i,:) = pc;
                    end
                else
                    r = corr(y,x,'type',option.dist,'rows','pairwise');
                end                
                r = single(r);                
                tmp_out(sub,:,:) = single(r);
            else % run each model seperately
                for m = 1:nmods
                    % Remove effect of time
                    if option.remtime
                        warning off
                        sprintf('......Removing effect of time......')
                        times=option.timevector;
                        drift_reg = [];
                        for n = 1:option.ndrifts
                            % Drift models
                            drift_reg = [drift_reg, (pdist(times').^n)'];
                        end
                        % Remove NaNs
                        x = [];
                        x = [model(:,m), meg_data];
                        drift_reg(any(isnan(x(:,[1 round(size(x,2)/2)])),2),:) = [];
                        x(any(isnan(x(:,[1 round(size(x,2)/2)])),2),:) = [];
                        for i = 2:size(x,2)
                            [~,~,r] = regress(x(:,i),[ones(size(x,1),1) drift_reg]);
                            x(:,i) = r;
                        end % i
                        clear drift_reg r i
                        warning on
                    else
                        % Remove NaNs
                        x = [];
                        x = [model(:,m), meg_data]; % combine models and data in same matrix
                        x(any(isnan(x(:,[1 round(size(x,2)/2)])),2),:) = [];
                    end
                    
                    % correlate
                    sprintf('......model correlations......')
                    tmpr = corr(double(x(:,1)),double(x(:,2:end)),'type',option.dist,'rows','pairwise');
                    tmp_out(sub,m,:) = single(tmpr); clear tmpr
                end
            end
        end
    end
    
    for m = 1:nmods
        rsa_out.(mod_names{m}) = squeeze(tmp_out(:,m,:,:));  % store RSA timecourses
    end
    clear tmp_out
    
    options = optionsfile(s,1);
    option = setfield(option,'rfxdir',options.rfxdir);
    clear options        
    rfxdir = option.rfxdir;
    if ~exist(rfxdir)
        mkdir(rfxdir)
    end                
    
    % Save output
    if option.doTF && option.doTFbands==0
        if option.doPhase || option.doPhasePower
            outname = [option.rfxdir option.tf_pre 'phz_' option.rsafront option.masknic{mask} option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'];
            save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out
        else
            outname = [option.rfxdir option.tf_pre option.rsafront option.masknic{mask} option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'];
            save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out
        end
    elseif option.doTFbands
        if option.doTWs
            if option.doPhase || option.doPhasePower
                outname = [option.rfxdir option.tf_pre 'bands_phz_' option.rsafront option.masknic{mask} option.midname '.mat'];
                save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out
            else
                outname = [option.rfxdir option.tf_pre 'bands_' option.rsafront option.masknic{mask} option.midname '.mat'];
                save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out
            end
        else
            if option.doPhase || option.doPhasePower
                outname = [option.rfxdir option.tf_pre 'bands_phz_' option.rsafront option.masknic{mask} option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'];
                save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out
            else
                outname = [option.rfxdir option.tf_pre 'bands_' option.rsafront option.masknic{mask} option.midname num2str(option.tw*option.tfstep) 'ms_sTW.mat'];
                save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out
            end
        end
    else
        outname = [option.rfxdir option.rsafront option.masknic{mask} option.midname num2str(option.tw*option.srate) 'ms_sTW.mat'];
        save(outname, 'rsa_out', 'mod_names','option'); clear rsa_out
    end               
    
end % masknic
