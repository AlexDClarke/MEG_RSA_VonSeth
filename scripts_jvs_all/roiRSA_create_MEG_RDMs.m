function roiRSA_create_MEG_RDMs(option,ss,s)

% Creates MEG-based ROI RDMs from MEG source localised ROI data
%
% Alex (10/2014) Feb 2015

%warning off

% Begin
for sub = 1:length(option.subs)        

    % Load list of trials to reject
    load([option.datadir option.sub_beg option.subs{sub} '/' option.trials option.subs{sub} '.mat']);
    reject = tmprej; clear tmprej
    
    trials = 1:size(ss,2);
            
    cd([option.datadir option.sub_beg option.subs{sub} option.subdir]);
    
    for mask = 1:length(option.masknic)
        
        sprintf('......Subject %s, Region %s......', option.subs{sub},num2str(mask))
                
        % Need for files extracted with tw as one value, and now want to
        % use another without re-extracting
        options = optionsfile(s,1);
        option = setfield(option,'masknic',options.masknic);
        option = setfield(option,'ROI_coords',options.ROI_coords);
        option = setfield(option,'tw',options.tw); clear options
        
        %% Setup matfile for outputs
        if option.doTF
            option.tw = round((option.tw*option.srate)/option.tfstep);
            if option.tw > 1
                if mod(option.tw,2)  % must be even number
                    option.tw = option.tw+1;
                end
            end
            if option.doavg
                if option.doPhase
                    rdms = matfile([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                else
                    rdms = matfile([option.tf_pre option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                end
            else
                if option.doPhase
                    rdms = matfile([option.tf_pre 'phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                else
                    rdms = matfile([option.tf_pre option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                end
            end
        elseif option.doTFbands
            if option.doTWs
                option.tw = 1;
            else
                option.tw = round((option.tw*option.srate)/option.tfstep);
            end
            if option.tw > 1
                if mod(option.tw,2)  % must be even number
                    option.tw = option.tw+1;
                end
            end
            if option.doavg
                if option.doPhase
                    rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                else
                    rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                end
            else
                if option.doPhase
                    if option.doTWs
                        rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'TW_RDMs.mat'],'Writable',true);
                    else
                        rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                    end
                else
                    if option.doTWs
                        rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'TW_RDMs.mat'],'Writable',true);
                    else
                        rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
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
                rdms = matfile([option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat'],'Writable',true);
            else
                rdms = matfile([option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.srate) 'ms_sTW.mat'],'Writable',true);
            end
        end                       
        
        %% Read ROI data
        if option.doTF || option.doTFbands
            if option.doPhase
                infile = [option.tf_pre 'phz_' option.masknic{mask} option.midname option.front option.subs{sub} '.mat'];
                load(infile);
            else
                infile = [option.tf_pre option.masknic{mask} option.midname option.front option.subs{sub} '.mat'];
                load(infile);
            end
        else
            infile = [option.masknic{mask} option.midname option.front option.subs{sub} '.mat'];
            load(infile);
        end    
        
        % Need for files extracted with tw as one value, and now want to
        % use another without re-extracting
        options = optionsfile(s,0);
        option = setfield(option,'masknic',options.masknic);
        option = setfield(option,'masknic',options.masknic);
        option = setfield(option,'ROI_coords',options.ROI_coords);
        option = setfield(option,'subs',options.subs);
        option = setfield(option,'suborder',options.suborder);
        option = setfield(option,'tw',options.tw); clear options        
        
        % Create RDMs over time
        %% TF version
        if option.doTF || option.doTFbands
            
            if option.doTWs
            option.tw = 1;  % tw needs to be 1 if we'll do TW averaging later
            else
            option.tw = round((option.tw*option.srate)/option.tfstep);
            if option.tw > 1
                if mod(option.tw,2)  % must be even number
                    option.tw = option.tw+1;
                end
            end
            end
            
            % Read all data first
            ao = D(:,:,:,trials);
            at = length(D(1,1,:,1)); % times
            af = length(D(1,:,1,1)); % freqs      
            clear D
            
            if option.tw > 1
                for time = (option.tw/2+1):at-(option.tw/2)
                    if option.doavg
                        a = mean(ao(:,:,(time-option.tw/2):(time+option.tw/2),:),3);  % vertices x freq x time x trials
                    else
                        a = ao(:,:,(time-option.tw/2):(time+option.tw/2),:);  % vertices x freq x time x trials
                    end
                    a = reshape(permute(a,[2 4 1 3]),af,length(trials),[]); % freqs x trials x [time x vertices]
                    
                    % RDMs for each freq
                    for f = 1:af
                        
                        st_rdm = squareform(pdist(squeeze(a(f,:,:)),option.dist1));
                        
                        % set incorrect/rejected as NaNs
                        st_rdm(find(reject(trials)),:) = NaN;
                        st_rdm(:,find(reject(trials))) = NaN;
                        
                        % set outlier trials as NaNs
                        b=[];
                        for i=1:length(trials)
                            if or(~isempty(find(a(f,i,:) > mean2(a(f,:,:))+2*std2(a(f,:,:)))),(~isempty(find(a(f,i,:) < mean2(a(f,:,:))-2*std2(a(f,:,:))))))
                                b(i) = 1;
                            end
                        end
                        st_rdm(find(b),:) = NaN;
                        st_rdm(:,find(b)) = NaN;
                        
                        % collect data
                        temp_o(1,1,:,:) = single(st_rdm);
                        rdms.ROI_RDMs(time,f,1:size(st_rdm,1),1:size(st_rdm,2)) = temp_o;
                        clear st_rdm temp_o
                    end
                    clear a
                end
                rdms.ROI_RDMs(at-(option.tw/2):at,1:af,1:length(rdms.ROI_RDMs(1,1,:,1)),1:length(rdms.ROI_RDMs(1,1,:,1))) = NaN;
            else
                for time = 1:at
                    a = ao(:,:,time,:);  % vertices x freqs x time x trials
                    a = reshape(permute(a,[2 4 1 3]),af,length(trials),[]); % freqs, trials x [time x vertices]
                    
                    % RDMs for each freq
                    for f = 1:af
                        st_rdm = squareform(pdist(squeeze(a(f,:,:)),option.dist1));
                        
                        % set incorrect/rejected as NaNs
                        st_rdm(find(reject(trials)),:) = NaN;
                        st_rdm(:,find(reject(trials))) = NaN;
                        
                        % set outlier trials as NaNs
                        b=[];
                        for i=1:length(trials)
                            if or(~isempty(find(a(f,i,:) > mean2(a(f,:,:))+2*std2(a(f,:,:)))),(~isempty(find(a(f,i,:) < mean2(a(f,:,:))-2*std2(a(f,:,:))))))
                                b(i) = 1;
                            end
                        end
                        st_rdm(find(b),:) = NaN;
                        st_rdm(:,find(b)) = NaN;
                        
                        % collect data
                        temp_o(1,1,:,:) = single(st_rdm);
                        rdms.ROI_RDMs(time,f,1:size(st_rdm,1),1:size(st_rdm,2)) = temp_o;
                        clear st_rdm temp_o
                    end % f
                    clear a
                end % time
            end
            clear ao
            
            %% Time/frequency averaing of similarities
            if option.doTFbands                
                % Fband averaging
                for band = 1:length(option.fs)                    
                    
                    f = option.fs{band};
                    f = find(((option.fsc>=f(1))&(option.fsc<=f(end)))); % index of freqs in band                    
                    
                    if option.doTWs % Time window averaging
                        for tw = 1:length(option.timew)
                            t = find(((option.tfepoch>=option.timew{tw}(1))&(option.tfepoch<=option.timew{tw}(end))));
                            ROI_RDMs_TW(tw,band,1:length(rdms.ROI_RDMs(1,1,:,1)),1:length(rdms.ROI_RDMs(1,1,1,:))) = mean(mean(rdms.ROI_RDMs(t,f,1:length(rdms.ROI_RDMs(1,1,:,1)),1:length(rdms.ROI_RDMs(1,1,1,:))),1),2);
                        end
                    else
                        ROI_RDMs_TW(1:length(rdms.ROI_RDMs(1,1,1,:)),band,1:length(rdms.ROI_RDMs(1,1,:,1)),1:length(rdms.ROI_RDMs(1,1,:,1))) = mean(rdms.ROI_RDMs(1:length(rdms.ROI_RDMs(:,1,1,1)),f,1:length(rdms.ROI_RDMs(1,1,:,1)),1:length(rdms.ROI_RDMs(1,1,1,:))),2);
                    end                    
                end
                
                % Delete old matfile and make new TFbands one (ineligant way to get around replacing data in matfile)
                if option.doavg
                    if option.doPhase
                        delete([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                        rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                    else
                        delete([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                        rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                    end
                else
                    if option.doPhase
                        if option.doTWs
                            delete([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'TW_RDMs.mat']);
                            rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'TW_RDMs.mat'],'Writable',true);
                        else
                            delete([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                            rdms = matfile([option.tf_pre 'bands_phz_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                        end
                    else
                        if option.doTWs
                            delete([option.tf_pre 'bands_' option.masknic{mask} option.midname 'TW_RDMs.mat']);
                            rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'TW_RDMs.mat'],'Writable',true);
                        else
                            delete([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat']);
                            rdms = matfile([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'Writable',true);
                        end
                    end
                end
                
                rdms.ROI_RDMs(1:tw,1:band,1:length(ROI_RDMs_TW(1,1,:,1)),1:length(ROI_RDMs_TW(1,1,1,:))) = ROI_RDMs_TW(1:tw,1:band,1:length(ROI_RDMs_TW(1,1,:,1)),1:length(ROI_RDMs_TW(1,1,1,:)));
                clear ROI_RDMs_TW
            end 

        %% Do normal
        else
            if option.tw > 1
                if mod(option.tw,2)  % must be even number
                    option.tw = option.tw+1;
                end
            end
            
            % Read all data first
            ao = D(:,:,trials);
            at = length(D(1,1,:,1)); % times
            clear D
            
            if option.tw > 1                                
                for time = (option.tw/2+1):at-(option.tw/2)
                    if option.doavg
                        a = mean(ao(:,(time-option.tw/2):(time+option.tw/2),:),2);  % vertices x time x trials
                    else
                        a = ao(:,(time-option.tw/2):(time+option.tw/2),:);  % vertices x time x trials
                    end
                    a = reshape(permute(a,[3 1 2]),length(trials),[]); % trials x [time x vertices]
                    st_rdm = squareform(pdist(a,option.dist1));
                    
                    % set incorrect/rejected as NaNs
                    st_rdm(find(reject(trials)),:) = NaN;
                    st_rdm(:,find(reject(trials))) = NaN;
                    
                    % set outlier trials as NaNs
                    b=[];
                    for i=1:length(trials)
                        if or(~isempty(find(a(i,:) > mean2(a)+2*std2(a))),(~isempty(find(a(i,:) < mean2(a)-2*std2(a)))))
                            b(i) = 1;
                        end
                    end
                    st_rdm(find(b),:) = NaN;
                    st_rdm(:,find(b)) = NaN;
                    
                    % collect data
                    temp_o(1,:,:) = single(st_rdm);
                    rdms.ROI_RDMs(time,1:size(st_rdm,1),1:size(st_rdm,2)) = temp_o;
                    clear a st_rdm temp_o
                end
                rdms.ROI_RDMs(at-(option.tw/2):at,1:length(rdms.ROI_RDMs(1,1,:,1)),1:length(rdms.ROI_RDMs(1,1,:,1))) = NaN;
            else
                for time = 1:at
                    a = ao(:,time,:);  % vertices x time x trials
                    a = reshape(permute(a,[3 1 2]),length(trials),[]); % trials x [time x vertices]
                    st_rdm = squareform(pdist(a,option.dist1));
                    
                    % set incorrect/rejected as NaNs
                    st_rdm(find(reject(trials)),:) = NaN;
                    st_rdm(:,find(reject(trials))) = NaN;
                    
                    % set outlier trials as NaNs
                    b=[];
                    for i=1:length(trials)
                        if or(~isempty(find(a(i,:) > mean2(a)+2*std2(a))),(~isempty(find(a(i,:) < mean2(a)-2*std2(a)))))
                            b(i) = 1;
                        end
                    end
                    st_rdm(find(b),:) = NaN;
                    st_rdm(:,find(b)) = NaN;
                    
                    % collect data
                    temp_o(1,:,:) = single(st_rdm);
                    rdms.ROI_RDMs(time,1:size(st_rdm,1),1:size(st_rdm,2)) = temp_o;
                    clear a st_rdm temp_o
                end
            end
        end % doTF       
        clear ao
        
        %% Add info to saved matfiles
        rdms.option = option;
        if option.doTF
            rdms.freqs = option.fsc;
        elseif option.doTFbands
            rdms.freqs = option.fs;
        end                        
        
    end % ROI    
end % sub
clear all