function roiRSA_create_MEG_RDMs(option,ss,s)

% Creates MEG-based ROI RDMs from MEG source localised ROI data
%
% Alex (10/2014) Feb 2015

%warning off

% Begin
for sub = 1:length(option.subs)

    if isfield(option,'trials')
    trials = importdata([option.trials option.suborder{sub} '.txt']);
    end   
    
    % Get ordering of selected trials
    for i = 1:length(ss)
        temptrials(i,1) = find(trials == ss(1,i));
    end
    trials = temptrials; clear temptrials
    
    cd([option.datadir option.sub_beg option.subs{sub} option.subdir]);
    
    for mask = 1:length(option.masknic)
        
        sprintf('......Subject %s, Region %s......', option.subs{sub},num2str(mask))
        
        % Read ROI data, initialise RDM matrix
        if option.doTF || option.doTFbands
            infile = [option.tf_pre option.masknic{mask} option.midname option.front option.subs{sub} '.mat'];
            load(infile);
            if isfield(option,'tf_pre2') && (option.doTFbands == 0)
                infile = [option.tf_pre2 option.masknic{mask} option.midname option.front option.subs{sub} '.mat'];
                load(infile);
            end
            if ~isfield(option,'trials')
                trials = 1:size(D,4);
            end                        
            if option.doTF
                ROI_RDMs = single(zeros(length(D(1,1,:,1)),length(D(1,:,1,1)),length(trials)/(length(trials)/24),length(trials)/(length(trials)/24)));
            elseif option.doTFbands
                ROI_RDMs = single(zeros(length(D(1,1,:,1)),length(option.fs),length(trials)/(length(trials)/24),length(trials)/(length(trials)/24)));
            end
        else
            infile = [option.masknic{mask} option.midname option.front option.subs{sub} '.mat'];
            load(infile);
            if ~isfield(option,'trials')
                trials = 1:size(D,3);
            end 
            ROI_RDMs = single(zeros(length(D(1,:,1)),length(trials)/(length(trials)/24),length(trials)/(length(trials)/24)));
        end
        
        % Need for files extracted with tw as one value, and now want to
        % use another without re-extracting
        options = optionsfile(s,1);
        option = setfield(option,'tw',options.tw);
        option = setfield(option,'s',s); clear options
        
        % Create RDMs over time
        %% TF version
        if option.doTF
                        
            option.tw = round((option.tw*option.srate)/option.tfstep);
            
            % Read all data first
            ao = D(:,:,:,trials);
            if isfield(option,'tf_pre2') && (option.doTFbands == 0)
                bo = D2(:,:,:,trials);  
            end
            if option.tw > 1                                
                for time = (option.tw/2+1):length(D(1,1,:,1))-(option.tw/2)                    
                    if option.doavg
                        a = mean(ao(:,:,(time-option.tw/2):(time+option.tw/2),:),3);  % vertices x freq x time x trials
                        if isfield(option,'tf_pre2')
                            b = mean(bo(:,:,(time-option.tw/2):(time+option.tw/2),:),3);  % vertices x freq x time x trials    
                        end
                    else
                        a = ao(:,:,(time-option.tw/2):(time+option.tw/2),:);  % vertices x freq x time x trials
                        if isfield(option,'tf_pre2')
                            b = bo(:,:,(time-option.tw/2):(time+option.tw/2),:);  % vertices x freq x time x trials
                        end
                    end                    
                    a = reshape(permute(a,[2 4 1 3]),length(D(1,:,1,1)),length(D(1,1,1,trials)),[]); % freqs x trials x [time x vertices]
                    if isfield(option,'tf_pre2')
                        b = reshape(permute(b,[2 4 1 3]),length(D(1,:,1,1)),length(D(1,1,1,trials)),[]); % freqs x trials x [time x vertices]    
                    end
                    
                    % RDMs for each freq
                    for f = 1:length(D(1,:,1,1))                                               
                        if isfield(option,'tf_pre2')
                            st_rdm = squareform(pdist([zscore(squeeze(a(f,:,:))) zscore(squeeze(b(f,:,:)))],'correlation'));
                        else
                            st_rdm = squareform(pdist(squeeze(a(f,:,:)),'correlation'));                        
                        end
                        
%                         % set incorrect/rejected as NaNs
%                         reject = abs(correct_trials-1);
%                         st_rdm(find(reject(trials)),:) = NaN;
%                         st_rdm(:,find(reject(trials))) = NaN;
                        
                        % collect data
                        temp_rdms = [];
                        for i = 1:(length(trials)/24) % for the multisession averaging
                            temp_rdms(:,:,i) = st_rdm(24*(i-1)+1:24*i,24*(i-1)+1:24*i);
                        end
                        st_rdm = mean(temp_rdms,3); clear temp_rdms i                        
                        ROI_RDMs(time,f,:,:) = single(st_rdm);
                        clear st_rdm
                    end
                    clear a
                end
            else
                for time = 1:length(D(1,1,:,1))
                    a = ao(:,:,time,:);  % vertices x freqs x time x trials
                    a = reshape(permute(a,[2 4 1 3]),length(D(1,:,1,1)),length(D(1,1,1,trials)),[]); % freqs, trials x [time x vertices]
                    if isfield(option,'tf_pre2')
                        b = bo(:,:,time,:);  % vertices x freqs x time x trials
                        b = reshape(permute(b,[2 4 1 3]),length(D(1,:,1,1)),length(D(1,1,1,trials)),[]); % freqs, trials x [time x vertices]
                    end
                    
                    % RDMs for each freq
                    for f = 1:length(D(1,:,1,1))                                                
                        if isfield(option,'tf_pre2')
                            st_rdm = squareform(pdist([squeeze(a(f,:,:)) squeeze(b(f,:,:))],'correlation'));
                        else
                            st_rdm = squareform(pdist(squeeze(a(f,:,:)),'correlation'));                        
                        end
                        
%                         % set incorrect/rejected as NaNs
%                         reject = abs(correct_trials-1);
%                         st_rdm(find(reject(trials)),:) = NaN;
%                         st_rdm(:,find(reject(trials))) = NaN;
                        
                        % collect data
                        temp_rdms = [];
                        for i = 1:(length(trials)/24) % for the multisession averaging
                            temp_rdms(:,:,i) = st_rdm(24*(i-1)+1:24*i,24*(i-1)+1:24*i);
                        end
                        st_rdm = mean(temp_rdms,3); clear temp_rdms i                        
                        ROI_RDMs(time,f,:,:) = single(st_rdm);
                        clear st_rdm
                    end
                    clear a
                end
            end
            
            %% TF bands version
        elseif option.doTFbands
            
            % Read all data first
            ao = D(:,:,:,trials);
            
            for band = 1:length(option.fs)
                sprintf('..... Frequency band %s .....', num2str(band))
                
                f = option.fs{band};
                if option.tw > 1
                    for time = (option.tw/2+1):length(D(1,1,:,1))-(option.tw/2)
                        if option.doavg
                            a = mean(mean(ao(:,f,(time-option.tw/2):(time+option.tw/2),:),3),2);  % vertices x freq x time x trials
                        else
                            a = mean(ao(:,f,(time-option.tw/2):(time+option.tw/2),:),2);  % vertices x freq x time x trials
                        end
                        a = reshape(permute(a,[4 1 2 3]),length(D(1,1,1,trials)),[]); % trials x [time x vertices x freqs]
                        st_rdm = squareform(pdist(a,'correlation'));
                        
%                         % set incorrect/rejected as NaNs
%                         reject = abs(correct_trials-1);
%                         st_rdm(find(reject(trials)),:) = NaN;
%                         st_rdm(:,find(reject(trials))) = NaN;
                        
                        % collect data
                        temp_rdms = [];
                        for i = 1:(length(trials)/24) % for the multisession averaging
                            temp_rdms(:,:,i) = st_rdm(24*(i-1)+1:24*i,24*(i-1)+1:24*i);
                        end
                        st_rdm = mean(temp_rdms,3); clear temp_rdms i                        
                        ROI_RDMs(time,band,:,:) = single(st_rdm);
                        clear a st_rdm
                    end
                else
                    for time = 1:length(D(1,1,:,1))
                        a = mean(ao(:,f,time,:),2);  % vertices x freqs x time x trials
                        a = reshape(permute(a,[4 1 2 3]),length(D(1,1,1,trials)),[]); % trials x [time x vertices x freqs]
                        st_rdm = squareform(pdist(a,'correlation'));
                        
%                         % set incorrect/rejected as NaNs
%                         reject = abs(correct_trials-1);
%                         st_rdm(find(reject(trials)),:) = NaN;
%                         st_rdm(:,find(reject(trials))) = NaN;
                        
                        % collect data
                        temp_rdms = [];
                        for i = 1:(length(trials)/24) % for the multisession averaging
                            temp_rdms(:,:,i) = st_rdm(24*(i-1)+1:24*i,24*(i-1)+1:24*i);
                        end
                        st_rdm = mean(temp_rdms,3); clear temp_rdms i
                        ROI_RDMs(time,band,:,:) = single(st_rdm);
                        clear a st_rdm
                    end
                end
            end
            %% Do normal
        else
            
            % Read all data first           
            ao = D(:,:,trials);
            
            if option.tw > 1
                
                if mod(option.tw,2)  % must be even number
                    option.tw = option.tw+1;
                end
                
                for time = (option.tw/2+1):length(D(1,:,1))-(option.tw/2)
                    if option.doavg
                        a = mean(ao(:,(time-option.tw/2):(time+option.tw/2),:),2);  % vertices x time x trials
                    else
                        a = ao(:,(time-option.tw/2):(time+option.tw/2),:);  % vertices x time x trials
                    end
                    a = reshape(permute(a,[3 1 2]),length(D(1,1,trials)),[]); % trials x [time x vertices]
                    st_rdm = squareform(pdist(a,'correlation'));
                    
%                     % set incorrect/rejected as NaNs
%                     reject = abs(correct_trials-1);
%                     st_rdm(find(reject(trials)),:) = NaN;
%                     st_rdm(:,find(reject(trials))) = NaN;
                    
                    % collect data
                    temp_rdms = [];
                    for i = 1:(length(trials)/24) % for the multisession averaging
                        temp_rdms(:,:,i) = st_rdm(24*(i-1)+1:24*i,24*(i-1)+1:24*i);
                    end
                    st_rdm = mean(temp_rdms,3); clear temp_rdms i
                    ROI_RDMs(time,:,:) = single(st_rdm);
                    clear a st_rdm
                end
            else
                for time = 1:length(D(1,:,1))
                    a = ao(:,time,:);  % vertices x time x trials
                    a = reshape(permute(a,[3 1 2]),length(D(1,1,trials)),[]); % trials x [time x vertices]
                    st_rdm = squareform(pdist(a,'correlation'));
                    
%                     % set incorrect/rejected as NaNs
%                     reject = abs(correct_trials-1);
%                     st_rdm(find(reject(trials)),:) = NaN;
%                     st_rdm(:,find(reject(trials))) = NaN;
                    
                    % collect data
                    temp_rdms = [];
                    for i = 1:(length(trials)/24) % for the multisession averaging
                        temp_rdms(:,:,i) = st_rdm(24*(i-1)+1:24*i,24*(i-1)+1:24*i);
                    end
                    st_rdm = mean(temp_rdms,3); clear temp_rdms i
                    ROI_RDMs(time,:,:) = single(st_rdm);
                    clear a st_rdm
                end
            end
        end % doTF
        
        clear ao
        
        %% Save data
        if option.doTF
            freqs = option.fsc;
            if option.doavg
                if isfield(option,'tf_pre2')
                    save([option.tf_pre option.tf_pre2 option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'ROI_RDMs','freqs','option');
                else
                    save([option.tf_pre option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'ROI_RDMs','freqs','option');
                end
            else
                if isfield(option,'tf_pre2')
                    save([option.tf_pre option.tf_pre2 option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'ROI_RDMs','freqs','option');
                else
                    save([option.tf_pre option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw*option.tfstep) 'ms_sTW.mat'],'ROI_RDMs','freqs','option');
                end
            end
            
        elseif option.doTFbands
            for i = 1:length(option.fs)
                freqs(i,1) = D.frequencies(option.fs{i}(1));
                freqs(i,2) = D.frequencies(option.fs{i}(end));
            end
            if option.doavg
                save([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw) 'ms_sTW.mat'],'ROI_RDMs','freqs','option');
            else
                save([option.tf_pre 'bands_' option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw) 'ms_sTW.mat'],'ROI_RDMs','freqs','option');
            end
            
        else
            if option.doavg
                save([option.masknic{mask} option.midname 'spatio_RDMs_' num2str(option.tw) 'ms_sTW.mat'],'ROI_RDMs','option');
            else
                save([option.masknic{mask} option.midname 'spatiotemporal_RDMs_' num2str(option.tw) 'ms_sTW.mat'],'ROI_RDMs','option');
            end
        end
        
        clear ROI_RDMs
        
    end % ROI
    
end % sub

clear all
