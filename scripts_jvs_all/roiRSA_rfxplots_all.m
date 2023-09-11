function roiRSA_rfxplots_all(option,s,resetoption)

% Creates timecourse plots for RFX analysis with all models on one plot

% %twsz = {[84/2:114/2 354/2:414/2] [108/2:154/2] [66/2:140/2] [130/2:602/2 606/2:646/2] [258/2:400/2 490/2:538/2 594/2:802/2]};
% twsz = {[88/2:140/2 176/2:218/2 380/2:442/2 576/2:612/2 626/2:670/2]
% [116/2:166/2]
% [74/2:178/2 236/2:388/2 420/2:472/2 536/2:582/2 604/2:648/2]
% [74/2:716/2]
% [272/2:802]};
% twsz = {[]
% []
% []
% [242/2:408/2 500/2:536/2]
% [382/2:416/2 522/2:552/2 600/2:802/2]};
% twsz = {[92/2:154/2 334/2:438/2]
% [102/2:216/2]
% [68/2:400/2]
% [80/2:756/2]
% [146/2:802/2]};
% twsz = {[]
% []
% []
% [168/2:324/2 346/2:418/2]
% [384/2:462/2]};
% twsz = {[92/2:150/2]
% [74/2:200/2]
% [74/2:268/2 288/2:478/2]
% [78/2:754/2]
% [146/2:416/2 554/2:632/2 642/2:802/2]};
% twsz = {[]
% []
% []
% []
% [440/2:496/2 600/2:656/2]};
%twsz = {[] [] [] [] []};
twsz = {[46/2:798/2] [46/2:798/2] [48/2:798/2] [42/2:798/2] [86/2:798/2]};
%twsz = {[] [] [] [] [438/2:482/2]};

%cd(option.rfxdir);
load([option.var_dir option.models '.mat']);  % Load to get names
alpha = option.alpha;
leg={};

if ~or(option.doTF,option.doTFbands)
    %% Standard timecourses
    for mask = 1:length(option.masknic)

        masknam = [option.masknic{mask}];
        if option.tw > 1
            if mod(option.tw,2)  % must be even number
                option.tw = option.tw+1;
            end
        end

        % Load data
        load([option.rsafront masknam option.midname num2str(option.tw*option.srate) 'ms_sTW.mat'],'rsa_out','option');

        mod_names = fieldnames(Models);
        nmods = size(mod_names,1);

        if resetoption == 1
            option = optionsfile(s);
        end

        figure; hold on
        map = colormap('lines');
        map = map([1 6 5 3 4],:) .* 0.95;
%       map = brewermap(5,'*YlGnBu');
        for modl = 1:nmods

            % Calculate mean, st, t-stats
            x = rsa_out.(mod_names{modl});
            x(isnan(x)) = 0;
            x(option.sub_out,:) = [];
            
            x = smoothdata(x,2,'gaussian',12);

            av = mean(x);
            rsa_std = std(x,[],1);

            % Spearmans plot
%           Color = [1 0 0];
            Color = map(modl,:);
            if modl == 1
                plot(zeros(1,length(av)),'LineWidth',1,'Color',[0 0 0])
            end
            se = rsa_std(1,:)/sqrt(size(x,1));
            h1 = boundedline([1:length(av)],av,se,'alpha','cmap',Color);
            set(h1,'LineWidth',3);   
            
            sig = nan(1,length(av));
            tws = twsz{modl};
            if ~isempty(tws)
%               sig(tws+500) = -0.005 - (modl*0.0005);
%               sig(tws+500) = -0.0015 - (modl*0.00015);
%               sig(tws+500) = -0.0004 - (modl*0.00005);
                sig(tws+500) = -0.01 - (modl*0.002);
                plot(sig,'LineWidth',5,'Color',map(modl,:));
            end
            
        end % modl

        %leg = {'','',mod_names{1},'',mod_names{2},'',mod_names{3},'',mod_names{4},'',mod_names{5},''};
        %legend(leg);
        %legend('boxoff')
        set(gca,'XTick',[1:option.jump/option.srate:option.epoch_length/option.srate]);
        set(gca,'XTickLabel',[round(option.start):option.jump:round(option.stop)]);
        set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')
%       ylim([-0.01 0.015])
%       ylim([-0.003 0.004])
%       ylim([-0.003 0.006])
%       ylim([-0.00075 0.0015])
        ylim([-0.03 0.1])
        set(get(gca,'YLabel'),'String','Similarity', 'fontweight','b')
        xlim(option.range)
        hold off
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff',[masknam 'all_mods_correlation']);
        close

    end % mask
end

if option.doTF && option.doTFbands==0
end

if option.doTFbands
end
