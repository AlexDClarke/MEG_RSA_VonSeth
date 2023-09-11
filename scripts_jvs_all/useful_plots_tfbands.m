
load('RS_timecourse_RATLsource_60ms_sTW.mat')

band = 4;

%in = [1:10, 12:15];
in = [1 2 4:8 10 12:15];

evc(1,:) = mean(rsa_out.alexnet_convpca_1(in,band,:),1);
evc(2,:) = mean(rsa_out.alexnet_convpca_25(in,band,:),1);
evc(3,:) = mean(rsa_out.alexnet_convpca_67(in,band,:),1);
evc(4,:) = mean(rsa_out.mikenet_pca_ticks24(in,band,:),1);
evc(5,:) = mean(rsa_out.mikenet_pca_ticks59(in,band,:),1);
evc(6,:) = mean(rsa_out.mikenet_pca_ticks1020(in,band,:),1);

figure;
plot(zeros(1,length(evc(1,:))),'LineWidth',1,'Color',[0 0 0])  
hold on
set(gca,'ColorOrder',[0 0 0; 0 0 1;.0 .3 .5;1 0 0;.75 .3 0;.6 .6 0],'NextPlot','replacechildren');
plot(evc','LineWidth',1.8); hold on
plot(zeros(1,length(evc(1,:))),'LineWidth',1,'Color',[0 0 0])  
hold off
set(gca,'XTick',[1:option.jump/option.tfstep:(option.epoch_length-1)/option.tfstep]);
set(gca,'XTickLabel',[-option.baseline_length:option.jump:(option.epoch_length-option.baseline_length-1)]);
set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')
xlim([1 length(option.tfepoch)]);
xlim([13 76]);
%legend('alex1', 'alex25', 'alex67', 'mike24', 'mike59', 'mike1020')
ylim([-5e-3 0.02])
