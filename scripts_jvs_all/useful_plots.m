
load('RS_timecourse_RATLsource_60ms_sTW.mat')

in = [2 4 5 7:15];

evc(1,:) = mean(rsa_out.alexnet_convpca_1(in,:),1);
evc(2,:) = mean(rsa_out.alexnet_convpca_25(in,:),1);
evc(3,:) = mean(rsa_out.alexnet_convpca_67(in,:),1);
evc(4,:) = mean(rsa_out.mikenet_pca_ticks26(in,:),1);
evc(5,:) = mean(rsa_out.mikenet_pca_ticks720(in,:),1);


figure;
set(gca,'ColorOrder',[0 0 0; 0 0 1;.0 .3 .5;1 0 0;.75 .3 0;.6 .6 0],'NextPlot','replacechildren');
plot(evc','LineWidth',2);
set(gca,'XTick',[1:option.jump/option.srate:option.epoch_length/option.srate]);
set(gca,'XTickLabel',[round(option.start):option.jump:round(option.stop)]);
set(get(gca,'XLabel'),'String','Time (ms)', 'fontweight','b')
set(get(gca,'YLabel'),'String','Similarity', 'fontweight','b')
xlim([501 1150])
%hold on
%plot(nc_out.lower,'--');
%hold off
legend('alex1', 'alex25', 'alex67', 'mike26', 'mike720')
