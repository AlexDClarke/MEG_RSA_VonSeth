
%% Create models based on PCA/SVD of model activations

part = {'conv1 2' 'conv2 2' 'conv3 2' 'conv4 2' 'conv5 2' 'fc6 2' 'fc7 2'};
for p = 1:length(part)
    M = single(load(['concept_X_node_activations_' part{p} '.txt']));
    [~,N,~,~,ex] = pca(M);
    N = N .* repmat(ex',[size(N,1) 1]);
    alexnet_pca(p,:,:) = N;
    
    MM = squareform(pdist(N,'correlation'));
    fullrdm_conv_5_2 = MM(logical(are_alexnet627_in_302set),logical(are_alexnet627_in_302set));
    aaa(alexnet_order_to_c302_order,alexnet_order_to_c302_order) = fullrdm_conv_5_2;
    aaa([221,266],:) = NaN;
    aaa(:,[221,266]) = NaN;
    aab = aaa(c302_to_catorder,c302_to_catorder);
    Models.(['alexnet_convpca_' num2str(p)]) = aab;
    clear N aaa aab
    
end

%% combine layers after PCA
%load('alexnet_layer_activation_pca_627items.mat')
alexnet_pca = permute(alexnet_pca,[2 3 1]);

l25 = alexnet_pca(:,:,2:5);
l25 = reshape(l25,627,[]);
MM = squareform(pdist(l25,'correlation'));
fullrdm_conv_5_2 = MM(logical(are_alexnet627_in_302set),logical(are_alexnet627_in_302set));
aaa(alexnet_order_to_c302_order,alexnet_order_to_c302_order) = fullrdm_conv_5_2;
aaa([221,266],:) = NaN;
aaa(:,[221,266]) = NaN;
aab = aaa(c302_to_catorder,c302_to_catorder);
Models.alexnet_convpca_25 = aab;

l25 = alexnet_pca(:,:,6:7);
l25 = reshape(l25,627,[]);
MM = squareform(pdist(l25,'correlation'));
fullrdm_conv_5_2 = MM(logical(are_alexnet627_in_302set),logical(are_alexnet627_in_302set));
aaa(alexnet_order_to_c302_order,alexnet_order_to_c302_order) = fullrdm_conv_5_2;
aaa([221,266],:) = NaN;
aaa(:,[221,266]) = NaN;
aab = aaa(c302_to_catorder,c302_to_catorder);
Models.alexnet_convpca_67 = aab;

l25 = alexnet_pca(:,:,:);
l25 = reshape(l25,627,[]);
MM = squareform(pdist(l25,'correlation'));
fullrdm_conv_5_2 = MM(logical(are_alexnet627_in_302set),logical(are_alexnet627_in_302set));
aaa(alexnet_order_to_c302_order,alexnet_order_to_c302_order) = fullrdm_conv_5_2;
aaa([221,266],:) = NaN;
aaa(:,[221,266]) = NaN;
aab = aaa(c302_to_catorder,c302_to_catorder);
Models.alexnet_convpca_all = aab;

% c(:,1) = vectorizeRDM(Models.alexnet_convpca_1);
% c(:,2) = vectorizeRDM(Models.alexnet_convpca_2);
% c(:,3) = vectorizeRDM(Models.alexnet_convpca_3);
% c(:,4) = vectorizeRDM(Models.alexnet_convpca_4);
% c(:,5) = vectorizeRDM(Models.alexnet_convpca_5);
% c(:,6) = vectorizeRDM(Models.alexnet_convpca_6);
% c(:,7) = vectorizeRDM(Models.alexnet_convpca_7);

Models2.alexnet_convpca_1 = Models.alexnet_convpca_1(logical(catorder_in_reducedset),logical(catorder_in_reducedset));
Models2.alexnet_convpca_2 = Models.alexnet_convpca_2(logical(catorder_in_reducedset),logical(catorder_in_reducedset));
Models2.alexnet_convpca_3 = Models.alexnet_convpca_3(logical(catorder_in_reducedset),logical(catorder_in_reducedset));
Models2.alexnet_convpca_4 = Models.alexnet_convpca_4(logical(catorder_in_reducedset),logical(catorder_in_reducedset));
Models2.alexnet_convpca_5 = Models.alexnet_convpca_5(logical(catorder_in_reducedset),logical(catorder_in_reducedset));
Models2.alexnet_convpca_6 = Models.alexnet_convpca_6(logical(catorder_in_reducedset),logical(catorder_in_reducedset));
Models2.alexnet_convpca_7 = Models.alexnet_convpca_7(logical(catorder_in_reducedset),logical(catorder_in_reducedset));
Models2.alexnet_convpca_25= Models.alexnet_convpca_25(logical(catorder_in_reducedset),logical(catorder_in_reducedset));
Models2.alexnet_convpca_67= Models.alexnet_convpca_67(logical(catorder_in_reducedset),logical(catorder_in_reducedset));
Models2.alexnet_convpca_all=Models.alexnet_convpca_all(logical(catorder_in_reducedset),logical(catorder_in_reducedset));
