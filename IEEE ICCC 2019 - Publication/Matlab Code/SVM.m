rng(0)
data = data(randperm(end),:);
data = data(randperm(end),:);
data = data(randperm(end),:);
data = data(randperm(end),:);

total_examples = size(data,1);

X = data(1:end,1:end-1);
y = data(1:end,end);

test_threshold = ceil(0.10*total_examples);
rand_num = randperm(total_examples);
X_train = X(rand_num(1:total_examples-test_threshold),:);
y_train = y(rand_num(1:total_examples-test_threshold),:);
X_test = X(rand_num(end-test_threshold+1:end),:);
y_test = y(rand_num(end-test_threshold+1:end),:);

c = cvpartition(y_train,'k',5);

opts = statset('display','iter');
fun = @(train_data,train_labels,test_data,test_labels)...
    sum(predict(fitcsvm(train_data,train_labels,'KernelFunction','rbf'),test_data) ~= test_labels);

[fs,history] = sequentialfs(fun,X_train,y_train,'cv',c,'options',opts,'nfeatures',2,'direction','backward');

X_train_w_best_features = X_train(:,fs);
Mdl = fitcsvm(X_train_w_best_features,y_train,'KernelFunction','rbf',...
    'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus','Showplots',true));

X_test_w_best_features = X_test(:,fs);
accuracy = sum(predict(Mdl,X_test_w_best_features) == y_test)/length(y_test) *100













% rng(0)
% data = data(randperm(end),:);
% data = data(randperm(end),:);
% data = data(randperm(end),:);
% data = data(randperm(end),:);
% 
% 
% 
% X = data(1:5000,1:end-1);
% y = data(1:5000,end);
% 
% rand_num = randperm(5000);
% X_train = X(rand_num(1:4000),:);
% y_train = y(rand_num(1:4000),:);
% X_test = X(rand_num(4001:5000),:);
% y_test = y(rand_num(4001:5000),:);
% 
% c = cvpartition(y_train,'k',5);
% 
% opts = statset('display','iter');
% fun = @(train_data,train_labels,test_data,test_labels)...
%     sum(predict(fitcsvm(train_data,train_labels,'KernelFunction','rbf'),test_data) ~= test_labels);
% 
% [fs,history] = sequentialfs(fun,X_train,y_train,'cv',c,'options',opts,'nfeatures',12,'direction','backward');
% 
% X_train_w_best_features = X_train(:,fs);
% Mdl = fitcsvm(X_train_w_best_features,y_train,'KernelFunction','rbf',...
%     'OptimizeHyperparameters','auto',...
%     'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
%     'expected-improvement-plus','Showplots',true));
% 
% X_test_w_best_features = X_test(:,fs);
% accuracy = sum(predict(Mdl,X_test_w_best_features) == y_test)/length(y_test) *100



% load fisheriris.mat
% species_num = grp2idx(species);
% 
% rng(0)
% data = data(randperm(end),:);
% data = data(randperm(end),:);
% data = data(randperm(end),:);
% data = data(randperm(end),:);
% 
% 
% 
% X = randn(100,10);
% X = randn(5000,30);
% 
% X(:,[1,3,5,7]) = meas(1:100,:);
% X(:,[1,3,5,7,9,11,13,15,17,19,21,23]) = data(1:5000,1:end-1);
% X = data(1:5000,1:end-1);
% 
% y = species_num(1:100,:);
% y = data(1:5000,end);
% 
% rand_num = randperm(100);
% X_train = X(rand_num(1:80),:);
% y_train = y(rand_num(1:80),:);
% rand_num = randperm(5000);
% X_train = X(rand_num(1:4000),:);
% y_train = y(rand_num(1:4000),:);
% 
% 
% X_test = X(rand_num(81:end),:);
% y_test = y(rand_num(81:end),:);
% X_test = X(rand_num(4001:end),:);
% y_test = y(rand_num(4001:end),:);
% 
% c = cvpartition(y_train,'k',5);
% 
% opts = statset('display','iter');
% fun = @(train_data,train_labels,test_data,test_labels)...
%     sum(predict(fitcsvm(train_data,train_labels,'KernelFunction','rbf'),test_data) ~= test_labels);
% 
% [fs,history] = sequentialfs(fun,X_train,y_train,'cv',c,'options',opts,'nfeatures',2);
% [fs,history] = sequentialfs(fun,X_train,y_train,'cv',c,'options',opts,'nfeatures',12);
% 
% X_train_w_best_features = X_train(:,fs);
% Mdl = fitcsvm(X_train_w_best_features,y_train,'KernelFunction','rbf',...
%     'OptimizeHyperparameters','auto',...
%     'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
%     'expected-improvement-plus','Showplots',true));
% 
% X_test_w_best_features = X_test(:,fs);
% accuracy = sum(predict(Mdl,X_test_w_best_features) == y_test)/length(y_test) *100;
% 
% % hyperplane ??
% 
% figure;
% hgscatter = gscatter(X_train_w_best_features(:,1),X_train_w_best_features(:,2),y_train);
% hold on;
% h_sv=plot(Mdl.SupportVectors(:,1),Mdl.SupportVectors(:,2),'ko','markersize',8);
% 
% 
% % test set? data? ?? ??? ????.
% 
% gscatter(X_test_w_best_features(:,1),X_test_w_best_features(:,2),y_test,'rb','xx')
% 
% % decision plane
% XLIMs = get(gca,'xlim');
% YLIMs = get(gca,'ylim');
% [xi,yi] = meshgrid([XLIMs(1):0.01:XLIMs(2)],[YLIMs(1):0.01:YLIMs(2)]);
% dd = [xi(:), yi(:)];
% pred_mesh = predict(Mdl, dd);
% redcolor = [1, 0.8, 0.8];
% bluecolor = [0.8, 0.8, 1];
% pos = find(pred_mesh == 1);
% h1 = plot(dd(pos,1), dd(pos,2),'s','color',redcolor,'Markersize',5,'MarkerEdgeColor',redcolor,'MarkerFaceColor',redcolor);
% pos = find(pred_mesh == 2);
% h2 = plot(dd(pos,1), dd(pos,2),'s','color',bluecolor,'Markersize',5,'MarkerEdgeColor',bluecolor,'MarkerFaceColor',bluecolor);
% uistack(h1,'bottom');
% uistack(h2,'bottom');
% legend([hgscatter;h_sv],{'setosa','versicolor','support vectors'})