% NN_Shihao was obtained when rng(0) was set for the first time and then it
% was disabled in the 2nd run. the threshold was set to 'if tot_err(i) <= 0.2470'

rng(0);
load('25-Mar-NN-data-100m.mat','data','original','Sample')
% load('data1019_d.mat','data','original','Sample')
% load('data1015_d.mat','data','original','Sample')
data = data(randperm(end),:);
data = data(randperm(end),:);
data = data(randperm(end),:);
data = data(randperm(end),:);
data = data(randperm(end),:);
data = data(randperm(end),:);

data = data(randperm(end),:);
data = data(randperm(end),:);
data = data(randperm(end),:);
data = data(randperm(end),:);
data = data(randperm(end),:);
data = data(randperm(end),:);
data = data(randperm(end),:);
data = data(randperm(end),:);

% data=[];Sample=[];original=[];

inp = [];
targ = [];
inputs = [];
targets = [];
inp = data(:, 1:6);   % Change the columns number as required
targ = data(:,7);     % Change the columns number as required
net = feedforwardnet; %([10 10]);
% net = cascadeforwardnet;
% net = feedforwardnet([10 10]);
net.trainFcn = 'trainlm';
net.layers{1}.transferFcn = 'tansig'; % default is tansig
% net.layers{1}= leakyReluLayer(0.01);

rows = size(inp,1);
tot_err=[];
count=1;

for i=3:1:rows-1
    inputs  = inp(1:i+1,:);
    targets = targ(1:i+1,:);
    [net,parameters] = train(net,inputs',targets'); % 'parameters' contains all of the information concerning the training of the network.
    Output = net(Sample')';
    predicted = Output;
    predicted = round(predicted);
    predicted(predicted<0) = 0;
    predicted(predicted>1) = 1;
    tot_err(1)=0.5;
    tot_err(2)=0.5;
%     tot_err(3)=0.5;
    tot_err(i) = cmp2(original,predicted);
    plot(tot_err,'r')
%     ylim([0 0.6]);
    hold on
    total_err_record(1)=0.5;
    total_err_record(2)=0.5;
%     total_err_record(3)=0.5;
    tot_err(i)
    total_err_record(i)=tot_err(i);
%     if i>15 && tot_err(i) <= 0.001 % NN_Shihao Setting
%         %     if tot_err(i) <= 0.1          % NN_Ihsan_Setting
%         disp ('found!!!')
%         break
%     end
%     net_last=net;
end

% for i=4:2:rows-1
%     inputs  = inp(1:i,:);
%     targets = targ(1:i,:);
%     [net,parameters] = train(net,inputs',targets'); % 'parameters' contains all of the information concerning the training of the network.
%     Output = net(Sample')';
%     predicted = Output;
%     predicted= round(predicted);
%     predicted(predicted<0) = 0;
%     predicted(predicted>1) = 1;
%     
%     if mod(i,2)==0
%         tot_err(count) = cmp2(original,predicted);
%     end
%     plot(tot_err,'k')
%     hold on
% 
%     total_err_record(count)=tot_err(count);
%     total_err_record(end)
%     if count>15 && tot_err(count) < 0.0005
%         disp ('found!!!')
%         break
%     end
%     count=count+1;
% end