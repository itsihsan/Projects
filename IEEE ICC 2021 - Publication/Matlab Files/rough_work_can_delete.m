% Test_noise=readmatrix('3N_Test_half.csv','Range','G1:K15000');
% test_noise_rss = Test_noise(:,1:3);
% test_noise_coordinates = Test_noise(:,4:5);
% clear Test_noise
% Train_noise=readmatrix('3N_Train_half.csv','Range','G1:K236001');
% train_noise_rss = Train_noise(:,1:3);
% train_noise_coordinates = Train_noise(:,4:5);
% clear Train_noise

% Test_noise=readmatrix('5N_Test_half_lognormal.csv','Range','K1:Q15000');
% test_noise_rss = Test_noise(:,1:5);
% test_noise_coordinates = Test_noise(:,6:7);
% clear Test_noise
% Train_noise=readmatrix('5N_Train_half_lognormal.csv','Range','K1:Q236001');
% train_noise_rss = Train_noise(:,1:5);
% train_noise_coordinates = Train_noise(:,6:7);
% clear Train_noise

for i=1:length(test_noise_rss)

%     fprintf('\n Example number : ')
    disp(i)
    for j=1:length(train_noise_rss)
        mse(j,1) = norm(train_noise_rss(j,:) - test_noise_rss(i,:));
%         diff_man(j,1) = sqrt( (train_noise_rss(j,1)-test_noise_rss(i,1))^2 + (train_noise_rss(j,2)-test_noise_rss(i,2))^2 + (train_noise_rss(j,3)-test_noise_rss(i,3))^2 + (train_noise_rss(j,4)-test_noise_rss(i,4))^2 + (train_noise_rss(j,5)-test_noise_rss(i,5))^2 );
    end
    [value(1,i), indice(1,i)] = min(mse);
    estimated_position(i,:) = train_noise_coordinates(indice(1,i),:);
    error(i) = norm(test_noise_coordinates(i,:)-estimated_position(i,:));
%     estimated_position(end,:)
%     error(end)
    fprintf('The x-coordinates is = %.1f , y-coordinate is = %.1f, and mse is = %.4f',estimated_position(end,1),estimated_position(end,2),error(end))
    j=[];
    mse=[];
%     diff_man=[];
end


