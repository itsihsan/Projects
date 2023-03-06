% function [verifier, claimed_location , attack_location, U, Delta_U, R, D] = Optimum_DRSS(distant_constraint)
% function [data] = fingerprint_temp_nn_data(samples)

speedOfLight = 3e8;
verifier =  [100,100 ; 100,400 ; 400,100; 400,400 ; 250,250];          % 5 verifiers with users around them
% verifier =  [70,20 ; 100,180 ; 200,50; 125,120 ; 250,250];          % 5 verifiers with users around them

N = length(verifier); % Number of verifiers
fp_x = 0:500/1000:500;
fp_y = 0:500/1000:500;
% fp_x = 0:1:250;
% fp_y = 0:1:250;
% fp_x = 0:3/3:3;
% fp_y = 0:3/3:3;
% standard_deviation_ToA_noise = 1.6667e-08; % equivalent to 5 meters.
standard_deviation_ToA_noise = 3.3333e-08; % equivalent to 10 meters.
Data=[];

for x=0:1:length(fp_x)-1
    for y=0:1:length(fp_y)-1
        for n=1:N
            d2t(n) = sqrt((fp_x(x+1) - verifier(n,1))^2 + (fp_y(y+1) - verifier(n,2))^2);
            ToA(n) = d2t(n)/speedOfLight;
        end
%         RSS_noise = noise(standard_deviation_rss_noise, N)
%         RSS = RSS + RSS_noise
        U(y+1,:) = [fp_x(x+1), fp_y(y+1), ToA];
        d2t=[];
        ToA=[];
    end
    Data = [Data;U];
    U=[];
end
ToA_noise = [zeros(length(fp_x)*length(fp_y) , 2), (standard_deviation_ToA_noise * randn((length(fp_x)*length(fp_y)), N))];
Data_noise = Data + ToA_noise;

% Data export to Excel
filename = 'ToA_data_noise.xlsx';
writematrix(Data_noise,filename,'Sheet',1,'Range','A1:G1002001')