% function [verifier, claimed_location , attack_location, U, Delta_U, R, D] = Optimum_DRSS(distant_constraint)
% function [data] = fingerprint_temp_nn_data(samples)

pathloss_constant = 3;
standard_deviation_rss_noise = 0; % MAKE SURE TO CHANGE THIS PARAMETER

transmit_power = 30; % in dBm
do = 1; % Reference Distance
lembda = (3*10^8)/(9*10^9); % lembda = c/f for 900 MHz band
p = transmit_power + 20 * log10(lembda/(4*pi*do));

% verifier =  [100,100 ; 100,400 ; 400,100; 400,400 ; 250,250];          % 5 verifiers with users around them
verifier =  [50,50 ; 50,200 ; 200,50; 200,200 ; 125,125];          % 5 verifiers with users around them
% verifier =  [50,50 ; 200,50; 200,200 ; 125,125]; % 4 verifiers with users around them
% verifier =  [50,50 ; 200,200 ; 125,125]; % 3 verifiers with users around them
% verifier =  [0,0 ; 0,200 ; 200,0; 200,200]; % Rob request for extra figures

N = length(verifier); % Number of verifiers
fp_x = 40:40:160;
fp_y = 40:40:160;


Data=[];

for x=0:1:length(fp_x)-1
    for y=0:1:length(fp_y)-1
        for n=1:N
            d2t(n) = sqrt((fp_x(x+1) - verifier(n,1))^2 + (fp_y(y+1) - verifier(n,2))^2);
            RSS(n) = p - 10 * pathloss_constant * log10(d2t(n)/do);
        end
        U(y+1,:) = [fp_x(x+1), fp_y(y+1), RSS];
        d2t=[];
        RSS=[];
    end
    Data = [Data;U];
    U=[];
end
RSS_noise = [zeros(length(fp_x)*length(fp_y) , 2), (standard_deviation_rss_noise * randn((length(fp_x)*length(fp_y)), N))];
Data_noise = Data + RSS_noise;

% % Data export to Excel - Quarter
% filename = '3N_rawdata_noise_half.xlsx';
% writematrix(Data_noise,filename,'Sheet',1,'Range','A1:E251001')

% % Data export to Excel - Quarter
% filename = '3N_rawdata_noise_quarter.xlsx';
% writematrix(Data_noise,filename,'Sheet',1,'Range','A1:E1002001')