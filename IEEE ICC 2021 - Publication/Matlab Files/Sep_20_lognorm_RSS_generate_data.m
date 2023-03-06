pathloss_constant = 3;
standard_deviation_shadow = 0; % MAKE SURE TO CHANGE IT

transmit_power = 30; % in dBm
do = 1; % Reference Distance
lembda = (3*10^8)/(9*10^9); % lembda = c/f for 900 MHz band
p_d0 = transmit_power + (10 * pathloss_constant * log10(lembda/(4*pi*do)));

% verifier =  [50,50 ; 50,200 ; 200,50; 200,200 ; 125,125];          % 5 verifiers with users around them
% fp_x = 0:0.5:250;
% fp_y = 0:0.5:250;

% % % Rob's query
% % verifier =  [0,0 ; 0,200 ; 200,0; 200,200];
verifier =  [0,0 ; 0,200 ; 200,0; 200,200 ; 100,100];
fp_x = 0:0.5:200;
fp_y = 0:0.5:200;

N = length(verifier); % Number of verifiers

Data=[];

for x=0:1:length(fp_x)-1
    for y=0:1:length(fp_y)-1
        for n=1:N
            shadow = standard_deviation_shadow * randn(1, 1);
            di(n) = sqrt((fp_x(x+1) - verifier(n,1))^2 + (fp_y(y+1) - verifier(n,2))^2);
            P_di(n) = transmit_power - (p_d0 + (10 * pathloss_constant * log10(di(n)/do)) + shadow);
            shadow = [];
        end
        U(y+1,:) = [fp_x(x+1), fp_y(y+1), P_di];
        di=[];
        P_di=[];
    end
    Data = [Data;U];
    U=[];
end
% RSS_noise = [zeros(length(fp_x)*length(fp_y) , 2), (standard_deviation_shadow * randn((length(fp_x)*length(fp_y)), N))];
% Data_noise = Data + RSS_noise;

% % Data export to Excel - Quarter
% filename = '3N_rawdata_noise_half.xlsx';
% writematrix(Data_noise,filename,'Sheet',1,'Range','A1:E251001')

% % Data export to Excel - Quarter
% filename = '3N_rawdata_noise_quarter.xlsx';
% writematrix(Data_noise,filename,'Sheet',1,'Range','A1:E1002001')