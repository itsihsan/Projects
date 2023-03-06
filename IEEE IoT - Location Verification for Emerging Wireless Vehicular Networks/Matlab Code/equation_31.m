% biases thermal 400 and NLoS 200 ns

% U = [U_Good;U_Good]; V = [V_Bad;V_Bad]; Y = [Y_Good;Y_Bad];
% clear U_Good; clear V_Bad; clear Y_Bad; clear Y_Good;
% clear;clc
% load 'test.mat'
% rows = size(Y,1);
% columns = size(Y,2);
% Same paramteres for a user in Ho and H1 Version
% U=U_Good(:,1);
% V=V_Bad(:,1);
% U=U_Good(1,:);
% V=V_Bad(1,:);
% U=U_Good(1,1);
% V=V_Bad(1,1);
V=V_Bad; U=U_Good;
rows = size(U,1);
columns = size(U,2);

Y1 = 0;
Y2 = 10000;
a_D10 = sort([Y1,Y2]);
bin_Y = linspace(a_D10(1),a_D10(2),5001);
delta_Y = (abs(Y1)+abs(Y2))/5000;

for n = 1:rows
    for m = 1:columns
        
        Ho_first_num  = exp(-((U(n,m) - bin_Y).^2 ./(2*(var_i + var_t))));
        Ho_first_den  = sqrt(2*pi*(var_i + var_t));
        Ho_second_num = std_i*(U(n,m)- bin_Y);
        Ho_second_den = std_t*sqrt(2*(var_i + var_t));
        Ho = (Ho_first_num ./ Ho_first_den) .* erfc(Ho_second_num ./ Ho_second_den);
        %         Ho(n,:) = Ho_;
        
        H1_first_num  = exp(-((V(n,m) - bin_Y).^2 ./(2*(var_i + var_t))));
        H1_first_den  = sqrt(2*pi*(var_i + var_t));
        H1_second_num = std_i*(V(n,m)- bin_Y);
        H1_second_den = std_t*sqrt(2*(var_i + var_t));
        H1 = (H1_first_num ./ H1_first_den) .* erfc(H1_second_num ./ H1_second_den);
        %         H1(n,:) = H1_;
        
        D_10_V1 = H1 .* log10(H1./Ho);
        D_10(n,m) =    delta_Y * sum(D_10_V1);
        
        D_01_V1 = Ho .* log10(Ho./H1);
        D_01(n,m) = delta_Y * sum(D_01_V1);
    end
end
D_10_V2 = mean(D_10,1)
D_10_Final = sum(D_10_V2)
% D_10_Final = mean(D_10_V2)

D_01_V2 = mean(D_01,1)
D_01_Final = sum(D_01_V2)
% D_01_Final = mean(D_01_V2)

Total_Error_New = 0.5 * (1 - min( sqrt(0.5*D_10_Final),sqrt(0.5*D_01_Final)))


