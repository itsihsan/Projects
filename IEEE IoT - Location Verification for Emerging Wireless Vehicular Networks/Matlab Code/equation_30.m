%  Implementation of equaion (25) for D10

rows = size(Y_Bad,1);
columns = size(Y_Bad,2);
%  rows = 20000;
%  columns = 1;
 
for nn = 1:rows
    for mm = 1:columns
        
%         a_D10 = [];
%         D_10 = [];
%         bin_Y_Bad = [];
%         delta_Y_Bad = [];
        
        a_D10 = sort([0,Y_Bad(nn,mm)]);
        bin_Y_Bad = linspace(a_D10(1),a_D10(2),100000);
        delta_Y_Bad = abs(Y_Bad(nn,mm))/100001;
        
        D_10_V1 =     ([exp(-((V_Bad(nn,mm) - bin_Y_Bad).^2/(2*(var_i + var_t))))]/...
                      [sqrt(2*pi*(var_i + var_t))]) .* ...
                      (erfc([std_i*(V_Bad(nn,mm)- bin_Y_Bad)]/...
                      [std_t*sqrt(2*(var_i + var_t))]))...
                      .*...
                      log10([([exp(-((V_Bad(nn,mm)- bin_Y_Bad).^2/(2*(var_i + var_t))))]/...
                      [sqrt(2*pi*(var_i + var_t))]) .* ...
                      (erfc([std_i*(V_Bad(nn,mm)- bin_Y_Bad)]/...
                      [std_t*sqrt(2*(var_i + var_t))]))] ./ ...
                      H_o(nn,mm));
        
        D_10(nn,mm) = delta_Y_Bad * sum(D_10_V1);
    end
end
D_10_V2 = mean(D_10,1)
D_10 = sum(D_10_V2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Implementation of equaion (26) for D01

% rows = size(Y_Bad,1);
% columns = size(Y_Bad,2);

for n = 1:rows
    for m = 1:columns
        
%         a_D01 = [];
%         D_01 = [];
%         bin_Y_Good = [];
%         delta_Y_Good = [];
        
        a_D01 = sort([0,Y_Good(n,m)]);
        bin_Y_Good = linspace(a_D01(1),a_D01(2),100000);
        delta_Y_Good = abs(Y_Good(n,m))/100001;
        
        D_01_V1 = ([exp(-((U_Good(n,m) - bin_Y_Good).^2/(2*(var_i + var_t))))]/...
                      [sqrt(2*pi*(var_i + var_t))]) .* ...
                      (erfc([std_i*(U_Good(n,m) - bin_Y_Good)]/...
                      [std_t*sqrt(2*(var_i + var_t))]))...
                      .*...
                      log10([([exp(-((U_Good(n,m) - bin_Y_Good).^2/(2*(var_i + var_t))))]/...
                      [sqrt(2*pi*(var_i + var_t))]) .* ...
                      (erfc([std_i*(U_Good(n,m) - bin_Y_Good)]/...
                      [std_t*sqrt(2*(var_i + var_t))]))] ./ ...
                      H_1(n,m));
        
        D_01(n,m) = delta_Y_Good * sum(D_01_V1);
    end
end
D_01_V2 = mean(D_01,1)
D_01 = sum(D_01_V2)

Total_Error_New = 0.5 * (1 - min( sqrt(0.5*D_10),sqrt(0.5*D_01)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The below code is for checking a M X 1 matrix only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U_Good = U_Good(:,1);
% Y_Good = Y_Good(:,1);
% Y_Bad = Y_Bad(:,1);
% V_Bad = V_Bad(:,1);
% H_1 = H_1(:,1);
% H_o = H_o(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for nn=1:length(Y_Bad)
%     a_D10 = sort([0,Y_Bad(nn)]);
%     bin_Y_Bad = linspace(a_D10(1),a_D10(2),100000);
%     delta_Y_Bad = abs(Y_Bad(nn))/100001;
%     
%     D_10_New   = ([exp(-((V_Bad(nn) - bin_Y_Bad).^2/(2*(var_i + var_t))))]/...
%         [sqrt(2*pi*(var_i + var_t))]) .* ...
%         (erfc([std_i*(V_Bad(nn) - bin_Y_Bad)]/...
%         [std_t*sqrt(2*(var_i + var_t))]))...
%         .*...
%         log10([([exp(-((V_Bad(nn) - bin_Y_Bad).^2/(2*(var_i + var_t))))]/...
%         [sqrt(2*pi*(var_i + var_t))]) .* ...
%         (erfc([std_i*(V_Bad(nn) - bin_Y_Bad)]/...
%         [std_t*sqrt(2*(var_i + var_t))]))] ./ ...
%         H_o(nn));
%     
%     D_10_Final(nn) = delta_Y_Bad * sum(D_10_New);
% end
% D_10_Final_mean_2 = mean(D_10_Final)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for nn=1:length(Y_Good)
%     a_D01 = sort([0,Y_Good(nn)]);
%     bin_Y_Good = linspace(a_D01(1),a_D01(2),100000);
%     delta_Y_Good = abs(Y_Good(nn))/100001;
%     
%     D_01_New   = ([exp(-((U_Good(nn) - bin_Y_Good).^2/(2*(var_i + var_t))))]/...
%         [sqrt(2*pi*(var_i + var_t))]) .* ...
%         (erfc([std_i*(U_Good(nn) - bin_Y_Good)]/...
%         [std_t*sqrt(2*(var_i + var_t))]))...
%         .*...
%         log10([([exp(-((U_Good(nn) - bin_Y_Good).^2/(2*(var_i + var_t))))]/...
%         [sqrt(2*pi*(var_i + var_t))]) .* ...
%         (erfc([std_i*(U_Good(nn) - bin_Y_Good)]/...
%         [std_t*sqrt(2*(var_i + var_t))]))] ./ ...
%         H_o(nn));
%     
%     D_01_Final(nn) = delta_Y_Good * sum(D_01_New);
% end
% D_01_Final_mean = mean(D_01_Final)
% 

