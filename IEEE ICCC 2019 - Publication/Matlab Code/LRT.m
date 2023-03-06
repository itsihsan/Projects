function [U,Y_Good,V,Y_Bad,falsepositive,detection,LRT_value,LRT_value_1] = LRT (d_u, d_v, d_y, d_y_m, R)
%New function [U,Y_Good,V,Y_Bad,falsepositive,detection,LRT_value,LRT_value_1] = LRT (d_u, d_v, d_y, d_y_m, R)

% U_Good = toa_Array_pre(:,1:6);
% Y_Good = toa_Array_pre(:,13:18);
% V_Bad = toa_Array_pre(:,25:30);
% Y_Bad = toa_Array_pre(:,31:36);
% 
% U_Good = [] ;Y_Good = [] ; V_Bad = []; Y_Bad = [];
% [toa_Array] = TOA_new(6,1.2);
U_Good1 = toa_Array_LRT(:,1:6);
Y_Good1 = toa_Array_LRT(:,13:18);
V_Bad1  = toa_Array_LRT(:,26:31);
Y_Bad1  = toa_Array_LRT(:,32:37);
R = (0.01^2)*eye(6)
% U_Good = [] ;Y_Good = [] ; V_Bad = []; Y_Bad = [];
% clear;clc


% Likelihood Calculation
%------------------------------------

% STEP1: Put DiffRSS from bad guy (at optimal location) into two PDFs
%------------------------------------
%Old ones
%pdf_value_m = gaussian_pdf(Y_B,V,R);
%pdf_value_c = gaussian_pdf(Y_B,U,R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% delete the below commands as they are not part of the file.
clear;clc
R = (var_t)*eye(3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdf_value_m = gaussian_pdf(Y_Bad,V_Bad,R);
pdf_value_c = gaussian_pdf(Y_Bad,U_Good,R);

% pdf_value_m = lognormal_pdf(Y_Bad,V_Bad,R);
% pdf_value_c = lognormal_pdf(Y_Bad,U_Good,R);
% LRT Decision Rule: find out the bad guy.
%------------------------------------
for i=1:363
    LRT_value(i,:) = pdf_value_m(i,:) / pdf_value_c(i,:);
        if LRT_value(i,:) >= 1
        detection(i,:) = 1;  % Malicious User is Detected correctly
    else
        detection(i,:) = 0;  % Malicious User is Detected incorrectly
    end
end
% STEP2: Put DiffRSS from good guy (at claimed location) into two PDFs
%------------------------------------

pdf_value_m_1 = gaussian_pdf(Y_Good,V_Bad,R);
pdf_value_c_1 = gaussian_pdf(Y_Good,U_Good,R);

% pdf_value_m_1 = lognormal_pdf(Y_Good,V_Bad,R);
% pdf_value_c_1 = lognormal_pdf(Y_Good,U_Good,R);
% LRT Decision Rule: Confirm a good guy.
%------------------------------------
for j=1:364
    LRT_value_1(j,:) = pdf_value_m_1(j,:)/pdf_value_c_1(j,:);
    if LRT_value_1(j,:) > 1
%         if LRT_value_1(j,:) >= 1 (ziqing old)
        falsepositive(j,:) = 1; % Genuine User is Detected incorrectly
    else
        falsepositive(j,:) = 0; % Genuine User is Detected correctly
    end
end

end
