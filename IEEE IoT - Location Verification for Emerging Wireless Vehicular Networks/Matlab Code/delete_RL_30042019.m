% This file can be deleted as it is just to serve the temporary purpose of
% calculating revised LRT with the RL Arpil 30 2019
R = (var_t)*eye(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdf_value_m = gaussian_pdf(Y_Bad,V_Bad,R);
pdf_value_c = gaussian_pdf(Y_Bad,U_Good,R);

% LRT Decision Rule: find out the bad guy.
%------------------------------------
for i=1:5000
    LRT_value(i,:) = pdf_value_m(i,:) / pdf_value_c(i,:); %H1/H0
    if LRT_value(i,:) >= 0.25
        detection(i,:) = 1;  % Malicious User is Detected correctly
    else
        detection(i,:) = 0;  % Malicious User is Detected incorrectly
    end
end
% STEP2: Put DiffRSS from good guy (at claimed location) into two PDFs
%------------------------------------
%Old ones
% pdf_value_m_1 = gaussian_pdf(Y_G,V,R);
% pdf_value_c_1 = gaussian_pdf(Y_G,U,R);

pdf_value_m_1 = gaussian_pdf(Y_Good,V_Bad,R);
pdf_value_c_1 = gaussian_pdf(Y_Good,U_Good,R);
% LRT Decision Rule: Confirm a good guy.
%------------------------------------
for j=1:5000
    LRT_value_1(j,:) = pdf_value_m_1(j,:)/pdf_value_c_1(j,:); %H1/H0
    if LRT_value_1(j,:) > 0.25
%         if LRT_value_1(j,:) >= 1 (ziqing old)
        falsepositive(j,:) = 1; % Genuine User is Detected incorrectly
    else
        falsepositive(j,:) = 0; % Genuine User is Detected correctly
    end
end

