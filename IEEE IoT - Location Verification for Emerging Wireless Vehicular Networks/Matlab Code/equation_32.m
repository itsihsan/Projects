%********************************************************
% Revised Version for calculation total error 26-Sep-2018
% % This file calculates the Actual Total Error.
%********************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version is based on normalization where the individual probabilities are not normalized first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%********************************************************
% Revised Version for calculation total error 26-Sep-2018
%********************************************************
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
bin_Y = linspace(a_D10(1),a_D10(2),10001);
delta_Y = ( a_D10(2)-a_D10(1) )/10000;
total_bins = length(bin_Y);

for n = 1:rows
    for k = 1:total_bins
        for m = 1:columns
            
            Ho_first_num(k,m)  = exp(-((U(n,m) - bin_Y(k)).^2 ./(2*(var_i + var_t))));
            Ho_first_den  = sqrt(2*pi*(var_i + var_t));
            Ho_second_num(k,m) = std_i*(U(n,m)- bin_Y(k));
            Ho_second_den = std_t*sqrt(2*(var_i + var_t));
%             Ho(k,m) = (Ho_first_num(k,m) ./ Ho_first_den) .* erfc(Ho_second_num(k,m) ./ Ho_second_den);
            Ho(k,m) = delta_Y *((Ho_first_num(k,m) ./ Ho_first_den) .* erfc(Ho_second_num(k,m) ./ Ho_second_den));
            
            H1_first_num(k,m)  = exp(-((V(n,m) - bin_Y(k)).^2 ./(2*(var_i + var_t))));
            H1_first_den  = sqrt(2*pi*(var_i + var_t));
            H1_second_num(k,m) = std_i*(V(n,m)- bin_Y(k));
            H1_second_den = std_t*sqrt(2*(var_i + var_t));
%             H1(k,m) = (H1_first_num(k,m) ./ H1_first_den) .* erfc(H1_second_num(k,m) ./ H1_second_den);
            H1(k,m) = delta_Y *((H1_first_num(k,m) ./ H1_first_den) .* erfc(H1_second_num(k,m) ./ H1_second_den));
            
        end
    end
    
    P_bold = prod(Ho,2); % dimension 2 is product of a row
    Q_bold = prod(H1,2);
    Normalize_P_bold = sum(P_bold);
    Normalize_Q_bold = sum(Q_bold);
    P = P_bold/Normalize_P_bold;
    Q = Q_bold/Normalize_Q_bold;
%     difference = 0.5 * delta_Y * (sum(abs(P - Q)));
    difference = 0.5 * (sum(abs(P - Q)));

    TE(n) = 0.5*(1-difference);
    diff(n) = difference;
    difference = [];
end

TR_diff = mean(diff)
TE_mean = mean(TE)



% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % This version is based on normalization where the individual probabilities
% % % are normalized first, there product is normalized afterwards as well.
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % % U=U_Good(:,1);
% % % V=V_Bad(:,1);
% % % U=U_Good(1,:);
% % % V=V_Bad(1,:);
% % % U=U_Good(1,1);
% % % V=V_Bad(1,1);
% %
% % V=V_Bad; U=U_Good;
% % rows = size(U,1);
% % columns = size(U,2);
% %
% % Y1 = 0;
% % Y2 = 4000;
% % a_D10 = sort([Y1,Y2]);
% % bin_Y = linspace(a_D10(1),a_D10(2),1001);
% % delta_Y = (abs(Y1)+abs(Y2))/1000;
% % total_bins = length(bin_Y);
% %
% % for n = 1:rows
% %     for k = 1:total_bins
% %         for m = 1:columns
% %
% %             Ho_first_num(k,m)  = exp(-((U(n,m) - bin_Y(k)).^2 ./(2*(var_i + var_t))));
% %             Ho_first_den  = sqrt(2*pi*(var_i + var_t));
% %             Ho_second_num(k,m) = std_i*(U(n,m)- bin_Y(k));
% %             Ho_second_den = std_t*sqrt(2*(var_i + var_t));
% %             Ho(k,m) = (Ho_first_num(k,m) ./ Ho_first_den) .* erfc(Ho_second_num(k,m) ./ Ho_second_den);
% %
% %             H1_first_num(k,m)  = exp(-((V(n,m) - bin_Y(k)).^2 ./(2*(var_i + var_t))));
% %             H1_first_den  = sqrt(2*pi*(var_i + var_t));
% %             H1_second_num(k,m) = std_i*(V(n,m)- bin_Y(k));
% %             H1_second_den = std_t*sqrt(2*(var_i + var_t));
% %             H1(k,m) = (H1_first_num(k,m) ./ H1_first_den) .* erfc(H1_second_num(k,m) ./ H1_second_den);
% %
% %         end
% %     end
% %
% %     sum_Ho = sum(Ho,1); % dimension 1 is sum of a column
% %     sum_H1 = sum(H1,1); % dimension 1 is sum of a column
% %     Normalize_Ho = repmat(sum_Ho,[total_bins,1]);
% %     Normalize_H1 = repmat(sum_H1,[total_bins,1]);
% %     P = Ho./Normalize_Ho;
% %     Q = H1 ./Normalize_H1;
% %     P_bold = prod(P,2); % dimension 2 is product of a row
% %     Q_bold = prod(Q,2);
% %     Normalize_P_bold = sum(P_bold,1);
% %     Normalize_Q_bold = sum(Q_bold,1);
% %     Nor_P = P_bold/Normalize_P_bold;
% %     Nor_Q = Q_bold/Normalize_Q_bold;
% %     difference = 0.5*(sum(abs(Nor_P - Nor_Q)));
% %     diff(n) = difference;
% %     TE(n) = 0.5*(1-difference);
% % end
% % Total_Error = mean(TE)
% % TR_diff = mean(diff)
