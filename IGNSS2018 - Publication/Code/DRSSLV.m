function [falsepositive, detection,LRT_value,LRT_value_1] = DRSSLV (verifier_pre,verifier,claimed_location,optimal_attack_location,R,D)
% DRSSLV - Description
%
% Syntax: [falsepositive, detection] = DRSSLV(verifier_pre,verifier,claimed_location,optimal_attack_location,R,D)
% 
%------------------------------------
% Observation Model
N = length(verifier_pre);
% Calculate RSS a[1,1].-[1;1]t N different verifiers (From the optimal attacking position)
mean_rss_from_optimal_r = channel_model (optimal_attack_location,verifier_pre);
rss_from_optimal_r = add_shadowing (mean_rss_from_optimal_r,R);

% Calculate DRSS at 6 different verifiers (From the optimal attacking position)
for a = 1:N-1
        d_y_m(a) = rss_from_optimal_r(a) - rss_from_optimal_r(N);
end

% Calculate RSS at 6 different verifiers (From the claimed position, assuming there's a good guy there)
mean_rss_from_claim_r = channel_model(claimed_location,verifier_pre);
rss_from_claimed_r = add_shadowing (mean_rss_from_claim_r,R);

% Calculate DRSS at 6 different verifiers (From the claimed position, assuming there's a good guy there)
for b = 1:N-1
        d_y(b) = rss_from_claimed_r(b) - rss_from_claimed_r(N);
end


% Calculate some mean DRSS vectors which may be used in Location Verification
mean_rss_from_claim_t = channel_model(claimed_location,verifier);
mean_rss_from_optimal_t = channel_model(optimal_attack_location,verifier);

for m = 1:N-1
        d_u(m) = mean_rss_from_claim_t(m) - mean_rss_from_claim_t(N);
end

for n = 1:N-1
        d_v(n) = mean_rss_from_optimal_t(n) - mean_rss_from_optimal_t(N);
end

% Likelihood Calculation
%------------------------------------

% STEP1: Put DiffRSS from bad guy (at optimal location) into two PDFs
%------------------------------------
pdf_value_m = mv_gaussian_pdf(d_y_m,d_v,D);   
pdf_value_c = mv_gaussian_pdf(d_y_m,d_u,D);

% LRT Decision Rule: find out the bad guy.
%------------------------------------
LRT_value = pdf_value_m / pdf_value_c;
if LRT_value >= 1          
    detection = 1;
else
    detection = 0;
end

% STEP2: Put DiffRSS from good guy (at claimed location) into two PDFs
%------------------------------------
pdf_value_m_1 = mv_gaussian_pdf(d_y,d_v,D);
pdf_value_c_1 = mv_gaussian_pdf(d_y,d_u,D);

% LRT Decision Rule: Confirm a good guy.
%------------------------------------
LRT_value_1 = pdf_value_c_1 / pdf_value_m_1;
if LRT_value_1 <= 1
    falsepositive = 1;
else
    falsepositive = 0;
end

%d_y_m
%d_y
%d_u
%d_v
%mean_rss_from_optimal_t

end