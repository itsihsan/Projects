function [tot_err] = cmp2(original,predicted)
% function [tot_err,Count11,Count01,Count00,Count10] = cmp2(original,predicted)

% function [tot_err,r11,r01,r00,r10,Count11,Count01,Count00,Count10] = cmp2(original,predicted)

original=round(original);
predicted=round(predicted);

Count00=0;
Count11=0;
Count10=0;
Count01=0;

total_bad_users  = sum(original(:) == 1);
total_good_users = sum(original(:) == 0);

tt = size(original,1);

for i = 1:1:tt
    if original(i)==0 && predicted(i)==0
    Count00 = Count00 + 1; % no. of both 0, good-good
    end
end

for i = 1:1:tt
    if original(i)==1 && predicted(i)==1
    Count11 = Count11 + 1; % no. of both 1, bad-bad
    end
end

for i = 1:1:tt
    if original(i)==0 && predicted(i)==1
    Count01 = Count01 + 1; % good-bad
    end
end

for i = 1:1:tt
    if original(i)==1 && predicted(i)==0
    Count10 = Count10 + 1; % bad-good
    end
end

r11 = Count11 / total_bad_users; % detection rate - Changes made for Po
r01 = Count01 / total_good_users; %false positive rate - Changes made for Po
r00 = Count00 / total_good_users;
r10 = Count10 / total_bad_users;
% r11 = 2* Count11 / tt; % detection rate
% r01 = 2* Count01 / tt; %false positive rate
% r00 = 2* Count00 / tt;
% r10 = 2* Count10 / tt;
tot_err = (0.5 * r01) + 0.5 * ( 1 - r11 );

end