 tt=1
 tick = 0
 false_detection_record=[];
 detection_record=[];
  for i = 1 : 2000

            [verifier_pre, verifier, claimed_location , optimal_attack_location, U, Delta_U, R, D] = Optimum_DRSS(100);

            [falsepositive, detection] = DRSSLV (verifier_pre,verifier,claimed_location,optimal_attack_location,R,D);

            false_detection_record = [false_detection_record , falsepositive];
            detection_record = [detection_record , detection];
            tick = tick + 1 % iteration indicator
        end 

    False_Positive_Rate = numel(find(false_detection_record == 1))/numel(false_detection_record);
    Detection_Rate = numel(find(detection_record == 1))/numel(detection_record);
    a (tt) = False_Positive_Rate;
    b (tt) = Detection_Rate;
    c (tt) = 0.5 * a (tt) + 0.5 * (1 - b (tt));
end