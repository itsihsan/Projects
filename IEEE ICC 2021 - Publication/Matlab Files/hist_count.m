edges = 1:5:300;
edges(1)=0;
% edges = 1:10:500;
% nn03 = NN_Samples(:,1);
% nn04 = NN_Samples(:,2);
% nn05 = NN_Samples(:,3);
% nn06 = NN_Samples(:,4);
% nn07 = NN_Samples(:,5);
% nn08 = NN_Samples(:,6);
% nn09 = NN_Samples(:,7);
% nn10 = NN_Samples(:,8);

nn03 = NN_Samples(:,1);
nn04 = NN_Samples(:,2);
nn05 = NN_Samples(:,3);
nn06 = NN_Samples(:,4);
nn07 = NN_Samples(:,5);
nn08 = NN_Samples(:,6);
nn09 = NN_Samples(:,7);
nn10 = NN_Samples(:,8);
nn11 = NN_Samples(:,9);
nn12 = NN_Samples(:,10);
nn13 = NN_Samples(:,11);
nn14 = NN_Samples(:,12);
nn15 = NN_Samples(:,13);
nn16 = NN_Samples(:,14);
nn17 = NN_Samples(:,15);
nn18 = NN_Samples(:,16);
nn19 = NN_Samples(:,17);
nn20 = NN_Samples(:,18);
nn30 = NN_Samples(:,19);
nn40 = NN_Samples(:,20);
nn50 = NN_Samples(:,21);


% hist_count_nn02 = histcounts(nn02,edges);
hist_count_nn03 = histcounts(nn03,edges);
hist_count_nn04 = histcounts(nn04,edges);
hist_count_nn05 = histcounts(nn05,edges);
hist_count_nn06 = histcounts(nn06,edges);
hist_count_nn07 = histcounts(nn07,edges);
hist_count_nn08 = histcounts(nn08,edges);
hist_count_nn09 = histcounts(nn09,edges);
hist_count_nn10 = histcounts(nn10,edges);
hist_count_nn11 = histcounts(nn11,edges);
hist_count_nn12 = histcounts(nn12,edges);
hist_count_nn13 = histcounts(nn13,edges);
hist_count_nn14 = histcounts(nn14,edges);
hist_count_nn15 = histcounts(nn15,edges);
hist_count_nn16 = histcounts(nn16,edges);
hist_count_nn17 = histcounts(nn17,edges);
hist_count_nn18 = histcounts(nn18,edges);
hist_count_nn19 = histcounts(nn19,edges);
hist_count_nn20 = histcounts(nn20,edges);
hist_count_nn30 = histcounts(nn30,edges);
hist_count_nn40 = histcounts(nn40,edges);
hist_count_nn50 = histcounts(nn50,edges);
% hist_count_fp = histcounts(fp,edges);

hist_count_nn= [hist_count_nn03; hist_count_nn04; 
hist_count_nn05; hist_count_nn06; hist_count_nn07;
hist_count_nn08; hist_count_nn09; hist_count_nn10; 
hist_count_nn11; hist_count_nn12; hist_count_nn13; hist_count_nn14; 
hist_count_nn15; hist_count_nn16; hist_count_nn17; hist_count_nn18; 
hist_count_nn19; hist_count_nn20; hist_count_nn30; hist_count_nn40; 
hist_count_nn50];
