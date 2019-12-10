%% how many genes change expression and do/don't change location
DD = '~/Develop/DiGiovanni_DiStefano_FC/Data/';
M = readtable( [ DD 'FoldChangeExpression__GenesThatMove_gt5.txt']);
T = readtable( [ DD 'FoldChangeExpression__AllGenes.txt']);
S = T( ~ismember(T.ORF,M.ORF),:);

%%


fprintf('%d\tgenes don''t move >= 5%%.\n' , height(S) ); 
fprintf('%d\tgenes move >= 5%%.\n' , height(M) );

p = 0.001 ; 

fprintf('p-value threshold = %0.03f\n' , p )

fprintf('%d\tgenes change expression and don''t move >= 5%%.\n' , sum(S.p<p) ); 
fprintf('%d\tgenes change expression and move >= 5%%.\n' , sum(M.p<p) );


fprintf('%0.02f%%\tof the genes that move >= 5%% change expression.\n' , 100*mean(M.p<p) ); 
fprintf('%0.02f%%\tof the genes that don''t move >= 5%% change expression.\n' , 100*mean(S.p<p) ); 
