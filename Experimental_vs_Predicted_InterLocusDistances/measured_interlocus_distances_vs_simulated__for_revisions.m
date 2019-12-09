% Figure 6. Validation of polymer models by live and fixed cell microscopy. 
% Fig 6E : simulations vs measured

FIGNAME = '~/Downloads/Figure6E__Simulated_vs_Measured_interlocus_distances__TRP1-LYS4' ; 
% load data
cd( '~/Develop/DiGiovanni_DiStefano_FC/Experimental_vs_Predicted_InterLocusDistances' )
S = readtable('simulations.txt','FileType','text');
S.X = round(S.X);
S = sortrows(S,{'X','Y'},'ascend');
S.whichpoint = categorical( repmat( {'low' 'mean' 'high'}' , 11 , 1) );
genotypes = upper( {'WT' 'FC(IV:XII)cen4' 'FC(IV:XII)cen12' 'FC(IV:XV)cen4' 'FC(IV:XV)cen15' ...
    'FC(IV:XV:V)cen4' 'FC(IV:XV:V)cen5' 'FC(IV:XV:XVI)cen4' 'FC(IV:XV:XVI)cen16' 'FC(IV:XV:V:VII)cen4'...
    'FC(IV:XV:V:VII)cen7'} )  ; 

S.genotypes = categorical( reshape( repmat(genotypes,3,1) , height(S) ,1) ) ;
S.genotype = cellfun(@(X)regexprep( X,'[:()]','_'),cellstr(S.genotypes) ,'UniformOutput',false) ;
S.genotype2 = cellstr( S.genotype);
S = sortrows(S,'genotype2','descend')

%% load experimental data
E = readtable('Results FC TRP1-LYS4 distances.xlsx');
E = stack(E,E.Properties.VariableNames) ; 
E.Properties.VariableNames = {'genotype' 'distance'} ;
E.genotype = cellfun(@(X)regexprep( X,'.*YMM\d+_',''),cellstr(E.genotype) ,'UniformOutput',false) ;
G = grpstats( E , 'genotype' , {'mean' 'median' 'std' 'sem'} , 'DataVars','distance');
G.genotype2 = cellstr( G.genotype);
G = sortrows(G,'genotype2','descend')

%%

figure; 
subplot(2,1,1)
errorbar( 1:11 , S.Y( S.whichpoint=='mean') ,  S.Y( S.whichpoint=='mean')-S.Y( S.whichpoint=='low') , S.Y( S.whichpoint=='high')- S.Y( S.whichpoint=='mean') ,'-ok','MarkerFaceColor',[.7 .7 .7])
axis tight; 
title('simulated')
set(gca,'xtick',1:11);
set(gca,'xticklabel',regexprep( cellstr( S.genotype(S.whichpoint=='mean')) , '_' ,' ') );
xlim([0.5 11.5])

subplot(2,1,2)
errorbar( 1:11 , G.mean_distance , G.sem_distance ,'-ok','MarkerFaceColor',[.7 .7 .7])
title('experimental')
set(gca,'xtick',1:11);
axis tight; 
set(gca,'xticklabel',regexprep( cellstr(G.genotype) , '_' ,' ') )
xlim([0.5 11.5])

[c,p] = corr( S.Y(S.whichpoint=='mean') ,  G.mean_distance , 'type','Spearman')
[c,p] = corr( S.Y(S.whichpoint=='mean') ,  G.mean_distance , 'type','Pearson')

%%
% G = experimental ;    S = simulated
fh = figure('units','centimeters','position',[5 5 8 8]);
hold on ;
herrorbar( S.Y( S.whichpoint=='mean') , G.mean_distance ,  S.Y( S.whichpoint=='mean')-S.Y( S.whichpoint=='low') , S.Y( S.whichpoint=='high')- S.Y( S.whichpoint=='mean') , '.k') ;
errorbar(  S.Y( S.whichpoint=='mean') , G.mean_distance ,  G.std_distance , G.std_distance ,'ok') ;
ylabel('Experimentally measured distance')
xlabel('Simulated inter-locus distance')
gscatter( S.Y( S.whichpoint=='mean') , G.mean_distance ,   regexprep( cellstr(G.genotype) , '_' ,' ') , parula(height(G)) ,'.',30)
print( '-dpng' ,[ FIGNAME  '-1' ] , '-r300') ; 
legend('off');
print( '-dpng' ,[ FIGNAME  '-2' ] , '-r300') ; 
close ; 