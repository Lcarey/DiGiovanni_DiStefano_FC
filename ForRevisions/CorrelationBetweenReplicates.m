%% Plot correlations between biological replicates and across strains
% for Genetics revision
% LBC November 2017

%% load data
DATADIR = '~/Develop/DiGiovanni_DiStefano_FC/Data/kallisto/' ;
fl = dir( [ DATADIR '*.tsv'] );
exprmat = NaN(0);
rep_vect = NaN(0);
strain_vect = NaN(0);
for I = 1:numel(fl)
    repN = str2double( regexprep( regexprep( fl(I).name , '.tsv' ,'') , '.*rep' , '') )  ; 
    strainN = str2double( regexprep( regexprep( fl(I).name , '_rep.*' ,'') , '.*_' , '') ) ;
    T = readtable( [fl(I).folder filesep fl(I).name]  , 'FileType','text');
    exprmat  =  [ exprmat T.tpm] ;
    rep_vect = [ rep_vect repN ];
    strain_vect = [ strain_vect strainN ] ;
end

%% calculate correlations between biological replicates
unique_strains = unique(strain_vect);
cSpearman_between_bio_reps = NaN(0); 
cPearson_between_bio_reps = NaN(0); 
for ui = 1:numel(unique_strains)
    idx_this_strain = strain_vect == unique_strains(ui) ; 
    [c,p] = corr( exprmat( : , idx_this_strain) , 'Type' , 'Spearman') ;
    c = tril(c,-1) ;
    c = c(c>0);
    cSpearman_between_bio_reps = vertcat( cSpearman_between_bio_reps , c) ; 
    
    
    [c,p] = corr( log2( 1+ exprmat( : , idx_this_strain)) , 'Type' , 'Pearson') ;
    c = tril(c,-1) ;
    c = c(c>0);
    cPearson_between_bio_reps = vertcat( cPearson_between_bio_reps , c) ;     
end

% correlations among all experiments
[c,p] = corr( log2( 1+ exprmat ) , 'Type' , 'Pearson') ;
c = tril(c,-4) ;
c = c(c>0);

fh = figure('units','centimeters','position',[5 5 14 8]);
subplot(1,2,1)
hold on ; 
histogram( cPearson_between_bio_reps , 0.8:0.02:1 , 'Normalization','Probability') ;
histogram( c , 0.8:0.02:1 , 'Normalization','Probability') ;
%xlabel('Correlation between experiments')
ylabel('Fraction of pairs of experiments')
legend( {'bioglocal replicates' 'between strains'} ,'location','NorthOutside')
set(gca,'xtick',0:0.05:1)
xlim([0.8 1])

subplot(1,2,2)
hold on ; 
[f,x] = ecdf( cPearson_between_bio_reps ) ;
plot(x,f,'-','LineWidth',3);
[f,x] = ecdf( c ) ;
plot(x,f,'-','LineWidth',3);
ylabel('Cummulative fraction of pairs')
legend( {'bioglocal replicates' 'between strains'} ,'location','NorthOutside')
[~,p,kstat] = kstest2( c , cPearson_between_bio_reps) 
set(gca,'xtick',0:0.05:1)
xlim([0.8 1])

print( '-dpng' , '~/Downloads/correlations_between_biological_replicates_and_strains' , '-r600') ;
close ;