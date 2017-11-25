%addpath(genpath('~/Develop/Matlab/'));
addpath(genpath('~/Develop/DiGiovanni_DiStefano_FC/'));
WD = '~/Develop/DiGiovanni_DiStefano_FC/' ;
RNA = 'Data/kallisto/' ;
LOC = 'Data/percentage_in_periphery_per_gene/' ;
cd(WD)
%% which genes deleted in each strain?
% run Mendoza17_ID_Deleted_ORFs_From_Annotation.m
 [ A , SGD ] = IDDeletedORFsFromAnnotation() ;

%% load data
strains = {'409' '524' '527' '1138' '1228' '1379' '1387' '1380' '1388' '1788' '1793'} ;
s = struct();
PP = table();
renormalization_vector = NaN(0);
eps = 0.1 ;
for SI = 1:numel(strains)
    L  = readtable( [ WD LOC 'percentage_in_periphery_per_gene_for_' strains{SI} '.txt'] ,'FileType','text','ReadVariableNames',false) ;
    L.Properties.VariableNames = {'target_id' 'PctPeriphery'} ;
    E1 =  readtable( [ WD RNA 'abundance_' strains{SI} '_rep1.tsv'  ] ,'FileType','text','ReadVariableNames',true,'Delimiter','\t'); E1.TPM1 = E1.tpm ; E1.Cnts1 = E1.est_counts ;
    E2 =  readtable( [ WD RNA 'abundance_' strains{SI} '_rep2.tsv'  ] ,'FileType','text','ReadVariableNames',true,'Delimiter','\t'); E2.TPM2 = E2.tpm ; E2.Cnts2 = E2.est_counts ;
    E3 =  readtable( [ WD RNA 'abundance_' strains{SI} '_rep3.tsv'  ] ,'FileType','text','ReadVariableNames',true,'Delimiter','\t'); E3.TPM3 = E3.tpm ; E3.Cnts3 = E3.est_counts ;
    E4 =  readtable( [ WD RNA 'abundance_' strains{SI} '_rep4.tsv'  ] ,'FileType','text','ReadVariableNames',true,'Delimiter','\t'); E4.TPM4 = E4.tpm ; E4.Cnts4 = E4.est_counts ;
    
    T = [ E1(:,{'target_id' 'Cnts1' 'TPM1'})   E2(:, {'Cnts2' 'TPM2'})  E3(:, {'Cnts3' 'TPM3'}) E4(:, {'Cnts4' 'TPM4'})];
    T.MedianExpr = median( [T.TPM1 T.TPM2 T.TPM3 T.TPM4] , 2);
    T.MedianCounts = median( [ T.Cnts1 T.Cnts2 T.Cnts3 T.Cnts4] , 2);
    T = innerjoin( T , L ) ;
    
    s(SI).T = T ;
    s(SI).strain = strains{SI} ;
    PP.( ['Expr_' strains{SI}]) = T.MedianExpr + eps ;
    PP.( ['PP_' strains{SI}]) = T.PctPeriphery + eps ;
    PP.( ['Counts_' strains{SI}]) = T.MedianCounts ;
    PP.target_id = T.target_id ;
    PP.chr = cellfun(@(X)double(X(2))  , PP.target_id) - 64 ; % to numeric Chr
    PP.pos = cellfun( @(X)str2double(X(4:6)) , PP.target_id) ;
end
PP = innerjoin(PP , SGD(:,{'ORF','GENE' 'Start' 'End'}) , 'LeftKey','target_id','RightKey','ORF');

%% for each gene, calc the distance to the start or end of chromosome
%  use whichever distance is smaller
G = grpstats(SGD , 'Chr' , 'max', 'DataVars' ,{'Start' 'End'});
G.max=max(G.max_Start,G.max_End);
PP.nt_to_closest_end = NaN( height(PP) , 1);
for I = 1:height(PP)
    min_from_left = min( PP.Start(I), PP.End(I));
    min_from_right = G.max( str2double(G.Chr) == PP.chr(I)) - min_from_left ; 
    PP.nt_to_closest_end(I) = min(min_from_left ,  min_from_right);
end

%% set deleted genes to NaN in PP
for I = 2:height(A)
    vn =  [ 'Expr_' num2str(A.ID(I)) ] ;
    expression_vect = PP.(vn);
    expression_vect( ismember(PP.target_id , A.deleted_orfs{I})) = NaN ; 
    PP.(vn) =  expression_vect ;
end
%% calculate 'correct' WT expression value by averaging across strains that don't change %P
PP.Expr = NaN( height(PP) , 1);
vn = PP.Properties.VariableNames;
vnPPidx = find( contains( vn  , 'PP_') );
THRESH = 2 ; 
for I = 1:height(PP)
    pp = table2array( PP(I,vnPPidx)) ;
    ppidx_doesnt_change = find( abs( pp(1) - pp) <= THRESH ) ;
    expr_idx_pp_doesnt_change = vnPPidx(ppidx_doesnt_change)-1 ; 
    PP.Expr(I) = nanmedian( table2array(PP(I,expr_idx_pp_doesnt_change)));
end

%% calc Expression Fold Change and %P diff
for I = 2:height(A)
    vnE =  [ 'Expr_' num2str(A.ID(I)) ] ;
    vnP =  [ 'PP_' num2str(A.ID(I)) ] ;
    PP.([ 'Pdiff_' num2str(A.ID(I)) ]) = PP.PP_409 - PP.(vnP) ; 
    PP.([ 'Efc_' num2str(A.ID(I)) ]) = log2( PP.(vnE)  ./ PP.Expr  );
end
save('PP.mat','PP' ,'A');
