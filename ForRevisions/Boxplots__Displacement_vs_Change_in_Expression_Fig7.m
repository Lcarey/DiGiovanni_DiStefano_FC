%  Figure 7. Gene displacement away from the nuclear periphery correlates with increased expression. 
% boxplots of displacement vs change in expression for all genes,
% subtelomeric genes, and non-subtelomeric genes
%
% LBC November 2019


%% load data
BB = readtable('~/CareyLab/ExternalData/Brauer08/TableS1.xls');
ESR_UP_list = BB.ORF( strcmp(BB.ESR,'up'));
ESR_DN_list = BB.ORF( strcmp(BB.ESR,'down'));
GR_UP_list = BB.ORF(strcmp(BB.GrowthRateResponse_1_5SD,'up'));
GR_DN_list = BB.ORF(strcmp(BB.GrowthRateResponse_1_5SD,'down'));

cd('~/Develop/DiGiovanni_DiStefano_FC/RNASeqAnalysis/');
load('PP.mat');
SGD = dataset2table( loadSGDFeatures() );
%high_expressed_idx =  PP.Expr_409 > prctile(PP.Expr_409 , 50) ; 
%PP = PP( high_expressed_idx , :) ; 
% Plot some figures
FIGNAME = '~/Downloads/FoldChange_vs_Pdiff';

% calculate the distance to the centromere for each gene
SGD.KB2CEN = NaN( height(SGD),1);
for I = 1:height(SGD)
    cen_pos = SGD.Start( strcmp(SGD.Chr , SGD.Chr{I})  & strcmp(SGD.TYPE,'centromere'));
    if ~isempty( cen_pos )
        SGD.KB2CEN(I) = abs(SGD.Start(I) - cen_pos) ./ 1000 ;
    end
end
SGD = SGD( strcmp(SGD.TYPE,'ORF'),:);
PP = innerjoin( SGD( :  ,  {'ORF' 'KB2CEN'}) , PP , 'LeftKey','ORF','RightKey','target_id' ); 
%% Arms with deletions

    
X = NaN(0);
Y = NaN(0);
dt = NaN(0);
dCEN = NaN(0);
orfs = cell(0);
ID = NaN(0);
idx_on_arm_with_deletion = NaN(0);
expr_in_409 = NaN(0);
for I = 2:height(A)
    vnE =  [ 'Efc_' num2str(A.ID(I)) ] ;
    vnP =  [ 'Pdiff_' num2str(A.ID(I)) ] ;
    expr = PP.(vnE) ; 
    expr_in_409 = vertcat( expr_in_409 , PP.Expr_409 ) ;
    gene_chr_arm_ids = cellfun( @(X)X(1:3) , PP.ORF ,'UniformOutput',false );
    these_arms_have_tel_deletions = unique(gene_chr_arm_ids(isnan(expr))) ; 
    idx_on_arm_with_deletion = vertcat( idx_on_arm_with_deletion , ismember(gene_chr_arm_ids,these_arms_have_tel_deletions) ); 
    X = vertcat( X , PP.(vnP));
    Y = vertcat( Y , expr );
    dt = vertcat( dt , PP.nt_to_closest_end ) ; 
    dCEN = vertcat( dCEN , PP.KB2CEN );
    ID = vertcat( ID , repmat(A.ID(I) , numel(PP.(vnP)),1) );
    orfs = vertcat( orfs , PP.ORF);
end
idx_on_arm_with_deletion = logical(idx_on_arm_with_deletion) ; 
%% group by %PP
X(X>30) = 31 ; 
GROUPS = round(X/5)*5 ; 
[ug,n]=count_unique(GROUPS) ;
keep_groups = ug(n>=10);
Y = Y( ismember(GROUPS,keep_groups));
dt = dt(ismember(GROUPS,keep_groups));
dCEN = dCEN(ismember(GROUPS,keep_groups));
GROUPS = GROUPS( ismember(GROUPS,keep_groups));
%% optionally, remove genes that are part of the ESR or Growth Rate Response
DISCARD_FLAG  = false ; 

if DISCARD_FLAG
idx_to_discard_esr = ismember(orfs,ESR_UP_list) | ismember(orfs,ESR_DN_list) ; 
idx_to_discard_gr_response = ismember(orfs,GR_DN_list) | ismember(orfs,GR_UP_list) ; 
idx_to_keep = ~idx_to_discard_esr  & ~idx_to_discard_gr_response ;

orfs = orfs(idx_to_keep) ;
X = X(idx_to_keep) ;
Y = Y(idx_to_keep) ;
ID = ID(idx_to_keep) ;
GROUPS = GROUPS(idx_to_keep) ; 
end
%% Main figures for text
idx = dt > 50*1000 ; 
%idx = true(numel(dt),1);

fh = figure('units','centimeters','position',[5 5 12 7]);
hold on ;
line([-1 100] , [nanmedian(Y) nanmedian(Y)],'LineStyle','--','Color',[.7 .7 .7])
bh = boxplot( Y(idx) , GROUPS(idx) ,'notch','on','Symbol','');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),[.7 .7 .7],'FaceAlpha',.5);
end
for I = 1:size(bh,2)
    set(bh(6,I),'LineWidth',3)
end
ylim([-0.75 0.75])
ylabel('Fold change in expression')
xlabel('Decrease in % peripheral')
xlim([1.5 max(xlim)])
%title( 'all genes in all FC strains')

[~,tbl,stats] = anova1( Y(idx) , GROUPS(idx) ,'off');
[c,m,h,nms] = multcompare(stats,'Display','off') ;
p = c(c(:,1)==2 , 6) ; 
for I = 1:numel(p)
    text( I+1.7 , 0.75 , sprintf( '%0.0e\n' , p(I)) );% , 'BackgroundColor' , 'white' );
end

print('-dpsc2', [ FIGNAME sprintf('_%d',sum(idx)) ],'-append');
close; 

