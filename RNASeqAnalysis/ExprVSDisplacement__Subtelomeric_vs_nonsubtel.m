
load('PP.mat');
%high_expressed_idx =  PP.Expr_409 > prctile(PP.Expr_409 , 50) ; 
%PP = PP( high_expressed_idx , :) ; 
% Plot some figures
figname = '~/Downloads/FoldChange_vs_Pdiff.eps';

%% Arms with deletions

    
X = NaN(0);
Y = NaN(0);
dt = NaN(0);
orfs = cell(0);
ID = NaN(0);
idx_on_arm_with_deletion = NaN(0);
expr_in_409 = NaN(0);
for I = 2:height(A)
    vnE =  [ 'Efc_' num2str(A.ID(I)) ] ;
    vnP =  [ 'Pdiff_' num2str(A.ID(I)) ] ;
    expr = PP.(vnE) ; 
    expr_in_409 = vertcat( expr_in_409 , PP.Expr_409 ) ;
    gene_chr_arm_ids = cellfun( @(X)X(1:3) , PP.target_id ,'UniformOutput',false );
    these_arms_have_tel_deletions = unique(gene_chr_arm_ids(isnan(expr))) ; 
    idx_on_arm_with_deletion = vertcat( idx_on_arm_with_deletion , ismember(gene_chr_arm_ids,these_arms_have_tel_deletions) ); 
    X = vertcat( X , PP.(vnP));
    Y = vertcat( Y , expr );
    dt = vertcat( dt , PP.nt_to_closest_end ) ; 
    ID = vertcat( ID , repmat(A.ID(I) , numel(PP.(vnP)),1) );
    orfs = vertcat( orfs , PP.target_id);
end
idx_on_arm_with_deletion = logical(idx_on_arm_with_deletion) ; 

% group by %PP
X(X>30) = 31 ; 
GROUPS = round(X/5)*5 ; 
[ug,n]=count_unique(GROUPS) ;
keep_groups = ug(n>=10);
Y = Y( ismember(GROUPS,keep_groups));
dt = dt( ismember(GROUPS,keep_groups));
GROUPS = GROUPS( ismember(GROUPS,keep_groups));

% remove non-subtelomeric genes
idx_to_keep = dt > (33 * 1e3) ; %100kb
Y = Y(idx_to_keep);
GROUPS = GROUPS(idx_to_keep);


% Main figures for text

fh = figure('units','centimeters','position',[5 5 12 7]);
hold on ;
line([-1 100] , [nanmedian(Y) nanmedian(Y)],'LineStyle','--','Color',[.7 .7 .7])
bh = boxplot( Y , GROUPS ,'notch','on','Symbol','');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),[.7 .7 .7],'FaceAlpha',.5);
end
for I = 1:size(bh,2)
    set(bh(6,I),'LineWidth',3)
end
ylim([-0.75 0.75])
ylabel('Fold change in expression')
xlabel('Decrease in time spent in the nuclear periphery')
xlim([1.5 max(xlim)])
title([ 'non-subtelomeric (>33kb) # genes = ' num2str(numel(Y))])
%print('-dpsc2',figname,'-append');
%close; 