%% Figure 9. For genes further from the telomere, the predicted change in the time a 
% gene spends in the nuclear periphery is a better predictor of changes in gene expression. 
% Using data from all FC strains, we selected the subset of genes on chromosome arms that
% underwent fusion, and calculated the fold-change in expression (relative to wild-type), 
% the change in % peripheral, and the distance to the former telomere. Each point shows the 
% fold difference between (the ability of changes in % peripheral to predict expression) and
% (the ability of the log distance to the telomere to predict expression).  Each value is the 
% log2(r2%peripheral / r2dist-to-tel) for the set of genes that are < X kb from the former 
% telomere. Gene sets in which changes in % peripheral are better predictors of changes in
% expression ( log2(r2%peripheral / r2dist-to-tel) > 0 ) are colored blue. 
%
% LBC November 2019

%% load data
cd('~/Develop/DiGiovanni_DiStefano_FC/RNASeqAnalysis/');
load('PP.mat');
SGD = dataset2table( loadSGDFeatures() );
%high_expressed_idx =  PP.Expr_409 > prctile(PP.Expr_409 , 50) ; 
%PP = PP( high_expressed_idx , :) ; 
% Plot some figures
FIGNAME = '~/Downloads/Dist2Tel__vs__PctPeriphery';

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

% group by %PP
X(X>30) = 31 ; 
GROUPS = round(X/5)*5 ; 
[ug,n]=count_unique(GROUPS) ;
keep_groups = ug(n>=10);
Y = Y( ismember(GROUPS,keep_groups));
dt = dt(ismember(GROUPS,keep_groups));
dCEN = dCEN(ismember(GROUPS,keep_groups));
GROUPS = GROUPS( ismember(GROUPS,keep_groups));


%% Higher correlation w/fold change in expr, %P or dist-to-tel
xl  = [10000:50:1e5 ]  ;
c1 = NaN( 1 , numel(xl) );
c2 = NaN( 1 , numel(xl) );
c3 = NaN( 1 , numel(xl) );
c4 = NaN( 1 , numel(xl) );
c5 = NaN( 1 , numel(xl) );
p = NaN( 2 , numel(xl) );
corrtype = 'Pearson' ;
warning('off')
parfor I = 1:numel(xl)
  %  expr_to_keep_idx = expr_in_409 > prctile( expr_in_409 , 75 ) ;
    idx = dt < xl(I) ;  % threshold < X
    idx = (dt < (xl(I) + 10*1000 )) & (dt > (xl(I) - 10*1000 )) ;  % threshold +/- 10kb (moving window)
    idx = idx & idx_on_arm_with_deletion ; % & expr_to_keep_idx ;
    if sum(idx)>10
            [ c3(I) ] = corr( Y(idx) , dt(idx) ,'rows','complete','Type',corrtype) ;
            [ c4(I) ] = corr( Y(idx) , log10(dt(idx)) ,'rows','complete','Type',corrtype) ;
            [ c5(I) ] = corr( Y(idx) , X(idx) ,'rows','complete','Type',corrtype) ;
        ft = fittype( 'poly1' );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Robust = 'Bisquare';
        [xData, yData] = prepareCurveData( zscore(X(idx)), Y(idx) );
        if numel(xData)>100
            [fitresultPP, gofPP] = fit( xData, yData , ft, opts );
            [xData, yData] = prepareCurveData( zscore( ( dt(idx))) , Y(idx) );
            [fitresultD2T, gofD2T] = fit( xData, yData, ft, opts );
            c1(I) = gofD2T.rsquare ;
            c2(I) = gofPP.rsquare ;
        end
    end
 %   fprintf( '%d\t%d\t%d\t%d\t%0.02f\t%d-%d\n' , I , xl(I) , sum(idx) , numel(xData) , c1(I) , (xl(I) - 10*1000 ) , (xl(I) + 10*1000 )   ) ;
    if mod(I,50)==0 , fprintf('%d %0.0f%%\n' , I , (I/numel(xl))*100) , end  ;
end
c = vertcat(c1,c2);
%c2 = c ; c2(p<0.01) = NaN ;


% FIGURES %%%
fh = figure('units','centimeters','position',[5 5 10 9]) ; 
hold on ;
idx = c(1,:) > 0 & c(2,:) > 0  ;
yyy = log2(c(2,idx) ./ c(1,idx)) ; 
yyy(yyy<-0.2)=-0.2; yyy(yyy>0.2)=0.2;
xxx = xl(idx)./1000 ; 
plot(xxx,yyy,'.-r','LineWidth',3);
plot(xxx(yyy>0),yyy(yyy>0),'.-b','LineWidth',3)

set(gca','xscale','log')
ylim([-0.2 0.2])
set(gca,'xtick',[0:5:40 50 60 75 100 125 150 200] )
line(xlim,[0 0],'LineStyle','--','Color',[.7 .7 .7])
ylabel('log2( % peripheral / dist-to-telomere )')
xlim([21 100])
set(gca,'ytick',-1:0.1:1 )
%xlabel('Using only genes < Xkb from deteleted telomeres')
xlabel('Using only genes Xkb from deteleted telomeres')
print( '-dpng' , FIGNAME , '-r600') ; 
close ; 

% show actual r2 values
fh = figure('units','centimeters','position',[5 5 10 9]) ; 
hold on ;
plot(xl./1000 , c2 , '.-','LineWidth',2,'DisplayName','% peripheral');
plot(xl./1000 , c1 , '.-','LineWidth',2,'DisplayName','Dist-to-deleted-tel');
set(gca','xscale','log')
legend('location','nw')
%xlabel('Genes < Xkb from deteleted telomeres')\
xlabel('Genes Xkb from deteleted telomeres')
ylabel('Ability to predict changes in expression (r^2)')
set(gca,'xtick',[0:5:40 50 60 75 100 125 150 200] )
xlim([21 100])

print( '-dpng' , [ FIGNAME '_inset'] , '-r600') ; 
close ; 