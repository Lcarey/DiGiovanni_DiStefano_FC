BB = readtable('~/Data/Brauer08/TableS1.xls');
ESR_UP_list = BB.ORF( strcmp(BB.ESR,'up'));
ESR_DN_list = BB.ORF( strcmp(BB.ESR,'down'));
GR_UP_list = BB.ORF(strcmp(BB.GrowthRateResponse_1_5SD,'up'));
GR_DN_list = BB.ORF(strcmp(BB.GrowthRateResponse_1_5SD,'down'));
%%
cd('~/Develop/DiGiovanni_DiStefano_FC/RNASeqAnalysis/');
load('PP.mat');
SGD = dataset2table( loadSGDFeatures() );
%high_expressed_idx =  PP.Expr_409 > prctile(PP.Expr_409 , 50) ; 
%PP = PP( high_expressed_idx , :) ; 
% Plot some figures
figname = '~/Downloads/FoldChange_vs_Pdiff.eps';

%% calculate the distance to the centromere for each gene
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
title([ '# genes = ' num2str(numel(Y))])
%print('-dpsc2',figname,'-append');
%close; 


%% Xaxis is dist from telomere

GROUPS = round( ( dt./10000 )) * 10 ;

fh = figure('units','centimeters','position',[5 5 12 7]);
hold on ;
line([-1 100] , [nanmedian(Y) nanmedian(Y)],'LineStyle','--','Color',[.7 .7 .7])
bh = boxplot( Y(idx_on_arm_with_deletion) , GROUPS(idx_on_arm_with_deletion) ,'notch','on','Symbol','');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),[.7 .7 .7],'FaceAlpha',.5);
end
for I = 1:size(bh,2)
    set(bh(6,I),'LineWidth',3)
end
ylim([-0.75 0.75])
xlim([1.1 20])
ylabel('Fold change in expression')
xlabel('kb from the deleted telomere')
%xlim([1.5 max(xlim)])
title([ '# genes = ' num2str(sum(idx_on_arm_with_deletion))])
%print('-dpsc2',figname,'-append');
%close; 


%% Higher correlation w/fold change in expr, %P or dist-to-tel
xl  = 1:5:1e5 ;
c = NaN( 2 , numel(xl) );
p = NaN( 2 , numel(xl) );
corrtype = 'Spearman' ;
warning('off')
for I = 1:numel(xl)
  %  expr_to_keep_idx = expr_in_409 > prctile( expr_in_409 , 75 ) ;
    idx = dt < xl(I) ;
    idx = idx & idx_on_arm_with_deletion ; % & expr_to_keep_idx ;
    if sum(idx)>100
        %    [ c(2,I) , p(1,I)] = corr( Y(idx) , dt(idx) ,'rows','complete','Type',corrtype) ;
        %    [ c(3,I) , p(2,I)] = corr( Y(idx) , log10(dt(idx)) ,'rows','complete','Type',corrtype) ;
        %    [ c(1,I) , p(3,I)] = corr( Y(idx) , X(idx) ,'rows','complete','Type',corrtype) ;
        ft = fittype( 'poly1' );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Robust = 'Bisquare';
        [xData, yData] = prepareCurveData( zscore(X(idx)), Y(idx) );
        if numel(xData)>100
            [fitresultPP, gofPP] = fit( xData, yData , ft, opts );
            [xData, yData] = prepareCurveData( zscore(dt(idx)), Y(idx) );
            [fitresultD2T, gofD2T] = fit( xData, yData, ft, opts );
            c(1,I) = gofD2T.rsquare ;
            c(2,I) = gofPP.rsquare ;
        end
    end
    if mod(I,50)==0 , fprintf('%d %0.0f%%\n' , I , (I/numel(xl))*100) , end  ;
end
%c2 = c ; c2(p<0.01) = NaN ;
%%
fh = figure('units','centimeters','position',[5 5 10 8]) ; 
hold on ;
idx = c(1,:) > 0 & c(2,:) > 0  ;
yyy = log2(c(2,idx) ./ c(1,idx)) ; 
yyy(yyy<-0.2)=-0.2;
xxx = xl(idx)./1000 ; 
plot(xxx,yyy,'.-r','LineWidth',3);
plot(xxx(yyy>0 & xxx>30),yyy(yyy>0& xxx>30),'.-b','LineWidth',3)

set(gca','xscale','log')
ylim([-0.2 0.2])
set(gca,'xtick',[0:5:40 50 60 75 100 125 150 200] )
line(xlim,[0 0],'LineStyle','--','Color',[.7 .7 .7])
ylabel('log2( % peripheral / dist-to-telomere )')
xlim([21 100])
set(gca,'ytick',-1:0.1:1 )
xlabel('Using only genes < Xkb from deteleted telomeres')

%%
Yo = Y;
Yo(Yo>2) = 2;
Yo(Yo<-2)=-2;
fh = figure('units','centimeters','position',[5 5 12 7]);
hold on; 
sh = scatter( X , Yo , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
xlim([-1 max(X)])
ylim([-2.1 2.1])
ylabel('Fold change in expression')
xlabel('Decrease in time spent in the nuclear periphery')
line([-1 100] , [nanmedian(Y) nanmedian(Y)],'LineStyle','--','Color',[.7 .7 .7])
title([ 'N = ' num2str(numel(Y))])

print('-dpsc2',figname,'-append');
close; 

grps = round(X/5)*5 ; 
[p,tbl,stats] = kruskalwallis(Y,grps,'off');
[c,m,h,nms] = multcompare(stats,'Display','off');
pvals = squareform(c(:,6)) ;
pvals(: , 2)
close ; 

%% show that this is the same result for all strains
Q=table();
Q.ID = ID; 
Q.X = X ; 
Q.Y = Y;
Q.GRP = NaN(height(Q),1);
Q.GRP( Q.X>-1 & Q.X<1 ) = 0 ;
Q.GRP( Q.X>10 ) = 1 ;

G = grpstats( Q , {'ID' 'GRP'} ,{'mean' 'median' 'std'},'Datavars','Y');

fh = figure('units','centimeters','position',[5 5 12 7]);
bar(reshape(G.median_Y,2,[])');
legend({'no change in location' '>10% less peripheral'},'location','nw')
ylim([0 0.3])
ylabel('Fold change in expression')
set(gca,'xticklabel',A.ID(2:end))

set(gca,'ytick',arrayfun( @(X) round(log2(X/100)*100)/100 , 100:5:200))
ylabel('% increase in expression')
set(gca,'yticklabel',[0:5:100])
print('-dpsc2',figname,'-append');

%% plot area around STL1

%% ChrIV 1,507,200 - 1,516,800 : Displacement of away from the NE in FC strains
orfs_of_interest = PP.target_id( PP.chr == 4 & PP.Start > 1503000 & PP.Start <  1517000) ;
% orfs_of_interest = PP.target_id( PP.chr == 4 & PP.Start > 1507200 & PP.Start <  1516800) ;
ally = NaN(0);
idx = ismember(orfs , orfs_of_interest);
syms = 'o^phsvd<>+';
fh = figure('units','centimeters','position',[5 5 12 7]);
hold on ;
usids = unique(ID)
for I = 1:numel(usids)
    sidx = ID==usids(I);
    usids(I)
    gh = gscatter(  X(idx & sidx) ,  Y(idx & sidx) , orfs(idx & sidx) , parula(6) , syms(I));
    arrayfun(@(X)set(X,'MarkerFaceColor',get(X,'Color')) , gh);
end
axis tight;
lh = line( xlim , [0 0],'LineStyle','--','Color',[.7 .7 .7]) ; 
set(get(get(lh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylabel('Fold change in expression');
xlabel('Decrease in time spent in the nuclear periphery')

for I = orfs_of_interest'
    edata = Y( strcmp(orfs , char(I))) ;
    [~,p] = ttest( edata );
    fprintf('%s\t%0.05f\n' , char(I) , p)
end

figure; hold on;
for I = 1:numel(usids)
    plot(1,1,'ok','Marker',syms(I),'DisplayName',A.genotype{I+1}) ;
end
legend('location','best')


%% For reviewer: 
% 11.	Figure 7A: I would suggest fitting and plotting a loess regression or other solutions 
%  for showing the trend in the data, as otherwise it looks like a general shift up-ward of the data, 
%  which may be confounded by normalization issues
%
[X,o] = sort(X);
Y = Y(o);
yy100 = smooth( X , Y , 100 ,  'rloess' );
Xmod = X ; 
Xmod(Xmod>25) = 25 ; 


fh = figure('units','centimeters','position',[5 5 12 7]);
hold on; 
sh = scatter( Xmod , Y , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);

xlim([-1 max(Xmod)])
ylim([-2.1 2.1])
ylabel('Log_2 fold change in expression')
xlabel('Decrease in time spent in the nuclear periphery (predicted)')
line([-1 100] , [nanmedian(Y) nanmedian(Y)],'LineStyle','--','Color',[.7 .7 .7])

plot( Xmod , yy100 ,'-r','LineWidth',1)

ylim([-1 1])
print('-dpsc2','~/Downloads/NewFig7AForRevision.eps');
close;
