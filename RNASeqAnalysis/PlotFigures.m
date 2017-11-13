%% Plot some figures
figname = 'FoldChange_vs_Pdiff.eps';
load('PP.mat');

X = NaN(0);
Y = NaN(0);
orfs = cell(0);
ID = NaN(0);
for I = 2:height(A)
    vnE =  [ 'Efc_' num2str(A.ID(I)) ] ;
    vnP =  [ 'Pdiff_' num2str(A.ID(I)) ] ;
    X = vertcat( X , PP.(vnP));
    Y = vertcat( Y , PP.(vnE));
    ID = vertcat( ID , repmat(A.ID(I) , numel(PP.(vnP)),1) );
    orfs = vertcat( orfs , PP.target_id);
end

%%

fh = figure('units','centimeters','position',[5 5 12 7]);
hold on ;
line([-1 100] , [nanmedian(Y) nanmedian(Y)],'LineStyle','--','Color',[.7 .7 .7])
X(X>30) = 31 ; 
bh = boxplot( Y , round(X/5)*5,'notch','on','Symbol','');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),[.7 .7 .7],'FaceAlpha',.5);
end
for I = 1:length(bh)
    set(bh(6,I),'LineWidth',3)
end
ylim([-0.75 0.75])
ylabel('Fold change in expression')
xlabel('Decrease in time spent in the nuclear periphery')
xlim([1.5 max(xlim)])
print('-dpsc2',figname,'-append');
close; 


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
print('-dpsc2',figname,'-append');
close; 

grps = round(X/5)*5 ; 
[p,tbl,stats] = kruskalwallis(Y,grps,'off');
[c,m,h,nms] = multcompare(stats,'Display','off');
pvals = squareform(c(:,6)) ;
pvals(: , 2)

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