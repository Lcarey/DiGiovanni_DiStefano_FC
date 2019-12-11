%% Plot some figures
%   how many genes increase expression? How many decrease
FIGNAME = 'FoldChange_vs_Pdiff';
WD = '~/Develop/DiGiovanni_DiStefano_FC/' ; 
load([ WD 'RNASeqAnalysis/PP.mat' ] );
T = table();
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
T.X = X;
T.FoldChangeExpression = Y;
T.ORF = orfs;
T.ID = ID;

G = grpstats( T( T.X>5 ,:) , 'ORF' , {'mean' 'median'} ,'DataVars','FoldChangeExpression');
%G.Properties.VariableNames{strcmp(G.Properties.VariableNames,'nanmean_FoldChangeExpression')} = 'mean_FoldChangeExpression' ; 
G.Properties.VariableNames{2} = 'N_FC_Strains' ; 

G.p = NaN(height(G),1);
for I = 1:height(G)
    [~, G.p(I) ] = ttest( T.FoldChangeExpression( strcmp(T.ORF , G.ORF{I})));
end
writetable( G , [WD '/Data/FoldChangeExpression__GenesThatMove_gt5.txt' ] ) ; 


Gall = grpstats( T  , 'ORF' , {'mean' 'median'} ,'DataVars','FoldChangeExpression');
Gall.Properties.VariableNames{2} = 'N_FC_Strains' ; 

Gall.p = NaN(height(Gall),1);
for I = 1:height(Gall)
    [~, Gall.p(I) ] = ttest( T.FoldChangeExpression( strcmp(T.ORF , Gall.ORF{I})));
end
writetable( Gall , '~/Develop/DiGiovanni_DiStefano_FC/Data/FoldChangeExpression__AllGenes.txt' ) ; 

%%
fh = figure('units','centimeters','position',[5 5 5 7]);
data = NaN(0);
data(3) = sum(G.p>0.05);
data(1) = sum(G.p<0.05 & G.mean_FoldChangeExpression < 0);
data(2) = sum(G.p<0.05 & G.mean_FoldChangeExpression > 0);
bh = bar( 1:3 , data ,'FaceColor', [.8 .8 .8] ) ; 
ylabel('# of genes that move > 5%')
set(gca,'xticklabel',{'< 0' '> 0' 'none'})
xlabel('Change in expression')

xlim([.5 3.5])
%title('Genes that change location > 5%  (FDR=10%)')
text( 2.6 , 40 , string(data(3)) ,'FontSize',15, 'FontWeight','bold')
text( 1.7 , 40 , string(data(2)) ,'FontSize',15 , 'Color' ,'b', 'FontWeight','bold')
text( 0.7 , 40 , string(data(1)) ,'FontSize',15 , 'Color' ,[0.1137 0.6941 0] , 'FontWeight','bold')

ylim([0 410])

%% volcano plot of fold-change vs p-value
x = G.mean_FoldChangeExpression; x(x<-2)=-2 ; x(x>2)=2 ;
y = -1 * log10(G.p) ; y(y>4)=4 ;
fh = figure('units','centimeters','position',[5 5 6.5 7]);
hold on ; 
plot( x(x>0) , y(x>0) , '.','Color', [0 .7 .7],'MarkerSize',10);
plot( x(G.p<0.05 & x>0) , y(G.p<0.05 & x>0) , '.','Color', 'b','MarkerSize',15);
plot( x(x<0) , y(x<0) , '.','Color', [0.3 .8 0.4],'MarkerSize',10);
plot( x(G.p<0.05 & x<0) , y(G.p<0.05 & x<0) , '.','Color', [0.1137 0.6941 0] ,'MarkerSize',15);
grid on 
set(gca,'xtick',-5:5)
line( [ 0 0 ] , ylim , 'LineStyle','-','Color','k')
xlim( [-2.1 2.1])
ylim([0 4.1])
xlabel('Fold change in expression')
ylabel('p-value')