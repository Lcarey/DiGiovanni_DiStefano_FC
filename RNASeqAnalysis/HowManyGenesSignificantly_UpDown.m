%% Plot some figures
figname = 'FoldChange_vs_Pdiff.eps';
load('PP.mat');
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
T.Y = Y;
T.ORF = orfs;
T.ID = ID;
%%
G = grpstats( T( T.X>5 ,:) , 'ORF' , {'mean' 'median'} ,'DataVars','Y');
G.p = NaN(height(G),1);
for I = 1:height(G)
    [~, G.p(I) ] = ttest( T.Y( strcmp(T.ORF , G.ORF{I})));
end

fh = figure('units','centimeters','position',[5 5 12 7]);
data = NaN(0);
data(3) = sum(G.p>0.05);
data(1) = sum(G.p<0.05 & G.mean_Y < 0);
data(2) = sum(G.p<0.05 & G.mean_Y > 0);
bh = bar(data ,'FaceColor', [.7 .7 .7] ) ; 
ylabel('# of genes that move > 5%')
set(gca,'xticklabel',{'fold change < 0' 'fold change >0' 'no change'})
xlim([.5 3.5])
%title('Genes that change location > 5%  (FDR=10%)')
text( [0.9 1.85 2.85] , [200 200 200] , string(data) ,'FontSize',15)
ylim([0 410])