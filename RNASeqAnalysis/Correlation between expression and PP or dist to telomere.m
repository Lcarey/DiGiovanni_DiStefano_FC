%% Calculate correlation between expression & distance to the telomere in WT cells
load PP
Y1 = log10(PP.Expr) ; % .Expr has the RNAseq expression measurements normalized across all experiments for which the gene isn't predicted to change expr
%Y2 = log10(PP.Expr_409) ; % only WT data, probably slightly lower quality. corr = 0.99 between the two of them

Xlin = PP.nt_to_closest_end ; 
Xlog = log10( Xlin ) ; 
Xpp  = PP.PP_409 ; % predicted % peripheral in 409


T = table;
T.Y = Y1 ; 
T.Xlog = Xlog ; 
T.PP = Xpp ;

%% randomly sample X% of genes
rmat = NaN( 1000 , 6);
idx_top_25 = Y1 > prctile(Y1,75); 
idx_bottom_25 = Y1 < prctile(Y1,25); 
idx_mid_50 = Y1 >= prctile(Y1,25) & Y1 <= prctile(Y1,75); 

for I = 1:size(rmat,1)
    idx = random('Uniform',0,1,numel(Y1),1) <= 0.90 ; 
    rmat(I,1) = corr(Y1(idx & idx_bottom_25) , Xlog(idx & idx_bottom_25));
    rmat(I,2) = corr(Y1(idx & idx_bottom_25) , Xpp(idx & idx_bottom_25));
    rmat(I,3) = corr(Y1(idx & idx_mid_50) , Xlog(idx & idx_mid_50));
    rmat(I,4) = corr(Y1(idx & idx_mid_50) , Xpp(idx & idx_mid_50));
    rmat(I,5) = corr(Y1(idx & idx_top_25) , Xlog(idx & idx_top_25));
    rmat(I,6) = corr(Y1(idx & idx_top_25) , Xpp(idx & idx_top_25));
    rmat(I,7) = corr(Y1(idx) , Xlog(idx));
    rmat(I,8) = corr(Y1(idx) , Xpp(idx));
end
rmat = abs(rmat) ; 

%%
x1 = 10 ; 
x2 = 90 ; 
R=1 ; C = 4 ; 

figure; 
subplot(R,C,4)
boxplot(rmat(:,1:2),'symbol','','notch','on','Widths',0.8)
ylim( [ prctile(reshape(rmat(:,1:2),[],1),x1) prctile(reshape(rmat(:,1:2),[],1),x2) ])
set(gca,'ytick',[ min(get(gca,'ytick')) max(get(gca,'ytick')) ] )
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),'b','FaceAlpha',.5); 
patch(get(h(2),'XData'),get(h(2),'YData'),'r','FaceAlpha',.5); 
set(gca,'xtick',[])
title('bottom 25% by expression')


subplot(R,C,2)
boxplot(rmat(:,3:4),'symbol','','notch','on','Widths',0.8)
ylim( [ prctile(reshape(rmat(:,3:4),[],1),x1) prctile(reshape(rmat(:,3:4),[],1),x2) ])
set(gca,'ytick',[ min(get(gca,'ytick')) max(get(gca,'ytick')) ] )
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),'b','FaceAlpha',.5); 
patch(get(h(2),'XData'),get(h(2),'YData'),'r','FaceAlpha',.5); 
set(gca,'xtick',[])
title('mid 50% by expression')


subplot(R,C,3)
boxplot(rmat(:,5:6),'symbol','','notch','on','Widths',0.8)
ylim( [ prctile(reshape(rmat(:,5:6),[],1),x1) prctile(reshape(rmat(:,5:6),[],1),x2) ])
set(gca,'ytick',[ min(get(gca,'ytick')) max(get(gca,'ytick')) ] )
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),'b','FaceAlpha',.5); 
patch(get(h(2),'XData'),get(h(2),'YData'),'r','FaceAlpha',.5); 
set(gca,'xtick',[])
title('top 25% by expression')

subplot(R,C,1)
boxplot(rmat(:,7:8),'symbol','','notch','on','Widths',0.8)
ylim( [ prctile(reshape(rmat(:,7:8),[],1),x1) prctile(reshape(rmat(:,7:8),[],1),x2) ])
set(gca,'ytick',[ min(get(gca,'ytick')) max(get(gca,'ytick')) ] )
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),'b','FaceAlpha',.5); 
patch(get(h(2),'XData'),get(h(2),'YData'),'r','FaceAlpha',.5); 
set(gca,'xtick',[])
title('all genes')

%% shadedErrorBar plots
T = sortrows(T,'PP','descend');
X = movmedian( T.PP , 50);
Y = movmedian( T.Y , 50); 
S = movstd( T.Y , 1000);
fh = figure('units','centimeters','position',[5 5 5 5]) ;
shadedErrorBar(X,Y,S,{'b-','markerfacecolor','b'});    
axis tight; 
ylim([-0.4 2])
set(gca,'ytick',0:2);
xlabel('% peripheral (predicted)')
ylabel('Expression (log10(TPM))')
corr(T.PP,T.Y)

T = sortrows(T,'Xlog','descend');
X = movmedian( T.Xlog , 50);
Y = movmedian( T.Y , 50); 
S = movstd( T.Y , 1000);
fh = figure('units','centimeters','position',[5 5 5 5]) ;
shadedErrorBar(X,Y,S,{'r-','markerfacecolor','r'});  
axis tight ; 
ylim([-0.4 2]);
set(gca,'ytick',0:2);
set(gca,'xdir','rev')
set(gca,'xtick',2:6) % nt from end
set(gca,'xticklabels',[.1 1 10 100 1e3]) % nt from end
xlabel('kb from telomere')
ylabel('Expression (log10(TPM))')
corr(T.Xlog,T.Y)
