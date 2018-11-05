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

%% randomly sample X% of genes using correlation
rmat = NaN( 1000 , 6);
corr_type = 'Pearson' ; 
idx_top_25 = Y1 > prctile(Y1,75); 
idx_bottom_25 = Y1 < prctile(Y1,25); 
idx_mid_50 = Y1 >= prctile(Y1,25) & Y1 <= prctile(Y1,75); 

for I = 1:size(rmat,1)
    idx_d2t = Xlin > (30*1000) &  Xlin < (1e4*1000); 
    idx = random('Uniform',0,1,numel(Y1),1) <= 0.9 &  idx_d2t ; 
    rmat(I,1) = corr(Y1(idx & idx_bottom_25) , Xlog(idx & idx_bottom_25),'Type',corr_type);
    rmat(I,2) = corr(Y1(idx & idx_bottom_25) , Xpp(idx & idx_bottom_25),'Type',corr_type);
    rmat(I,3) = corr(Y1(idx & idx_mid_50) , Xlog(idx & idx_mid_50),'Type',corr_type);
    rmat(I,4) = corr(Y1(idx & idx_mid_50) , Xpp(idx & idx_mid_50),'Type',corr_type);
    rmat(I,5) = corr(Y1(idx & idx_top_25) , Xlog(idx & idx_top_25),'Type',corr_type);
    rmat(I,6) = corr(Y1(idx & idx_top_25) , Xpp(idx & idx_top_25),'Type',corr_type);
    rmat(I,7) = corr(Y1(idx) , Xlog(idx),'Type',corr_type);
    rmat(I,8) = corr(Y1(idx) , Xpp(idx) ,'Type',corr_type);
end
rmat = abs(rmat) ; 

% %% randomly sample using robust linear fit
% warning('off');
% ft = fittype( 'poly1' );
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% opts.Robust = 'Bisquare';
% rmat = NaN( 1000 , 6);
% idx_top_25 = Y1 > prctile(Y1,75); 
% idx_bottom_25 = Y1 < prctile(Y1,25); 
% idx_mid_50 = Y1 >= prctile(Y1,25) & Y1 <= prctile(Y1,75); 
% 
% for I = 1:size(rmat,1)
%     if mod(I,10)==0, fprintf('.');end;if mod(I,100)==0, fprintf('\n');end;
%         
%     idx = random('Uniform',0,1,numel(Y1),1) <= 0.90  ;
% 
%     [xData, yData] = prepareCurveData( Xlog(idx & idx_bottom_25) , Y1(idx & idx_bottom_25) );
%     [~, gof] = fit( xData, yData , ft, opts );
%     rmat(I,1) = gof.rsquare ; 
%     
%     [xData, yData] = prepareCurveData( Xpp(idx & idx_bottom_25) , Y1(idx & idx_bottom_25) );
%     [~, gof] = fit( xData, yData , ft, opts );
%     rmat(I,2) = gof.rsquare ; 
% 
%     [xData, yData] = prepareCurveData( Xlog(idx & idx_mid_50) , Y1(idx & idx_mid_50) );
%     [~, gof] = fit( xData, yData , ft, opts );
%     rmat(I,3) = gof.rsquare ; 
%     
%     [xData, yData] = prepareCurveData( Xpp(idx & idx_mid_50) , Y1(idx & idx_mid_50) );
%     [~, gof] = fit( xData, yData , ft, opts );
%     rmat(I,4) = gof.rsquare ; 
% 
%     [xData, yData] = prepareCurveData( Xlog(idx & idx_top_25) , Y1(idx & idx_top_25) );
%     [~, gof] = fit( xData, yData , ft, opts );
%     rmat(I,5) = gof.rsquare ; 
%     
%     [xData, yData] = prepareCurveData( Xpp(idx & idx_top_25) , Y1(idx & idx_top_25) );
%     [~, gof] = fit( xData, yData , ft, opts );
%     rmat(I,6) = gof.rsquare ; 
%         
%     [xData, yData] = prepareCurveData( Xlog(idx) , Y1(idx) );
%     [~, gof] = fit( xData, yData , ft, opts );
%     rmat(I,7) = gof.rsquare ; 
%     
%     [xData, yData] = prepareCurveData( Xpp(idx) , Y1(idx) );
%     [~, gof] = fit( xData, yData , ft, opts );
%     rmat(I,8) = gof.rsquare ; 
% end
% 

%% boxplots of correlation between expression & %P vs D2T in WT cells
x1 = 7.5 ; x1p = 0.05  ; 
x2 = 99 ; 
R=1 ; C = 4 ; 

fh = figure('units','centimeters','position',[5 5 24 6]);

subplot(R,C,4)
boxplot(rmat(:,1:2),'symbol','','notch','on','Widths',0.8)
ylim( [ prctile(reshape(rmat(:,1:2),[],1),x1) x1p+prctile(reshape(rmat(:,1:2),[],1),x1) ])
%set(gca,'ytick',[ min(get(gca,'ytick')) max(get(gca,'ytick')) ] )
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),'b','FaceAlpha',.5); 
patch(get(h(2),'XData'),get(h(2),'YData'),'r','FaceAlpha',.5); 
set(gca,'xtick',[])
title('bottom 25% by expr.')


subplot(R,C,2)
boxplot(rmat(:,3:4),'symbol','','notch','on','Widths',0.8)
ylim( [ prctile(reshape(rmat(:,3:4),[],1),x1) x1p+prctile(reshape(rmat(:,3:4),[],1),x1) ])
%set(gca,'ytick',[ min(get(gca,'ytick')) max(get(gca,'ytick')) ] )
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),'b','FaceAlpha',.5); 
patch(get(h(2),'XData'),get(h(2),'YData'),'r','FaceAlpha',.5); 
set(gca,'xtick',[])
title('mid 50% by expr.')


subplot(R,C,3)
boxplot(rmat(:,5:6),'symbol','','notch','on','Widths',0.8)
ylim( [ prctile(reshape(rmat(:,5:6),[],1),x1) x1p+prctile(reshape(rmat(:,5:6),[],1),x1) ])
%set(gca,'ytick',[ min(get(gca,'ytick')) max(get(gca,'ytick')) ] )
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),'b','FaceAlpha',.5); 
patch(get(h(2),'XData'),get(h(2),'YData'),'r','FaceAlpha',.5); 
set(gca,'xtick',[])
title('top 25% by expr.')

subplot(R,C,1)
boxplot(rmat(:,7:8),'symbol','','notch','on','Widths',0.8)
ylim( [ prctile(reshape(rmat(:,7:8),[],1),x1) x1p+prctile(reshape(rmat(:,7:8),[],1),x1) ])
%set(gca,'ytick',[ min(get(gca,'ytick')) max(get(gca,'ytick')) ] )
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
