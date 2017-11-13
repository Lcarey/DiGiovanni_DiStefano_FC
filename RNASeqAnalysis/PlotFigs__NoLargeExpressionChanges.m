% Plot some figures
figname = 'NoLargeExpressionChangesFigures.eps';
delete(figname);
load('PP.mat');

% Show that we don't have large expression changes
%typical X vs Y
for I = 1:numel(usids)
    fh = figure('units','centimeters','position',[5 5 7 7 ]);
    Y = log10( 0.1 + PP.(['Expr_' num2str(usids(I)) ]) );
    X = log10( 0.1 + PP.Expr_409 );
    X = X(~isnan(Y)); Y = Y(~isnan(Y));
    dscatter(X,Y)
%    plot( X , Y , 'ok','MarkerFaceColor',[.7 .7 .7])
    xlabel('WT expression (log_{10}(TPM))');
    ylabel('FC expression (log_{10}(TPM))');
    title( num2str(usids(I)) )
    set(gca,'xtick',0:1:10)
    set(gca,'ytick',0:1:10)
    xlim([0 4]);ylim(xlim)
    line([0 5],[0 5],'Color','k')
    line([0 5],[0.5 5.5],'Color',[.7 .7 .7],'LineStyle','--')
    line([0.5 5.5],[0 5],'Color',[.7 .7 .7],'LineStyle','--')
    line([0 5],[1 6],'Color',[.7 .7 .7])    
    line([1 6],[0 5],'Color',[.7 .7 .7])
    print('-dpsc2',figname,'-append');
    close;
end
 
for I = 1:numel(usids)
    fh = figure('units','centimeters','position',[5 5 7 7 ]);
    Y = log2(  ( 0.1 + PP.(['Expr_' num2str(usids(I))] ))  ./ ( 0.1 + PP.Expr_409 ) );
    idx = ~isnan(Y) & ~isinf(Y);
    Y(Y<-2) = -2; Y(Y>2)=2 ; 
    dscatter(PP.PP_409(idx),Y(idx))
    ylim([-2.1 2.1])
    xlabel('% peripheral in WT');
    title( num2str(usids(I)) )
    set(gca,'xtick',0:25:100)
    set(gca,'ytick',-10:0.5:10)
    line(xlim,[0 0],'Color','k')
    line(xlim,[0.5 0.5],'Color',[.7 .7 .7],'LineStyle','--')
    line(xlim,[-.5 -.5],'Color',[.7 .7 .7],'LineStyle','--')
    line(xlim,[1 1],'Color',[.7 .7 .7])    
    line(xlim,[1 1],'Color',[.7 .7 .7])    
    ylabel('Fold change in expression');    
    print('-dpsc2',figname,'-append');
    close;
end

