% Plot some figures
% load data
figbasename = '~/Downloads/NoLargeExpressionChangesFigures_';
load('PP.mat');
IDs = PP.Properties.VariableNames(regexpcmp(PP.Properties.VariableNames,'Efc_')) ; 
IDs = regexprep( IDs , 'Efc_' ,'');
usids = unique(IDs);
%%

% Show that we don't have large expression changes
%typical X vs Y
for I = 1:numel(usids)
    fh = figure('units','centimeters','position',[5 5 7 7 ]);
    Y = log2( 0.1 + PP.(['Expr_' (usids{I}) ]) );
    X = log2( 0.1 + PP.Expr_409 );
    X = X(~isnan(Y)); Y = Y(~isnan(Y));
    dscatter(X,Y)
%    plot( X , Y , 'ok','MarkerFaceColor',[.7 .7 .7])
    xlabel('WT expression (log_{2}(TPM))');
    ylabel('FC expression (log_{2}(TPM))');
    title( usids{I} )
    set(gca,'xtick',0:1:10)
    set(gca,'ytick',0:1:10)
    xlim([0 12]);ylim(xlim)
    line([0 100],[0 100],'Color','k')
     line([0 100],[1 101],'Color','r')    
    line([1 101],[0 100],'Color','r')
    print('-dpng',[figbasename num2str(I) '.png'] , '-r300' );
    close;
end
 
%%
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

