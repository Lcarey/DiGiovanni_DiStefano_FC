% replot measured vs predicted distances for Figure 6C final revision
%  Please make the panels in Figure 6C larger (it's difficult so read the legend and numbers of the axes).
%  Also, in LYS4-NE and TRP1-NE, please add at least one more value to each axis so that it's that that it's a linear axis. 
% load data
cd('~/Develop/DiGiovanni_DiStefano_FC/Figure6C_Measured_vs_Predicted_distances');
fl = dir('*.txt');
clrs = cbrewer('qual','Set1',5);
for flI = 1:numel(fl)
    fl(flI).txt = regexprep( fl(flI).name , 'scatter_plot_' ,'') ;
    fl(flI).txt = regexprep( fl(flI).txt , '_distance.txt' ,'') ;
    fl(flI).txt = regexprep( fl(flI).txt , '_' ,'??') ;
    T = readtable( fl(flI).name , 'ReadVariableNames',false);
    fh = figure('units','centimeters','position',[5 5 8 8]);
    hold on ;
    line([0 2],[0 2],'LineStyle','--','Color',[.7 .7 .7]);
    for I = 1:height(T)
        h = errorbar(  T.Var2(I) , T.Var4(I) , T.Var5(I),T.Var5(I) , T.Var3(I), T.Var3(I) ,'Color',clrs(I,:) ,'LineWidth',2 ,'Marker','.','MarkerSize',15);
    end
    title(fl(flI).txt)
    xlabel('Predicted distance (\mum)')
    ylabel('Measured distance (\mum)')
    xlim([0 max(T.Var2+T.Var3)])
    ylim([0 max(T.Var4+T.Var5)])
    set(gca,'xtick',0:0.25:5)
    set(gca,'ytick',0:0.25:5)
    legend( vertcat('xy' , cellfun(@(X)X((end-1):end) , T.Var1,'UniformOutput',false)) , 'Location','nw')
    print('-dpng' , [ fl(flI).txt ] , '-r300' );
    close ; 
end
    