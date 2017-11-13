%addpath(genpath('~/Develop/Matlab/'));
WD = '~/Develop/DiGiovanni_DiStefano_FC/' ;
RNA = 'Data/kallisto/' ;
LOC = 'Data/percentage_in_periphery_per_gene/' ;
cd(WD)
%% which genes deleted in each strain?
% run Mendoza17_ID_Deleted_ORFs_From_Annotation.m
IDDeletedORFsFromAnnotation

%% load data
strains = {'409' '524' '527' '1138' '1228' '1379' '1387' '1380' '1388' '1788' '1793'} ;
%figure;  hold on ;
s = struct();
PP = table();
renormalization_vector = NaN(0);
eps = 0.1 ;
for SI = 1:numel(strains)
    L  = readtable( [ WD LOC 'percentage_in_periphery_per_gene_for_' strains{SI} '.txt'] ,'FileType','text','ReadVariableNames',false) ;
    L.Properties.VariableNames = {'target_id' 'PctPeriphery'} ;
    E1 =  readtable( [ WD RNA 'abundance_' strains{SI} '_rep1.tsv'  ] ,'FileType','text','ReadVariableNames',true,'Delimiter','\t'); E1.TPM1 = E1.tpm ; E1.Cnts1 = E1.est_counts ;
    E2 =  readtable( [ WD RNA 'abundance_' strains{SI} '_rep2.tsv'  ] ,'FileType','text','ReadVariableNames',true,'Delimiter','\t'); E2.TPM2 = E2.tpm ; E2.Cnts2 = E2.est_counts ;
    E3 =  readtable( [ WD RNA 'abundance_' strains{SI} '_rep3.tsv'  ] ,'FileType','text','ReadVariableNames',true,'Delimiter','\t'); E3.TPM3 = E3.tpm ; E3.Cnts3 = E3.est_counts ;
    E4 =  readtable( [ WD RNA 'abundance_' strains{SI} '_rep4.tsv'  ] ,'FileType','text','ReadVariableNames',true,'Delimiter','\t'); E4.TPM4 = E4.tpm ; E4.Cnts4 = E4.est_counts ;
    
    T = [ E1(:,{'target_id' 'Cnts1' 'TPM1'})   E2(:, {'Cnts2' 'TPM2'})  E3(:, {'Cnts3' 'TPM3'}) E4(:, {'Cnts4' 'TPM4'})];
    T.MedianExpr = median( [T.TPM1 T.TPM2 T.TPM3 T.TPM4] , 2);
    T.MedianCounts = median( [ T.Cnts1 T.Cnts2 T.Cnts3 T.Cnts4] , 2);
    T = innerjoin( T , L ) ;
    
    s(SI).T = T ;
    s(SI).strain = strains{SI} ;
    %    G = grpstats( T , 'PctPeriphery' , 'mean' ,'DataVars','MedianExpr');
    %    plot(G.PctPeriphery , G.mean_MedianExpr ,'DisplayName',strains{SI} ) ;
    PP.( ['Expr_' strains{SI}]) = T.MedianExpr + eps ;
    PP.( ['PP_' strains{SI}]) = T.PctPeriphery + eps ;
    PP.( ['Counts_' strains{SI}]) = T.MedianCounts ;
    PP.target_id = T.target_id ;
    PP.chr = cellfun(@(X)double(X(2))  , PP.target_id) - 64 ; % to numeric Chr
    PP.pos = cellfun( @(X)str2double(X(4:6)) , PP.target_id) ;
end
PP = innerjoin(PP , SGD(:,{'ORF','GENE' 'Start' 'End'}) , 'LeftKey','target_id','RightKey','ORF');

% distance from start or end of chromosome
%  use whichever distance is smaller
G = grpstats(SGD , 'Chr' , 'max', 'DataVars' ,{'Start' 'End'});
G.max=max(G.max_Start,G.max_End);
PP.nt_to_closest_end = NaN( height(PP) , 1);
for I = 1:height(PP)
    min_from_left = min( PP.Start(I), PP.End(I));
    min_from_right = G.max( str2double(G.Chr) == PP.chr(I)) - min_from_left ; 
    PP.nt_to_closest_end(I) = min(min_from_left ,  min_from_right);
end

% set deleted genes to NaN in PP
for I = 2:height(A)
    vn =  [ 'Expr_' num2str(A.ID(I)) ] ;
    expression_vect = PP.(vn);
    expression_vect( ismember(PP.target_id , A.deleted_orfs{I})) = NaN ; 
    PP.(vn) =  expression_vect ;
end
%% calculate 'correct' WT expression value by averaging across strains that don't change %P
PP.Expr = NaN( height(PP) , 1);
vn = PP.Properties.VariableNames;
vnPPidx = find( regexpcmp(vn , '^PP_'));
THRESH = 2 ; 
for I = 1:height(PP)
    pp = table2array( PP(I,vnPPidx)) ;
    ppidx_doesnt_change = find( abs( pp(1) - pp) <= THRESH ) ;
    expr_idx_pp_doesnt_change = vnPPidx(ppidx_doesnt_change)-1 ; 
    PP.Expr(I) = nanmedian( table2array(PP(I,expr_idx_pp_doesnt_change)));
end

%% calc Expression Fold Change and %P diff
for I = 2:height(A)
    vnE =  [ 'Expr_' num2str(A.ID(I)) ] ;
    vnP =  [ 'PP_' num2str(A.ID(I)) ] ;
    PP.([ 'Pdiff_' num2str(A.ID(I)) ]) = PP.PP_409 - PP.(vnP) ; 
    PP.([ 'Efc_' num2str(A.ID(I)) ]) = log2( PP.(vnE)  ./ PP.Expr  );
end

% Plot fold change vs %P diff
figname = 'FoldChange_vs_Pdiff.eps';
X = NaN(0);
Y = NaN(0);
orfs = cell(0);
for I = 2:height(A)
    vnE =  [ 'Efc_' num2str(A.ID(I)) ] ;
    vnP =  [ 'Pdiff_' num2str(A.ID(I)) ] ;
    X = vertcat( X , PP.(vnP));
    Y = vertcat( Y , PP.(vnE));
    orfs = vertcat( orfs , PP.target_id);
end

% extremes
idx = X>25 & ~isnan(Y);
EXTREMES = table();
EXTREMES.orf = orfs(idx) ;
EXTREMES.FC = Y(idx);
EXTREMES.X = X(idx);
EXTREMES = sortrows(EXTREMES,'orf')
%
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
ylim([-0.5 0.5])
ylabel('\Delta expression ( log_2(FC/WT) )')
xlabel('Increase in time spent in the periphery')
xlim([1.5 max(xlim)])
print('-dpsc2',figname,'-append');
close; 

grps = round(X/5)*5 ; 
[p,t,st] = anova1(Y,grps,'off');
[c,m,h,nms] = multcompare(st,'display','on');
for I = unique(grps)'
    [~,p] = ttest2( Y(grps==0) , Y(grps==I))
end

Yo = Y;
Yo(Yo>2) = 2;
Yo(Yo<-2)=-2;
fh = figure('units','centimeters','position',[5 5 12 7]);
hold on; 
sh = scatter( X , Yo , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
xlim([-1 max(X)])
ylim([-2.1 2.1])
ylabel('\Delta expression ( log_2(FC/WT) )')
xlabel('Increase in time spent in the periphery')
line([-1 100] , [nanmedian(Y) nanmedian(Y)],'LineStyle','--','Color',[.7 .7 .7])
print('-dpsc2',figname,'-append');
close; 



%%
xl = -2:0.1:2 ; 
figure; hold on; 
histogram(Y(grps==0),xl,'Normalization','Probability');
histogram(Y(grps==1),xl,'Normalization','Probability');
legend({'-1 - 1' '5 - 15'})
[~,p0] = ttest2(Y(grps==0) , Y(grps==1))
[~,p1] = ttest2(Y(grps==0) , Y(grps==-1))
[~,p2] = ttest2(Y(grps==0) , Y(grps==2))
[~,p3] = ttest2(Y(grps==0) , Y(grps==3))

title(sprintf( 'p = %0.02e' , p) )
xlabel('Expression fold change')
ylabel('Fraction of genes')

% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % %% Take the avg of each fusion. 
% % % % %   november 2017
% % % % T = table();
% % % % T.target_id = PP.target_id ; 
% % % % T.Expr_409 = PP.Expr_409 ;
% % % % T.PP_409 = PP.PP_409 ;
% % % % T.GENE = PP.GENE ;
% % % % ugd = unique(A.N_gd); ugd = ugd(ugd>0);
% % % % unique_genotypes = cell(0);
% % % % unique_genotypes_l = cell(0);
% % % % for I = 1:numel(ugd)
% % % %     Atmp = A( A.N_gd == ugd(I) , :);
% % % %     sids = Atmp.ID ; 
% % % %     E  = mean( [ PP.(['Expr_',num2str(sids(1))])  PP.(['Expr_',num2str(sids(2))])] , 2);
% % % %     P  = mean( [ PP.(['PP_',num2str(sids(1))])  PP.(['PP_',num2str(sids(2))])] , 2);
% % % %     
% % % %     idx = ismember( PP.target_id , Atmp.orfs_deleted{1});
% % % %     E(idx) = NaN; 
% % % %     P(idx) = NaN;
% % % %     FCid_l = regexprep( A.genotype{ find(A.N_gd==ugd(I),1)} , 'cen.*','');
% % % %     FCid = regexprep( FCid_l,'[-()]','') ;
% % % %     T.(['Expr_' FCid]) = E ; 
% % % %     T.(['PP_' FCid]) = P ; 
% % % %     T.(['Efc_' FCid]) = log2( E ./ T.Expr_409) ;
% % % %     T.(['Epc_' FCid]) = ((E - T.Expr_409) ./ T.Expr_409) * 100 ;
% % % %     T.(['Ppc_' FCid]) = ((T.PP_409 - P) ./ T.PP_409) * 100 ;
% % % %     T.(['Pdiff_' FCid]) = (T.PP_409 - P)  ;
% % % % 
% % % %     unique_genotypes{end+1} = FCid ; 
% % % %     unique_genotypes_l{end+1} = FCid_l ; 
% % % % end
% % % % 
% % % % for I = 1:height(T) , if isempty(T.GENE{I}) , T.GENE{I} = T.target_id{I}; end ; end ; 
% % % % 
% % % % %% Show that we don't have large expression changes
% % % % % typical X vs Y
% % % % figname = 'no_major_expression_changes.eps' ;
% % % % delete(figname)
% % % % for I = 1:numel(unique_genotypes)
% % % %     fh = figure('units','centimeters','position',[5 5 7 7 ]);
% % % %     Y = log10( 0.1 + T.(['Expr_' unique_genotypes{I}]) );
% % % %     X = log10( 0.1 + T.Expr_409 );
% % % %     X = X(~isnan(Y)); Y = Y(~isnan(Y));
% % % %     dscatter(X,Y)
% % % %     %plot( X , Y , 'ok','MarkerFaceColor',[.7 .7 .7])
% % % %     xlabel('WT expression (log_{10}(TPM))');
% % % %     ylabel('FC expression (log_{10}(TPM))');
% % % %     title( unique_genotypes_l{I} )
% % % %     set(gca,'xtick',0:1:10)
% % % %     set(gca,'ytick',0:1:10)
% % % %     xlim([0 4]);ylim(xlim)
% % % %     line([0 5],[0 5],'Color','k')
% % % %     line([0 5],[0.5 5.5],'Color',[.7 .7 .7],'LineStyle','--')
% % % %     line([0.5 5.5],[0 5],'Color',[.7 .7 .7],'LineStyle','--')
% % % %     line([0 5],[1 6],'Color',[.7 .7 .7])    
% % % %     line([1 6],[0 5],'Color',[.7 .7 .7])
% % % %     print('-dpsc2',figname,'-append');
% % % %     close;
% % % % end
% % % %  
% % % % for I = 1:numel(unique_genotypes)
% % % %     fh = figure('units','centimeters','position',[5 5 7 7 ]);
% % % %     Y = log2(  ( 0.1 + T.(['Expr_' unique_genotypes{I}]) )  ./ ( 0.1 + T.Expr_409 ) );
% % % %     idx = ~isnan(Y) & ~isinf(Y);
% % % %     Y(Y<-2) = -2; Y(Y>2)=2 ; 
% % % %     dscatter(T.PP_409(idx),Y(idx))
% % % %     ylim([-2.1 2.1])
% % % %     ylabel('\Delta expression log_2(FC/WT) TPM');
% % % %     xlabel('% peripheral in WT');
% % % %     title( unique_genotypes_l{I} )
% % % %     set(gca,'xtick',0:10:100)
% % % %     set(gca,'ytick',-10:0.5:10)
% % % %     line(xlim,[0 0],'Color','k')
% % % %     line(xlim,[0.5 0.5],'Color',[.7 .7 .7],'LineStyle','--')
% % % %     line(xlim,[-.5 -.5],'Color',[.7 .7 .7],'LineStyle','--')
% % % %     line(xlim,[1 1],'Color',[.7 .7 .7])    
% % % %     line(xlim,[1 1],'Color',[.7 .7 .7])    
% % % %     
% % % %     print('-dpsc2',figname,'-append');
% % % %     close;
% % % % end
% % % % 
% % % % 
% % % % %% ChrIV 1,507,200 - 1,516,800 : Displacement of away from the NE in FC strains
% % % % orfs = PP.target_id( PP.chr == 4 & PP.Start > 1503000 & PP.Start <  1519000) ;
% % % % %orfs = PP.target_id( PP.chr == 4 & PP.Start > 1507200 & PP.Start <  1516800) ;
% % % % ally = NaN(0);
% % % % idx = ismember(T.target_id , orfs);
% % % % syms = 'o+phsv';
% % % % fh = figure('units','centimeters','position',[5 5 12 12 ]);
% % % % hold on; 
% % % % for I = 1:numel(unique_genotypes)
% % % %     Y = T.(['Efc_' unique_genotypes{I}]) ;
% % % %     X = T.(['Pdiff_' unique_genotypes{I}]) ; 
% % % %     gh = gscatter(  X(idx) ,  Y(idx) , T.GENE(idx) , parula(numel(orfs)) , syms(I));
% % % %     arrayfun(@(X)set(X,'MarkerFaceColor',get(X,'Color')) , gh);
% % % %     ally = vertcat( ally ,  Y(idx));
% % % % end
% % % % lh = line( xlim , [0 0],'LineStyle','--','Color',[.7 .7 .7]) ; 
% % % % set(get(get(lh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% % % % ylabel('\Delta expression log_2(FC/WT) TPM');
% % % % xlabel('\Delta % peripheral (WT - FC)');
% % % % 
% % % % %%
% % % % figure; hold on;
% % % % for I = 1:numel(unique_genotypes_l)
% % % %     plot(1,1,'ok','Marker',syms(I),'DisplayName',unique_genotypes_l{I});
% % % % end
% % % % legend('location','best')
% % % % %%
% % % % vn = PP.Properties.VariableNames;
% % % % vn = vn(regexpcmp(vn,'^PP'));
% % % % rmat = NaN(numel(vn));
% % % % for I = 1:numel(vn)
% % % %     for J = 1:numel(vn)
% % % %         rmat(I,J) = corr( log10( PP.(vn{I})+0.1) , log10( PP.(vn{J})+0.1) ,'rows','complete');
% % % %     end
% % % % end
% % % % figure; 
% % % % imagesc(rmat,[0.95 1])
% % % % set(gca,'xticklabel',A.genotype)
% % % % xticklabel_rotate([],45)
% % % % title('correlation in %Peripherial')
% % % % % % % % %% 
% % % % % % % % PPC1 = NaN(0);
% % % % % % % % PPC0 = NaN(0);
% % % % % % % % for I = 1:numel(unique_genotypes)
% % % % % % % %     X = T.(['Pdiff_' unique_genotypes{I} ]) ;
% % % % % % % %     Y = T.(['Epc_' unique_genotypes{I}]) ;
% % % % % % % %     
% % % % % % % %     PC1 = Y(X>5 & X<15) ; PPC1 = vertcat(PC1,PPC1) ; 
% % % % % % % %     PC0 = Y(X>-1 & X<1) ; PPC0 = vertcat(PC0,PPC0) ; 
% % % % % % % %    % figure ; hold on; 
% % % % % % % %    % histogram(PC0, -100:5:100 ,'Normalization','Probability');
% % % % % % % %   %  histogram(PC1, -100:5:100,'Normalization','Probability');
% % % % % % % %     plot(median(PC0),median(PC1),'ok')
% % % % % % % %  %   [~,p]=ttest2(PC1,PC0);
% % % % % % % % %    title(  sprintf( '%s %0.03f' , unique_genotypes{I} , p ) )
% % % % % % % % end
% % % % %%  For each fusion, plot the difference in expression 
% % % % %  for genes that don't change vs genes that do change 
% % % % %fh = figure('units','centimeters','position',[5 5 12 12 ]);
% % % % %hold on ; 
% % % % Xo = NaN(0);
% % % % Yo = NaN(0);
% % % % Eo = NaN(0);
% % % % Elr = NaN(0);
% % % % Do = NaN(0);
% % % % Go = NaN(0);
% % % % data = NaN(0);
% % % % genotypes = cell(0);
% % % % MAXPCTCHANGE = 20 ; 
% % % % for SI = 2:numel(strains)
% % % %         genotype = A.genotype{ A.ID == str2double(strains{SI})};
% % % %         genotypes{end+1} = genotype ; 
% % % %         orfs_deleted = A.orfs_deleted{ A.ID== str2double(strains{SI}) } ;
% % % %         idx = ~ismember(PP.target_id, orfs_deleted) ;
% % % %         tmpY = ((PP.( ['Expr_' strains{SI}]) - PP.Expr_409 ) ./ PP.Expr_409 ) .* 100 ;
% % % %         tmpYlr = log2(PP.( ['Expr_' strains{SI}]) ./ PP.Expr_409 )   ;
% % % %         Elr = [ Elr tmpYlr(idx)'];
% % % %         tmpX = (PP.PP_409 - PP.( ['PP_' strains{SI}]) )  ;
% % % %         Yo = [Yo  tmpY(idx)'];
% % % %         Xo = [Xo  tmpX(idx)'] ;
% % % %         Do = [Do  PP.nt_to_closest_end(idx)'] ;
% % % %         Eo = [Eo  PP.Expr_409(idx)'] ;
% % % %         Go = [Go PP.chr(idx)'];
% % % %         PC0  = Yo(Xo>-1 & Xo<1); 
% % % %         PC1 = Yo(Xo>5 & Xo<MAXPCTCHANGE) ; 
% % % %         %plot(median(PC0),median(PC1),'.k','MarkerSize',0.1);
% % % %         %text( median(PC0) , median(PC1) , genotype);
% % % %       %  data(end+1 , 1 ) = median(PC1) - median(PC0);
% % % %         pc0b = (bootstrp( 1e1 , @median , PC0)); 
% % % %         pc1b = (bootstrp( 1e1 , @median , PC1)); 
% % % %         data(end+1 , 1) = mean( (pc1b ./ pc0b) );
% % % %         data(end , 2) = std( (pc1b ./ pc0b) ) ; 
% % % % 
% % % % end
% % % % %xlim([0 12])
% % % % %ylim(xlim)
% % % % %line(xlim,xlim)
% % % % %xlabel('median \Delta expression for genes that don''t change location')
% % % % %ylabel('median \Delta expression when 5-15% less peripheral')
% % % % 
% % % % 
% % % % 
% % % % %%
% % % % % bar graph of differences
% % % % figname = ['effect_moving_NOchange_to_5' 'div' num2str(MAXPCTCHANGE) '%_change.eps' ] ; 
% % % % delete(figname);
% % % % % fh = figure('units','centimeters','position',[5 5 10 25 ]);
% % % % % hold on ;
% % % % % bar( data(:,1) ,'FaceColor',[.7 .7 .7])
% % % % % %is_negerrb = data(:,1) - data(:,2) ;
% % % % % %negerrb = data(:,2);
% % % % % %negerrb(is_negerrb<0) = 0 ;
% % % % % errorbar( 1:nrows(data) , data(:,1) , data(:,2) ,'ok');
% % % % % set(gca,'xtick',1:10); 
% % % % % set(gca,'xticklabel',genotypes) ;
% % % % % rotateticklabel(gca,10,45);
% % % % % print('-dpsc2',figname,'-append');
% % % % % close; 
% % % % 
% % % % 
% % % % fh = figure('units','centimeters','position',[5 5 10 10 ]);
% % % % hold on ;
% % % % bar( data(:,1) ,'FaceColor',[.7 .7 .7])
% % % % errorbar( 1:nrows(data) , data(:,1) , data(:,2) ,'ok');
% % % % set(gca,'xtick',1:10); 
% % % % set(gca,'xticklabel',[]) ;
% % % % % line( xlim , [1 1] , 'LineStyle','--' ,'Color',[.5 .5 .5])
% % % % ylim([ 0 max(ylim)])
% % % % print('-dpsc2',figname,'-append');
% % % % close; 
% % % % %%
% % % % fh = figure('units','centimeters','position',[5 5 12 7 ]);
% % % % hold on ;
% % % % dsX = Xo(~isinf(Elr) & ~isnan(Elr))' ;
% % % % dsY = Elr(~isinf(Elr) & ~isnan(Elr))' ; 
% % % % dsY(dsY>2)=2; dsY(dsY<-2)=-2;
% % % % dscatter( dsX , dsY )
% % % % axis tight ;
% % % % grid on ; 
% % % % line( xlim , [0 0] ,'LineStyle','--')
% % % % xlabel('\Delta in % peripheral (WT-FC)')
% % % % ylabel('\Delta expression (log_2(FC/WT))')
% % % % %%
% % % % figure; hold on; 
% % % % scatter(Xo,Yo,7,log10(Eo),'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);  
% % % % grid on ;    
% % % % %line( xlim , [ median(Yo(:)) median(Yo(:)) ])
% % % % %line( xlim , [ mean(Yo(:)) mean(Yo(:)) ])
% % % % set(gca,'ColorOrder',[0 0 0 ; 0 0 0 ; 1 0 0 ; 1 0 0 ; 0 1 0 ; 0 1 0 ; 0 0 1 ; 0 0 1 ; 1 1 0 ; 1 1 0])
% % % % [X,o] = sort(Xo); Y = Yo(o);
% % % % plot( X  , movmedian( Y  , round(numel(X)/500) ) ,'LineWidth',1)
% % % % xlabel('\Delta % Peripheral (WT-mut)')
% % % %  ylabel('\DeltaExpression (log2(WT/mut))')
% % % %   title('only >100 kb, log10(dist)')
% % % %   title('colored by log10(Expr WT), threshold at 100 TPM')
% % % % colorbar
% % % % %ylim([-2.1 2.1])
% % % % 
% % % % %% high expression and high DeltaX genes
% % % % high_delta_x_high_expr_genes = unique(Go( Xo>20 ));
% % % % PP( ismember( PP.target_id,high_delta_x_high_expr_genes) ,:)
% % % % 
% % % % %%
% % % % figure; 
% % % % boxplot( Yo , round(Xo/5)*5 ,'notch','on')
% % % % grid on ;    
% % % % line( xlim , [ median(Y(:)) median(Y(:)) ])
% % % % line( xlim , [ mean(Y(:)) mean(Y(:)) ])
% % % % xlabel('\Delta % Peripheral (WT-mut)')
% % % % ylabel('\DeltaExpression (log2(WT/mut))')
% % % % ylim([-0.5 0.5])
% % % % %% 
% % % % %  plot each gene once
% % % % X = mean(Xo,2) ;
% % % % Y = mean(Yo,2) ;
% % % % [X,o] = sort(X);
% % % % Y = Y(o);
% % % % figure; hold on; 
% % % % plot(X,Y,'.','Color',[.7 .7 .7])    
% % % % grid on ;    
% % % % line( xlim , [ median(Y) median(Y) ])
% % % % line( xlim , [ mean(Y) mean(Y) ])
% % % % plot( X , movmedian( Y , round(numel(X)/100) ) ,'LineWidth',4)
% % % % 
% % % % xlabel('\DeltaPeripheral (WT-mut)')
% % % %  ylabel('\DeltaExpression (log2(WT/mut))')
% % % %  title('each gene once')
% % % %  ylim([-0.5 0.5])
% % % %  
% % % % % plot each fusion once
% % % % X = [ mean( [ Xo(:,1) Xo(:,2)] , 2)     mean( [ Xo(:,3) Xo(:,4)] , 2)   mean( [Xo(:,5) Xo(:,6)] , 2)     mean( [ Xo(:,7) Xo(:,8)] , 2) mean( [ Xo(:,9) Xo(:,10)] , 2)]; 
% % % % Y = [ mean( [ Yo(:,1) Yo(:,2)] , 2)     mean( [ Yo(:,3) Yo(:,4)] , 2)   mean( [Yo(:,5) Yo(:,6)] , 2)     mean( [ Yo(:,7) Yo(:,8)] , 2) mean( [ Yo(:,9) Yo(:,10)] , 2)]; 
% % % % X = X(:);
% % % % Y = Y(:);      
% % % % [X,o] = sort(X);
% % % % Y = Y(o);
% % % % figure; hold on; 
% % % % plot(X,Y,'.','Color',[.7 .7 .7])    
% % % % grid on ;    
% % % % line( xlim , [ median(Y) median(Y) ])
% % % % line( xlim , [ mean(Y) mean(Y) ])
% % % % plot( X , movmedian( Y , round(numel(X)/100) ) ,'LineWidth',4)
% % % % 
% % % % xlabel('\DeltaPeripheral (WT-mut)')
% % % %  ylabel('\DeltaExpression (log2(WT/mut))')
% % % %   title('each fusion once')
% % % %  ylim([-0.5 0.5])
% % % % 
% % % % %% October 25, 2017
% % % % 
% % % % orfs_deleted = A.orfs_deleted{ A.ID==527} ;
% % % % idx = ismember(PP.target_id, orfs_deleted) ;
% % % % %idxHIGHexpr = PP.Counts_409>10 | PP.Counts_524>10 | PP.Counts_527>10; 
% % % %  X1 =  log2( PP.PP_409 ./  PP.PP_524) ;
% % % %  X2 =  log2( PP.PP_409 ./  PP.PP_527) ;
% % % %  X1 =  ( PP.PP_409 -  PP.PP_524) ;
% % % %  X2 =  ( PP.PP_409 -  PP.PP_527) ;
% % % % 
% % % %  Y1 =  log2( PP.Expr_409 ./  PP.Expr_524) ;
% % % %  Y2 =  log2( PP.Expr_409 ./  PP.Expr_527) ;
% % % % 
% % % %  X = mean( [X1 X2] , 2) ;
% % % %  Y = mean( [Y1 Y2] , 2) ;
% % % % 
% % % %  
% % % %  figure ; 
% % % %  hold on ;
% % % %  plot( X1 , Y1 ,'.r')
% % % %  plot( X2 , Y2 ,'.b')
% % % % %plot( X , Y , '.k','MarkerSize',10)
% % % %  plot( X1(idx) , Y1(idx) ,'om')
% % % %  plot( X2(idx) , Y2(idx) ,'oc')
% % % %  
% % % %  
% % % % % plot( X(idxHIGHexpr) , Y(idxHIGHexpr) ,'ok')
% % % %  
% % % % % errorbar(  0 , median(Y(X<1 & X>-1)) , std(Y(X<1 & X>-1)) ,'Color','k','LineWidth',2);
% % % %  xlabel('\Delta%Peripheral (WT - mut)')
% % % %  ylabel('\DeltaExpression (log2(WT/mut))')
% % % %  
% % % % grid on 
% % % % 
% % % % 
% % % % of_interest = find(X>20)
% % % % for I = 1:numel(of_interest)
% % % %     text( X(of_interest(I)) , Y(of_interest(I)) , PP.target_id{of_interest(I)});
% % % % end
% % % % 
% % % % %%
% % % % figure; hold on; 
% % % % ecdf(Y1( X1 > -0.5 & X1 < 0.5)  )
% % % % ecdf(Y1( X1 > 0.5 & X1 < 2)  )
% % % % ecdf(Y1( X1 > 2)  )
% % % % legend({'no change' '0.5 - 2' '>2'})
% % % % title('525')
% % % % 
% % % % figure; hold on; 
% % % % ecdf(Y2( X2 > -0.5 & X2 < 0.5)  )
% % % % ecdf(Y2( X2 > 0.5 & X2 < 2)  )
% % % % ecdf(Y2( X2 > 2)  )
% % % % legend({'no change' '0.5 - 2' '>2'})
% % % % title('527')
% % % % 
% % % % 
% % % % %%
% % % % X = NaN(0) ;
% % % % Y = NaN(0) ;
% % % % SI = 1;
% % % % eps = 0.001 ;
% % % % high_expr = true(0);
% % % % for SJ = 2:numel(strains)
% % % %     E1 = (eps + PP.( ['Expr_' strains{SI}]))  ;
% % % %     E2 = (eps + PP.( ['Expr_' strains{SJ}]))  ;
% % % %     E2 = E2 ./ sum(E2) ; 
% % % %     E1 = E1 ./ sum(E1) ;
% % % %     EXPR_I_Over_J = log2( E1 ./ E2 ) ;
% % % %     
% % % %     PP_I_Minus_J = PP.( ['PP_' strains{SI}]) -   PP.( ['PP_' strains{SJ}]) ;
% % % %     
% % % %     X = vertcat( X , PP_I_Minus_J  )  ;
% % % %     Y = vertcat( Y , EXPR_I_Over_J )  ;
% % % %     high_expr = vertcat( high_expr , (PP.( ['Counts_' strains{SJ}]) > 100) & (PP.( ['Counts_' strains{SI}]) > 100)  )  ;
% % % % end
% % % % 
% % % % [~,o] = sort(X);
% % % % X = X(o); 
% % % % Y = Y(o) ; 
% % % % high_expr = high_expr(o);
% % % % 
% % % % figure; hold on ;
% % % % idx = high_expr & X>0 & abs(Y)< 1000 ;
% % % % idx = high_expr ; 
% % % % plot(X,Y,'.','Color',[.7 .7 .7]);
% % % % 
% % % % [xData, yData] = prepareCurveData( X(idx), Y(idx) );
% % % % ft = fittype( 'poly1' );
% % % % opts = fitoptions( 'Method', 'LinearLeastSquares' );
% % % % opts.Robust = 'Bisquare';
% % % % [fitresult, gof] = fit( xData, yData, ft, opts );
% % % %  boxplot( yData , round(xData/5)*5 ,'notch','on' )
% % % % plot( fitresult  );
% % % % 
% % % % ci=confint(fitresult)
% % % % title( sprintf('slope=%0.04f(+/-%0.04f) R^2=%0.02f' , fitresult.p1, abs(ci(1,1)-fitresult.p1) , gof.rsquare ) );
% % % % ylim([-0.5 0.5])
% % % % grid on ;
% % % % %xlim([0.5 7.5])
% % % % 
% % % % plot(xData , movmedian( yData , 1000),'.k');
% % % % %%
% % % % rsq = NaN(0) ;
% % % % slp = NaN(0) ;
% % % % mx  = NaN(0) ;
% % % % c   = NaN(0) ;
% % % % p   = NaN(0) ;
% % % % uchr = unique(PP.chr);
% % % % for SJ = 2:numel(strains)
% % % %     Y = log2( (eps + PP.( ['Expr_' strains{1}])) ./ (eps + PP.( ['Expr_' strains{SJ}])) ) ;
% % % %     X = PP.( ['PP_' strains{1}]) -   PP.( ['PP_' strains{SJ}]) ;
% % % %     high_expr = (PP.( ['Counts_' strains{SJ}]) > 10) & (PP.( ['Counts_' strains{SI}]) > 10)    ;
% % % %     for I = 1:numel(uchr)
% % % %         idx = high_expr & X > 0 & PP.chr==uchr(I) ;
% % % %         try
% % % %         [xData, yData] = prepareCurveData( X(idx), Y(idx) );
% % % %         ft = fittype( 'poly1' );
% % % %         opts = fitoptions( 'Method', 'LinearLeastSquares' );
% % % %         opts.Robust = 'Bisquare';
% % % %         [fitresult, gof] = fit( xData, yData, ft, opts );
% % % %         [c(end+1),p(end+1)]=corr(xData,yData) ;
% % % %         %  plot( fitresult, xData, yData );
% % % %         rsq(end+1) = gof.rsquare ;
% % % %         slp(end+1) = fitresult.p1 ; 
% % % %         mx(end+1) = mean(xData) ;
% % % %         catch
% % % %         end
% % % %     end
% % % % end
% % % % 
% % % % %%  %% Genes that change localization a lot are very lowly expressed %% 
% % % % fh = figure('units','centimeters','position',[5 5 8 8]);
% % % % hold on ;
% % % % [f,x] = ecdf(log10(Eo)); plot(x,f,'-k','LineWidth',2)
% % % % [f,x] = ecdf(log10(Eo(Xo>10 & Xo <20))) ;  plot(x,f,'-','LineWidth',2)
% % % % [f,x] = ecdf(log10(Eo(Xo>20))) ; plot(x,f,'-','LineWidth',2)
% % % % legend({'all' '\Delta%P>10 & \Delta%P<20' '\Delta%P>20'}, 'location','se')
% % % % xlabel('Expression in WT (TPM)')
% % % % ylabel('Fraction of genes')
% % % % ylabel('% of genes')
% % % % xlim([0 3])
% % % % set(gca,'ytick',0:.1:2)
% % % % set(gca,'xtick',log10([1 2 5 10 25 50 100 250 500 1000]))
% % % % set(gca,'xticklabel', [1 2 5 10 25 50 100 250 500 1000] )
% % % % grid on;
% % % % 
% % % % %% new idea: for each gene, get expression for all strains where %PP doesn't change from 409
% % % % od = cell(0);
% % % % for I = 1:height(A)
% % % %     od = unique( vertcat( A.orfs_deleted{I}' , od) ) ;
% % % % end
% % % % 
% % % % R = table();
% % % % vn = PP.Properties.VariableNames ;
% % % % expr_vars = vn( regexpcmp( vn , 'Expr_'));
% % % % pp_vars = vn( regexpcmp( vn , 'PP_'));
% % % % 
% % % % R.target_id = unique(PP.target_id);
% % % % R.expr = cell( height(R) , 1);
% % % % R.pp = cell( height(R) , 1);
% % % % R.c = NaN( height(R) , 1);
% % % % R.gene = cell( height(R) , 1);
% % % % 
% % % % for I = 1:numel(R.target_id)
% % % %     idx = strcmp( PP.target_id , R.target_id{I});
% % % %     R.expr{I} = table2array(PP( idx , expr_vars)) ;
% % % %     R.pp{I} = table2array(PP( idx , pp_vars)) ;
% % % %     R.c(I) = corr( R.expr{I}' , R.pp{I}');
% % % %     R.gene{I} = PP.GENE{ idx } ;
% % % % 
% % % % end
% % % %     
% % % %  %%
% % % %  R.pprange = cellfun(@range, R.pp);
% % % %  R.exprrange = cellfun(@range, R.expr);
% % % %  R.exprrange_log2 = cellfun(@range, cellfun(@log2,R.expr,'UniformOutput',false) );
% % % % 
% % % %  R = sortrows( R , 'pprange' , 'descend') ;
% % % %  idx  = find(~isnan(R.c) & R.pprange > 5 & ~ismember(R.target_id,od)  )  ;
% % % %  figname = 'all_nondeleted_genes_that_change_loc.eps';
% % % %  delete(figname);
% % % %  G = cellfun(@str2double , regexprep( pp_vars , 'PP_' ,''))  ;
% % % %  clrs = 'kbbrrmgmgcc';
% % % % for I = idx'
% % % %     fh = figure('units','centimeters','position',[5 5 10 7 ]); hold on ;grid on;
% % % %     %plot(R.pp{I},R.expr{I},'ok','MarkerFaceColor',[.7 .7 .7])
% % % %     gscatter( R.pp{I}' , R.expr{I}' , G  , clrs ,'o');
% % % %     grid on;
% % % %     ylabel('Expression (TPM)');
% % % %     xlabel('% peripheral');
% % % %     legend('location','EastOutside')
% % % %     if range(ylim) < 4
% % % %         ylim(  [ min(ylim)-2  max(ylim)+2]);
% % % %     end
% % % %    title( sprintf('%s %s c=%0.02f' , R.target_id{I} , R.gene{I} , corr( R.pp{I}' , R.expr{I}' )) )
% % % %     print('-dpsc2',figname,'-append')
% % % %   
% % % %     close;
% % % %     
% % % % end
% % % %   
% % % %  %% split at %PP == 10
% % % % THRESH = 20 ; 
% % % % idx  = find(  ~isnan(R.c) & R.pprange > THRESH & ~ismember(R.target_id,od)  )  ;
% % % % R.mean_expr_high_pp_10  = NaN(height(R),1);
% % % % R.mean_expr_low_pp_10   = NaN(height(R),1);
% % % % R.ttest2_low_high_pp_10 = NaN(height(R),1);
% % % % R.mean_cidx = NaN( height(R) , 1);
% % % % for I = idx'
% % % %     cidx = kmeans( R.pp{I}'  , 2); 
% % % %     low_pp_cidx = 1 ;
% % % %     if (mean(R.pp{I}(cidx==1)) > mean(R.pp{I}(cidx==2)))
% % % %         low_pp_cidx = 2 ;
% % % %     end
% % % %     exp_pp_low  = R.expr{I}(cidx==low_pp_cidx) ;
% % % %     exp_pp_high = R.expr{I}(cidx~=low_pp_cidx) ;
% % % %     R.mean_expr_high_pp_10(I)  = mean(exp_pp_high);
% % % %     R.mean_expr_low_pp_10(I)  = mean(exp_pp_low);
% % % %     [~,R.ttest2_low_high_pp_10(I)]=ttest2( exp_pp_high , exp_pp_low) ;
% % % %     R.mean_cidx(I) = mean(cidx);
% % % % end
% % % % %
% % % % fh = figure('units','centimeters','position',[5 5 8 8 ]);hold on; 
% % % % histogram( R.c( R.ttest2_low_high_pp_10>0.05) , -1:.1:1)
% % % % histogram( R.c( R.ttest2_low_high_pp_10<0.05) , -1:.1:1)
% % % % legend({'expression does not change' 'expr changes w/%p' },'location','ne')
% % % % xlabel('Correlation between %P and expr')
% % % % ylabel('# of genes')
% % % % a = sum(R.ttest2_low_high_pp_10>0.05 ) / numel(idx) * 100 ;
% % % % b = sum(R.ttest2_low_high_pp_10<0.05 & R.c<0)  / numel(idx) * 100 ;
% % % % c = sum(R.ttest2_low_high_pp_10<0.05 & R.c>0) / numel(idx) * 100 ;
% % % % title(['\Delta' sprintf('%%P>=%d  %0.0f%%,%0.0f%%,%0.0f%%' , THRESH , a,b,c)] )
% % % % 
% % % % %%
% % % % 
