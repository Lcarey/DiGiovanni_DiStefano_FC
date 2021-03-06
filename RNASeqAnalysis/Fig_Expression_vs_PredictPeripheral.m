%% Figure for expression vs % peripheral
load PP
SGD = loadSGDFeatures();
Q = innerjoin(PP,dataset2table(SGD(strcmp(SGD.TYPE,'ORF'),:)),'LeftKey','target_id','RightKey','ORF');

geneproperties = readtable('~/Data/Yeast/protein_properties.tab','FileType','text');
Q = innerjoin( Q, geneproperties,'LeftKey','target_id','RightKey','ORF');
Q.ORF = Q.target_id ; 
Q.log_nt_to_closest_end = log10(Q.nt_to_closest_end) ; 
%%
figname = 'PredictExpression.eps' ;
delete(figname);
fh = figure('units','centimeters','position',[5 5 7 7 ]);
sh = scatter( PP.PP_409 , PP.Expr_409 , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
set(gca,'yscale','log')
ylabel('Expression (TPM)')
xlabel('% peripheral (predicted)')
axis tight;
[cP,~] = corr(  PP.PP_409 , log10(PP.Expr_409) ,'type','pearson')
[cS,~] = corr(  PP.PP_409 , log10(PP.Expr_409) ,'type','spearman')
[cK,~] = corr(  PP.PP_409 , log10(PP.Expr_409) ,'type','kendall')
title( sprintf('Pearson corr = %0.04f' , cP ))
print('-dpsc2',figname,'-append');
close;

fh = figure('units','centimeters','position',[5 5 7 7 ]);
sh = scatter( log10(PP.nt_to_closest_end ./ 1000) , PP.Expr_409 , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
set(gca,'yscale','log')
ylabel('Expression (TPM)')
xlabel('kb to telomere')
axis tight;
[c,p] = corr(  log10(PP.nt_to_closest_end ./ 1000)  , log10(PP.Expr_409) );
set(gca,'xdir','reverse')
set(gca,'xtick',[-1 0 1 2]); set(gca,'xticklabels',[1 10 100 1000]);
title( sprintf('Pearson corr = %0.04f' , c ))
print('-dpsc2',figname,'-append');
close;
%%
%X = zscore([ Q.CAI  Q.CodonBias Q.FOPScore   Q.ProteinLength Q.PI Q.PP_409 log(Q.nt_to_closest_end) ]) ; 
%Y = (Q.Expr_409)  ; 
K = 10 
C = cvpartition( height(Q) ,'KFold',K) ; 
ydist = 'poisson' ;
for cvI = 1:K
    idx_train = training(C,cvI); 
    idx_test  = test(C,cvI); 
    mdlDISTlin = fitglm( Q, 'PredictorVars',{'CAI' 'CodonBias' 'FOPScore' 'ProteinLength' 'PI' 'nt_to_closest_end'} , 'ResponseVar','Expr_409');
    
    
    
    
ydist = 'poisson' ;
mdlstr = 'y ~ x1 + x2 + x3 + x6 +x7' ;
mdl = fitglm(X,Y,mdlstr ,'Distribution',ydist);
mdl1 = fitglm(X,Y,[ mdlstr ' + x4 '] ,'Distribution',ydist);
mdlNT = addTerms(mdl , 'x5');
mdlPP = addTerms(mdl , 'x4');

t1 = prctile(Y,75) ; %t1 = 100 ; 
t2 = prctile(Y,50) ; %t2 = 50 ; 
t3 = prctile(Y,25) ; %t3 = 10 ; 

mdlLOW = fitglm(X(Y>t1,:),Y(Y>t1),mdlstr,'Distribution',ydist);
mdlNTLOW = addTerms(mdlLOW , 'x5');
mdlPPLOW = addTerms(mdlLOW , 'x4');

mdlVLOW = fitglm(X(Y<t1 & Y>t3,:),Y(Y<t1 & Y>t3),mdlstr,'Distribution',ydist);
mdlNTVLOW = addTerms(mdlVLOW , 'x5');
mdlPPVLOW = addTerms(mdlVLOW , 'x4');

mdlVVLOW = fitglm(X(Y<t3,:),Y(Y<t3),mdlstr,'Distribution',ydist);
mdlNTVVLOW = addTerms(mdlVVLOW , 'x5');
mdlPPVVLOW = addTerms(mdlVVLOW , 'x4');


%mdlVVVLOW = fitglm(X(Y<t4,:),Y(Y<t4),mdlstr,'Distribution',ydist);
%mdlNTVVVLOW = addTerms(mdlVVVLOW , 'x5');
%mdlPPVVVLOW = addTerms(mdlVVVLOW , 'x4');

%fprintf('%0.04f\t%0.04f\t%0.04f\n' , mdl.Rsquared.Ordinary , mdlPP.Rsquared.Ordinary , mdl1.Rsquared.Ordinary );
% %
% fh = figure('units','centimeters','position',[5 5 7 7 ]);
% sh = scatter( mdlPP.predict , PP.Expr_409 , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% ylabel('Expression (TPM)')
% xlabel('Predicted expression (TPM)')
% axis tight;
% set(gca,'xtick',logspace(-1,6,8))
% set(gca,'ytick',logspace(-1,6,8))
% xlim([1 2e4])
% ylim(xlim)
% title( sprintf('r^2 = %0.05f' , mdlPP.Rsquared.Ordinary) ) 
%
data = [ mdlPP.Rsquared.Ordinary-mdl.Rsquared.Ordinary mdlNT.Rsquared.Ordinary-mdl.Rsquared.Ordinary ; ...
        mdlPPLOW.Rsquared.Ordinary-mdlLOW.Rsquared.Ordinary mdlNTLOW.Rsquared.Ordinary-mdlLOW.Rsquared.Ordinary ; ...
        mdlPPVLOW.Rsquared.Ordinary-mdlVLOW.Rsquared.Ordinary mdlNTVLOW.Rsquared.Ordinary-mdlVLOW.Rsquared.Ordinary ; ...
        mdlPPVVLOW.Rsquared.Ordinary-mdlVVLOW.Rsquared.Ordinary mdlNTVVLOW.Rsquared.Ordinary-mdlVVLOW.Rsquared.Ordinary ; ...
    ];

fh = figure('units','centimeters','position',[5 5 7 7 ]);
%bar([2 1] , data(1:2,:));
bar( [4 3 2 1] , data(1:4,:))
ylim([0 max(ylim)])
ylabel('Increase in r^2')
%set(gca,'xticklabel',{'<100 TPM' 'all genes' })
%set(gca,'xticklabel',{'all' '<100 TPM' '<50' '<10 TPM'})
legend({'% peripheral' 'kb to telomere'})
set(gca,'xtick',[])
%% compare + / - model w/low expr
%%
fh = figure('units','centimeters','position',[5 5 7 7 ]);
sh = scatter( mdlLOW.predict , mdlLOW.Variables.y , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
set(gca,'yscale','log')
set(gca,'xscale','log')
ylabel('Expression (TPM)')
xlabel('Predicted expression (TPM)')
axis tight;
set(gca,'xtick',[0 1 2 5 10 15 20 30 40])
set(gca,'ytick',[0 1 2 5 10 15 20 30 40])
xlim([2 75])
ylim(xlim)
title( sprintf('r^2 = %0.05f' , mdlLOW.Rsquared.Ordinary) ) 

fh = figure('units','centimeters','position',[5 5 7 7 ]);
sh = scatter( mdlPPLOW.predict , mdlPPLOW.Variables.y , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
set(gca,'yscale','log')
set(gca,'xscale','log')
ylabel('Expression (TPM)')
xlabel('Predicted expression (TPM)')
axis tight;
set(gca,'xtick',[0 1 2 5 10 15 20 30 40])
set(gca,'ytick',[0 1 2 5 10 15 20 30 40])
xlim([2 75])
ylim(xlim)
title( sprintf('r^2 = %0.05f' , mdlPPLOW.Rsquared.Ordinary) ) 

%% % peripheral vs nt from end
%%
fh = figure('units','centimeters','position',[5 5 7 7 ]);
sh = scatter( PP.nt_to_closest_end./1000 , PP.PP_409  , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
set(gca,'yscale','lin')
set(gca,'xscale','log')
axis tight;
ylim([0 100])
ylabel('% peripheral')
xlabel('kb from telomere')
set(gca,'xtick',logspace(-1,6,8))
set(gca,'xticklabel',[.1 1 10 100 1000])
set(gca,'ytick',0:25:100)
title( sprintf('Pearson corr = %0.05f' , corr(PP.nt_to_closest_end./1000 , PP.PP_409) )) 
%%