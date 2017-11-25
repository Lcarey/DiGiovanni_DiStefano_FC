%% Figure for expression vs % peripheral
load PP
%%
geneproperties = readtable('~/Google Drive/CareyLab/ExternalData/Yeast/protein_properties.tab','FileType','text');
%%
fh = figure('units','centimeters','position',[5 5 7 7 ]);
sh = scatter( PP.PP_409 , PP.Expr_409 , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
set(gca,'yscale','log')
ylabel('Expression dTPM)')
xlabel('% peripheral (predicted)')
axis tight;
[c,p] = corr(  PP.PP_409 , log10(PP.Expr_409) )

%%
fh = figure('units','centimeters','position',[5 5 7 7 ]);
sh = scatter( log10(PP.nt_to_closest_end ./ 1000) , PP.Expr_409 , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
set(gca,'yscale','log')
ylabel('Expression dTPM)')
xlabel('kb to telomere (log10)')
axis tight;
[c,p] = corr(  log10(PP.nt_to_closest_end ./ 1000)  , log10(PP.Expr_409) )

%%
X = zscore([ Q.CAI  Q.CodonBias Q.FOPScore  Q.PP_409 Q.nt_to_closest_end Q.ProteinLength Q.PI ]) ; 
Y = (Q.Expr_409)  ; 
ydist = 'poisson' ;
mdlstr = 'y ~ x1 + x2 + x3 + x6 +x7' ;
mdl = fitglm(X,Y,mdlstr ,'Distribution',ydist);
mdl1 = fitglm(X,Y,[ mdlstr ' + x4 '] ,'Distribution',ydist);
mdlNT = addTerms(mdl , 'x5');
mdlPP = addTerms(mdl , 'x4');

t1 = prctile(Y,75) ; t1 = 100 ; 
t2 = prctile(Y,50) ; t2 = 50 ; 
t3 = prctile(Y,25) ; t3 = 10 ; 

mdlLOW = fitglm(X(Y<t1,:),Y(Y<t1),mdlstr,'Distribution',ydist);
mdlNTLOW = addTerms(mdlLOW , 'x5');
mdlPPLOW = addTerms(mdlLOW , 'x4');

mdlVLOW = fitglm(X(Y<t2,:),Y(Y<t2),mdlstr,'Distribution',ydist);
mdlNTVLOW = addTerms(mdlVLOW , 'x5');
mdlPPVLOW = addTerms(mdlVLOW , 'x4');

mdlVVLOW = fitglm(X(Y<t3,:),Y(Y<t3),mdlstr,'Distribution',ydist);
mdlNTVVLOW = addTerms(mdlVVLOW , 'x5');
mdlPPVVLOW = addTerms(mdlVVLOW , 'x4');


mdlVVVLOW = fitglm(X(Y<t4,:),Y(Y<t4),mdlstr,'Distribution',ydist);
mdlNTVVVLOW = addTerms(mdlVVVLOW , 'x5');
mdlPPVVVLOW = addTerms(mdlVVVLOW , 'x4');

fprintf('%0.04f\t%0.04f\t%0.04f\n' , mdl.Rsquared.Ordinary , mdlPP.Rsquared.Ordinary , mdl1.Rsquared.Ordinary );
%%
fh = figure('units','centimeters','position',[5 5 7 7 ]);
sh = scatter( mdlPP.predict , PP.Expr_409 , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
set(gca,'yscale','log')
set(gca,'xscale','log')
ylabel('Expression (TPM)')
xlabel('Predicted expression (TPM)')
axis tight;
set(gca,'xtick',logspace(-1,6,8))
set(gca,'ytick',logspace(-1,6,8))
xlim([1 2e4])
ylim(xlim)
title( sprintf('r^2 = %0.05f' , mdlPP.Rsquared.Ordinary) ) 
%%
data = [ mdlPP.Rsquared.Ordinary-mdl.Rsquared.Ordinary mdlNT.Rsquared.Ordinary-mdl.Rsquared.Ordinary ; ...
        mdlPPLOW.Rsquared.Ordinary-mdlLOW.Rsquared.Ordinary mdlNTLOW.Rsquared.Ordinary-mdlLOW.Rsquared.Ordinary ; ...
        mdlPPVLOW.Rsquared.Ordinary-mdlVLOW.Rsquared.Ordinary mdlNTVLOW.Rsquared.Ordinary-mdlVLOW.Rsquared.Ordinary ; ...
        mdlPPVVLOW.Rsquared.Ordinary-mdlVVLOW.Rsquared.Ordinary mdlNTVVLOW.Rsquared.Ordinary-mdlVVLOW.Rsquared.Ordinary ; ...
    ];

fh = figure('units','centimeters','position',[5 5 7 7 ]);
bar([2 1] , data(1:2,:))
ylim([0 max(ylim)])
ylabel('Increase in r^2')
set(gca,'xticklabel',{'<100 TPM' 'all genes' })
%set(gca,'xticklabel',{'all' '<100 TPM' '<50' '<10 TPM'})
legend({'% peripheral' 'kb to telomere'})

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
title( sprintf('r^2 = %0.05f' , corr(PP.nt_to_closest_end./1000 , PP.PP_409) )) 
%%