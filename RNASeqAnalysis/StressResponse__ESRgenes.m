%% load data
BB = readtable('~/Data/Brauer08/TableS1.xls');
ESR_UP_list = BB.ORF( strcmp(BB.ESR,'up'));
ESR_DN_list = BB.ORF( strcmp(BB.ESR,'down'));
GR_UP_list = BB.ORF(strcmp(BB.GrowthRateResponse_1_5SD,'up'));
GR_DN_list = BB.ORF(strcmp(BB.GrowthRateResponse_1_5SD,'down'));


load('PP.mat');
%PP = PP( PP.nt_to_closest_end > 5e4 , : );
% PP = PP( PP.nt_to_closest_end > 5e4 , : ); % remove subtelomeric

vn = PP.Properties.VariableNames ; 
exp_idx = find(regexpcmp( vn , 'Expr_') ); 
clrs = parula(numel(exp_idx)) ;
expr_wt = PP.Expr_409 ; 
%%
fh = figure('units','centimeters','position',[5 5 9 9 ]); 
hold on ;
for I = 2:numel(exp_idx)
    expr_data = table2array( PP(:,exp_idx(I)));
    data = log2( (expr_data+1)  ./ (expr_wt+1));
    data(data<-1)=-1 ; data(data>1)=1;
    [f,x] = ecdf( data( ismember( PP.target_id , GR_UP_list)));
    plot(x , f , 'Color',clrs(I,:),'DisplayName', regexprep(vn{exp_idx(I)},'Expr_',''),'LineWidth',2);
end
grid on ;
legend('location','nw')
title('GR UP regulated')
lh = line([0 0 ],ylim,'Color','k','LineStyle','--','LineWidth',2,'HandleVisibility','off');
lh = line(xlim,[0.5 0.5],'Color','k','LineStyle','--','LineWidth',2,'HandleVisibility','off');
xlim([-1 1])
xlabel('log2( FC_{strain} / WT )')
ylabel('Fraction of genes')

fh = figure('units','centimeters','position',[5 5 9 9 ]); 
hold on ;
for I = 2:numel(exp_idx)
    expr_data = table2array( PP(:,exp_idx(I)));
    data = log2( (expr_data+1)  ./ (expr_wt+1));
    data(data<-1)=-1 ; data(data>1)=1;
    [f,x] = ecdf( data( ismember( PP.target_id , GR_DN_list)));
    plot(x , f , 'Color',clrs(I,:),'DisplayName', regexprep(vn{exp_idx(I)},'Expr_',''),'LineWidth',2);
end
grid on ;
legend('location','nw')
title('GR DOWN regulated')
lh = line([0 0 ],ylim,'Color','k','LineStyle','--','LineWidth',2,'HandleVisibility','off');
lh = line(xlim,[0.5 0.5],'Color','k','LineStyle','--','LineWidth',2,'HandleVisibility','off');

xlim([-1 1])
xlabel('log2( FC_{strain} / WT )')
ylabel('Fraction of genes')