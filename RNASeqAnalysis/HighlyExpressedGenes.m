
%%
load('PP.mat');
PP = PP( PP.Expr_409 , : );
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
    [f,x] = ecdf( data( ismember( PP.target_id , ESR_UP_list)));
    plot(x , f , 'Color',clrs(I,:),'DisplayName', regexprep(vn{exp_idx(I)},'Expr_',''),'LineWidth',2);
end
grid on ;
legend('location','nw')
title(' ESR UP regulated')
line([0 0 ],ylim,'Color','k','LineStyle','--','LineWidth',2)
line(xlim,[0.5 0.5],'Color','k','LineStyle','--','LineWidth',2)
xlim([-1 1])
xlabel('log2( STRAIN / WT )')
ylabel('Fraction of genes')

fh = figure('units','centimeters','position',[5 5 9 9 ]); 
hold on ;
for I = 2:numel(exp_idx)
    expr_data = table2array( PP(:,exp_idx(I)));
    data = log2( (expr_data+1)  ./ (expr_wt+1));
    data(data<-1)=-1 ; data(data>1)=1;
    [f,x] = ecdf( data( ismember( PP.target_id , ESR_DN_list)));
    plot(x , f , 'Color',clrs(I,:),'DisplayName', regexprep(vn{exp_idx(I)},'Expr_',''),'LineWidth',2);
end
grid on ;
legend('location','nw')
title(' ESR DOWN regulated')
line([0 0 ],ylim,'Color','k','LineStyle','--','LineWidth',2)
line(xlim,[0.5 0.5],'Color','k','LineStyle','--','LineWidth',2)

xlim([-1 1])
xlabel('log2( STRAIN / WT )')
ylabel('Fraction of genes')