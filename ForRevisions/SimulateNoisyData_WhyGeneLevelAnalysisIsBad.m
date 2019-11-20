%% simulate measurement error in which the trend is detected perfectly, 
% but no individual point is significantly different from 0
% 
% Why analysis by position or predicted location is more powerful than single-gene analysis. 
%     For each of positions (~genes) along X the real value (grey circles) is y=0.1*x. 
%     We then simulated measurements (~expression data, black points) with normally 
%        distributed measurement noise. Only two positions (green asterisk) have measured 
%       values statistically significant from zero. However, fitting an equation (red +) 
%       to all the measurements shows the true value at all positions. 
% 
% LBC November 2019

FIGNAME = '~/Downloads/example.png' ; 
delete(FIGNAME) ; 
for J = 1:100

X = 0:10 ; 
Y = X * 0.1 ; 

Nsamples = 10 ; 
Y_Noise = random( 'normal' , 0 , 0.75 , Nsamples , numel(X) ) ;
X_Noise = random( 'normal' , 0 , 0.1 , Nsamples , numel(X) ) ;
p = NaN( numel(X),1);
c = 0 ; 
fh = figure('units','centimeters','position',[50 50 6 4]);
hold on ; 
line( [-1 12] , [0 0 ], 'LineStyle' ,'--' ,'Color',[.8 .8 .8] )

plot( X , Y , 'o' , 'Color', [.5 .5 .5] ,'MarkerFaceColor',[.5 .5 .5] ,'DisplayName' , 'True Value')
all_data = NaN(0,2);
    for I = 1:numel(X)
    y = Y(I) + Y_Noise(:,I) ; ;
    x = X(I) + X_Noise(:,I) ;
    plot( x , y , '.k')
    all_data = vertcat( all_data , [ x  y ] );
    [~,p(I)] = ttest( y );
    
  %  if (p(I)) <0.05 
  %      text( X(I) , 2 , '*' ,'Color','b')
  %  end
    if (p(I)*numel(X)) <0.05 
        c = c + 1; 
        text( X(I) , 2 , '*' ,'Color','g')
    end
    
end

ylim( [ -0.9 2.2 ] )
xlim( [-0.5 10.5])
set(gca,'xtick',0:10)
%grid on ; 
%set(gca,'ytick',-1:10)


% Set up fittype and options.
ft = fittype( 'poly1' );
[fitresult, gof] = fit( all_data(:,1), all_data(:,2), ft );
%plot(fitresult)
legend('off')
plot( X , feval(fitresult,X) , '+r');

if mean(abs(Y - feval(fitresult,X)')) < 0.03 & c < 3 
    print('-dpng',FIGNAME,'-r600');
end

close ;

end