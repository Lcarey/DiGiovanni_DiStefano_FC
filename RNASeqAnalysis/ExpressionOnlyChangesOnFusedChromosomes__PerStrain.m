fh = figure('units','centimeters','position',[5 5 5 20 ]);
R = 6 ; 
C = 2 ; 

idx = PP.chr == 4 | PP.chr == 15 |PP.chr == 16 ; 
subplot(R,C,1) ; hold on ; 
ecdf(PP.Pdiff_1388(idx))
%ecdf(PP.Pdiff_1388)
ecdf(PP.Pdiff_1388(~idx))
subplot(R,C,2) ; hold on ; 
ecdf(PP.Efc_1138(idx))
%ecdf(PP.Efc_1138)
ecdf(PP.Efc_1138(~idx))

idx = PP.chr == 4 | PP.chr == 15 |PP.chr == 15 |PP.chr == 7 ; 
subplot(R,C,3) ; hold on ; 
ecdf(PP.Pdiff_1788(idx))
%ecdf(PP.Pdiff_1788)
ecdf(PP.Pdiff_1788(~idx))
subplot(R,C,4) ; hold on ; 
ecdf(PP.Efc_1788(idx))
%ecdf(PP.Efc_1788)
ecdf(PP.Efc_1788(~idx))

idx = PP.chr == 4 | PP.chr == 15 |PP.chr == 5 ; 
subplot(R,C,5) ; hold on ; 
ecdf(PP.Pdiff_1379(idx))
%ecdf(PP.Pdiff_1379)
ecdf(PP.Pdiff_1379(~idx))
subplot(R,C,6) ; hold on ; 
ecdf(PP.Efc_1379(idx))
%ecdf(PP.Efc_1379)
ecdf(PP.Efc_1379(~idx))

idx = PP.chr == 4 | PP.chr == 12 ; 
subplot(R,C,7) ; hold on ; 
ecdf(PP.Pdiff_524(idx))
%ecdf(PP.Pdiff_524)
ecdf(PP.Pdiff_524(~idx))
subplot(R,C,8) ; hold on ; 
ecdf(PP.Efc_524(idx))
%ecdf(PP.Efc_524)
ecdf(PP.Efc_524(~idx))

idx = PP.chr == 4 | PP.chr == 15 | PP.chr == 5  ; 
subplot(R,C,9) ; hold on ; 
ecdf(PP.Pdiff_1387(idx))
%ecdf(PP.Pdiff_1387)
ecdf(PP.Pdiff_1387(~idx))
subplot(R,C,10) ; hold on ; 
ecdf(PP.Efc_1387(idx))
%ecdf(PP.Efc_1387)
ecdf(PP.Efc_1387(~idx))

clrs = vertcat( [ 0 0 0 ; winter(5) ] );
subplot(R,C,11) ; hold on ; 
for I = 1:5
    [f,x] = ecdf( random('normal' , 0 , (I-0.9) , 1e4 , 1));
    plot(x,f,'-','Color',clrs(I,:))
end
subplot(R,C,12) ; hold on ; 
for I = 1:5
    [f,x] = ecdf( random('normal' , 0 , (I-0.9)./200 , 1e4 , 1));
    plot(x,f,'-','Color',clrs(I,:))
end


for I = 1:R*C
    subplot(R,C,I);
    ylabel('');
    xlabel('')
    set(gca,'ytick',[0 1])
    xlim([-0.025 0.025])
end
for I = 1:2:R*C
    subplot(R,C,I)
    xlim([-2 5])
end

%% try w/bar plots
fh = figure('units','centimeters','position',[5 5 5 20 ]);
T = 5 ;
E = -0.1 ; 
idx = PP.chr == 4 | PP.chr == 15 |PP.chr == 16 ; 
subplot(R,1,1) ; hold on ; 
data = [ mean((PP.Pdiff_1388(idx)) > T)  mean((PP.Pdiff_1388(~idx)) > T) ;...
    mean((PP.Efc_1138(idx))<E)  mean((PP.Efc_1138(~idx))<E) ] ; 
bar(data)

subplot(R,1,2)
idx = PP.chr == 4 | PP.chr == 15 |PP.chr == 15 |PP.chr == 7 ; 
data = [ mean((PP.Pdiff_1788(idx)) > T)  mean((PP.Pdiff_1788(~idx)) > T) ;...
    mean((PP.Efc_1788(idx))<E)  mean((PP.Efc_1788(~idx))<E) ] ; 
bar(data)

idx = PP.chr == 4 | PP.chr == 15 |PP.chr == 5 ; 
subplot(R,1,3) ; hold on ; 
data = [ mean((PP.Pdiff_1379(idx)) > T)  mean((PP.Pdiff_1379(~idx)) > T) ;...
    mean((PP.Efc_1379(idx))<E)  mean((PP.Efc_1379(~idx))<E) ] ; 
bar(data)

