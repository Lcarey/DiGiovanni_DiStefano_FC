WORKDIR = '/Users/lcarey/Develop/DiGiovanni_DiStefano_FC/Data/ModelPredictions/';
fn = 'genomic_distance_from_telomere_vs_displacement_in_FC_strains_wrt_WT_with_chromosome_position.txt'  ;
% fn = 'genomic_distance_from_closest_chr_arm_end_vs_displacement_in_FC_strains_wrt_WT_with_chromosome_position.txt'  ;

T = readtable( [WORKDIR fn]);
T = sortrows(T,'DisplacementFromNEwrtWT','descend');


T.chr = str2double(regexprep( regexprep(T.Locus,'^chr','') , '_.*','')) ; 
T.chrROM = num2roman(T.chr);


% which are subtelomeric?
T.subtelomeric50 = T.GenDistToTEL < 50*1000 ;
T.subtelomeric100 = T.GenDistToTEL < 100*1000 ;
T.subtelomeric75 = T.GenDistToTEL < 75*1000 ;

%% load CEN positions
SGD = dataset2table( loadSGDFeatures() );
SGD = SGD( strcmp(SGD.TYPE,'centromere') , :);
SGD.Chr = str2double(SGD.Chr);

% identify midpoint of each locus, and dist 2 centromere
T.chrpos = NaN(height(T),1);
for I = 1:height(T)
    metadata = regexp( T.Locus{I} , '_' , 'split');
    T.chrpos(I) = mean(str2double(metadata(3:4))) ; 
    T.KB2CEN(I) = abs( T.chrpos(I) - SGD.Start( SGD.Chr==T.chr(I)) ) ./ 1000 ;
end

% which are peri-centromatic ?
T.pericen25 = T.KB2CEN < 25 ;
T.pericen50 = T.KB2CEN < 50 ;
T.pericen75 = T.KB2CEN < 75 ;
T.pericen100 = T.KB2CEN < 100 ;

%% we don't really care about Donor vs Recipient. Just Fused vs non-fused
idx = find( ~strcmp(T.Fused_NoFused,'NoFused')) ; 
for I = 1:numel(idx)
    T.Fused_NoFused{idx(I)} = 'Fused';
end


 
%%
    %gscatter(  T.GenDistToTEL./1000 , T.DisplacementFromNEwrtWT , { T.Terminal_InternalTEL T.Fused_NoFused} )
figbasename = '~/Downloads/Fig5B_' ; 
T.ID = regexprep( strcat( T.Terminal_InternalTEL, ' and  ' , ' ; ' ,T.Fused_NoFused) , ';' , '') ; 
uids = unique(T.ID);
clrs = [ 0 0 1 ; 0 1 0 ;1 0 0]
for I = 1:numel(uids)

    fh = figure('units','centimeters','position',[5 5 8 8 ]); 
    hold on ;
    plot( T.GenDistToTEL./1000 , T.DisplacementFromNEwrtWT , 'o' ,'Color',[.5 .5 .5]  , 'MarkerFaceColor' , [.7 .7 .7] ) ; 
 
    idx = strcmp(T.ID,uids{I});
    plot( T.GenDistToTEL(idx)./1000 , T.DisplacementFromNEwrtWT(idx) , 'o' ,'Color',[.1 .1 .1] , 'MarkerFaceColor' , clrs(I,:) ) ; 
    
%     plot( T.GenDistToTEL(idx & T.pericen75)./1000 , T.DisplacementFromNEwrtWT(idx& T.pericen75) , '.c' ) ; 
%     plot( T.GenDistToTEL(idx & T.pericen50)./1000 , T.DisplacementFromNEwrtWT(idx& T.pericen50) , '.g' ) ; 
%     plot( T.GenDistToTEL(idx & T.pericen25)./1000 , T.DisplacementFromNEwrtWT(idx& T.pericen25) , '.r' ) ; 

    
xlabel('Distance from the telomere (kb)');
ylabel('Displacement from N.E. (nm)')
set(gca,'xscale','log')
xlim([0 1000])
ylim([-50 250])
%title(uids{I});
print('-dpng', [figbasename num2str(I) '_' uids{I} ] , '-r300');
close;
end

    fh = figure('units','centimeters','position',[5 5 14 8 ]); 
    hold on ;
 %   idx = regexpcmp(T.ID,'d');  
    plot( T.GenDistToTEL./1000 , T.DisplacementFromNEwrtWT , 'o','Color',[.5 .5 .5] , 'MarkerFaceColor' , [.7 .7 .7] ) ; 
 %   plot( T.GenDistToTEL( idx&T.pericen75)./1000 , T.DisplacementFromNEwrtWT(idx&T.pericen75) , '.c' ) ; 
 %   plot( T.GenDistToTEL( idx&T.pericen50)./1000 , T.DisplacementFromNEwrtWT(idx&T.pericen50) , '.g' ) ; 
 %   plot( T.GenDistToTEL( idx&T.pericen25)./1000 , T.DisplacementFromNEwrtWT(idx&T.pericen25) , '.r' ) ; 
    scatter( T.GenDistToTEL(T.pericen100)./1000 , T.DisplacementFromNEwrtWT(T.pericen100) , 15 , T.KB2CEN(T.pericen100) ,'filled')
    
xlabel('Distance from the telomere (kb)');
ylabel('Displacement from N.E. (nm)')
set(gca,'xscale','log')
xlim([0 1000])
ylim([-50 250])
colorbar('Direction','reverse')
print('-dpng', [figbasename 'pcen' ] , '-r300');
close;


%%   fix bugs in code that generated the .txt file
% no longer necessary. file fixed
% % % fix bug that wrongly assigns fused / not fused
% % for I = 1:height(T)
% %     if ~isempty(regexp( T.x_Fusion{I} , ['_' T.chrROM{I} '_'])) %matches
% %         T.Fused_NoFused{I} = 'Fused' ;
% %     end
% % end
% % 
% % % $ grep InternalTEL genomic_distance_from_telomere_vs_displacement_in_FC_strains_wrt_WT_with_chromosome_position.txt|grep NoFused|head
% % % #Fusion Locus Fused/NoFused Arm Telomere Terminal/InternalTEL GenDistToTEL DisplacementFromNEwrtWT
% % % LC_IV_XV_V_VII_cen4 chr5_1_0_9600 NoFused armL TELL5 InternalTEL 4800 207.85
% % % LC_IV_XV_V_VII_cen4 chr5_2_9600_19200 NoFused armL TELL5 InternalTEL 14400 189.965
% % % LC_IV_XV_V_VII_cen4 chr5_3_19200_28800 NoFused armL TELL5 InternalTEL 24000 172.175
% % 
% % % can't be InternalTEL & NoFused. Impossible. 
% % % change all to Donor
% % idx = find(strcmp(T.Terminal_InternalTEL,'InternalTEL') & strcmp(T.Fused_NoFused,'NoFused')) ; 
% % for I = 1:numel(idx)
% %     T.Fused_NoFused{idx(I)} = 'Donor';
% % end
