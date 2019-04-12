WORKDIR = '/Users/lcarey/Develop/DiGiovanni_DiStefano_FC/Data/ModelPredictions/';

T = readtable( [WORKDIR 'scatter_plot_distance_from_chr_arm_end_vs_displacement_in_FC.txt']);
%% Figure 5. Loci near the end of fused chromosomes are predicted to be displaced away from the nuclear periphery. (A) The predicted displacement with respect to the nuclear envelope for loci in fused (blue) and non-fused (orange) chromosomes. (B) The predicted displacement with respect to the nuclear envelope for loci as a function of distance from the end of the chromosome in wild-type cells. Loci that are closer to the ends of the chromosomes exhibit a greater change away from the nuclear envelope. 
figure ; 
plot(T.Var1./1000,T.Var2,'ok')
set(gca,'xscale','log')
xlim([0 600])
ylim([-50 250])