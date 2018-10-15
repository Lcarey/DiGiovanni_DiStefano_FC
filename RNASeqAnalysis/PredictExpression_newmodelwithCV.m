

%% Figure for expression vs % peripheral
load PP
SGD = loadSGDFeatures();
Q = innerjoin(PP,dataset2table(SGD(strcmp(SGD.TYPE,'ORF'),:)),'LeftKey','target_id','RightKey','ORF');

geneproperties = readtable('~/Data/Yeast/protein_properties.tab','FileType','text');
Q = innerjoin( Q, geneproperties,'LeftKey','target_id','RightKey','ORF');
Q.ORF = Q.target_id ;
Q.log_nt_to_closest_end = log10(Q.nt_to_closest_end) ;
Q = Q( : ,     { 'Expr_409'  'CAI' 'CodonBias' 'FOPScore' 'ProteinLength' 'PI' 'nt_to_closest_end' 'PP_409'} ) ;
Q.log2EXPR409 = log2( Q.Expr_409) ;
Q.log10_kb_to_closest_end = log10( Q.nt_to_closest_end ./ 1000 ) ;
clear 'SGD' 'PP' 'geneproperties' ;

%%
T = table();
vn = {'nt_to_closest_end' 'log10_kb_to_closest_end' 'PP_409'} ;
N = 3 ;
for pIvn = 1:numel(vn)
    rvect = NaN(N,1) ;
    for pJ = 1:N
        inputTable = Q;
        predictorNames = {'CAI' 'CodonBias' 'FOPScore' 'ProteinLength' 'PI' vn{pIvn} };
        
        predictors = inputTable(:, predictorNames);
        response = inputTable.log2EXPR409;
        
        % Train a egression model
        % This code specifies all the model options and trains the model.
        template = templateTree(...
            'MinLeafSize', 10);
        regressionEnsemble = fitrensemble(...
            predictors, ...
            response, ...
            'Method', 'Bag', ...
            'NumLearningCycles', 30, ...
            'Learners', template);
        
        % Create the result struct with predict function
        predictorExtractionFcn = @(t) t(:, predictorNames);
        ensemblePredictFcn = @(x) predict(regressionEnsemble, x);
        trainedModel.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));
        trainedModel.RegressionEnsemble = regressionEnsemble;
        
        % Perform cross-validation
        partitionedModel = crossval(trainedModel.RegressionEnsemble, 'KFold', 5);
  
        rvect(pJ) = corr(response, kfoldPredict(partitionedModel))^2  ;
    end
    T.(vn{pIvn}) = rvect
end