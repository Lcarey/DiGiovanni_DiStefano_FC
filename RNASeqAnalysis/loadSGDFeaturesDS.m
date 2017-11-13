function SGDFeatures = loadSGDFeaturesDS(SGDFeatureFile)
%% SGDFeatures = loadSGDFeaturesDS(SGDFeatureFile)
% Find or download SGD_Features.tab
% LBC 

URL = 'https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab' ;
SystematicNameField = 'WellContentAlias';
SGDTableFieldName = 'ORF';
SGDVarNames = {'SGDID','TYPE','VERIFIED','ORF','GENE','ALIAS','Chromosome','loc','Chr','Start','End','Strand','GeneticPosCM','moddate',    'adddate','description'};


% default SGD Features file
if ~exist('SGDFeatureFile','var')
    SGDFeatureFile = '~/Data/Gene/Description/Yeast/Sgd/SGD_features.tab'   ;
end

% look for file
if exist(SGDFeatureFile,'file')
    fprintf('Found  %s , loading...\n' , SGDFeatureFile );
    SGDFeatures = dataset('file',SGDFeatureFile,'ReadVarNames',true);
elseif exist('SGD_features.tab','file' ) % check current working dir
    disp(['Found SGD_features.tab, loading...']);
    SGDFeatures = dataset('file','SGD_features.tab','ReadVarNames','off');
else
    disp('downloading SGD Features File')
    a = urlwrite(URL,'SGD_features.tab');
    if (~ a)
        error('SGD_features.tab doesnt exist and couldnt download.');
    end
    SGDFeatures = dataset('file','SGD_features.tab','ReadVarNames',true);
    if (isempty(SGDFeatures) )
        error('Could download but couldnt read SGD_features.tab');
    end
end
SGDFeatures = set(SGDFeatures,'VarNames',SGDVarNames);

end
