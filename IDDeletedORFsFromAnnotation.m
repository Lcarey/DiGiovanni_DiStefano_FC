function [ A , SGD ] = IDDeletedORFsFromAnnotation()
%% [ A , SGD ] = IDDeletedORFsFromAnnotation();
% 
% parse strain annotation tab file and id ORFs that are deleted in FC
% strains
%
% LBC

%%  load SGD features .tab file
SGD = loadSGDFeaturesDS();
SGD = dataset2table( SGD( strcmp(SGD.TYPE,'ORF'),:) ) ;
SGD = SGD( : , {'ORF' 'GENE' 'Chr' 'Start' 'End'} );
SGD.max_pos = max( [SGD.Start SGD.End], [] , 2) ; 
SGD.min_pos = min( [SGD.Start SGD.End], [] , 2) ; 
A = readtable('strain_annotation_regions.tab','FileType','text','Delimiter','\t');

%% fill in gene names
for I = 1:height(SGD)
    if isempty(SGD.GENE{I})
        SGD.GENE{I} = SGD.ORF{I};
    end
end
%% pick out ORFs that are deleted in each strain
A.deleted_orfs = cell( height(A), 1);
A.deleted_genes = cell( height(A), 1);
for I = 2:height(A)
   deleted_regions = regexp( A.regions_deleted{I} , ',' ,'split');
   for dr = deleted_regions
       chr = regexprep( char(dr) , ':.*','') ;
       first_base_del =  str2double( regexprep( regexprep( dr , '.*:','') , '-.*',''));
       last_base_del  =  str2double( regexprep( dr , '.*-','')); 
       idx = find(strcmp(SGD.Chr,chr) & (...
           (SGD.min_pos >= first_base_del & SGD.min_pos <= last_base_del) | ...
           (SGD.max_pos >= first_base_del & SGD.max_pos <= last_base_del) ) ) ; 
       A.deleted_orfs{I} = vertcat( A.deleted_orfs{I} , SGD.ORF(idx) );
       A.deleted_genes{I} = vertcat( A.deleted_genes{I} , SGD.GENE(idx) );
   end
end

%%
