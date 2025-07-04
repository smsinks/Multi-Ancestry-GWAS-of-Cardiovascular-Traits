%% SNPs Associated w/ Heart Function Among AFR and EUR

clc; clear; close all force

try
    cd '/Users/sinkala/Documents/MATLAB/UKBioBank/cardiovascularNewAnalysis'
    
    % add the path to the GWAS results for Africans and Europeans
    addpath(['/Users/sinkala/Documents/MATLAB/UKBioBank/' ...
        'cardiovascularNewAnalysis/eur_afr_gwas'])
    
    onHPC_Cluster = false;
catch
    addpath('/scratch/snkmus003/ukbiobank')
    cd '/scratch/snkmus003/ukbiobank/heartPaper'
   
    % specificy that we are working on the cluster
    onHPC_Cluster = true;
    numOfCores = 10;
end

%% Run the analysis to collect the data from the LDLink 

warning('off')
if ~exist('LDLinkData','dir')
    mkdir('LDLinkData')
end
addpath('/cardiovascularNewAnalysis/LDLinkData')

% process the data 
if ~exist('gwasCatalog_heartAFR.txt','file')

    % laod the AFR data if it exists
    curAFRdata = readtable(['FUMA_multitrait_cardio_AFR_results/', ...
        'IndSigSNPs.txt'] ) ;

    % here are all the variables and an empty table
    allVariants = curAFRdata.rsID ;

    gwasCat_AFR = [] ;
    % ****************************************************************

    % run a parrel for loop here
    for kk = 1:length(allVariants)

        fprintf(['\nDownloading LD data from LDlink.com ',  ...
            'for %s #%d of %d\n'],  allVariants{kk}, kk, ...
            length(allVariants) )

        % check that the file already exists
        filename = [allVariants{kk},'.txt'] ;

        % here are the curl options
        curl_URL = sprintf(['curl -k -H "Content-Type: application/json" '...
            '-X POST -d ''{"snps": "%s", "pop": "AFR", "r2_d": "r2", ' ...
            '"r2_d_threshold": "0.1", "window": "500000", ' ...
            '"genome_build": "grch37"}'' ''https://ldlink.nih.gov/' ...
            'LDlinkRest/ldtrait?token=97fb8311b229'' -o %s.txt'], ...
            allVariants{kk}, allVariants{kk});

        % now download the data
        if ~exist(filename,'file')
            try
                system(curl_URL)
            catch
                fprintf('\nFailed to DOWNLOAD file for %s\n', ...
                    allVariants{kk} )
                continue ;
            end
        end

        % Load the data into a MATLAB table
        try
            curLinkData = readtable(filename,'ReadVariableNames',true) ;
            
            % check if the file has beeen in correctly
            if ~any(contains(curLinkData.Properties.VariableNames,  ...
                    'NoEntriesInTheGWASCatalogAreIdentified' ))
                if ~any(contains(curLinkData.Properties.VariableNames, ...
                        'GWASTrait'))
                    % read the file in a different way
                    curLinkData = readtable(filename,'Delimiter', '\t') ;
                end
            end
            
        catch
            fprintf('\nFailed to READ file for %s\n', allVariants{kk})
            continue
        end

        % remove the file after reading i
        system(['mv ',filename,' LDLinkData/'])

        % check that data are okay
        if any(contains(curLinkData.Properties.VariableNames,  ...
                'NoEntriesInTheGWASCatalogAreIdentified' ))
            curAFRdata.Novelty(kk) = {'Novel'} ;
            curAFRdata.theTraits(kk) = {'None'} ;
        elseif any(contains(curLinkData.Properties.VariableNames,  ...
                ' GWASTrait'))
            curAFRdata.Novelty(kk) = {'Has Trait'} ;
            curAFRdata.theTraits(kk) = cellsr( strjoin( ...
                curLinkData.GWASTrait,';') ) ;
        end

        if any(contains(curLinkData.Properties.VariableNames,'error'))
            continue
        end

        % clean up the variable names and the data
        curLinkData.Position = extractAfter( ...
            curLinkData.Position_GRCh37_ ,':');
        curLinkData.Chrom = extractAfter( extractBefore( ...
            curLinkData.Position_GRCh37_ ,':') , 'chr') ;
        curLinkData.Position = str2double(curLinkData.Position);
        curLinkData.Chrom = str2double( curLinkData.Chrom) ;
        curLinkData.Position_GRCh37_ = [] ;
        
        % convert the risk allele to cell str
        if iscell( curLinkData.RiskAllele)
            curLinkData.RiskAllele = str2double(curLinkData.RiskAllele) ;
        end
        
        if iscell( curLinkData.EffectSize_95_CI_)
            curLinkData.EffectSize_95_CI_ = ...
                str2double(curLinkData.EffectSize_95_CI_) ;
        end

        % remove the rows where the index SNP and the Linked SNP
        % remove the data with no linked snp
        gwasCat_AFR = [gwasCat_AFR; curLinkData] ;

    end

    % save the data to a file
    writetable(gwasCat_AFR,'gwasCatalog_heartAFR.txt');
    writetable(curAFRdata,'gwasCatalog_AFR_Trait_IndSNPs.txt')
    afrCatalog = curAFRdata; 
else
    gwasCat_AFR = readtable('gwasCatalog_heartAFR.txt')
    afrCatalog = readtable('gwasCatalog_AFR_Trait_IndSNPs.txt')
end

% ********************************************************************
% ******************* Do the same for the EUR groups  ****************
if ~exist('gwasCatalog_heartEUR.txt','file')

    % laod the AFR data if it exists
    curEURdata = readtable(['FUMA_multitrait_cardio_EUR_results/', ...
        'IndSigSNPs.txt'] ) ;

    % here are all the variables and an empty table
    allVariants = curEURdata.rsID ;

    gwasCat_EUR = [] ;
    % ****************************************************************

    % run a parrel for loop here
    for kk = 1:length(allVariants)

        fprintf(['\nDownloading LD data from LDlink.com ',  ...
            'for %s #%d of %d\n'],  allVariants{kk}, kk, ...
            length(allVariants) )

        % check that the file already exists
        filename = [allVariants{kk},'.txt'] ;

        % here are the curl options
        curl_URL = sprintf(['curl -k -H "Content-Type: application/json" '...
            '-X POST -d ''{"snps": "%s", "pop": "EUR", "r2_d": "r2", ' ...
            '"r2_d_threshold": "0.1", "window": "500000", ' ...
            '"genome_build": "grch37"}'' ''https://ldlink.nih.gov/' ...
            'LDlinkRest/ldtrait?token=97fb8311b229'' -o %s.txt'], ...
            allVariants{kk}, allVariants{kk});

        % now download the data
        if ~exist(filename,'file')
            try
                system(curl_URL)
            catch
                fprintf('\nFailed to DOWNLOAD file for %s\n', ...
                    allVariants{kk} )
                continue ;
            end
        end

        % Load the data into a MATLAB table
        try
            curLinkData = readtable(filename,'ReadVariableNames',true) ;
            
            % check if the file has beeen in correctly
            if ~any(contains(curLinkData.Properties.VariableNames,  ...
                    'NoEntriesInTheGWASCatalogAreIdentified' ))
                if ~any(contains(curLinkData.Properties.VariableNames, ...
                        'GWASTrait'))
                    % read the file in a different way
                    curLinkData = readtable(filename,'Delimiter', '\t') ;
                end
            end
        catch
            fprintf('\nFailed to READ file for %s\n', allVariants{kk})
            continue
        end

        % remove the file after reading i
        system(['mv ',filename,' LDLinkData/'])

        % check that data are okay
        if any(contains(curLinkData.Properties.VariableNames,  ...
                'NoEntriesInTheGWASCatalogAreIdentified' ))
            curEURdata.Novelty(kk) = {'Novel'} ;
            curEURdata.theTraits(kk) = {'None'} ;
        elseif any(contains(curLinkData.Properties.VariableNames,  ...
                ' GWASTrait'))
            curEURdata.Novelty(kk) = {'Has Trait'} ;
            curEURdata.theTraits(kk) = cellsr( strjoin( ...
                curLinkData.GWASTrait,';') ) ;
        end

        % check that data are okay
        if any(contains(curLinkData.Properties.VariableNames,'error')) ...
                || any(contains(curLinkData.(1){1},'error'))
            continue
        end

        % clean up the variable names and the data
        curLinkData.Position = extractAfter( ...
            curLinkData.Position_GRCh37_ ,':');
        curLinkData.Chrom = extractAfter( extractBefore( ...
            curLinkData.Position_GRCh37_ ,':') , 'chr') ;
        curLinkData.Position = str2double(curLinkData.Position);
        curLinkData.Chrom = str2double( curLinkData.Chrom) ;
        curLinkData.Position_GRCh37_ = [] ;
        
        % convert the risk allele to cell str
        if iscell( curLinkData.RiskAllele)
            curLinkData.RiskAllele = str2double(curLinkData.RiskAllele) ;
        end
        
        if iscell( curLinkData.EffectSize_95_CI_)
            curLinkData.EffectSize_95_CI_ = ...
                str2double(curLinkData.EffectSize_95_CI_) ;
        end

        % remove the rows where the index SNP and the Linked SNP
        % remove the data with no linked snp
        gwasCat_EUR = [gwasCat_EUR; curLinkData] ;

    end

    % save the data to a file
    writetable(gwasCat_EUR,'gwasCatalog_heartEUR.txt');
    writetable(curEURdata,'gwasCatalog_EUR_Trait_IndSNPs.txt')
    eurCatalog = curEURdata ;
else
    gwasCat_EUR = readtable('gwasCatalog_heartEUR.txt')
    eurCatalog = readtable('gwasCatalog_EUR_Trait_IndSNPs.txt')
end