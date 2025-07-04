% Define traits and populations to be processed
traits = {'DiastolicBP', 'MaxHeartRate', 'PulseRate', 'SystolicBP'};
populations = {'AFR', 'EUR'};

% load the variant matrix and get the high confidence SNPs and from the
% data and only SNPs that have a SNP ids 
high_freq_high_quality_variant = readtable( ...
    ['newVariantMetric/',...
    'high_quality_variants_with_high_freq.txt'])  ;

% Loop through each population
for p = 1:length(populations)
    
    % here is the population 
    population = populations{p};

    % Loop through each trait
    for t = 1:length(traits)
        trait = traits{t};
        
        % Construct the base filename
        baseFilename = sprintf('%s_%s', trait, population);

        % Define path to the input file
        inputFile = sprintf('%s.txt', baseFilename);

        % Read the entire table once
        opts = detectImportOptions(inputFile);
        opts.Delimiter = ',';
        dataTable = readtable(inputFile, opts);

        % Calculate Z-scores for the entire dataset
        dataTable.Z = dataTable.Beta ./ dataTable.SE;
        
        % return only the high frequency and high quality variants from the
        % data
        dataTable = dataTable( ismember( dataTable.SNP, ...
           high_freq_high_quality_variant.rsid), :) ;
       
        % Loop through each chromosome (1-22, adjust if necessary)
        parfor chr = 1:22
            
            % print the process
            fprintf('Processing %s - %s - Chromosome %d\n',...
                trait, population, chr);

            % Index rows belonging to the current chromosome
            chrData = dataTable(dataTable.CHR == chr, :);

            % Prepare output table with required columns
            outputTable = chrData(:, {'SNP', 'BP', 'A1', 'A2', 'Z'});
            outputTable.Properties.VariableNames{'SNP'} = 'rsID';
            outputTable.Properties.VariableNames{'BP'} = 'pos';

            % Define output filename for the current chromosome
            outFilename = sprintf('jass_%s/z_UKB_%s_chr%d.txt', ...
                upper(population), upper(trait), chr);

            % Write the chromosome-specific table to a file
            writetable(outputTable, outFilename, 'Delimiter', '\t', ...
                'WriteVariableNames', true);
        end

        fprintf('Data preparation for %s - %s completed.\n', ...
            trait, population);
    end
end

fprintf('All data preparation tasks completed.\n');

