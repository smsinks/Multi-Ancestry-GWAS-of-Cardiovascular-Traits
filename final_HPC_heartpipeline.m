%% SNPs Associated w/ Heart Function Among AFR and EUR

clc; clear; close all force

try
    cd '/Users/sinkala/Documents/MATLAB/UKBioBank/cardiovascularNewAnalysis'
    
    % add the path to the GWAS results for Africans and Europeans
    addpath(['/Users/sinkala/Documents/MATLAB/UKBioBank/' ...
        'Manuscript Heart/eur_afr_gwas'])
    
    onHPC_Cluster = false;
catch
    addpath('/scratch/snkmus003/ukbiobank')
    cd '/scratch/snkmus003/ukbiobank/heartPaper'
   
    % specificy that we are working on the cluster
    onHPC_Cluster = true;
    numOfCores = 4;
end

% the colors to use the for two groups
groupColors = [ 0.47,0.67,0.19; 0.85,0.33,0.10  ] ;

%% Load the UK Biobank data

fprintf('\n Now loading the data \n')
data = readtable("hpc_ukbiobank_heart.csv");
data.EthnicGroupBW = categorical(data.EthnicGroupBW);

% % get a smaller datasets for the cbio cluster
% myVars = {'PulseRate_AutomatedReading', ...
%     'SystolicBloodPressure_AutomatedReading',...
%     'DiastolicBloodPressure_AutomatedReading',...
%     'MaximumHeartRateDuringFitnessTest'};
% getVars = [{'EthnicGroupBW','BodyMassIndex_BMI_','AgeAtRecruitment',...
%     'Weight','WaistCircumference','StandingHeight'},myVars] ;
% data = data(:,getVars);

% change the variable names if they are not AFR and EUR
if ~any( ismember( data.EthnicGroupBW, {'EUR'}) )
    data.EthnicGroupBW = renamecats(data.EthnicGroupBW, ...
        categories(data.EthnicGroupBW) ,{'AFR','EUR'}) ;
end

%% Run the T-tests: comparing between AFR and White

% % the colors to use the for two groups
% groupColors = [ 0.47,0.67,0.19; 0.85,0.33,0.10  ] ;

% get the variables
myVars = {'PulseRate_AutomatedReading', ...
    'SystolicBloodPressure_AutomatedReading',...
    'DiastolicBloodPressure_AutomatedReading',...
    'MaximumHeartRateDuringFitnessTest'};

figTitles = {'Pulse Rate','Systolic Blood Pressure',...
    'Diastolic Blood Pressure','Max Heart Rate During Fitness Test'};

% set up the unit and the short titles of the variable names
varUnits = {'bpm','mmHg','mmHg','bpm'} ;

% shortTitles = {'Pulse Rate','Max HR During Fitness Test','Diastolic BP'} ;
shortTitles = {'Pulse Rate','Systolic BP','Diastolic BP',...
    'Max Heart Rate'} ;

% preallocated the t-test table
ttestTable = table('Size',[1,7],'VariableTypes',...
    {'string','double','double','double','double','double','double'}, ...
    'VariableNames',{'Variable','meanAFR','meanWhite','tValue', ...
    'lowerBound','upperBound','pValue'});

% run the test in a loop
for ii = 1:length( myVars)
    
    % get the current var name
    curVarName = myVars{ii} ;
    
    % get the data in the current column and fill the outliers using linear
    % interpolation
    curData =  filloutliers( data.(curVarName),'linear') ;

    % perform a ttest with unequal variable assumed
    [~,p,ci,stats] = ttest2(curData(data.EthnicGroupBW =='AFR'),...
        curData(data.EthnicGroupBW == 'EUR'),'Vartype','unequal') ;

    % get the mean values for the two groups
    meanAFR = mean(curData(data.EthnicGroupBW == 'AFR'),'omitnan') ;
    meanEUR = mean(curData(data.EthnicGroupBW == 'EUR'),'omitnan') ;

    % add to the table
    ttestTable.Variable(ii) = figTitles(ii) ;
    ttestTable.meanAFR(ii) = meanAFR;
    ttestTable.meanWhite(ii) = meanEUR ;
    ttestTable.lowerBound(ii) = ci(1);
    ttestTable.upperBound(ii) = ci(2);
    ttestTable.pValue(ii) = p ;
    ttestTable.tValue(ii) = stats.tstat;
    
    % plot the figure
    colourBoxPlot(curData,data.EthnicGroupBW,groupColors,true)
    hold on 
    
    set(gca,'FontSize',15)
    
   % put the p value on the plot depending on its value
    if p < 0.0001
        text(0.5,0.9,['p ', convertPValue2SuperScript(p)],...
            'Units','normalized','FontWeight','bold','FontSize',14, ...
            'HorizontalAlignment','center')
    else
        text(0.5,0.9,sprintf('p = %0.4f',p),...
            'Units','normalized','FontWeight','bold','FontSize',14, ...
            'HorizontalAlignment','center')
    end
    
    % add the title to the plot and a y-axis label
    if ii == 3 % for peak expiratory follow
        ylabel([shortTitles{ii},' (L/min)'])
    else
        ylabel([shortTitles{ii},' (L)'])
    end
    title(shortTitles{ii},'FontSize',16,'FontWeight','bold')
    hold off

end 

hold off 

% also get the summary statistics of the variable of interest 
varStats = readtable('stats_ukbiobank_variables.xlsx');
varStats = varStats(ismember(varStats.Variable , myVars), :) ;

writetable(ttestTable,'Supplementary Data 1.xlsx','Sheet','Welch Results')
writetable(varStats,'Supplementary Data 1.xlsx','Sheet','Variable Stats')

clear aa ans ci curData curVarName groups ii meanAFR meanEUR ...
    p plotNames stats xValues plotName

%% Reviewer Comment -- unnecessary 

% Keep only the rows we need
want = {'mean','std','range','nnz'};
vs = varStats(ismember(varStats.Stats,want), :);

% Pivot to wide format
T = unstack(vs,"African","Stats");   % AFR
T.Properties.VariableNames(2:end) = ...
    strcat(T.Properties.VariableNames(2:end),'_AFR');

T2 = unstack(vs,"White","Stats");    % EUR
T2.Properties.VariableNames(2:end) = ...
    strcat(T2.Properties.VariableNames(2:end),'_EUR');

% Join and rename columns
out = join(T,T2,"Keys","Variable");
out = renamevars(out,["mean_AFR","std_AFR","range_AFR","nnz_AFR", ...
                      "mean_EUR","std_EUR","range_EUR","nnz_EUR"], ...
                     ["Mean AFR","SD AFR","Range AFR","N AFR", ...
                      "Mean EUR","SD EUR","Range EUR","N EUR"]);

writetable(out,"Supplementary Data 1.xlsx","Sheet","Variable Statistics");

%% 0.  Select only the rows we need ------------
want     = {'mean','std','range','nnz'};          % stats to keep
vs       = varStats(ismember(varStats.Stats,want), :);   % filter rows

% 1.  Build AFR wide table
% ---------------------------------------------------
vsAFR    = vs(:, {'Variable','Stats','African'}); % drop the White column
T_AFR    = unstack(vsAFR, 'African', 'Stats', ...
                   'GroupingVariables', 'Variable');

% Rename columns to friendly labels
T_AFR.Properties.VariableNames(2:end) = ...
    {'Mean','N','Range','SD'};
T_AFR = addvars( T_AFR , repmat({'AFR'}, height(T_AFR), 1) ,  ...
    'NewVariableNames',{'Ethnicity'}) ;

% 2.  Build EUR wide table
% ---------------------------------------------------
vsEUR    = vs(:, {'Variable','Stats','White'});   % drop the African column
T_EUR    = unstack(vsEUR, 'White', 'Stats', ...
                   'GroupingVariables', 'Variable');

T_EUR.Properties.VariableNames(2:end) = ...
    {'Mean','N','Range','SD'};
T_EUR = addvars( T_EUR , repmat({'EUR'}, height(T_EUR), 1) ,  ...
    'NewVariableNames',{'Ethnicity'}) ;

% 3.  Join the two tables
% ----------------------------------------------------
out = [T_AFR; T_EUR];

out.Variable = replace( out.Variable,  ...
    {'PulseRate_AutomatedReading', ...
    'DiastolicBloodPressure_AutomatedReading',...
    'MaximumHeartRateDuringFitnessTest', ...
    'SystolicBloodPressure_AutomatedReading'}, ...
    {'Pulse Rate','Diastolic BP','Systolic BP','Max Heart Rate'} ) ;

% Optional: put columns in logical order
out = movevars(out, ...
      {'Ethnicity','Mean','SD','Range','N'},'After', 'Variable');
 

% 4.  Write to Excel
% ---------------------------------------------------------
writetable(out, 'Supplementary Data 1.xlsx', ...
           'Sheet', 'Variable Statistics');

%% Correlation between the cardiac measure 

% MATLAB Code for Group-wise Correlation Analysis and Export
% This script assumes your data is already loaded into a table variable named 'data'.
% It calculates correlations separately for specified groups and saves the
% results in a long-format table to an Excel file.

% --- Step 1: Configuration ---
% Define the measures to analyze (long names from the 'data' table).
% You can change or add to this list.
selectedMeasures = {'PulseRate_AutomatedReading', ...
                    'SystolicBloodPressure_AutomatedReading', ...
                    'DiastolicBloodPressure_AutomatedReading', ...
                    'MaximumHeartRateDuringFitnessTest'};

% Define corresponding short names for the final report.
% Ensure this list has the same number of items as 'selectedMeasures'.
shortNames = {'Pulse Rate', 'Systolic BP', 'Diastolic BP', 'Max Heart Rate'};

% Define the ethnic groups to analyze from the 'EthnicGroupBW' column.
groupsToAnalyze = {'AFR', 'EUR'};

% --- Step 2: Pre-allocate the results table ---
% An empty table is created to store the results from the loops.
finalTable = table('Size', [0, 4], ...
                   'VariableTypes', ...
                   {'string', 'string', 'double', 'double'}, ...
                   'VariableNames', ...
                   {'Measure', 'Ethnicity', 'R_squared', 'P_value'});

% --- Step 3: Loop through groups and pairs of measures ---
% The outer loop iterates through each ethnic group.
for i = 1:numel(groupsToAnalyze)
    currentGroup = groupsToAnalyze{i};
    
    % Filter the main data table for the current group.
    groupData = data(data.EthnicGroupBW == currentGroup, :);
    
    % Check if there is enough data for the analysis.
    if height(groupData) < 2
        warning('Skipping group %s: not enough data rows for analysis.',...
            currentGroup);
        continue;
    end
    
    % The nested loops iterate through each unique pair of measures.
    % The second loop starts from j+1 to avoid duplicate pairs and self-comparisons.
    numMeasures = numel(selectedMeasures);
    for j = 1:numMeasures
        for k = j + 1:numMeasures
            
            % Get the long names of the variables for data extraction.
            measure1_long = selectedMeasures{j};
            measure2_long = selectedMeasures{k};
            
            % Extract the two columns of data for the current pair.
            pairData = groupData{:, {measure1_long, measure2_long}};
            
            % Calculate Pearson correlation coefficients and p-values.
            % 'rows','pairwise' handles any missing NaN values.
            [R, P] = corr(pairData, 'Type','Pearson' ,'Rows','complete');
            
            % Extract the r-value and p-value from the 2x2 matrices.
            % The value at (1,2) is the correlation between the two variables.
            % A check handles cases where correlation could not be computed (e.g., all NaN data).
            if all(isnan(R(:)))
                r_val = NaN;
                p_val = NaN;
            else
                r_val = R(1, 2);
                p_val = P(1, 2);
            end
            
            % Calculate R-squared.
            r_squared = r_val^2;
            
            % Create a descriptive name for the variable pair using the short names.
            variablePairName = sprintf('%s vs. %s', shortNames{j}, shortNames{k});
            
            % Create a new row as a table, which is easy to concatenate.
            newRow = {variablePairName, currentGroup, r_squared, p_val};
            
            % Append the new row to the final results table.
            finalTable = [finalTable; newRow];
        end
    end
end

% --- Step 4: Display the final table in the MATLAB Command Window ---
if isempty(finalTable)
    error('No results were generated. Please check group names and data content.');
end
disp('Combined Phenotypic Correlation Results:');
disp(finalTable);

% --- Step 5: Save the final table to an Excel file ---
fileName = 'SupplementaryData1.xlsx';
sheetName = 'Phenotypic Correlations';
try
    writetable(finalTable, fileName, 'Sheet', sheetName, 'WriteRowNames', false);
    fprintf('\nPhenotypic correlations successfully written to "%s" (sheet: %s).\n', fileName, sheetName);
catch ME
    fprintf('\nError: Could not write to file "%s".\n', fileName);
    fprintf('Please ensure you have write permissions and that the file is not open elsewhere.\n');
    fprintf('Original error message: %s\n', ME.message);
end

%%

% --- Step 6: Generate and Save Heatmap Visualization ---
fprintf('\nGenerating heatmap of overall correlations...\n');

% Define the very short names for the heatmap labels.
heatmapNames = {'Pulse', 'SBP', 'DBP', 'MHR'};

% Calculate the overall correlation matrix from the entire dataset (all groups).
% Round the results to 2 decimal places for cleaner display on the heatmap.
overallCor = round(corrcoef(data{:, selectedMeasures}, 'rows', 'pairwise'), 2);

% Create the heatmap.
figure; % Create a new figure window
h = heatmap(heatmapNames, heatmapNames, overallCor);

% Customize the heatmap appearance.
h.Title = 'Overall Correlation Matrix of Phenotypic Measures';
h.XLabel = 'Measures';
h.YLabel = 'Measures';
% h.Colormap = parula; % e.g., parula, jet


%% Do the Pulse Rate plot against BMI percetile

% A combined one and then one separate for AFR and white
bmiData = data(:,[{'EthnicGroupBW','BodyMassIndex_BMI_', ...
    'AgeAtRecruitment','StandingHeight'},myVars]) ;

% make the variable names easier to work with 
bmiData.Properties.VariableNames = {'Race','BMI','Age','Height',...
    'Pulse Rate','Systolic BP','Diastolic BP','Max Heart Rate'};

% remove the Nan values from teh data 
bmiData(any(isnan(bmiData{:,2:end}),2),:) = [] ;

% set up the counfounder names and make a new short title and weather to
% save the figure
myConfounders = {'BMI','Age','Height'} ;
toSave = false ;

for kk = 1:length(myConfounders)
    
    % get the current confounder and the confounder variable in the table 
    curConf = myConfounders{kk} ;
    confVar = [lower(curConf),'Perc'] ;
    
    % Find the 40th and 60th percentiles of the elements of X.
    confounderPerc = prctile(bmiData.(curConf),[5:5:95]) ;
    
    % add the tenth percentile
    locPercentile = bmiData.(curConf) <= confounderPerc(1) ;
    bmiData.(confVar)(locPercentile) = 5;
    
    % preallocate the percentile to start with
    percVar = 10 ;
    
    % now add a new variable to the table called BMI percentile
    for jj = 1:length(confounderPerc)
        
        % get the location of the percentiles
        locPercentile = bmiData.(curConf) >= confounderPerc(jj) ;
        
        % add to the table the percentile
        bmiData.(confVar)(locPercentile) = percVar ;
        
        % increase by 10
        percVar = percVar + 5 ;
        
    end
    
    % convert to a categorical array
    bmiData.(confVar) = categorical(bmiData.(confVar)) ;
    
    % make two separate data set
    bmiDataBlack = bmiData(bmiData.Race == 'AFR', :) ;
    bmiDataWhite = bmiData(bmiData.Race == 'EUR', :) ;
    
    % get a random samples of the EUR equal to the black
    s = RandStream('dsfmt19937') ;
    smallSample = randsample(s,height(bmiDataWhite),height(bmiDataBlack));
    bmiDataWhite = bmiDataWhite(smallSample,:);

    % get the group stats
    plotDataBlack = grpstats(bmiDataBlack,confVar,{'mean','sem'}, ...
        'DataVars',shortTitles);
    plotDataWhite = grpstats(bmiDataWhite,confVar,{'mean','sem'}, ...
        'DataVars',shortTitles);
    
    % plot the variables in a loop
    for ii = 1:length(shortTitles)
        
        % plot the figure for both AFR and white
        figure()
        errorbar(plotDataBlack.(confVar), ...
            plotDataBlack.(['mean_',shortTitles{ii}]), ...
            plotDataBlack.(['sem_', shortTitles{ii}]),'-o',...
            'Color',groupColors(2,:),...
            'MarkerSize',10,...
            'MarkerEdgeColor',groupColors(2,:),...
            'MarkerFaceColor',groupColors(2,:), ...
            'LineWidth',2)
        
        hold on
        
        errorbar(plotDataWhite.(confVar), ...
            plotDataWhite.(['mean_',shortTitles{ii}]), ...
            plotDataWhite.(['sem_', shortTitles{ii}]),'-o',...
             'Color',groupColors(1,:),...
            'MarkerSize',10,...
            'MarkerEdgeColor',groupColors(1,:),...
            'MarkerFaceColor',groupColors(1,:),...
            'LineWidth',2)
        
        % annotate the plot
        set(gca,'FontSize',16,'LineWidth',1.5,'Box','off')
        xlabel([curConf,' Percentile'])
        ylabel([shortTitles{ii}, ' (', varUnits{ii}, ')'])
        legend({'African','EUR'},'Location','best')
        
        % change the postion of the legend from height 
        if strcmp(curConf,'Height')
           legend({'AFR','EUR'},'Location','SouthEast')
        end
        
        % add the title
        title(figTitles{ii},'FontSize',16,'FontWeight','bold')
        
        hold off
        
        % save the figure
        if toSave == true
            % create Figure Name to filenames
            name = [ replace(figTitles{ii},'/','-'), ' ', curConf] ;
            
            % set the properties of the figure and save it
            set(gcf,'paperunits','centimeters','paperposition',[0 0 30 10])
            print(name,'-dpng','-r300','-vector') % -r600
        end
        
        % ==============================================================
        % Ensure trendStats exists only once
        % ==============================================================
        if ~exist('trendStats', 'var')
            trendStats = table(); %create an empty table on first encounter
        end
        
        % here are the statistical test 
        trait  = shortTitles{ii};   % e.g., 'Pulse Rate'
        pctVar = confVar;           % percentile variable (e.g., 'bmiperc')
        
        % run the trend tests for AFR and EUR (function shown earlier) 
        afrRow = runTrendTests(bmiDataBlack, pctVar, trait, 'AFR', curConf);
        eurRow = runTrendTests(bmiDataWhite, pctVar, trait, 'EUR', curConf);
        
        % Append to the growing results table
        trendStats = [trendStats; afrRow; eurRow];
        
    end
end

% Save or inspect at the end
disp(trendStats)
writetable(trendStats, 'TrendStatistics.xlsx')

% ***************** also add some scatter plots *******************
figure()
gscatter(bmiData.('Pulse Rate'),bmiData.BMI, bmiData.Race,...
    groupColors,'..',10)
set(gca,'LineWidth',1,'Box','off','FontWeight','bold')
xlabel('Pulse Rate') ;
ylabel('BMI') ;
title('Correlation Between Pulse Rate and BMI','FontSize',14',...
    'FontWeight','bold')
hold on

% add a reference line and the correlation coefficent
addReferenceLineToPlot(bmiData.('Pulse Rate'),bmiData.BMI)
legend({'AFR','EUR'} ,'Location','Best')

% remove some of the variables 
clear ii locPercentile bmiPerc percVar plotDataBlack plotDataWhite ...
    smallSample ageDataWhite ageDataBlack curConf kk jj bmiShortTitles ...
    bmiVarUnits bmiFigTitles s

%% Correlation Between the variable of interest

% a figure for SBP and DBP
figure()
gscatter(data.DiastolicBloodPressure_AutomatedReading,...
    data.SystolicBloodPressure_AutomatedReading,...
    data.EthnicGroupBW , groupColors, '..',10)
set(gca,'LineWidth',1.5,'Box','off','FontSize',15)
xlabel('Diastolic BP') ;
ylabel('Systolic BP') ;
title('Correlation Between Pulse Rate and Max Heart Rate',...
    'FontSize',14','FontWeight','bold')
hold on

% add a reference line and the correlation coefficent
addReferenceLineToPlot(data.DiastolicBloodPressure_AutomatedReading,...
    data.SystolicBloodPressure_AutomatedReading)
legend({'Blakcs','EUR'} ,'Location','Best')

% ===================== Bin Scatter Plot =====================
% get the variable to plot here 
plotVars = {'BodyMassIndex_BMI_','Weight','WaistCircumference',...
    'SystolicBloodPressure_AutomatedReading'} ;
yTitles = {'BMI','Weight','Waist Circumference','Systolic BP'} ;
    
% Create a binned scatter plot of some random data points.
for ii = 1:length(plotVars)
    
    figure()
    dbp = data.DiastolicBloodPressure_AutomatedReading ;
    curVar = data.(plotVars{ii}) ;
    binscatter(dbp,curVar,'NumBins',[80 150],'ShowEmptyBins','off')
    colormap(gca,'turbo')
    ax = gca;
    ax.ColorScale = 'log';
    hold on
    
    % edit some plot properties.
    xlabel("Diastolic BP (mmHg)") ;
    ylabel(yTitles{ii}) ;
    set(gca,'FontSize',16,'LineWidth',1.5 ,'Box','off')
    
    % calculate the pearson's linear correation for complete rows
    [r2 , pR ] = corr(dbp,curVar, 'Type','Pearson' ,'Rows','complete');
    if pR == 0
        text(0.5, 0.95, strcat( sprintf("R = %0.2f,  P ", r2), ...
            convertPValue2SuperScript(pR) ) ,'FontWeight','bold', ...
            'FontSize',14 ,'Units','normalize', ...
            'HorizontalAlignment','Center')
    else
        text(0.5, 0.95, strcat( sprintf("R = %0.2f,  P = ", r2), ...
            convertPValue2SuperScript(pR) ) ,'FontWeight','bold', ...
            'FontSize',14 ,'Units','normalize', ...
            'HorizontalAlignment','Center')
    end
    
    % add a reference line and the correlation coefficent
    % add a reference line to the plot
    curVar = fillmissing(curVar , 'nearest') ;
    dbp = fillmissing(dbp, 'nearest') ;
    X = [ones(length(dbp),1) dbp] ;
    b1 = X\curVar ;     yCalc1 = X*b1 ;
    plot(dbp,yCalc1,'LineWidth',2,'Color','k')
    
    % save the figures
    saveas(gcf,['bin_scatter_',plotVars{ii},'_.fig'],'fig')
    
end

hold off

clear bmiVar dbp ax r2 pR X yCalc1 yTitles manifest curVar confVar ...
    confounderPerc b1 ans

%% *********************** REVIEVER COMMENT 1 *********************

% the colors to use the for two groups
if ~exist('groupColours','var')
    groupColors = [ 0.47,0.67,0.19; 0.85,0.33,0.10  ] ;
end

% Load candidate SNPs for AFR and EUR datasets
afr_candidate_snps = ...
    readtable(['FUMA_multitrait_cardio_AFR_results/', 'snps.txt']);
eur_candidate_snps = ...
    readtable(['FUMA_multitrait_cardio_EUR_results/', 'snps.txt']);

% Assuming that 'LeadSNP' is a column in your table that identifies the
% lead SNP of each LD block If this column does not exist, you will need to
% determine how to identify lead SNPs based on your dataset

% For AFR dataset
% Group by LeadSNP and find the maximum CADD score within each group
[afr_groups, afr_group_ids] = findgroups(afr_candidate_snps.IndSigSNP);
afr_maxCADD = splitapply(@max, afr_candidate_snps.CADD, afr_groups);

% For EUR dataset
% Group by LeadSNP and find the maximum CADD score within each group
[eur_groups, eur_group_ids] = findgroups(eur_candidate_snps.IndSigSNP);
eur_maxCADD = splitapply(@max, eur_candidate_snps.CADD, eur_groups);

% here is the CADD threshold 
caddThreshold = 10; 

% Find the fraction of groups (LD blocks) with max CADD score above
% threshold
afr_deleterious_fraction = sum(afr_maxCADD > caddThreshold) ...
    / numel(afr_maxCADD);
eur_deleterious_fraction = sum(eur_maxCADD > caddThreshold) ...
    / numel(eur_maxCADD);

% Display the results
fprintf(['Fraction of LD blocks with deleterious lead SNPs',...
    ' in AFR dataset: %.2f%%\n'], afr_deleterious_fraction * 100);
fprintf(['Fraction of LD blocks with deleterious lead SNPs ',...
    'in EUR dataset: %.2f%%\n'], eur_deleterious_fraction * 100);

% Plotting the fractions
% Create an array of the fractions
deleterious_fractions = [afr_deleterious_fraction,...
    eur_deleterious_fraction];

% Plot histograms of CADD scores and a bar plot of deleterious fractions in
% a tiled layout

% Set up the figure and tiled layout
figure();
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'loose');

% here are the letters 
theLetters = 'a':'j';

% Bar plot for deleterious fractions
nexttile;
bar(deleterious_fractions, 'FaceColor', 'flat', 'CData', groupColors);
set(gca, 'xticklabel', {'AFR', 'EUR'}, 'LineWidth', 1, 'Box', 'off', ...
    'FontSize', 14);
title('Fraction of Deleterious SNPs','FontSize',16);
ylabel('Fraction');
% Display the fractions above the bars
for i = 1:length(deleterious_fractions)
    text(i, deleterious_fractions(i), ...
        sprintf('%.0f%%', deleterious_fractions(i) * 100),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% add a letter to the figure and remove the letter from the array
text(-0.15, 1 ,theLetters(1),'Units','normalized',...
    'FontWeight','bold','FontSize',30)
theLetters(1) = [] ;


% here is the second plot 
nexttile;
histogram(afr_candidate_snps.CADD, 'BinWidth', 1, 'FaceColor', ...
    groupColors(1,:), 'EdgeColor', 'none');
title('Distribution in AFR Dataset','FontSize',16);
set(gca, 'LineWidth', 1, 'Box', 'off','FontSize',14);
xlabel('CADD Score');
ylabel('Frequency');
xlim([0, max(afr_candidate_snps.CADD)+1]); % Adjust X-axis to cover the range of CADD scores

% add a letter to the figure and remove the letter from the array
text(-0.15, 1 ,theLetters(1),'Units','normalized',...
    'FontWeight','bold','FontSize',30)
theLetters(1) = [] ;

% Histogram for EUR dataset CADD scores
nexttile;
histogram(eur_candidate_snps.CADD, 'BinWidth', 1, ...
    'FaceColor', groupColors(2,:), 'EdgeColor', 'none');
title('Distribution in EUR Dataset','FontSize',16);
set(gca, 'LineWidth', 1, 'Box', 'off', 'FontSize', 14);
xlabel('CADD Score');
ylabel('Frequency');
xlim([0, max(eur_candidate_snps.CADD)+1]); % Adjust X-axis to cover the range of CADD scores
% ax = gca ; ax.YAxis.Exponent = 0;

% add a letter to the figure and remove the letter from the array
text(-0.15, 1 ,theLetters(1),'Units','normalized',...
    'FontWeight','bold','FontSize',30)

%% Do the enrichment analysis results 

% Define colors for AFR and EUR
groupColors = [ 0.47,0.67,0.19; 0.85,0.33,0.10 ];

% Define the enrichment result files (use consistent naming with
% "_table.txt")
enrichrFiles = { ...
    'GO_Biological_Process_2023', ...
    'GWAS_Catalog_2023', ...
    'MGI_Mammalian_Phenotype_Level_4_2021', ...
    'OMIM_Expanded', ...
    'PhenGenI_Association_2021', ...
    'UK_Biobank_GWAS_v1', ...
};

% Define human-readable plot titles
plotNames = { ...
    'GO - Biological Process', ...
    'GWAS Catalog', ...
    'MGI Mammalian Phenotype', ...
    'PhenGenI Association', ...
    'UK Biobank GWAS', ...
};

% Define groups to compare
myGroups = {'AFR','EUR'};

% Folder where enrichment files are located
myPath = ['/Users/sinkala/Documents/MATLAB/UKBioBank/',...
    'cardiovascularNewAnalysis/enrichR_results'];

% Loop through each enrichment file
for ii = 1:length(enrichrFiles)
    
    % Set up the figure
    figure();
    tiledlayout(1,2,'Padding','compact');
    set(gcf,'Position',[700,500, 800,400]);
    
    % Loop through AFR and EUR results
    for jj = 1:length(myGroups)
        
        fprintf('\nProcessing %s for %s group...\n', ...
            plotNames{ii}, myGroups{jj});
        
        % Construct file path
        fileName = fullfile(myPath, [ lower(myGroups{jj}), '_', ...
            enrichrFiles{ii}, '_enrichment.txt'] ) ;
        
        % Check if file exists
        if ~isfile(fileName)
            warning('File not found: %s. Skipping...', fileName);
            continue;
        end
        
        % Load enrichment results
        curEnrich = readtable(fileName, 'Format', 'auto') ;
        
        % skip when there is no data 
        if isempty( curEnrich)
            continue
        end
        
        % Clean term names
        curEnrich.Term = strtrim(regexprep(curEnrich.Term, ...
            {'\(+\w*','\:+\w*',')','\R-HSA-+\w*', ...
            'Homo sapiens','\w* raw'},''));
        
        % Plot the figure
        nexttile
        Enrichr_One_bar_plot_logP(curEnrich, myGroups{jj}, ...
            groupColors(jj,:));
        title(myGroups{jj}, 'FontWeight','bold', 'FontSize', 12);
        
        % Write to Excel
        writetable(curEnrich, 'Supplementary Data 4.xlsx', ...
            'Sheet', [myGroups{jj},'_',plotNames{ii}]);
    end
    
    % Add a global title to the figure
    sgtitle(plotNames{ii}, 'FontSize', 14, 'FontWeight', 'bold');
    
    % Save the plot
    figOut = fullfile(myPath, [ plotNames{ii} , '.png'] );
    saveas(gcf, figOut, 'png');
    
    close(gcf);
end

clear ii jj curEnrich Terms

%% Process the gwas results from the uk biobank into two set 

% one whould be for africans and the other from europeans
ukGwasFilesNames = {'PulseRate.txt';'SystolicBP.txt';'DiastolicBP.txt';...
    'MaxHeartRate.txt'} ;

% change my vars to the compact form
myVars = extractBefore(ukGwasFilesNames,'.txt') ;

if ~exist('ukb_cardio_manifest.txt','file')
    % get the file names from the manifest files
    manifest = readtable('Pan-UK Biobank phenotype manifest 2.xlsx');
    manifest.trait_type = categorical(manifest.trait_type);
    manifest.phenocode = str2double(manifest.phenocode);
    
    % also get the data for continuous traits
    manifest = manifest(manifest.trait_type == "continuous",:) ;
    
    % get only the files of the variables of interest
    manifest = manifest( ismember(manifest.description, ...
        {'Pulse rate, automated reading',...
        'Diastolic blood pressure, automated reading',...
        'Systolic blood pressure, automated reading',...                               }
        'Maximum heart rate during fitness test'}), :) ;
    
    % save the processed manifest
    writetable(manifest , 'ukb_cardio_manifest.txt')
else
    manifest = readtable('ukb_cardio_manifest.txt') ;
end

% get the file names 
ukGwasFiles = replace(manifest.filename,'tsv.bgz','txt') ;

% check if the variant qc materic file exist 
if ~exist('variant_qc_metrics.txt','file')
    variantqc = readtable('full_variant_qc_metrics.txt') ;
    variantqc = variantqc(:,{'chrom','pos','ref','alt',...
        'rsid','nearest_genes'}) ;
    writetable(variantqc,'variant_qc_metrics.txt') ;
end

% create the files if they do not exits
if ~exist('jass_results/fuma_multitrait_cardio_EUR.txt','file')
    
    % load the qc matrics of the data and match the first variable name to
    % that of the gwas results
    fprintf('\nLoading the variant metrics data \n')
    variantqc = readtable('variant_qc_metrics.txt') ;
    variantqc.Properties.VariableNames([1,5]) = {'chr','ID'};
    variantqc.alt = [];
    
    % parfor loop over the data
    for ii = 1:length(ukGwasFilesNames)
        
        % see if the processed file exists
        if exist(ukGwasFiles{ii},'file')
            % load the clean data from the .txt file
            fprintf('\nLoading the processed GWAS data\n')
            curResults = readtable(ukGwasFiles{ii}) ;
        else
            % load the gwas results from the tsv file
            curResults = readtable( ...
                replace( ukGwasFiles{ii},'txt','tsv') , ...
                'FileType','Text') ;
            
            % process the file it it contains more than the required
            % columns
            if width(curResults) > 8
                fprintf('\nRemoving the unrequired variables from GWAS\n')
                
                % turn only the information of EUR and AFR groups
                curResults = curResults(:,...
                    {'chr','pos','ref','alt','pval_AFR','pval_EUR',...
                    'beta_AFR','beta_EUR','se_AFR','se_EUR'});
                
                % save to a text file
                fprintf('\n Saving the gwas file with less columns\n')
                writetable(curResults,ukGwasFiles{ii}) ;
            end
        end
        
        % merge the current results with the gwas data 
        fprintf('\nMerging with the variant metric data\n')
        curResults = innerjoin(curResults,variantqc, ...
            'Keys',{'chr','pos','ref'});
        
        % *******************************************************
        % convert the p-values to the actual p-values and not the
        % natural log values and remove the p-values that are NaN
        curResults.pval_EUR = exp(curResults.pval_EUR);
        curResults.pval_AFR = exp(curResults.pval_AFR);
        
        % *******************************************************
        
        % ********************* FUMA ******************************
        % save the gwas files for the FUMA analysis
        fuma_AFR = removevars(curResults, ...
            {'pval_EUR','beta_EUR','se_EUR','nearest_genes'} );
        fuma_EUR = removevars(curResults, ...
            {'pval_AFR','beta_AFR','se_AFR','nearest_genes'} );
        % ***********************************************************
        
        % get only the required coloumns for the analysis
        curResults = removevars( curResults ,  ...
            {'beta_AFR','beta_EUR','se_AFR','se_EUR'} ) ;
        
        % print something to the screen
        fprintf('\nCreating new GWAS results for %s: #%d\n',...
            extractBefore(ukGwasFilesNames{ii},'.txt'), ii)
        
        % break it into two files
        gwasAFR = removevars(curResults,'pval_EUR');
        gwasEUR = removevars(curResults,'pval_AFR');
        
        % move the p-value to the end of the file 
        gwasEUR = movevars(gwasEUR ,'pval_EUR','After','nearest_genes');
        gwasAFR = movevars(gwasAFR ,'pval_AFR','After','nearest_genes');
        
        % change the variable names
        gwasAFR.Properties.VariableNames(end) = "P" ;
        gwasEUR.Properties.VariableNames(end) = "P" ;
        
        % change the variable names to upper case to match those required
        % for plotting
        gwasAFR.Properties.VariableNames = upper( ...
            gwasAFR.Properties.VariableNames) ;
        gwasEUR.Properties.VariableNames = upper( ...
            gwasEUR.Properties.VariableNames) ;
        
        % save to two text files
        fprintf('\nSaving the gwas data for %s\n', shortTitles{ii})
        writetable(gwasAFR,['gwas_', ...
            replace(ukGwasFilesNames{ii},'.txt','_AFR.txt')] );
        writetable(gwasEUR,['gwas_', ...
            replace(ukGwasFilesNames{ii},'.txt','_EUR.txt')] );
        
        % ************** also do the same for the fuma data ********
        % let see the head of the file
        fprintf('\n The size of AFR FUMA %s data is %d %d\n', ...
            shortTitles{ii}, height(fuma_AFR) , width(fuma_AFR) )
        fprintf('\n The size of EUR FUMA %s data is %d %d\n', ...
            shortTitles{ii}, height(fuma_EUR) , width(fuma_EUR) )
        
        % change the names to upper case
        fuma_AFR.Properties.VariableNames = upper( ...
            fuma_AFR.Properties.VariableNames ) ;
        fuma_EUR.Properties.VariableNames = upper( ...
            fuma_EUR.Properties.VariableNames ) ;
        
        % replace the variable names to match those of FUMA 
        fuma_EUR.Properties.VariableNames = ...
            {'CHR','BP','A2','A1','P','Beta','SE','SNP'} ;
        fuma_AFR.Properties.VariableNames = ...
            {'CHR','BP','A2','A1','P','Beta','SE','SNP'} ;
        
        % rearrange the variable in the table
        fuma_EUR = fuma_EUR(:,{'CHR','BP','SNP','P','A1','A2','Beta','SE'});
        fuma_AFR = fuma_AFR(:,{'CHR','BP','SNP','P','A1','A2','Beta','SE'});
        
        % **************** save the metal input data *****************
        fprintf('\nSaving the METAL data data %s %d\n', ...
            extractBefore(ukGwasFilesNames{ii},'.txt'), ii)
        writetable( fuma_AFR,['jass_results/', ...
            replace(ukGwasFilesNames{ii},'.txt','_AFR.txt')]);
        writetable( fuma_EUR,['jass_results/', ...
            replace(ukGwasFilesNames{ii},'.txt','_EUR.txt')] );
        % ************************************************************
        
        % remove the p-values that NaN in the data so that the data is
        % small enough for FUMA analysis and those with p-values greater
        % than 0.1
        fuma_EUR(isnan(fuma_EUR.P),:) = [] ;
        fuma_AFR(isnan(fuma_AFR.P),:) = [] ;
        fuma_EUR(fuma_EUR.P > 0.1,:) = [] ;
        fuma_AFR(fuma_AFR.P > 0.1,:) = [] ;
        
        % save the files
        fprintf('\nSaving the FUMA data %s %d\n', ...
            extractBefore(ukGwasFilesNames{ii},'.txt'), ii)
        writetable( fuma_AFR,['fuma_input/fuma_', ...
            replace(ukGwasFilesNames{ii},'.txt','_AFR.txt')]);
        writetable( fuma_EUR,['fuma_input/fuma_', ...
            replace(ukGwasFilesNames{ii},'.txt','_EUR.txt')] );
        
    end
end

% % ************** stop the code from excuting after processing the files 
% fprintf('\nI have stopped the code on line 525\n')
% return

% % save the GWAS resutls to suppementary File 1 
% sig_pulseRate_AFR = readtable('jass_results/all_sig_pulseRate_AFR.txt') ;
% sig_pulseRate_AFR.Properties.VariableNames = ...
%     {'CHR','BP','SNP','P','A1','A2','Beta','SE'} ;
% 
% % save the file to a supplementary fiel 
% writetable( sig_pulseRate_AFR, 'Supplementary File 1.xlsx','Sheet', ...
%     'AFR Pulse Rate SNPs') ; 

%% load the UK bio bank that has SNP frequencies

% create the data for the analysis 
if ~exist('coloc_eur.txt','file')
    
    % process the SNP freq data
    fprintf('\nProcessing SNP colocalisation data \n')
    
    % change this the snps that are present in the UK Biobank from the data
    % which also included the imputed snps check if the variant qc materic
    % file exist
    if ~exist('variantqc','var')
        % load the data
        variantqc = readtable('afr_eur_variant_qc_metrics.txt') ;
        
        % return only the required table variables
        variantqc = variantqc(:,{'rsid','af_EUR','af_AFR'});
    end
    
    % read in the data of the current condition in
    curAFR = readtable('JASS/afr_zscore_mtag.txt','FileType','text') ;
    curEUR = readtable('JASS/eur_zscore_mtag.txt','FileType','text') ;
    
    % throw in an assertion 
    assert( all( strcmp(curAFR.SNP_ID, curEUR.SNP_ID ) ))
    
    % merge the current results with the gwas data
    fprintf('\nMerging with the variant metric data\n')
    curAFR = innerjoin(curAFR, variantqc, 'LeftKeys', {'SNP_ID'} ,...
        'RightKeys', {'rsid'});
    curEUR = innerjoin(curEUR, variantqc, 'LeftKeys', {'SNP_ID'} ,...
        'RightKeys', {'rsid'});
    
    % remove the allele frequence from the data and clean up the variable
    % names 
    curAFR = removevars(curAFR,{'af_EUR'});
    curEUR = removevars(curEUR,{'af_AFR'});
    curAFR.Properties.VariableNames(end) = "freq" ;
    curEUR.Properties.VariableNames(end) = "freq" ;
    curAFR.Properties.VariableNames(1) = "SNP";
    curEUR.Properties.VariableNames(1) = "SNP";
    
    % now covert the zscore to beta values and se values 
    % Calculate SE using Z-score and P-value
    curAFR.SE = abs(curAFR.Zscore) ./ norminv( ...
        1 - curAFR.P_value / 2, 0, 1);
    curEUR.SE = abs(curEUR.Zscore) ./ norminv( ...
        1 - curEUR.P_value / 2, 0, 1);
    
    % Estimate beta using Z-score and SE
    curAFR.beta = curAFR.Zscore .* curAFR.SE;
    curEUR.beta = curEUR.Zscore .* curEUR.SE;
    
    % Estimate varbeta using SE
    curAFR.varbeta = curAFR.SE.^2;
    curEUR.varbeta = curEUR.SE.^2;
    
    % throw in an assertion 
    assert( all( strcmp(curAFR.SNP_ID, curEUR.SNP_ID ) ))
    
    % save the data for the GWAS results with SE and Beta estimate
    fprintf('\nSaving Coloc Data \n')
    writetable(curAFR,'JASS/afr_jass_fuma_input.txt')  
    writetable(curEUR,'JASS/eur_jass_fuma_input.txt')  
    
    % here is the header of the file 
    head(curAFR)
    
    % get the snps significant at the suggestive p-value threshold
    toKeep = curAFR.P_value > 0.01 | curEUR.P_value > 0.01 ;
    curAFR(toKeep, :) = [];
    curEUR(toKeep, :) = [];
    
    % Your dataset should now have the columns "SNP", "beta", "varbeta",
    % "Pvalue", and "freq"
    curAFR = curAFR(:, {'SNP','beta','varbeta', 'P_value', 'freq','SE'});
    curEUR = curEUR(:, {'SNP','beta','varbeta', 'P_value', 'freq','SE'});
    
    % replace some variables names
    curAFR.Properties.VariableNames(4) = "Pvalue";
    curEUR.Properties.VariableNames(4) = "Pvalue";
   
    % save the data to a csv file
    fprintf('\nSaving Coloc Data \n')
    writetable(curAFR,'coloc_afr.txt')  
    writetable(curEUR,'coloc_eur.txt')  
    
else
    % read the data for colocalisation 
    coloc_AFR = readtable('coloc_afr.txt' , ...
        'Format', '%s%f%f%f%f%f', 'Delimiter', ',') ;   
    coloc_EUR = readtable('coloc_eur.txt',...
         'Format', '%s%f%f%f%f%f', 'Delimiter', ',') ;
    
     % replace some variables names
     coloc_AFR.Properties.VariableNames(4) = "p";
     coloc_EUR.Properties.VariableNames(4) = "p";
     
%      % add the SE from the jass input comment out this code because it
%      % will cause bugs in future
%      coloc_AFR = innerjoin( coloc_AFR, jassrep_AFR(:,{'SNP','SE'}) );
%      coloc_EUR = innerjoin( coloc_EUR, jassrep_EUR(:,{'SNP','SE'}) );
      
     % save the data to a csv file
     fprintf('\nSaving Coloc Data \n')
     writetable(coloc_AFR,'coloc_afr.txt')
     writetable(coloc_EUR,'coloc_eur.txt')
end

% % ************** stop the code from excuting after processing the files
% fprintf('\nI have stopped the code on line 744\n')
% return

% %% check for SNPs for are present in both groups 
% 
% % Regarding your second question, an "Inf" Z-score essentially means that
% % the observed effect size is extremely large and suggests that there might
% % be an error in the data. In practice, this might occur if the standard
% % error (SE) is very small, leading to an extremely large Z-score when the
% % beta coefficient is divided by the SE.
% % 
% % To convert a Z-score back to a beta coefficient, you would need to know
% % the standard error (SE) of the beta coefficient. The formula to calculate
% % a Z-score from a beta coefficient and its SE is Z = beta / SE. Therefore,
% % to get the beta coefficient from a Z-score, you can rearrange this
% % formula to beta = Z * SE.
% % 
% % However, in the case of an infinite Z-score, you can't calculate a
% % meaningful beta coefficient or standard error. This is because an
% % infinite Z-score suggests that the observed effect is infinitely many
% % standard errors away from zero, which isn't really possible in real data.
% % In this case, you might want to check for errors in your data or consider
% % removing that SNP from your analysis.
% 
% % read the data for colocalisation
% coloc_AFR = readtable('coloc_afr.txt' , ...
%     'Format', '%s%f%f%f%f%f', 'Delimiter', ',') ;
% coloc_EUR = readtable('coloc_eur.txt',...
%     'Format', '%s%f%f%f%f%f', 'Delimiter', ',') ;
% 
% % get the coloc data and change the variable names so that the data works
% % with colocalisation analysis
% coloc_afr_input = coloc_AFR(:, {'SNP','beta','SE','p'}) ;
% coloc_eur_input = coloc_EUR(:, {'SNP','beta','SE','p'} );
% 
% coloc_afr_input.Properties.VariableNames = {'snp', 'beta', 'se', 'pval'} ;
% coloc_eur_input .Properties.VariableNames = {'snp', 'beta', 'se', 'pval'} ;
% 
% % remove the rows that have all NaN in one datasets 
% toGo =  all( isnan(coloc_eur_input{:,2:end}) , 2) | ...
%     all( isnan(coloc_afr_input{:,2:end}) , 2)  ;
% coloc_eur_input(toGo, :) = [] ;
% coloc_afr_input(toGo, :) = [] ;
% 
% % also remove the rows with beta values = 0 , these are coming from the
% % infinite Z-scores data 
% toGo = coloc_eur_input.beta == 0 | coloc_afr_input.beta == 0  ;
% coloc_eur_input(toGo, :) = [] ;
% coloc_afr_input(toGo, :) = [] ;
% 
% % Prior probabilities for H0, H1, H2, H3, H4
% prior = [0.01, 0.84, 0.1, 0.025, 0.025];
% 
% coloc_results = ...
%     colocalization_analysis(coloc_afr_input, coloc_eur_input, prior) ;
% % 
% % % perform a more complex analysis
% % N1 = 2885;
% % N2 = 419353 ;
% % coloc_advanced = colocalization_analysis_advanced(...
% %     coloc_afr_input, coloc_eur_input, prior, N1, N2) ;

%% Plot the Beta Estimate between AFR and EUR 

% read the data for colocalisation
coloc_AFR = readtable('coloc_afr.txt' ,'Format', '%s%f%f%f%f%f', ...
    'Delimiter', ',') ;
coloc_EUR = readtable('coloc_eur.txt','Format', '%s%f%f%f%f%f',...
    'Delimiter', ',') ;
     
% get the data plotting and clean the variable names
curData = innerjoin( coloc_AFR, coloc_EUR ,'Key',{'SNP'}) ;
curData.Properties.VariableNames = replace( ...
    curData.Properties.VariableNames, 'coloc_', '') ;

% get only the candidate SNPs 
afr_candidate_snps = readtable(['FUMA_multitrait_cardio_AFR_results/',...
    'snps.txt']);
eur_candidate_snps = readtable(['FUMA_multitrait_cardio_EUR_results/',...
    'snps.txt']);

% now here are all the candidate SNPs 
all_candidate_snps = [afr_candidate_snps.rsID ; ...
    afr_candidate_snps.IndSigSNP ; eur_candidate_snps.rsID ; ...
    eur_candidate_snps.IndSigSNP ] ;

% % return only the candidate snps 
% curData = curData( ismember( curData.SNP, all_candidate_snps) , :) ;

% get the beta values for africans and 0 z-score results these are coming
% from the infinite Z-scores data
curData(isnan(curData.beta_EUR) | isnan(curData.beta_AFR), :) = [] ;
curData(curData.beta_EUR == 0 | curData.beta_AFR == 0, :) = [] ;

% remove the snps that are same from the data 
[~ , theUnique ] = unique(curData.SNP,'first') ;
curData = curData(theUnique, :) ;

% curData(curData.beta_AFR > 0 & curData.beta_EUR < 0, :) = [] ;
% curData(curData.beta_AFR < 0 & curData.beta_EUR > 0, :) = [] ;

%%

curData = sortrows(curData,'p_AFR', 'ascend') ;

% the colors to use the for two groups
groupColors = [ 0.47,0.67,0.19; 0.85,0.33,0.10  ] ;

figure()
% plot the bar graphs using a tiled layout
tiledlayout(2,1,'padding','compact');

% the letters 
theLetters = ['a','b','c'];

% here are the groups 
theGroups = {'AFR','EUR'} ;

% plot the figure of the beta estimate in a loop 
for ii = 1:2
   
    % here are the tiles 
    nexttile
    
    % sort the data according the to the pvalues in EUR
    curData = sortrows(curData, ['p_',theGroups{ii}], 'ascend') ;

    % get a smaller datasets
    numOfVariants = 200 ;
    plotData = curData(1:numOfVariants,:) ;
    
    % plot the figure for both AFR and EUR
    errorbar(1:height(plotData), plotData.beta_EUR,  ...
        plotData.SE_EUR ,'o','vertical',...
        'MarkerSize',10,...
        'MarkerEdgeColor',groupColors(1,:),...
        'MarkerFaceColor',groupColors(1,:), ...
        'LineWidth',2)
    
    hold on
    
    errorbar(1:height(plotData), plotData.beta_AFR,  ...
        plotData.SE_AFR ,'o','vertical',...
        'MarkerSize',10,...
        'MarkerEdgeColor',groupColors(2,:),...
        'MarkerFaceColor',groupColors(2,:),...
        'LineWidth',2)
    
%     % ******************* Here are variants of positive effects ********
%     % get the variants which are have oppositve effect in the two
%     % populaiton 
%     loc_opposite = (plotData.beta_AFR > 0 & plotData.beta_EUR < 0) ...
%         | (plotData.beta_AFR < 0 & plotData.beta_EUR > 0) ;
%     oppositeEffects = plotData( loc_opposite , :) ;
%     
%     % plot the variants of opppostive effects
%     errorbar(1:height(oppositeEffects), oppositeEffects.beta_AFR,  ...
%         oppositeEffects.SE_AFR ,'o','vertical',...
%         'MarkerSize',10,...
%         'MarkerEdgeColor',groupColors(2,:),...
%         'MarkerFaceColor',[1, 0.1 0.1 ],...
%         'LineWidth',2)
%     
%     errorbar(1:height(oppositeEffects), oppositeEffects.beta_EUR,  ...
%         oppositeEffects.SE_EUR ,'o','vertical',...
%         'MarkerSize',10,...
%         'MarkerEdgeColor',groupColors(1,:),...
%         'MarkerFaceColor',[ 1, 0.1 0.1 ], ...
%         'LineWidth',2)
%     % *****************************************************************
    
    % annotate the plot
    set(gca,'FontSize',14,'LineWidth',1.5,'Box','off')
    xlabel('Genetic Variants')
    ylabel('Beta Values')
    legend({'EUR','AFR'},'Location','NorthEast')
    % add the title
    title( sprintf('Top %d most significant SNPs in %s', ...
        numOfVariants, theGroups{ii}), 'FontSize',15)
    
    % add a letter to the figure and remove the letter from the array
    text(-0.05, 1 ,theLetters(1),'Units','normalized',...
        'FontWeight','bold','FontSize',24)
    theLetters(1) = [] ;
end

writetable(curData,'Supplementary Data 2.xlsx','Sheet','Beta Comparison')

clear loc_opposite ii oppositeEffects

%% Get some Very Important SNPs from the Data 

% Identify SNPs with consistent effect direction
sameDir = (curData.beta_AFR > 0 & curData.beta_EUR > 0) | ...
          (curData.beta_AFR < 0 & curData.beta_EUR < 0);
consistentSnps = curData(sameDir, :);

% Display summary
fprintf('\nTotal SNPs with consistent effect direction: %d\n', ...
    height(consistentSnps));

% Find SNPs with strongest significance (smallest p-values in both groups)
sortedConsistent = sortrows(consistentSnps, {'p_AFR','p_EUR'}, ...
    {'ascend','ascend'});

% Print top SNPs
topConsistent = sortedConsistent(1:5,:);
disp('Top 5 consistent-effect SNPs (by p-value):');
disp(topConsistent(:, {'SNP','beta_AFR','p_AFR','beta_EUR', ...
    'p_EUR','freq_AFR','freq_EUR'}));

% Identify the SNP with largest beta_AFR
[~, idxMaxAFR] = max(abs(consistentSnps.beta_AFR));
snpMaxAFR = consistentSnps(idxMaxAFR,:);

fprintf('\nSNP with largest beta in AFR:\n');
disp(snpMaxAFR(:, {'SNP','beta_AFR','p_AFR','beta_EUR','p_EUR', ...
    'freq_AFR','freq_EUR'}));

% Identify the SNP with largest beta_EUR
[~, idxMaxEUR] = max(abs(consistentSnps.beta_EUR));
snpMaxEUR = consistentSnps(idxMaxEUR,:);

fprintf('\nSNP with largest beta in EUR:\n');
disp(snpMaxEUR(:, {'SNP','beta_EUR','p_EUR','beta_AFR','p_AFR', ...
    'freq_AFR','freq_EUR'}));

% Save to a new worksheet
writetable(sortedConsistent, 'Supplementary Data 2.xlsx', ...
    'Sheet','Consistent SNP Effects');

%% printf something to screen 

allPos = sum( curData.beta_AFR > 0 & curData.beta_EUR > 0) ;
fprintf('\nThe number of with positive Beta in EUR and AFR is %d\n',allPos)

allNeg = sum( curData.beta_AFR < 0 & curData.beta_EUR < 0) ;
fprintf('\nThe number of with negative Beta in EUR and AFR is %d\n',allNeg)
    
oppositeBeta = height(curData) - (allNeg+allPos) ;
fprintf('\nThe number of with opposite Beta in EUR and AFR is %d\n',...
    oppositeBeta)

% get the proportion of SNPs with a common effect size 
perc_same_effect = (allNeg+allPos)/height(curData) ;
fprintf(['\nThe percentage of SNPs with same direction beta ', ...
    'in EUR and AFR is %d\n'],perc_same_effect)

clear causalEUR causalAFR ii curData ans theLetters 

%% Try another replication strategy 

% % read the data for colocalisation
% coloc_AFR = readtable('coloc_afr.txt' ,'Format', '%s%f%f%f%f', ...
%     'Delimiter', ',') ;
% coloc_EUR = readtable('coloc_eur.txt','Format', '%s%f%f%f%f',...
%     'Delimiter', ',') ;
%      
% % get the data plotting and clean the variable names
% replicationData = innerjoin( coloc_AFR, coloc_EUR ,'Key',{'SNP'}) ;
% replicationData.Properties.VariableNames = replace( ...
%     replicationData.Properties.VariableNames, 'coloc_', '') ;
% 
% % get only the candidate SNPs 
% afr_lead_snps = readtable(['FUMA_multitrait_cardio_AFR_results/',...
%     'leadSNPs.txt']);
% eur_lead_snps = readtable(['FUMA_multitrait_cardio_EUR_results/',...
%     'leadSNPs.txt']);
% 
% % now here are all the candidate SNPs 
% all_candidate_snps = [afr_lead_snps.rsID ] % ; eur_lead_snps.rsID ] ;
% 
% % return only the candidate snps 
% replicationData = replicationData( ismember( replicationData.SNP, ...
%     all_candidate_snps) , :) ;

%% Get the Top SNPs in EUR and AFR 

% read the data for colocalisation
coloc_AFR = readtable('coloc_afr.txt' ,'Format', '%s%f%f%f%f%f', ...
    'Delimiter', ',') ;
coloc_EUR = readtable('coloc_eur.txt','Format', '%s%f%f%f%f%f',...
    'Delimiter', ',') ;
     
% get the data plotting and clean the variable names
curData = innerjoin( coloc_AFR, coloc_EUR ,'Key',{'SNP'}) ;
curData.Properties.VariableNames = replace( ...
    curData.Properties.VariableNames, 'coloc_', '') ;

% get the beta values for africans and 0 z-score results these are coming
% from the infinite Z-scores data
% curData(isnan(curData.beta_EUR) | isnan(curData.beta_AFR), :) = [] ;

curData = sortrows(curData,'p_AFR', 'ascend') ;
head(curData)

curData = sortrows(curData,'p_EUR', 'ascend') ;
head(curData)

clear coloc_AFR coloc_EUR curData

%% DownLoad the Results of Know Variants from LDLINK

warning('off')
if ~exist('LDLinkDataAFR','dir')
    mkdir('LDLinkDataAFR')
end
addpath('/Manuscript Heart/LDLinkDataAFR')
addpath('/Manuscript Heart/LDLinkDataEUR')
addpath('cardiovascularNewAnalysis/LDLinkDataAFR')
addpath('cardiovascularNewAnalysis/LDLinkDataEUR')

% DO ANOTHER ANALYSIS USING ALL THE POPULATIONS 

% process the data 
if ~exist('gwasCatalog_heartAFR.txt','file')

    % laod the AFR data if it exists
    curAFRdata = readtable(['FUMA_multitrait_cardio_AFR_results/', ...
        'leadSNPs.txt'] ) ;

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
                % get the data 
                system(curl_URL)
                
                % remove the file after reading i
                system(['mv ',filename,' LDLinkDataAFR'])
            catch
                fprintf('\nFailed to DOWNLOAD file for %s\n', ...
                    allVariants{kk} )
                continue ;
            end
        end

        % Load the data into a MATLAB table
        try
            curLinkData = readtable( ['LDLinkDataAFR/', filename],...
                'ReadVariableNames',true) ;
            
            % check if the file has been in correctly
            if ~any(contains(curLinkData.Properties.VariableNames,  ...
                    'NoEntriesInTheGWASCatalogAreIdentified' ))
                if ~any(contains(curLinkData.Properties.VariableNames, ...
                        'GWASTrait'))
                    % read the file in a different way
                    curLinkData = readtable( ...
                        ['LDLinkDataAFR/', filename],...
                        'Delimiter', '\t') ;
                end
            end
            
        catch
            fprintf('\nFailed to READ file for %s\n', allVariants{kk})
            continue
        end

        % check that data are okay
        if any(contains(curLinkData.Properties.VariableNames,  ...
                'NoEntriesInTheGWASCatalogAreIdentified' ))
            curAFRdata.Novelty(kk) = {'Novel'} ;
            curAFRdata.theTraits(kk) = {'None'} ;
        elseif any(contains(curLinkData.Properties.VariableNames,  ...
                'GWASTrait'))
            curAFRdata.Novelty(kk) = {'Has Trait'} ;
            curAFRdata.theTraits(kk) = cellstr( strjoin( ...
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
    writetable(curAFRdata,'gwasCatalog_AFR_Trait_leadSNPs.txt')
else
    gwasCat_AFR = readtable('gwasCatalog_heartAFR.txt');
end

%% ********************************************************************
% ******************* Do the same for the EUR groups  ****************
if ~exist('gwasCatalog_heartEUR.txt','file')

    % laod the AFR data if it exists
    curEURdata = readtable(['FUMA_multitrait_cardio_EUR_results/', ...
        'leadSNPs.txt'] ) ;

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
        if ~exist(['LDLinkDataEUR/',filename],'file')
            try
                system(curl_URL)
            catch
                fprintf('\nFailed to DOWNLOAD file for %s\n', ...
                    allVariants{kk} )
                continue ;
            end
        end
        
        % move the file 
        if exist(filename,'file')
            % remove the file after reading i
            system(['mv ',filename,' LDLinkDataEUR'])
        end

        % Load the data into a MATLAB table
        try
            curLinkData = readtable( ['LDLinkDataEUR/',filename],...
                'ReadVariableNames',true) ;
            
            % check if the file has beeen in correctly
            if ~any(contains(curLinkData.Properties.VariableNames,  ...
                    'NoEntriesInTheGWASCatalogAreIdentified' ))
                if ~any(contains(curLinkData.Properties.VariableNames, ...
                        'GWASTrait'))
                    % read the file in a different way
                    curLinkData = readtable(...
                        ['LDLinkDataEUR/',filename],...
                        'Delimiter', '\t') ;
                end
            end
        catch
            fprintf('\nFailed to READ file for %s\n', allVariants{kk})
            continue
        end

        % check that data are okay
        if any(contains(curLinkData.Properties.VariableNames,  ...
                'NoEntriesInTheGWASCatalogAreIdentified' ))
            curEURdata.Novelty(kk) = {'Novel'} ;
            curEURdata.theTraits(kk) = {'None'} ;
        elseif any(contains(curLinkData.Properties.VariableNames,  ...
                'GWASTrait'))
            curEURdata.Novelty(kk) = {'Has Trait'} ;
            curEURdata.theTraits(kk) = cellstr( strjoin( ...
                curLinkData.GWASTrait,';') ) ;
        end

        % check that data are okay
        if any(contains(curLinkData.Properties.VariableNames,'error')) ...
                || any(contains(curLinkData.(1){1},'error')) ...
                || any(contains(curLinkData.(1),'Gateway Timeout') )  
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

        % add the results to the growing table 
        gwasCat_EUR = [gwasCat_EUR; curLinkData] ;

    end

    % save the data to a file
    writetable(gwasCat_EUR,'gwasCatalog_heartEUR.txt');
    writetable(curEURdata,'gwasCatalog_EUR_Trait_leadSNPs.txt')
else 
    gwasCat_EUR = readtable('gwasCatalog_heartEUR.txt') ;
end

% ************ # exit the code here ***************

%% Noventy of Variants 

fprintf('\nHere is the information about the Significant AFR SNPs\n')

% here are the types of traits that I will use to anonate the snps 
% from the reported snps we the those associated with heart disease and
% related phenotypes
HeartAssocSnps = {'BMI','length','Body mass index','Body size',...
    'Hip circumference','weight','Waist','LDL cholesterol','HDL',...
    'kidney disease','Low density lipoprotein cholesterol levels',...
    'Glomerular filtration rate','Chronic kidney disease',...
    'Body fat percentage'} ;

% here are the phenotypes reported in gwas studies 
HeartReportedSnps = {'Systolic blood','Pulse pressure','Heart rate', ...
     'Electrocardio','Diastolic blood','Cardiac structure and function',...
     'heart','Blood pressure','pulse','Cardiac','cardio',...
     'hypertension','Ischemic stroke','Coronary artery','atrial', ...
     'Myocardial infarction','Cardiovascular disease',...
     'Myocardial','Coronary heart disease','aortic','aorta'} ;
 
% process the GWAS Catalog data return only only cardiovascular trait
% studies 
if ~exist('gwasCatalog_database_heart.txt','file')
    
    % load the GWAS catalog data 
    gwasHeart = readtable('gwas_catalog_v1.0-studies_r2025-03-08.tsv', ...
        'FileType','text');
    
    % return only the required studies
    gwasHeart = gwasHeart( contains( gwasHeart.DISEASE_TRAIT,  ...
        [HeartReportedSnps, HeartAssocSnps], 'IgnoreCase',true ), :) ;
    
    % save the files
    writetable( gwasHeart, 'gwasCatalog_database_heart.txt')
else
    gwasHeart = readtable('gwasCatalog_database_heart.txt');
end
 
% Next, we aimed to identify the novel and known variants among the
% predicted causal variants associated with cardiovascular function.

% load the Tag SNPs by Ldlink analysis and the FUMA analysis 
afrCatalog = readtable('gwasCatalog_heartAFR.txt') ;
afr_fuma_Catalog = readtable( ...
    'FUMA_multitrait_cardio_AFR_results/gwascatalog.txt') ; 

% it turns out that the LD LINK analysis show more variants associated with
% cardiovascular traits than that FUMA analysis 
fprintf('\nThe number of SNPs in Catalog using LDLink is %d\n', ...
    length(unique(afrCatalog.Query)) )
fprintf('\nThe number of SNPs in Catalog using FUMA is %d\n', ...
    length(unique( afr_fuma_Catalog.IndSigSNP)) )

% % check that all the snps in FUMA are present in LD Link 
% assert(all( ismember( afr_fuma_Catalog.IndSigSNP, afrCatalog.Query ) ) )

% here is the other table showing novel traits and know traits 
afr_trait = readtable('gwasCatalog_AFR_Trait_leadSNPs.txt') ;

% add the trait type to the table 
afr_trait = addvars( afr_trait, afr_trait.Novelty,  ...
    'NewVariableNames',{'NoveltyClass'}, 'After','Novelty') ;

% reannoate those variant 
afr_trait.NoveltyClass( contains( ...
    afr_trait.theTraits ,HeartAssocSnps,'IgnoreCase',true )) = ...
    {'Cardiovascular Associated Trait'}  ;
afr_trait.NoveltyClass( contains( ...
    afr_trait.theTraits ,HeartReportedSnps,'IgnoreCase',true )) = ...
    {'Cardiovascular Trait'}  ;
afr_trait.NoveltyClass(ismember(afr_trait.NoveltyClass ,'Has Trait')) = ...
    {'Novel'} ;

% find the number of snps in each class 
afr_trait.NoveltyClass = categorical(afr_trait.NoveltyClass) ;
fprintf('\n Here is a summary of hte traits\n')
summary(afr_trait.NoveltyClass)
length(unique(afr_trait.rsID))

head(afr_trait)

% save the results to a supplementary file 
writetable( afrCatalog, 'Supplementary Data 3.xlsx', ...
    'Sheet','AFR Known Catalog Traits') ;
writetable( afr_trait, 'Supplementary Data 3.xlsx', ...
    'Sheet','AFR IndSNP Traits') ;

%% ****************** Find the Novel SNPs in EUR **********************

fprintf('\nHere is the information about the Significant EUR SNPs\n')

% load the Tag SNPs by Ldlink analysis and the FUMA analysis 
eurCatalog = readtable('gwasCatalog_heartEUR.txt') ;
eur_fuma_Catalog = readtable( ...
    'FUMA_multitrait_cardio_EUR_results/gwascatalog.txt') ; 

% it turns out that he LD LINK analysis show more variants associated with
% cardiovascular traits than that FUMA analysis 
fprintf('\nThe number of EUR SNPs in Catalog using LDLink is %d\n', ...
    length(unique(eurCatalog.Query)) )
fprintf('\nThe number of EUR SNPs in Catalog using FUMA is %d\n', ...
    length(unique( eur_fuma_Catalog.IndSigSNP)) )

% check that all the snps in FUMA are present in LD Link 
% assert(all( ismember( eur_fuma_Catalog.IndSigSNP,eurCatalog.Query)))

% here is the other table showing novel traits and know traits 
eur_trait = readtable('gwasCatalog_EUR_Trait_leadSNPs.txt') ;
summary(categorical(eur_trait.Novelty))

% ********************* Here is the difference **********************
% I turn out that some of the SNPs that are in eurCatalog are not present
% in teh eur_fuma_Catalog because the assertion failed. I have to add those
% SNPs to the table

[~,locDiff] = setdiff( eur_fuma_Catalog.IndSigSNP, eurCatalog.Query ) ;
diff_catalog = eur_fuma_Catalog(locDiff, :) ;

% throw in an assertion 
assert( ~all(ismember(diff_catalog.IndSigSNP,eurCatalog.Query)))

% run the analyis in a loop 
for ii = 1:height(diff_catalog)
    
    if rem(ii,20) == 0
        fprintf('\nRunning analysis for SNP #%d of %d\n', ii, ...
            height(diff_catalog) ) ;
    end
    
    % check that the current SNP in present in the data
    if any( ismember(eur_trait.rsID, diff_catalog.IndSigSNP(ii) ) )
        
        %get the current snp from teh data
        curSnp = diff_catalog.IndSigSNP(ii) ;
        
        % get the location of the current variant in the data 
        [~,locSnp] = intersect( eur_trait.rsID, ...
            diff_catalog.IndSigSNP(ii), 'Stable') ;
    
        % add the current trait to the table
        eur_trait.Novelty(locSnp) = {'Has Trait'} ;
        
        % do it different for empty rows and those that are not empty
        if isempty(eur_trait.theTraits{locSnp})
            eur_trait.theTraits(locSnp) = diff_catalog.Trait(ii)  ;
        else
            % for those that dont have a trait
            if strcmp( eur_trait.theTraits(locSnp) ,'None')
                
                % add the trait
                eur_trait.theTraits(locSnp) = diff_catalog.Trait(ii)  ;
            else
                % add the trait to the other traits
                eur_trait.theTraits(locSnp) = strcat( ...
                    eur_trait.theTraits(locSnp),';', ...
                    diff_catalog.Trait(ii) ) ;
            end
        end
    end   
end

% process data for the variants that do not have standard SNP Ids 
for ii = 1:height(eur_trait)
    
    if rem(ii,20) == 0
        fprintf('\nRunning analysis for in Catalog SNP #%d of %d\n',ii, ...
            height(eur_trait) ) ;
    end
    
    % for the snps that dont have rsID
    if ~contains(eur_trait.rsID(ii) ,'rs') || ...
            isempty(eur_trait.Novelty{ii})
        
        % get the IdSigSNPs for the SNP 
        curIndSNP = split( eur_trait.IndSigSNPs(ii) ,';') ;
        
        % find those Snps in the GWAS catalog for heart function and add
        % the phenotype in type in 
        curCatalog = gwasHeart( ismember( gwasHeart.SNPS , curIndSNP), :);
        
        % now let get the traits
        if isempty(curCatalog)
            % we have now confirmed that the snps has no trait in the GWAS
            % catalogy
            eur_trait.Novelty(ii) = {'Novel'} ;
            eur_trait.theTraits(ii) = {'None'} ;
            
            continue
        end
        
        % now let get the traits associated with the SNPs 
        curCatalogTrait = curCatalog.DISEASE_TRAIT  ;  
        eur_trait.theTraits(ii) = cellstr( strjoin(curCatalogTrait,';')) ;
    end 
end

% ********************* Here is the difference *********************

% add the trait type to the table 
eur_trait = addvars( eur_trait, eur_trait.Novelty,  ...
    'NewVariableNames',{'NoveltyClass'}, 'After','Novelty') ;

% reannoate those variant 
eur_trait.NoveltyClass( contains( ...
    eur_trait.theTraits ,HeartAssocSnps,'IgnoreCase',true )) = ...
    {'Cardiovascular Associated Trait'}  ;
eur_trait.NoveltyClass( contains( ...
    eur_trait.theTraits ,HeartReportedSnps,'IgnoreCase',true )) = ...
    {'Cardiovascular Trait'}  ;
eur_trait.NoveltyClass(ismember(eur_trait.NoveltyClass,'Has Trait')) = ...
    {'Novel'} ;

% find the number of snps in each class 
eur_trait.NoveltyClass = categorical(eur_trait.NoveltyClass) ;
fprintf('\nHere is a summary of the traits in EUR\n')
summary(eur_trait.NoveltyClass)
length(unique(eur_trait.rsID))

% get the top Novel variants 
eur_novel = eur_trait( eur_trait.NoveltyClass == "Novel", :) ;
eur_novel = sortrows(eur_novel,'p','ascend') ;
fprintf('\nhere are the top EUR novel variants\') 
eur_novel(1:5,{'rsID','p','NoveltyClass'})

% get the top Known variaints 
eur_known = eur_trait( eur_trait.NoveltyClass == "Cardiovascular Trait",:);
eur_known = sortrows(eur_known,'p','ascend') ;
fprintf('\nhere are the top EUR Cardiovascular Trait variants\') 
eur_known(1:5,{'rsID','p','NoveltyClass'})

% get the top associated variants 
eur_assoc = eur_trait( eur_trait.NoveltyClass ==  ...
    "Cardiovascular Associated Trait",:);
eur_assoc = sortrows(eur_assoc,'p','ascend') ;
fprintf('\nhere are the top EUR Cardiovascular Associated Trait snps\') 
eur_assoc(1:5,{'rsID','p','NoveltyClass'})

% save the results to a supplementary file 
writetable(eurCatalog,'Supplementary Data 3.xlsx', ...
    'Sheet','EUR Known Catalog Traits') ;
writetable(eur_trait,'Supplementary Data 3.xlsx', ...
    'Sheet','EUR IndSNP Traits') ;

clear curCatalogTrait ii curCatalog curIndSNP

%% Create the Multitrait GWAS FUMA input 

% one whould be for africans and the other from europeans
if ~exist('jass_results/fuma_multitrait_cardio_AFR.txt','file')
    
    % print something to the screen
    fprintf('\nCreating new FUMA results for the multitrait analysis\n')
 
    % try reading the files that do not exist
    try
        % these are the old files 
        fuma_AFR = readtable( ...
            'jass_results/afr_multitrait_results_1.tbl', ...
            'FileType','text');
        fuma_EUR = readtable( ...
            'jass_results/eur_multitrait_results_1.tbl', ...
            'FileType','text');
    
    % these are the variable names 
    % SNP_ID Allele1 Allele2 Weight	Zscore	P-value	Direction
    
    % remove the p-values that NaN in the data so that the data is
    % small enough for FUMA analysis and those with p-values greater
    % than 0.1
    fuma_EUR(isnan(fuma_EUR.P_value),:) = [] ;
    fuma_AFR(isnan(fuma_AFR.P_value),:) = [] ;
    fuma_EUR(fuma_EUR.P_value > 0.1,:) = [] ;
    fuma_AFR(fuma_AFR.P_value > 0.1,:) = [] ;
    
    % get the unique rows
    fuma_AFR = unique(fuma_AFR) ;
    fuma_EUR = unique(fuma_EUR) ;
    
    % and change the last value for the p-values to P
    % SNP_ID,Allele1,Allele2,Weight,Zscore,P_value,Direction
    fuma_AFR.Properties.VariableNames(6) = "P" ;
    fuma_EUR.Properties.VariableNames(6) = "P" ;
    fuma_AFR.Properties.VariableNames(1) = "ID" ;
    fuma_EUR.Properties.VariableNames(1) = "ID" ;
    
    % remove the variable that are not required 
    fuma_AFR = removevars(fuma_AFR,{'Allele1','Allele2'}) ;
    fuma_EUR = removevars(fuma_EUR,{'Allele1','Allele2'}) ;
    
    % aso save the manhattan data with the followng column names 
    % [T.CHR,T.BP,T.P]; this will have to be added to the file 
    chromPosAFR = readtable('gwas_DiastolicBP_AFR.txt') ;
    chromPosEUR = readtable('gwas_DiastolicBP_EUR.txt') ;
    
    % return only the required variables 
    chromPosAFR = chromPosAFR(:,{'CHR','POS','REF','ALT','ID',...
        'NEAREST_GENES'}) ;
    chromPosEUR = chromPosEUR(:,{'CHR','POS','REF','ALT','ID',...
        'NEAREST_GENES'}) ;
    
    % join the two table 
    fuma_AFR = innerjoin(chromPosAFR, fuma_AFR, 'Keys', {'ID'}) ;
    fuma_EUR = innerjoin(chromPosEUR, fuma_EUR, 'Keys', {'ID'}) ;
    
    % read the data that is already processed 
    catch
        
        % Load the AFR and EUR FUMA inputs
        fprintf('\nReading the JASS results \n')
        fuma_AFR = readtable('JASS/afr_jass_fuma_input.txt');
        fuma_EUR = readtable('JASS/eur_jass_fuma_input.txt');
        
        % Step 4: Rename variables
        fuma_AFR = renamevars( fuma_AFR,{'P_value','chrom'}, {'P','CHR'}) ;
        fuma_EUR = renamevars( fuma_EUR,{'P_value','chrom'}, {'P','CHR'}) ;
        
        % Step 3: Remove duplicate rows (based on all columns)
        fuma_AFR = unique(fuma_AFR);
        fuma_EUR = unique(fuma_EUR);
        
        % Step 4: Rename variables
        fuma_AFR.Properties.VariableNames{'SNP'} = 'ID';
        fuma_EUR.Properties.VariableNames{'SNP'} = 'ID';
        
        % Step 5: Remove columns not needed (REF in this case) If you had
        % Allele1/Allele2 before, you'd drop those. Now we remove 'REF'.
        fuma_AFR = removevars(fuma_AFR, {'ref'});
        fuma_EUR = removevars(fuma_EUR, {'ref'});
        
        % Step 2: Read only the required columns: rsid (SNP ID) and
        % nearest_genes
        fprintf('\nNow reading the variant matrix file\n')
        geneInfo = readtable('JASS/variant_qc_sig_snps_in_jass.txt');
        
        % Step 3: Rename the rsid column to match 'ID' in your fuma tables
        geneInfo.Properties.VariableNames{'rsid'} = 'ID';
        
        % return only the unique genes 
        [~ , theUnique ] = unique(geneInfo.ID ) ;
        geneInfo = geneInfo(theUnique, :) ;
        
        % Step 4: Join to add gene names
        fuma_AFR_joined = innerjoin( ...
            geneInfo(:, {'ID','ref','alt','nearest_genes'}), fuma_AFR,  ...
             'Keys', 'ID');
         
        fuma_EUR_joined = innerjoin( ...
               geneInfo(:, {'ID','ref','alt','nearest_genes'}), ...
               fuma_EUR,  'Keys', 'ID');
        
%         % Step 5: Assertion — make sure no SNPs were lost
%         assert(height(fuma_AFR_joined) == height(fuma_AFR), ...
%             'Some SNPs were lost during AFR gene name join');
%         
%         assert(height(fuma_EUR_joined) == height(fuma_EUR), ...
%             'Some SNPs were lost during EUR gene name join');
        
        % Step 6 (Optional): Replace the original tables
        fuma_AFR = fuma_AFR_joined;
        fuma_EUR = fuma_EUR_joined;
        
        % convert the variable names to up case and then arrange them in a
        % particular order 
        fuma_AFR.Properties.VariableNames = upper( ...
            fuma_AFR.Properties.VariableNames) ;
        fuma_EUR.Properties.VariableNames = upper( ...
            fuma_EUR.Properties.VariableNames) ;
        
        % arrange them in a particular order
        fuma_AFR = fuma_AFR(:, {'CHR','POS','REF','ALT','ID',...
            'NEAREST_GENES','P','SE','BETA','ZSCORE'}) ;
        fuma_EUR = fuma_EUR(:, {'CHR','POS','REF','ALT','ID',...
            'NEAREST_GENES','P','SE','BETA' ,'ZSCORE'}) ;
        
        % change the last variable names 
        fuma_AFR.Properties.VariableNames(end) = "Zscore" ;
        fuma_EUR.Properties.VariableNames(end) = "Zscore" ;
       
        clear fuma_AFR_joined fuma_EUR_joined qcFile geneInfo
        
    end
    
%     % I need all the data for replacation creation - save the files
%     fprintf('\nSaving the FUMA data for multitrait analysis\n')
%     writetable(fuma_AFR,'jass_results/all_multitrait_cardio_AFR.txt');
%     writetable(fuma_EUR,'jass_results/all_multitrait_cardio_EUR.txt' );
    
    % Step 1: Remove rows with NaN p-values
    fuma_AFR(isnan(fuma_AFR.P), :) = [];
    fuma_EUR(isnan(fuma_EUR.P), :) = [];
    
    % Step 2: Remove rows with P > 0.1
    fuma_AFR(fuma_AFR.P > 0.1, :) = [];
    fuma_EUR(fuma_EUR.P > 0.1, :) = [];
        
    % save the files
    fprintf('\nSaving the FUMA data for multitrait analysis\n')
    writetable(fuma_AFR,'jass_results/fuma_multitrait_cardio_AFR.txt');
    writetable(fuma_EUR,'jass_results/fuma_multitrait_cardio_EUR.txt' );
    
end

clear chromPosAFR chromPosEUR ii

%% Create the second jass input for replacition 

if ~exist('jass_results/jass_rep_afr_data.txt','file')
    
    % save the files
    fprintf('\nload jass data for jass replication\n')
    jassrep_AFR = readtable( ...
        'jass_results/fuma_multitrait_cardio_AFR.txt');
    jassrep_EUR = readtable(...
        'jass_results/fuma_multitrait_cardio_EUR.txt');
    
    % process the data for another meta analysis with jass
    % Calculate SE using Z-score and P-value
    
    % check first if the data contains the Beta values
    if any(strcmp( jassrep_AFR.Properties.VariableNames, 'BETA'))
        % change the variable name case 
        jassrep_AFR = renamevars(jassrep_AFR, "BETA", "beta") ;
        jassrep_EUR = renamevars(jassrep_EUR, "BETA", "beta") ;
        
    else % computer the SE when they are missing 
        jassrep_AFR.SE = abs(jassrep_AFR.Zscore) ./ norminv( ...
            1 - jassrep_AFR.P / 2, 0, 1);
        jassrep_EUR.SE = abs(jassrep_EUR.Zscore) ./ norminv( ...
            1 - jassrep_EUR.P / 2, 0, 1);
        
        % Estimate beta using Z-score and SE
        jassrep_AFR.beta = jassrep_AFR.Zscore .* jassrep_AFR.SE ;
        jassrep_EUR.beta = jassrep_EUR.Zscore .* jassrep_EUR.SE ;
    end
    
    % jass input requried
    
    % # === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
    % SEPARATOR COMMA
    % DEFAULT 5978
    % MARKER SNP
    % ALLELE A2 A1 % A2 if the reference allele
    % EFFECT Beta
    % PVALUE P
    % STDERR SE
    % PROCESS DiastolicBP_AFR.txt
    %
    % % jass results have this
    % CHR
    % POS
    % REF
    % ALT
    % ID
    % NEAREST_GENES
    % Weight
    % Zscore
    % P
    % Direction
    
    % get only the required variables CHR,BP,SNP,P,A1,A2,Beta,SE
    jassrep_AFR = jassrep_AFR(:, ...
        {'CHR','POS','ID','P','ALT','REF','beta','SE'}) ;
    jassrep_EUR = jassrep_EUR(:, ...
        {'CHR','POS','ID','P','ALT','REF','beta','SE'}) ;
    
    % change the variable names
    jassrep_AFR.Properties.VariableNames = ...
        {'CHR','BP','SNP','P','A1','A2','Beta','SE'} ;
    jassrep_EUR.Properties.VariableNames = ...
        {'CHR','BP','SNP','P','A1','A2','Beta','SE'} ;
    
    % save the file to text file in jass folder
    writetable( jassrep_AFR, 'jass_results/jass_rep_afr_data.txt')
    writetable( jassrep_EUR, 'jass_results/jass_rep_eur_data.txt')
    
end

%% Make the Manhattan Plot for the GWAS data 

% here are the snps to exclude

% check for the file in the directory
if ~exist('jass_results/fuma_multitrait_cardio_AFR_ManhattanPlot.fig',...
        'file')
    
    fprintf('\n Plotting the Manhattan plots\n')
    
    % plot the manhattam plots in a loop - I use the FUMA which has
    % fewer snps so that the plot data point can be fewer
    figure()
    ManhattanPlot('jass_results/fuma_multitrait_cardio_AFR.txt',...
        'title','Cardiovascular traits - AFR')
    
    figure()
    ManhattanPlot('jass_results/fuma_multitrait_cardio_EUR.txt',...
        'title','Cardiovascular traits - EUR')
end

%% Here is teh combined plots -- which is better 

% here are teh groups 
myGroups = {'EUR','AFR'};

% here are the letters
theLetters = 'a':'j';

% Create a new figure and tiled layout
figure('Position', [100, 100, 1400, 800]);
t = tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

for jj = 1:length(myGroups)
    
    % Current figure file
    curFigTitle = ['jass_results/fuma_multitrait_cardio_', ...
        myGroups{jj},'_ManhattanPlot.fig'];
    
    % Open the figure invisibly
    tempFig = openfig(curFigTitle, 'invisible');
    
    % Get the axes handle from the opened figure
    axOld = findall(tempFig, 'type', 'axes');
    
    % Create a subplot for this group
    nexttile(t, jj);
    axNew = copyobj(axOld, gcf);
    
    % Adjust the copied axis to fit in the subplot tile
    set(axNew, 'Position', get(gca, 'Position'));
    delete(gca); % Remove the dummy axis
    
    % Custom formatting
    set(axNew,'FontSize',16,'LineWidth',1.5,'Box','off');
    title(axNew, ['Cardiovascular traits - ',myGroups{jj}],...
        'FontSize',16,'FontWeight','normal', ...
        'Units','normalized','Position',[0.5, 0.92]);
    
    % add a letter to the figure and remove the letter from the array
    text(-0.03, 1 ,theLetters(jj),'Units','normalized',...
        'FontWeight','bold','FontSize',28)
    theLetters(1) = [] ;
    
    xlabel(axNew, '');
    close(tempFig);  % Close the temporary figure
end

% Save combined plot
set(gcf, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 50 25])
print(gcf, 'jass_results/ManhattanPlot_Combined.png', '-dpng', '-r300');

%% Compare SNPs frequency with SNPs Assocations

% Some SNPs are associated with FEV1, FVC and PEF. Are they found at the
% same or different frequencies among EUR and black

% load the UK bio bank that has SNP frequencies
if exist('jass_results/snp_freq_heart.csv','file')
    
    % load the SNP freq data that is processed
    fprintf('\nLoading the processed SNP freq data \n')
    snpFreq = readtable('jass_results/snp_freq_heart.csv') ;
    
else
    % process the SNP freq data
    fprintf('\nProcessing SNP freq data \n')
    
    % change this the snps that are present in the UK Biobank from the data
    % which also included the imputed snps
    snpFreq = readtable('JASS/variant_qc_full_sig_snps_in_jass.txt');
    
    % return only the required table variables
    snpFreq = snpFreq(:,{'chrom','pos','alt','ref','rsid',...
        'af_EUR','af_AFR'});
    
    % change the variables names to those present in the snp comparison
    % table
    snpFreq.Properties.VariableNames = {'Chrom','Position',...
        'Alternative','Reference','Variant','EUR','AFR'};
    
    % read in the data of the current condition in
    curAFR = readtable('jass_results/fuma_multitrait_cardio_AFR.txt',...
        'FileType','text') ;
    curEUR = readtable('jass_results/fuma_multitrait_cardio_EUR.txt',...
        'FileType','text') ;
    
    % get the snps significant at the suggestive p-value threshold
    curAFR(curAFR.P > 5e-8, :) = [];
    curEUR(curEUR.P > 5e-8, :) = [];
    
    % add those condition to the table where the SNPs is true
    locAFR = ismember(snpFreq.Variant, curAFR.ID) ;
    locEUR = ismember(snpFreq.Variant, curEUR.ID) ;
    
    % add the variable to the table
    snpFreq.('gwas_P_AFR') = locAFR ;
    snpFreq.('gwas_P_EUR') = locEUR ;
    
    % return only the SNPs that are significant in the GWAS of the three
    % variables
    % I dont know why indexing by variable
    targetVars = {'gwas_P_EUR','gwas_P_AFR'};
    
    snpFreq = snpFreq( any(snpFreq{:,{'gwas_P_EUR','gwas_P_AFR'}},2), :);
    
    % add a new variable at the end of the column to show which SNPs
    % related to EUR and those related to AFR
    snpFreq.SigEUR = any(snpFreq{:,{'gwas_P_EUR'}},2) ;
    snpFreq.SigAFR = any(snpFreq{:,{'gwas_P_AFR'}},2) ;
    
    % get instance where the gwas is significant in EUR
    snpFreq.SigGroup(snpFreq.SigEUR==true) = {'EUR'} ;
    snpFreq.SigGroup(snpFreq.SigAFR==true) = {'AFR'} ;
    
    % get the instances where both snps are significant
    bothSig = all([snpFreq.SigEUR,snpFreq.SigAFR],2) ;
    snpFreq.SigGroup(bothSig)= {'Both'};
    
    % sort the rows of the data
    snpFreq = sortrows(snpFreq,'EUR','descend');
    
    % save the data to a csv file
    fprintf('\nSaving the SNP freq data to a csv files\n')
    writetable(snpFreq,'jass_results/snp_freq_heart.csv')
    
end

% also add for both groups
snpFreq.SigGroup  = categorical(snpFreq.SigGroup );

% disply the signficant snps for each group
summary(snpFreq.SigGroup)

clear locAFR locEUR ii curAFR curEUR

% % ************** stop the code from excuting after processing the files 
% fprintf('\nI have stopped the code on line 716\n')
% return

%% FUMA analysis 

% I have put the process files for FUMA analysis here 
% /scratch/snkmus003/ukbiobank/heartPaper/fuma_input
warning('off','all');

% This part of script get statistc of the common snps, the lead snps, and
% the independent snps, and also the genes within which of the snps are
% located based on the magma analysis. I will also need to plot scatter
% plot location for which the lead variant in AFR are unique from those in
% EUR if such are found.

if ~exist('jass_results/snpFreq_lead_snps.csv','var')

    fprintf('\nProcessing the FUMA results \n')
    
    % load the snpfreq data
    if ~exist('snpFreq','var')
        % here is the snp freq data
        snpFreq = readtable('jass_results/snp_freq_heart.csv') ;
        % also add for both groups
        snpFreq.SigGroup  = categorical(snpFreq.SigGroup );
    end
    
    % preallocate the table for the significant snps
    numSigSnps = array2table(myVars,'VariableNames',{'Trait'} ) ;
    
    % pre allocate the table for all the lead snp
    all_indSig_snp = [] ;
    all_magma_genes = [] ;
    
    % plot venn diagram for the common snps while considering the LD of the
    % snps from the FUMA results
    
    % get the number of SNP ids found to be significant
    numSigSnps.all_AFR(1) = sum(snpFreq.SigAFR) ;
    numSigSnps.all_EUR(1) = sum(snpFreq.SigEUR) ;
    
    % read the fuma results of the SNPs for the current variable
    curEURdata = readtable(['FUMA_multitrait_cardio_EUR_results/', ...
        'IndSigSNPs.txt'] ) ;
    
    % laod the AFR data if it exists
    try
        curAFRdata = readtable(['FUMA_multitrait_cardio_AFR_results/', ...
            'IndSigSNPs.txt'] ) ;
    catch
        % prepare the variable types for the current data
        varTypes = cell(1,width(curEURdata)) ;
        for vv = 1:width(curEURdata )
            varTypes{vv} = class(curEURdata.(vv)) ;
        end
        curAFRdata = table('Size',[0,width(curEURdata)], ...
            'VariableNames',curEURdata.Properties.VariableNames, ...
            'VariableTypes',varTypes);
    end
    
    % add the number of lead snps
    numSigSnps.IndSigSNPs(1) = length(unique( curAFRdata.rsID)) ;
    numSigSnps.IndSigSNPs(1) = length(unique( curEURdata.rsID)) ;
   
    % ***** I HAVE CONSIDERED THE LEAD SNPS AS DISCOVERIES ****
    % sometimes the is no african data so I have to use a try catch
    % statemeent
    % read the fuma results of the SNPs for the current variable
    curEURdata = readtable(['FUMA_multitrait_cardio_EUR_results/', ...
        'leadSNPs.txt'] ) ;
    % laod the AFR data if it exists
    try
        curAFRdata = readtable(['FUMA_multitrait_cardio_AFR_results/',...
            'leadSNPs.txt'] ) ;
    catch
        % prepare the variable types for the current data
        varTypes = cell(1,width(curEURdata)) ;
        for vv = 1:width(curEURdata )
            varTypes{vv} = class(curEURdata.(vv)) ;
        end
        curAFRdata = table('Size',[0,width(curEURdata)], ...
            'VariableNames',curEURdata.Properties.VariableNames, ...
            'VariableTypes',varTypes);
    end
   
    % add the number of lead snps
    numSigSnps.leadSNPs_AFR(1) = length(unique( curAFRdata.rsID)) ;
    numSigSnps.leadSNPs_EUR(1) = length(unique( curEURdata.rsID)) ;
    
    % load the LD data from for the snps
    eurLD = readtable(['FUMA_multitrait_cardio_EUR_results/', ...
        'ld.txt'] ,'FileType','Text','Delimiter','tab') ;
    
    % prepare the variable types for the current data
    try
        afrLD = readtable(['FUMA_multitrait_cardio_AFR_results/', ...
            'ld.txt'] ,'FileType','Text','Delimiter','tab') ;
    catch
        varTypes = cell(1,width(eurLD)) ;
        for vv = 1:width( eurLD)
            varTypes{vv} = class(eurLD.(vv)) ;
        end
        afrLD = table('Size',[0,width(eurLD)], ...
            'VariableNames',eurLD.Properties.VariableNames, ...
            'VariableTypes',varTypes);
    end
    
    % put the LD data together and remove the rows with the same data
    ldData = unique([afrLD; eurLD]) ;
    ldData( strcmp(ldData.SNP1, ldData.SNP2), :)  = [] ;
    
    % ********************* LD Cross Ref Analysis ********************
    
    % add variable to the table for the top ld variant and R square value
    curAFRdata = addvars( curAFRdata,  ...
        repmat( {'none'}, height(curAFRdata),1), ...
        zeros(height(curAFRdata),1), 'NewVariableNames' , ...
        {'top_ldsnp','ld_r2'}  );
    curEURdata = addvars( curEURdata,  ...
        repmat( {'none'}, height(curEURdata),1), ...
        zeros(height(curEURdata),1), 'NewVariableNames' , ...
        {'top_ldsnp','ld_r2'}  );
    
    % *********************** LD Consideration *********************
    % loop over the variants in AFR and see they occur in EUR groups
    for jj = 1:height(curAFRdata)
        
        % get the current variant
        curSnpId = curAFRdata.uniqID(jj) ;
        
        % do a snps id to snp id comparison for before doing the LD
        % maping
        if any(ismember(curAFRdata.rsID(jj), curEURdata.rsID))
            
            % for the african group
            curAFRdata.top_ldsnp(jj) = curAFRdata.rsID(jj);
            curAFRdata.ld_r2(jj) = 1 ;
            
            % for the european group
            [~,Bi] = intersect(curEURdata.rsID, curAFRdata.rsID);
            curEURdata.top_ldsnp(Bi) = curAFRdata.rsID(jj) ;
            curEURdata.ld_r2(Bi) = 1 ;
            
            % go back to the top of the loop
            continue
        end
        
    end
    
    % ****************************************************************
    try
        % check if the data has already been processed
        ldLinkData = readtable('heart_ldlink_data_AFR.txt') ;
        ldLinkData = movevars(ldLinkData,{'IndexSNP'},'Before',1) ;
    catch
        % here are all the variables and an empty table 
        allVariants = [curEURdata.rsID ;curAFRdata.rsID ] ;
        ldLinkData = [] ;
        
        % beging a parallel works 
        if onHPC_Cluster == true
            % work on the entire datasets
            parpool('local',numOfCores)
        end
        
        % run a parrel for loop here 
        for kk = 1:length(allVariants)
            
            % here is the file name
            filename = ['lead_snp_files_AFR/', allVariants{kk},'.txt'] ;
            
            % file the file does not exist on the cluster
            if ~exist(filename,'file')
                
                fprintf(['\nDownloading LD data from LDlink.com ',  ...
                    'for %s #%d of %d\n'],  allVariants{kk}, kk, ...
                    length(allVariants) )
                
                % here are the curl options
                curl_URL = sprintf("curl -k -X GET " + ...
                    "'https://ldlink.nci.nih.gov/LDlinkRest/" + ...
                    "ldproxy?var=%s&pop=AFR&r2_d=r2&window=500000&" + ...
                    "genome_build=grch37&token=97fb8311b229'" + ...
                    " -o %s", allVariants{kk}, filename ) ;
                
                % get the current variant from only from online curSnp
                url = sprintf(['https://ldlink.nci.nih.gov/LDlinkRest/',...
                    'ldproxy?var=%s&pop=EUR&r2_d=r2&window=500000&',...
                    'genome_build=grch37&token=97fb8311b229'], ...
                    allVariants{kk}) ;
                
                % here are the options
                options = weboptions('Timeout',30) ;
                
                % try to download the data five times
                try
                    % use curl to down load the data
                    system(curl_URL)
                catch
                    try
                        % use curl to down load the data
                        websave( filename, url, options) ;
                    catch
                        fprintf('\nFailed to DOWNLOAD file for %s\n', ...
                            allVariants{kk} )
                        continue ;
                    end
                end 
            end
            
            % Load the data into a MATLAB table
            try
                fprintf(['\nRunning LD data for varaints',  ...
                    ' %s #%d of %d\n'],  allVariants{kk}, kk, ...
                    length(allVariants) )
                % read the data from the save file 
                curLinkData = readtable(filename) ;
            catch
                fprintf('\nFailed to READ file for %s\n', allVariants{kk} )
                continue
            end
           
            
            % check that data are okay
            if any(contains(curLinkData.Properties.VariableNames,'error'))
                continue
            end
            if isempty(curLinkData)
                fprintf('\nFile NOT okay for %s\n', allVariants{kk})
                continue
            end
            
            % get only the variant with R square > 0.1, and get only the
            % SNP Ids and the Rsqurae values, and the coodinates
            try
                curLinkData = curLinkData(curLinkData.R2 > 0.1, ...
                    {'RS_Number','Coord','R2'}) ;
            catch
                fprintf('\nFile NOT okay for %s\n', allVariants{kk})
                continue
            end
            
            % change the vairable names and add the Index SNP to the data
            curLinkData.Properties.VariableNames([1,3]) =  ...
                {'LinkedSNP','Rsquare'} ;
            curLinkData = addvars( curLinkData, ...
                repmat( allVariants(kk), height(curLinkData), 1),  ...
                'NewVariableNames', {'IndexSNP'} , 'Before', 1) ;
            
            % clean up the variable names and the data
            curLinkData.Position = extractAfter(curLinkData.Coord,':');
            curLinkData.Chrom = extractAfter( extractBefore( ...
                curLinkData.Coord,':') , 'chr') ;
            curLinkData.Position = str2double(curLinkData.Position);
            curLinkData.Chrom = str2double( curLinkData.Chrom) ;
            curLinkData.Coord = [] ;
            
            % remove the rowws where the index SNP and the Linked SNP
            curLinkData( strcmp( ....
                curLinkData.IndexSNP, curLinkData.LinkedSNP), :)  = [] ;
            
            % if no data remains
            if isempty(curLinkData)
                fprintf('\nNo variants in LD with %s\n', allVariants{kk})
                continue
            end
            
            % let see the header of the file 
            if rem(kk,100) == 0
                head(curLinkData)
            end
            
            % remove the data with no linked snp
            ldLinkData = [ldLinkData; curLinkData] ;          
        end
        
        % let see how the data looks like
        disp( ldLinkData(1:10,:) )
        
        % save the data to a file 
        writetable(ldLinkData,'heart_ldlink_data_AFR.txt');
    end
    
    % #################### Call the LD function as well ###############
    
    % specify the rCutOff and windowSize
    wndwSize = 500000 ;
    rCutOff = 0.1 ;
    
    % this function does more LD pruning of the variants associated
    % with the traits in AFR and EUR
    [~, curEURdata] = ld_crossref_cardio( ldLinkData, curEURdata, ...
        rCutOff, wndwSize) ;
    [~, curAFRdata] = ld_crossref_cardio( ldLinkData, curAFRdata, ...
        rCutOff, wndwSize) ;
    
    % #################################################################
    
    % ***************** Combined the results into one table ***********
    curEURdata = addvars( curEURdata, ...
        repmat({'EUR'}, height(curEURdata), 1), ...
        'NewVariableNames',{'Ethnicity'} ,'Before',1) ;
    
    curAFRdata = addvars( curAFRdata, ...
        repmat({'AFR'}, height(curAFRdata), 1), ...
        'NewVariableNames',{'Ethnicity'} ,'Before',1) ;
    
    % here are all the snps put into table
    all_indSig_snp = [all_indSig_snp; [curEURdata; curAFRdata] ]  ;
    
    % reduct the snpFreq data to return only the lead snps from the FUMA
    % analysis
    snpFreq = snpFreq(ismember(snpFreq.Variant, all_indSig_snp.rsID), :) ;
    
    % save the new snpFreq data to a csv file
    writetable(snpFreq,'jass_results/snpFreq_lead_snps.csv') ;
    
    % get the current afrSnp with ld snps consdiers
    afrSnps_ld = curAFRdata.rsID;
    locLdsnps = ~ismember( curAFRdata.top_ldsnp, 'none') ;
    afrSnps_ld(locLdsnps) = curAFRdata.top_ldsnp(locLdsnps) ;
    
    % plot the venn diagram of SNPs unique in each population
    commonSnps = intersect(curEURdata.rsID, afrSnps_ld) ;
    curEurSnps = setdiff(curEURdata.rsID, afrSnps_ld);
    curAFRSnps = setdiff(afrSnps_ld,curEURdata.rsID);
    
    % plot a venn diagram of the snps draw the venn diagram for the
    % SNPs that are intersecting
    figure()
    venn([200,200],75,'EdgeColor','black')
    hold on
    
    % add text to the figure: get numbers of intersects to add to the
    % plots
    myNumText(3) = length(curAFRSnps) ;
    myNumText(2) = length(commonSnps);
    myNumText(1) = length(curEurSnps) ;
    
    % here are text positions and the titles of the plots
    textPos = [0.25, 0.47 , 0.70] ;
    for jj = 1:length(textPos)
        text(textPos(jj),0.5,num2str(myNumText(jj)),...
            'Units','normalized',...
            'FontWeight','bold','FontSize',16,...
            'HorizontalAlignment','center')
    end
    
    % add the titles to the circle
    myGroups = {'EUR','AFR'};
    textPos = [0.33 , 0.66] ;
    for jj = 1:length(textPos)
        text(textPos(jj), 0.98 ,[myGroups{jj}],'Units','normalized',...
            'FontWeight','bold','FontSize',14,...
            'HorizontalAlignment','center')
    end
    hold off
    
    % save the figure
    saveas(gcf,'jass_results/Venn_Diagram_of_SNPs_05.fig','fig')
    
    % #################################################################
    % ############## Venn Diagram of Gene Based Test ##################
    
    % read the fuma results of the SNPs for the current variable
    curEURgenes = readtable( ...
        ['FUMA_multitrait_cardio_EUR_results/', ...
        'magma.genes.out'] ,'FileType','Text','Delimiter','tab') ;
    % load the AFR data if it exists
    curAFRgenes = readtable( ...
        ['FUMA_multitrait_cardio_AFR_results/', ...
        'magma.genes.out'],'FileType','Text','Delimiter','tab' );
    
    % adjust the p-values for multiple comparison
   
    % Input SNPs were mapped to 18514 protein coding genes. Genome wide
    % significance (red dashed line in the plot) was defined at P =
    % 0.05/18514 = 2.701e-6
    afrCutoff = 0.05/18514 ;
    eurCutoff = 0.05/18514;
    
    % return only the signficant genes at the p-values of less than
    % 0.05
    curAFRgenes = curAFRgenes( curAFRgenes.P < afrCutoff, ...
        {'CHR','NSNPS','P','SYMBOL'}) ;
    curEURgenes = curEURgenes( curEURgenes.P < eurCutoff, ...
        {'CHR','NSNPS','P','SYMBOL'}) ;
    
    % combine the results in a table
    curEURgenes = addvars( curEURgenes, ...
        repmat({'EUR'}, height(curEURgenes), 1), ...
        'NewVariableNames',{'Ethnicity'} ,'Before',1) ;
    
    curAFRgenes = addvars( curAFRgenes, ...
        repmat({'AFR'}, height(curAFRgenes), 1), ...
        'NewVariableNames',{'Ethnicity'} ,'Before',1) ;
    
    % here are all the snps put into table
    all_magma_genes = [all_magma_genes; [curEURgenes; curAFRgenes] ] ;
    
    % plot the venn diagrams for the genes plot the venn diagram of
    % SNPs unique in each population
    commonGenes = unique( intersect(...
        curEURgenes.SYMBOL,curAFRgenes.SYMBOL));
    onlyEurGenes = unique( setdiff(...
        curEURgenes.SYMBOL,curAFRgenes.SYMBOL));
    onlyAFRgenes = unique( setdiff(...
        curAFRgenes.SYMBOL,curEURgenes.SYMBOL));
    
    % plot a venn diagram of the snps draw the venn diagram for the
    % SNPs that are intersecting
    figure()
    venn([200,200],75,'EdgeColor','black')
    hold on
    
    % add text to the figure: get numbers of intersects to add to the
    % plots
    myNumText(3) = length(onlyAFRgenes) ;
    myNumText(2) = length(commonGenes);
    myNumText(1) = length(onlyEurGenes) ;
    
    % here are text positions and the titles of the plots
    textPos = [0.25, 0.47 , 0.70] ;
    for jj = 1:length(textPos)
        text(textPos(jj),0.5,num2str(myNumText(jj)),'Units', ...
            'normalized','FontWeight','bold','FontSize',16,...
            'HorizontalAlignment','center')
    end
    
    % add the titles to the circle
    myGroups = {'EUR','AFR'};
    textPos = [0.33 , 0.66] ;
    for jj = 1:length(textPos)
        text(textPos(jj), 0.98 ,[myGroups{jj}], ...
            'Units','normalized','FontWeight','bold','FontSize',14,...
            'HorizontalAlignment','center')
    end
    hold off
    
    % save the figure
    saveas(gcf,'jass_results/Venn_Diagram_of_Genes.fig','fig')
    
    writetable(all_magma_genes,'jass_results/all_magma_genes.csv')
    
    % also save the list of genes for the phenotype to an excel fiel 
    writetable( curEURgenes,'jass_results/each_magma_genes.xlsx', ...
        'Sheet','EUR - MAGMA')
    writetable( curAFRgenes,'jass_results/each_magma_genes.xlsx', ...
        'Sheet','AFR - MAGMA')
     
    % ##############################################################
    % plot the final venn diagram for when all all the snps put
    % together
    
    % save the results to an Excel files
    writetable(all_indSig_snp,'jass_results/all_lead_snps.csv')
    
else
    % laod the processed data
    snpFreqLead = readtable('jass_results/snpFreq_lead_snps.csv');
end

%% ***************** PROCESS THE LD DATA FOR THE EUR *******************
% *********************************************************************

try
    % check if the data has already been processed
    ldLinkData = readtable('heart_ldlink_data_EUR.txt') ;
    ldLinkData = movevars(ldLinkData,{'IndexSNP'},'Before',1) ;
catch
    
    % ***** I HAVE CONSIDERED THE LEAD SNPS AS DISCOVERIES ****
    % sometimes the is no african data so I have to use a try catch
    % statemeent
    % read the fuma results of the SNPs for the current variable
    curEURdata = readtable(['FUMA_multitrait_cardio_EUR_results/', ...
        'leadSNPs.txt'] ) ;
    curAFRdata = readtable(['FUMA_multitrait_cardio_AFR_results/', ...
        'leadSNPs.txt'] ) ;
    
    % here are all the variables and an empty table
    allVariants = unique([curEURdata.rsID ;curAFRdata.rsID ]);
    ldLinkData = [] ;
    
    % start the parallel pool
    try
        % beging a parallel works
        if onHPC_Cluster == true
            % work on the entire datasets
            parpool('local',numOfCores)
        else
            error('\nThis section must be on a high memory computer\n')
        end
    catch
        % catch the error for the if statement
        fprintf('\nParallel pool already running\n')
    end
    
    % run a parrel for loop here
    for kk = 1:length(allVariants)
        
        % here is the file name
        filename = ['lead_snp_files_EUR/',allVariants{kk},'.txt'] ;
        
        % file the file does not exist on the cluster
        if ~exist(filename,'file')
            
            fprintf(['\nDownloading LD data from LDlink.com ',  ...
                'for %s #%d of %d\n'],  allVariants{kk}, kk, ...
                length(allVariants) )
            
            % here are the curl options
            curl_URL = sprintf("curl -k -X GET " + ...
                "'https://ldlink.nci.nih.gov/LDlinkRest/" + ...
                "ldproxy?var=%s&pop=EUR&r2_d=r2&window=500000&" + ...
                "genome_build=grch37&token=97fb8311b229'" + ...
                " -o %s.txt", allVariants{kk},  filename ) ;
            
            % get the current variant from only from online curSnp
            url = sprintf(['https://ldlink.nci.nih.gov/LDlinkRest/',...
                'ldproxy?var=%s&pop=EUR&r2_d=r2&window=500000&',...
                'genome_build=grch37&token=97fb8311b229'], allVariants{kk}) ;
            
            % here are the options
            options = weboptions('Timeout',30) ;
            
            % try to download the data five times
            try
                % use curl to down load the data
                system(curl_URL)
            catch
                try
                    % use curl to down load the data
                    websave(filename,url,options);
                catch
                    fprintf('\nFailed to DOWNLOAD file for %s\n', ...
                        allVariants{kk} )
                    continue ;
                end
            end
        else
            fprintf('\nThe File Exist: for %s #%d of %d\n', ...
                allVariants{kk}, kk, length(allVariants) )
        end
    end
    
    % run a parrel for loop here
    for kk = 1:length(allVariants)
        
        % here is the file name
        filename = ['lead_snp_files_EUR/',allVariants{kk},'.txt']  ;
        
        % Load the data into a MATLAB table
        try
            fprintf(['\nRunning LD data for varaints',  ...
                ' %s #%d of %d\n'],  allVariants{kk}, kk, ...
                length(allVariants) )
            % read the data from the save file
            curLinkData = readtable(filename) ;
        catch
            fprintf('\nFailed to READ file for %s\n', allVariants{kk})
            continue
        end
          
        % check that data are okay
        if any(contains(curLinkData.Properties.VariableNames,'error'))
            continue
        end
        if isempty(curLinkData)
            fprintf('\nFile NOT okay for %s\n', allVariants{kk})
            continue
        end
        
        % get only the variant with R square > 0.1, and get only the
        % SNP Ids and the Rsqurae values, and the coodinates
        try
            curLinkData = curLinkData(curLinkData.R2 > 0.1, ...
                {'RS_Number','Coord','R2'}) ;
        catch
            fprintf('\nFile NOT okay for %s\n', allVariants{kk})
            continue
        end
        
        % change the vairable names and add the Index SNP to the data
        curLinkData.Properties.VariableNames([1,3]) =  ...
            {'LinkedSNP','Rsquare'} ;
        curLinkData = addvars( curLinkData, ...
            repmat( allVariants(kk), height(curLinkData), 1),  ...
            'NewVariableNames', {'IndexSNP'} , 'Before', 1) ;
        
        % clean up the variable names and the data
        curLinkData.Position = extractAfter(curLinkData.Coord,':');
        curLinkData.Chrom = extractAfter( extractBefore( ...
            curLinkData.Coord,':') , 'chr') ;
        curLinkData.Position = str2double(curLinkData.Position);
        curLinkData.Chrom = str2double( curLinkData.Chrom) ;
        curLinkData.Coord = [] ;
        
        % remove the rowws where the index SNP and the Linked SNP
        curLinkData( strcmp( ....
            curLinkData.IndexSNP, curLinkData.LinkedSNP), :)  = [] ;
        
        % if no data remains
        if isempty(curLinkData)
            fprintf('\nNo variants in LD with %s\n', allVariants{kk})
            continue
        end
        
        % let see the header of the file
        if rem(kk,100) == 0
            head(curLinkData)
        end
        
        % remove the data with no linked snp
        ldLinkData = [ldLinkData; curLinkData] ;
    end
    
    % let see how the data looks like
    disp( ldLinkData(1:10,:) )
    
    % save the data to a file
    writetable(ldLinkData,'heart_ldlink_data_EUR.txt');
    
    clear afr_snp_in_ld eur_snp_in_ld jj ii myGroups textPos jj ...
        curTitle myNumText pics2_ld_data cur_pics2 curAFR_ld ...
        cur_AFR_data curAFRsnps curEURdata curEurSnps curLd curSnpld ...
        curTitle curVariants
end

% % ************** stop the code from excuting after processing the files
% fprintf('\nI have stopped the code on line 1214\n')
% return

clear ldData eurLD afrLD

%% Create the Coloc Multitrait Cardio Input Files for Replications

if ~exist('jass_results/replic_multitrait_cardio_AFR.txt','file')
    
    % create the tmp file
    if ~exist('jass_results/tmp_multitrait_cardio_AFR.txt','file')
        
        % here is the filter variant matrxi
        fprintf('\nLoading variantQc data for sig snps \n')
        sigVariantqc = readtable("variant_qc_sig_snps_in_jass.txt") ;
        
        % here are the zscore score data 
        fprintf('\nLoading JASS zcore data for sig snps \n')
        zscoreAFR = readtable("JASS/afr_zscore_mtag.txt") ;
        zscoreEUR = readtable("JASS/eur_zscore_mtag.txt") ;
        
        % get an intersection of the snps 
        commonSnps = [ intersect(zscoreAFR.SNP_ID, sigVariantqc.rsid) ; ...
            intersect(zscoreEUR.SNP_ID, sigVariantqc.rsid) ];
        
        % return those common snps in variantqc 
        sigVariantqc = sigVariantqc( ismember(sigVariantqc.rsid , ...
            commonSnps), :) ;
        
        % now join them 
        fprintf('\nProcessing JASS zcore data for sig snps \n')
        zscoreAFR = outerjoin(sigVariantqc , zscoreAFR,  ...
            'LeftKey','rsid','RightKey','SNP_ID','MergeKey',true, ...
            'Type','left' ) ;
        zscoreEUR = outerjoin(sigVariantqc, zscoreEUR, ...
            'LeftKey','rsid','RightKey','SNP_ID','MergeKey',true, ...
            'Type','left') ;
        
        % create a new variable name 
        zscoreAFR.Properties.VariableNames{'rsid_SNP_ID'} = 'SNP';
        zscoreEUR.Properties.VariableNames{'rsid_SNP_ID'} = 'SNP';
        
        % Throw in an assertion to ensure both files contain the same SNPs
        assert( all( strcmp( zscoreAFR.SNP, zscoreEUR.SNP) ), ...
            'SNP IDs do not match between AFR and EUR files');
        
        % remove the common NaN rows and common rows with p-values greater
        % than 0.01 Remove rows with NaN in the p-value column
        rowsToRemove = ...
            isnan(zscoreAFR.P_value) & isnan(zscoreEUR.P_value) | ...
            zscoreAFR.P_value > 0.01 & zscoreEUR.P_value > 0.01 ;
        
        % Remove rows with p-value greater than 0.1 (suggestive threshold)
        zscoreAFR(rowsToRemove,:) = [] ;
        zscoreAFR(rowsToRemove,:) = [] ;
        
        % Throw in an assertion to ensure both files contain the same SNPs
        assert(all(strcmp(zscoreAFR.SNP, zscoreEUR.SNP)), ...
            'SNP IDs do not match between AFR and EUR files');
        
        % Read the temporary files into MATLAB tables
        writetable(zscoreAFR,'jass_results/tmp_multitrait_cardio_AFR.txt');
        writetable(zscoreEUR,'jass_results/tmp_multitrait_cardio_EUR.txt');
        
    end
    
    fprintf('\nCreating Coloc Multitrait Cardio Input Files\n');
    
    % Read the temporary files into MATLAB tables
    replicAFR = readtable('jass_results/tmp_multitrait_cardio_AFR.txt');
    replicEUR = readtable('jass_results/tmp_multitrait_cardio_EUR.txt');
    
    % Throw in an assertion to ensure both files contain the same SNPs
    assert(all(strcmp(replicAFR.SNP, replicEUR.SNP)), ...
        'SNP IDs do not match between AFR and EUR files');
    
    % Now convert the Z-score to SE, beta, and varbeta for each table.
    replicAFR.SE = abs(replicAFR.Zscore) ./ norminv(1 - replicAFR.P_value / 2, 0, 1);
    replicEUR.SE = abs(replicEUR.Zscore) ./ norminv(1 - replicEUR.P_value / 2, 0, 1);
    
    replicAFR.beta = replicAFR.Zscore .* replicAFR.SE;
    replicEUR.beta = replicEUR.Zscore .* replicEUR.SE;
    
    replicAFR.varbeta = replicAFR.SE.^2;
    replicEUR.varbeta = replicEUR.SE.^2;
    
    % covert all variable names to upper case 
    replicAFR.Properties.VariableNames = upper( ...
          replicAFR.Properties.VariableNames) ;
    replicEUR.Properties.VariableNames = upper( ...
          replicEUR.Properties.VariableNames) ;
      
    % rename some variables 
    replicAFR = renamevars( replicAFR, {'CHROM','P_VALUE'}, {'CHR','P'}) ;
    replicEUR = renamevars( replicEUR, {'CHROM','P_VALUE'}, {'CHR','P'}) ;
    
    % rename some variable
    replicAFR = replicAFR(:, {'CHR','POS','REF','ALT','SNP', ...
          'NEAREST_GENES','P','SE','BETA','VARBETA','ZSCORE'}) ;
    replicEUR = replicEUR(:, {'CHR','POS','REF','ALT','SNP', ...
          'NEAREST_GENES','P','SE','BETA','VARBETA','ZSCORE'}) ;

    % (Optional) Reorder or select columns if needed; here we keep:
    % SNP, beta, varbeta, p, freq, SE
    replicAFR.Properties.VariableNames(end) = "Zscore" ;
    replicEUR.Properties.VariableNames(end) = "Zscore" ;
    
    % Throw in an assertion to ensure both files contain the same SNPs
    assert(all(strcmp(replicAFR.SNP, replicEUR.SNP)), ...
        'SNP IDs do not match between AFR and EUR files');
    
    % Save the processed files
    writetable(replicAFR,'jass_results/replic_multitrait_cardio_AFR.txt');
    writetable(replicEUR,'jass_results/replic_multitrait_cardio_EUR.txt');
    
else
    fprintf('\nColoc Multitrait Cardio Input Files already exist.\n');
end

%% Prepare the file for variant replication

% process the file
if ~exist('jass_results/gwas_all_EUR_4_variants_replication.txt','file')
    
    % Print status
    fprintf('\nCreating the file for replication analysis\n')
    
    % =============== AFR DATA ===============
    
    % 1) Read the new AFR association file
    %    (e.g., fuma_multitrait_cardio_AFR.txt or afr_jass_fuma_input.txt)
    gwas_all_AFR = readtable( ...
        'jass_results/replic_multitrait_cardio_AFR.txt','FileType', 'text');
    
    % 2) Remove duplicate rows
    gwas_all_AFR = unique(gwas_all_AFR);
    
    % 3) Rename columns using renamevars for a clean implementation.
    %    For example, if the file contains 'SNP_ID' or 'P-value', rename
    %    them to 'ID' and 'P', respectively.
    if any(strcmp(gwas_all_AFR.Properties.VariableNames, 'SNP'))
        gwas_all_AFR = renamevars(gwas_all_AFR, 'SNP', 'ID');
    end
    if any(strcmp(gwas_all_AFR.Properties.VariableNames, 'P-value'))
        gwas_all_AFR = renamevars(gwas_all_AFR, 'P-value', 'P');
    end
    
    % 4) Remove variables not required (if they exist)
    varsToRemove = {'Allele1','Allele2'};
    varsToRemove = varsToRemove(ismember(varsToRemove, ...
        gwas_all_AFR.Properties.VariableNames));
    if ~isempty(varsToRemove)
        gwas_all_AFR = removevars(gwas_all_AFR, varsToRemove);
    end
    
    % 5) Read the variant QC file (which replaces the old
    % gwas_DiastolicBP_AFR.txt)
    chromPos = readtable('JASS/variant_qc_sig_snps_in_jass.txt');
    [~, theUnique ] = unique(chromPos.rsid) ;
    chromPos = chromPos(theUnique,:) ;
    
    % 6) Keep only the relevant columns and rename them for joining
    %    The QC file originally has: chrom, pos, ref, alt, rsid,
    %    nearest_genes
    chromPos = chromPos(:, ...
        {'chrom','pos','ref','alt','rsid','nearest_genes'});
    chromPos = renamevars(chromPos, ...
        {'chrom','pos','ref','alt','rsid','nearest_genes'}, ...
        {'CHR','POS','REF','ALT','ID','NEAREST_GENES'});
    
    % 7) Join the two tables on 'ID'
    gwas_all_AFR = innerjoin( chromPos, ... 
        gwas_all_AFR(:, {'ID','P','SE','BETA','Zscore'}), ...
        'Keys', {'ID'} ,'MergeKey',true, 'Type','lef');
    
    % 8) Save the final AFR replication file
    writetable(gwas_all_AFR, ...
        'jass_results/gwas_all_AFR_4_variants_replication.txt');
    
    % =============== EUR DATA ===============
    
    % 1) Read the new EUR association file
    %    (e.g., fuma_multitrait_cardio_EUR.txt or eur_jass_fuma_input.txt)
    gwas_all_EUR = readtable( ...
        'jass_results/all_multitrait_cardio_EUR.txt', 'FileType', 'text');
    
    % 2) Remove duplicate rows
    gwas_all_EUR = unique(gwas_all_EUR);
    
    % 3) Rename columns using renamevars
    if any(strcmp(gwas_all_EUR.Properties.VariableNames, 'SNP_ID'))
        gwas_all_EUR = renamevars(gwas_all_EUR, 'SNP_ID', 'ID');
    end
    if any(strcmp(gwas_all_EUR.Properties.VariableNames, 'P-value'))
        gwas_all_EUR = renamevars(gwas_all_EUR, 'P-value', 'P');
    end
    
    % 4) Remove variables not required (if present)
    varsToRemove = {'Allele1','Allele2'};
    varsToRemove = varsToRemove(ismember(varsToRemove, ...
        gwas_all_EUR.Properties.VariableNames));
    if ~isempty(varsToRemove)
        gwas_all_EUR = removevars(gwas_all_EUR, varsToRemove);
    end
    
    % 7) Join the two tables on 'ID'
    gwas_all_EUR = outerjoin( chromPos, ...
        gwas_all_EUR(:, {'ID','P','SE','BETA','Zscore'}), ...
        'Keys', {'ID'} ,'MergeKey',true ,'Type','left' );
    
    % 8) Save the final EUR replication file
    writetable(gwas_all_EUR, ...
        'jass_results/gwas_all_EUR_4_variants_replication.txt');
    
     % Throw in an assertion to ensure both files contain the same SNPs
    assert(all(strcmp(gwas_all_EUR.ID, gwas_all_AFR.ID)), ...
        'SNP IDs do not match between AFR and EUR files');
    
    % Clean up workspace variables
    clear chromPosAFR chromPosEUR gwas_all_EUR gwas_all_AFR
    
end

%% *********** Replicate the EUR Sig variants in AFR group ************* 

% check if the variants have been replicated in AFR from EUR
if ~exist('replicated_snps_in_afr.txt','file')
    
    % check that the ld data exist
    if ~exist('gwasStats_afr','var')
        
        % check if the data has already been processed
        ldLinkData = readtable('heart_ldlink_data_AFR.txt') ;
        ldLinkData = movevars(ldLinkData,{'IndexSNP'},'Before',1) ;
        
        % get the afr only snps and change the variable names for all the
        % snps and see which among those can be replicated
        gwasStats_afr = readtable(...
            'jass_results/gwas_all_AFR_4_variants_replication.txt', ...
            'FileType','text') ;
        loc_rsid = find( ismember( ...
            gwasStats_afr.Properties.VariableNames,'ID')) ;
        gwasStats_afr.Properties.VariableNames(loc_rsid) = "rsid" ;
        
    end
    
    % here are the EUR sig snps that need to be replace
    eur_id_snps = readtable(['FUMA_multitrait_cardio_EUR_results/', ...
        'leadSNPs.txt'] ) ;
    
    % multily the independent sig snps by the LD snps 
    
    % ******** try to use bash for this -- it will be faster **********
    % THE APPROACH IS SIMPLE, I GET THE ldLinkData, JOIN IT TO THE 
    % gwasStats_afr AND THEN JOIN THE eur_id_snps TO THAT. THIS SHOULD GIVE
    % ME THE LIST OF SNPS IN EUR, AND THE SNPS IN LD WITH THOSE IN THE 
    % ldLinkData, IF THOSE SNPS ARE SIGNFICANT IN TEH gwasStats_afr data
    
    % here are the replicated snps
    replicated_snps_afr  = [] ;
    
    % start the parallel pool 
    try
        % beging a parallel works
        if onHPC_Cluster == true
            % work on the entire datasets
            parpool('local',numOfCores)
        else
            error('\nThis section must be on a high memory computer\n')
        end
    catch
        % catch the error for the if statement 
        fprintf('\nParallel pool already running\n')
    end
    
    % let check how many variables we have in the workspace 
    whos
    
    % find the variants in the intergrated data are associated with
    % the traits
    parfor jj = 1:height(eur_id_snps)
        
        % print something to the screen
        if rem(jj,20) == 0
            fprintf('\nRunning analysis for %s #%d of %d\n', ...
                eur_id_snps.rsID{jj}, jj, height(eur_id_snps) )
        end
        
        % here is the snps to be used for replication
        theSnp = eur_id_snps.rsID(jj) ;
        
        % replicate the results in the EUR group
        cur_rep = replicateVariant(theSnp, ldLinkData, gwasStats_afr);
        
        % put in the same table
        replicated_snps_afr = [replicated_snps_afr ;cur_rep] ;
    end
    
    % some linked snps are '.' - remove those
    replicated_snps_afr(ismember(replicated_snps_afr.IndexSNP,'.'),:) = [];
    
    % save the file
    writetable(replicated_snps_afr,'replicated_snps_in_afr.txt');
    
    clear gwasStats_afr
    
    % check that the snp data dose not exist for only the african variants
end

if ~exist('replicated_snps_in_AFR_only.txt','file')
    
    % load the file
    replicated_snps_afr = readtable('replicated_snps_in_afr.txt') ;
    
    % ******************* replicate clean up ************************
    % change some variables names in the replicated snps data for the
    % AFR groups
    
    % the "IndexSNP" is the the column for snp Ids in AFR grop that
    % came from the GWAS summary statisitics data that have been
    % replicated
    loc_rsid = find( ismember( ...
        replicated_snps_afr.Properties.VariableNames,'IndexSNP')) ;
    replicated_snps_afr.Properties.VariableNames(loc_rsid) = ...
        "rsID_indSNP_AFR" ;
    
    % the "P" is the the column for replication P-values in AFR group
    % that was calculated using the "replicateVariant.m" function
    loc_rsid = find( ismember( ...
        replicated_snps_afr.Properties.VariableNames,'P')) ;
    replicated_snps_afr.Properties.VariableNames(loc_rsid) = ...
        "rep_p_in_AFR" ;
    
    % clean up the variable names. The "rsid" is the snps that came
    % from the snps in LD with the Index snps
    loc_rsid = find( ismember( ...
        replicated_snps_afr.Properties.VariableNames,'rsid')) ;
    replicated_snps_afr.Properties.VariableNames(loc_rsid) = "LinkedSNP" ;
    replicated_snps_afr = movevars( replicated_snps_afr, ...
        "rsID_indSNP_AFR" , "Before", 1) ;
    replicated_snps_afr = addvars( replicated_snps_afr, ...
        replicated_snps_afr.LinkedSNP, ...
        'NewVariableNames',...
        "LinkedSNP_in_AFR",'After',"rsID_indSNP_AFR") ;
    
    % **************************************************************************
    
    % check how many AFR variants are replicated from the original
    % variants and plot a a new venn diagram - I have to load the EUR
    % data and see how many snps are presents
    fprintf('\n Loading the EUR signficant variants for replication\n')
    eur_id_snps = readtable(['FUMA_multitrait_cardio_EUR_results/', ...
        'leadSNPs.txt'] ) ;
    loc_rsid = find( ismember( ...
        eur_id_snps.Properties.VariableNames,'rsID')) ;
    eur_id_snps.Properties.VariableNames(loc_rsid) = "rsID_indSNP_EUR";
    eur_id_snps.Properties.VariableNames(7) = "gwas_p_EUR" ;
    
    % *****************CONSIDER THE SNPS IN LD*************
    % **************************************************************
    % now let see how many of the EUR snps have been replicated
    replicated_snps_afr = innerjoin( ...
        eur_id_snps(:,{'rsID_indSNP_EUR','gwas_p_EUR'}), ...
        replicated_snps_afr, 'LeftKeys',{'rsID_indSNP_EUR'}, ...
        'RightKeys',{'LinkedSNP'}) ;
    
    % add the original p-values in AFR group from the GWAS summary
    % statistics
    fprintf('\n Loading all AFR GWAS STATS for replication\n')
    all_afr_gwas_stats = readtable(...
        'jass_results/gwas_all_AFR_4_variants_replication.txt', ...
        'FileType','text') ;
    loc_rsid = find( ismember( ...
        all_afr_gwas_stats.Properties.VariableNames,'P')) ;
    all_afr_gwas_stats.Properties.VariableNames(loc_rsid) = ...
        "orignal_p_in_AFR";
    
    % join with the replicated snps
    replicated_snps_afr = innerjoin( replicated_snps_afr ,  ...
        all_afr_gwas_stats(:, {'ID','orignal_p_in_AFR'}) , ...
        'LeftKeys',{'rsID_indSNP_AFR'}, 'RightKeys',{'ID'}) ;
    
    replicated_snps_afr = movevars( replicated_snps_afr,  ...
        {'LinkedSNP_in_AFR','Rsquare','orignal_p_in_AFR',...
        'rep_p_in_AFR'} , 'After',{'rsID_indSNP_AFR'} );
    
    % get only the variants that are significants
    replicated_snps_afr = replicated_snps_afr( ...
        replicated_snps_afr.rep_p_in_AFR < 0.1,:);
    
    % save to a file
    writetable(replicated_snps_afr,'replicated_snps_in_AFR_only.txt') ;
    
    % get the data for the venn diagram
    vennEURdata = readtable(['FUMA_multitrait_cardio_EUR_results/',...
        'leadSNPs.txt'] ) ;
    
    % replace the EUR snps Ids with the AFR snps ID in teh vennEURdata
    [~,Ai,Bi] = intersect(vennEURdata.rsID,  ...
        replicated_snps_afr.rsID_indSNP_EUR,'Stable') ;
    vennEURdata.rsID(Ai) = replicated_snps_afr.rsID_indSNP_AFR(Bi) ;
    
    % plot the venn diagram of SNPs unique in each population
    commonSnps = intersect( vennEURdata.rsID, ...
        replicated_snps_afr.rsID_indSNP_AFR);
    curEurSnps = setdiff( vennEURdata.rsID, ...
        replicated_snps_afr.rsID_indSNP_AFR);
    curAFRSnps = setdiff(replicated_snps_afr.rsID_indSNP_AFR, ...
        vennEURdata.rsID);
    
    % plot a venn diagram of the snps draw the venn diagram for the
    % SNPs that are intersecting
    figure()
    venn([200,200],75,'EdgeColor','black')
    hold on
    
    % add text to the figure: get numbers of intersects to add to the
    % plots
    myNumText(3) = length(curAFRSnps) ;
    myNumText(2) = length(commonSnps);
    myNumText(1) = length(curEurSnps) ;
    
    % here are text positions and the titles of the plots
    textPos = [0.25, 0.47 , 0.70] ;
    for jj = 1:length(textPos)
        text(textPos(jj),0.5,num2str(myNumText(jj)),...
            'Units','normalized',...
            'FontWeight','bold','FontSize',16,...
            'HorizontalAlignment','center')
    end
    
    % add the titles to the circle
    myGroups = {'EUR','AFR'};
    textPos = [0.33 , 0.66] ;
    for jj = 1:length(textPos)
        text(textPos(jj), 0.98 ,[myGroups{jj}],'Units','normalized',...
            'FontWeight','bold','FontSize',14,...
            'HorizontalAlignment','center')
    end
    hold off
    
    % save the figure
    saveas(gcf,'jass_results/Venn_Diagram_of_repAFR_SNPs.fig','fig')
    
    clear all_afr_gwas_stats vennEURdata
else
    % load the replciated variants
    replicated_snps_afr = readtable('replicated_snps_in_AFR_only.txt') ;
end

%% ********* Replicate the AFR variants in EUR group **************

% check if the variants have been replicated in EUR from EUR
if ~exist('replicated_snps_in_EUR.txt','file')
    
    % check that the ld data exist
    if ~exist('gwasStats_EUR','var')
        
        % check if the data has already been processed
        ldLinkData = readtable('heart_ldlink_data_EUR.txt') ;
        ldLinkData = movevars(ldLinkData,{'IndexSNP'},'Before',1) ;
        
        % get the EUR only snps and change the variable names for all the
        % snps and see which among those can be replicated
        gwasStats_EUR = readtable(...
            'jass_results/gwas_all_EUR_4_variants_replication.txt', ...
            'FileType','text') ;
        loc_rsid = find( ismember( ...
            gwasStats_EUR.Properties.VariableNames,'ID')) ;
        gwasStats_EUR.Properties.VariableNames(loc_rsid) = "rsid" ;
        
    end
    
    % here are the AFR sig snps that need to be replace
    afr_id_snps = readtable(['FUMA_multitrait_cardio_AFR_results/', ...
        'leadSNPs.txt'] ) ;
    
    % here are the replicated snps
    replicated_snps_eur  = [] ;
    
    % start the parallel pool 
    try
        % beging a parallel works
        if onHPC_Cluster == true
            % work on the entire datasets
            parpool('local',numOfCores)
        else
            error('\nThis section must be on a high memory computer\n')
        end
    catch
        % catch the error for the if statement 
        fprintf('\nParallel pool already running\n')
    end
    
    % let see what there 
    whos
    
    % find the variants in the intergrated data are associated with
    % the traits
    parfor jj = 1:height(afr_id_snps)
        
        % print something to the screen
        if rem(jj,20) == 0
            fprintf('\nRunning analysis for %s #%d of %d\n', ...
                afr_id_snps.rsID{jj}, jj, height(afr_id_snps) )
        end
        
        % here is the snps to be used for replication
        theSnp = afr_id_snps.rsID(jj) ;
        
        % replicate the results in the EUR group
        cur_rep = replicateVariant(theSnp,ldLinkData,gwasStats_EUR);
        
        % put in the same table
        replicated_snps_eur = [replicated_snps_eur ; cur_rep] ;
    end
    
    % some linked snps are '.' - remove those
    replicated_snps_eur( ...
        ismember( replicated_snps_eur.IndexSNP,'.'), :)  = [] ;
    
    % save the file
    writetable(replicated_snps_eur,'replicated_snps_in_EUR.txt');
    
    clear gwasStats_EUR
 
    % check that the snp data dose not exist for only the EURican variants
end


if ~exist('replicated_snps_in_EUR_only.txt','file')
    
    % load the file
    replicated_snps_eur = readtable('replicated_snps_in_EUR.txt') ;
    
    % ******************* replicate clean up ************************
    % change some variables names in the replicated snps data for the
    % EUR groups
    
    % the "IndexSNP" is the the column for snp Ids in EUR group that
    % came from the GWAS summary statisitics data that have been
    % replicated
    loc_rsid = find( ismember( ...
        replicated_snps_eur.Properties.VariableNames,'IndexSNP')) ;
    replicated_snps_eur.Properties.VariableNames(loc_rsid) = ...
        "rsID_indSNP_EUR" ;
    
    % the "P" is the the column for replication P-values in EUR group
    % that was calculated using the "replicateVariant.m" function
    loc_rsid = find( ismember( ...
        replicated_snps_eur.Properties.VariableNames,'P')) ;
    replicated_snps_eur.Properties.VariableNames(loc_rsid) = ...
        "rep_p_in_EUR" ;
    
    % clean up the variable names. The "rsid" is the snps that came
    % from the snps in LD with the Index snps
    loc_rsid = find( ismember( ...
        replicated_snps_eur.Properties.VariableNames,'rsid')) ;
    replicated_snps_eur.Properties.VariableNames(loc_rsid) = ...
        "LinkedSNP" ;
    replicated_snps_eur = movevars( replicated_snps_eur, ...
        "rsID_indSNP_EUR" , "Before", 1) ;
    replicated_snps_eur = addvars( replicated_snps_eur, ...
        replicated_snps_eur.LinkedSNP, ...
        'NewVariableNames',...
        "LinkedSNP_in_EUR",'After',"rsID_indSNP_EUR") ;
    
    % ***************************************************************
    
    % check how many EUR variants are replicated from the original
    % variants and plot a a new venn diagram - I have to load the EUR
    % data and see how many snps are presents
    fprintf('\nLoading the EUR signficant variants for replication\n')
    afr_id_snps = readtable(['FUMA_multitrait_cardio_AFR_results/', ...
        'leadSNPs.txt'] ) ;
    loc_rsid = find( ismember( ...
        afr_id_snps.Properties.VariableNames,'rsID')) ;
    afr_id_snps.Properties.VariableNames(loc_rsid) = "rsID_indSNP_AFR";
    afr_id_snps.Properties.VariableNames(7) = "gwas_p_AFR" ;
    
    % *****************CONSIDER THE SNPS IN LD*************!!!!!!!!!
    % **************************************************************
    % now let see how many of the EUR snps have been replicated
    replicated_snps_eur = innerjoin( ...
        afr_id_snps(:,{'rsID_indSNP_AFR','gwas_p_AFR'}), ...
        replicated_snps_eur, 'LeftKeys',{'rsID_indSNP_AFR'}, ...
        'RightKeys',{'LinkedSNP'}) ;
    
    % add the original p-values in EUR group from the GWAS summary
    % statistics
    fprintf('\n Loading all EUR GWAS STATS for replication\n')
    all_EUR_gwas_stats = readtable(...
        'jass_results/gwas_all_EUR_4_variants_replication.txt', ...
        'FileType','text') ;
    loc_rsid = find( ismember( ...
        all_EUR_gwas_stats.Properties.VariableNames,'P')) ;
    all_EUR_gwas_stats.Properties.VariableNames(loc_rsid) = ...
        "orignal_p_in_EUR";
    
    % join with the replicated snps
    replicated_snps_eur = innerjoin( replicated_snps_eur ,  ...
        all_EUR_gwas_stats(:, {'ID','orignal_p_in_EUR'}) , ...
        'LeftKeys',{'rsID_indSNP_EUR'}, 'RightKeys',{'ID'}) ;
    
    replicated_snps_eur = movevars( replicated_snps_eur,  ...
        {'LinkedSNP_in_EUR','Rsquare','orignal_p_in_EUR',...
        'rep_p_in_EUR'} , 'After',{'rsID_indSNP_EUR'} );
    
    % get only the variants that are significants
    replicated_snps_eur = replicated_snps_eur( ...
        replicated_snps_eur.rep_p_in_EUR < 0.1,:);
    
    % save to a file
    writetable(replicated_snps_eur,'replicated_snps_in_EUR_only.txt');
    
    % get the data for the venn diagram
    vennAFRdata = readtable(['FUMA_multitrait_cardio_AFR_results/',...
        'leadSNPs.txt'] ) ;
    
    % replace the EUR snps Ids with the EUR snps ID in teh vennEURdata
    [~,Ai,Bi] = intersect(vennAFRdata.rsID,  ...
        replicated_snps_eur.rsID_indSNP_AFR,'Stable') ;
    vennAFRdata.rsID(Ai) = replicated_snps_eur.rsID_indSNP_EUR(Bi) ;
    
    % plot the venn diagram of SNPs unique in each population
    commonSnps = intersect(vennAFRdata.rsID, ...
        replicated_snps_eur.rsID_indSNP_EUR);
    curEurSnps = setdiff(vennAFRdata.rsID, ...
        replicated_snps_eur.rsID_indSNP_EUR);
    curEURSnps = setdiff(replicated_snps_eur.rsID_indSNP_EUR, ...
        vennAFRdata.rsID);
    
    % plot a venn diagram of the snps draw the venn diagram for the
    % SNPs that are intersecting
    figure()
    venn([200,200],75,'EdgeColor','black')
    hold on
    
    % add text to the figure: get numbers of intersects to add to the
    % plots
    myNumText(3) = length(curEURSnps) ;
    myNumText(2) = length(commonSnps);
    myNumText(1) = length(curEurSnps) ;
    
    % here are text positions and the titles of the plots
    textPos = [0.25, 0.47 , 0.70] ;
    for jj = 1:length(textPos)
        text(textPos(jj),0.5,num2str(myNumText(jj)),...
            'Units','normalized',...
            'FontWeight','bold','FontSize',16,...
            'HorizontalAlignment','center')
    end
    
    % add the titles to the circle
    myGroups = {'AFR','EUR'};
    textPos = [0.33 , 0.66] ;
    for jj = 1:length(textPos)
        text(textPos(jj), 0.98 ,[myGroups{jj}],'Units','normalized',...
            'FontWeight','bold','FontSize',14,...
            'HorizontalAlignment','center')
    end
    hold off
    
    % save the figure
    saveas(gcf,'jass_results/Venn_Diagram_of_repEUR_SNPs.fig','fig')
    
    clear all_EUR_gwas_stats vennEURdata
else
    % load the replciated variants
    replicated_snps_eur = readtable('replicated_snps_in_EUR_only.txt');
end

% % ************** stop the code from excuting after processing the files 
% fprintf('\nI have stopped the code on line 1890\n')
% return

%% Check out many variants are replicated in each groups 

% here are the replicated variants
try
    % load the replicated snps
    replicated_snps_afr = readtable('replicated_snps_in_AFR_only.txt') ;
    replicated_snps_eur = readtable('replicated_snps_in_EUR_only.txt') ;
    
    % reduce the data to make get only the required variable s
    replicated_snps_afr = replicated_snps_afr(:,[1:7,8,12]) ;
    
    % reduce the data to make get only the required variable s
    replicated_snps_eur = replicated_snps_eur(:,[1:7,8,12] ) ;
    
    writetable(replicated_snps_afr, 'Supplementary Data 2.xlsx','Sheet', ...
        'Replicated SNPs in AFR')
    writetable(replicated_snps_eur, 'Supplementary Data 2.xlsx','Sheet', ...
        'Replicated SNPs in EUR')
catch
    
    % load the snps from a excel file
    replicated_snps_afr = readtable('Supplementary Data 2.xlsx','Sheet', ...
        'AFR LeadSNPs Rep In EUR') ;
    replicated_snps_eur = readtable('Supplementary Data 2.xlsx','Sheet', ...
        'EUR LeadSNPs Rep In AFR') ;
    
end

% load the snp freq data 
snpFreq = readtable('jass_results/snpFreq_all.txt');


%% --- 1. Pick the most significant replication for each lead SNP -------------

% AFR table: sort by LeadSNP_AFR then ascending LinkedSNP_EUR_Pvalue
topRep_AFR = sortrows(replicated_snps_afr, ...
            {'leadSNP_AFR','LinkedSNP_EUR_Pvalue'}, {'ascend','ascend'});

% EUR table: sort by LeadSNP_EUR then ascending LinkedSNP_AFR_Pvalue
topRep_EUR = sortrows(replicated_snps_eur, ...
            {'leadSNP_EUR','LinkedSNP_AFR_Pvalue'}, {'ascend','ascend'});

% --- 3. Outer join with frequency/location table ---
afrJoined = outerjoin(topRep_AFR, snpFreq, 'LeftKey','leadSNP_AFR', ...
    'RightKey','Variant','MergeKeys',true, 'Type','left');
eurJoined = outerjoin(topRep_EUR, snpFreq,'LeftKey','leadSNP_EUR', ... 
    'RightKey','Variant','MergeKeys',true, 'Type','left');
                  
% --- 5. Drop redundant columns and rename frequency variables ---
% AFR table
afrJoined.Properties.VariableNames{'EUR'} = 'EURfreq';
afrJoined.Properties.VariableNames{'AFR'} = 'AFRfreq';

% EUR table
eurJoined.Properties.VariableNames{'EUR'} = 'EURfreq';
eurJoined.Properties.VariableNames{'AFR'} = 'AFRfreq';

% --- 4. Save to workbook ---
writetable(afrJoined, 'Supplementary Data 2.xlsx', ...
           'Sheet','AFR LeadSNPs Rep In EUR','WriteMode','overwritesheet');

writetable(eurJoined, 'Supplementary Data 2.xlsx', ...
           'Sheet','EUR LeadSNPs Rep In AFR','WriteMode','overwritesheet');







%% Let find the top signficant snps between the gruops 

% Let check how many of those snps only signficant in one group or that
% could not be replicated had nan p-values
AFR = readtable('jass_results/afr_jass_with_beta.txt',...
    'Format', '%s%f%f%f%f%f%f', ...
    'Delimiter', ',') ; 
EUR = readtable('jass_results/eur_jass_with_beta.txt', ...
    'Format', '%s%f%f%f%f%f%f', ...
    'Delimiter', ',') ;

% here is the combined dataset 
allBoth = innerjoin(AFR, EUR, 'Key','SNP') ;

% here the indipendent snps in AFR group
afr_lead_snps = readtable(['FUMA_multitrait_cardio_AFR_results/', ...
            'leadSNPs.txt'] ) ;
        
% here are the replicated variants
replicated_snps_afr = readtable('replicated_snps_in_AFR_only.txt') ;

% the distinct loci in AFR group 
distinctAFR = afr_lead_snps( ~ismember( afr_lead_snps.rsID , ...
    [ replicated_snps_afr.rsID_indSNP_AFR; ...
    replicated_snps_afr.LinkedSNP_in_AFR]), :)  ;
distinctAFR = sortrows(distinctAFR, 'p','ascend') ;

% reduce the data to make get only the required variable s
replicated_snps_afr = replicated_snps_afr(:,1:7 ) ;

% here the indipendent snps in AFR group
eur_lead_snps = readtable(['FUMA_multitrait_cardio_EUR_results/', ...
            'leadSNPs.txt'] ) ;
        
% here are the replicated variants
replicated_snps_eur = readtable('replicated_snps_in_EUR_only.txt') ;

% reduce the data to make get only the required variable s
replicated_snps_eur = replicated_snps_eur(:,1:7 ) ;
        
% the distinct loci in AFR group 
distinctEUR = eur_lead_snps( ~ismember( eur_lead_snps.rsID , ...
    [ replicated_snps_eur.rsID_indSNP_EUR; ...
    replicated_snps_eur.LinkedSNP_in_EUR]), :)  ;
distinctEUR = sortrows(distinctEUR, 'p','ascend') ;

% now put together the top 5 snps in each group 
top5_afr = distinctAFR(1:5,:)  ;
top5_eur = distinctEUR(1:5,:)  ;

theTopSNPs = [top5_afr.rsID; top5_afr.rsID] ;
% theTopSNPs = allBoth( ismember( allBoth.SNP, theTopSNPs), :) ;

% get all the common variants between AFR and EUR 
commonVariant = [ replicated_snps_afr{:, ...
    {'rsID_indSNP_EUR','rsID_indSNP_AFR','LinkedSNP_in_AFR'}} ;
    replicated_snps_eur{:, ...
    {'rsID_indSNP_AFR','rsID_indSNP_EUR','LinkedSNP_in_EUR'} } ] ;
commonVariant = unique( commonVariant(:) )

% save the significant snps in EUR and AFR to the supplementary Data 
writetable( afr_lead_snps ,'Supplementary Data 2.xlsx', ...
    'Sheet', 'AFR - Lead SNPs') 
writetable( eur_lead_snps ,'Supplementary Data 2.xlsx', ...
    'Sheet', 'EUR - Lead SNPs') 


% here the indipendent snps in AFR group
afr_loci = readtable(['FUMA_multitrait_cardio_AFR_results/', ...
            'GenomicRiskLoci.txt'] ) ;
% here the indipendent snps in AFR group
eur_loci = readtable(['FUMA_multitrait_cardio_AFR_results/', ...
            'GenomicRiskLoci.txt'] ) ;
        
% save the significant snps in EUR and AFR to the supplementary Data 
writetable( afr_loci ,'Supplementary Data 2.xlsx', ...
    'Sheet', 'AFR - Genomic Risk Loci') 
writetable( eur_loci ,'Supplementary Data 2.xlsx', ...
    'Sheet', 'EUR - Genomic Risk Loci') 

clear afr_loci eur_loci

%% Now many among the AFR distict variants have nan p-values in EUR 

% ********************************************************************
% THIS IS THE CORRECT ANALYSIS NOT THAT IN THE PREVIOUS SECTION 

% here are the snps signficant in both groups 
afrSigAll = allBoth( allBoth.P_value_AFR < 5e-8 & ...
    ( allBoth.P_value_EUR < 5e-8 | isnan(allBoth.P_value_EUR ) ), :) ;
eurSigAll = allBoth( allBoth.P_value_EUR < 5e-8 & ...
    ( allBoth.P_value_AFR < 5e-8 | isnan(allBoth.P_value_AFR ) ), :) ;

% print something to the screen 
fprintf('\nThe number of Independent Sig SNPs in AFR = %d\n', ...
    height(afr_lead_snps))
fprintf('\nThe number of Independent Sig SNPs in EUR = %d\n', ...
    height(eur_lead_snps))

% here are the number of replicated snps in each population 
fprintf('\nThe number of replicated SNP in AFR = %d\n', ...
    length(unique( replicated_snps_afr.rsID_indSNP_EUR)))
fprintf('\nThe number of replicated SNPs in EUR = %d\n',...
    length(unique( replicated_snps_eur.rsID_indSNP_AFR)))

% print to the screen the number of SNPs with NaN pvalue in AFR
fprintf('\nThe number of AFR sig SNPs with NAN p in EUR is %d\n', ...
    sum( isnan( allBoth.P_value_EUR( ...
    ismember(allBoth.SNP, afr_lead_snps.rsID)))) )
fprintf('\nThe number of EUR sig SNPs with NAN p in AFR is %d\n', ...
    sum( isnan( allBoth.P_value_AFR( ...
    ismember( allBoth.SNP, eur_lead_snps.rsID) ) ) ) )

% save the results to supplementary data 
nanEUR = allBoth( isnan( allBoth.P_value_EUR( ...
    ismember(allBoth.SNP, afr_lead_snps.rsID))), : ) ;
nanAFR = allBoth( isnan( allBoth.P_value_AFR( ...
    ismember( allBoth.SNP, eur_lead_snps.rsID))), :) ;

% This is the correct common variant variable 
commonVariant = [ unique( replicated_snps_eur.rsID_indSNP_AFR) ; ...
   unique( replicated_snps_afr.rsID_indSNP_EUR) ]

% now here are the distinct snps that could not be replicated in either
% groups
fprintf(['\nThe number of common variants between AFR and EUR\n ', ...
    'account for the SNPs in LD with SNPs signficant in each group\n', ...
    ' is %d\n'], length(commonVariant) ) 

% now lets find the number of snps unique to each group 
% the distinct loci in AFR group 
distinctAFR = afr_lead_snps( ~ismember( afr_lead_snps.rsID , ...
    commonVariant), :)  ;
distinctAFR = unique(distinctAFR);
distinctAFR = sortrows( distinctAFR, 'p','ascend') ;

% the distinct loci in AFR group 
distinctEUR = eur_lead_snps( ~ismember(eur_lead_snps.rsID,commonVariant), :); 
distinctEUR = unique(distinctEUR);
distinctEUR = sortrows(distinctEUR, 'p','ascend') ;

% here are the number of SNPs unique to each group 
fprintf('\nThe number of AFR specific SNPs = %d\n',height(distinctAFR))
fprintf('\nThe number of EUR specific SNPs = %d\n',height(distinctEUR))

% plot a venn diagram of the snps draw the venn diagram for the
% SNPs that are intersecting
figure()
venn([200,200],75,'EdgeColor','black')
hold on

% add text to the figure: get numbers of intersects to add to the
% plots
myNumText(1) = height(distinctAFR) ;
myNumText(2) = length(commonVariant); % 79
myNumText(3) = height(distinctEUR) ;

% here are text positions and the titles of the plots
textPos = [0.25, 0.47 , 0.70] ;
for jj = 1:length(textPos)
    text(textPos(jj),0.5,num2str(myNumText(jj)),...
        'Units','normalized',...
        'FontWeight','bold','FontSize',16,...
        'HorizontalAlignment','center')
end

% add the titles to the circle
myGroups = {'AFR','EUR'};
textPos = [0.33 , 0.66] ;
for jj = 1:length(textPos)
    text(textPos(jj), 0.98 ,[myGroups{jj}],'Units','normalized',...
        'FontWeight','bold','FontSize',14,...
        'HorizontalAlignment','center')
end
hold off

% save the figure
saveas(gcf,'jass_results/Venn_Diagram_of_CommonSNP_Real.fig','fig')

% here are the snps signficant in both groups 
afrSig = allBoth( allBoth.P_value_AFR < 5e-8 , :) ;
eurSig = allBoth( allBoth.P_value_EUR < 5e-8 , :) ;

% get the afr and eur sig snps are are part of the independence snsp 
afrSig = afrSig( ismember(afrSig.SNP, afr_lead_snps.rsID), :) ;
eurSig = eurSig( ismember(eurSig.SNP, eur_lead_snps.rsID), :) ;

% here are the top 5 snps: these will be obained from the distinctAFR and
% distinctEUR table and then extracted from the both Sig SNPs table 
% here are the fuma results of the snps variannt 
afrSig = sortrows(afrSig,'P_value_AFR','ascend') ;
eurSig = sortrows(eurSig,'P_value_EUR','ascend') ;
top5both = allBoth( ismember( allBoth.SNP, ...
    [afrSig.SNP(1:5) ; eurSig.SNP(1:5)] ), :) ;

% add the gene names to the table and positions to the data 
posData = [afr_lead_snps ;eur_lead_snps ] ;
top5both = innerjoin( posData(:,{'rsID','chr','pos'}),...
    top5both, 'LeftKey',{'rsID'},'RightKey',{'SNP'})  

% also add the genes symbols to the data 
eurGenes = readtable('FUMA_multitrait_cardio_EUR_results/snps.txt');
afrGenes = readtable('FUMA_multitrait_cardio_AFR_results/snps.txt') ;
eurGenes = unique( eurGenes(:,{'IndSigSNP','rsID','nearestGene','dist'})) ;
afrGenes = unique( afrGenes(:,{'IndSigSNP','rsID','nearestGene','dist'} ));
theGenes = [ eurGenes; afrGenes] ;
top5both.Properties.VariableNames(1) = {'rsID'} ;

% add the genes to the table 
for ii = 1:height(top5both)
    
   % get the genes
   curGenes = theGenes( ...
       ismember( theGenes.IndSigSNP, top5both.rsID(ii) ) ,:) ;
   
   % sort the table to return the nearest genes 
   curGenes = sortrows( curGenes, 'dist','ascend') ;
   
   % add the first gene to the table 
   top5both.nearestGene(ii) = curGenes.nearestGene(1) ; 
    
end

% remove the variable that are not reqirued 
top5both = removevars( top5both, ...
    {'SE_AFR','beta_AFR','varbeta_AFR','SE_EUR','beta_EUR','varbeta_EUR'});
top5both = movevars( top5both, 'nearestGene','After','pos') 

% save the data to a file
writetable(top5both ,'Table1.xlsx')

clear posData clear 

% % ************** stop the code from excuting after processing the files
% fprintf('\nI have stopped the code on line 744\n')
% return

%% load the UK bio bank that has SNP frequencies

% create the data for the analysis 
if exist('tableau_cardio_plots.csv','file')
    
    % load the SNP freq data that is processed
    fprintf('\nLoading the processed SNP freq data for tableau \n')
    tableauData = readtable('tableau_cardio_plots.csv') ;
    
else
    % process the SNP freq data
    fprintf('\nProcessing SNP freq data \n')
    
    % change this the snps that are present in the UK Biobank from the data
    % which also included the imputed snps
    tableauData = readtable('full_variant_qc_metrics.txt');
    
    % return only the required table variables
    tableauData = tableauData(:,{'chrom','pos','alt','ref','rsid',...
        'af_EUR','af_AFR'});
    
    % change the variables names to those present in the snp comparison
    % table
    tableauData.Properties.VariableNames = {'Chrom','Position',...
        'Alternative','Reference','Variant','EUR','AFR'};
    
    % read in the data of the current condition in
    curAFR = readtable('jass_results/afr_multitrait_results_1.tbl',...
        'FileType','text') ;
    curEUR = readtable('jass_results/eur_multitrait_results_1.tbl',...
        'FileType','text') ;
    
    % throw in an assertion 
    assert( all( strcmp( curAFR.SNP_ID , curEUR.SNP_ID ) ))
    
    % get all the SNPs that are signicant in AFR group only and reutl only
    % those for from the data 
    % ********************************************************************
    snpFreqSig = readtable('jass_results/snpFreq_indSig_snps.csv') ;
    
    % merge with EUR and AFR original data 
    toAddGWAS = curAFR(:, {'SNP_ID','Zscore','P_value','Direction'});
    toAddGWAS.Properties.VariableNames(2:end) = strcat( ...
         'AFR_', toAddGWAS.Properties.VariableNames(2:end) ) ;
     
    % add that to the table
    snpFreqSig = innerjoin( snpFreqSig, toAddGWAS, 'LeftKeys', ...
        {'Variant'},'RightKey',{'SNP_ID'}) ;
    
    % merge with EUR and AFR original data
    toAddGWAS = curEUR(:, {'SNP_ID','Zscore','P_value','Direction'});
    toAddGWAS.Properties.VariableNames(2:end) = strcat( ...
        'EUR_', toAddGWAS.Properties.VariableNames(2:end) ) ;
    
    % add that to the table
    snpFreqSig = innerjoin( snpFreqSig, toAddGWAS, 'LeftKeys', ...
        {'Variant'},'RightKey',{'SNP_ID'} ) ;
   
    % save to a file
    writetable(snpFreqSig,'jass_results/snpFreq_indSig_snps_with_ps.csv') 
    
    % *******************************************************************
    
    % get the snps significant at the suggestive p-value threshold
    toKeep = curAFR.P_value > 0.01 | curEUR.P_value > 0.01 ;
    curAFR(toKeep, :) = [];
    curEUR(toKeep, :) = [];
    
    % add those condition to the table where the SNPs is true
    locAFR = ismember(tableauData.Variant, curAFR.SNP_ID) ;
    locEUR = ismember(tableauData.Variant, curEUR.SNP_ID) ;
    
    % add the variable to the table
    tableauData.('gwas_P_AFR') = locAFR ;
    tableauData.('gwas_P_EUR') = locEUR ;
    
    % add a new variable at the end of the column to show which SNPs
    % related to EUR and those related to AFR
    tableauData.SigEUR = any(tableauData{:,{'gwas_P_EUR'}},2) ;
    tableauData.SigAFR = any(tableauData{:,{'gwas_P_AFR'}},2) ;
    
    % get instance where the gwas is significant in EUR
    tableauData.SigGroup(tableauData.SigEUR==true) = {'EUR'} ;
    tableauData.SigGroup(tableauData.SigAFR==true) = {'AFR'} ;
    
    % get the instances where both snps are significant
    bothSig = all([tableauData.SigEUR,tableauData.SigAFR],2) ;
    tableauData.SigGroup(bothSig)= {'Both'};
    
    % sort the rows of the data
    tableauData = sortrows(tableauData,'EUR','descend');
    
    % save the data to a csv file
    fprintf('\nSaving the tableau data for plotting\n')
    writetable(tableauData,'tableau_cardio_plots.csv')  
end

% % ************** stop the code from excuting after processing the files
% fprintf('\nI have stopped the code on line 744\n')
% return



%% Produce Locus Zoom plots for all AFR sig SNPs traits 

% the plot show have AFR and bottom plot have EUR, the points must be
% coloured based on the LD of the variants to the lead variants and this
% information is present the ld_link_data. I have also created a function
% using chat GPT to create regional plots

% here is the gwas data with only the required variable returns
gwas_data = readtable('jass_results/fuma_multitrait_cardio_AFR.txt') ;
gwas_data = gwas_data(:,{'ID','CHR','POS','P'} ) ;
gwas_data.Properties.VariableNames = ["SNP_ID", "CHR", "POS", "P"];

% get the recombination rate and return only the required variable 
recomb_rate_data = readtable('genetic_map_GRCh37.txt') ;
recomb_rate_data = recomb_rate_data(:,{'Position_bp_','Rate_cM_Mb_'}) ;
recomb_rate_data.Properties.VariableNames = ["pos","rate"] ;

% get the data for the genes position
gene_pos = readtable('mart_export_with_TSS.txt'); 
gene_pos = gene_pos(:, {'GeneStart_bp_','GeneEnd_bp_', ...
    'Chromosome_scaffoldName','GeneName'} ) ; 
gene_pos.Properties.VariableNames =  ...
    {'start_pos','end_pos','chr','gene_name'} ;

% here the indipendent snps in AFR group
afr_lead_snps = readtable(['FUMA_multitrait_cardio_AFR_results/', ...
            'IndSigSNPs.txt'] ) ;
        
% % try to use LD data for AFR
% ldLinkData,'heart_ldlink_data_AFR.txt'
ld_data2 = readtable('heart_ldlink_data_AFR.txt') ;

%% Functional Impact of SNPs 

% plot the enrichment scores of the data
theGroups = {'AFR','EUR'} ;

figure()
tiledlayout(2,1,'TileSpacing','compact')

% here are the letter 
theLetters = ['a','b'];

% plot the bargraph in a loop 
for jj = 1:length(theGroups)
    
    % Load your data
    curData = readtable( ...
        sprintf('FUMA_multitrait_cardio_%s_results/annov.stats.txt',  ...
        theGroups{jj} ) ) ;
    
    fprintf('\nHere the functional enrichment data for %s\n',theGroups{jj})
    curData = sortrows( curData, 'prop','descend')
    
    % Get the annotations, proportions, enrichment scores, and p-values
    annot = curData.annot;
    props = curData.prop;
    enrich = curData.enrichment+0.1;
    p_vals = curData.fisher_P;
    
    % Convert enrichment scores to log2 scale
    log_enrich = log2(enrich) ;
    
    % Generate a colormap based on the enrichment scores
    colors = cbrewer('div','RdBu',length(enrich)) ;
    [~, sort_order] = sortrows(log_enrich,'descend') ;
    
    % here is the next tile 
    nexttile
    
    % Create the bar graph
    hBar = bar(props);
    hBar.FaceColor = 'flat';
    
    % change the colours of the bars
    for ii = 1:length(enrich)
        hBar.CData(sort_order(ii),:) = colors(ii,:);
    end
    
    % Label the axes
    xlabel('Annotation');
    ylabel('Proportion');
    set(gca,'FontSize',12,'LineWidth',1,'Box','off','TickDir','out')
    title(theGroups{jj})
    
    % Rotate x-axis labels to prevent them from overlapping
    xticks(1:length(annot));
    xticklabels( replace(annot , '_',' ') ) ;
    xtickangle(45);
     
    % Add significance indicators on top of the bars
    for i = 1:length(p_vals)
        % adjust as needed to position the indicators above the bars
        y = props(i) + 0.005;
        if p_vals(i) < 0.001
            text(i, y, '***','HorizontalAlignment','center','FontSize',18);
        elseif p_vals(i) < 0.01
            text(i, y, '**', 'HorizontalAlignment','center','FontSize',18);
        elseif p_vals(i) < 0.05
            text(i, y, '*', 'HorizontalAlignment','center','FontSize',18);
        end
    end
    
    % add a letter to the figure and remove the letter from the array
    text(-0.10, 1 ,theLetters(1),'Units','normalized',...
        'FontWeight','bold','FontSize',24)
    theLetters(1) = [] ;
     
end

% % Save the figure
% saveas(gcf, 'bar_graph.png');
%%

% Loop over the groups to create pie charts
for jj = 1:length(theGroups)
    
    % Start a new figure for the pie charts
    figure()

    % Load your data
    curData = readtable( ...
        sprintf('FUMA_multitrait_cardio_%s_results/annov.stats.txt',  ...
        theGroups{jj} ) ) ;
    
    fprintf('\nHere the functional enrichment data for %s\n',theGroups{jj})
    curData = sortrows( curData, 'prop','descend');
    
    % Get the annotations, proportions, enrichment scores, and p-values
    annot = replace(curData.annot , '_',' ');
    props = curData.prop;
    enrich = curData.enrichment+0.1;
    p_vals = curData.fisher_P;
    
    % Convert enrichment scores to log2 scale
    log_enrich = log2(enrich) ;
    
    % Generate a colormap based on the enrichment scores
    colors = cbrewer('div','RdBu',length(enrich)) ;
    [~, sort_order] = sortrows(log_enrich,'descend') ;
        
    % Create the pie chart
    % find which categories have proportions less than 0.05
    exploded = props < 0.05;
    % use the 'colors' array for the pie chart colors
    hPie = pie(props, exploded, annot);
    
    % Apply the colors to the pie chart
    for k = 1:length(props)
        hPie(2*sort_order(k)-1).FaceColor = colors(k,:);
    end
    
    % Remove labels from exploded parts
    exploded_indices = find(exploded);
    for k = 1:length(exploded_indices)
        hPie(2*(exploded_indices(k))).String = '';
    end
    
    % Label the axes and the title
    title(theGroups{jj})
    
    % Save the figure
    saveas(gcf, [theGroups{jj} '_pie_graphs.png']);
    
end


%% Here the SNPs that are signficance -- SNP Freq analysis

% load the snp freq data 
snpFreq = readtable('jass_results/snpFreq_indSig_snps.csv');
% here are the signficant snp freq data 

fprintf('\The is of the snp frequency data is %d\n', size(snpFreq) )

fprintf('\nHere are the SNPs signicant between groups\n')
summary(categorical(snpFreq.SigGroup))

% here are the first 20 rows on the data 
snpFreq(1:20,:)

afrSnpFreq = snpFreq( ismember( snpFreq.SigGroup , 'AFR'), :) 


%% Magma Gene Analysis

% plot the graph for mapped genes
eur_mapped_genes = readtable(['FUMA_multitrait_cardio_EUR_results/',...
    'genes.txt']) ;
afr_mapped_genes = readtable(['FUMA_multitrait_cardio_AFR_results/',...
    'genes.txt']) ;

% combine the results in a table
eur_mapped_genes = addvars( eur_mapped_genes, ...
    repmat({'EUR'}, height(eur_mapped_genes), 1), ...
    'NewVariableNames',{'Ethnicity'} ,'Before',1) ;
size(eur_mapped_genes)

afr_mapped_genes = addvars( afr_mapped_genes, ...
    repmat({'AFR'}, height(afr_mapped_genes), 1), ...
    'NewVariableNames',{'Ethnicity'} ,'Before',1) ;
size(afr_mapped_genes)

% plot the venn diagrams for the genes plot the venn diagram of
% SNPs unique in each population
commonGenes = unique( intersect(...
    eur_mapped_genes.symbol, afr_mapped_genes.symbol));
onlyEurMappedGenes = unique( setdiff(...
    eur_mapped_genes.symbol, afr_mapped_genes.symbol));
onlyAFRMappedgenes = unique( setdiff(...
    afr_mapped_genes.symbol, eur_mapped_genes.symbol));

% plot a venn diagram of the snps draw the venn diagram for the
% SNPs that are intersecting
figure()
venn([200,200],75,'EdgeColor','black')
hold on

% add text to the figure: get numbers of intersects to add to the
% plots
myNumText(3) = length(onlyAFRMappedgenes) 
myNumText(2) = length(commonGenes)
myNumText(1) = length(onlyEurMappedGenes) 

% here are text positions and the titles of the plots
textPos = [0.25, 0.47 , 0.70] ;
for jj = 1:length(textPos)
    text(textPos(jj),0.5,num2str(myNumText(jj)),'Units', ...
        'normalized','FontWeight','bold','FontSize',16,...
        'HorizontalAlignment','center')
end

% add the titles to the circle
myGroups = {'EUR','AFR'};
textPos = [0.33 , 0.66] ;
for jj = 1:length(textPos)
    text(textPos(jj), 0.98 ,[myGroups{jj}], ...
        'Units','normalized','FontWeight','bold','FontSize',14,...
        'HorizontalAlignment','center')
end
hold off

% save the figure
saveas(gcf,'jass_results/Venn_Diagram_of_MappedGenes.fig','fig')

% also save the list of genes for the phenotype to an excel fiel
writetable( eur_mapped_genes,'jass_results/each_mapped_genes.xlsx', ...
    'Sheet','EUR')
writetable( afr_mapped_genes,'jass_results/each_mapped_genes.xlsx', ...
    'Sheet','AFR')

% also save the list of genes for the phenotype to an excel fiel
writetable( eur_mapped_genes,'Supplementary Data 5.xlsx', ...
    'Sheet','Mapped Genes - EUR')
writetable( afr_mapped_genes,'Supplementary Data 5.xlsx', ...
    'Sheet','Mapped Genes - AFR') 

%%

for ii = 25 % height(gwas_data)
    
    % get the current snp ID 
    SNP_ID = afr_lead_snps.rsID{ii} ;
    
    % load the LD data for the current snp and get only the required
    % columns and name them appropriately for the funciton
    try
        LD_data = readtable( ['lead_snp_files_AFR/',SNP_ID,'.txt'] ) ;
        LD_data = addvars( LD_data(:,{'RS_Number','R2','Coord'}) , ...
            repmat( afr_lead_snps.rsID(ii), height(LD_data), 1) ,  ...
            'Before',1,'NewVariableNames',{'SNP_ID1'}) ;   
    catch
        continue
    end
   
    % add the postion to the data and remove the unrequried and rename the
    % variable in the table 
    LD_data.POS = str2double( extractAfter( LD_data.Coord,':') ) ;
    LD_data.Coord = [] ;
    LD_data.Properties.VariableNames = ["SNP_ID1","SNP_ID2","r2","POS"];
    
%     LD_data = ld_data2( ismember( ld_data2.IndexSNP, SNP_ID), 1:4) ;
%     LD_data.Properties.VariableNames = ["SNP_ID1", "SNP_ID2", "r2", "POS"];
%     LD_data( ismember( LD_data.SNP_ID1, '.') |  ...
%         ismember( LD_data.SNP_ID2, '.'), :) = [] ;
    
    % now call the regional plot function 
    plot_regional(gwas_data, recomb_rate_data, gene_pos, LD_data, SNP_ID)
    
end

%% Report the indendent SNPs and in each groups 

% one whould be for africans and the other from europeans
ukGwasFilesNames = {'PulseRate.txt';'SystolicBP.txt';'DiastolicBP.txt';...
    'MaxHeartRate.txt'} ;

% change my vars to the compact form
myVars = extractBefore(ukGwasFilesNames,'.txt') ;

% here are the the groups
theGroups = {'AFR','EUR'} ;

% here are the results for one traits 
oneTraitResults = table('Size',[0,3], ...
    'VariableTypes',{'cell','cell','double'}, ...
    'VariableNames',{'Group','Measure','NumIndSnps'} ) ;

for jj = 1:length( theGroups)
    
    for ii = 1:length(myVars)
        
        % here is the current folder
        curFolder = sprintf('FUMA_%s_%s', myVars{ii}, theGroups{jj}) ;
        
        % here are the current results
        try
            % load the file
            curResults = readtable(['Old_Results/fuma_results/', ...
                curFolder, '/IndSigSNPs.txt' ])  ;
            
            % here are the current snps
            curSnps = addvars( array2table( ...
                [ theGroups(jj), myVars(ii) ]) , ...
                sum(curResults.p < 5e-8)  ) ;
            curSnps.Properties.VariableNames = ...
                {'Group','Measure','NumIndSnps'} ;
        catch 
            % here are the current snps
            curSnps = addvars(array2table( [theGroups(jj), myVars(ii)] )...
                ,0  ) ;
            curSnps.Properties.VariableNames = ...
                {'Group','Measure','NumIndSnps'} ;
        end
        
        % get count the number of Snps involved
        oneTraitResults = [ oneTraitResults; curSnps] ;
        
        % remove the variable 
        curSnps = [] ;
    end
end

%% here are the unique snps 

% here the indipendent snps in AFR group
afr_lead_snps = readtable(['FUMA_multitrait_cardio_AFR_results/', ...
            'IndSigSNPs.txt'] ) ;
        
% here are the replicated variants
replicated_snps_afr = readtable('replicated_snps_in_AFR_only.txt') ;

% reduce the data to make get only the required variable s
replicated_snps_afr = replicated_snps_afr(:,1:7 ) ;
        
% the distinct loci in AFR group 
distinctAFR = afr_lead_snps( ~ismember( afr_lead_snps.rsID , ...
    [ replicated_snps_afr.rsID_indSNP_AFR; ...
    replicated_snps_afr.LinkedSNP_in_AFR]), :)  ;

distinctAFR = sortrows( distinctAFR, 'p','ascend') ;

% here the indipendent snps in AFR group
eur_lead_snps = readtable(['FUMA_multitrait_cardio_EUR_results/', ...
            'IndSigSNPs.txt'] ) ;
        
% here are the replicated variants
replicated_snps_eur = readtable('replicated_snps_in_EUR_only.txt') ;

% reduce the data to make get only the required variable s
replicated_snps_eur = replicated_snps_eur(:,1:7 ) ;
        
% the distinct loci in AFR group 
distinctEUR = eur_lead_snps( ~ismember( eur_lead_snps.rsID , ...
    [ replicated_snps_eur.rsID_indSNP_EUR; ...
    replicated_snps_eur.LinkedSNP_in_EUR]), :)  ;

distinctEUR = sortrows(distinctEUR, 'p','ascend') ;

% now put together the top 5 snps in each group 
top5_afr = distinctAFR(1:5,:)  ;
top5_eur = distinctEUR(1:5,:)  ;

theTopSNPs = [top5_afr.rsID; top5_afr.rsID] ;

%% Create a table for the Top5 Snps in each group 

% I processed the usign teh extractTop10Snps.sh script
if ~exist('jass_results/eur_multitrait_topresults.txt','file')
    % process the data 
    afr_multi = readtable('jass_results/fuma_multitrait_cardio_AFR.txt',...
          'FileType','text' );
    eur_multi = readtable('jass_results/fuma_multitrait_cardio_EUR.txt',...
          'FileType','text');
    afr_multi = afr_multi(ismember(afr_multi.ID, theTopSNPs), :);
    eur_multi = eur_multi(ismember(eur_multi.ID, theTopSNPs), :);
    
    % save the data
    writetable(afr_multi,'jass_results/afr_multitrait_topresults.txt')
    writetable(eur_multi,'jass_results/eur_multitrait_topresults.txt')
end

% load the data 
eur_top5 = readtable('jass_results/eur_multitrait_topresults.txt') ;
afr_top5 = readtable('jass_results/afr_multitrait_topresults.txt') ;


%% New compare the enrichement analysis between AFR and EUR 

% the colors to use the for two groups
groupColors = [ 0.47,0.67,0.19; 0.85,0.33,0.10  ] ;

myPath = ['/Users/sinkala/Documents/MATLAB/UKBioBank/', ...
    'cardiovascularNewAnalysis/'];

% set up the genes 
enrichrVars = {'ClinVar_2019_table','DisGeNET_table',...
    'UK_Biobank_GWAS_v1_table','WikiPathway_2021_Human_table'} ;

% set up the plot names
plotNames = {'ClinVar','DisGeNET','UK Biobank GWAS','WikiPathway'} ;

warning('off')

% plot the graphs in a loop 
for ii = 1:length(enrichrVars)
    
    % % load the data
    enrichAFR = readtable([myPath, enrichrVars{ii}, '-African.txt'],...
        'Format','auto'); 
    enrichEUR = readtable([myPath, enrichrVars{ii}, '-European.txt'],...
        'Format','auto'); 
    
    fprintf('\nThese are the top 5 signficant pathways for %s\n', ...
        plotNames{ii} );
    
    % let see the head of these files 
    enrichAFR_top5 = enrichAFR(1:5,1:8)
    enrichEUR_top5 = enrichEUR(1:5,1:8)
    
    % clean up the enriched term names
    enrichEUR.Term = strtrim( regexprep(enrichEUR.Term, ...
        {'\(+\w*','\:+\w*',')','\R-HSA-+\w*','Homo sapiens','\w* raw'},''));
    enrichAFR.Term = strtrim( regexprep(enrichAFR.Term, ...
        {'\(+\w*','\:+\w*',')','\R-HSA-+\w*','Homo sapiens','\w* raw'},''));
    
    % create a plot and also produce the go biological data to use for the
    % venn diagrams
   enrichrPlot_pValue(enrichAFR, enrichEUR,plotNames{ii} ,...
        flipud(groupColors),{'AFR','EUR'} ) ;
    
    % save to a supplementary file 
    writetable(enrichEUR,'Supplementary Data 5.xlsx',...
        'Sheet',[plotNames{ii},'-EUR'])
    writetable(enrichAFR,'Supplementary Data 5.xlsx',...
        'Sheet',[plotNames{ii},'-AFR'])
    
end

clear enrichEUR  enrichAFR ii plotNames enrichVars

     
%% Apply a chi square test to the SNP frequencies

% load the SNP frequency data 
snpFreq = readtable('jass_results/snpFreq_lead_snps.csv') ;

% preallocate the chi-square table 
chiTable = snpFreq(:,{'Chrom','Position','Alternative','Reference',...
    'Variant','EUR','AFR','SigGroup'}) ;

% run the data in a loop
for ii = 1:height(chiTable)

    % print something to the screen
    if rem(ii,50) == 0
        fprintf('\nRunning chi-square test # %d of %d Variables \n',...
            ii, height(chiTable))
    end
    
    % convert the variables to a categorical array and if it
    % fails then continue
    curData = chiTable(ii,{'AFR','EUR'}) ;
    curData = curData{1,:}' ;
    
    % convert the nan to zero
    if any(isnan(curData))
        curData(isnan(curData)) = 0;
    end
    
    % get the data for the current fisher exact test 
    fisherData = table('Size',[2,2],'VariableTypes',{'double','double'},...
        'VariableNames',{'snp','Nosnp'} ) ;
    
    % add the snp to the fisher exact test table and delete the count
    % variable that is not need for the fisher exact test
    fisherData.snp = round([383471;5978].*curData)  ;
    fisherData.Nosnp = round([383471;5978] - fisherData.snp) ;
    
    % now perform the chi-square test
    [~, p, stats ] = fishertest(fisherData);
    
    % add the results to the table
    % add the variable name to the table
    chiTable.pValue(ii) = p ; 
    chiTable.OddRatio(ii) = stats.OddsRatio ; 
    chiTable.LowerBound(ii) = stats.ConfidenceInterval(1);
    chiTable.UpperBound(ii) = stats.ConfidenceInterval(2);
    
end  

% add the variable names to the table and then remove the rows with missing
% results and add the adjusted p value to the table 
chiResults = addvars(chiTable, ...
    mafdr(chiTable.pValue,'BHFDR',true),'After','pValue', ...
    'NewVariableNames','FDR') ;
chiResults = sortrows(chiResults,'pValue','ascend');

% get only the unique snps
[~,theUnique ]= unique(chiResults.Variant) ;
chiResults = chiResults(theUnique, :)  ;

% save the results to excel 
writetable(chiResults,'jass_results/snp_comparisons_imputed.csv')
writetable(chiResults, 'Supplementary Data 3.xlsx', ...
    'Sheet','SNP Frequency Comp') 

clear ii curData tbl chi2 chiPvalue labels ci p stats meanAfricans  ...
    meanWhites chiTable ttestTable chiData xValues groups plotName ...
    curCol locFreq africaChiTable

%% Get the SNPs Whose Frequencies varies the most 

% load the chi square results
chiResults  = readtable('jass_results/snp_comparisons_imputed.csv');

% get only the sig SNPs
mostSigSNPfreq = chiResults(:,{'Variant','EUR','AFR','SigGroup',...
    'pValue','FDR','OddRatio','LowerBound','UpperBound'}) ;

% get the unique snps
[~,theUnique] = unique(mostSigSNPfreq.Variant);
mostSigSNPfreq = mostSigSNPfreq(theUnique, :);

% find the difference in the snp freq
mostSigSNPfreq = addvars(mostSigSNPfreq, ...
    mostSigSNPfreq.EUR - mostSigSNPfreq.AFR,'NewVariableNames',...
    "freqDiff",'After','AFR') ;

% return the variants that are present in the data 
% ********************************************************************
eurGenes = readtable('FUMA_multitrait_cardio_EUR_results/snps.txt');
afrGenes = readtable('FUMA_multitrait_cardio_AFR_results/snps.txt') ;
eurGenes = unique( eurGenes(:,{'IndSigSNP','rsID'})) ;
afrGenes = unique( afrGenes(:,{'IndSigSNP','rsID'} ));
mySigSnps = [eurGenes.IndSigSNP ; eurGenes.rsID;  ...
    afrGenes.IndSigSNP ; afrGenes.rsID] ;

% now subset teh table 
mostSigSNPfreq = mostSigSNPfreq( ismember( mostSigSNPfreq.Variant , ...
    mySigSnps), :) ;

% ********************************************************************

% sort the table 
mostSigSNPfreq = sortrows(mostSigSNPfreq,'freqDiff','descend');

% plot two figures of the differences get the data and sort according the
% frequency of the snps in EUR 
plotData = mostSigSNPfreq([1:8,end-8:end],{'Variant','EUR','AFR',...
    'freqDiff','pValue'});
plotData.Variant = categorical(plotData.Variant) ;
plotData = sortrows(plotData,'freqDiff','ascend');

% the colors to use the for two groups
groupColors = [ 0.47,0.67,0.19; 0.85,0.33,0.10  ] ;

% produce a bar graph
figure()
bar1 = bar(plotData.Variant, plotData{:,{'EUR','AFR'}});

% change the properties of the bars graphs
set(bar1(2),'FaceColor',groupColors(1,:),'EdgeColor',groupColors(1,:));
set(bar1(1),'FaceColor',groupColors(2,:),'EdgeColor',groupColors(2,:));

set(gca,'FontSize',14,'LineWidth',2,'Box','off','TickDir','out')
ylabel('SNP Frequency')
title('SNP Frequency Comparison','FontSize',16,'FontWeight','bold')
legend({'EUR','AFR'})

% save the figure
saveas(gcf,'SNP_frequency_comparison.fig','fig')

clear plotData theUnique

%% Find the SNP that vary the most between the Groups 

chiResults  = readtable('snp_comparisons_imputed.csv');
snps_func_impact = chiResults  ;

%%

snps_func_impact = addvars( snps_func_impact, ...
    snps_func_impact.AFR./snps_func_impact.EUR , ...
    'After','AFR','NewVariableNames',{'freqDiv'});
snps_func_impact = sortrows(snps_func_impact, 'freqDiv','descend') ;


%% ***************** Some Internal Function *********************

function pSuperScript = convertPValue2SuperScript(p)
    % converts the number to scientific superscript for printing on a
    % figure
    pS = num2str(p) ;

    % get the first number
    firstNumbers = extractBefore(pS,'e') ;

    % check if there is a decimal place. then only get the first 4 numbers
    if contains( firstNumbers  ,'.')
        firstNumbers = firstNumbers(1:4) ;
    end

    % get the correctly formated p value
    pSuperScript = sprintf('= %s x 10^{%d}', firstNumbers, ...
        str2double(extractAfter(pS, 'e') )) ;

    % if the p value is large
    if p > 0.0001
        pSuperScript = sprintf('= %0.4f', p) ;
    elseif p == 0
        pSuperScript = sprintf('< 1 x 10^{%d}', -350) ;
    end

end


function addReferenceLineToPlot(xVariable, yVariable)
% This function add a reference line to the plot and the R square value

% replace the missing values with mean
xVariable = fillmissing(xVariable,"linear");
yVariable = fillmissing(yVariable,"linear");

% preform a regression calculation and add it to the plot
X = [ones(length(xVariable),1) xVariable] ;
b1 = X\yVariable    ;     yCalc1 = X*b1 ;
plot(xVariable,yCalc1,'LineWidth',1.5,'Color','k')

% calculate the pearson's linear correation  and add it to the plot
[r2 , pR ] = corr(xVariable,yVariable, 'Type','Pearson',...
    'Rows','complete');
if pR < 0.0001
    text(0.4, 0.9 , strcat( sprintf( ...
        "R = %0.2f, P ", r2), convertPValue2SuperScript(pR)  ), ...
        'Units','normalized','FontSize',11,'FontWeight','bold')
else
    text(0.4, 0.9 ,sprintf('R = %0.2f, P = %0.4f', r2, pR), ...
        'Units','normalized','FontSize',11,'FontWeight','bold')
end
end

%% ====================== Internal Functions =================

function createLegendInternal(yPoint, xStart, legendLabels , plotColors,...
    myLgdTitle , fontSizes ,rectAndTextBox)

% specificy the y values starts and mode of actions for the drugs
% yPoint = 0.830 ; xStart = 0.1023 ;
xStartText = xStart + 0.01 ;
yPointTitle = yPoint + 0.03 ;

% specify the font size to be used in the plot
if ~exist('fontSizes','var')
    fontSizes = [10, 12] ;
end

% specifiy the rectangle and text length
if ~exist('rectAndTextBox','var')
    rectAndTextBox = [0.018 ,0.12] ;
end

% check for errors
if ~isnumeric(yPoint) || ~isnumeric(xStart)
    error('Both yPoint and xStarts should be numeric values')
elseif yPoint > 1 || xStart > 1
    error('Both yPoint and xStarts should be less than 1')
elseif ~isnumeric(plotColors)
    error('plot Color should be numeric')
end

if size(plotColors,1) ~= length(legendLabels)
    error('There should be a color for each legend names')
end

if iscategorical( legendLabels)
    legendLabels = categories(legendLabels);
end

for ii = 1:length(legendLabels)
    % add the legend color
    annotation('rectangle',[xStart yPoint rectAndTextBox(1) 0.023],...
        'EdgeColor', plotColors(ii,:), ...
        'FaceColor', plotColors(ii,:));
    
    % add the legend text
    annotation('textbox',[xStartText yPoint rectAndTextBox(2) 0.0230],...
        'String',legendLabels{ii},'FontSize',fontSizes(1),...
        'FontName','Helvetica Neue','FitBoxToText','off',...
        'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1] , ...
        'VerticalAlignment','middle','FontWeight','normal')
    
    % move the y point down
    yPoint = yPoint - 0.03 ;
end

% add the title
annotation('textbox',[xStart yPointTitle rectAndTextBox(2) 0.0230],...
    'String', myLgdTitle,'FontSize',fontSizes(2),...
    'FontName','Helvetica Neue','FitBoxToText','off',...
    'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1] , ...
    'VerticalAlignment','middle','FontWeight','bold',...
    'HorizontalAlignment','left');

end

% **************
%% The Internal Function 

function ManhattanPlotGrey( filename, varargin )

%   This function takes a GWAS output file  and plots a Manhattan Plot.

% ARGUMENTS:
%   sex: defaults to 0, set to 1 to include sex chromosomes

%   sig: defaults to 5e-8, significance threshold for horizontal line

%   vert: defaults to 0, set to 1 for tick labels have chr in them and are
%   written upwards

%   labels: defaults to [-1,-1]. Set to [x,y] to label the top SNP on each
%   locus with p<x. Locus defined by windows of y base pairs.

%   outfile: defaults to filename of input file, set to something else to
%   change name of output file

%   title: defaults to 'Manhattan Plot'. Fairly self-explanatory

%   save: defaults to 1, set to 0 to disable saving to save time

%   format: defaults to PLINK, options {PLINK,BOLT-LMM,SAIGE} This is used
%   to identify the correct column headings. If using anything else, rename
%   the header line so that CHR, BP, and P are the headers for chromosome,
%   base pair and p value, and the default PLINK should catch it

%   Usage:
%   ManhattanPlot('gwas.assoc.fisher',varargin) will take the association
%   analysis in 'gwas.assoc.fisher, generate a Manhattan Plot, and store it
%   in gwas_ManhattanPlot.png, which should be publication-ready, and
%   gwas_ManhattanPlot.fig for minor readjustments in MATLAB's GUI. This is
%   fine as .fisher is a PLINK format file. For a SAIGE output (and
%   significance threshold 5e-9) with sex chromosomes shown, and the top
%   hit for each SNP within 1 mb and below p=1e-6 labelled, use
%   ManhattanPlot('gwas.saige','format','SAIGE','sig',5e-9,'sex',1,
%   'labels',[1e-6,1000000])

%   Tested using assoc, assoc.logistic and assoc.fisher files generated by
%   Plink 1.7. Also tested using BOLT-LMM and SAIGE.
%
%   Harry Green (2019) Genetics of Complex Traits, University of Exeter

% Reading in optional arguments and input file

p = inputParser;
defaultsex = 0;
defaultvert = 0;
defaultsig = 5e-8;
defaultsave = 1;
defaultoutfile = filename;
defaulttitle= 'Manhattan Plot';
defaultlabels = [-1, -1];
defaultformat = 'PLINK';
expectedformat = {'PLINK','BOLT-LMM','SAIGE'};

addRequired(p,'filename');
% addOptional(p,'outfile',defaultoutfile,@ischar);
addOptional(p,'sex',defaultsex,@isnumeric);
addOptional(p,'vert',defaultvert,@isnumeric);
addOptional(p,'labels',defaultlabels,@isnumeric);
addOptional(p,'sig',defaultsig,@isnumeric);
addOptional(p,'save',defaultsave,@isnumeric);
addParameter(p,'format',defaultformat,...
    @(x) any(validatestring(x,expectedformat)));
addParameter(p,'outfile',defaultoutfile,@ischar);
addParameter(p,'title',defaulttitle,@ischar);
addParameter(p,'plotColor',@isnumeric);
% addParameter(p,'plotNumber',defaultplotNumber,@isnumeric);

parse(p,filename,varargin{:});

filename=p.Results.filename;
sex=p.Results.sex;
vert=p.Results.vert;
sig=p.Results.sig;
save=p.Results.save;
format=p.Results.format;
outfile=p.Results.outfile;
labels=p.Results.labels;
plottitle=p.Results.title;
plotColor=p.Results.plotColor ;
% plotNumber=p.Results.plotNumber ;

% ################ Check if the input is a table for a file name #########

if ischar(filename)
    % MATLAB refuses to read file formats like .assoc, so first rename as
    % .txt
    copyfile(filename, strcat(filename,'.txt'));
    opts = detectImportOptions(strcat(filename,'.txt'),'NumHeaderLines',0);
    
    % the columns of T will be the headers of the assoc file
    T = readtable(strcat(filename,'.txt'),opts);
    
    %delete unwanted .txt file that we didn't want anyway
    delete(strcat(filename,'.txt'))
elseif istable(filename)
    T = filename ;
end


% File format defintitions
% This section uses the default column headings from PLINK, BOLT-LMM and
% SAIGE. Easy to modify for other software packages

% check the CHR variable exist then change to PLINK format
if ~ismember(T.Properties.VariableNames,'CHR')
    T.Properties.VariableNames([1,2]) = {'CHR','BP'} ;
end

if strcmp(format,'PLINK')
    % These are the only useful columns for a manhattan plot, chromosome,
    % base pair, p value
    tab2=[T.CHR,T.BP,T.P]; 
end

if strcmp(format,'BOLT-LMM')
    % These are the only useful columns for a manhattan plot, chromosome,
    % base pair, p value
    tab2=[T.CHR,T.BP,T.P_BOLT_LMM]; 
end

if strcmp(format,'SAIGE')
    % These are the only useful columns for a manhattan plot, chromosome,
    % base pair, p value
    tab2=[T.CHR,T.POS,T.p_value]; 
end

%

if sex==0
    tab2=tab2(tab2(:,1)<23,:); %remove chr23 if not wanting sex chromosomes
end

% variable to track base pairs, this helps the plotting function to know
% where to start plotting the next chromosome
bptrack=0; 
tab2(tab2(:,3)==0,:)=[];

% get the signficant loci based on the pvalues
tabSig = tab2(tab2(:,3) < sig, :) ;

% get the  non significant loci
tab2 = tab2(tab2(:,3) >= sig, :) ;

% % get the plot
% figure(plotNumber)

% produce the plot in a loop
for i=1:22+sex
    hold on
    
    %a scatterplot. On the x axis, the base pair number + bptrack, which
    %starts at 0. On the y axis, - log 10 p
    plot( tab2(tab2(:,1)==i,2)+max(bptrack),...
        -log10(tab2(tab2(:,1)==i,3)), '.' , ...
        'MarkerEdgeColor',[0.5 0.5 0.5],...
        'MarkerFaceColor',[0.5 0.5 0.5],...
        'MarkerSize',10);
    
    % produce a scatter plot of the significant variants
    plot( tabSig(tabSig(:,1)==i,2)+max(bptrack),...
        -log10(tabSig(tabSig(:,1)==i,3)), '.' , ...
        'MarkerEdgeColor',plotColor,...
        'MarkerFaceColor',plotColor,...
        'MarkerSize',10);
    
    %this updates bptrack by adding the highest base pair number of the
    %chromosome. At the end, we should get the total number of base pairs
    %in the human genome. All values of bptrack are stored in a vector.
    %They're useful later
    bptrack=[bptrack,max(bptrack)+max(max(tab2(tab2(:,1)==i,:)))]; 
end

%if strongest hit is 1e-60, plot window goes up to 1e-61
if isempty(tabSig)
    %and down to the highest p value.
    ylimit=max(max(-log10(tab2(:,3)))+1);
    ylimitmin=min(min(-log10(tab2(:,3))));
else
    %and down to the highest p value.
    ylimit=max(max(-log10(tabSig(:,3)))+1);
    ylimitmin=min(min(-log10(tabSig(:,3))));
end

%genome wide significant line, uses the sig optional argument which
%defaults to 5e-8
plot([0,max(bptrack)],-log10([sig,sig]),'k--') 

% ylabel('$-\log_{10}(p)$','Interpreter','latex')
ylabel('-log_1_0(p)')
xlim([0,max(bptrack)])

% find the previous ylimit
curYlim = get(gca, 'YLim');

% compared with the new yLim
if max(ylimit,-log10(sig/5)) > curYlim
    %y axis will always go to the significance threshold/5 at least
    ylim([0,max(ylimit,-log10(sig/5))]) % floor(ylimitmin)
end

%this calculates the moving average of bptrack and uses them as chromosome
%label markers. This puts them in the middle of each chromosome.
M=movmean(bptrack,2); 
xticks(M(2:end));

% Rotation section
% this section of the code changes the x axis label orientation.

if vert ==0
    xlabel('Chromosome')% ,'Interpreter','latex')
    xticklabels( 1:23 );
end
if vert ==1
    xtickangle(90)
    xticklabels( {'chr1','chr2','chr3','chr4','chr5','chr6','chr7',...
        'chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15',...
        'chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrXY'});
end

% Annotation section
% this section of the code annotates the top SNP on each locus

% Neither of these should be negative. Defaults to -1 -1
if ~(labels(1)==-1||labels(2)==-1) 
    labellim=labels(1);
    labelprox=labels(2);
    
    %as nothing else will be used for labelling, the reduced table is
    %easier to manage
    tab3=tab2(tab2(:,3)<labellim,:); 
    
    %easier if they're in ascending order of p value;
    [~,index]=sort(tab3(:,3));
    tab3=tab3(index,:);
    
    %this becomes a list of SNPs that have been labelled. It starts with a
    %dummy SNP to give the first SNP something to compare with. This will
    %always fail so the first label always gets added
    labelledsnps=[0,0,0]; 
    
    for i=1:size(tab3,1)
        test=abs(tab3(i,:)-labelledsnps);
        %if there are no snps on the same chromosome and within 1 MB
        if sum(test(:,1)==0&test(:,2)<labelprox)==0 
            
            %this plots a square over the corresponding dot
            plot(bptrack(tab3(i,1))+tab3(i,2),-log10(tab3(i,3)),'ks',...
                'MarkerSize',8,'LineWidth',1) 
            
            %this puts together a strong of chrx:y
            labels=strcat('chr',string(tab3(i,1)),':',...
                string(bptrack(tab3(i,1))+tab3(i,2))); 
            
            %this plots the label above and to the left of the dot, shifted
            %up by 0.05 to give some space
            text(bptrack(tab3(i,1))+tab3(i,2),-log10(tab3(i,3))+0.05,...
                char(labels),'VerticalAlignment','bottom',...
                'HorizontalAlignment','left') %,'Interpreter','latex') 
            labelledsnps=[labelledsnps;tab3(i,:)];
        end
    end
end
% finalising file

% takes the output file name up to the first decimal point
a = strsplit(outfile,'.'); 

% my code 
set(gca,'FontSize',12,'LineWidth',0.5)

box on
title(plottitle,'fontsize',14)
% title(plottitle,'Interpreter','latex','fontsize',14)

if save~=0
    name=[a{1},'_ManhattanPlot.fig']; %adds _ManhattanPlot to filenames
    savefig(name);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 30 10])
    print([a{1},'_ManhattanPlot.png'],'-dpng','-r300', '-vector') % -r600
end



end


%% Another internal function 

function varargout = venn(varargin)

%VENN   Plot 2- or 3- circle area-proportional Venn diagram
%
%  venn(A, I)
%  venn(Z)
%  venn(..., F)
%  venn(..., 'ErrMinMode', MODE)
%  H = venn(...)
%  [H, S] = venn(...)
%  [H, S] = venn(..., 'Plot', 'off')
%  S = venn(..., 'Plot', 'off')
%  [...] = venn(..., P1, V1, P2, V2, ...) 
%
%venn(A, I) by itself plots circles with total areas A, and intersection
%area(s) I. For two-circle venn diagrams, A is a two element vector of circle 
%areas [c1 c2] and I is a scalar specifying the area of intersection between 
%them. For three-circle venn diagrams, A is a three element vector [c1 c2 c3], 
%and I is a four element vector [i12 i13 i23 i123], specifiying the 
%two-circle intersection areas i12, i13, i23, and the three-circle
%intersection i123.
%
%venn(Z) plots a Venn diagram with zone areas specified by the vector Z. 
%For a 2-circle venn diagram, Z is a three element vector [z1 z2 z12]
%For a 3-circle venn, Z is a 7 element vector [z1 z2 z3 z12 z13 z23 z123]
%
%venn(..., F) specifies optional optimization options. VENN uses FMINBND to
%locate optimum pair-wise circle distances, and FMINSEARCH to optimize
%overall three-circle alignment. F is a structure with fields specifying
%optimization options for these functions. F may be a two-element array of
%structures, in which case the first structure is used for FMINBND
%function calls, and the second structure is used for FMINSEARCH function
%calls.
%
%venn(..., 'ErrMinMode', MODE)
%Used for 3-circle venn diagrams only. MODE can be 'TotalError' (default), 
%'None', or 'ChowRodgers'. When ErrMinMode is 'None', the positions and 
%sizes of the three circles are fixed by their pairwise-intersections, 
%which means there may be a large amount of error in the area of the three-
%circle intersection. Specifying ErrMinMode as 'TotalError' attempts to 
%minimize the total error in all four intersection zones. The area of the 
%three circles are kept constant in proportion to their populations. The 
%'ChowRodgers' mode uses the the method proposed by Chow and Rodgers 
%[Ref. 1] to draw 'nice' three-circle venn diagrams which appear more 
%visually representative of the desired areas, although the actual areas of 
%the circles are allowed to deviate from requested values.
%
%H = venn(...) returns a two- or three- element vector to the patches 
%representing the circles. 
%
%[H, S] = venn(...) returns a structure containing descriptive values
%computed for the requested venn diagram. S is a structure with the
%following fields, where C is the number of circles (N = 2 or 3), Z is
%the number of zones (Z = 3 or 7), and I is the number of intersection 
%areas (1 or 4)
%
% Radius            C-element vector of circle radii
%
% Position          C*2 array of circle centers
%
% ZoneCentroid      Z*2 array of zone centroids (Can be used for labeling)
%
% CirclePop         C-element vector of supplied circle populations. 
%                   (I.e., the 'true' circle areas)
%
% CircleArea        C-element of actual circle areas
%
% CircleAreaError   = (CircleArea-CirclePop)/CirclePop
%
% IntersectPop      I-element vector of supplied intersection populations
%                   (I.e., the 'true' intersection areas)
%
% IntersectArea     I-element vector of actual intersection areas
%
% IntersectError    = (IntersectArea-IntersectPop)/IntersectPop
%
% ZonePop           Z-element vector of supplied zone populations. (I.e.
%                   'true' zone areas
%
% ZoneArea          Z-element vector of actual zone areas.
%
% ZoneAreaError     = (ZoneArea-ZonePop)/ZonePop
% 
%
%[H, S] = venn(..., 'Plot', 'off')
%S = venn(..., 'Plot', 'off')

% Returns a structure of computed values, without plotting the diagram.
% This which can be useful when S is used to draw custom venn diagrams or
% for exporting venn diagram data to another application. When Plot is set
% to off, the handles vector H is returned as an empty array.
% Alternatively, the command

% S = venn(..., 'Plot', 'off) will return only the output structure.

%[...] = venn(..., P1, V1, P2, V2, ...) 

%Specifies additional patch settings in standard Matlab parameter/value
%pair syntax. Parameters can be any valid patch parameter. Values for patch
%parameters can either be single values, or a cell array of length
%LENGTH(A), in which case each value in the cell array is applied to the
%corresponding circle in A.

%Examples
%
%   %Plot a simple 2-circle venn diagram with custom patch properties
%   figure, axis equal, axis off
%   A = [300 200]; I = 150;
%   venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black')
%
%   %Compare ErrMinModes
%   A = [350 300 275]; I = [100 80 60 40];
%   figure
%   subplot(1,3,1), h1 = venn(A,I,'ErrMinMode','None');
%   axis image,  title ('No 3-Circle Error Minimization')
%   subplot(1,3,2), h2 = venn(A,I,'ErrMinMode','TotalError');
%   axis image,  title ('Total Error Mode')
%   subplot(1,3,3), h3 = venn(A,I,'ErrMinMode','ChowRodgers');
%   axis image, title ('Chow-Rodgers Mode')
%   set([h1 h2], 'FaceAlpha', 0.6)
%
%   %Using the same areas as above, display the error optimization at each 
%   iteration. Get the output structure.
%   F = struct('Display', 'iter');
%   [H,S] = venn(A,I,F,'ErrMinMode','ChowRodgers','FaceAlpha', 0.6);
%
%   %Now label each zone 
%   for i = 1:7
%       text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), ['Zone ' num2str(i)])
%   end
%
%See also patch, bar, optimset, fminbdn, fminsearch
%
%Copyright (C) 2008 Darik Gamble, University of Waterloo.
%dgamble@engmail.uwaterloo.ca
%
%References
%1. S Chow and P Rodgers. Extended Abstract: Constructing Area-Proportional
%   Venn and Euler Diagrams with Three Circles. Presented at Euler Diagrams 
%   Workshop 2005. Paris. Available online: 
%   http://www.cs.kent.ac.uk/pubs/2005/2354/content.pdf
%
%2. S Chow and F Ruskey. Drawing Area-Proportional Venn and Euler Diagrams. 
%   Lecture Notes in Computer Science. 2004. 2912: 466-477. Springer-Verlag. 
%   Available online: http://www.springerlink.com/content/rxhtlmqav45gc84q/
%
%3. MP Fewell. Area of Common Overlap of Three Circles. Australian Government 
%   Department of Defence. Defence Technology and Science Organisation. 2006. 
%   DSTO-TN-0722. Available online:
%   http://dspace.dsto.defence.gov.au/dspace/bitstream/1947/4551/4/DSTO-TN-0722.PR.pdf


%Variable overview
%   A0, A   Desired and actual circle areas
%               A = [A1 A2] or [A1 A2 A3]
%   I0, I   Desired and actual intersection areas
%               I = I12 or [I12 I13 I23 I123]
%   Z0, Z   Desired and actual zone areas
%               Z = [Z1 Z2 Z12] or [Z1 Z2 Z3 Z12 Z13 Z23 Z123]
%   x, y    Circle centers
%               x = [x1 x2] or [x1 x2 x3]
%   r       Circle radii
%               r = [r1 r2] or [r1 r2 r3]
%   d       Pair-wise distances between circles
%               d = d12 or [d12 d13 d23]
    


    %Parse input arguments and preallocate settings
    [A0, I0, Z0, nCirc, fminOpts, vennOpts, patchOpts] = ...
        parseArgsIn(varargin);
    [d, x, y, A, I, Z] = preallocVectors (nCirc);
    zoneCentroids = []; %Will only be calculated if needed
    
    %Circle Radii
    r = sqrt(A0/pi);

    %Determine distance between first circle pair    
    d(1) = circPairDist(r(1), r(2), I0(1), fminOpts(1));
    
    %Position of second circle is now known
    x(2) = d(1); 
    
    %First intersection area 
    I(1) = areaIntersect2Circ(r(1), r(2), d(1));
    
    if nCirc==3
        %Pairwise distances for remaining pairs 1&3 and 2&3
        d(2) = circPairDist(r(1), r(3), I0(2), fminOpts(1)); %d13
        d(3) = circPairDist(r(2), r(3), I0(3), fminOpts(1)); %d23

        %Check triangle inequality
        srtD = sort(d);
        if ~(srtD(end)<(srtD(1)+srtD(2)))
            error('venn:triangleInequality', ...
                'Triangle inequality not satisfied')
        end

        %Guess the initial position of the third circle using the law of
        %cosines
        alpha = acos( (d(1)^2 + d(2)^2 - d(3)^2)  / (2 * d(1) * d(2)) );
        x(3) = d(2)*cos(alpha);
        y(3) = d(2)*sin(alpha);

        %Each pair-wise intersection fixes the distance between each pair
        %of circles, so technically there are no degrees of freedom left in
        %which to adjust the three-circle intersection. We can either try
        %moving the third circle around to minimize the total error, or
        %apply Chow-Rodgers 
        
        switch vennOpts.ErrMinMode
            case 'TotalError'
                %Minimize total intersection area error by moving the third
                %circle
                pos = fminsearch(@threeCircleAreaError, [x(3) y(3)],...
                    fminOpts(2));
                x(3) = pos(1);
                y(3) = pos(2);
            case 'ChowRodgers'
                %note that doChowRodgersSearch updates x and y in this
                %workspace as a nested fcn
                doChowRodgersSearch;
        end

        %Make sure everything is 'up to date' after optimization
        update3CircleData;
        
    end
    
    % Are we supposed to plot?
    if vennOpts.Plot
        if isempty(vennOpts.Parent)
            vennOpts.Parent = gca;
        end
        hVenn = drawCircles(vennOpts.Parent, x, y, r, ...
            patchOpts.Parameters, patchOpts.Values);
    else
        hVenn = [];
    end
    
    %Only determine zone centroids if they're needed 
    %Needed for output structure 
    nOut = nargout;
    if (nOut==1 && ~vennOpts.Plot) || nOut==2
        if nCirc == 2
            %Need to calculate new areas
            A = A0; %Areas never change for 2-circle venn
            Z = calcZoneAreas(2, A, I);
            zoneCentroids = zoneCentroids2(d, r, Z);
        else
            zoneCentroids = zoneCentroids3(x, y, d, r, Z);
        end
    end
        
    %Figure out output arguments
    if nOut==1
        if vennOpts.Plot
            varargout{1} = hVenn;
        else
            varargout{1} = getOutputStruct;
        end
    elseif nOut==2
        varargout{1} = hVenn;
        varargout{2} = getOutputStruct;
    end
        
    
    
    function err = threeCircleAreaError (pos)
        
        x3 = pos(1);
        y3 = pos(2);
        
        %Calculate distances
        d(2) = sqrt(x3^2 + y3^2); %d13
        d(3) = sqrt((x3-d(1))^2 + y3^2); %d23
        
        %Calculate intersections
        %Note: we're only moving the third circle, so I12 is not changing
        I(2:3) = areaIntersect2Circ (r(1:2), r([3 3]), d(2:3)); %I13 and I23
        I(4) = areaIntersect3Circ (r, d); %I123
        
        %Replace 0 (no intersection) with infinite error
        I(I==0) = Inf;
        
        %Error
        err = sum(abs((I-I0)./I0));
        
    end


    function doChowRodgersSearch
        
        %Adapted from Ref. [1]
        
        %Initialize an index matrix to select all 7choose2 zone pairs (21 pairs)
        idx = nchoosek(1:7, 2);
                
        %Which zone-zone pairs are considered equal?
        %Zones within 10% of each other considered equal
        zonePairAreas0 = Z0(idx);
        
        %Percent difference in population between the two members of a pair
        ar0 = 2*abs(zonePairAreas0(:,1)-zonePairAreas0(:,2))./sum(zonePairAreas0, 2)*100;
        eqPairCutoff = 10;  
        pairIsEq = ar0<=eqPairCutoff;
        
        %Calculate allowable range for pairs of zones considered unequal
        if any(~pairIsEq)
            %Sort zone areas
            [zUneqAreas0, zUneqAreasSrtIdx] = sort(zonePairAreas0(~pairIsEq,:), 2);
            
            %Make a real index array out of the inconvenient index sort returns
            n = sum(~pairIsEq);
            zUneqAreasSrtIdx = sub2ind([n,2], [1:n; 1:n]', zUneqAreasSrtIdx);
            
            %rp = (largepopulation/smallpopulation)-1
            rp = zUneqAreas0(:,2)./zUneqAreas0(:,1)-1;
            rpMin = 1 + 0.3*rp;
            rpMax = 1 + 2*rp;
        end
        
        %Preallocate zone error vector
        zoneErr = zeros(1,21); 

        %Initialize independent parameters to search over
        guessParams = [r(1) x(2) r(2) x(3) y(3) r(3)];
        
        %Search!
        pp = fminsearch(@chowRodgersErr, guessParams, fminOpts(2));  
        
        [r(1) x(2) r(2) x(3) y(3) r(3)] = deal(pp(1), pp(2), pp(3), pp(4), pp(5), pp(6));
        
        
        function err = chowRodgersErr (p)
            
            %params = [x2 r2 x3 y3 r3]
            [r(1), x(2), r(2), x(3), y(3), r(3)] = deal(p(1), p(2), p(3), p(4), p(5), p(6));
             
            %After changing x2, r2, x3, y3, and r3, update circle areas,
            %distances, intersection areas, zone areas
            update3CircleData;

            if any(pairIsEq)
                %For zone pairs considered equal, error is equal to square of the
                %distance beyond the cutoff; 0 within cutoff
                zAreas = Z(idx(pairIsEq,:));
                ar = 2*abs(zAreas(:,1)-zAreas(:,2))./sum(zAreas, 2)*100;
                isWithinRange = ar<eqPairCutoff;
                ar(isWithinRange) = 0;
                ar(~isWithinRange) = ar(~isWithinRange) - eqPairCutoff;

                %Amplify error for equal zones with unequal areas
                eqZoneUneqAreaErrorGain = 10;
                ar(~isWithinRange) = ar(~isWithinRange)*eqZoneUneqAreaErrorGain;

                zoneErr(pairIsEq) = ar.^2;
            end

            if any(~pairIsEq)
                %For zone pairs considered unequal, error is equal to square of
                %the distance from the allowable range of rp

                %rp = (largepopulation/smallpopulation)-1
                zUneqPairAreas = Z(idx(~pairIsEq,:));
                
                %Sort based on the population sizes (determined by parent
                %function doChowRodgersSearch)
                zUneqPairAreas = zUneqPairAreas(zUneqAreasSrtIdx);
                rp = zUneqPairAreas(:,2)./zUneqPairAreas(:,1)-1;

                lessThanMin = rp<rpMin;
                moreThanMax = rp>rpMax;
                rp(~lessThanMin & ~moreThanMax) = 0;
                
                %Determine how far out of range errors are
                rp(lessThanMin) = rp(lessThanMin) - rpMin(lessThanMin);
                rp(moreThanMax) = rp(moreThanMax) - rpMax(moreThanMax);          

                %Consider the case where rp < rpMin to be more
                %erroneous than the case where rp > rpMax 
                tooSmallErrorGain = 10;
                rp(lessThanMin) = rp(lessThanMin)*tooSmallErrorGain;

                zoneErr(~pairIsEq) = rp.^2;
            end
            
            %Total error
            err = sum(zoneErr);
            
        end %chowRodgersErr
        
    end %doChowRodgersSearch

    function update3CircleData
        
        %Circle areas
        A = pi*r.^2;

        %Calculate distances
        d(1) = abs(x(2)); %d12
        d(2) = sqrt(x(3)^2 + y(3)^2); %d13
        d(3) = sqrt((x(3)-d(1))^2 + y(3)^2); %d23

        %Calculate actual intersection areas
        I(1:3) = areaIntersect2Circ (r([1 1 2]), r([2 3 3]), d); %I12, I13, I23
        I(4) = areaIntersect3Circ (r, d); %I123

        %Calculate actual zone areas
        Z = calcZoneAreas(3, A, I);

    end

    function S = getOutputStruct
               
        S = struct(...
            'Radius'                ,r                      ,...        
            'Position'              ,[x' y']                ,...
            'ZoneCentroid'          ,zoneCentroids          ,...
            'CirclePop'             ,A0                     ,...
            'CircleArea'            ,A                      ,...
            'CircleAreaError'       ,(A-A0)./A0             ,...
            'IntersectPop'          ,I0                     ,...
            'IntersectArea'         ,I                      ,...
            'IntersectError'        ,(I-I0)./I0             ,...
            'ZonePop'               ,Z0                     ,...
            'ZoneArea'              ,Z                      ,...
            'ZoneAreaError'         ,(Z-Z0)./Z0             );  
        end

end %venn

        
function D = circPairDist (rA, rB, I, opts)
    %Returns an estimate of the distance between two circles with radii rA and
    %rB with area of intersection I
    %opts is a structure of FMINBND search options
    D = fminbnd(@areadiff, 0, rA+rB, opts);
    function dA = areadiff (d)
        intersectArea = areaIntersect2Circ (rA, rB, d);
        dA = abs(I-intersectArea)/I;
    end
end

function hCirc = drawCircles(hParent, xc, yc, r, P, V,c)

    hAx = ancestor(hParent, 'axes');
    nextplot = get(hAx, 'NextPlot');
    
    %P and V are cell arrays of patch parameter/values
    xc = xc(:); yc = yc(:);     %Circle centers
    r = r(:);                   %Radii
    n = length(r);              
    
    %Independent parameter
    dt = 0.05;
    t = 0:dt:2*pi;

    % Origin centered circle coordinates
    X = r*cos(t);
    Y = r*sin(t);
    
    hCirc = zeros(1,n);
    % c = [0.00,0.45,0.74; 0.85,0.33,0.10; 0.47,0.67,0.19];
    % c = [0,0.450,0.740;0.640,0.080,0.18] ; %blue and orange-red
    c = [ 0.85,0.33,0.10; 0.47,0.67,0.19  ]; % green and orange-brown
    
    % c = {'r', 'g', 'b'};                        %default colors
    fa = {0.6, 0.6, 0.6};                       %default face alpha
    tag = {'Circle1', 'Circle2', 'Circle3'}; 	%default tag
    
    for i = 1:n
        xx = X(i,:)+xc(i);  
        yy = Y(i,:)+yc(i);
        
        if iscell(c)
            % if the color is cell array
            hCirc(i) = patch (xx, yy, c{i}, 'FaceAlpha', fa{i},...
                'Parent', hParent, 'Tag', tag{i});
        else % or my colors
            hCirc(i) = patch (xx, yy, c(i,:), 'FaceAlpha', fa{i},...
                'Parent', hParent, 'Tag', tag{i} ,'LineWidth',1);
        end
        
        if i==1
            set(hAx, 'NextPlot', 'add');
        end
    end
    set(hAx, 'NextPlot', nextplot);
    
    % make the axis invisible 
    set(findobj(gcf, 'type','axes'), 'Visible','off')

%     % Custom patch parameter values
%     if ~isempty(P)
% 
%         c = cellfun(@iscell, V);
% 
%         %Scalar parameter values -- apply to all circles
%         if any(~c)
%             set(hCirc, {P{~c}}, {V{~c}});
%         end
% 
%         %Parameters values with one value per circle
%         if any(c)
%             %Make sure all vals are column cell arrays
%             V = cellfun(@(val) (val(:)), V(c), 'UniformOutput', false);
%             set(hCirc, {P{c}}, [V{:}])
%         end
%     end
    
end %plotCircles

 
function A = areaIntersect2Circ (r1, r2, d)
    %Area of Intersection of 2 Circles
    %Taken from [2]
    
    alpha = 2*acos( (d.^2 + r1.^2 - r2.^2)./(2*r1.*d) );
    beta  = 2*acos( (d.^2 + r2.^2 - r1.^2)./(2*r2.*d) );
    
    A =    0.5 * r1.^2 .* (alpha - sin(alpha)) ...  
         + 0.5 * r2.^2 .* (beta - sin(beta));
    
end

function [A, x, y, c, trngArea] = areaIntersect3Circ (r, d)
    %Area of common intersection of three circles
    %This algorithm is taken from [3]. 
    %   Symbol    Meaning
    %     T         theta
    %     p         prime
    %     pp        double prime
        
    %[r1 r2 r3] = deal(r(1), r(2), r(3));
    %[d12 d13 d23] = deal(d(1), d(2), d(3));

    %Intersection points
    [x,y,sinTp,cosTp] = intersect3C (r,d);
    
    if any(isnan(x)) || any(isnan(y))
        A = 0;
        %No three circle intersection
        return
    end
    
    %Step 6. Use the coordinates of the intersection points to calculate the chord lengths c1,
    %c2, c3:
    i1 = [1 1 2];
    i2 = [2 3 3];
    c = sqrt((x(i1)-x(i2)).^2 + (y(i1)-y(i2)).^2)';

    %Step 7: Check whether more than half of circle 3 is included in the circular triangle, so
    %as to choose the correct expression for the area
    lhs = d(2) * sinTp;
    rhs = y(2) + (y(3) - y(2))/(x(3) - x(2))*(d(2)*cosTp - x(2));
    if lhs < rhs
        sign = [-1 -1 1];
    else
        sign = [-1 -1 -1];
    end
    
    %Calculate the area of the three circular segments.
    ca = r.^2.*asin(c/2./r) + sign.*c/4.*sqrt(4*r.^2 - c.^2);

    trngArea = 1/4 * sqrt( (c(1)+c(2)+c(3))*(c(2)+c(3)-c(1))*(c(1)+c(3)-c(2))*(c(1)+c(2)-c(3)) );
    A = trngArea + sum(ca);
    
end

function [x, y, sinTp, cosTp] = intersect3C (r, d)
    %Calculate the points of intersection of three circles
    %Adapted from Ref. [3]
    
    %d = [d12 d13 d23]
    %x = [x12; x13; x23]
    %y = [y12; y13; y23]

    %   Symbol    Meaning
    %     T         theta
    %     p         prime
    %     pp        double prime
    
    x = zeros(3,1);
    y = zeros(3,1);
     
    %Step 1. Check whether circles 1 and 2 intersect by testing d(1)
    if ~( ((r(1)-r(2))<d(1)) && (d(1)<(r(1)+r(2))) )
        %x = NaN; y = NaN;
        %bigfix: no returned values for sinTp, cosTp
        [x, y, sinTp, cosTp] = deal(NaN);
        return
    end

    %Step 2. Calculate the coordinates of the relevant intersection point of circles 1 and 2:
    x(1) = (r(1)^2 - r(2)^2 + d(1)^2)/(2*d(1));
    y(1) = 0.5/d(1) * sqrt( 2*d(1)^2*(r(1)^2 + r(2)^2) - (r(1)^2 - r(2)^2)^2 - d(1)^4 );

    %Step 3. Calculate the values of the sines and cosines of the angles tp and tpp:
    cosTp  =  (d(1)^2 + d(2)^2 - d(3)^2) / (2 * d(1) * d(2));
    cosTpp = -(d(1)^2 + d(3)^2 - d(2)^2) / (2 * d(1) * d(3));
    sinTp  =  (sqrt(1 - cosTp^2));
    sinTpp =  (sqrt(1 - cosTpp^2));

    %Step 4. Check that circle 3 is placed so as to form a circular triangle.
    cond1 = (x(1) - d(2)*cosTp)^2 + (y(1) - d(2)*sinTp)^2 < r(3)^2;
    cond2 = (x(1) - d(2)*cosTp)^2 + (y(1) + d(2)*sinTp)^2 > r(3)^2;
    if  ~(cond1 && cond2)
        x = NaN; y = NaN;
        return
    end

    %Step 5: Calculate the values of the coordinates of the relevant intersection points involving
    %circle 3
    xp13  =  (r(1)^2 - r(3)^2 + d(2)^2) / (2 * d(2));
    %yp13  = -0.5 / d(2) * sqrt( 2 * d(2)^2 * (r(2)^2 + r(3)^2) - (r(1)^2 - r(3)^2)^2 - d(2)^4 );
    yp13  = -0.5 / d(2) * sqrt( 2 * d(2)^2 * (r(1)^2 + r(3)^2) - (r(1)^2 - r(3)^2)^2 - d(2)^4 );

    x(2)   =  xp13*cosTp - yp13*sinTp;
    y(2)   =  xp13*sinTp + yp13*cosTp;

    xpp23 =  (r(2)^2 - r(3)^2 + d(3)^2) / (2 * d(3));
    ypp23 =  0.5 / d(3) * sqrt( 2 * d(3)^2 * (r(2)^2 + r(3)^2) - (r(2)^2 - r(3)^2)^2 - d(3)^4 );

    x(3) = xpp23*cosTpp - ypp23*sinTpp + d(1);
    y(3) = xpp23*sinTpp + ypp23*cosTpp;

end



function z = calcZoneAreas(nCircles, a, i)
    
    %Uses simple set addition and subtraction to calculate the zone areas
    %with circle areas a and intersection areas i

    if nCircles==2
        %a = [A1 A2]
        %i = I12
        %z = [A1-I12, A2-I12, I12]
        z = [a(1)-i, a(2)-i, i];
    elseif nCircles==3
        %a = [A1  A2  A3]
        %i = [I12 I13 I23 I123]
        %z = [A1-I12-I13+I123, A2-I12-I23+I123, A3-I13-I23+I123, ...
        %     I12-I123, I13-I123, I23-I123, I123];
        z = [a(1)-i(1)-i(2)+i(4), a(2)-i(1)-i(3)+i(4), a(3)-i(2)-i(3)+i(4), ...
                i(1)-i(4), i(2)-i(4), i(3)-i(4), i(4)];
    else
        error('')
        %This error gets caught earlier in the stack w. better error msgs
    end
end

function [Cx, Cy, aiz] = centroid2CI (x, y, r)

    %Finds the centroid of the area of intersection of two circles.
    %Vectorized to find centroids for multiple circle pairs
    %x, y, and r are nCirclePairs*2 arrays
    %Cx and Cy are nCirclePairs*1 vectors

    %Centroid of the area of intersection of two circles
    n = size(x,1);
    xic = zeros(n,2);
    az = zeros(n,2);
    
    dx = x(:,2)-x(:,1);
    dy = y(:,2)-y(:,1);
    d = sqrt(dx.^2 + dy.^2);
    
    %Translate the circles so the first is at (0,0) and the second is at (0,d)
    %By symmetry, all centroids are located on the x-axis.
    %The two circles intersect at (xp, yp) and (xp, -yp)
    xp = 0.5*(r(:,1).^2 - r(:,2).^2 + d.^2)./d;

    %Split the inner zone in two
    %Right side (Area enclosed by circle 1 and the line (xp,yp) (xp,-yp)
    %Angle (xp,yp) (X1,Y1) (xp,-yp)
    alpha = 2*acos(xp./r(:,1));
    %Area and centroid of the right side of the inner zone
    [xic(:,1) az(:,1)] = circleChordVals (r(:,1), alpha);
    %Angle (xp,yp) (X2,Y2) (xp,-yp)
    alpha = 2*acos((d-xp)./r(:,2));
    %Area and centroid of the left side of the inner zone
    [xic(:,2) az(:,2)] = circleChordVals (r(:,2), alpha);
    xic(:,2) = d - xic(:,2);

    %Thus the overall centroid  & area of the inner zone
    aiz = sum(az,2);
    Cx = sum(az.*xic,2)./aiz;
    
    %Now translate the centroid back based on the original positions of the
    %circles
    theta = atan2(dy, dx);
    Cy = Cx.*sin(theta) + y(:,1);
    Cx = Cx.*cos(theta) + x(:,1);
    
end

function centroidPos = zoneCentroids2 (d, r, Z)
    
    centroidPos = zeros(3,2);
    
    %Find the centroids of the three zones in a 2-circle venn diagram
    %By symmetry, all centroids are located on the x-axis.
    %First, find the x-location of the middle (intersection) zone centroid
    
    %Centroid of the inner zone
    centroidPos(3,1) = centroid2CI([0 d], [0 0], r);
    
    %Now, the centroid of the left-most zone is equal to the centroid of
    %the first circle (0,0) minus the centroid of the inner zone
    centroidPos(1,1) = -centroidPos(3,1)*Z(3)/Z(1);
    
    %Similarly for the right-most zone; the second circle has centroid at x=d
    centroidPos(2,1) = (d*(Z(2)+Z(3)) - centroidPos(3,1)*Z(3))/Z(2);
    
end

function centroidPos = zoneCentroids3 (x0, y0, d, r, Z)

    Z = Z(:);
        
    %Get area, points of intersection, and chord lengths
    [act, xi, yi, c, atr] = areaIntersect3Circ (r, d);
    atr = atr(:);
    r = r(:);
    
    %Area and centroid of the triangle within the circular triangle is
    xtr = sum(xi/3); 
    ytr = sum(yi/3);
        
    %Now find the centroids of the three segments surrounding the triangle
    i = [1 2; 1 3; 2 3]; 
    xi = xi(i); yi = yi(i);
    [xcs, ycs, acs] = circSegProps (r(:), x0(:), y0(:), xi, yi, c(:));
    
    %Overall centroid of the circular triangle
    xct = (xtr*atr + sum(xcs.*acs))/act;
    yct = (ytr*atr + sum(ycs.*acs))/act;
    
    %Now calculate the centroids of the three two-pair intersection zones
    %(Zones 12 13 23)
    %Entire zone centroid/areas

    %x, y, and r are nCirclePairs*2 arrays
    %Cx and Cy are nCirclePairs*1 vectors
    i = [1 2; 1 3; 2 3];
    [x2c, y2c, a2c] = centroid2CI (x0(i), y0(i), r(i));
    
    %Minus the three-circle intersection zone
    xZI2C = (x2c.*a2c - xct*act)./(a2c-act);
    yZI2C = (y2c.*a2c - yct*act)./(a2c-act);
    
    x0 = x0(:);
    y0 = y0(:);
    
    %Finally, the centroids of the three circles minus the intersection
    %areas
    i1 = [4 4 5]; i2 = [5 6 6];
    j1 = [1 1 2]; j2 = [2 3 3];
    x1C = (x0*pi.*r.^2 - xZI2C(j1).*Z(i1) - xZI2C(j2).*Z(i2) - xct*act)./Z(1:3);
    y1C = (y0*pi.*r.^2 - yZI2C(j1).*Z(i1) - yZI2C(j2).*Z(i2) - yct*act)./Z(1:3);
    
    %Combine and return
    centroidPos = [x1C y1C; xZI2C yZI2C; xct yct];
end


function [x, a] = circleChordVals (r, alpha)
    %For a circle centered at (0,0), with angle alpha from the x-axis to the 
    %intersection of the circle to a vertical chord, find the x-centroid and
    %area of the region enclosed between the chord and the edge of the circle
    %adapted from http://mathworld.wolfram.com/CircularSegment.html
    a = r.^2/2.*(alpha-sin(alpha));                         %Area
    x = 4.*r/3 .* sin(alpha/2).^3 ./ (alpha-sin(alpha));    %Centroid
end

function [xc, yc, area] = circSegProps (r, x0, y0, x, y, c)

    %Translate circle to (0,0)
    x = x-[x0 x0];
    y = y-[y0 y0];

    %Angle subtended by chord
    alpha = 2*asin(0.5*c./r);
       
    %adapted from http://mathworld.wolfram.com/CircularSegment.html
    area = r.^2/2.*(alpha-sin(alpha));                         %Area
    d   = 4.*r/3 .* sin(alpha/2).^3 ./ (alpha-sin(alpha));    %Centroid
   
    %Perpindicular bisector of the chord
    m = -(x(:,2)-x(:,1))./(y(:,2)-y(:,1));
    
    %angle of bisector
    theta = atan(m);
    
    %centroids
    xc = d.*cos(theta);
    yc = d.*sin(theta);
    
    %Make sure we're on the correct side
    %Point of intersection of the perp. bisector and the circle perimeter
    xb = (x(:,1)+x(:,2))/2;
    xc(xb<0) = xc(xb<0)*-1;
    yc(xb<0) = yc(xb<0)*-1;
    
    %Translate back
    xc = xc + x0;
    yc = yc + y0;
end


function [A0, I0, Z0, nCircles, fminOpts, vennOpts, patchOpts] = parseArgsIn (args)

    [A0, I0, Z0] = deal([]);
    nIn = length(args);
    badArgs = false;
    
    %Get the easy cases out of the way
    if nIn == 0
        badArgs = true;
    elseif nIn == 1
        %venn(Z)
        Z0 = args{1};
        nIn = 0;
    elseif nIn == 2
        if isnumeric(args{2})
            %venn (A,I)
            [A0, I0] = deal(args{1:2});
            nIn = 0;
        else
            %venn (Z, F)
            Z0 = args{1};
            args = args(2);
            nIn = 1;
        end
    else
        %Find the first non-numeric input arg
        i = find(~cellfun(@isnumeric, args), 1);
        if i == 2
            %venn(Z, ....)
            Z0 = args{1};
        elseif i == 3
            %venn(A, I, ...)
            [A0, I0] = deal(args{1:2});
        else
            badArgs = true;
        end
        nIn = nIn - i + 1;
        args = args(i:end);
    end
    
    if badArgs
        error('venn:parseInputArgs:unrecognizedSyntax', 'Unrecogized input syntax')
    end
    try
        [A0, I0, Z0] = parseInputAreas (A0, I0, Z0);
    catch
        error('venn:parseArgsIn:parseInputAreas', 'Incorrect size(s) for area vector(s)')
    end
    nCircles = length(A0);
    nZones = length(Z0);
             
    %Any arguments left?
    if nIn > 0 
        
        if isstruct(args{1})
            %FMIN search options
            f = args{1};
           
            nIn = nIn - 1;
            if nIn>0, args = args(2:end); end

            if length(f) == 1
                %Just double up
                fminOpts = [f f];
            elseif length(f) == 2
                %ok
                fminOpts = f;
            else
                error('venn:parseArgsIn', 'FMINOPTS must be a 1 or 2 element structure array.')
            end
        else
            %Use defaults
            fminOpts = [optimset('fminbnd'), optimset('fminsearch')];
        end
    else
        %Use defaults
        fminOpts = [optimset('fminbnd'), optimset('fminsearch')];
    end

    %If there's an even number of args in remaining
    if nIn>0 
        if mod(nIn, 2)==0
            %Parameter/Value pairs
            p = args(1:2:end);
            v = args(2:2:end);
            [vennOpts, patchOpts] = parsePVPairs (p, v, nZones);
        else
            error('venn:parseArgsIn', 'Parameter/Value options must come in pairs')
        end
    else
        vennOpts = defaultVennOptions;
        patchOpts = struct('Parameters', [], 'Values', []);
    end

end %parseArgsIn

function [vennOpts, patchOpts] = parsePVPairs (p, v, nZones)

    p = lower(p);

    %Break up P/V list into Venn parameters and patch parameters
    vennParamNames = {'plot', 'errminmode', 'parent'};
    [isVennParam, idx] = ismember(p, vennParamNames);
    idx = idx(isVennParam);
    %vennParams = p(isVennParam);
    vennVals = v(isVennParam);
    
    %First do Patch options
    patchOpts.Parameters = p(~isVennParam);
    patchOpts.Values = v(~isVennParam);
        
    %Now do Venn options
    vennOpts = defaultVennOptions;
        
    %PLOT
    i = find(idx==1, 1);
    if i
        plot = lower(vennVals{i});
        if islogical(plot)
            vennOpts.Plot = plot;
        else
            if ischar(plot) && any(strcmp(plot, {'on', 'off'}))
                vennOpts.Plot = strcmp(plot, 'on');
            else
                error('venn:parsePVPairs', 'Plot must be ''on'', ''off'', or a logical value.')
            end
        end
    end
    
    %ERRMINMODE
    i = find(idx==2, 1);
    if i
        mode = lower(vennVals{i});
        okModes = {'None', 'TotalError', 'ChowRodgers'};
        [isOkMode, modeIdx] = ismember(mode, lower(okModes));
        if isOkMode                
            vennOpts.ErrMinMode = okModes{modeIdx};
        else
            error('venn:parsePVPairs', 'ErrMinMode must be None, TotalError, or ChowRodgers')
        end
    end

    %PARENT
    i = find(idx==5, 1);
    if i
        h = v{i};
        if length(h)==1 && ishandle(h) 
            vennOpts.Parent = h;
        else
            error('venn:parsePVPairs', 'Parent must be a valid scalar handle')
        end
    end
    
       
end %parsePVPairs

function [A0, I0, Z0] = parseInputAreas (A0, I0, Z0)

    %Switch to row vectors
    A0 = A0(:)';
    I0 = I0(:)';
    Z0 = Z0(:)';

    if isempty(Z0)
        %A0 and I0 supplied
        
        Z0 = calcZoneAreas (length(A0), A0, I0);
    else
        %Z0 supplied
        switch length(Z0)
            case 3
                A0 = Z0(1:2)+Z0(3);
                I0 = Z0(3);
            case 7
                A0 = Z0(1:3)+Z0([4 4 5])+Z0([5 6 6])+Z0(7);
                I0 = [Z0(4:6)+Z0(7) Z0(7)];
            otherwise
                error('')
        end
    end
end

function vennOpts = defaultVennOptions 
    
    vennOpts = struct(...
        'Plot'          ,true               ,...
        'Labels'        ,[]                 ,...
        'PopLabels'     ,false              ,...
        'DrawLabels'    ,false              ,...
        'Parent'        ,[]                 ,...
        'Offset'        ,[0 0]              ,...
        'ErrMinMode'    ,'TotalError'       );
    
end

function [d, x, y, A, I, Z] = preallocVectors (nCirc)

    %Initialize position vectors
    x = zeros(1, nCirc);
    y = zeros(1, nCirc);

    if nCirc==2
        d = 0;
        I = 0;    
        A = zeros(1,2);
        Z = zeros(1,3);

    else %nCirc==3
        d = zeros(1,3);
        I = zeros(1,4);
        A = zeros(1,3);
        Z = zeros(1,7);
    end
end
    
%% Another internal function 

% Helper: run trend tests for a trait
% ---------------------------------
function rowOut = runTrendTests(tbl, xName, yName, raceLabel, confLabel)
    % xName : percentile variable (numeric)
    % yName : outcome (numeric)
    % raceLabel : 'AFR' or 'EUR'
    % confLabel : 'BMI', 'Age', etc. (string)
    %
    % Returns a table row with statistics
    
    x  = double(tbl.(xName));  % numeric vector of percentiles
    y  = tbl.(yName);          % outcome
    
    % Linear regression slope and p-value
    lmFit   = fitlm(x, y);
    slope   = lmFit.Coefficients{'x1','Estimate'};
    pSlope  = lmFit.Coefficients{'x1','pValue'};
    
    % Spearman
    [rho, pSpearman] = corr(x, y, 'Type', 'Spearman');
    
    % Construct a one-row table
    rowOut = table({raceLabel}, {confLabel}, {yName}, slope, pSlope, rho, pSpearman, ...
                   'VariableNames', {'Ancestry','Confounder','Trait', ...
                                     'OLS_Slope','OLS_p','SpearmanRho','Spearman_p'});
end

%% Another internal function 

function colourBoxPlot(plotData,groups, color, includeScatter)

% set the color to the box plots
if nargin == 2 || isempty(color)
    rng(6);
    color = rand(length(unique(groups )),3) ;
end

% plot the data
figure()
boxplot(plotData,groups,'Color', flipud(color) ,'Symbol','k+', ...
    'OutlierSize',5) ;

% set some figure properties and add title ot the figure
set(gca,'FontSize',14,'LineWidth',1.5,'Box','off')

% set the line width of the box plots
set(findobj(gca,'type','line'),'linew',2)
set(findobj(gca,'Tag','Lower Whisker'),'LineStyle','-')
set(findobj(gca,'Tag','Upper Whisker'),'LineStyle','-')

% set the color of the box plots
h4 = findobj(gca,'Tag','Box') ;
for kk=1:length(h4)
    patch(get(h4(kk),'XData'),get(h4(kk),'YData'),...
        color(kk,:),'FaceAlpha', 0.3,'LineStyle','-');
end

% add a scatter plot if we that is true
if includeScatter
   
    % get the unique groups 
    uniqueGroups = unique(groups) ;

    % add the scatter plots
    hold on
    
    % add scatter plot to the box plots
    groupSc1 = plotData(groups == uniqueGroups(1,1));
    groupSc2 = plotData(groups == uniqueGroups(2,1));
    
    % set the size the size of the marker
    makerSz = 30;
     
    % get only a smaller subset if the data points are too many
    if length(groupSc1) > 600
        groupSc1 = randsample(groupSc1,600) ;
        
        % change the markerSize 
        makerSz = 15 ;
    end
    
    if length(groupSc2) > 600
        groupSc2 = randsample(groupSc2,600) ;
    end
    
    x = ones(length(groupSc1)).*(1+(rand(length(groupSc1))-0.5)/5) ;
    x1 = ones(length(groupSc2)).*(1+(rand(length(groupSc2))-0.5)/10);
    
    % here is the first scatter plot
    scatter(x(:,1), groupSc1, makerSz, color(2,:),'filled', ...
        'MarkerFaceColor',color(2,:),'Marker','o','MarkerFaceAlpha',0.8)
    hold on
    scatter(x1(:,2).*2, groupSc2, makerSz, color(1,:),'filled', ...
        'MarkerFaceColor',color(1,:),'Marker','o','MarkerFaceAlpha',0.8)
    
    hold off
end

end

%% Another internal function 

function ManhattanPlot( filename, varargin )

%   This function takes a GWAS output file  and plots a Manhattan Plot.

%%  ARGUMENTS:
%   sex: defaults to 0, set to 1 to include sex chromosomes

%   sig: defaults to 5e-8, significance threshold for horizontal line

%   vert: defaults to 0, set to 1 for tick labels have chr in them and are
%   written upwards

%   labels: defaults to [-1,-1]. Set to [x,y] to label the top SNP on each
%   locus with p<x. Locus defined by windows of y base pairs.

%   outfile: defaults to filename of input file, set to something else to
%   change name of output file

%   title: defaults to 'Manhattan Plot'. Fairly self-explanatory

%   save: defaults to 1, set to 0 to disable saving to save time

%   format: defaults to PLINK, options {PLINK,BOLT-LMM,SAIGE} This is used
%   to identify the correct column headings. If using anything else, rename
%   the header line so that CHR, BP, and P are the headers for chromosome,
%   base pair and p value, and the default PLINK should catch it

%%   Usage:
%   ManhattanPlot('gwas.assoc.fisher',varargin) will take the association
%   analysis in 'gwas.assoc.fisher, generate a Manhattan Plot, and store it
%   in gwas_ManhattanPlot.png, which should be publication-ready, and
%   gwas_ManhattanPlot.fig for minor readjustments in MATLAB's GUI. This is
%   fine as .fisher is a PLINK format file. For a SAIGE output (and
%   significance threshold 5e-9) with sex chromosomes shown, and the top
%   hit for each SNP within 1 mb and below p=1e-6 labelled, use
%   ManhattanPlot('gwas.saige','format','SAIGE','sig',5e-9,'sex',1,
%   'labels',[1e-6,1000000])

%   Tested using assoc, assoc.logistic and assoc.fisher files generated by
%   Plink 1.7. Also tested using BOLT-LMM and SAIGE.
%
%   Harry Green (2019) Genetics of Complex Traits, University of Exeter

%% Reading in optional arguments and input file

p = inputParser;
defaultsex = 0;
defaultvert = 0;
defaultsig = 5e-8;
defaultsave = 1;
defaultoutfile = filename;
defaulttitle= 'Manhattan Plot';
defaultlabels = [-1, -1];
defaultformat = 'PLINK';
expectedformat = {'PLINK','BOLT-LMM','SAIGE'};

addRequired(p,'filename');
% addOptional(p,'outfile',defaultoutfile,@ischar);
addOptional(p,'sex',defaultsex,@isnumeric);
addOptional(p,'vert',defaultvert,@isnumeric);
addOptional(p,'labels',defaultlabels,@isnumeric);
addOptional(p,'sig',defaultsig,@isnumeric);
addOptional(p,'save',defaultsave,@isnumeric);
addParameter(p,'format',defaultformat,...
    @(x) any(validatestring(x,expectedformat)));
addParameter(p,'outfile',defaultoutfile,@ischar);
addParameter(p,'title',defaulttitle,@ischar);

parse(p,filename,varargin{:});

filename=p.Results.filename;
sex=p.Results.sex;
vert=p.Results.vert;
sig=p.Results.sig;
save=p.Results.save;
format=p.Results.format;
outfile=p.Results.outfile;
labels=p.Results.labels;
plottitle=p.Results.title;

% ################ Check if the input is a table for a file name #########

if ischar(filename)
    % MATLAB refuses to read file formats like .assoc, so first rename as
    % .txt
    try
        copyfile(filename, strcat(filename,'.txt'));
        opts = detectImportOptions(strcat(filename,'.txt'),...
            'NumHeaderLines',0);
        % the columns of T will be the headers of the assoc file
        T = readtable(strcat(filename,'.txt'),opts);
        
        %delete unwanted .txt file that we didn't want anyway
        delete(strcat(filename,'.txt'))
    catch
        T = readtable(filename);
    end
elseif istable(filename)
    T = filename ;
end

%% File format defintitions
% This section uses the default column headings from PLINK, BOLT-LMM and
% SAIGE. Easy to modify for other software packages

% check the CHR variable exist then change to PLINK format
if ~ismember(T.Properties.VariableNames,'CHR')
    T.Properties.VariableNames([1,2]) = {'CHR','BP'} ;
end

if strcmp(format,'PLINK')
    % These are the only useful columns for a manhattan plot, chromosome,
    % base pair, p value
    try 
        tab2=[T.CHR,T.BP,T.P]; 
    catch
        % some time this does not work
        locBP = find(ismember(T.Properties.VariableNames,'POS'), true);
        T.Properties.VariableNames(locBP) = "BP";
        tab2=[T.CHR,T.BP,T.P]; 
    end
end

if strcmp(format,'BOLT-LMM')
    % These are the only useful columns for a manhattan plot, chromosome,
    % base pair, p value
    tab2=[T.CHR,T.BP,T.P_BOLT_LMM]; 
end

if strcmp(format,'SAIGE')
    % These are the only useful columns for a manhattan plot, chromosome,
    % base pair, p value
    tab2=[T.CHR,T.POS,T.p_value]; 
end

%%

if sex==0
    tab2=tab2(tab2(:,1)<23,:); %remove chr23 if not wanting sex chromosomes
end

% variable to track base pairs, this helps the plotting function to know
% where to start plotting the next chromosome
bptrack=0; 
tab2(tab2(:,3)==0,:)=[];
for i=1:22+sex
    hold on
    %a scatterplot. On the x axis, the base pair number + bptrack, which
    %starts at 0. On the y axis, - log 10 p
    plot(tab2(tab2(:,1)==i,2)+max(bptrack), ...
        -log10(tab2(tab2(:,1)==i,3)), ...
        '.','MarkerSize',10); 
    
    %this updates bptrack by adding the highest base pair number of the
    %chromosome. At the end, we should get the total number of base pairs
    %in the human genome. All values of bptrack are stored in a vector.
    %They're useful later
    bptrack=[bptrack,max(bptrack)+max(max(tab2(tab2(:,1)==i,:)))]; 
end

%if strongest hit is 1e-60, plot window goes up to 1e-61
ylimit=max(max(-log10(tab2(:,3)))+1); 
ylimitmin=min(min(-log10(tab2(:,3)))); %and down to the highest p value.

%genome wide significant line, uses the sig optional argument which
%defaults to 5e-8
plot([0,max(bptrack)],-log10([sig,sig]),'k--') 

% ylabel('$-\log_{10}(p)$','Interpreter','latex')
ylabel('-log_1_0(p)')
xlim([0,max(bptrack)])
%y axis will always go to the significance threshold/5 at least
ylim([floor(ylimitmin),max(ylimit,-log10(sig/5))])

%this calculates the moving average of bptrack and uses them as chromosome
%label markers. This puts them in the middle of each chromosome.
M=movmean(bptrack,2); 
xticks(M(2:end));

%% Rotation section
% this section of the code changes the x axis label orientation.

if vert ==0
    xlabel('Chromosome')% ,'Interpreter','latex')
    xticklabels( 1:23 );
end
if vert ==1
    xtickangle(90)
    xticklabels( {'chr1','chr2','chr3','chr4','chr5','chr6','chr7',...
        'chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15',...
        'chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrXY'});
end

%% Annotation section
% this section of the code annotates the top SNP on each locus

% Neither of these should be negative. Defaults to -1 -1
if ~(labels(1)==-1||labels(2)==-1) 
    labellim=labels(1);
    labelprox=labels(2);
    
    %as nothing else will be used for labelling, the reduced table is
    %easier to manage
    tab3=tab2(tab2(:,3)<labellim,:); 
    
    %easier if they're in ascending order of p value;
    [~,index]=sort(tab3(:,3));
    tab3=tab3(index,:);
    
    %this becomes a list of SNPs that have been labelled. It starts with a
    %dummy SNP to give the first SNP something to compare with. This will
    %always fail so the first label always gets added
    labelledsnps=[0,0,0]; 
    
    for i=1:size(tab3,1)
        test=abs(tab3(i,:)-labelledsnps);
        %if there are no snps on the same chromosome and within 1 MB
        if sum(test(:,1)==0&test(:,2)<labelprox)==0 
            
            %this plots a square over the corresponding dot
            plot(bptrack(tab3(i,1))+tab3(i,2),-log10(tab3(i,3)),'ks',...
                'MarkerSize',10,'LineWidth',1) 
            
            %this puts together a strong of chrx:y
            labels=strcat('chr',string(tab3(i,1)),':',...
                string(bptrack(tab3(i,1))+tab3(i,2))); 
            
            %this plots the label above and to the left of the dot, shifted
            %up by 0.05 to give some space
            text(bptrack(tab3(i,1))+tab3(i,2),-log10(tab3(i,3))+0.05,...
                char(labels),'VerticalAlignment','bottom',...
                'HorizontalAlignment','left') %,'Interpreter','latex') 
            labelledsnps=[labelledsnps;tab3(i,:)];
        end
    end
end
%% finalising file

% takes the output file name up to the first decimal point
a = strsplit(outfile,'.'); 

% my code 
set(gca,'FontSize',12,'LineWidth',1)

box on
title(plottitle,'fontsize',14)
% title(plottitle,'Interpreter','latex','fontsize',14)

if save~=0
    name=[a{1},'_ManhattanPlot.fig']; %adds _ManhattanPlot to filenames
    savefig(name);
    set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 30 10])
    print([a{1},'_ManhattanPlot.png'],'-dpng','-r300', '-painters') % -r600
end



end
        
% ******************* Another Internal function ******************

function exitPos = ultraBars(inData, colors, rowNames, ...
    horizonBars , nextAxis)

% Input:
% inData: a matrix and vector of plotting data
% colors: colors for each unique value of inData
% rowNames: a cell array for name for each row
% horizonBar: specifies whether to add a horizontal bar to each plot or not
 
% get the number of unique numbers and remove the 0 which means not plot
uniqueVars = unique(inData) ;
uniqueVars(uniqueVars == 0) = [] ;

% create the legend variable if the inData is a row vector
if size(inData,1) == 1
    lgdVar = split( num2str(uniqueVars) )';
end

% get the number of colours to plot
if ~exist('colors','var') || isempty(colors)
    colors = rand(length(uniqueVars), 3);
end

% check the are is a row name for each row in inData
if size(inData,1) > 1 % only valid for matrix data
    if size(inData,1) ~= size(rowNames,1)
        error('row names should contain a name for each row in plotting data')
    end
end

% check if the orizontal bars are required on the plot
if ~exist('horizonBars','var')
    horizonBars = false;
end

% check for color errors 
if length(uniqueVars) > size(colors,1) 
    error('A color must be specified for each value in Data')
end

% plot the figure for row vectors 
% I have separated the two peaces of code because one does not need a box
% around the data
if size(inData,1) == 1
    % figure()
    % axes('position',[0.1300,0.11,0.7,0.04]);
    % set the bar width for large case
    if size(inData,2) < 500
        barwidth = 0.9 ;
    else
        barwidth = 1 ;
    end
    
    % now plot the data
    for ii = 1:length(uniqueVars)
        % get the data for that values and plot it in that color
        plotData = double(ismember(inData,uniqueVars(ii) )) ;
        bar(plotData,'FaceColor',colors(ii,:),...
            'EdgeColor',[1 1 1] ,'BarWidth',barwidth) ;
        hold on
    end
    set(gca,'GridColor',[1 1 1], 'XLim', [0.5 size(inData,2)+0.5 ], ...
        'XColor',[1 1 1] ,'YColor',[1 1 1],'XTick',[] , 'YTick',[],...
        'FontWeight','bold')
    % add a legend to the bar using the legendflex function
   
% if the is more thhan one row in the dataset
else
    % make a plot of multiple bar chart of top of each other
    % add the subtype plots to the heatmap using a loop
    % initialise some variables
    global plotTime
    figure(plotTime*100)
    % set(gcf,'position',[100,50,800,600])
        
    yInitial = 0.05; yPosSaved = yInitial;
    ySize = 0.02; increaseby = 0.022; % 0.44 0.4
    % change the size the bar if plotTime == 3
    if plotTime == 3
        ySize = 0.05; increaseby = 0.055;
    end
    xEndPos = 0.7 ; % 0.7750
    % loop over the rows and ascend by column
    for jj = 1:size(inData,1) % begin plotting from the bottem
        % define the next axis for the follow up plots
        if ~exist('nextAxis','var') || isempty(nextAxis)
            axes('position',[0.1300,yInitial,xEndPos,ySize]);
        elseif exist('nextAxis','var') && jj > 1
            axes('position',[0.1300,yInitial,xEndPos,ySize]);   
        else
            axes('position',nextAxis);
            yInitial = nextAxis(2) ;
            ySize = nextAxis(4) ; xEndPos = nextAxis(3) ;
            yPosSaved = yInitial ;
        end         
        for ii = 1:numel(uniqueVars) 
            plotData = double(ismember(inData(jj,:),uniqueVars(ii) )) ;
            bar(plotData,'FaceColor',colors(ii,:),'EdgeColor',[1 1 1] ,...
                'BarWidth',0.9) ;   
            hold on
            
            % add the name of the genes to the left of heatmap
            if exist('rowNames','var')
                dim = [0.02 yInitial 0.11 increaseby];
                annotation('textbox',dim,'String',rowNames{jj},...
                'FitBoxToText','on','FontSize',10,'EdgeColor','none',...
                    'HorizontalAlignment','right','FontWeight','bold',...
                    'VerticalAlignment','middle');
            end
        end
        % change the plot properties
        set(gca,'GridColor',[1 1 1], 'XLim', [0.5  size(inData,2)+0.5],...
            'XColor',[1 1 1] ,'YColor',[1 1 1],'YTickLabel',[],...
            'XTickLabel',[],'FontWeight','bold','YTick',[],'XTick',[])
        % increase the value to change the colors and plot positions
        yInitial = yInitial + increaseby;
    end
    % add a grey box to the plot usign the annotation function
    if plotTime ~= 3 % dont add the box if plotTime is == 3
    dim = [0.1300, yPosSaved, xEndPos, increaseby*size(inData,1)];
    annotation('rectangle',dim ,'Color',[0.5, 0.5, 0.5])
    end
    hold off 
   
    % plot the horizontal bars if they are required
    % prellocate the bar data size
    barhData = zeros(size(inData,1),numel(uniqueVars)) ;
    if horizonBars == true
        axes('position',[xEndPos+0.137, yPosSaved, 0.12, ...
            increaseby*size(inData,1)+0.001 ]);
        for kk = 1:numel(uniqueVars) 
            barhData(:,kk) = sum(inData == uniqueVars(kk),2) ;   
        end
        bar1 = barh(barhData,'stacked','BarWidth',0.85) ;
        % make sure there are no colors and spaces between the axis and the
        % first and last bar
        set(gca,'GridColor',[1 1 1], 'YLim', [0.5 size(inData,1)+0.52 ], ...
            'XColor',[0.5 0.5 0.5] ,'YColor',[1 1 1],'FontSize',10,...
            'YTick',[],'FontWeight','bold','Box','off','TickDir', 'out',...
            'LineWidth',1,'XAxisLocation','origin')
        
        % annoate the bar graph
        for ii = 1:size(barhData,2)
            set(bar1(ii),'FaceColor',colors(ii,:))
        end
        
    end
    exitPos = [0.1300,yInitial,xEndPos, ySize] ;
end

end

% ************************ Another Internal Function ******************

% function for LD cross referecing when comparing variants
function [allLd, snpData, numOfLDsnps ] = ld_crossref(ldData, snpData, ...
    rSquareCutOff , wndwSize)

% the function has five (5) inputs
% #1 LD table for UKB or 1000G
% #2 snpData - are the results of a gwas study for causal snps
% #3 allSnpsData - all the snps not only the causal variants - NOT INCLUDED
% #4 rSquareCutOff - values for considering LD - e.g., 0.5
% #5 window size - up or downstream the causal variant

% get the current ld data and the curnspl data
curLddata = ldData(ldData.Rsquare >= rSquareCutOff, :) ;
allLd = [] ;
numOfLDsnps = 0 ;
seenVariants = {'Dummy'};

% here are the data
for ii = 1:height(snpData)

    % print something to the screen
    if rem(ii,100) == 0
        fprintf('\nrunning analysis for variant #%d of %d\n', ...
            ii, height(snpData) )
    end

    % get the current snp
    curSnp = snpData.Variant(ii) ;

    % check if the variants has already been seen
    if ~isempty(seenVariants)
        if any( ismember(seenVariants,curSnp) )
            continue
        end
    end

    % add the current variant to those we have already seen
    seenVariants = unique( [seenVariants;curSnp] );

    % get the snps in close proximity with the current snp -- these are
    % snps which are within 500KB of the lead snp that is predicted
    % caused -- THESE SHOULD BE TAKEN FROM ALL THE SNPS AVAILABLE IN
    % THE DATA
    closeSnps = snpData( ismember( snpData.Position , ...
        snpData.Position(ii)-wndwSize:snpData.Position(ii)+wndwSize) & ...
        snpData.Chrom == snpData.Chrom(ii) , : ) ;

    % remove the current snp from the table
    closeSnps(ismember(closeSnps.Variant, curSnp), :) = [] ;

    % check if the table has data
    if isempty(closeSnps)
        continue
    end

    % find if any of those snps are in LD with the lead snp AND ALL THE
    % OTHER SNPS THAT ARE SIGNFICANT WITHIN A PARTICULAR REGION
    curLd = curLddata( ismember(curLddata.IndexSNP,curSnp), :) ;

    % find the snps are in strong LD
    curLd = curLd(ismember( curLd.LinkedSNP ,closeSnps.Variant), :) ;
    curLd = unique(curLd) ;

    % check if the table has data
    if isempty(curLd)
        continue
    end

    % get the variants in ld
    ldVariants = unique( [curLd.LinkedSNP] );

    % add the current variant to those we have already seen
    seenVariants = unique( [seenVariants;ldVariants] );

    % create the current ld struct
    allLd = [allLd ;curLd] ;

    % get the location of those viariants
    locVariants = ismember(snpData.Variant, ldVariants);

    % print some thing the screen
    numOfLDsnps = numOfLDsnps + sum(locVariants) -1  ;

    % change those snps in the snpLd table
    snpData.Variant(locVariants) = curSnp ;
    snpData.Chrom(locVariants) = snpData.Chrom(ii);
    snpData.Position(locVariants) = snpData.Position(ii) ;
    try
        snpData.Alternative(locVariants) = snpData.Alternative(ii);
    catch
    end
    try
        snpData.Reference(locVariants) = snpData.Reference(ii);
    catch
    end

%     % *******************************************************************
%     % brute force approach as well
% 
%     % get the snps with 250kb of the lead snp
%     closeSnps = snpData( ismember( snpData.Position , ...
%         snpData.Position(ii)-500000:snpData.Position(ii)+500000) & ...
%         snpData.Chrom == snpData.Chrom(ii) , : ) ;
% 
%     % covert those snps to the lead snp
%     locVariants = ismember(snpData.Variant,closeSnps.Variant);
%     snpData.Variant(locVariants) = curSnp; 
%     snpData.Chrom(locVariants) = snpData.Chrom(ii);
%     snpData.Position(locVariants) = snpData.Position(ii) ;
%     try
%         snpData.Alternative(locVariants) = snpData.Alternative(ii);
%     catch
%     end
%     try
%         snpData.Reference(locVariants) = snpData.Reference(ii);
%     catch
%     end
%     % ******************************************************************
% 
%     % print some thing the screen
%     locVariants =  ismember(snpData.Variant,closeSnps.Variant) ;
%     numOfLDsnps = numOfLDsnps + sum(locVariants) -1  ;
    fprintf('\nThe number of LD variant identified is %d\n',numOfLDsnps)

end
end
