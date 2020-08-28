%%% Creates plots of USA's COVID cases, because why not
%%% 
%%% Ideas for future enhancements:
%%% 1. More sophisticated model of Prevalence Ratio that makes it a 
%%%    function of both test positivity rate and number of tests per capita
%%%    and allows the [a b c] parameters to vary over time: 
%%%    https://covid19-projections.com/estimating-true-infections/
%%%
%%% 2. County level and Nationwide plots that incorporate Prevalance Ratio
%%% 
%%% 3. Incorporating estimates of unreported COVID deaths: 
%%%    https://weinbergerlab.github.io/excess_pi_covid/
%%% 
%%% 4. Quantififying the amount of correlation between shutdowns/mask
%%%    mandates/etc and number of new infections and/or R_t. This seems
%%%    hard...and FUN!
%%%


clearvars; close all;

answer = questdlg('Create new figures? This takes longer to run',  ...
    '', 'Yes', 'Just State Level', 'No', 'Yes');

switch answer
    case 'Yes'
        makeFigsFlag = true;
        makeCountyFigsFlag = true;
    case 'No'
        makeFigsFlag = false;
        makeCountyFigsFlag = false;
    case 'Just State Level'
        makeFigsFlag = true;
        makeCountyFigsFlag = false;
    otherwise
        return
end

%Load some data. Can download new data from these websites:
url = 'https://usafactsstatic.blob.core.windows.net/public/data/covid-19/';
url2 = 'https://covidtracking.com/api/v1/states/';
url3 = 'https://d14wlfuexuxgcm.cloudfront.net/covid/';
theRoot = ['https://raw.githubusercontent.com/youyanggu/' ...
    'covid19_projections/master/r_values/'];
restOfName = '_r_values_us.csv';
popFile = 'covid_county_population_usafacts.csv';
caseFile = 'covid_confirmed_usafacts.csv';
deathFile = 'covid_deaths_usafacts.csv';
testFile = 'daily.csv';
RTfile = 'rt.csv';
web([url popFile], '-browser');
web([url caseFile], '-browser');
web([url deathFile], '-browser');
web([url2 testFile], '-browser');
web([url3 RTfile], '-browser');
wait_a_sec = msgbox('Downloading some data. Wait a sec then click OK');
uiwait(wait_a_sec);

dataLoc = 'C:\Users\schei\Downloads\';
mainDir = 'C:\Users\schei\OneDrive\Documents\MATLAB\COVID\';
newLoc = [mainDir 'data\'];

%Put deaths in perspective with other disastrous events
benghazi = 4;
titanic = 1517;
katrina = 1836;
pearl = 2467;
nine11 = 2996;
gettysburg = 3155;
lung = 135720 / 12;
flu = sum(1000 * [37 12 43 38 51 23 38 61 34.157]) / 9;
heart = 647000 / 12;
hiroshima = 75000;

%Extract the pertinent data from COVID-Projections RT model
infectionLength = 15;  %Average duration in days of COVID-19
daysPerMonth = [31 29 31 30 31 30 31 31 30 31 30 31]; %2020 has a Leap Day
dayOfYear = 145;  %For offsetting "days since 2020 began" in plots
firstMonth = 5;  %First file predicting R_t is May 24th
firstDay = 24;
newFormatMonthDay = [7 22];  %First file predicting R_t in different format
commasPerLine = 5;

%Need to make the covid_projections.com R_t file string
currMonth = firstMonth;
currDay = firstDay;
iter = 1;
newFormatFlag = false;
notDoneFlag = true;
while notDoneFlag
    if daysPerMonth(currMonth) < currDay
        currDay = 1;  %We have moved to a new month
        currMonth = currMonth + 1;
        if currMonth == 13  %A new year! script needs updating
            keyboard
        end
    end
    
   %See if we have reached the new file format on July 22
   if currMonth >= newFormatMonthDay(1)
       if currDay >= newFormatMonthDay(2)
           newFormatFlag = true;
           commasPerLine = 7;
       end
   end
    
   %Build the date string for the file we will be reading
    dayString = num2str(currDay);
    if length(dayString) < 2  %must be in days 1-9, add 0 to string
        dayString = ['0' dayString];
    end
    
    monthString = num2str(currMonth);
    if length(monthString) < 2  %must be in months 1-9, add 0 to string
        monthString = ['0' monthString];
    end
    
    dateString = ['2020-' monthString '-' dayString];
    url4 = [theRoot dateString restOfName];
    
    %Use a try/catch in case we are trying to read a file that doesn't
    %exist yet
    try
        RTstring = webread(url4);
    catch  %No more files to read in!
        notDoneFlag = false;
        continue;
    end
    
    commas = strfind(RTstring, ',');
    
    %Grab the R_t data
    if newFormatFlag
        for nextLine = commasPerLine:commasPerLine: ...
                length(commas) - commasPerLine
            lineNum = nextLine / commasPerLine;
            nextCommaInd = commas(nextLine+1);
            if ~exist('RTstateArray')  %Only need to grab states list once
                stateString = RTstring(nextCommaInd-2 : nextCommaInd-1);
                RTstateArray{lineNum} = stateString;
            end
            lastCommaInd = commas(nextLine + commasPerLine);
            theRTs{lineNum}(iter) = str2num(RTstring(lastCommaInd+1 : ...
                lastCommaInd+4));
        end
    else  %Old format with different R_t column
        for nextLine = commasPerLine:commasPerLine: ...
                length(commas) - commasPerLine
            lineNum = nextLine / commasPerLine;
            nextCommaInd = commas(nextLine+1);
            stateString = RTstring(nextCommaInd-2 : nextCommaInd-1);
            RTstateArray{lineNum} = stateString;
            lastCommaInd = commas(nextLine + commasPerLine);
            theRTs{lineNum}(iter) = str2num(RTstring(lastCommaInd-4 : ...
                lastCommaInd-1));
        end
    end
    
    iter = iter + 1;
    currDay = currDay + 1;
end

%Read other data now that it's downloaded
population = readmatrix([dataLoc popFile]);
locations = importdata([dataLoc popFile]);
cases = readmatrix([dataLoc caseFile]);
deaths = readmatrix([dataLoc deathFile]);
tests = readmatrix([dataLoc testFile]);
testsStates = importdata([dataLoc testFile]);
RTlive = importdata([dataLoc RTfile]);

%Move data to our MATLAB work area
movefile([dataLoc popFile], [newLoc popFile]);
movefile([dataLoc caseFile], [newLoc caseFile]);
movefile([dataLoc deathFile], [newLoc deathFile]);
movefile([dataLoc testFile], [newLoc testFile]);
movefile([dataLoc RTfile], [newLoc RTfile]);

%remove date header row from cases/deaths and county/state
%columns from population
offsetDays = 21;  %Because first date in some spreadsheets is January 22nd
cases = cases(2:end, :);
deaths = deaths(2:end, :);
population(:, 2:3) = [];
counties = locations.textdata(2:end, 2);
states = locations.textdata(2:end, 3);
testsStates = flipud(testsStates.textdata(2:end, 2));
RTliveStates = RTlive.textdata(2:end,2);
RTliveValues = RTlive.data(:,2);


%See if county identifier columns match
if ~isequal(population(:,1), cases(:,1), deaths(:,1))
    %Data needs to be inspected more carefully
    for iter = 1:length(population(:,1))
        %Find where there is a mismatch
        if ~isequal(population(iter,1), cases(iter,1), deaths(iter,1))
            %Sometimes the population file records statewide as a 1 instead
            %of a zero. If both cases/deaths match, make population match
            %too.
            if isequal(cases(iter,1), deaths(iter,1))
                population(iter,1) = cases(iter,1);
            else
                keyboard;  %Manual inspection required
            end
        end
    end
end

%Remove extraneous columns now that data is aligned
cases = cases(:, 5:end);
deaths = deaths(:, 5:end);
population = population(:, 2);
tests = flipud(tests(:, [1 3 4]));


%Have to address the weird New Jersey artifact on July 17th of excess
%deaths. Choosing to just ignore it for now and copy the data from July
%16th instead.
theDiff = deaths(1807:1827, 158) - deaths(1807:1827, 157);
deaths(1807:1827, 158:end) = deaths(1807:1827, 158:end) - theDiff;


%Calculate per million multiplier and fix the divide by zeros
multiplier = 10^6 ./ population;
fixIndices = isinf(multiplier);
multiplier(fixIndices) = 0;


%Initialize arrays to store all the data and a couple counters
dimensions = size(cases);
casesPerDay = zeros(dimensions(1), dimensions(2));
deathsPerDay = zeros(dimensions(1), dimensions(2));
cases7day = zeros(dimensions(1), dimensions(2));
deaths7day = zeros(dimensions(1), dimensions(2));
cases7dayPerMil = zeros(dimensions(1), dimensions(2));
deaths7dayPerMil = zeros(dimensions(1), dimensions(2));
casesPerDayStatewide = zeros(51, dimensions(2));
deathsPerDayStatewide = zeros(51, dimensions(2));
cases7dayStatewide = zeros(51, dimensions(2));
deaths7dayStatewide = zeros(51, dimensions(2));
cases7dayPerMilStatewide = zeros(51, dimensions(2));
deaths7dayPerMilStatewide = zeros(51, dimensions(2));
deathsStatewideCumul = zeros(51, dimensions(2));
estimatedPercentInfected = zeros(51, 1);
previousState = states{1};
stateOrder = cell(51,1);
stateBeginIndex = 1;
currState = 0;


%Create 7 day averages for cases and deaths per million population
%Have to loop through each county
for currCounty = 1:dimensions(1)
    %Set up a waitbar
    if exist('wait', 'var')
        waitbar(currCounty/dimensions(1), wait);
    else
        wait = waitbar(currCounty/dimensions(1));
    end
    
    %cases/deaths are both cumulative, so need to figure out per day stats
    casesPerDay(currCounty,:) = ...
        [cases(currCounty,1) diff(cases(currCounty,:))];
    deathsPerDay(currCounty,:) = ...
        [deaths(currCounty,1) diff(deaths(currCounty,:))];
    
    %Have to loop through each day
    for currDay = 1:dimensions(2)
        %Average everything if haven't hit 7th day yet
        if currDay <= 6
            cases7day(currCounty, currDay) = ...
                sum(casesPerDay(currCounty, 1:currDay) / currDay);
            deaths7day(currCounty, currDay) = ...
                sum(deathsPerDay(currCounty, 1:currDay) / currDay);
            
        else   %otherwise just average most recent 7 days
            cases7day(currCounty, currDay) = ...
                sum(casesPerDay(currCounty, currDay-6:currDay) / 7);
            deaths7day(currCounty, currDay) = ...
                sum(deathsPerDay(currCounty, currDay-6:currDay) / 7);
        end
    end
    
    %Create statewide plot if we've iterated to a new state
    newState = ~matches(previousState, states{currCounty});
    if (newState || currCounty == dimensions(1)) && makeFigsFlag
        currState = currState + 1;
        stateOrder{currState} = previousState;
        if stateBeginIndex ~= (currCounty-1)  %This means we aren't DC
            casesPerDayStatewide(currState, :) = ...
                sum(casesPerDay(stateBeginIndex:currCounty-1, :));
            deathsPerDayStatewide(currState, :) = ...
                sum(deathsPerDay(stateBeginIndex:currCounty-1, :));
            cases7dayStatewide(currState, :) = ...
                sum(cases7day(stateBeginIndex:currCounty-1, :));
            deaths7dayStatewide(currState, :) = ...
                sum(deaths7day(stateBeginIndex:currCounty-1, :));
            cases7dayPerMilStatewide(currState, :) = ...
                sum(cases7dayPerMil(stateBeginIndex:currCounty-1, :));
            deaths7dayPerMilStatewide(currState, :) = ...
                sum(deaths7dayPerMil(stateBeginIndex:currCounty-1, :));
            deathsStatewideCumul(currState,:) = ...
                sum(deaths(stateBeginIndex:currCounty-1,:));
            statePop = sum(population(stateBeginIndex:currCounty-1));
            multiplierStatewide = 10^6 / statePop;
        else  %Otherwise we're in Washington DC
            casesPerDayStatewide(currState, :) = ...
                casesPerDay(stateBeginIndex:currCounty-1, :);
            deathsPerDayStatewide(currState, :) = ...
                deathsPerDay(stateBeginIndex:currCounty-1, :);
            cases7dayStatewide(currState, :) = ...
                cases7day(stateBeginIndex:currCounty-1, :);
            deaths7dayStatewide(currState, :) = ...
                deaths7day(stateBeginIndex:currCounty-1, :);
            cases7dayPerMilStatewide(currState, :) = ...
                cases7dayPerMil(stateBeginIndex:currCounty-1, :);
            deaths7dayPerMilStatewide(currState, :) = ...
                deaths7dayPerMil(stateBeginIndex:currCounty-1, :);
            deathsStatewideCumul(currState,:) = ...
                deaths(stateBeginIndex:currCounty-1,:);
            statePop = sum(population(stateBeginIndex:currCounty-1));
            multiplierStatewide = 10^6 / statePop;
        end
        
        
        %Calculate the cumulative cases and deaths per million population
        casesPerMilStatewide = round(multiplierStatewide * ...
            sum(casesPerDayStatewide(currState,:)), -2);
        deathPerMilStatewide = round(multiplierStatewide * ...
            sum(deathsPerDayStatewide(currState,:)), 1);
        
        
        %Calculate test positivity rate
        testInds = matches(testsStates, previousState);
        currTests = tests(testInds, :);
        %Remove rows without a tally for negative tests
        noNegativesInds = isnan(currTests(:,3));
        currTests(noNegativesInds,:) = [];
        positivityRate = currTests(:,2) ./ ...
            (currTests(:,2) + currTests(:,3));
        positivityRate(isnan(positivityRate)) = 0;  %Replace divide-by-0s
        
        %Have to loop through each day to calculate 7 day average
        positivityRate7day = zeros(length(positivityRate), 1);
        dayNumber = zeros(length(positivityRate), 1);
        confirmedCases = zeros(length(positivityRate), 1);
        for currentDay = 1:length(positivityRate)
            %Average everything if haven't hit 7th day yet
            if currentDay <= 6
                positivityRate7day(currentDay) = ...
                    sum(positivityRate(1:currentDay) / currentDay);
            else   %otherwise just average most recent 7 days
                positivityRate7day(currentDay) = ...
                    sum(positivityRate(currentDay-6:currentDay) / 7);
            end
            
            %Figure out how many days from January 21st we are:
            monthDay = currTests(currentDay,1) - 20200000;
            if monthDay > 1231 
                keyboard;  %Need to adjust code to accept non-2020 years!
            end
            
            numMonths = round(monthDay/100);
            if numMonths > 1  %Need to add previous months days
                dayNumber(currentDay) = mod(monthDay, 100) + ...
                    sum(daysPerMonth(1:numMonths-1));
            else
                dayNumber(currentDay) = mod(monthDay, 100);
            end
            
            %Grab confirmed cases 7 day average for this day and state
            theDay = dayNumber(currentDay) - offsetDays;
            if theDay <= size(cases7dayStatewide, 2)
                confirmedCases(currentDay) = ...
                    cases7dayStatewide(currState, theDay);
            else  %Don't have new versions of the other spreadsheets
                %Just assume same as previous day
                confirmedCases(currentDay) = ...
                    cases7dayStatewide(currState, theDay-1);
            end
        end
        
        %Calculate estimate of TRUE number of cases as a function of 7 day
        %average of test positivity rate, as described here:
        %https://covid19-projections.com/estimating-true-infections
        MayInd = find(dayNumber > sum(daysPerMonth(1:4)), 1);
        a1 = 16;
        b1 = 0.5;
        c1 = 2.5;
        prevalenceRatioEarly = a1 .* ...
            (positivityRate7day(1:(MayInd-1))) .^ b1 + c1;
        a2 = 10;
        b2 = 0.4;
        c2 = 2.5;
        prevalenceRatioLater = a2 .* ...
            (positivityRate7day(MayInd:end)) .^ b2 + c2;
        prevalenceRatio = [prevalenceRatioEarly; prevalenceRatioLater];
        
        estimatedCases = prevalenceRatio .* confirmedCases;
        estimatedPercentInfected(currState) = ...
            round(100 * (sum(estimatedCases) / statePop), 1);
        
        %Calculate estimate of percent currently infected
        estimatedCurrentlyInfected = zeros(length(estimatedCases), 1);
        
        for dayIter = 1:length(estimatedCases)
            if dayIter <= infectionLength
                estimatedCurrentlyInfected(dayIter) = 100 * ...
                    (sum(estimatedCases(1:dayIter) / statePop));
            else
                estimatedCurrentlyInfected(dayIter) = 100 * ...
                    (sum(estimatedCases(dayIter-infectionLength:dayIter)...
                    / statePop));
            end
        end
        disp([previousState ': ' ...
            num2str(round(estimatedCurrentlyInfected(end), 2)) ...
            '% Currently Infected']);
        
        %Grab the R_t data for this state
        casesInds = strcmp(previousState, RTliveStates);
        RTcases = RTliveValues(casesInds);
        deathsInd = strcmp(previousState, RTstateArray);
        RTdeaths = theRTs{deathsInd};
        
        theFig = tiledlayout(2, 2 , 'TileSpacing', 'Compact');
        title(theFig, ['COVID-19 Stats for the Great State of ' ...
            previousState], 'FontSize', 16);
        
        
        nexttile;
        grid on; hold on; box on;
        plot((1:currDay)+offsetDays, casesPerDayStatewide(currState,:), ...
            'Color', [0.9290 0.6940 0.1250]);
        plot((1:currDay)+offsetDays, cases7dayStatewide(currState,:), ...
            'LineWidth', 4, 'Color', [0.9290 0.6940 0.1250]);
        plot(dayNumber, estimatedCases, 'LineWidth', 4, ...
            'Color', [0.4940 0.1840 0.5560]);
        xlabel('Days Since 2020 Began');
        ylabel('New COVID-19 Cases');
        legend(['New Confirmed Cases: ' num2str(casesPerMilStatewide) ...
            ' per million cumulatively'], ...
            '7 Day Case Average', ['Estimated 7 Day Avg: ' ...
            num2str(estimatedPercentInfected(currState)) ...
            '% of population have had COVID-19'], 'Location', 'northwest');
        xLims = get(gca, 'XLim');
        
        nexttile;
        grid on; hold on; box on;
        yyaxis left;
        plot(dayNumber, estimatedCurrentlyInfected, 'LineWidth', 4, ...
            'Color', [0.4940 0.1840 0.5560]);
        infectedNow = estimatedCurrentlyInfected(end);
        ylabel('Estimated % Infected at the Time', ...
            'Color', [0.4940 0.1840 0.5560]);
        yyaxis right;
        plot((1:currDay)+offsetDays, deathsStatewideCumul(currState,:), ...
            'LineWidth', 4, 'Color', [0.8500 0.3250 0.0980]);
        xlabel('Days Since 2020 Began');
        ylabel('# of Deaths', 'Color', [0.8500 0.3250 0.0980]);
        legend(['% of Population Infected: Currently ' ...
            num2str(round(infectedNow,2)) '%'], ...
            ['Cumulative Deaths: ' num2str(deathPerMilStatewide) ...
            ' per million so far'], 'Location', 'northwest');
        xlim(xLims);

        nexttile;
        grid on; hold on; box on;
        xlabel('Days Since 2020 Began');
        RTcasesDays = dayNumber(end)-length(RTcases)+1 : dayNumber(end);
        plot(RTcasesDays, RTcases, 'LineWidth', 4, ...
            'Color', [0.9290 0.6940 0.1250]);
        RTdeathsDays = dayNumber(end)-length(RTdeaths)+1 : dayNumber(end);
        plot(RTdeathsDays, RTdeaths, 'LineWidth', 4, ...
            'Color', [0.8500 0.3250 0.0980]);
        RTcasesNow = round(RTcases(end), 2);
        RTdeathsNow = round(RTdeaths(end), 2);
        ylabel('Effective Reproduction Number (R_t)');
        legend(['R_t Estimate from Cases: Currently ' ...
            num2str(RTcasesNow)],['R_t Estimate from Deaths: Currently '...
            num2str(RTdeathsNow)], 'Location', 'northeast');
        xlim(xLims);
        
        nexttile;
        grid on; hold on; box on;
        theCDF = 100 * (1 - binocdf(0, 1:1000, infectedNow/100));
        firstOver90 = find(theCDF > 90, 1);
        if ~isempty(firstOver90)  %Want to plot up to 90% chance
            plot(1:firstOver90, theCDF(1:firstOver90), 'LineWidth', 4);
        else
            plot(1:1000, theCDF, 'LineWidth', 4);
        end
        xlabel(['# of People Interacted With (Assumes ' ...
            num2str(round(infectedNow,2)) '% Infected)']);
        ylabel('% Chance You Interacted With Infected Person');
        
        %Make figs fullsize and save figs and jpgs
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gcf, 'Visible', 'on');
        exportgraphics(gcf, ...
            [mainDir '\jpgs\' previousState '\(STATEWIDE).jpg']);
        savefig([mainDir '\figs\' previousState '\(STATEWIDE).fig']);

        close all;
        
        %Update beginning index for the new state
        stateBeginIndex = currCounty;
        previousState = states{currCounty};
    end
    
    
    %Calculate the cumulative cases and deaths per million population
    casesPerMil = round(multiplier(currCounty) * ...
        sum(casesPerDay(currCounty,:)), -2);
    deathPerMil = round(multiplier(currCounty) * ...
        sum(deathsPerDay(currCounty,:)), 1);
    
    %Generate and save county plots. Don't create plots if "Unallocated"
    unallocatedIndex = strfind(counties{currCounty}, 'Unallocated');
    if ~isempty(unallocatedIndex) || ~makeFigsFlag || ~makeCountyFigsFlag
        continue
    end
    figure('Visible', 'off'); grid on; hold on;
    plot((1:currDay)+offsetDays, casesPerDay(currCounty,:), ...
        'Color', [0.9290 0.6940 0.1250]);
    plot((1:currDay)+offsetDays, deaths(currCounty,:), 'LineWidth', 4, ...
        'Color', [0.8500 0.3250 0.0980]);
    plot((1:currDay)+offsetDays, cases7day(currCounty,:), ...
        'LineWidth', 4, 'Color', [0.9290 0.6940 0.1250]);
    xlabel('Days Since 2020 Began');
    ylabel('New COVID-19 Cases');
    title(['COVID-19 Stats For ' counties{currCounty} ', ' ...
        states{currCounty}]);
    legend(['New Confirmed Cases: ' num2str(casesPerMil) ...
        ' per million'], ['Cumulative Deaths: ' num2str(deathPerMil) ...
        ' per million'], '7 Day Case Average', 'Location', 'northwest');
    saveString = strrep(counties{currCounty}, ' ', '');
    saveString = strrep(saveString, '.', '');
    
    %Don't need to save a County fig for Washington DC
    if ~matches('DC', states{currCounty})
        set(gcf, 'Visible', 'on');
        savefig([mainDir '\figs\' states{currCounty} '\' saveString]);
        exportgraphics(gcf, ...
            [mainDir '\jpgs\' states{currCounty} '\' saveString '.jpg']);
    end
    
    close all;
end

if makeFigsFlag
    %Create national plots
    figure('Visible', 'off'); grid on; hold on;
    plot((1:currDay)+offsetDays, sum(casesPerDay), ...
        'Color', [0.9290 0.6940 0.1250]);
    plot((1:currDay)+offsetDays, sum(deaths), 'LineWidth', 4, ...
        'Color', [0.8500 0.3250 0.0980]);
    plot((1:currDay)+offsetDays, sum(cases7day), 'LineWidth', 4, ...
        'Color', [0.9290 0.6940 0.1250]);
    xlabel('Days Since 2020 Began');
    ylabel('New COVID-19 Cases');
    title('COVID-19 Stats for the United States of America');
    multiplierNation = 10^6 / sum(population);
    casesPerMilNation = round(multiplierNation * sum(cases(:,end)), -2);
    deathPerMilNation = round(multiplierNation * sum(deaths(:,end)));
    legend(['New Confirmed Cases: ' num2str(casesPerMilNation) ...
        ' per million'], ['Cumulative Deaths: ' ...
        num2str(deathPerMilNation) ' per million'], ...
        '7 Day Case Average', 'Location', 'northwest');
    set(gcf, 'Visible', 'on');
    savefig([mainDir '\figs\NATIONWIDE_CASES.fig']);
    exportgraphics(gcf, [mainDir '\jpgs\NATIONWIDE_CASES.jpg']);
    
    figure('Visible', 'off'); grid on; hold on;
    plot((1:currDay)+offsetDays, sum(deathsPerDay), ...
        'Color', [0.8500 0.3250 0.0980]);
    plot((1:currDay)+offsetDays, sum(deaths7day), 'LineWidth', 4, ...
        'Color', [0.8500 0.3250 0.0980]);
    xlabel('Days Since 2020 Began');
    ylabel('New COVID-19 Deaths');
    title('COVID-19 Stats for the United States of America');
    legend(['COVID-19 Deaths: ' num2str(deathPerMilNation) ...
        ' per million'], '7 Day Average of Deaths', ...
        'Location', 'northwest');
    set(gcf, 'Visible', 'on');
    savefig([mainDir '\figs\NATIONWIDE_DEATHS.fig']);
    exportgraphics(gcf, [mainDir '\jpgs\NATIONWIDE_DEATHS.jpg']);
    
    %Put deaths in perspective with other disastrous events  
    figure('Visible', 'off'); grid on; hold on;
    plot((1:currDay)+offsetDays, sum(deaths) / benghazi, 'LineWidth', 3);
    plot((1:currDay)+offsetDays, sum(deaths) / titanic, 'LineWidth', 3);
    plot((1:currDay)+offsetDays, sum(deaths) / katrina, 'LineWidth', 3);
    plot((1:currDay)+offsetDays, sum(deaths) / pearl, 'LineWidth', 3);
    plot((1:currDay)+offsetDays, sum(deaths) / nine11, 'LineWidth', 3);
    plot((1:currDay)+offsetDays, sum(deaths) / gettysburg, 'LineWidth', 3);
    plot((1:currDay)+offsetDays, sum(deaths) / lung, 'LineWidth', 3);
    plot((1:currDay)+offsetDays, sum(deaths) / flu, 'LineWidth', 3);
    plot((1:currDay)+offsetDays, sum(deaths) / hiroshima, 'LineWidth', 3);
    numBenghazis = floor(max(sum(deaths) / benghazi));
    numTitanics = round(max(sum(deaths) / titanic));
    numKatrinas = round(max(sum(deaths) / katrina));
    numPearls = round(max(sum(deaths) / pearl));
    numNines = round(max(sum(deaths) / nine11));
    numGettys = round(max(sum(deaths) / gettysburg));
    numLungs = round(max(sum(deaths) / lung));
    numFluYears = round(max(sum(deaths) / flu), 1);
    numHiroshimas = round(max(sum(deaths) / hiroshima), 1);
    
    xlabel('Days Since 2020 Began');
    ylabel('Number of Events');
    ylim([0 round(max(sum(deaths) / titanic), -1)]);
    title('COVID-19 Deaths Comparison for the US so Far');
    legend(['Benghazis: ' num2str(numBenghazis)], ...
        ['Titanics: ' num2str(numTitanics)], ...
        ['Hurricane Katrinas: ' num2str(numKatrinas)], ...
        ['Pearl Harbors: ' num2str(numPearls)], ...
        ['9/11s: ' num2str(numNines)], ...
        ['Gettysburgs: ' num2str(numGettys)], ...
        ['Months of US Lung Cancer Deaths: ' num2str(numLungs)], ...
        ['Years of US Flu Deaths: ' num2str(numFluYears)], ...
        ['Hiroshimas: ' num2str(numHiroshimas)], ...
        'Location', 'northwest');
    set(gcf, 'Visible', 'on');
    savefig([mainDir '\figs\NATIONWIDE_DEATHS_COMPARISON.fig']);
    exportgraphics(gcf, ...
        [mainDir '\jpgs\NATIONWIDE_DEATHS_COMPARISON.jpg']);
end

%Print some stats to the command line for the worst impacted states
indices = find(estimatedPercentInfected >= 15);
for ii = 1:length(indices)
    disp([stateOrder{indices(ii)} ': ' ...
        num2str(estimatedPercentInfected(indices(ii))) ...
        '% Total Population Infected']);
end

