%% Script to call function script to calculate sumax 
%##############################################################################################################
% Script and variables need to be adapted in sections: 'load data'; 'assign variable names'; 'configure run specifics'

% required input data:
% - P & date: daily precipitation with corresponding date vector (time series, time in rows, average per outlet in columns)
% - Ep: long term monthly averaged daily potential evaporation  (time series, same length as precipitation series and with corresponding dates, same units as precipitation)
% - Q: long term averaged daily discharge (per outlet averaged or actual values in row vector, same units and length as precipitation series)
% - imax: interception capacity (row vector with averages per outlet, same units as precipitation)
% - Kcap: capilary rise (same time unit as preciptation (eg. probably "per day", row vector with values per outlet, zeros is no capilary rise)
% - ndvi: ndvi/evi data, only for some options (long term monthly averages, same length as precipitation series)
% - dorm: summer and winter dormancy months (12 values per outlet, following

% Written by Tanja Euser (t.euser@tudelft.nl), May 2015 
% Last adapted by Tanja Euser, May 2015
%##############################################################################################################

%% load data
% This section needs to be written specifically for each case. 


%% Assign variable names
% In this section all variables are filled, all values can be case specific
% and therefore needs to be checked for correctness!
date = dateVar;
P = precipVar; 
Ep = evapVar;
Q = Qvar;
imax = [1.5 1.5 1.5]; 
ndvi = ndviVar;
dorm = [0 0 0];       %months in which dormancy of vegetation occurs (summer or winter)
Kcap = [0 0 0];
Treturn = [1.1 2 3 4 5 10 15 20 30 40];

%% configure run specifics
% In this section the specifications of the sumax calculation are defined
monthS = 4;                     %month to start the analysis (best to choose a month in which catchment is saturated
smthSpan = 1;                   %if 1, than precipitation is not smoothed,otherwise it is smoothed with this window
EaMethod = 'flatLine';          % method to apply seasonal cycle to Ea (options: flatLine/gradual/gradual_ndvi)
extremeValueDist = 'Gumbel';    % extreme value distributed used to calculate sumax for different return periods (options: Gumbel/Weibul)
dormFormat = 'month';           % method to apply dormancy to Ea (options: precipForm/month)
correctRC_P = 0;                % correction of precip data in case runoff coefficient is larger than threshold   
correctRC_Q = 0;                % correction of discharge data in case runoff coefficient is larger than threshold   
thresRC = 0.95;                 % threshold for runoff coefficient, only relevant in case of correction
folder = 'D:\TEuser\Onderzoek\Promotie\data processing\generic functions\sumax_calculation\';                 % folder to save output data
fileName = ['sumaxCalculation_',datestr(date(1)),'-',datestr(date(end)),'.mat'];                              % file name for output  

%% check runoff coefficient (P > Q)
%In this section it is checked whether the precipitation is larger than the
%discharge. If wished, the precipitation or discharge can be corrected in
%order to get a runoff coefficient smaller than the specified threshold
Ptot = zeros(length(P(1,:)),1);
Eptot = zeros(length(P(1,:)),1);
Qtot = zeros(length(P(1,:)),1);
if or(correctRC_P == 1, correctRC_Q == 1)
    for i = 1:length(P(1,:))
        ind = find(P(:,i) >= 0 & Q(:,i) >= 0); 
        Ptot(i) = sum(P(ind,i));
        Qtot(i) = sum(Q(ind,i));
        Eptot(i) = sum(Ep(ind,i));
    
        rcT = Qtot(i) / Ptot(i); 
        if rcT > thresRC
            % correction of precipitation data
            if correctRC_P == 1
                P(:,i) = P(:,i) * rcT / thresRC;
                Ptot(i) = sum(P(ind,i));
            % correction of discharge data
            elseif correctRC_Q == 1
                Q(:,i) = Q(:,i) * thresRC / rcT;
                Qtot(i) = sum(Q(ind,i));
            end
        end
    end
    RC = Qtot ./ Ptot;
else
    RC = Qtot ./ Ptot;
end

%% calculating sumax values - on continuous basis
% In this section the script is called to calculate sumax and the values
% are combined into 1 matrix 
sumaxYear = zeros(floor(length(P(:,1))/365),length(P(1,:)));
soilDefit = zeros(size(P));
Psmth = zeros(size(P));
Ea = zeros(size(P));
EpAve = zeros(size(P));
for i = 1:length(P(1,:))
    disp(['station ',num2str(i)])
    [Temp1,Temp2,Temp3,Temp4] = calcSumax_endlessReservoir(date,P(:,i),Q(:,i),Ep(:,i),imax(i),smthSpan,monthS,EaMethod,dorm,dormFormat,ndvi(:,i),Kcap(i),extremeValueDist);
    sumaxYear(:,i) = Temp1;                 % maximum required sumax per year
    Psmth(:,i) = Temp2;                     % smoothed precipitation
    Ea(:,i) = Temp3;                        % estimated actual evaporation
    soilDefit(:,i) = Temp4(2:end);          % simulated soil deficit   
end

%% calculate sumax correspondig to drought for specific return periods
% In this section the required storage (sumax) for droughts with specified
% return periods are calculated

%Gumbel
if strcmp(extremeValueDist,'Gumbel')
sumax = zeros(length(sumaxYear(1,:)),length(Treturn));
for i = 1:length(sumaxYear(1,:))
    for n = 1:length(Treturn)
        ind = find(isfinite(sumaxYear(:,i)));
        sumax(i,n) = calcGumbel(sumaxYear(ind,i),Treturn(n));
    end
end
%Weibull
elseif strcmp(extremeValueDist,'Weibul')
    disp('This option is not programmed yet')
end
    
%% write data to file
% In this section the calculated values are written to a .mat file,
% together with the input data and run specifics
f0 = 'date'; v0 = date;
f1 = 'precipitation'; v1 = P;
f2 = 'potentialEvaporation'; v2 = Ep;
f3 = 'discharge'; v3 = Q;
f4 = 'interceptionStorage'; v4 = imax;
f5 = 'NDVI'; v5 = ndvi;
f6 = 'smoothingSpan'; v6 = smthSpan;
f7 = 'startMonthAnalysis';  v7 = monthS;
f8 = 'periodsWithVegetationDormancy';  v8 = dorm;
f9 = 'capillaryRise';  v9 = Kcap;
f10 = 'returnPeriodsForExtremeValuesDistribution';  v10 = Treturn;
f11 = 'EaMethodForSeasonality';  v11 = EaMethod;
f12 = 'extremeValueDistribution';  v12 = extremeValueDist;
f13 = 'correctionAppliedForRC_P';  v13 = correctRC_P;
f14 = 'correctAppliedForRC_Q';  v14 = correctRC_Q; 
f15 = 'thresholdAppliedForRCcorrection';  v15 = thresRC;
runSpecifics = struct(f0,v0,f1,v1,f2,v2,f3,v3,f4,v4,f5,v5,f6,v6,f7,v7,f8,v8,f9,v9,f10,v10,f11,v11,f12,v12,f13,v13,f14,v14,f15,v15);

save([folder,fileName],'runSpecifics','sumax','RC','sumaxYear','soilDefit')



