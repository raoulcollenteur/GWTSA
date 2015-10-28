% function to calculate sumax values based on precipitation deficits

function[sumax_Ymax, Psmth, Ea, soilDefit] = calcSumax_endlessReservoir(date,P,Q,Ep,imax,smthSpan,startM,EaMethod,dorm,dormFormat,ndvi,Kcap,extremeValueDist)
% date = datenumbers corresponding to P and Q data; column vector
% P = precipitation [L/T]; column vector
% Q = observed discharge [L/T]; column vector
% Ep = calculated potential evaporation [L/T]; column vector
% imax = interception evaporation [L]
% smthSpan = number of timestep to smooth over [T]
% startM = starting month --> catchment being saturated
% EaMethod = method to calculate Ea estimate (flatLine/gradual/ndvi)
% dorm = indication for dormancy of vegetation: 0/1 in precip format, or months in which dormancy occurs
% dormFormat = format of dormancy vector (precipFormat/month)
% ndvi = NDVI [-]; column vector
% Kcap = coefficient for capillary rise [L/T]

%% preparing input data
% convert date to datevector
if length(date(1,:)) == 1
    dateV = datevec(date);
else
    dateV = data;
end

% substracting interception from precipitation and potential evaporation via an interception reservoir 
Icep = zeros(length(P(:,1)),1); Si = 0; PminI = zeros(length(P(:,1)),1); EpminI = zeros(length(P(:,1)),1);
for i = 1:length(P(:,1))
    Si = Si + P(i);
    PminI(i) = max(Si-imax,0); Si = Si - PminI(i);      % effective rainfall        
    Icep(i) = min(Si,Ep(i)); Si = Si - Icep(i);         % interception evaporation   
    EpminI(i) = Ep(i) - Icep(i);                        % potential transpiration
end
ind = find(P >= 0 & Q >= 0);
PminI_filt = PminI(ind);                               %including interception excluding years with missing data in P/Q, for calculation RC   
EpminI_filt = EpminI(ind);
PminI(isnan(P)) = NaN;

% estimating actual evaporation
Ptot = sum(PminI_filt); 
Qtot = sum(Q(ind));
Eptot = sum(EpminI_filt);
EaTot = min(Ptot - Qtot, Eptot);
EaAve = EaTot / length(PminI_filt);
if EaAve < 0
    disp('warning: Ea negative')
    EaAve = 0;
end
if strcmp(EaMethod,'gradual')
    Ea = getEaSeasonal_gradual(EpminI,EaAve);
elseif strcmp(EaMethod,'flatLine')
    Ea = getEaSeasonal_flatLine(EpminI,EaAve);
elseif strcmp(EaMethod,'gradual_NDVI')
    Ea = getEaSeasonal_gradual_NDVI(EpminI,EaAve,ndvi);
end

%account for dormancy
if strcmp(dormFormat,'month')
    Ea(ismember(dateV(:,2),dorm)) = 0;
elseif strcmp(dormFormat,'precipFormat')
    Ea(dorm == 1) = 0;
end

%smoothing rainfall data
Psmth = smooth(PminI,smthSpan);

%% keeping track of soil moisture deficit
soilDefit = zeros(length(Psmth)+1,1);
for i = 1:length(Psmth)
    soilDefit(i+1) = min(soilDefit(i) + Psmth(i) - Ea(i) + (soilDefit(i) * Kcap * -1), 0);           %calculate cummulative of sum P - Ea with maximum of zero, so only accounting "negative storage"
end
soilDefit(isnan(Psmth)) = NaN;

%% calculating max sumax per hydrological year
if strcmp(extremeValueDist, 'Gumbel')
indD = find(dateV(:, 2) == startM & dateV(:, 3) == 1);               %find all starts of the hyd year
sumax_Ymax = zeros(length(indD),1);                                  %set length of sumax vector   

for n = 1:length(indD)
    try
       X = min(soilDefit(indD(n):indD(n+1)-1,1)) ;
       sumax_Ymax(n) = X * -1;                  
       if  soilDefit(indD(n+1)) ~= 0
            indO1 = find(soilDefit(indD(n):indD(n+1)-1,1) == 0, 1, 'last');
            indO2 = find(soilDefit(indD(n+1):indD(n+2)-1,1) == 0, 1, 'first');
            if length(soilDefit(indD(n):indD(n+1)-1,1)) - indO1 > indO2             %dry period mainly occurs in yr 1
                indD(n+1) = indD(n+1) + indO2;
            elseif length(soilDefit(indD(n):indD(n+1)-1,1)) - indO1 < indO2             %dry period mainly occurs in yr 2
                indD(n+1) = indD(n) + indO1;
                X = min(soilDefit(indD(n):indD(n+1)-1,1)) ;
                sumax_Ymax(n) = X * -1;
            end    
       end
    catch
       disp(['year ', num2str(n), ' has been skipped']) 
       sumax_Ymax(n) = NaN;
    end
end

elseif strcmp(extremeValueDist, 'Weibull')
    disp('This is not programmed yet')
end