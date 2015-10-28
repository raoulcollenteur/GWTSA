% function to apply a seasonal signal to the longterm average evaporation,
% based on a gradual increase of actual evaporation during summer

function[Ea] = getEaSeasonal_gradual(EpAve,EaAveC)

%Ep = daily/monthly potential evaporation series
%EaAveC = longterm average actual evaporation for same timestep as Ep
%date = datenumbers corresponding with Ep data
%Ea = actual evaporation time series with seasonal pattern

EaAve = ones(size(EpAve)) * EaAveC;            %get timeseries with constant average Ea

a = min(EpAve,EaAve);
A = sum(EaAve - a);
B = sum(EpAve - a);
Ea = A ./ B .* (EpAve - a) + a;