% function to apply a seasonal signal to the longterm average evaporation,
% based on a constant actual evaporation during summer and corrected with
% NDVI or EVI data. Estimates of actual evaporation are corrected with
% potential evaporation

function [Ea] = getEaSeasonal_gradual_ndvi(EpAve,EaAveC,ndvi)

%Ep = monthly potential evaporation series
%EaAveC = longterm average actual evaporation for same timestep as Ep
%ndvi = ndvi/evi longterm monthly average values

EaAve = ones(size(EpAve)) * EaAveC;            %get timeseries with constant average Ea

EaTemp = EaAve .* ndvi ./ mean(ndvi);
a = min(EpAve,EaTemp);
A = sum(EaTemp - a);
c = ndvi;
c((EpAve - a) == 0) = 0;
C = sum(c);
Ea = A ./ C .* (c - a) + a;
