% function to apply a seasonal signal to the longterm average evaporation,
% based on a constant actual evaporation during summer

function[Ea] = getEaSeasonal_flatLine(EpAve,EaAveC)

%Ep = daily/monthly potential evaporation series
%EaAveC = longterm average actual evaporation for same timestep as Ep
%date = datenumbers corresponding with Ep data
%Ea = actual evaporation time series with seasonal pattern

EaAve = ones(size(EpAve)) * EaAveC;            %get timeseries with constant average Ea

a = min(EpAve,EaAve);
EaMiss = EaAve - a;
EpMiss = EpAve - a;
x = EaAveC / 8;
i = 2;
err = max(EpAve);
while abs(err) > 0.01 * x
    x(i) = x(i-1) + err/(length(EpAve)*0.8);
    addEa = min(EpMiss, x(i));
    A = sum(EaMiss);
    B = sum(addEa);
    err = A - B;
    
    i = i+1;
    disp(num2str(i-1))
    if EpMiss == addEa
        break
    end
end

try
    Ea = min(EpAve, EaAve + addEa);
catch
    disp('addEa was not defined')
    Ea = zeros(size(EaAve)) * NaN;
end
