% funtion to calculate Gumbel distribution for a specific gauge and
% returning extremes with different return periods

function[Xreturn] = calcGumbel(X,Treturn)

% X = (unsorted) data series of maxima
% Treturn = vector with required return periods
% Xreturn = vector with values corresponding to Treturn

Xsort = sort(X,'descend');

i = 1:1:length(Xsort);
T = (length(Xsort) + 1) ./ i;
oneT = ones(size(T));
y = -log(-log(oneT - oneT ./ T));

z = polyfit(y,Xsort',1);

Xreturn = zeros(length(Treturn),1);
for i = 1:length(Treturn)
    yreturn = -log(-log(1 - 1 / Treturn(i)));
    Xreturn(i) = z(1) * yreturn + z(2);
end
