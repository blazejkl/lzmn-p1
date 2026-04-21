function [w, wp, wpp] = horner_pochodne(wspolczynniki, z) 
    w = wspolczynniki(1) * ones(size(z));
    wp = zeros(size(z));
    wpp = zeros(size(z));
    
    for i = 2:length(wspolczynniki)
        wpp = wpp .* z + 2 * wp;
        wp = wp .* z + w;
        w = w .* z + wspolczynniki(i);
    end
end

% komentarz

% funkcja przyjmuje wektor współczynników wielomianu (w kolejności od
% najwyższej potęgi do wyrazu wolnego) oraz konkretne punkty w jego
% dziedzinie i zwraca wartości wielomianu i jego dwóch pierwszych
% pochodnych w tych punktach wykorzystując schemat Hornera

% w pierwszym etapie, przed pętlą, do w przypisujemy współczynnik przy
% najwyższej potędze, do wp i wpp 0 (dla każdego punktu)

% no i potem iterujemy po pozostałych współczynnikach
% dla w: mnożymy w przez z i dodajemy kolejny współczynnik - wprost wzór
% schematu Hornera
% dla wp: mnożymy wp przez z i dodajemy w (pochodna ze współczynnika -
% stałej jest równa zero, zostaje więc pochodna iloczynu w i z - wp * z +
% w)
% dla wpp: mnożymy wpp przez z i dodajemy 2 * wp (różniczkujemy wp i znowu 
% z pochodnej iloczynu i sumy otrzymujemy wzó®tarz do
% khalleyatarzdohalleya