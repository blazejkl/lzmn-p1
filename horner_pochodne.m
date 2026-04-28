function [w, wp, wpp] = horner_pochodne(wspolczynniki, z) 
% Funkcja oblicza wartosc wielomianu oraz jego pierwszej i drugiej pochodnej
% za pomoca zwektoryzowanego algorytmu Hornera.
% Wejscie:
%   wspolczynniki - wektor wspolczynnikow wielomianu (od najwyzszej potegi)
%   z - punkt lub macierz punktow ewaluacji
% Wyjscie:
%   w, wp, wpp - wartosc funkcji, pierwszej i drugiej pochodnej
    w = wspolczynniki(1) * ones(size(z));
    wp = zeros(size(z));
    wpp = zeros(size(z));
    
    for i = 2:length(wspolczynniki)
        wpp = wpp .* z + 2 * wp;
        wp = wp .* z + w;
        w = w .* z + wspolczynniki(i);
    end
end
