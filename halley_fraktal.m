function wynik = halley_fraktal(wspolczynniki, rozdzielczosc, maks_iteracji, tolerancja, zakres)
% Funkcja wyznacza zera wielomianu metoda Halleya dla siatki punktow.
% Wejscie:
%   wspolczynniki - wektor wspolczynnikow wielomianu
%   rozdzielczosc - liczba punktow siatki w jednym wymiarze (N x N)
%   maks_iteracji - limit iteracji algorytmu
%   tolerancja - kryterium stopu dla wartosci wielomianu (reziduum)
%   zakres - wektor graniczny siatki [Re_min, Re_max, Im_min, Im_max]
% Wyjscie:
%   wynik - struktura z danymi siatki, iteracjami oraz znalezionymi zerami
    os_rzeczywista = linspace(zakres(1), zakres(2), rozdzielczosc);
    os_urojona = linspace(zakres(3), zakres(4), rozdzielczosc);
    [X, Y] = meshgrid(os_rzeczywista, os_urojona);
    Z = X + 1i * Y;
    iteracje = zeros(size(Z));

    for k = 1:maks_iteracji
        [w, wp, wpp] = horner_pochodne(wspolczynniki, Z);
        mianownik = 2 .* wp.^2 - w .* wpp;
        maska = abs(w) > tolerancja & abs(mianownik) > 1e-12;

        if ~any(maska, 'all')
            break;
        end

        Z(maska) = Z(maska) - (2 .* w(maska) .* wp(maska)) ./ mianownik(maska);
        iteracje(maska) = iteracje(maska) + 1;
    end

    [w, ~, ~] = horner_pochodne(wspolczynniki, Z);
    maska_zbiezne = abs(w) <= tolerancja;

    Z_koncowe = Z(maska_zbiezne);
    Z_zaokraglone = round(1000 * real(Z_koncowe)) / 1000 + 1i * round(1000 * imag(Z_koncowe)) / 1000;
    zera = unique(Z_zaokraglone).';

    mapa_zer = zeros(size(Z));
    for k = 1:length(zera)
        mapa_zer(maska_zbiezne & abs(Z - zera(k)) < 1e-3) = k;
    end

    wynik.os_rzeczywista = os_rzeczywista;
    wynik.os_urojona = os_urojona;
    wynik.iteracje = iteracje;
    wynik.mapa_zer = mapa_zer;
    wynik.zera = zera;
end
