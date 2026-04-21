function wynik = halley_fraktal(wspolczynniki, rozdzielczosc, maks_iteracji, tolerancja, zakres)
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

% komentarz

% funkcja przyjmuje:
% - wspolczynniki wielomianu (od najwyzszej potegi do wyrazu wolnego)
% - rozdzielczosc siatki
% - maksymalna liczbe iteracji metody halleya
% - tolerancje, czyli kiedy uznajemy ze jestesmy juz dostatecznie blisko zera
% - zakres rysowania na plaszczyznie zespolonej [xmin xmax ymin ymax]

% funkcja zwraca strukture wynik, w ktorej sa:
% - os_rzeczywista i os_urojona, czyli osie do rysowania wykresu
% - iteracje, czyli ile krokow halleya bylo potrzebne w kazdym punkcie
% - mapa_zer, czyli do ktorego zera zbiegl dany punkt startowy
% - zera, czyli znalezione przyblizenia zer wielomianu

% najpierw tworzymy prostokatna siatke punktow na plaszczyznie zespolonej
% linspace robi rowno rozlozone punkty na osi rzeczywistej i urojonej
% meshgrid zamienia to w pelna siatke punktow
% potem skladamy z tego liczby zespolone z = x + iy

% iteracje to macierz zer o tym samym rozmiarze co siatka
% bedziemy tam zapisywac ile iteracji wykonalo sie dla kazdego punktu startowego

% potem idzie glowna petla po kolejnych krokach metody halleya
% nie lecimy punkt po punkcie osobno, tylko robimy wszystko naraz na calej macierzy Z

% w kazdym kroku liczymy:
% - w(Z), czyli wartosc wielomianu
% - wp(Z), czyli pierwsza pochodna
% - wpp(Z), czyli druga pochodna
% za pomoca funkcji horner_pochodne

% dalej liczymy mianownik ze wzoru halleya:
% 2 * wp^2 - w * wpp

% maska to tablica logiczna true/false
% true jest tam, gdzie punkt nadal trzeba iterowac:
% - bo abs(w) > tolerancja, czyli jeszcze nie doszlismy do zera
% - i jednoczesnie mianownik nie jest za maly, zeby nie dzielic przez prawie 0

% if ~any(maska, 'all') znaczy:
% jesli nigdzie nie ma juz punktu do poprawiania, to przerywamy petle
% ~ to negacja, czyli "nie"

% potem wykonujemy wlasciwy krok metody halleya:
% z_{k+1} = z_k - (2 * w * wp) / (2 * wp^2 - w * wpp)
% ale tylko dla tych punktow, gdzie maska = true

% przy okazji zwiekszamy licznik iteracji dla tych punktow o 1

% po zakonczeniu petli jeszcze raz liczymy w w punktach koncowych
% i sprawdzamy ktore punkty faktycznie spelniaja warunek abs(w) <= tolerancja

% potem bierzemy tylko te punkty zbiezne i probujemy z nich odczytac same zera
% poniewaz numerycznie punkty bardzo bliskie temu samemu zeru nie sa idealnie rowne,
% to zaokraglamy je do 3 miejsc po przecinku
% potem unique usuwa powtorzenia i zostaja rozne znalezione zera

% mapa_zer to macierz, w ktorej wpisujemy numer zera dla kazdego punktu siatki
% dzieki temu mozna potem narysowac obszary przyciagania zer

% na koncu wszystko pakujemy do struktury wynik, zeby latwo bylo to potem wykorzystac
% w skrypcie glownym do wyswietlenia zer i rysowania wykresow