wspolczynniki = [12, -3, 5, -1];

wynik = halley_fraktal(wspolczynniki, 800, 40, 1e-6, [-2 2 -2 2]);

disp('Znalezione zera:')
disp(wynik.zera.')

figure;
imagesc(wynik.os_rzeczywista, wynik.os_urojona, wynik.iteracje);
axis xy equal tight;
colorbar;
title('Liczba iteracji metody Halleya');
xlabel('Re(z)');
ylabel('Im(z)');

figure;
imagesc(wynik.os_rzeczywista, wynik.os_urojona, wynik.mapa_zer);
axis xy equal tight;
colorbar;
title('Obszary przyciagania zer');
xlabel('Re(z)');
ylabel('Im(z)');

% komentarz

% tutaj podajemy wspolczynniki wielomianu i wywolujemy glowna funkcje
% wynik zawiera znalezione zera, liczbe iteracji oraz mape obszarow przyciagania

% najpierw wypisujemy znalezione zera w command window

% potem rysujemy 2 wykresy:
% - liczbe iteracji metody halleya w kazdym punkcie startowym
% - obszary przyciagania zer, czyli do ktorego zera zbiega dany punkt

% imagesc sluzy do pokazania tych macierzy jako obraz