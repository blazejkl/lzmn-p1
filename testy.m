a = [2 -3 0 1];

disp('Test 1: horner - punkt rzeczywisty')
z = 2;
[w, wp, wpp] = horner_pochodne(a, z);
assert(abs(w - polyval(a, z)) < 1e-12);
assert(abs(wp - polyval(polyder(a), z)) < 1e-12);
assert(abs(wpp - polyval(polyder(polyder(a)), z)) < 1e-12);

disp('Test 2: horner - punkt zespolony')
z = 1 + 2i;
[w, wp, wpp] = horner_pochodne(a, z);
assert(abs(w - polyval(a, z)) < 1e-12);
assert(abs(wp - polyval(polyder(a), z)) < 1e-12);
assert(abs(wpp - polyval(polyder(polyder(a)), z)) < 1e-12);

disp('Test 3: horner - wektor punktow')
z = [-2 -1 0 1 2];
[w, wp, wpp] = horner_pochodne(a, z);
assert(max(abs(w - polyval(a, z))) < 1e-12);
assert(max(abs(wp - polyval(polyder(a), z))) < 1e-12);
assert(max(abs(wpp - polyval(polyder(polyder(a)), z))) < 1e-12);

disp('Test 4: horner - wektor zespolony')
z = [-1-1i, -1+i, 1-1i, 1+i];
[w, wp, wpp] = horner_pochodne(a, z);
assert(max(abs(w - polyval(a, z))) < 1e-12);
assert(max(abs(wp - polyval(polyder(a), z))) < 1e-12);
assert(max(abs(wpp - polyval(polyder(polyder(a)), z))) < 1e-12);

disp('Test 5: horner - wielomian staly')
a_stale = [5];
z = [0 1 2 3+4i];
[w, wp, wpp] = horner_pochodne(a_stale, z);
assert(all(w == 5));
assert(all(wp == 0));
assert(all(wpp == 0));

disp('Test 6: halley - z^4 - 1 daje 4 zera')
a1 = [1 0 0 0 -1];
wynik1 = halley_fraktal(a1, 300, 40, 1e-6, [-2 2 -2 2]);
zera1 = roots(a1);
assert(length(wynik1.zera) == 4);

disp('Test 7: halley - znalezione zera sa blisko roots dla z^4 - 1')
for k = 1:length(zera1)
    assert(min(abs(wynik1.zera - zera1(k))) < 1e-2);
end

disp('Test 8: halley - residuale male dla z^4 - 1')
for k = 1:length(wynik1.zera)
    [w, ~, ~] = horner_pochodne(a1, wynik1.zera(k));
    assert(abs(w) < 1e-4);
end

disp('Test 9: halley - z^3 - 1 daje 3 zera')
a2 = [1 0 0 -1];
wynik2 = halley_fraktal(a2, 300, 40, 1e-6, [-2 2 -2 2]);
zera2 = roots(a2);
assert(length(wynik2.zera) == 3);

disp('Test 10: halley - znalezione zera sa blisko roots dla z^3 - 1')
for k = 1:length(zera2)
    assert(min(abs(wynik2.zera - zera2(k))) < 1e-2);
end

disp('Test 11: halley - z^2 + 1 daje 2 zera')
a3 = [1 0 1];
wynik3 = halley_fraktal(a3, 300, 40, 1e-6, [-2 2 -2 2]);
zera3 = roots(a3);
assert(length(wynik3.zera) == 2);

disp('Test 12: halley - znalezione zera sa blisko roots dla z^2 + 1')
for k = 1:length(zera3)
    assert(min(abs(wynik3.zera - zera3(k))) < 1e-2);
end

disp('Wszystkie testy zaliczone.')

% komentarz

% tutaj podajemy wspolczynniki wielomianu i wywolujemy glowna funkcje
% wynik zawiera znalezione zera, liczbe iteracji oraz mape obszarow przyciagania

% najpierw wypisujemy znalezione zera w command window

% potem rysujemy 2 wykresy:
% - liczbe iteracji metody halleya w kazdym punkcie startowym
% - obszary przyciagania zer, czyli do ktorego zera zbiega dany punkt

% imagesc sluzy do pokazania tych macierzy jako obraz