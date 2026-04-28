a = [2 -3 0 1];

disp('Test 1: horner - punkt rzeczywisty')
z = 2;
[w, wp, wpp] = horner_pochodne(a, z);
assert(abs(w - polyval(a, z)) < 1e-12);

disp('Test 2: horner - wektor zespolony')
z = [-1-1i, -1+i, 1-1i, 1+i];
[w, wp, wpp] = horner_pochodne(a, z);
assert(max(abs(w - polyval(a, z))) < 1e-12);

disp('Test 3: horner - wielomian staly')
a_stale = [5];
z = [0 1 2 3+4i];
[w, wp, wpp] = horner_pochodne(a_stale, z);
assert(all(w == 5) && all(wp == 0) && all(wpp == 0));

a1 = [1 0 0 0 -1];
disp('Test 4: halley - z^4 - 1 daje 4 zera')
wynik1 = halley_fraktal(a1, 100, 40, 1e-6, [-2 2 -2 2]);
assert(length(wynik1.zera) == 4);

disp('Test 5: halley - residuale male dla z^4 - 1')
for k = 1:length(wynik1.zera)
    [w, ~, ~] = horner_pochodne(a1, wynik1.zera(k));
    assert(abs(w) < 1e-4);
end

disp('Test 6: halley - z^3 - 1 daje 3 zera')
a2 = [1 0 0 -1];
wynik2 = halley_fraktal(a2, 100, 40, 1e-6, [-2 2 -2 2]);
assert(length(wynik2.zera) == 3);

% testy parametryczne

disp('Test 7: Wplyw tolerancji na liczbe iteracji (1e-4 vs 1e-14)')
wynik_tol_slaba = halley_fraktal(a1, 100, 40, 1e-4, [-2 2 -2 2]);
wynik_tol_ostra = halley_fraktal(a1, 100, 40, 1e-14, [-2 2 -2 2]);
sr_iter_slaba = mean(wynik_tol_slaba.iteracje, 'all');
sr_iter_ostra = mean(wynik_tol_ostra.iteracje, 'all');
assert(sr_iter_ostra >= sr_iter_slaba);
assert(length(wynik_tol_slaba.zera) == 4 && length(wynik_tol_ostra.zera) == 4);

disp('Test 8: Wplyw rozdzielczosci na wymiary siatki wyników')
wynik_res_mala = halley_fraktal(a1, 50, 40, 1e-6, [-2 2 -2 2]);
wynik_res_duza = halley_fraktal(a1, 150, 40, 1e-6, [-2 2 -2 2]);
assert(size(wynik_res_mala.iteracje, 1) == 50);
assert(size(wynik_res_duza.iteracje, 1) == 150);

disp('Test 9: Zbieznosc liniowa dla pierwiastkow wielokrotnych')
a_pojedynczy = [1 0 -1]; % z^2 - 1 (zera pojedyncze)
a_wielokrotny = [1 -2 1]; % (z-1)^2 (zero podwojne)

wynik_poj = halley_fraktal(a_pojedynczy, 50, 40, 1e-6, [-2 2 -2 2]);
wynik_wiel = halley_fraktal(a_wielokrotny, 50, 40, 1e-6, [-2 2 -2 2]);

assert(mean(wynik_wiel.iteracje, 'all') > mean(wynik_poj.iteracje, 'all'));

disp('Wszystkie testy zaliczone.')

% tabelaryczne zestawienie wyników i metryk
disp(' ');
disp('--- TABELARYCZNE ZESTAWIENIE WYNIKÓW (Halley vs roots) ---');

wielomiany = {
    [1 0 0 0 -1], ...       % 1. Standardowy (z^4 - 1)
    [1 0 1], ...            % 2. Zera tylko urojone (z^2 + 1)
    [1 0 -2 2], ...         % 3. Zera mieszane (z^3 - 2z + 2)
    [1 -2 1], ...           % 4. Pierwiastek wielokrotny (z-1)^2
    [1000 500 -200 100], ...% 5. Duze wspolczynniki
    [5]                     % 6. Brak rozwiazan (wielomian staly w(z)=5)
};

nazwy = {
    'z^4 - 1 (zera proste)',
    'z^2 + 1 (zera urojone)',
    'z^3 - 2z + 2 (zera mieszane)',
    '(z - 1)^2 (zero wielokrotne, spadek zbieznosci)',
    '1000z^3 + 500z^2 - 200z + 100 (duze wspolczynniki)',
    'w(z) = 5 (brak zer - przypadek bez wyniku)'
};

for i = 1:length(wielomiany)
    a_tab = wielomiany{i};
    fprintf('\nPrzyklad %d: %s\n', i, nazwy{i});
    
    wynik_tab = halley_fraktal(a_tab, 100, 40, 1e-10, [-2 2 -2 2]);
    zera_ref = roots(a_tab);
    
    % Sprawdzenie przypadku nr 6 (gdzie algorytm nie zwraca wyniku)
    if isempty(wynik_tab.zera)
        disp(' -> Algorytm nie zwrocil zadnego wyniku. Uzasadnienie:');
        disp(' -> Funkcja stala ma pochodne rowne 0. Mianownik we wzorze Halleya');
        disp(' -> wynosi 0, wiec bezpiecznik w kodzie przerywa iteracje.');
        continue;
    end

    % Rysowanie tabelki dla przykladow, ktore maja zera
    fprintf('%-6s | %-22s | %-12s | %-12s\n', 'Zero', 'Wynik Halley', 'Blad Abs.', 'Reziduum |w|');
    fprintf(repmat('-', 1, 64)); fprintf('\n');

    for k = 1:length(wynik_tab.zera)
        z_h = wynik_tab.zera(k);
        
        if isempty(zera_ref)
            blad = NaN;
        else
            [blad, ~] = min(abs(z_h - zera_ref)); 
        end
        
        [w_val, ~, ~] = horner_pochodne(a_tab, z_h);
        reziduum = abs(w_val); 
        
        fprintf('Z_%-2d  | %9.5f %+.5fi | %12.2e | %12.2e\n', k, real(z_h), imag(z_h), blad, reziduum);
    end
end
disp(' ');

% graficzna analiza parametrów wejściowych

a_wykres = [1 0 0 0 -1];
tolerancje = [1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14];
srednie_iteracje = zeros(size(tolerancje));

for i = 1:length(tolerancje)
    w_temp = halley_fraktal(a_wykres, 50, 40, tolerancje(i), [-2 2 -2 2]);
    srednie_iteracje(i) = mean(w_temp.iteracje, 'all');
end

figure('Name', 'Analiza parametrow wejsciowych - Metoda Halleya');
semilogx(tolerancje, srednie_iteracje, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
set(gca, 'XDir', 'reverse');
grid on;
title('Wplyw tolerancji (\epsilon) na sredni koszt obliczeniowy');
xlabel('Tolerancja bledu (skala logarytmiczna)');
ylabel('Srednia liczba iteracji algorytmu na siatce');
