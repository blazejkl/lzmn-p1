wspolczynniki = [1 -1 1 -1];

wynik = halley_fraktal(wspolczynniki, 2000, 40, 1e-6, [-2 2 -2 2]);

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