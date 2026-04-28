function [w, wp, wpp] = horner_pochodne_petla(wspolczynniki, z)
    % Niezwektoryzowana wersja schematu Hornera.
    % Oblicza wartosci element po elemencie za pomoca petli.
    w = zeros(size(z));
    wp = zeros(size(z));
    wpp = zeros(size(z));
    
    for k = 1:numel(z)
        akt_z = z(k);
        akt_w = wspolczynniki(1);
        akt_wp = 0;
        akt_wpp = 0;
        
        for i = 2:length(wspolczynniki)
            akt_wpp = akt_wpp * akt_z + 2 * akt_wp;
            akt_wp = akt_wp * akt_z + akt_w;
            akt_w = akt_w * akt_z + wspolczynniki(i);
        end
        
        w(k) = akt_w;
        wp(k) = akt_wp;
        wpp(k) = akt_wpp;
    end
end