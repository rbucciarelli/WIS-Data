% guide to directions

%% defining directions

%%% "CARTESIAN" - 0/360 at east, 90 at north, 180 at west, 270 at south
cart = [1:360];

%%% "POLAR"  - 0 at east, 90 at north, 180/-180 at west, -90 at south
pola = [1:180 -179:0];

%%% "OCEANOGRAPHIC" - 90 at east, 0/360 at north, 270 at west, 180 at south
%%% angles indicate direction FROM
ocea = [89:-1:1 360:-1:90];

%%% "METEOROLOGICAL" - 90 at east, 0/360 at north, 270 at west, 180 at south
%%% angles indicate direction TOWARDS
mete = [89:-1:1 360:-1:90];
% ocdir = metdir-180;
% ocdir(ocdir<0) = ocdir(ocdir<0) + 360;

figure(1); clf; hold on 
plot(cart, '.', 'DisplayName', 'cartesian')
plot(pola, '.', 'DisplayName', 'polar')
plot(ocea, '.', 'DisplayName', 'oceanographic')
% plot(mete, '.', 'DisplayName', 'meteorological')
yticks([-180:45:360]);
grid on; box on;
xticks([0:90:360])
set(gca, 'XTickLabel', {'E', 'N', 'W', 'S'})
legend;

%% conversions

%%% CARTESIAN TO POLAR
cart2pola = cart;
cart2pola(cart2pola>180) = cart2pola(cart2pola>180)-360;
% figure(2); clf; hold on; plot(cart2pol, '.'); plot(pola, '.')
% yticks([-180:45:360]); grid on; box on; xticks([0:90:360]); set(gca, 'XTickLabel', {'E', 'N', 'W', 'S'})

%%% POLAR TO CARTESIAN
pola2cart = pola;
pola2cart(pola<0) = pola2cart(pola<0)+360;
% figure(2); clf; hold on; plot(pola2cart, '.'); plot(cart, '.')
% yticks([-180:45:360]); grid on; box on; xticks([0:90:360]); set(gca, 'XTickLabel', {'E', 'N', 'W', 'S'})


%%% CARTESIAN TO OCEANOGRAPHIC
cart2ocea = cart;
cart2ocea = 450-cart2ocea;
cart2ocea(cart2ocea>360) = cart2ocea(cart2ocea>360)-360;
% figure(2); clf; hold on; plot(cart2ocea, '.'); plot(ocea, '.')
% yticks([-180:45:360]); grid on; box on; xticks([0:90:360]); set(gca, 'XTickLabel', {'E', 'N', 'W', 'S'})

%%% OCEANOGRAPHIC TO CARTESIAN
ocea2cart = ocea;
ocea2cart = 450-ocea2cart;
ocea2cart(ocea2cart>360) = ocea2cart(ocea2cart>360)-360;
% figure(2); clf; hold on; plot(ocea2cart, '.'); plot(cart, '.')
% yticks([-180:45:360]); grid on; box on; xticks([0:90:360]); set(gca, 'XTickLabel', {'E', 'N', 'W', 'S'})


%%% METEOROLOGICAL TO OCEANOGRAPHIC (and vis versa?)
mete2ocea = mete;
mete2ocea = mete2ocea-180;
mete2ocea(mete2ocea<0) = mete2ocea(mete2ocea<0)+360;
% figure(2); clf; hold on; plot(mete2ocea, '.'); plot(ocea, '.')
% yticks([-180:45:360]); grid on; box on; xticks([0:90:360]); set(gca, 'XTickLabel', {'E', 'N', 'W', 'S'})


