% data la matrice del datacube il codice dovrà eseguire DTFT lungo le tre
% direzioni della matrice
clear
clc

fastTimeIndex = 100;
slowTimeIndex = 10;
spatialIndex = 2;

% Number of bins for the FFT computation, for each dimension. Should be a power of 2.
bins_fastTime = 2^nextpow2(500);
bins_slowTime = 2^nextpow2(500);
bins_spatial = 2^nextpow2(500);

datacube = zeros(fastTimeIndex,slowTimeIndex,spatialIndex);

FFTcube = fftn(datacube,[bins_fastTime bins_slowTime bins_spatial]);