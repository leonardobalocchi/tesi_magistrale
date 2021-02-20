% data la matrice del datacube il codice dovrà eseguire DTFT lungo le tre
% direzioni della matrice
clear
clc

fastTimeIndex = 10
slowTimeIndex = 10
spatialIndex = 10

binsNumber = 1024 % Number of bins for the FFT computation. Should be a power of 2.

datacube = zeros(fastTimeIndex,slowTimeIndex,spatialIndex)

FFTcube_slow = zeros(fastTimeIndex,spatialIndex,binsNumber); % to fill with datacube's FFT through slow time axis
FFTcube_fast = zeros(slowTimeIndex,spatialIndex,binsNumber); % to fill with datacube's FFT through fast time axis
FFTcube_spatial = zeros(fastTimeIndex,slowTimeIndex,binsNumber); % to fill with datacube's FFT through spatial axis

% We don't have experimental/simulated data, so we can't fill the datacube yet.

% FFT of the slow time signals.
for i=1:fastTimeIndex

   for j=1:spatialIndex

       FFTcube_slow(i,j,:) = fft(datacube(i,:,j),binsNumber);

   end

end

% FFT of the fast time signals.
for i=1:slowTimeIndex

   for j=1:spatialIndex

       FFTcube_fast(i,j,:) = fft(datacube(:,i,j),binsNumber);

   end

end

% FFT of the "spatial" signal.
for i=1:fastTimeIndex

   for j=1:slowTimeIndex

       FFTcube_spatial(i,j,:) = fft(datacube(i,j,:),binsNumber);

   end

end
