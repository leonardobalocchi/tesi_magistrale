% data la matrice del datacube il codice dovrà eseguire DTFT lungo le tre
% direzioni della matrice
clear
clc
datacube = rand(10, 10, 10);

[N, K, M] = size(datacube);
% N fast time index
% K slow time index
% M spatial index

datacube_slow = zeros(N, K, M); % to fill with datacube's FFT through K axis
datacube_fast = zeros(N, K, M); % to fill with datacube's FFT through N axis
datacube_spatial = zeros(N, K, M); % to fill with datacube's FFt throgh M axis

for i=1:M
    
   for j=1:N
       
       datacube_slow(j, :, i) = fft(datacube(j, :, i));  
   
   end
   
end

for i=1:M
    
   for j=1:K
       
       datacube_fast(:, j, i) = fft(datacube(:, j, i));  
   
   end
   
end

for i=1:K
    
   for j=1:N
       
       datacube_spatial(j, i, :) = fft(datacube(j, i, :));  
   
   end
   
end


