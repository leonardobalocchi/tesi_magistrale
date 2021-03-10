%% originale

fs = 500e3;
sLFM = phased.LinearFMWaveform('SampleRate',fs,...
    'SweepBandwidth',200e3,...
    'PulseWidth',1e-3,'PRF',1e3);

lfmwav = step(sLFM);
nsamp = size(lfmwav,1);
t = [0:(nsamp-1)]/fs;
plot(t*1000,real(lfmwav))
xlabel('Time (millisec)')
ylabel('Amplitude')
grid

%% Waveform Specs
fs = 640e3;
sweep_bandwidth = 670e6;
sweep_slope = 21e12; %supponendo che in matlab lo slope sia in secondi e non in microsecondi, quindi devo convertire in secondi
pulse_width = sweep_bandwidth/sweep_slope;
prf = 3200;
sLFM = phased.LinearFMWaveform('SampleRate',fs,...
    'SweepBandwidth',sweep_bandwidth,...
    'PulseWidth',pulse_width,'PRF',prf);

lfmwav = step(sLFM);
nsamp = size(lfmwav,1);
t = [0:(nsamp-1)]/fs;
plot(t*1000,real(lfmwav),'o')
xlabel('Time (millisec)')
ylabel('Amplitude')
grid

