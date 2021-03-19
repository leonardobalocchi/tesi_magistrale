clc
close all
clearvars

centerFreq=77e9; % Hz
bandwidth=2e9; % Hz
sweepTime=500e-6; % s
pulseRepFreq=1/sweepTime; % Hz
lightSpeed=299792458; % m/s
wavelength=lightSpeed/centerFreq; % m
receiverChannels=8;
chirpCycles=20000;
fastTimeSamples=50000; % sample
fastTimeRate=fastTimeSamples/sweepTime; % sample/s

% Antenna definition
antenna=phased.ULA; % Uniform Linear Array
antenna.NumElements=receiverChannels;
antenna.Element=phased.CosineAntennaElement; 
spacing_tuning=0; % Use this parameter to tune the antenna radiation pattern.
antenna.ElementSpacing=wavelength/2+spacing_tuning;

figure('Name','Array geometry')
viewArray(antenna)
figure('Name','3D Directivity')
pattern(antenna,centerFreq,'Type','directivity','PropagationSpeed',lightSpeed)
% pattern(antenna,centerFreq,-180:180,0,'Type','directivity','PropagationSpeed',lightSpeed)
figure('Name','Directivity at 0° elevation')
pattern(antenna,centerFreq,-180:180,0,'Type','directivity','PropagationSpeed',lightSpeed)

%% Waveform specifications
waveformSysObj=phased.LinearFMWaveform;
waveformSysObj.SampleRate = fastTimeRate;
waveformSysObj.PRF=pulseRepFreq;
waveformSysObj.PulseWidth=sweepTime;
numSamples=fastTimeRate*sweepTime;
waveform=step(waveformSysObj);
timeAxis=(0:numSamples-1)/fastTimeRate;

figure('Name','Radar waveform')
plot(timeAxis,real(waveform),'.');title('Radar Waveform'); xlabel('Time (sec)'); ylabel('Amplitude')
prf = waveform.PRF;
