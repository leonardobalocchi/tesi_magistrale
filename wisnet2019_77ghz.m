%% System definition
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
fastTimeSamples=256; % sample
fastTimeRate=fastTimeSamples/sweepTime; % sample/s
tiltAngle=deg2rad(50);

% Absolute maximum values
maximumRange=['The maximum range is: ',num2str(lightSpeed*fastTimeSamples/4/bandwidth),' m'];
disp(maximumRange)
maximumSpacing=['The maximum spacing for antenna elements is: ',num2str(lightSpeed/2/centerFreq*1e3),' mm'];
disp(maximumSpacing)
maximumTargetVelocity=['The maximum target velocity is: ',num2str(lightSpeed/4/centerFreq/sweepTime/cos(tiltAngle)),' m/s'];
disp(maximumTargetVelocity)

%% Antenna definition
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
%%
figure('Name','Spectrogram')
spectrogram(waveform,'yaxis') % to be fixed

%% Transmitter
TX=phased.Transmitter('Gain',20);

%% Target Specs
targetRCS=0.1;
TgtModel=phased.RadarTarget;
tgtPos=[100*sqrt(3);100;0];              % Target at 20 km distance, 30 degree azimuth
tgtVel=[75*sqrt(3);75;0];                  % Radial velocity is 150 m/sec