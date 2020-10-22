% modelFreqanalFromGlm
%
%        $Id: computeFreqanal.m 1950 2010-12-18 10:12:48Z julien $
%      usage: [co, amp, ph] = computeFreqanal(tSeries,params)
%         by: julien besle
%       date: 2010-11-15
%    purpose: computes correlation analysis on timeseries array, called by FreqAnal, FreqAnalPlot and modelFreqanalFromGLM
%        e.g:
%
function [co, amp, co2, amp2,ptSeries] = computeFreqanal(tSeries,nCycles,detrend,spatialNormalization,trigonometricFunction,TR)

nFrames = size(tSeries,1);
cHz=nCycles/60;
% Set highpassPeriod
highpassPeriod = round(nFrames/nCycles);

% Remove dc, convert to percent, detrend, and spatial normalization
warnState = warning('query','MATLAB:divideByZero');
ptSeries = percentTSeries(tSeries,...
    'detrend', detrend,...
    'highpassPeriod', highpassPeriod,...
    'spatialNormalization', spatialNormalization,...
    'subtractMean', 'Yes',...
    'temporalNormalization', 'No');
warning(warnState.state,warnState.identifier);

% Compute Fourier transform
ft = fft(ptSeries);

ft = ft(1:1+fix(size(ft, 1)/2), :);
ampFT1 = 2*abs(ft)/nFrames;
L=size(tSeries,1);
Fs=1/TR;
f1 = Fs/2*linspace(0,1,size(ft,1));


% Compute co and amp (avoiding divide by zero)
% amp1 = ampFT1(nCycles+1,:);
% co1 = zeros(size(amp1),'single');
% sumAmp1 = sqrt(sum(ampFT1.^2));
% nonzeroIndices = find(sumAmp1 >0);
% co1(nonzeroIndices) = ampFT1(nCycles+1,nonzeroIndices) ./ sumAmp1(nonzeroIndices);

% Calculate phase:
% switch(trigonometricFunction)
%   case 'Sine'
%     % 1) add pi/2 so that it is in sine phase.
%     % 2) minus sign because sin(x-phi) is shifted to the right by phi.
%     ph = -(pi/2) - angle(ft(nCycles+1,:));   %
%   case 'Cosine'
%     ph = - angle(ft(nCycles+1,:));   
% end
% % 3) Add 2pi to any negative values so phases increase from 0 to 2pi.
% ph(ph<0) = ph(ph<0)+pi*2;

% snr = ampFT1(nCycles+1,:)./mean(ampFT1(round(2*size(ampFT1,1)/3:size(ampFT1,1)),:));

% keyboard
L=size(ptSeries,1);
NFFT=2^nextpow2(L);
Fs=1/TR;

ft = fft(ptSeries,NFFT);
f1=[0:size(ft,1)]*Fs/size(ft,1);

ft = ft(1:1+fix(size(ft, 1)/2), :);
f1 = f1(1:1+fix(size(f1, 2)/2));
ampFT=2*abs(ft)/nFrames;
amp = sum(ampFT(f1 < .5 & f1 >.008,:),1);
% co = zeros(size(amp),'single');
sumAmp = sqrt(sum(ampFT.^2));
% nonzeroIndices = find(sumAmp >0);
co=trapz(f1(f1 < .5 & f1 >.008),ampFT(f1 < .5 & f1 >.008,:),1)./trapz(f1,ampFT,1);
amp=trapz(f1(f1 < .5 & f1 >.008),ampFT(f1 < .5 & f1 >.008,:),1);
co2=trapz(f1(f1 < (cHz+.25) & f1 >(cHz-.25)),ampFT(f1 < (cHz+.25) & f1 >(cHz-.25),:),1)./trapz(f1,ampFT,1);
co2 = sum(ampFT(f1 < (cHz+.25) & f1 >(cHz-.25),:),1)./sum(ampFT,1);
amp2=trapz(f1(f1 < (cHz+.25) & f1 >(cHz-.25)),ampFT(f1 < (cHz+.25) & f1 >(cHz-.25),:),1);
% keyboard
% L=size(ptSeries,1);
% NFFT=2^nextpow2(L);
% f = Fs/2*linspace(0,1,NFFT/2+1);
% fftTSeries = fft(ptSeries,NFFT)/L;
% fftTSeries = fftTSeries(1:NFFT/2+1,:);
% ampFT_2 = 2*abs(fftTSeries);



% absft = abs(ft(1:round((length(ft)/2))));
% signalAmp = absft(nCycles+1);
% noiseFreq = round(2*length(absft)/3):length(absft);
% noiseAmp = mean(absft(noiseFreq));
% snr = signalAmp/noiseAmp;