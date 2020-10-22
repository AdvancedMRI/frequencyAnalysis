function freqAnalPlot(view,overlayNum,scan,x,y,s,roi)

% check arguments
if ~any(nargin == [1:7])
  help eventRelatedPlot
  return
end


% select the window to plot into
selectGraphWin;
global MLR;
fignum = MLR.graphFigure;

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','freqAnalPlot');

junkFrames = viewGet(view, 'junkFrames', scan);
nFrames = viewGet(view,'nFrames',scan);
framePeriod = viewGet(view,'framePeriod',scan);


tSeries = squeeze(loadTSeries(view,scan,s,[],x,y));
jSeries = tSeries(1:junkFrames);
tSeries = tSeries(junkFrames+1:junkFrames+nFrames);

xx=viewGet(view,'co',scan);
ncycles=xx.params.ncycles(scan);
cHz=ncycles/60;


time = linspace(framePeriod,nFrames*framePeriod,nFrames)';
Fs=1/framePeriod;
L=length(tSeries);
NFFT=2^nextpow2(L);
f = Fs/2*linspace(0,1,NFFT/2+1);

if length(tSeries)>100
    d=fdesign.bandpass('N,F3dB1,F3dB2',10,(cHz-.25),cHz+.25,Fs);
else
    d=fdesign.lowpass('Fp,Fst,Ap,Ast',2,4,2,60,Fs);
end
bpHd=design(d,'butter');
d=fdesign.lowpass('N,F3dB',10,0.5,Fs);

lpHd=design(d,'butter');
d=fdesign.bandpass('N,F3dB1,F3dB2',6,.1,0.5,Fs);
bp2Hd=design(d,'butter');

d=fdesign.lowpass('N,F3dB',10,4,Fs);
lpHd2=design(d,'butter');

nCols = 1;
if ~isempty(roi)
    % get the roi
    roi = roi{1};
    nCols = 2;
    % get time series for roi
    roi = loadROITSeries(view,roi,[],[],'keepNAN','true');
    roi.tSeriesMean = mean(roi.tSeries,1);
    roi.tSeriesMean = roi.tSeriesMean(junkFrames+1:junkFrames+nFrames);
    minTs=min(roi.tSeriesMean);maxTs=max(roi.tSeriesMean);meanTs=mean(roi.tSeriesMean);
    % plot mean time series
    subplot(2,nCols,2);
    plot(time,roi.tSeriesMean,'k-');
    hold on
    %filter time series
    bpTs=filtfilthd(bpHd,roi.tSeriesMean);
    maxBp=max(bpTs);minBp=(min(bpTs));meanBp=mean(bpTs);
    bp2Ts=filtfilthd(bp2Hd,roi.tSeriesMean);
    maxBp2=max(bp2Ts(50:end-50));minBp2=(min(bp2Ts(50:end-50)));meanBp2=mean(bp2Ts(50:end-50));
    lpTs=filtfilt(lpHd.sosMatrix,lpHd.ScaleValues,roi.tSeriesMean);
    lpTs2=filtfilt(lpHd2.sosMatrix,lpHd2.ScaleValues,roi.tSeriesMean);
    maxLp=max(lpTs(50:end-50));minLp=(min(lpTs(50:end-50)));meanLp=mean(lpTs(50:end-50));
    plot(time,lpTs,'b','linewidth',2);
    plot(time,lpTs2,'k','linewidth',2);
    plot(time,bpTs,'r','linewidth',2);
    line([min(time) max(time)],[0 0],'color','k','LineStyle','--')

    % load phys data and plot
    try load(stripext([viewGet(view,'groupname'),'/TSeries/',viewGet(view,'tSeriesFile')]))
        plot(time,respTS(s,junkFrames+1:junkFrames+nFrames),'g');%/max(abs(respTS(s,junkFrames+1:junkFrames+nFrames))),'g')
    catch
    end
    
    stimFile=viewGet(view,'stimfilename');
   
    
    if ~isempty(stimFile)   
        try 
            load(stimFile{1})        
        catch
            disp('No stimfile')
        end
        resp=acq.data((acq.data(1:end-1,4)-acq.data(2:end,4)>1),3);
        resp=resp/max(resp);
        [PKS,LOCS]= findpeaks(-resp,'MINPEAKDISTANCE',30,'MINPEAKHEIGHT',0.1);
        plot(time,resp(junkFrames+1:junkFrames+nFrames)+mean(PKS),'g','linewidth',2);
        plot(time(LOCS),-PKS+mean(PKS),'og')
        for i=1:length(LOCS)-1
            respAvg1(i,:)=resp(LOCS(i):LOCS(i)+round(mean(LOCS(2:end)-LOCS(1:end-1))));
            tSAvg(i,:)=roi.tSeriesMean(LOCS(i):LOCS(i)+round(mean(LOCS(2:end)-LOCS(1:end-1))));
            lpAvg(i,:)=lpTs(LOCS(i):LOCS(i)+round(mean(LOCS(2:end)-LOCS(1:end-1))));
        end
    end
    xlabel('Time (s)');
    ylabel('Velocity (cm/s)');
    %ylabel('Phase');
    axis tight;
    title(sprintf('Mean tSeries of %s, max= %f, min = %f mean = %f \n card: max = %f min = %f mean = %f \n lp resp: max = %f min = %f mean = %f \n bp resp: max = %f min = %f mean = %f',roi.name,maxTs,minTs,meanTs,maxBp,minBp,meanBp,maxLp,minLp,meanLp,maxBp2,minBp2,meanBp2));
    fprintf('\nMean tSeries of %s, max/min/mean %f\t%f\t%f \ncard: max/min/mean = %f\t%f\t%f \nlp resp: max/min/mean = %f\t%f\t%f\nbp resp: max/min/mean = %f\t%f\t%f\n',roi.name,maxTs,minTs,meanTs,maxBp,minBp,meanBp,maxLp,minLp,meanLp,maxBp2,minBp2,meanBp2);

    
    
    subplot(2,nCols,4);
   
    roi.fftTSeriesMean = fft(roi.tSeriesMean,NFFT)/L;
    %roi.fftTSeriesMean(1) = 0;
    plot(f,abs(roi.fftTSeriesMean(1:NFFT/2+1)),'k-');
    hold on, plot(f(f < .5 & f >.000),abs(roi.fftTSeriesMean(f < .5 & f >.000)),'b-','linewidth',2)
    hold on, plot(f(f < (cHz+.25) & f >(cHz-.25)),abs(roi.fftTSeriesMean(f < (cHz+.25) & f >(cHz-.25))),'r-','linewidth',2)

    title(sprintf('Mean FFT of %s',roi.name));
    title(sprintf('resp power = %f, cardiac power =  %f',sum(abs(roi.fftTSeriesMean(f < .5 & f >.008)))/sum(abs(roi.fftTSeriesMean(f < 4.5 & f >.008))),sum(abs(roi.fftTSeriesMean(f < (cHz+.25) & f >(cHz-.25))))/sum(abs(roi.fftTSeriesMean(f < 4.5 & f >.008)))));
    fprintf('resp power/card power = %f\t%f\n',sum(abs(roi.fftTSeriesMean(f < .5 & f >.008)))/sum(abs(roi.fftTSeriesMean(f < 4.5 & f >.008))),sum(abs(roi.fftTSeriesMean(f < (cHz+.25) & f >(cHz-.25))))/sum(abs(roi.fftTSeriesMean(f < 4.5 & f >.008))));

    xlabel('FFT components');
    ylabel('FFT of Velocity (cm/s)');
    axis tight;
    
%     try
%         respVar=load('respVar');
%         resp=respVar.respR{s,scan};
%         B = pinv([resp ones(size(resp,1),1)])*roi.tSeriesMean';
%         subplot(2,nCols,2);
%         hold on
%         plot(time,respVar.respR{s,scan}/max(respVar.respR{s,scan}),'b','linewidth',2)
%         subplot(2,nCols,4);
%         hold on
%         ft=fft([resp ones(size(resp,1),1)]*B,NFFT)/L;
%         ft=fft(respVar.respR{s,scan}/max(respVar.respR{s,scan}),NFFT)/L;
%         ft(1)=0;
%         plot(f,abs(ft(1:NFFT/2+1)),'b','linewidth',2)
%         axis('tight')
%     catch
%         %disp('Error in resp fitting')
%     end
%     try
%         respVar=load('respVar');
%         
%         card=respVar.resp2SR{s,scan};
%         B=[];
%         for i=1:size(card,1)
%             B(i,:)=pinv([card(i,:)' ones(size(card,2),1)])*roi.tSeriesMean';
%         end
%         [~,j]=max(B(:,1));
%         subplot(2,nCols,2);
%         hold on
% %         plot(time,[card(j,:)' ones(size(card,2),1)]*B(j,:)','r--')
%         subplot(2,nCols,4);
%         hold on
%         ft=fft([card(j,:)' ones(size(card,2),1)]*B(j,:)',NFFT)/L;
% %         ft(1)=0;
%         plot(f,abs(ft(1:NFFT/2+1)),'r--');
%     catch
%         %disp('Error in Cardiac fitting')
%     end
end

% get the overlays
titleStr = sprintf('Voxel: [%i, %i, %i]',x,y,s);
for iOverlay = 1:viewGet(view,'nOverlays');
  % get name
  overlayName = viewGet(view,'overlayName',iOverlay);
  % get value
  overlay = viewGet(view,'overlayData',scan,iOverlay);
  overlayValue = overlay(x,y,s);
  % add to titleStr
  titleStr = sprintf('%s %s=%f',titleStr,overlayName,overlayValue);
  if iOverlay == 3
    titleStr = sprintf('%s\n',titleStr);
  end
end

subplot(2,nCols,1);
plot(time,tSeries, 'k.-');
hold on
plot(time,filtfilthd(bpHd,tSeries),'r');
plot(time,filtfilthd(lpHd,tSeries),'b');
lpTS=filtfilthd(lpHd,tSeries);

try load(stripext([viewGet(view,'groupname'),'/TSeries/',viewGet(view,'tSeriesFile')]))
        plot(time,respTS(s,junkFrames+1:junkFrames+nFrames),'g');%/max(abs(respTS(s,junkFrames+1:junkFrames+nFrames))),'g')
catch
end
stimFile=viewGet(view,'stimfilename');
    try load(stimFile{1})
        resp=acq.data((acq.data(1:end-1,4)-acq.data(2:end,4)>1),3);
        resp=resp/max(resp);
        [PKS,LOCS]= findpeaks(-resp,'MINPEAKDISTANCE',30,'MINPEAKHEIGHT',0.1);
        plot(time,resp(junkFrames+1:junkFrames+nFrames)+mean(PKS),'g','linewidth',2);
        plot(time(LOCS),-PKS+mean(PKS),'og')
        %line([time(LOCS);time(LOCS)],repmat([-8;8],[1 length(LOCS)]),'color','g')
    catch
        disp('No stimfile')
    end

title(titleStr);
xlabel('Time (s)');
%ylabel('Velocity (cm/s)');
ylabel('Phase');
axis tight;
% draw borders between rund
concatInfo = viewGet(view,'concatInfo',scan);
if ~isempty(concatInfo)
  vline(concatInfo.runTransition(2:end,1)-1,'r-');
end

subplot(2,nCols,nCols+1);
fftTSeries = fft(tSeries,NFFT)/L;
% set mean to zero
% fftTSeries(1) = 0;
% plot it 
plot(f,abs(fftTSeries(1:NFFT/2+1)),'k.-');
hold on, plot(f(f < .5 & f >.000),abs(fftTSeries(f < .5 & f >.000)),'bo-')
hold on, plot(f(f < (cHz+.25) & f >(cHz-.25)),abs(fftTSeries(f < (cHz+.25) & f >(cHz-.25))),'ro-')
%hold on, plot(f(f < (cHz+.25) & f >(cHz-.25)),abs(fftTSeries(f < (cHz+.25) & f >(cHz-.25))),'ro-')

title(sprintf('Voxel: [%i %i %i]',x,y,s));
xlabel('FFT components');
ylabel('FFT of Velocity (cm/s)');
axis tight;


try
    respVar=load('respVar');
        % X=respVar.respSR{s,scan}';
%         resp=respVar.respSR{s,scan};
%         for i=1:size(resp,1)
%             B(i,:)=pinv([resp(i,:)' ones(size(resp,2),1)])*roi.tSeriesMean';
%         end
%         [~,j]=max(B(:,1));
        resp=respVar.respR{s,scan};
        B = pinv([resp ones(size(resp,1),1)])*tSeries;
        subplot(2,nCols,1);
        hold on
%         plot(time,[resp(j,:)' ones(size(resp,2),1)]*B(j,:)','b--')
        plot(time,[resp]);% ones(size(resp,1),1)]*B,'b--')

%         plot(time,respVar.respR{s,scan}*B,'k')
        subplot(2,nCols,nCols+1);
        hold on
%         ft=fft([resp(j,:)' ones(size(resp,2),1)]*B(j,:)',NFFT)/L;
        ft=fft([resp ones(size(resp,1),1)]*B,NFFT)/L;
%         ft(1)=0;
        plot(f,abs(ft(1:NFFT/2+1)),'b--');
catch
    %disp('Error in Resp fitting')
end
try
    respVar=load('respVar');
    
    % X=respVar.respSR{s,scan}';
    card=respVar.resp2SR{s,scan};
    B=[];
    for i=1:size(card,1)
        B(i,:)=pinv([card(i,:)' ones(size(card,2),1)])*tSeries;
    end
    [~,j]=max(B(:,1));
    subplot(2,nCols,1);
    hold on
    plot(time,[card(j,:)' ones(size(card,2),1)]*B(j,:)','r--')
    subplot(2,nCols,nCols+1);
    hold on
    ft=fft([card(j,:)' ones(size(card,2),1)]*B(j,:)',NFFT)/L;
%     ft(1)=0;
    plot(f,abs(ft(1:NFFT/2+1)),'r--');
catch
    %disp('Error in Cardiac fitting')
end

if ~isempty(roi)
    %setFont4CSF;
end

if ~isempty(stimFile) & ~isempty(roi)

figure,errorbar(time(1:size(respAvg1,2)),mean(respAvg1,1)+mean(PKS),std(respAvg1,[],1),'g')
hold on,errorbar(time(1:size(respAvg1,2)),mean(tSAvg,1),std(tSAvg,[],1),'k')
errorbar(time(1:size(respAvg1,2)),mean(lpAvg,1),std(lpAvg,[],1),'b')
line([min(time) max(time(1:size(respAvg1,2)))],[0 0],'color','k','LineStyle','--')

title(sprintf('RespAvg of %s, max= %f, min = %f mean = %f \n lp resp: max = %f min = %f mean = %f',roi.name,max(mean(tSAvg,1)),min(mean(tSAvg,1)),mean(mean(tSAvg,1)),max(mean(lpAvg,1)),min(mean(lpAvg,1)),mean(mean(lpAvg,1))));
end
    
return

