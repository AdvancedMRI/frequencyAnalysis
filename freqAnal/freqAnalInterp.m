function overlayImages = freqAnalInterp(thisView,scanNum,...
  baseDims,analysisNum,overlayCoords,interpMethod,interpExtrapVal,overlayList)
%
% freqAnalInterp: special case for freqAnal. Need to treat amp/ph as complex
% valued.

if ieNotDefined('overlayList')
  overlayList = 1:viewGet(thisView,'numberofoverlays',analysisNum);
end

% Initialize
overlayImages = zeros(size(overlayCoords,1), length(overlayList));

cOverlay=0;
for iOverlay = overlayList
  cOverlay = cOverlay+1;
  if iOverlay
    overlayData = viewGet(thisView,'overlayData',scanNum,iOverlay,analysisNum);
    %check if overlay is an amplitude or a phase
    % if it s the case, load both and interpolate the complex value
    % we assume that phase overlays always follow amplitude overlays
    % Doing a bit more of book keeping, we could avoid interpolating twice 
    % the phase and amplitude overlays, 
    % but I don't think it's worth the trouble
    overlayType = viewGet(thisView,'overlayType',iOverlay,analysisNum);
    switch(overlayType)
      case 'amp'
        if ~strcmp(viewGet(thisView,'overlayType',iOverlay+1,analysisNum),'ph')
%           mrWarnDlg(['(freqAnalInterp) Could not find phase associated to overlay' viewGet(thisView,'overlayName',iOverlay,analysisNum) ', using non-complex interpolation']);
        else
          amplitude = overlayData;
          phase = viewGet(thisView,'overlayData',scanNum,iOverlay+1,analysisNum);
          overlayData = amplitude .* exp(1i*phase);
        end
      case 'ph'
        phase = overlayData;
        if ~strcmp(viewGet(thisView,'overlayType',iOverlay-1,analysisNum),'amp')
          mrWarnDlg(['(freqAnalInterp) Could not find amplitude associated to overlay ' viewGet(thisView,'overlayName',iOverlay,analysisNum) ', using non-weighted phase interpolation']);
          overlayData = exp(1i*phase);
        else
          amplitude = viewGet(thisView,'overlayData',scanNum,iOverlay-1,analysisNum);
          overlayData = amplitude .* exp(1i*phase);
        end
    end

    if ~isempty(overlayData) && ~isempty(overlayCoords)
      % Extract the slice
      overlayImages(:,cOverlay) = interp3(overlayData,...
        overlayCoords(:,2),overlayCoords(:,1),overlayCoords(:,3),...
        interpMethod,interpExtrapVal);
      switch(overlayType)
        case 'amp'
          overlayImages(:,cOverlay) = abs(overlayImages(:,cOverlay));

        case 'ph'
          ang = angle(overlayImages(:,cOverlay));
          ang(ang < 0) = ang(ang < 0)+2*pi;
          overlayImages(:,cOverlay) = ang;
      end
    else
      overlayImages(:,cOverlay) = NaN;
    end
    
  else
    overlayImages(:,cOverlay) = NaN;
  end
end

