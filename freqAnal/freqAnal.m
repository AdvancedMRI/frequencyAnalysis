function [view params] = freqAnal(view,params,varargin)
%
% view = freqAnal(view,[params])
% 
% Loops throughs scans and slices, loads corresponding tSeries, computes
% correlation analysis from tSeries, and saves the resulting co, amp, and
% ph to the freqAnal.mat file along with the analysis parameters.
%
% Checks to see of co, amp, and ph are already loaded. If so, it uses the
% existing maps and recomputes as specified.
%
% params: Optional initial parameters. Default: user is prompted via
%    freqAnalGUI. If a freqAnal already exists and is loaded then the
%    freqAnalGUI is initialized with the existing parameters. Params must be
%    a structure with all of the following fields.
% groupName: group of scans that will be averaged.
%    Default: current group of view.
% recompute: vector of 1 and 0 specifying which scans to compute.
%    Default: all of the scans.
% ncycles: vector specifying number of cycles per scan.
%    Default: 1
% detrend: cell array of strings as in percentTSeries.
%    Default: 'None'
% spatialnorm: cell array of strings as in percentTSeries.
%    Default: 'None'
% tseriesfile: cell array of strings specifying tseries filenames. Or
%    'any' to allow any match with the tseries files. This is useful so
%    that you can run the analysis on multiple data sets using the same
%    params.
%
%
% Examples:
%
% params = freqAnalGUI('groupName','Raw');
% view = newView;
% view = freqAnal(view,params);
%
% view = freqAnal(view);
%
% To just get parameters
% [view params] = freqAnal(view,[],'justGetParams=1','defaultParams=1');
%
% djh, 5/2005, updated to mrLoadRet-4.0
% $Id: freqAnal.m 2684 2013-02-14 15:38:55Z julien $	

% check arguments
if ~any(nargin == [1 2 3 4 5])
  help freqAnal
  return
end

if ~isview(view)
    help freqAnal
    mrErrorDlg('(freqAnal) Invalid view.')
end

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end

% If freqAnal is loaded, then use it. Otherwise, viewGet returns [];
freqAnal = viewGet(view,'freqAnal');
co = viewGet(view,'co');
amp = viewGet(view,'amp');
co2 = viewGet(view,'co2');
amp2 = viewGet(view,'amp2');
% ph = viewGet(view,'ph');
% snr = viewGet(view,'snr');
if ~isempty(freqAnal)
    oldparams = freqAnal.params;
else
    oldparams = [];
end

% if we are just getting default parameters then
% get them by calling the reconcile function
if defaultParams
  params = freqAnalReconcileParams(viewGet(view,'groupName'),[]);
end

% Get analysis parameters from freqAnalGUI, using co.params if it exists
if ieNotDefined('params')
  if ~isempty(oldparams)
    % Initialize analysis parameters with previous values from loaded
    % coherence map and reconcile with current status of group.
    params = freqAnalGUI('groupName',viewGet(view,'groupName'),'params',oldparams);        
  else
    % Initialize analysis parameters with default values
    params = freqAnalGUI('groupName',viewGet(view,'groupName'));
  end
else
  % Reconcile params with current status of group and ensure that params
  % has the required fields.
  params = freqAnalReconcileParams(params.groupName,params);
end

% Abort if params empty
if ieNotDefined('params')
    mrMsgBox('freqAnal cancelled',1);
    return
end

% if just getting params then return
if justGetParams,return,end

% Change group, get nScans
groupName = params.groupName;
curGroup = viewGet(view,'currentGroup');
groupNum = viewGet(view,'groupNum',groupName);
if (groupNum ~= curGroup)
	mrWarnDlg(['Changing view to group: ',groupName]);
	view = viewSet(view,'currentGroup',groupNum);
end
nScans = viewGet(view,'nScans',groupNum);

% If co/amp/ph do not exist (above viewGet calls return [])...
% Intialize structure and initialize data to cell array of length nScans.
if isempty(co)
    co.name = 'co';
    co.function = 'freqAnal';
    co.reconcileFunction = 'freqAnalReconcileParams';
    co.mergeFunction = 'freqAnalMergeParams';
    co.data = cell(1,nScans);
end
if isempty(co2)
    co2.name = 'co2';
    co2.function = 'freqAnal';
    co2.reconcileFunction = 'freqAnalReconcileParams';
    co2.mergeFunction = 'freqAnalMergeParams';
    co2.data = cell(1,nScans);
end
if isempty(amp)
    amp.name = 'amp';
    amp.function = 'freqAnal';
    amp.reconcileFunction = 'freqAnalReconcileParams';
    amp.mergeFunction = 'freqAnalMergeParams';
    amp.data = cell(1,nScans);
end
if isempty(amp2)
    amp2.name = 'amp2';
    amp2.function = 'freqAnal';
    amp2.reconcileFunction = 'freqAnalReconcileParams';
    amp2.mergeFunction = 'freqAnalMergeParams';
    amp2.data = cell(1,nScans);
end
% if isempty(ph)
%     ph.name = 'ph';
%     ph.function = 'freqAnal';
%     ph.reconcileFunction = 'freqAnalReconcileParams';
%     ph.mergeFunction = 'freqAnalMergeParams';
%     ph.data = cell(1,nScans);
% end
% if isempty(snr)
%     snr.name = 'snr';
%     snr.function = 'freqAnal';
%     snr.reconcileFunction = 'freqAnalReconcileParams';
%     snr.mergeFunction = 'freqAnalMergeParams';
%     snr.data = cell(1,nScans);
% end

% Compute it
[co,amp,co2,amp2] = computeFreqAnalysis(view,params,co,amp,co2,amp2);

% Set params field, merging with previous params
if ~isempty(oldparams)
    params = freqAnalMergeParams(groupName,oldparams,params);
end
co.params = params;
amp.params = params;
co2.params = params;
amp2.params = params;
    
% Fill range fields, fixed values for co and ph
co.range = [0 .25];
co.clip = [0 1];
co2.range = [0 .25];
co2.clip = [0 1];

% Fill range field for amp
ampMin = realmax;
ampMax = 0;
amp2Min = realmax;
amp2Max = 0;
for scan=1:nScans
    if ~isempty(amp.data{scan})
        ampMin = min([ampMin min(amp.data{scan}(:))]);
        ampMax = max([ampMax max(amp.data{scan}(:))]);
    end
    if ~isempty(amp2.data{scan})
        amp2Min = min([amp2Min min(amp2.data{scan}(:))]);
        amp2Max = max([amp2Max max(amp2.data{scan}(:))]);
    end
end
if (ampMin <= ampMax)
    amp.range = [ampMin ampMax];
else
    % if amp data is empty, need to make sure min < max
    amp.range = [0 1];
end

if (amp2Min <= amp2Max)
    amp2.range = [amp2Min amp2Max];
else
    % if amp data is empty, need to make sure min < max
    amp2.range = [0 1];
end

% Fill colormap fields, keeping previous if already exists
if ~isfield(co,'colormap')
    co.colormap = jet(256);
end
if ~isfield(amp,'colormap')
    amp.colormap = hot(256);
end
if ~isfield(co2,'colormap')
    co2.colormap = jet(256);
end
if ~isfield(amp2,'colormap')
   amp2.colormap = hot(256);
end

% Set date field
dateString = datestr(now);
co.date = dateString;
amp.date = dateString;
co2.date = dateString;
amp2.date = dateString;

% Set interrogrator function
co.interrogator = 'freqAnalPlot';
amp.interrogator = 'freqAnalPlot';
co2.interrogator = 'freqAnalPlot';
amp2.interrogator = 'freqAnalPlot';

% Set groupName
co.groupName = params.groupName;
amp.groupName = params.groupName;
co2.groupName = params.groupName;
amp2.groupName = params.groupName;

% Install freqAnal in the view
freqAnal.name = 'freqAnal';  % This can be reset by editAnalysisGUI
freqAnal.type = 'freqAnal';
freqAnal.groupName = params.groupName;
freqAnal.function = 'freqAnal';
freqAnal.reconcileFunction = 'freqAnalReconcileParams';
freqAnal.mergeFunction = 'freqAnalMergeParams';
freqAnal.guiFunction = 'freqAnalGUI';
freqAnal.overlayInterpFunction = [];%'freqAnalInterp';
freqAnal.params = params;
freqAnal.date = dateString;
view = viewSet(view,'newanalysis',freqAnal);
view = viewSet(view,'newoverlay',co);
view = viewSet(view,'newoverlay',amp);
view = viewSet(view,'newoverlay',co2);
view = viewSet(view,'newoverlay',amp2);
if ~isempty(viewGet(view,'fignum'))
  refreshMLRDisplay(viewGet(view,'viewNum'));
end

% Save it
saveAnalysis(view,freqAnal.name);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [co,amp,co2,amp2] = computeFreqAnalysis(view,params,co,amp,co2,amp2)
% Required fields in params: 'recompute','ncycles','detrend','spatialnorm'

% Get scanList from params.recompute field
scanList = find(params.recompute(:));

disp('Computing freqAnal...');
warning('off','MATLAB:divideByZero');
for scanIndex=1:length(scanList)
  scanNum = scanList(scanIndex);
  waitHandle = mrWaitBar(0,['Computing Fequency Analysis for scan ' int2str(scanNum) ':']);

  % sliceDims: [ydim xdim] for single slice
  % volDims; [ydim xdim nslices] for single scan
  sliceDims = viewGet(view,'sliceDims',scanNum);
  volDims = viewGet(view,'dims',scanNum);

  % Initialize data with NaNs
  co.data{scanNum} = NaN*ones(volDims);
  amp.data{scanNum} = NaN*ones(volDims);
  co2.data{scanNum} = NaN*ones(volDims);
  amp2.data{scanNum} = NaN*ones(volDims);

  % check for freqAnal cycles set to 0
  if params.ncycles(scanList(scanIndex)) == 0
    mrWarnDlg(sprintf('(freqAnal:computeFrequencyAnalysis) !!! Scan %i has ncycles set to 0 - this needs to be set to how many cycles of the stimulus you had per scan. Skipping this scan !!!',scanList(scanIndex)));
    continue;
  end
  
  nslices = viewGet(view,'nslices',scanNum);
  for sliceNum = 1:nslices
    % Analysis parameters for this scan
    junkframes = viewGet(view,'junkframes',scanNum);
    nframes = viewGet(view,'nframes',scanNum);
    % Load tSeries
    tSeries = loadTSeries(view, scanNum, sliceNum);
    % Reshape the tSeries
    % ATTN: added reshapeTSeries function, since loadTSeries not longer reshapes when it loads -eli
    tSeries = reshapeTSeries(tSeries);
    
    % check that junkframes and nframes settings are ok
    if size(tSeries,1) < (junkframes+nframes)
      mrErrorDlg(sprintf('(freqAnal) Number of junkframes (%i) plus nframes (%i) should not be larger than number of volumes in scan %i',junkframes,nframes,size(tSeries,1)));
    end
    % Remove junkFrames
    tSeries = tSeries(junkframes+1:junkframes+nframes,:);
    %compute freqAnal
    TR=viewGet(view,'framePeriod',scanNum);
    [coSeries,ampSeries,co2Series,amp2Series,~] = computeFreqanal(tSeries,params.ncycles(scanNum),params.detrend{scanNum},params.spatialnorm{scanNum},params.trigonometricFunction{scanNum},TR);

    switch view.viewType
      case {'Volume'}
          co.data{scanNum}(:,:,sliceNum) = reshape(coSeries,sliceDims);
          amp.data{scanNum}(:,:,sliceNum) = reshape(ampSeries,sliceDims);
          co2.data{scanNum}(:,:,sliceNum) = reshape(co2Series,sliceDims);
          amp2.data{scanNum}(:,:,sliceNum) = reshape(amp2Series,sliceDims);
      case {'Surface'}
          co.data{scanNum} = coSeries;
          amp.data{scanNum} = ampSeries;
          co2.data{scanNum} = co2Series;
          amp2.data{scanNum} = amp2Series;
      case {'Flat'}
          co.data{scanNum}(:,:,sliceNum) = reshape(coSeries,sliceDims);
          amp.data{scanNum}(:,:,sliceNum) = reshape(ampSeries,sliceDims);
          co2.data{scanNum}(:,:,sliceNum) = reshape(co2Series,sliceDims);
           amp2.data{scanNum}(:,:,sliceNum) = reshape(amp2Series,sliceDims);
    end
    % Update waitbar
    mrWaitBar(sliceNum/nslices,waitHandle);
  end
  mrCloseDlg(waitHandle);
end
warning('on','MATLAB:divideByZero');
