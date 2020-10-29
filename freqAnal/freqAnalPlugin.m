% freqAnalPlugin.m
%
%        $Id: freqAnalPlugin.m 1969 2010-12-19 19:14:32Z julien $ 
%      usage: freqAnalPlugin(action,<thisView>)
%         by: alex beckett
%       date: 05/2013
%
function retval = freqAnalPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help freqAnalPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(freqAnalPlugin) Need a valid view to install plugin'));
  else
    %install menu Item
    mlrAdjustGUI(thisView,'add','menu','Frequency Analysis','/Analysis/Time Series Statistics','callback',@freqAnal_Callback);
    % Install default interrogators
    mlrAdjustGUI(thisView,'add','interrogator',{'freqAnalPlot'});
    
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Adds an item in Menu ''Analysis'' to do frequency analysis';
 otherwise
   disp(sprintf('(freqAnalPlugin) Unknown command %s',action));
end


% --------------------------------------------------------------------
function freqAnal_Callback(hObject, eventdata)

view = viewGet(getfield(guidata(hObject),'viewNum'),'view');
view = freqAnal(view);
