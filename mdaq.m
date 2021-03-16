function [out,out1,out2,out3,out4]=mdaq(varargin)

% Driver for use of MATLAB Data acquisition toolbox with SDT

obj=[];evt=[];
if nargin==0; CAM=''; Cam=''; carg=1;
elseif comstr(varargin{1},'reset');[CAM,Cam]=comstr(varargin{1},6);carg=2;
elseif nargin>2&&ischar(varargin{3})&&~ischar(varargin{1})
 obj=varargin{1};evt=varargin{2};
 [CAM,Cam]=comstr(varargin{3},1);carg=4;
else
 [CAM,Cam]=comstr(varargin{1},1);carg=2;
end

if comstr(Cam,'init');[CAM,Cam]=comstr(CAM,5);
 %% #Init -----------------------------------------------------------------1
 [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
  '-reset(#3#"To reset the displayed tab")' ...
  '-fig(10#%g#"Figure number to display time acquisitions")' ...
  ],{struct,CAM}); Cam=lower(CAM);
 [ci,cf]=iicom('dockid');
 UI=sdtroot('paramui',ci); PA=sdtroot('paramvh',ci);
 d=PA.Stack{'info','daqsdt'};
 figure(RO.fig);
 cingui('objset',RO.fig,{'@Dock',{'Name','Id','tile',ci}});figure(RO.fig);
 
 if comstr(Cam,'channels')
  %% #InitChannels : add inputs, define impact/accel channels...----------2
  % Tab.DAQ
  % AcqFcn
  % xxx : should be callback
  d = daq("ni"); d.UserData=struct; % When init create new instance
  stack_set(PA,'info','daqsdt',d); % Store in PA for exernal access
  addinput(d,"Dev1","ai0","IEPE");
  addinput(d,"Dev1","ai1","IEPE");
  addinput(d,"Dev1","ai2","IEPE");
  RA.ref='ai0';
  RA.refind=find(strcmpi({d.Channels.ID},RA.ref));
  
  d.UserData.Channels=RA;
  d.UserData.fig=RO.fig;
  % Store PA in d to easily access data from daq callbacks
  d.UserData.PA=PA; 
  d.UserData.ci=ci; 
  
 elseif comstr(Cam,'impact')
  %% #InitImpact : Open impact tab----------------------------------------2
  sdtroot('InitGenericImpact',UI,RO); % Generic init + field editions
  r1j=PA.Impact; r1=cinguj('ObjToStruct',r1j); R1=fe_def('cleanentry',r1);
  wire=cf.CStack{'SensDof',''};
  if ~isempty(wire)
   st=arrayfun(@(x) sprintf('%13.2f',x),wire.tdof(:,1),'uni',0);
   sdcedit(r1j,'SensDof','choices',st);
  end
  iicom(ci,'SetImpact',struct('FSamp',d.Rate));
  
 else; error('Init%s unknown',CAM);  
 end
elseif comstr(Cam,'set');[CAM,Cam]=comstr(CAM,4);
 %% #Set------------------------------------------------------------------1
 % Set Input parsing from sdtproject.m
 if isstruct(evt) && isfield(evt,'val')
  RO=fe_range('ToFlat',obj,evt); % get the input of command from obj and event
 else; RO=[];
 end
 if ~isempty(evt)&&isstruct(evt)&&isfield(evt,'CAM');
   [CAM,Cam]=comstr(evt.CAM,1); RO=evt.RO;
 elseif ~isempty(CAM)||isempty(RO) % generic get of SetFcn
   [RO,uo,CAM,Cam]=clean_get_uf('getuo',['SetStruct' CAM],obj,evt);
 else; RO=struct; % uo remain empty, CAM is specify outside
 end
 
 UI=sdtroot('paramui',obj);PA=sdtroot('paramvh',obj);
 d=PA.Stack{'info','daqsdt'};
 %% set_Commands   
 if comstr(Cam,'impact')
 %% #SetImpact-2
  [r1j,r1,st,PA,u0]=sdth.eMethods.defaultSet(UI,PA,RO,'Impact',...
   {'TriggerLevel','PreTrigger','Averages'});
  R1=fe_def('cleanentry',r1);
  j1=0;out='';
  while j1<length(st); j1=j1+1; val=RO.(st{j1});
   if strcmpi(st{j1},'ContinuousScan')
    %% #ContinuousScan : Continuously read and dispaly channels-----------3
    if d.Running; daqsdt('stop',d);
    else
     R1.AcqType='ContinuousScan';
     ReadRO(d,R1)
     start(d,"Continuous");
    end
   elseif strcmpi(st{j1},'FSamp')
    %% #FSamp : Define acquisition rate-----------------------------------3
    if val<d.RateLimit(1); val=d.RateLimit(1); end
    if val>d.RateLimit(2); val=d.RateLimit(2); end
    d.Rate=val;
    % Asked value can be changed to closest acceptable one->update display
    sdcedit(r1j,'FSamp',round(d.Rate*10)/10);
   elseif strcmpi(st{j1},'AcqTime')
    %% #AcqTime : Define acquisition time (linked with FreqRes)-----------3
    sdcedit(r1j,'AcqTime',round(val*1000)/1000);
    sdcedit(r1j,'FreqRes',round(1/val*1000)/1000);
   elseif strcmpi(st{j1},'FreqRes')
    %% #FreqRes : Define Frequency Resolution (linked with AcqTime)-------3
    sdcedit(r1j,'FreqRes',round(val*1000)/1000);
    sdcedit(r1j,'AcqTime',round(1/val*1000)/1000);   
   elseif strcmpi(st{j1},'SensDof')
    %% #SensDof-3
    if isnumeric(val) % Check if the number is in the sensdof list
     r2=str2num(r1.SensDof.get('choices'));
     i1=find(ismember(r2,val));
     if ~isempty(i1); val=i1; end
    end
    sdcedit(r1j,'SensDof',val);
    R1.SensDof=str2double(fe_def('cleanentry',r1.SensDof));
    cf=feplot(PA.FeplotFig,';');
    wire=cf.CStack{'SensDof',''};
    dof=unique(wire.tdof(:,1),'stable');
    i1=fe_c(dof,R1.SensDof,'ind');
    
    def=struct('def',zeros(size(dof)),'DOF',dof);
    def.def(i1)=1;
    
    cf.def=def;fecom(cf,';undefline;showdefearrow');
   elseif strcmpi(st{j1},'TriggerWait')
    %% #TriggerWait : wait for trigger condition to acquire data----------3
    if d.Running; daqsdt('stop',d); % AcqFcun=@mdaq; feval(AcqFcn,'stop',d)
    else
     R1.AcqType='TriggerWait';
     ReadRO(d,R1)
     start(d,"Continuous");
    end
   elseif strcmpi(st{j1},'Reject')
    %% #Reject : Remove last impact measurement---------------------------3
    frames=PA.Stack{'curve','frames'};
    if ~isempty(frames); frames(end)=[]; end
    stack_set(PA,'curve','frames',frames);
    TriggerPlot(d.UserData,frames);
   elseif strcmpi(st{j1},'Measured')
    %% #Measured : Number of measured impacts-----------------------------3
    sdcedit(r1j,'Measured','value',sprintf('%i/%i',val,R1.Averages));
   elseif strcmpi(st{j1},'Accept')
    %% #Accept : Accept current H1 estimation and store in Test---------3
    ci=d.UserData.ci;
    refind=d.UserData.Channels.refind;
    XF=ci.Stack{'Test'};
    
    % Generate H1
    frames=PA.Stack{'frames'};
    XF2=fe_curve(sprintf('h1h2 %i -stack',refind),frames,struct('Out','H1'));
    XF2.X{2}=str2num(R1.SensDof);
    XF2.X{3}=(1:size(XF2.Y,2))'+.99;
    XF2.Xlab(2:3)={'In' 'Out'};
    XF2.dof=[repmat(XF2.X{2},length(XF2.X{3}),1) XF2.X{3}];
    XF2.idopt=idopt;
    XF2.idopt.DataType='acc';
    % Merge with previous meas. points
    if ~isempty(XF)
     i1=fe_c(XF.X{2},str2num(R1.SensDof),'ind'); % Remove previous meas.
     if ~isempty(i1)
      i1=fe_c(XF.X{2},str2num(R1.SensDof),'ind',2);% Channels to keep
      XF=fe_def('subchcurve',XF,{'In',i1});
     end
     XF=fe_def('curvejoin',XF,XF2);
    else; XF=XF2;
    end
    iicom(ci,'curveinit','Test',XF); % Store in Test
    idcom(ci);
    % Flush frames + next SensDof
    stack_set(PA,'curve','frames',{});
    TriggerPlot(d.UserData,{});
    iicom(ci,'SetImpact',struct('SensDof',r1.SensDof.get('value')+1));
    %% #ImpactEnd -3
   else;error('Impact%s',st{j1});
   end
  end
  cingui('resize',UI.gf);
 %% #SetEnd -2
 else;error('Set%s',CAM);
 end
 % ua.JTable.repaint;
 return;

elseif comstr(Cam,'stop')
 %% #Stop : stop background reading---------------------------------------1
 d=varargin{carg};
 stop(d);flush(d);d.UserData.data=[];
 %% #End------------------------------------------------------------------1
elseif comstr(Cam,'@');out=eval(CAM);
elseif comstr(Cam,'cvs')
 out='$Revision: 1.00 $  $Date: 2020/01/13 09:32:24 $';
else; error('daqsdt %s unknown',CAM);
end
end

%% #SubFunc---------------------------------------------------------------1
function DaqRead(d,evt)
%% DaqRead : Backgroung reading (ContinuousScan,TriggerWait,...)----------2
[data,time] = read(d,d.ScansAvailableFcnCount,"OutputFormat","Matrix");
RO=d.UserData;
RA=RO.Read;
if strcmpi(RA.AcqType,'ContinuousScan')
 if isempty(RO.data); RO.data=zeros(RA.Samples,size(data,2));end
 data2=RO.data;
 data2(1:RA.RefreshSamples,:)=[];
 data2=[data2;data];
 RO.data=data2;
 
 for j1=1:size(data2,2)
  plot(RA.Axes(j1),RA.t,data2(:,j1));
  xlim(RA.Axes(j1),[0 RA.AcqTime]);drawnow;
 end
elseif strcmpi(RA.AcqType,'TriggerWait')
 if isempty(RO.data); RO.data=data; end % First read, store in pretrig band
 
 refind=RO.Channels.refind; % Indice of reference channel
 if any(data(:,refind)>RA.TriggerLevel)
  h=msgbox('Measuring...','modal');

  data=[RO.data;data]; % Concat pretrig band + current band
  i1=find(data(:,refind)>RA.TriggerLevel,1,'first');
  data=data(i1-RA.PreTriggerSamples:end,:); % Begin when trig condition is true 
  i2=RA.Samples-size(data,1); % Read following time samples
  data2=read(d,i2,"OutputFormat","Matrix");
  data2=[data;data2];
   
  frame=struct('X',{{RA.t (1:size(data2,2))'}},'Xlab',{{{'Time' 's' []} 'Channel'}},...
   'Y',data2);
  % Reorder Channels to have first the impact and then the accel
  ind=setdiff(1:size(data,2),refind);ind=[refind ind]; 
  frame=fe_def('subchcurve',frame,{'Channel',ind'});

  % Store with previous frames
  PA=RO.PA; ci=RO.ci;
  frames=PA.Stack{'curve','frames'};
  if isempty(frames); frames={frame}; 
  else; frames=[frames;{frame}];
  end
  stack_set(PA,'curve','frames',frames);
  
  TriggerPlot(RO,frames);
  daqsdt('stop',d); % Stop and flush data  
  
  start(d,"Continuous");
  close(h);
 end
else; sdtw('_IsGui','AcqType %s not handled',RO.AcqType); daqsdt('stop');
end
d.UserData=RO;
end

function ReadRO(d,RO)
%% ReadRO : Assign meas. settings used by DaqRead-------------------------2
% ReadRO(d,R1,'ContinuousScan');
 if ~isfield(RO,'AcqType'); error('Provide measurement AcqType'); end
 RA.AcqType=RO.AcqType;
 
 if strcmpi(RA.AcqType,'ContinuousScan')
  RA.AcqTime=1;
  RA.RefreshTime=.1;

  RA.Samples=ceil(d.Rate*RA.AcqTime);
  RA.RefreshSamples=ceil(d.Rate*RA.RefreshTime);
  RA.t=(1/d.Rate:1/d.Rate:1/d.Rate*RA.Samples)';
  
  % Create axes
  clf(d.UserData.fig);figure(d.UserData.fig)
  i1=length(d.Channels); % Number of channels
  for j1=1:i1
   RA.Axes(1,j1)=subplot(i1,1,j1);axis tight;
  end

  d.ScansAvailableFcnCount = RA.RefreshSamples;
  d.ScansAvailableFcn = daqsdt('@DaqRead');  
 elseif strcmpi(RA.AcqType,'TriggerWait')
  RA.AcqTime=RO.AcqTime;
  RA.PreTrigger=RO.PreTrigger;
  RA.TriggerLevel=RO.TriggerLevel;

  RA.Samples=ceil(d.Rate*RA.AcqTime);
  RA.PreTriggerSamples=ceil(d.Rate*RA.PreTrigger);
  RA.t=(1/d.Rate:1/d.Rate:1/d.Rate*RA.Samples)';
  
  % Create axes
  clf(d.UserData.fig);figure(d.UserData.fig)
  i1=length(d.Channels); % Number of channels
  for j1=1:i1
   for j2=1:3
    RA.Axes(j1,j2)=subplot(i1,3,(j1-1)*3+j2);axis tight;
   end
  end
  
  d.ScansAvailableFcnCount = RA.PreTriggerSamples;
  d.ScansAvailableFcn = daqsdt('@DaqRead');
 else; error('Measurement type "%s" unknown',st);  
 end

 d.UserData.Read=RA;
 d.UserData.data=[];% Initialize previous read data to empty
end

function TriggerPlot(RO,frames)
%% TriggerPlot : Refresh plots after each impact--------------------------2
 % TriggerPlot(d.UserData,PA.Stack{'curve','frames'})
 RA=RO.Read;
 refind=RO.Channels.refind; % Indice of reference channel
 % Update measured fields
 iicom(RO.ci,'SetImpact',struct('Measured',length(frames)));
 if isempty(frames)
  for j1=1:length(RA.Axes(:)); cla(RA.Axes(j1)); end
  return;
 end
 frame=frames{end}; % Frame is the last (current) meas.

 % Plot time measurements
 for j1=1:size(frame.Y,2) 
  ga=RA.Axes(j1,1);
  plot(ga,frame.X{1},frame.Y(:,j1));axis(ga,'tight');
  xlabel(ga,'Time [s]');
  ylabel(ga,sprintf('Channel %i',frame.X{2}(j1)));
 end
 % Compute current point spectra
 XF=fe_curve(sprintf('h1h2 %i -stack',refind),frame);
 Guu=stack_get(XF,'curve','Guu','get');
 ga=RA.Axes(1,2);
 plot(ga,Guu.X{1},Guu.Y);axis(ga,'tight');
 xlim(ga,[Guu.X{1}(2) Guu.X{1}(end)]);
 xlabel(ga,'Frequency [Hz]');
 ylabel(ga,sprintf('FFT Channel %i',refind));
 set(ga,'YScale','log');
 Gyy=stack_get(XF,'curve','Gyy','get');
 for j1=1:size(Gyy.Y,2)
  ga=RA.Axes(j1+1,2);
  plot(ga,Gyy.X{1},abs(Gyy.Y(:,j1)));axis(ga,'tight');
  xlim(ga,[Gyy.X{1}(2) Gyy.X{1}(end)]);
  xlabel(ga,'Frequency [Hz]');
  ylabel(ga,sprintf('FFT Channel %i',Gyy.X{2}(j1,1)));
  set(ga,'YScale','log');
 end
 % Compute estimators (all frames)
 XF=fe_curve(sprintf('h1h2 %i -stack',refind),frames);
 COH=stack_get(XF,'curve','COH','get');
 ga=RA.Axes(1,3);
 plot(ga,COH.X{1},COH.Y);
 xlim(ga,[COH.X{1}(2) COH.X{1}(end)]);ylim(ga,[0 1]);
 xlabel(ga,'Frequency [Hz]');
 ylabel(ga,'Coherence');
 H1=stack_get(XF,'curve','H1','get');
 for j1=1:size(H1.Y,2)
  ga=RA.Axes(j1+1,3);
  plot(ga,H1.X{1},abs(H1.Y(:,j1)));axis(ga,'tight');
  xlim(ga,[H1.X{1}(2) H1.X{1}(end)]);
  xlabel(ga,'Frequency [Hz]');
  ylabel(ga,sprintf('H1 Channel %i',H1.X{2}(j1,1)));
  set(ga,'YScale','log');
 end
 figure(get(ga,'Parent'));
 drawnow
end
