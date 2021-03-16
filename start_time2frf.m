function [out,out1,out2]=start_time2frf(varargin) %#ok<STOUT>

% Time signal acquisition example
%
% start_time2frf('h1h2 example') % Uses state-space and multiple signals 
%    to be averaged
% sdtweb('start_time2frf.m#hammer') % Illustrates a hammer test
% sdtweb('start_time2frf.m#fe_time') % Illustrates fe_time integration
% 
% See also : d_signal 


%       Etienne Balmes, Guillaume Martin
%       Copyright (c) 1990-2021 by SDTools, All Rights Reserved.
%       start_time2frf('cvs') for revision information

%#ok<*ASGLU,*NOSEM,*NASGU>

[CAM,Cam]=comstr(varargin{1},1);carg=2;

if comstr(Cam,'h1h2 example')

%% #Tuto H1H2_example --------------------------------------------------------------
   
% Step 0: Generate sample model in state space form
%  to be used for the example
 RT=struct('BlockSize',4096,'FSamp',300, ...
     'source',{{fe_curve('testeval double(t<=.005)')}}, ... % Impact
     'FrameGen','TimeFreq',...
     'TimeOutNoise',d_signal('@NoiseTimeAdd'),'TimeOutNoiseLevel',.1); 
 [fgen,model]=d_signal('ModelGartFESys -zeta .001',RT); %
 
 % fgen : curvemodel for experiment simulation
 
%% Step 1 : get on block of time data (normally done with data acquisition system)

 % each time fgen.GetData is called a new frame is generated using the
 % noise model
 ci=iiplot;
 iicom('curveinit',{'curve','Frame1',fgen.GetData;'curve','Frame2',fgen.GetData});
 iicom('ch2');% show accel response

%% Step 2: estimate the FRFs from time data (request H1)

 RT=struct('MaxFrame',1,'Out',{{'H1','Coh'}},'OutName',{{'Test','Coherence'}},'ci',ci);
 fgen.Source.NoiseLevel=0; fgen('H1H2',RT);
 iicom(';submagpha;showtest;');


%% Step 3 : identify mode shapes (see gartid, sdtweb('idrc'))

po=ii_pof(eig(fgen.Source.Model.sys.a)/2/pi,3,1);
po(abs(po(:,1))<1e-1,:)=[]; % No rigid body mode

ci.Stack{'curve','IdMain'}=struct('po',po);
idcom('estlocalpole');

if 1==2
 cf=feplot;feplot('initmodel-back',model);cf.sel='-Test';
 cf.def=r3{1,3};fecom('ch86') % Response near 86 Hz
end

%% Step 4 : now generate multiple frames with noise model and iiplot

 RT=struct('MaxFrame',3,'Out',{{'H1','Coh'}},'OutName',{{'Test','Coherence'}},'ci',ci);
 fgen.Source.NoiseLevel=1e-3;RT.MaxFrame=20;fgen('H1H2',RT)
 ci.os_('PiSubHCoh')
 
%% Step 5: now add windowing

 RT=struct('MaxFrame',3,'Out',{{'H1','Coh'}},'OutName',{{'Test','Coherence'}},'ci',ci);
 RT.Window='exponential 0 10 10 0';
 fgen.Source.NoiseLevel=2e-3;RT.MaxFrame=20;fgen('H1H2',RT)

 % ci.os_('PiSubHCoh')
 
 
%% Step 6: now use a long time signal
 RT=struct('BlockSize',4096*8,'FSamp',300, ...
     'source',{{fe_curve('testeval double(t<=.005)')}}, ... % Impact
     'FrameGen','TimeFreq',...
     'TimeOutNoise',d_signal('@NoiseTimeAdd'),'TimeOutNoiseLevel',.1); 
 'xxxgm noise'
 [fgen,model]=d_signal('ModelGartFESys -zeta .001',RT); %
 TI=fgen.GetData;% Long time buffer
 TI.Edit=struct('buflen',1000,'overlap',.9);
 
 RT=struct('BufTime',1,'Overlap',90,'fun',d_signal('@H1H2'))
 % xxx adapt  for H1H2  sdtweb ii_signal spectro.setCheck
 XFm=feval(ii_signal('@spectro'),'init',TI,RT)

 RT=struct('MaxFrame',3,'Out',{{'H1','Coh'}},'OutName',{{'Test','Coherence'}},'ci',ci);
 fgen.Source.NoiseLevel=1e-3;RT.MaxFrame=20;fgen('H1H2',RT)
 
 
%% EndTuto
%% ---------------------------------------------------------------

%% #fe_time --------------------------------------------------------------
elseif comstr(Cam,'scriptfe_time')

 model=start_time2frf('DemoGartfeTime');
 model=stack_set(model,'info','Rayleigh', ...
   [10 0   1e-5  0.0; ... % Elements of group 10 (masses)
     9 0   0.0  1e-3; ... % Elements of group 9 (springs)
     0 1   0.0  1e-4;      ... % Elements with MatId 1
     0 2   0.0  1e-4]);       % Elements with MatId 2
 opt=fe_time('Newmark .25 .5 0 1e-4 10e3');opt.NeedUVA=[0 1 0];
 
 model=fe_case(model,'setcurve','Force','TestRicker dt=.01 A=2');
 def=fe_time(opt,model);def.def=def.v; % Keep velocity only

 Sens=fe_case(model,'sens','Sensors');
 C1=fe_case('SensObserve',Sens,def);
 cf.def=fe_def('subdef',def,1:20:size(def.def,2));
 
 iicom('curveInit','Vel',C1);
 ci=iiplot;ii_mmif('FFT fmax =100',ci,'Vel');iicom('iixonly','fft(Vel)');

    
%% ---------------------------------------------------------------
elseif comstr(Cam,'timegen');[CAM,Cam]=comstr(CAM,11); 
%% #TimeGen : generate a time buffer given source and state space
 cur=varargin{carg};carg=carg+1;
 out=d_signal('RespFrame',varargin{2:end});
 

%% #DemoGartFE : using the gartfe mesh as example ----------------------------
elseif comstr(Cam,'demogartfe');[CAM,Cam]=comstr(CAM,11);

  CAM=['ModelGartFE' CAM];
  fprintf('Should use d_signal(''%s'')',CAM);
  if nargout==2;[out,out1]=d_signal(CAM,varargin{2:end});
  else; out=d_signal(CAM,varargin{2:end});
  end
  
%% ---------------------------------------------------------------------------
elseif comstr(Cam,'cvs');
 out=sdtcheck('revision');
 %out='$Revision: 531 $  $Date: 2020-12-16 21:52:35 +0100 (Wed, 16 Dec 2020) $';
else; error('%s unknown',CAM);
end
