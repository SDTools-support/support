function   out=d_sim(varargin)

%% Simulink examples 

[CAM,Cam]=comstr(varargin{1},1);carg=2;

if comstr(Cam,'init')
    
  ref=evalin('base','ref');
  if isfield(ref,'AnSampleTime') % dt,T0 set in Blockparameters
   r1=ref.AnSampleTime; if ~ischar(r1);r1=comstr(r1,-30);end
   set_param('sim_anim/anim_struct','SystemSampleTime',r1)
  end
  if isfield(ref,'cf')
   cf=ref.cf; % Actually declare the proper feplot
   figure(ref.cf.opt(1)); % Raise animation so that it can be viewed
  else cf=[];
  end
   
elseif comstr(Cam,'view');[CAM,Cam]=comstr(CAM,5);
 if comstr(Cam,'anim')
%% #ViewAnim_initializes the figure
 ref=u; 
 if isempty(ref);ref=evalin('base','ref');
 elseif isa(ref,'ss') % Sys,TR provided
  sys=ref; TR=varargin{1};mdl=varargin{2};
  ref=struct('a',sys.a,'b',sys.b,'c',eye(size(sys.a,1)/2,size(sys.a,1)), ...
    'd',zeros(size(sys.a,1)/2,size(sys.b,2)),'TR',TR,'mdl',mdl);
 % Adjust the size of feedback forces in the MUX based on b
 end
 set_param('sim_anim/Mux','Inputs',sprintf('[1 %i]',size(ref.b,2)-1)) 
 fprintf('Adjusted Fnl to dimension %i\n',size(ref.b,2));
 
 if ~isfield(ref,'cf'); 
 % Open the feplot figure and initialize for animation
  ref.cf=feplot(ref.mdl);
 end
 if nargin==3
 % display simout
 simout=varargin{1};
 if isempty(ref.cf);ref.cf=feplot;end
 ref.cf.def=struct('def',simout.Data', ...
    'DOF',(1:size(ref.c,1))'+.99, ...
    'data',simout.Time, 'TR',ref.TR,'fun',[0 4]); %#ok<STRNU>
 else
  % Initialize with a static response
  ref.cf.def=struct('def',ref.c*(ref.a\sum(ref.b,2)), ...
    'DOF',(1:size(ref.c,1))'+.99, ...
    'data',zeros(1), 'TR',ref.TR,'fun',[0 4]);
  fecom colordataevalx
  if nargout>0; out=ref; end
 end
 else;error('View%s',CAM);

 end
elseif comstr(Cam,'testubeam');[CAM,Cam]=comstr(CAM,10);
%% TestUbeam ubeam example
 mdl=demosdt('demo ubeam mix');cf=feplot;mdl=cf.mdl.GetData;
 mdl=fe_case(mdl,'reset','FixDof','base','z==0', ...
    'DofLoad','Ext',struct('def',1,'DOF',314.01), ...
    'DofLoad','Rel',struct('def',[1 -1;1 0],'DOF',[244.01;114.01]), ...
    'SensDof','Out',[244.01;114.01]);
 % uniform 1 % modal damping
 mdl=stack_rm(mdl,'info','RayLeigh');
 mdl=stack_set(mdl,'info','DefaultZeta',.01);

 if comstr(Cam,'ss')
  %% #ubeamss : use d_sim testubeamSS
  open_system sweep_ss.slx
  mdl=fe_case(mdl,'remove','Rel');
  [sys,TR] = fe2ss('free 6 10',mdl);  

  ref=struct('a',sys.a,'b',sys.b,'c',sys.c,'d',sys.d, ...
      'AnSampleTime',[.005 0]);

  r1=sim('sweep_ss', ...
      'SaveOutput','on','OutputSaveName','yout','SaveFormat', 'Dataset');
  
  C1s=r1.logsout{1}.Values;
  C1=resample(C1s,0:ref.AnSampleTime(1):4);
  C1=struct('X',{{C1.Time,sys.OutputName}},'Xlab',{{'Time','out'}}, 'Y',C1.Data)

  iicom('curveinit','x',C1); 

  %set_param('sweep_ss/fext','WaveForm','sine','Frequency','60');

 else
  [sys,TR] = fe2ss('free 6 10',mdl);  
  ref=d_sim('viewanim',sys,TR,mdl);
  ref.fnl.c=sys.b(size(sys.a,1)/2+1:end,1)'; % Relative motion
  ref.fnl.fun=@(y)1000*y.^3;
  ref.AnSampleTime=[.02 0];sdt_sim_anim('view',ref);
 
  % remove some inputs
  ref.b=ref.b(:,1:3); ref.d=ref.d(:,1:3);
  % Define a non-linearity
  ref.fnl.c=reshape(ref.b(:,2),[],2)';ref.fnl.c(size(ref.b,2)-1,1)=0;
  % Set animation step
 
  ref.AnSampleTime=[.005 0]; % Set sampling for animation
  set_param('sim_anim/fext','WaveForm','sine','Frequency','60');
 
  a=sim('sim_anim','StopTime',2);
 end
 assignin('base','ref',ref)

else
%%
    error('%s unkown',CAM)
end
 if 1==2
    
%% In ModelExplorer (block animation)
% coder.extrinsic('feplot')
% init : in BlockSet InitFcn callback
% 
a=get_param('sim_anim/anim_struct','DialogParameters')
a=get_param('sim_anim/anim_struct','ObjectParameters');
st=fieldnames(a);
for j1=1:size(st,1);b=a.(st{j1});st{j1,2}=b.Type;st{j1,3}=b.Enum;end

a(~cellfun('isempty',regexpi(a,'time')))


a=get_param('sim_anim/anim_struct','CompiledSampleTime')
a=get_param('sim_anim/anim_struct','SystemSampleTime')
% set_param('sim_anim/anim_struct','SystemSampleTime','.01')

get_param('sim_anim/StateSpace DispObs','SystemSampleTime') % dt,T0
get_param('sim_anim/anim_struct','SystemSampleTime',[.02 0]) % dt,T0
set_param('sim_anim/anim_struct','SystemSampleTime',[.01 0]) % dt,T0

%codegen -report anim_eb -args {1,ones(4,1),1}
%dbstack;keyboard
%t=varargin{1};
%if nargin==3; u=varargin{2};Scale=varargin{3};end
    
 end
    
end


