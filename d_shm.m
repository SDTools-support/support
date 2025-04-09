function [out,out1,out2]=d_shm(varargin); %#ok<*STOUT>

% D_SHM Parametric studies for SHM
%
%
% d_shm ScriptDocPiezo
%
% See <a href="matlab: sdtweb _taglist d_shm">TagList</a>
% Currently documented in PJE/tex/doc_beam_nl.tex
%

% Marc Rebillat, Mikhail Guskov, Eric Monteiro, ENSAM/PIMM
% Etienne Balmes, Jean-Philippe Bianchi, Guillaume Vermot des Roches, SDTools
%       Copyright (c) 1990-2025 by SDTools, All Rights Reserved.
%       For revision information use d_shm('cvs')

if nargin==0
  sdtweb _taglist d_shm % see structure of d_shm file
  return
end

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*DEFNU,*NOSEM>
if nargin<1; CAM=''; Cam='';
elseif ischar(varargin{1});
  obj=[];evt=[];[CAM,Cam]=comstr(varargin{1},1); carg=2;
else;
  obj=varargin{1};evt=varargin{2};[CAM,Cam]=comstr(varargin{3},1); carg=4;
  if isstruct(evt) && isfield(evt,'val')
    RO=fe_range('ToFlat',obj,evt); % get the input of command from obj and event
  else
  end
end

%% #Script -------------------------------------------------------------------
if comstr(Cam,'script');[CAM,Cam]=comstr(CAM,7);
  
  
  % Generic opening of project for all demos here
  RA=struct('ProjectWd',fullfile(sdtdef('TempDir'),'d_shm'),'root','base');
  sdtroot('setproject',RA);UI=sdtroot('ParamUI');
  
  % base script for integrated parametric studies of piezo transient simulations
  if comstr(Cam,'docpiezo');[CAM,Cam]=comstr(CAM,9);
    %% #ScriptDocPiezo
    
    
    wd=fullfile(sdtdef('TempDir'),'d_shm');sdtkey('mkdir',wd);
    RA=struct('ProjectWd',wd,'root','base');
    sdtroot('setproject',RA);UI=sdtroot('ParamUI');
    R1=struct('MeshCfg','docpiezo','SimuCfg','docpiezo','RunCfg','docpiezo');
    Range=d_shm('RangeCheck',R1);sdtroot('setrange',Range)
    Range=d_shm('RangeMeshCfg',Range,'docpiezo');
    Range.param.RunCfg.data{Range.val(1)}.save_def=1;
    
    R1=fe_range('loop',Range); % Perform the simulations
    cd(sdtroot('PARAM.Project.ProjectWd'));fe_range('dirscan',struct('wd',pwd));
    PA=sdtroot('paramvh')
    
    d_shm('CbMovie',z.model,z.def);
    
  elseif comstr(Cam,'nida13');[CAM,Cam]=comstr(CAM,7);
    %% #ScriptNida13
    % sdtweb nida13 scriptStiffener
    
    %% #Script M_parameters
    Range=d_shm('RangeCheck',struct('MeshCfg','PimmShell', ...
      'M_Def1_xc',[30 35],'AlphaR',0,'SimuCfg','PIMM', ...
      'RunCfg','docpiezo'));
    sdtroot('SetRange',Range);
    
    fe_range('loop',Range); % Perform the simulations
    
    
    
    %% #SVS shell volume shell config -2
    R1=struct('MeshCfg','PimmSVS','SimuCfg','PIMM','RunCfg','docpiezo');
    Range=d_shm('RangeCheck',R1); fe_range('loop',Range);
    sdtroot('SetRange',Range);
    
    
    if 1==2 % XXXeb
      
      Range=d_shm('RangeCheck',struct('MeshCfg','PimmShell', ...
        'M_Def1_xc',[30 35],'AlphaR',0,'SimuCfg','PIMM','RunCfg','docpiezo'));
      sdtroot('SetRange',Range);
      
      
      Range=d_shm('RangeCheck',struct('MeshCfg','PimmShellA',...
        'M_SVS_thick',25));
      fe_range('loop',Range); % Perform the simulations
      
      R1=struct('MeshCfg',{{'PimmShellA','PimmSVSA'}},'dt',[2.0000e-06; 3.0000e-06]);
      R1=fe_range('valSet',Range,R1);
      regexp('M_Act1_xc','M_(.*)_([^_]+)','tokens');ans{1} %#ok<NOANS>
      regexp('M_Act1_v2_xc','M_(.*)_([^_]+)','tokens');ans{1} %#ok<NOANS>
    end
    
    
    % Some other call examples: - - -
    % - this is kept but is not still working : to be implemented
    Range=d_shm('RangeCheck',struct('MeshCfg','PimmShellA',...
      'M_curve1_FullRow', ...
      {{'curve1','MeshDefCurve',struct('x0',[0.2 0.15 -.1],'dir',[1 0 0],'r',.1)}}));
    %  - ...
    if 1==1  % Position of actuator 1
      Range=d_shm('RangeM_Act1_xc',Range,[0.34 0.36]);
    else % single run
      Range=d_shm('RangeMeshCfg',Range,1);
    end
    
    %  define a base range with configuration values in step order
    %Range=struct('val',[1],'lab','MeshCfg');
    %Range=nida13('TimeRange',Range);
    
    %Range.val=1; Range.lab={'MeshCfg'}; Range.FileName={'@MeshCfg'};
    if nargout>0; out=Range;end
    
    %%
  elseif ~isempty(CAM); error('Script%s unknown',CAM);
  end
  
  %% #Study : some generic parametric studies ----------------------------------
elseif comstr(Cam,'study');[CAM,Cam]=comstr(CAM,6);
  
  %% #StudyRayleigh : some generic parametric studies --------------------------
  if comstr(Cam,'rayleigh');[CAM,Cam]=comstr(CAM,9);
    
    Range=d_shm('RangeCheck',RO.Range); % Build Range grid
    Range.FileName{1}=RO.root;
    sdtroot('SetRange',Range);
    sdtroot('SetProject',struct('ProjectWd',RO.wd));
    e_range('loop',Range); % Perform the simulations
    
    
  end
  %% #MeshCfg Step10 --------------------------------------------------------
  % should be moved to d_mesh('MeshCfg')
elseif comstr(Cam,'meshcfg');[CAM,Cam]=comstr(CAM,8);
  
  % take base range input
  if carg<=nargin; Range=varargin{carg};carg=carg+1;
  else; Range=d_shm('Range');end
  if isempty(CAM) && carg<=nargin;
    r1=varargin{carg}; carg=carg+1;
    [CAM,Cam]=comstr(r1,1);
  end
  MeshCfg={};
  % Step10 : mesh generation
  if ~isempty(Cam)
  elseif iscell(Range)&&size(Range,2)==2;
    MeshCfg=Range;Range=struct('param',struct);
  elseif carg>nargin
    % Default Mesh configurations
    out=d_shm('MeshCfg',Range,'pimmshella','pimmshellb','pimmshellacurve',...
      'pimmsvsa','pimmOmega','docpiezo');
    return;
  else;[CAM,Cam]=comstr(varargin{carg},1);carg=carg+1;
  end
  
  % build general fields :
  if ~isfield(Range,'param'); Range.param=struct; end
  if ~isfield(Range,'FileName'); Range.FileName={}; end
  if ~isfield(Range.param,'MeshCfg')
    Range.param.MeshCfg=struct('type','pop','value',1,'level',10,...
      'choices',{{}},'data',{{}},'SetFcn',{{@d_shm,'stepMesh'}},...
      'ShortFmt',1, ...
      'RepList', ... % 'TestEvtData',1,
      {{'M_(.*)_([^_]+)', ... % start with M_ and end with non_ parameter name
      struct('type',{'.','{}'},'subs',{'list',{'@token{1}',3}})
      }});
  end
  if exist('nida13','file');RunOpt.MatFun=@nida13;
  else; RunOpt.MatFun=@d_shm;
  end
  
  while 1 % loop on asked configs
    
    if comstr(Cam,'pimmshell');[CAM,Cam]=comstr(CAM,10);
      %% #MeshPIMMShell A|B --------------------------------------------------------2
      
      [RB,st,CAM]=cingui('paramedit -DoClean',[ ...
        'name("PIMMShell"#%s#"default name") ' ...
        'matlam("BmiNom"#%s#"laminated material name")' ...
        'curve(#3#"curve mode")' ...
        ],{struct,CAM});Cam=lower(CAM);
      RunOpt=sdth.sfield('AddMissing',RunOpt,RB);
      if RB.curve; RunOpt.curve=RB.curve; end
      
      RunOpt.Typ=Cam;
      
      mdl=struct('Node',[],'Elt',[],'pl',[],'il',[],'unit','MM');
      [mdl.pl,mdl.il,r2]=feval(RunOpt.MatFun,['matlam' RunOpt.matlam]);
      RN=struct; % Options for the current mesh config
      switch RunOpt.Typ
        case 'a' % #MeshPIMMShellA -3
          RN.list={'Name','Lam','shape'
            'Main_plate', mdl,struct('type','global','lx',400,'ly',300,'lc',8)
            'Pz1','BaseId1 +Noliac.NCE51.OD0TH1.in', ...
            struct('shape','circ','xc',27,'yc',75,'rc',10)
            'Pz2','BaseId1 +Noliac.NCE51.OD0TH1', ...
            struct('shape','circ','xc',200,'yc',75,'rc',10)
            'Pz3','BaseId1 +Noliac.NCE51.OD0TH1', ...
            struct('shape','circ','xc',373,'yc',75,'rc',10)
            'Pz4','BaseId1 +Noliac.NCE51.OD0TH1', ...
            struct('shape','circ','xc',113.5,'yc',225,'rc',10)
            'Pz5','BaseId1 +Noliac.NCE51.OD0TH1', ...
            struct('shape','circ','xc',286.5,'yc',225,'rc',10)
            'Def1','BaseId1', ...
            struct('shape','circ','xc',300,'yc',150,'rc',7)
            };
        case 'b' % #MeshPIMMShellB -3
          % Meriem nouvelle plaque
          RN.list={'Name','Lam','shape'
            'Main_plate', mdl,struct('type','global','lx',300,'ly',100,'lc',8)
            'Pz1','BaseId1 +Noliac.NCE51.OD0TH1.in', ...
            struct('shape','circ','xc',50,'yc',50,'rc',10)
            'Pz2','BaseId1 +Noliac.NCE51.OD0TH1', ...
            struct('shape','circ','xc',250,'yc',50,'rc',10)
            'Def1','BaseId1', ...
            struct('shape','circ','xc',150,'yc',50,'rc',7)
            };
        otherwise; error('PIMMShell %s',RunOpt.Typ);
      end
      RN.wv=7000; % approximate wave velocity
      RN.FeplotCf=2;
      if RunOpt.curve
        %% #MeshPIMMShellACurve -3
        RN.list(end+1,1:3)={'curve1','MeshDefCurve',...
          struct('x0',[200 150 -100],'dir',[1 0 0],'r',100)};
        RunOpt.name=[RunOpt.name '+C'];
      end
      RN.config='Epoxy'; % for legends
      MeshCfg={RunOpt.name RN};
      
      
      
      
    elseif comstr(Cam,'pimmomega');[CAM,Cam]=comstr(CAM,8);
      %% #MeshPIMMOmega : Monolitic with stiffener -------------------------------2
      
      % Coarse mesh with an Omega
      mdl.Node=[1 0 0 0   0  -50-22-78 0;
        2 0 0 0   0  -50-22    0;
        3 0 0 0   0  -50       0;
        4 0 0 0   0  +50       0;
        5 0 0 0   0  +50+22    0;
        6 0 0 0   0  +50+22+78 0;
        7 0 0 0   0  -37.5     61;
        8 0 0 0   0  +37.5     61];
      mdl.Elt=feutil('Objectbeamline 1 2 3 4 5 6 0 3 7 8 4',mdl.Node);
      mdl=feutil('Extrude 1 400 0 0',mdl);
      mdl.pl=[]; mdl.il=[]; mdl.unit='MM';
      try
        [mdl.pl,mdl.il,r2]=nida13('matlam BmiNom');
        mdl.Node(:,6)=mdl.Node(:,6)-min(mdl.Node(:,6));
        RN=d_shm('MeshCfg',[],'PimmShellA');RN=RN.param.MeshCfg.data{1};
        
        % Base plate but accept a coarse model of the base plate
        RN.list(2,:)={'Main_plate','mdl',struct('type','global','lc',5,'mdl',mdl)};
        MeshCfg={'PIMMOmega' RN};
      catch
        MeshCfg={'PIMMOmega' []};
      end
      
    elseif comstr(Cam,'pimmsvsa');[CAM,Cam]=comstr(CAM,9);
      %% #MeshPIMMSVS A ----------------------------------------------------------2
      
      RN=d_shm('MeshCfg',[],'PimmShellA');RN=RN.param.MeshCfg.data{1};
      [pl,il]=RunOpt.MatFun('matlam Bmi_CoreA');
      if size(pl,1)==2
        pl(:,1)=[1;10];
        RN.list{2,2}.pl=pl; % BMI core
        RN.list(end+1,1:3)={'SVS','', ...
          struct('height',25,'plcore',pl(2,:))};
        RN.wv=4500; % approximate wave velocity
        RN.config='Sandwich'; % for legends
        MeshCfg={'PimmSVSA' RN};
      else; MeshCfg={'PimmSVSA',[]};
      end
    elseif comstr(Cam,'docpiezo');[CAM,Cam]=comstr(CAM,9);
      %% #Meshdocpiezo -----------------------------------------------------------2
      
      % Start by defining properties of the underlying laminate
      mdl=struct('Node',[],'Elt',[], ... % empty model
        'pl', ... % composite layer property
        [1 fe_mat('m_elastic','SI',1) 42.5e9 .042 1490 3.35e9 .01], ...
        'il', ... % laminate definition (6 layers at 0,90,0,90,0,90)
        p_shell(['dbval 1 laminate    1 2.167e-4 0    1 2.167e-4 90 ' ...
        '1 2.167e-4 0 1 2.167e-4 90 1 2.167e-4 0 1 2.167e-4 90']), ...
        'unit','SI');
      RO=struct;
      RO.list={'Name','Lam','shape'
        'Main_plate', mdl,struct('global',1,'lx',.4,'ly',.3,'lc',.02);
        'Act1','BaseId1 +SmartM.MFC-P1.2814 -SmartM.MFC-P1.2814.in',struct('xc',.35,'yc',.25,'ang',30);
        'Sen2','BaseId1 +SmartM.MFC-P1.2814',struct('xc',.03,'yc',.05,'ang',30);
        'Sen3','BaseId1 +Noliac.NCE51.OD25TH1',struct('xc',.05,'yc',.25);
        };
      RO.FeplotCf=2; % store in feplot
      MeshCfg={'DocPiezo' RO};
      
    elseif comstr(Cam,'@');MeshCfg=eval(CAM(2:end));
      out1=struct('NeedInit',{MeshCfg(:,1)});
      
    elseif ~isempty(MeshCfg)
    else; error('Mesh%s unknown',CAM);
      %% #MeshEnd ----------------------------------------------------------------2
    end
    
    MeshCfg{2}.MeshCfg=MeshCfg{1};
    Range.param.MeshCfg=feval(fe_range('@popMerge'),Range,'MeshCfg',MeshCfg);
    if nargout>1; out1=struct('NeedInit',{MeshCfg});end
    if carg>nargin; break;end %Possibly multiple matches
    [CAM,Cam]=comstr(varargin{carg},1);carg=carg+1;
  end % loop on meshconfig
  
  st1='@MeshCfg';
  if ~any(strcmpi(Range.FileName,st1));Range.FileName{end+1}=st1;end
  
  out=Range;
  
  %% #SimuCfg : time integration and damping -----------------------------------
elseif comstr(Cam,'simucfg');[CAM,Cam]=comstr(CAM,8);
  
  % take base range input
  if carg<=nargin; Range=varargin{carg};carg=carg+1;
  else; Range=d_shm('Range');end
  if isempty(CAM) && carg<=nargin;
    r1=varargin{carg}; carg=carg+1;
    [CAM,Cam]=comstr(r1,1);
  end
  out1=[];
  % Step 20 : simulation configuration (TimeOpt)
  if isempty(Cam) % Default Mesh configurations
    out=d_shm('SimuCfg',Range,'Base','docpiezo');
    return;
  end
  
  % build general fields :
  if ~isfield(Range,'param'); Range.param=struct; end
  if ~isfield(Range.param,'SimuCfg')
    Range.param.SimuCfg=struct('type','pop','value',1,'level',20,...
      'choices',{{}},'data',{{}},'SetFcn',{{@nida13,'TimeOpt'}},...
      'ShortFmt',1);
  end
  
  while 1 % loop on asked configs
    
    if comstr(Cam,'base');[CAM,Cam]=comstr(CAM,5);
      %% #SimuBase : basic for PIMMa --------------------------------------------2
      RN=struct('dt',.3e-6,'tend',.4e-3,'Tcoef',1,...
        'AlphaR',0,'BetaR',2*.0025/200e3,'FeplotCf',2); %,'c_u',1
      %RN.AssembleCall='assemble -fetime Load'; % XXXeb ?
      %RN.FinalCleanupFcn='out.DOF=model.DOF;out.def=[Case.T*out.def+Case.TIn*[0 ft(1:end-1)'']];'; % XXXeb ?
      
      SimuCfg={'Base',RN};
    elseif comstr(Cam,'pimm');[CAM,Cam]=comstr(CAM,5);
      %% #SimuPIMM : basic for PIMMa --------------------------------------------2
      
      [RN,st,CAM]=cingui('paramedit -DoClean',[ ...
        'dt(.3e-6#%g#"time step") ' ...
        'Tend(.4e-3#%g#"end of simulation") ' ...
        'inf0(200#%g#"Input base frequency")' ...
        'Tcoef(1#%g#"xxx")' ...
        'AlphaR(0#%g#"Rayleigh coefficient on mass")' ...
        'BetaR(2.5e-8#%g#"Rayleigh coefficient on stiffness 2*.0025/200e3")' ...
        'BetaRVar(#%s#"callback to define a varying Rayleigh")' ...
        'zeta(0.0007#%g#"target damping ratio")' ...
        'FeplotCf(2#%g#"number of feplot figure")' ...
        'name(PIMM#%s#"configuration name")' ...
        'OSamp(0#%g#"output sampling frequency, 0=all")' ...
        'TVec(#%s#"function handle in string format to generate a non standard TimeVector")'...
        'tItSplit(10#%i#"number of time intervals to split dt, use with a TVec callback")' ...
        'Method("Newmark"#%s#"time scheme method")' ...
        ],{struct,CAM});Cam=lower(CAM);
      
      %RN.AssembleCall='assemble -fetime Load'; % XXXeb ?
      %RN.FinalCleanupFcn='out.DOF=model.DOF;out.def=[Case.T*out.def+Case.TIn*[0 ft(1:end-1)'']];'; % XXXeb ?
      
      SimuCfg={RN.name,RN};
      
    elseif comstr(Cam,'docpiezo');[CAM,Cam]=comstr(CAM,9);
      %% #SimuDocPiezo -----------------------------------------------------------2
      SimuCfg={'docpiezo',struct('dt',.3e-6,'NStep',200,...
        'AlphaR',0,'BetaR',2*.0025/200e3,'FeplotCf',2,...
        'AssembleCall','assemble -fetime Load')};
      Range.param.SimuCfg.SetFcn={@d_shm,'StepSimu'};
      
    else; error('SimuCfg%s unknown',CAM);
      %% #SimuEnd ----------------------------------------------------------------2
    end
    
    if isfield(Range,'SimuCfg')&&ischar(Range.SimuCfg);Range.SimuCfg=SimuCfg{1};end
    Range.param.SimuCfg=feval(fe_range('@popMerge'),Range,'SimuCfg',SimuCfg);
    
    if carg>nargin; break;end %Possibly multiple matches
    [CAM,Cam]=comstr(varargin{carg},1);carg=carg+1;
  end % loop on meshconfig
  if nargout>1; out1=struct('NeedInit',{SimuCfg});end
  
  st1='@SimuCfg';
  if ~any(strcmpi(Range.FileName,st1));Range.FileName{end+1}=st1;end
  
  out=Range;
  
  
  %% #RunCfg : perform simulation and post -------------------------------------
  % Step 30 : perform simulation and post
elseif comstr(Cam,'runcfg');[CAM,Cam]=comstr(CAM,7);
  
  if carg<=nargin; Range=varargin{carg};carg=carg+1;
  else; Range=d_shm('Range');end
  if isempty(CAM) && carg<=nargin;
    r1=varargin{carg}; carg=carg+1;  [CAM,Cam]=comstr(r1,1);
  end
  
  if isempty(Cam) % Default Mesh configurations
    out=d_shm('RunCfg',Range,'TransElec','docpiezo','docepiezo'); return;
  end
  
  % build general fields :
  if ~isfield(Range,'param'); Range.param=struct; end
  if ~isfield(Range.param,'RunCfg')
    Range.param.RunCfg=struct('type','pop','value',1,'level',30,...
      'choices',{{}},'data',{{}},'SetFcn',{{@nida13,'TimeElec-Run'}},...
      'ShortFmt',1);
  end
  
  while 1 % loop on asked configs
    
    if comstr(Cam,'transelec');[CAM,Cam]=comstr(CAM,10);
      %% #RunTransElec -----------------------------------------------------------2
      RunCfg={'TransElec',struct('cf',201,'save',1,'FeplotCf',[])};
      
    elseif comstr(Cam,'docepiezo');[CAM,Cam]=comstr(CAM,9);
      %% #Rundocepiezo : just elec ---------------------------------------2
      RunCfg={'docepiezo',struct('cf',201,'save',1,'OutElec',1,'FeplotCf',[])};
      Range.param.RunCfg.SetFcn={@d_shm,'StepRun'};
      
    elseif comstr(Cam,'docpiezo');[CAM,Cam]=comstr(CAM,9);
      %% #Rundocpiezo  -----------------------------------------------------------2
      RunCfg={'docpiezo',struct('cf',201,'save',1,'FeplotCf',2)};
      Range.param.RunCfg.SetFcn={@d_shm,'StepRun'};
      
    else; error('Run%s unknown',CAM);
      %% #RunEnd -----------------------------------------------------------------2
    end
    
    Range.param.RunCfg=feval(fe_range('@popMerge'),Range,'RunCfg',RunCfg);
    
    if carg>nargin; break;end %Possibly multiple matches
    [CAM,Cam]=comstr(varargin{carg},1);carg=carg+1;
  end % loop on RunCfg
  
  st1='@RunCfg';
  if ~any(strcmpi(Range.FileName,st1));Range.FileName{end+1}=st1;end
  out=Range;
  
  %% #MatLam : sample materials --------------------------------------------------
elseif comstr(Cam,'matlam');[CAM,Cam]=comstr(CAM,7);
  
  sdtw('_nb','Default d_shm laminate');
  out=[ 1 fe_mat('m_elastic','SI',6) 60e9 60e9 8e9 0.02 0.02 0.02 5e9 5e9 5e9 1550 ];
  out1=[ 1 fe_mat('p_shell','MM',2) -0.5 0 0 0 0 0 0  ...
    1 0.25 0 0 1 0.25 -45 0 1 0.25 45 0 1 0.25 0 0 ];
  out2=[];
  
  %% #Sens : generate sensors --------------------------------------------------
elseif comstr(Cam,'sens');[CAM,Cam]=comstr(CAM,5);
  
  model=varargin{carg};carg=carg+1;
  
  if comstr(Cam,'circle')
    Test=feutil('objectcircle .4 .4 0 .3 0 0 1 36');
    cut=fe_caseg('StressCut -selout -radius 3',Test,model);
    model=stack_set(model,'info','Circle',cut);
  end
  
  out=model;
  
  
  %% #Step -------------------------------------------------------------------
elseif comstr(Cam,'step');[CAM,Cam]=comstr(CAM,5);
  
  
  if ~isempty(evt)
    % Generic handling of step verifications, typically from fe_range('loop')
    Range=obj;
    if comstr(lower(evt.lab),'mesh')
      % M_ parameters in RO.list{1,3} with name based indexing
      RB=struct('NeedField',{{'MeshCmd',{@d_piezo,'MeshPlate'}}});
    else; RB=struct('NeedField',{{}});
    end
    
    RO=fe_range('ValEvtMerge',Range,evt); % Replace params of current level
    for j1=1:size(RB.NeedField,1) % Possibly check existence of fields
      if ~isfield(RO,RB.NeedField{j1,1});
        RO.(RB.NeedField{j1,1})=RB.NeedField{j1,2};
      end
    end
    % Low level input
  elseif carg<=nargin; RO=varargin{carg};carg=carg+1;
  else; RO=struct;
  end
  
  %% #StepMesh (level 10 uses MeshCfg) - - - - - - - - - - - - - - - - - - -
  if comstr(Cam,'mesh');[CAM,Cam]=comstr(CAM,5);
    
    % really mesh the plate
    RO.fname=fe_range('Fname -level 30',Range);
    
    out=feval(RO.MeshCmd{:},RO); % d_piezo('MeshPlate',RO);
    out=stack_set(out,'info','MeshInfo',RO); % store mesh RunOpt for flat param
    if isfield(RO,'FeplotCf') && ~isempty(RO.FeplotCf)
      cf=comgui('guifeplot -reset -project "SDT Root"',RO.FeplotCf);
      cf.model=out;
    else;RO.KeepRes=1;
    end
    
    if isfield(RO,'KeepRes') && ~isempty(RO.KeepRes)
      assignin('caller','Res',out);
    end
    
    %% #StepSimu : initialize data for simulation - - - - - - - - - - - - - - - -
    % generalized from nida13('timeOpt')
  elseif comstr(Cam,'simu');[CAM,Cam]=comstr(CAM,5);
    
    if ~isfield(RO,'FeplotCf') || isempty(RO.FeplotCf)
      error('This case is not handled. Must provide FeplotCf')
    else
      cf=clean_get_uf('feplotcf',RO.FeplotCf);
      model=cf.mdl;
    end
    % Build time opt:
    if isfield(RO,'Tcoef') && ~isempty(RO.Tcoef); RO.dt=RO.dt/RO.Tcoef;end
    if ~isfield(RO,'NStep')
      if isfield(RO,'tend') && ~isempty(RO.tend); RO.NStep=ceil(RO.tend/RO.dt);
      else; error('You must give NStep or Tend')
      end
    end
    %check method
    if any(cellfun(@(x)isfield(x,'NLdata'),model.Stack(:,3))); RO.Method='NLNewmark';end
    
    if isfield(RO,'Method')
      if strcmp(RO.Method,'NLNewmark')
        sdtkey('cvsnum>1.337;','fe_time')
        [opt,model]=d_fetime('TimeOpt NLNewmark',RO,model);
        if sdtkey('cvsnum','nl_spring')>=1.509
          opt.chandle='model.cv';
        end      
        % Rayleigh is better handled in opt.Rayleigh, then need to remove info from +c
        if isfield(RO,'BetaRVar')&&~isempty(RO.BetaRVar)
            warning('report example to Etienne');
          opt=feval(eval(RO.BetaRVar),opt,RO); % xxx should be done elsewhere
        end
      else
        warning('Does not support NL');
        opt=d_fetime(sprintf('TimeOpt Newmark .25 .5 0 %.15g %i',RO.dt,RO.NStep));
        % Rayleigh damping is in C matrix, handled by linear Newmark residual as c
      end
    end
    if isfield(RO,'TVec')&&~isempty(RO.TVec)
      if strncmpi(opt.Method,'nl',2); opt.JacobianUpdate=1; end % refact at dt change
      if isnumeric(RO.TVec); opt.TimeVector=RO.TVec;
      else % callback handling
        tvfcn=eval(RO.TVec);
        r1=tvfcn(RO);
        if isstruct(r1); opt.TimeVector=r1.TimeVector(:)'; else; opt.TimeVector=r1(:)'; end
      end
    end
    if isfield(RO,'OSamp')&&RO.OSamp;
      tout=0:1/RO.OSamp:RO.tend;
      opt.FinalCleanupFcn='d_shm(''CbFinalCleanup'');';
    else
      tout=[];
      opt.FinalCleanupFcn='d_shm(''CbFinalCleanup'');';
    end
    if any(cellfun(@(x)isfield(x,'NLdata'),model.Stack(:,3)))
      opt.AssembleCall=nl_spring('AssembleCall');
      opt.OutputFcn='of_time(''interp'',out,beta,gamma,Case.uva,a,tc-dt,tc,model.FNL)';
    else; opt.AssembleCall='assemble -fetime Load';
      opt.OutputFcn='of_time(''interp'',out,beta,gamma,Case.uva,a,tc-dt,tc)';
    end
    %opt.Rayleigh=[RO.AlphaR RO.BetaR]; % xxx handling depends on method
    if ~isempty(tout); opt.OutputFcn=tout; end
    if isfield(RO,'c_u')
      opt.c_u=RO.c_u; error('Report case should be handled by d_fetime');% xxx
    end
    opt=sdth.sfield('addselected',opt,RO,{'AssembleCall'});
    % defines piezo electrodes :
    r1=stack_get(model,'info','Electrodes','getdata');
    if isfield(RO,'Actuator') && ~isempty(RO.Actuator) % Defined as patch
      r1.data(:,2)=1; % OUT
      r1.data(RO.Actuator,2)=0; % IN
      model=stack_set(model,'info','Electrodes',r1);
    else; RO.Actuator=find(r1.data(:,2)==0);
    end
    p_piezo('electrodeinfo',model)
    model=p_piezo('electrode2case',model);  
    %add all PZT
    sens=struct('cta',speye(size(r1.data,1)),'DOF',r1.data(:,1)+.21, ...
        'lab',{cellfun(@(x) sprintf('pz%i',x),...
        num2cell((1:size(r1.data,1))'),'uni',0)}); 
    %  eval(iigui({'model','opt','sens'},'SetInBaseC')) % Allow debug
    model=fe_case(model,'SensDof','V_OUT',sens);
        
    %% Definition time variation of input
    if isfield(RO,'input'); %% Illustrated in lize16
      if ischar(RO.input);com=RO.input;
      elseif isfield(RO.input,'type');com=fe_curve(['test' RO.input.type],RO.input);com.name=RO.input.name;
      else;com=RO.input; %directly a curve
      end
    else
      if isfield(RO,'inA'); inA=RO.inA; else; inA=1; end
      if isfield(RO,'inn0'); inn0=RO.inn0; else; inn0=5; end
      if isfield(RO,'inf0'); inf0=RO.inf0; else; inf0=200e3; end
      com=sprintf('test coshan f0=%.15g n0=%.15g A=%.15g',inf0,inn0,inA);
    end
    if ischar(com);fprintf('In signal : %s\n',com)
    else;fprintf('In signal : %s\n',com.name)
    end
    if isempty(fe_case(model,'getdata','V_In'))
      r1=fe_case(model,'stack_get');
      r1=r1(ismember(lower(r1(:,1)),{'dofload','fsurf'}),:);
      model=fe_case(model,'setcurve',r1{1,2},'input',com);
    else
      model=fe_case(model,'setcurve','V_IN','input',com);
    end    
    %% Restart
    if isfield(RO,'restart')&&RO.restart>=0
     try;t=0:RO.dt:RO.tend;catch;keyboard;end
     C1=stack_get(model,'','input',3);ft=fe_curve('returny',C1,t);
     i1=find(ft~=0,1,'last')+5;
     if i1<length(ft);
       i2=find(t>=t(i1)+RO.restart,1,'first'); %u=du=d2u=d3u=0
       if isempty(i2);i2=i1+20;end;
       opt.SaveTimes=[t(i2) i2];
       wd=sdtroot('PARAM.Project.ProjectWd -safe');
       opt.SaveFile=fullfile(wd,'fe_time_save.mat');
     end;
    end    
    %% Save
    %r1=stack_get(model,'info','Electrodes','get');
    if isfield(RO,'OutElec')&&RO.OutElec
        opt.OutputInit='d_shm(''CbInitElecOutput'')';
    end

    if RO.AlphaR||RO.BetaR&&~isfield(opt,'Rayleigh');error('d_fetime problem');end
    model=stack_set(model,'info','TimeOpt',opt);
    model=stack_set(model,'info','SimuInfo',RO);% model is v_handle
    cf.mdl=model; %EM
    
    
    %% #StepRun - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  elseif comstr(Cam,'run');[CAM,Cam]=comstr(CAM,4);
    
    if ~isfield(RO,'FeplotCf') || isempty(RO.FeplotCf)
      r2=obj.param.SimuCfg.data{evt.val(strcmpi(obj.lab,'SimuCfg'))};
      RO.FeplotCf=r2.FeplotCf;
    end
    cf=clean_get_uf('feplotcf',RO.FeplotCf);model=cf.mdl.GetData;RO.cf=cf; 
    Range.FileName(strcmpi(Range.FileName,'@wdtest'))=[];
    RO.fname=fe_range('Fname -level 30',Range);
    if ~isfield(RO,'RunFcn'); RO.RunFcn={@DefaultRunFcn};end
    
    % RunFcn Method selected externally
    feval(RO.RunFcn{:});       
    % Store run step options
    model=stack_set(model,'info','RunInfo',RO);
    
    Res={def,C1,model};
    if ~isempty(cf)&&~isempty(fe_c(def.DOF,.21,'ind',2))
          cf.def=def; % xxx should be view
          fecom('ShowFiCevalA');fecom('scc1e-10');
    end
    % if isfield(RO,'ci')&&~isempty(RO.ci);
    C1=stack_set(C1,'info','Range',obj);% Save the Design point 
    ci=comgui('guiiiplot -project "SDT Root"');
    iicom(ci,'curveinit',C1.name,C1);
    try; iicom('chall');end

    if isfield(RO,'save') % #StepRunSave - - -3
      if ishandle(model); model=model.GetData; end
      RB=flatParam(struct('model',model));
      RB.Range=feutil('Rmfield',Range,'Res');
      % save list
      reslist={'C1';'model';'RO'};
      if isfield(RO,'save_def') && RO.save_def==1
        reslist=cat(1,reslist,{'def'});
      end
      wd=sdtroot('PARAM.Project.ProjectWd -safe');
      fname=fullfile(wd,RO.fname{1}); fname=[fname '_RES.mat'];
      RO=RB; save(fname,reslist{:})
      fprintf('Saved %s\n',fname);
    end
    assignin('caller','Res',Res); % Place result in range
    
    %% #StepEnd  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else
    error('Unknown Step%s',Cam)
  end
  
  
  %% #Range range params ---------------------------------
elseif comstr(Cam,'range'); [CAM,Cam]=comstr(CAM,6);
  
  Range=[];
  if carg<=nargin;val=varargin{carg};carg=carg+1; else;val=[];end
  if isfield(val,'param'); Range=val;
    if carg<=nargin; val=varargin{carg};carg=carg+1;end
  end
  if isempty(Range)&&carg<=nargin&&isfield(varargin{carg},'param')
    Range=varargin{carg};carg=carg+1;
  end
  if isempty(Cam)&&ischar(val);
    [CAM,Cam]=comstr(val,1);val=varargin{carg};carg=carg+1;
  end
  
  st1={'dt','s';'AlphaR','';'BetaR','';'Actuator','';'wdtest','';
    'inf0','Hz';'inn0','';'inA','V';'wv','mm/ms';'input','';
    'tend','s';'NStep',''}; % TimeOpt.Opt parameters
  i1=strcmpi(Cam,st1(:,1));
  if any(i1) % Simu parameters
    CAM=st1{i1,1};  unit=st1{i1,2};
    %% #Rangedt SimuCfg setting step 20 -2
    % d_shm('Rangedt',[1e-3 1e-4],Range)
    % d_shm('Rangeinput',{'I1',struct},Range);
    if iscell(val)
      out=struct('val',[],'lab',{{CAM}}, ...
        'param',struct(CAM,struct('type','pop','level',20)));
      if all(cellfun(@ischar,val(:))) || ...
          (size(val,2)==2&&all(cellfun(@ischar,val(:,1))))
        [out.param.(CAM),out.val]=feval(fe_range('@popMerge'),out,CAM,val);
      else; error('Not a valid val');
      end
    else
      out=struct('val',val(:),'lab',{{CAM}}, ...
        'param',struct(CAM, ...
        struct('type','double','level',20,...
        'Xlab',sprintf('%s [%s]',CAM,unit))),...
        'LabFcn',sprintf('sprintf(''%%.0f %s'',val);',unit));
    end
    
  elseif isempty(Cam)&&isfield(val,'NeedInit')
    %% #NeedInit : actually fill -2
    [out,out1]=d_shm(val.type,Range,val.NeedInit);
    % out1 should contain the selected values
    return;
    
  elseif ~isempty(regexpi(Cam,'^m_','once'))
    %% #RangeM_ : meshing parameters
    if iscell(val)
      out=struct('val',[],'lab',{{CAM}}, ...
        'param',struct(CAM,struct('type','pop','level',10)));
      [out.param.(CAM),out.val]=feval(fe_range('@popMerge'),out,CAM,val);
    else
      out=struct('val',val(:),'lab',{{CAM}}, ...
        'param',struct(CAM,struct('type','double','level',10)));
    end
    
  elseif ~isempty(regexp(Cam,'cfg$','once'))
    %% #RangeTypeCfg : support choices of
    % see more recent in sdtweb d_cor RangeTypeCfg
    if ischar(val); val={val}; end
    if iscell(val);
      [i1,r2]=ismember(lower(val),lower(Range.param.(CAM).choices));
      if ~all(i1); error('%s not found in %s',comstr(val(~i1),-30),CAM);end
    else; r2=val(:);
    end
    out=struct('val',r2(:),'lab',{{CAM}});
    
  elseif comstr(Cam,'check')||isempty(Cam)
    %% #RangeCheck : verify that all things needed are there
    
    st={'MeshCfg','SimuCfg','RunCfg'};
    if isempty(Range);Range=struct('param',[],'FileName',{{'SHM'}});
      Range.param.MainFcn={'d_shm','range'}; % Allow callback for param setting
    end
    for j1=1:length(st) % Verify base params
      if ~isfield(Range.param,st{j1});Range=d_shm(st{j1},Range); end
    end
    if isstruct(val); % Allow automated fill from struct
      val.param=Range.param; val.FileName=Range.FileName;
      Range=fe_range('grid -FileName',val);
    end
    %remove multiple healthy
    if isfield(Range.param,'input')
      i2=strcmpi(Range.lab,'M_Def1_Ecoef');
      if any(i2) %only Ecoef
        i1=strcmpi(Range.lab,'input');in1=unique(Range.val(:,i1));
        if length(in1)>1
          for j1=in1'
            in2= find(Range.val(:,i1)== j1 & Range.val(:,i2)== 1);
            in2(1)=[];Range.val(in2,:)=[];
          end
        end
      end
    else;dbstack;keyboard;
    end
    out=Range; return;
    
  else; error('''%s'' not a known range',CAM);
  end
  
  if ~isempty(Range) % Concatenate
    Range=fe_range('grid -level -replace=1',Range,out);
    if isfield(out,'lab')&&isfield(out,'val')&&~isempty(out.val)
      st1=sprintf('@%s',out.lab{1});
      if ~any(strcmpi(Range.FileName,st1));Range.FileName{end+1}=st1;end
    end
    out=Range;
  end
  if isfield(out,'FileName')
    out.FileName(strcmpi('@wdtest',out.FileName))=[];
    out.FileName(strcmpi('@wv',out.FileName))=[];
  end
  
  
  %% #Cb : callbacks ---------------------------------
elseif comstr(Cam,'cb'); [CAM,Cam]=comstr(CAM,3);
  
  if comstr(Cam,'matdef'); [CAM,Cam]=comstr(CAM,7);
    %% #cbMatDef model = d_shm('cbMatDef MatID ECoef',model)
    
    model=varargin{carg};carg=carg+1;
    
    opt=comstr(CAM,-1);
    pl=feutil(sprintf('getpl %i',opt(1)),model);
    st=fe_mat('typemstring',pl);
    switch lower(st)
      case 'm_elastic.1'; pl([3 6])=pl([3 6])*opt(2);
      case 'm_elastic.4'; i2=3:8;pl(i2)=pl(i2)*opt(2);
      case 'm_elastic.5'; i2=[3 4 6:8];pl(i2)=pl(i2)*opt(2);
      case 'm_elastic.6'; i2=[3:5 9:11]; pl(i2)=pl(i2)*opt(2);
      otherwise
        warning('%s not implemented',st);
    end
    out=feutil('setmat',model,pl);
    
    
  elseif comstr(Cam,'figdisppiezo'); [CAM,Cam]=comstr(CAM,13);
    %% #CBFigDispPiezo : display piezo output in Matlab Figure
    
    R1=struct('CurveFcn',{{'nida13','LoadRes'}},...
      'type','gridpos','IgnoreCol',0,'SortCol','auto');
    C1=fe_range(obj,evt,'Sel',R1);
    R2=struct;
    R2.objSet={'@title(1)',{'@toString',...
      {'@sprintf(''%s - Pz%i act'',uf.config,uf.Actuator)'}}};
    nida13('ViewCoSubElec-SmartLeg-keep-gf11-ImWrite',C1,R2);
    
  elseif comstr(Cam,'disppiezo'); [CAM,Cam]=comstr(CAM,10);
    %% #CBDispPiezo : display piezo output in iiplot
    RP=sdtroot('PARAM.Project');
    if isempty(RP); wd=pwd;
    else; wd=RP.ProjectWd;
    end
    
    [ind,ua]=feval(iimouse('@LinkedCh'),obj,evt,'tgetch row');
    st={};
    for j1=1:size(ua.row,1)
      C1=load(fullfile(wd,ua.row{j1,1}),'C1'); C1=C1.C1;
      C1.X{1}=C1.X{1}*1000; C1.Xlab{1}='Time [ms]';
      if ischar(C1.Xlab{1}) && strcmpi(C1.Xlab{1},'time')
        C1.Xlab{1}=fe_curve('datatypecell','time');
      end
      st(end+1,1:3)={'curve',ua.row{j1,1},C1}; %#ok<AGROW>
    end
    RB=fe_range('dirShortName',st(:,2));st(:,2)=RB.list;
    iicom('curveinit',st);
    
    
  elseif comstr(Cam,'dispdef'); [CAM,Cam]=comstr(CAM,8);
    %% #CBDispDef : display deformation in feplot
    'XXXJP TODO'
    keyboard
  elseif comstr(Cam,'movie'); [CAM,Cam]=comstr(CAM,6);
    %% #CBMovie : standard movie generation
    
    if isempty(Cam)
      model=varargin{carg};carg=carg+1;
      
      def=varargin{carg};carg=carg+1;
      
      if isa(model,'sdth');cf=model;cf.def=def;
      else
        cf=comgui('guifeplot');feplot(cf,model,def);
      end
      cf.sel(1)={'-linface','colordataevala -edgealpha0'};
      r1=jet(8);r1=reshape([r1';ones(3*9,size(r1,1))],3,[])';colormap(cf.ga,r1);
      fecom('colorscale instant')
      
      cf.sel(2)='seledgemat'; cf.o(2)='sel 2 ty4'; % undef
      fecom('scc1e-6')
      sdtweb('_link',sprintf('d_shm(''CBMovie %s'')',fullfile( ...
        sdtroot('param.Project.ProjectWd'),'Test.gif')));
    elseif ~isempty(strfind(Cam,'gif'))
      %% Actually generate
      r1=struct('FileProp',{{'Loopcount',2,'delayTime',.2}},'Reset',1);
      fecom('ch9');r1.FileName=CAM;fecom('animMovie',r1);
      % sdtweb nida13 echo.gif
    end
    
    
  elseif comstr(Cam,'initelecoutput'); [CAM,Cam]=comstr(CAM,8);
    %% #CBInitElecOutput :
    model=[];t=[]; opt=[];
    eval(iigui({'opt','model','Case','t'},'GetInCaller'))
    sens=fe_case(model,'sens');
    i1=fe_c(Case.DOF,sens.DOF(any(sens.cta)),'ind');
    %i1=int32(find(any(sens.cta)));
    out=struct('def',zeros(length(i1),length(t)),'data',t(:), ...
      'DOF',Case.DOF(i1),'fun',[0 4],'OutInd',int32(i1),'cur',zeros(1,3));
    %dbstack;keyboard;
    if ~isfield(opt,'OutputFcn')||~ischar(opt.OutputFcn)
      opt.OutputFcn='of_time(''interp'',out,beta,gamma,Case.uva,a,tc-dt,tc)';
    end
    eval(iigui({'opt','out'},'SetInCaller'))
    
  elseif comstr(Cam,'finalcleanup'); [CAM,Cam]=comstr(CAM,13);
   %% #CBFinalCleanupFcn :
   
   t=[];out=[];ft=[];model=[];
   eval(iigui({'model','out','Case','t','ft'},'MoveFromCaller'))
   %   %'out.DOF=model.DOF;out.def=[Case.T*out.def+Case.TIn*[0 ft(1:end-1)'']];';
   if isequal(t,out.data)
    u=[0 ft(1:end-1)''];
   else;u=of_time('lininterp',[t(:) [0 ft(1:end-1)']'],out.data(:),zeros(3,1))';
   end
   if ~isfield(Case,'mDOF');Case.mDOF=model.DOF;end
   if length(out.DOF)<length(Case.DOF) % Just selected DOFs (electric)
    c=Case.T(:,fe_c(Case.DOF,out.DOF,'ind'));
    i1=find(any(c,2)|any(Case.TIn,2));
    out.def=c(i1,:)*out.def+Case.TIn(i1,:)*u;
    out.DOF=Case.mDOF(i1);
   elseif size(out.def,1)==size(Case.T,2) % All DOFs 
    out.DOF=Case.mDOF;out.def=Case.T*out.def;    
    if isfield(Case,'TIn');out.def=out.def+Case.TIn*u;end   
    if isfield(out,'v');out.v=Case.T*out.v; if isfield(Case,'TIn'); out.v=out.v+Case.TIn*filter([1 -1],1,u);end;end   
   end
   eval(iigui({'out'},'SetInCaller'))
    
  elseif comstr(Cam,'finalmulti'); [CAM,Cam]=comstr(CAM,11);
    %% #CBFinalMulti 
    %copy from cbfinalcleanup: need to be cleaned!
    
    t=[];out=[];ft=[];model=[];Case=[];Load=[];
    eval(iigui({'model','out','Case','t','ft','Load'},'MoveFromCaller'))
    %   %'out.DOF=model.DOF;out.def=[Case.T*out.def+Case.TIn*[0 ft(1:end-1)'']];';
    
    if ~isfield(Case,'mDOF');Case.mDOF=model.DOF;end
    if ~isfield(Case,'TIn');eval(iigui({'out'},'SetInCaller'));return;end
    
    %add for multiples curves @EM: must be reviewed for multiples DOFSET (U&T)
    if isfield(Load,'bset')&&~isempty(Load.bset.def)
      %c1=fe_curve(Load.bset.curve,out.data);ft=c1.Y;
      ft=fe_curve('returny',model,Load.bset.curve{:},out.data);
    end
    if isequal(t,out.data)
      u=[0 ft(1:end-1)''];
    else;u=of_time('lininterp',[t(:) [0 ft(1:end-1)']'],out.data(:),zeros(3,1))';
    end
    
    if length(out.DOF)<length(Case.DOF) % Just selected DOFs (electric)
      c=Case.T(:,fe_c(Case.DOF,out.DOF,'ind'));
      i1=find(any(c,2)|any(Case.TIn,2));
      out.def=c(i1,:)*out.def+Case.TIn(i1,:)*u;
      out.DOF=Case.mDOF(i1);
      if isfield(out,'v');
        out.v=c(i1,:)*out.v;%+Case.TIn(i1,:)*filter([1 -1],1,u);
        out=rmfield(out,'OutInd');
      end
    elseif size(out.def,1)==size(Case.T,2)&&isfield(Case,'TIn'); % All DOFs
      out.DOF=Case.mDOF;
      out.def=Case.T*out.def+Case.TIn*u;
      if isfield(out,'v'); out.v=Case.T*out.v+Case.TIn*filter([1 -1],1,u);end
    end
    eval(iigui({'out'},'SetInCaller'))
    
  else; error('Cb%s unknown',CAM);
  end
  
  
  %% #Init : Init fcn ---------------------------------
elseif comstr(Cam,'init'); [CAM,Cam]=comstr(CAM,5);
  
  if comstr(Cam,'rvar'); [CAM,Cam]=comstr(CAM,5);
    %% #InitRVar - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    RO=varargin{carg}; carg=carg+1;
    st=CAM;
    if isfield(RO,'reset'); st=sprintf('%s',st); end
    if isfield(RO,'wd')
      sdtroot('SetProject',struct('ProjectWd',RO.wd,'UI','sdtroot'))
    elseif ~isfield(RO,'wd')&&isfield(RO,'val')&&isfield(RO,'FileName')
    end
    
    RO.MousePressed={'ContextCbk',{'tag','label','callback','level';
      'rvar_context','','',1; % contextmenu root
      'rvar_displayfig','Display piezo in matlab figure','d_shm(''CBFigDispPiezo'')',2;
      'rvar_comparetest','Compare simulated piezo signal to test signals','nida13(''LoadRes-plot'')',2;
      'rvar_displayiiplot','Display piezo in iiplot','d_shm(''CBDispPiezo'')',2;
      'rvar_displayfeplot','Display def in feplot','d_shm(''CBDispDef'')',2;
      'rvar_displaymovie','Generate movie','d_shm(''CBMovie'')',2;
      }};
    
    sdtroot(sprintf('InitRVar%s',st),'sdtroot',RO)
    
  else; error('Init%s unknown');
  end
  
elseif comstr(Cam,'delam');[CAM,Cam]=comstr(CAM,6);
  %% #delam
  if nargin>=carg;model=varargin{carg};carg=carg+1;else;return;end
  
  if comstr(Cam,'1d')
   %% #delam1D  -2
  %get defects
  RB=stack_get(model,'info','defect',3);prolay=vertcat(RB.list{:,3});
  in0=prolay(:,2)>1; in1=unique(prolay(:,2));in1(in1==0)=[];
  
  if any(in0)
    %get base plate (used & PZT)
    r1=p_piezo('electrodedof',model);ipzt=vertcat(r1{:,2}); npzt=feutil('getnode ',model,ipzt);
    nodes=feutilb('GetUsedNodes',model); mpid0=feutil('mpid',model);
    ind0=mpid0(:,2)>0;ipro=unique(mpid0(ind0,2));
    %get stack
    pro1=feutil(sprintf('GetIl %i',RB.ProId),model);pro2=pro1(10:end);
    %RO=feutil('getdd',[0 RB.ProId],model);RO.layer
    nlay=length(pro2)/4;  tlay=pro2(2:4:end);  %get thickness [matid h alpha 0]
    in2=[0;in1;nlay]+1; in3=[4*[0;in1]+9;length(pro1)]; hnew=zeros(length(in2)-1,1);    %new stacks
    for j1=1:size(hnew,1);hnew(j1)=sum(tlay(in2(j1):in2(j1+1)-1));end
    zlay=[0;cumsum(hnew)]-hnew(1)/2;  zref=(zlay(1:end-1)+zlay(2:end))/2;
    %duplicate base model
    mo1=model;mo1.Node=[];mo1.Elt=[];ofspro=ceil(max(model.il(:,1))/1e2)*1e2;
    RB.allpro=model.il(:,1)*ones(1,length(zref));    
    for j1=1:length(zref)
      %nodes, elt & pro
      mo2=struct('Node',nodes,'Elt',model.Elt);
      if j1==1;mo2.Node=model.Node;oldnode=nodes(:,1);end
      mo2.Node(:,7)=zref(j1);mpid=mpid0;ind=[1:9 in3(j1)+1:in3(j1+1)];
      if j1<length(zref);
        mo2.il=ones(length(ipro),1)*pro1(ind);mo2.il(:,1)=ipro+j1*ofspro; 
        mpid(ind0,1:2)=mpid(ind0,1:2)+j1*ofspro;  
        RB.allpro(:,j1)=RB.allpro(:,j1)+j1*ofspro;
      elseif j1==length(zref);
        mo2.il=pro1(ind);
      end       
      mo2.Elt=feutil('mpid',mo2,mpid);mo2.il(:, 3)= zlay(j1); % -hnew(j1)/2
      %merge
      if j1==1;mo1.Node=mo2.Node;mo1.Elt=mo2.Elt; mo1.il=mo2.il;
      else
        [mo1,in4]=feutil('addtest -noOri;',mo1,mo2);  oldnode(:,j1)=mo1.Node(in4,1);
        mo1.Elt=feutil('addelt',mo1.Elt,'rigid',[oldnode(:,j1-1:j1)  123456*ones(size(in4,1),1)]);
      end      
    end
    [~,ia]=ismember(ipro,RB.allpro(:,end));ia=setdiff(1:size(model.il,1),ia)';
    RB.allpro(ia,1:end-1)=0;
    %cut PZT/damage ProId
    in5=unique(mpid0(ind0 & mpid0(:,2)~=RB.ProId,2)); 
    nlayPZT=(in2(end)-in2(end-1));iend=size(mo1.il,2);
    for j1=in5'
      il=feutil(sprintf('getil %d',j1),model); m_function=fe_mat('typep',il(2));
      if comstr(m_function,'p_piezo')
        il(5)=nlayPZT;il0=feutil(sprintf('getil %d',il(3)),model);
        il0=il0(ind);il0(iend)=0;il0(3)=zlay(end-1);
        mo1.il=[mo1.il;il0;il(1:iend)];
      elseif comstr(m_function,'p_shell')
        il=il(ind);il(3)=zlay(end-1);mo1.il=[mo1.il;il(1:iend)];
        [~,ia]=ismember(j1,prolay(:,1));
        if ia
          nid=feutil(sprintf('findnode inelt{proid %d}',j1),mo1);
          [~,ib]=ismember(nid,oldnode(:,end));RB.list{ia,4}=oldnode(ib,:);
          ic=find(prolay(ia,2)>in2,1,'last'); if isempty(ic);ic=0;end
          RB.list{ia,4}(:,end+1)=ic;RB.list{ia,3}(3)=ic;
        end
      else dbstack;keyboard;
      end
    end
  else out=model;return;
  end
  %check intersection %reverse for xixi
  for j1=1:size(RB.list,1)-1
   r1=RB.list{size(RB.list,1)+1-j1,4};
   for j2=j1+1:size(RB.list,1)
    j3=size(RB.list,1)+1-j2;
    r2=RB.list{j3,4};RB.list{j3,4}=setdiff(r2,r1,'rows');
   end
   %mo1.Elt=feutil('addelt',mo1.Elt,'rigid',[r1 123456*ones(size(r1,1),1)]);
  end
  out=stack_set(mo1,'info','defect',RB); 
  
  elseif comstr(Cam,'3d')
   %% #delam3D  -2
   R1=stack_get(model,'info','defect',3);
   listmat=[0 0 0];matid=1;nldata=vertcat(R1.list{:,3});
   mo1=struct('Node',[],'Elt',[],'pl',[],'il',[]);
   for j1=1:size(model.il,1)
    p_func=fe_mat('typep',model.il(j1,2));
    if strcmpi(p_func,'p_shell')
     % if p_shell
     lam=reshape(model.il(j1,10:end),4,[])';lam(lam(:,1)==0,:)=[];
     %RO=feutil('getdd',[0 model.il(j1,1)],model);
     [curmat,~,imat]=uniquetol(lam(:,[1,3]),'byrows',1);
     [~,in1]=ismembertol(curmat,listmat(:,1:2),'byrows',1);in2=find(in1==0);
     if ~isempty(in2)
      %add new mat if needed
      listmat(end+(1:length(in2)),1:2)=curmat(in2,:);
      in1(in2)=size(listmat,1)+(-length(in2)+1:0);
      for j2=in2'
       curpl=feutil(sprintf('get pl%d',curmat(j2)),model);
       [m_func,m_unit,m_typ]=fe_mat('type',curpl(2));
       listmat(in1(j2),3)=matid;matid=matid+1;
       curpl(1)=listmat(in1(j2),3);
       if strcmpi(m_func,'m_elastic')
        %m_elastic
        if ismember(m_typ,[2,4,5]);dbstack;keyboard;end %incompatible type (shell)
        mo1=feutil('setmat',mo1,curpl);
        if curmat(j2,2)
         RO=feutil('getdd',[listmat(in1(j2),3) 0],mo1);
         bas=basis('rotate',[],sprintf('rz=%f;',listmat(in1(j2),2)),1);
         dd=m_elastic(sprintf('formulaPlAniso %d',listmat(in1(j2),3)),RO.dd,bas);
         dd(abs(dd)<1e-9)=0;[m_func2,~,m_typ2]=fe_mat('type',dd(2));
         dd(2)=fe_mat(m_func2,m_unit,m_typ2);mo1=feutil('setmat',mo1,dd);
        end
       elseif strcmpi(m_func,'m_piezo')
        [iok,ic]=ismember(curpl(3),listmat(:,1));
        if ~iok
         plmeca=feutil(sprintf('get pl%d',curpl(3)),model);
         listmat(end+1,[1 3])=[curpl(3) matid];matid=matid+1;
         curpl(3)=listmat(end,3);plmeca(1)=listmat(end,3);
         mo1=feutil('setmat',mo1,curpl);
         mo1=feutil('setmat',mo1,plmeca);
        else
         curpl(3)=listmat(ic,3);mo1=feutil('setmat',mo1,curpl);
        end
       else
        dbstack;keyboard %m_func
       end
      end
     end
    elseif strcmpi(p_func,'p_piezo')
     % if p_piezo
    else dbstack;keyboard %p_func
    end
    pref=[model.il(j1,1) fe_mat('p_solid','MM',1) 0 -3 0 1]; 
    mo1=feutil('setpro',mo1,pref);
   end
   %sort proid
   [~,ia]=sort(mo1.il(:,1));
   mo1.il=mo1.il(ia,:);
   %build geometry
   for j0=1:size(R1.allpro,2)
    mo3=struct('Node',[],'Elt',[],'pl',[],'il',[]);
    for j1=1:size(R1.allpro,1)
     mo2=struct('Node',model.Node,'Elt',[],'pl',[],'il',[]);
     [~,mo2.Elt]=feutil(sprintf('findelt proid %d',R1.allpro(j1,j0)),model);
     if ~isempty(mo2.Elt)
      %check type
      curil=model.il(model.il(:,1)==R1.allpro(j1,j0),:);
      [p_func,p_unit,p_typ]=fe_mat('typep',curil(2));
      if strcmpi(p_func,'p_shell')
      elseif strcmpi(p_func,'p_piezo')
       curil=model.il(model.il(:,1)==curil(3),:);
      else %p_typ
       keyboard;
      end
      %clean model
      mo2.Node=feutilb('GetUsedNodes',mo2);mo2=feutil('joinall',mo2);
      mo2=feutil('renumber -NoOri',mo2);
      lam=reshape(curil(10:end),4,[])';lam(lam(:,1)==0,:)=[];
      [~,imat]=ismembertol(lam(:,[1,3]),listmat(:,1:2),'byrows',1);
      if any(imat==0);dbstack;keyboard;end
      h1=curil(3)+[0;cumsum(lam(:,2))]; %RO=feutil('getdd',[0 curil(1)],model);RO.layer
      mo2=feutil('extrude 0 0 0 1',mo2,h1-mean(mo2.Node(1,7)));mpid=feutil('mpid',mo2);
      %assign mat
      for j2=1:size(lam,1)
       ie=feutil(sprintf('findelt innodes{z>=%f & z<=%f}',h1(j2:j2+1)),mo2);
       mpid(ie,1)=listmat(imat(j2),3);
      end
      mo2.Elt=feutil('mpid',mo2,mpid);
      %join layer
      if j1==1;mo3.Node=mo2.Node;mo3.Elt=mo2.Elt;
      else;mo3=feutil('addtest Merge-Edge -noorig;',mo3,mo2);
      end
     else %empty elt
      mo1.il(mo1.il(:,1)==R1.allpro(j1,j0),:)=[];
     end
    end
    %join layers
    if j0==1;mo1.Node=mo3.Node;mo1.Elt=mo3.Elt;
    else
     in1=find(nldata(:,3)==j0-1);
     if any(in1)
      [~,ir]=ismember(nldata(in1,1),R1.allpro(:,end));ipro=R1.allpro(ir,j0-1:j0);
      offset1=max(mo3.Node(:,7))*10;  allN=[]; allN0=[];
      %get nodes in previous layer
      for j2=1:size(ipro,1)
       [nid3,nodes3]=feutil(sprintf('findnode z==%f & inelt{proid %d}',min(mo3.Node(:,7)),ipro(j2,2)),mo3);
       [nodes1,ind1]=feutil('addnode', mo1.Node, nodes3);
       if size(nodes1,1)>size(mo1.Node,1);dbstack;keyboard;end %new nnodes!!
       allN0=[allN0;ind1];allN=[allN;nid3];R1.list{in1(j2),4}=[ind1 nid3];
      end
      %join
      mo3.Node(allN,7)=mo3.Node(allN,7)+offset1;
      [mo1,ind3]=feutil('addtest Merge-Edge -noorig;',mo1,mo3);
      mo1.Node(ind3(allN),7)=mo1.Node(allN0,7)+1e-12;%mo1.Node(ind3(allN),7)-offset1;
      %get nodes in next layer
      for j2=1:size(ipro,1)
       RB=R1.list{in1(j2),4};RB(:,2)=ind3(RB(:,2));
       %mo1.Elt=feutil('addelt',mo1.Elt,'rigid',[RB  123*ones(size(RB,1),1)]);
       R1.list{in1(j2),4}=RB;
      end
     else
      mo1=feutil('addtest Merge-Edge -noorig;',mo1,mo3);
     end
    end
   end
   mo1.Elt=feutilb('SeparateByMat',mo1.Elt);
   %check intersection
   for j1=1:size(R1.list,1)
    r1=R1.list{size(R1.list,1)+1-j1,4};
    for j2=j1+1:size(R1.list,1)
     j3=size(R1.list,1)+1-j2;
     r2=R1.list{j3,4};R1.list{j3,4}=setdiff(r2,r1,'rows');     
    end
    mo1.Elt=feutil('addelt',mo1.Elt,'rigid',[r1 123456*ones(size(r1,1),1)]);
   end   
   mo1=stack_set(mo1,'info','defect',R1);
   %electrode
   r1=p_piezo('electrodedof',model);RB=ones(size(r1,1),2);
   r1=str2double(regexp(r1(:,1),'\d*','once','match'));
   for j1=1:size(r1,1)
    r2=model.il(model.il(:,1)==r1(j1),[3 5]);
    m0=model.il(model.il(:,1)==r2(1),6+4*r2(2));
    m1=listmat(listmat(:,1)==m0,3);
    [nid,nodes]=feutil(sprintf('findnode inelt{proid %d & matid %d}',r1(j1),m1),mo1);
    [~,in1]=ismembertol(nodes(:,7),[min(nodes(:,7)),max(nodes(:,7))]);
    ntop=nid(in1==2);ndown=nid(in1==1);
    [mo1,InputDOF]=p_piezo(sprintf('ElectrodeMPC Neg%d -Ground',r1(j1)),mo1,nid(in1==1));
    [mo1,InputDOF(end+1,1)]=p_piezo(sprintf('ElectrodeMPC Pos%d',r1(j1)),mo1,nid(in1==2));
    RB(j1,1)=fix(InputDOF(end,1));
   end   
   %RB(1,2)=1; %actuator
   out=stack_set(mo1,'info','Electrodes',struct('data',RB,'ver',1));   
   
  else
   warning('no delam');%dbstack;keyboard
   out=model;
  end
  
elseif comstr(Cam,'view');[CAM,Cam]=comstr(CAM,5);
  %% #View
  if nargin>=carg;RO=varargin{carg};carg=carg+1;else;return;end
  if isfield(RO,'Node');cf=feplot(RO);
  elseif isfield(RO,'param')
    in1=strcmpi(RO.lab,'MeshCfg');i1=RO.val(1,in1);%get first only
    cf=feplot(d_piezo('meshplate',RO.param.MeshCfg.data{i1}));
  elseif isfield(RO,'list');cf=feplot(d_piezo('meshplate',RO));
  else;out=[];return;
  end
  fecom(';colordatapro;coloredgealpha.1;triaxon');out=cf;
  
elseif comstr(Cam,'dis');[CAM,Cam]=comstr(CAM,4);
  %% #DIs
  if nargin>=carg;Range=varargin{carg};carg=carg+1;else;Range=[];end
  
  if ~exist('compute_dis','file')
    warning('SHM@PIMM not available'); out=[];return;
  end
  
  if isempty(Range); out = compute_dis();
  else
    if isfield(Range.param,'input')
      i1=strcmpi(Range.lab,'input');i2=strcmpi(Range.lab,'M_Def1_Ecoef');
      if any(i2) %only Ecoef
        keyboard
      end
    end
  end
  
elseif comstr(Cam,'export');[CAM,Cam]=comstr(CAM,7);
  %% #Export  
   if nargin>=carg;Range=varargin{carg};carg=carg+1;else;warning('Missing argument Range');return;end
   if nargin>=carg;RO=varargin{carg};carg=carg+1;else;RO=[];end
   
   if isempty(RO);RO=struct('nb_rep',10,'snr_ratio',40:5:85);end
   
   %get structure
   i1=uniquetol(Range.val(:,strcmpi(Range.lab,'MeshCfg')));
   if length(i1)>1;dbstack;keyboard;end
   sname=Range.param.MeshCfg.choices{i1};
   
   %directory
   wd=sdtroot('PARAM.Project.ProjectWd -safe');
   RO.out_dir_name =fullfile(wd,[sname datestr(now,'_HH_MM_dd_mmm_yyyy')]) ;
   mkdir(RO.out_dir_name); % save([out_dir_name '/myRO'],'myRO') ;
   
   % Get current parameters
    % Act (Actuator) / Size (M_Def1_Ecoef) / Dx (M_Def1_xc) / Dy (M_Def1_yc) / lc / Plate
    % Run / Simu / Tend / dt (dt) / Amp / freq (inf0)
   inIN=strcmpi(Range.lab,'input');inac=strcmpi(Range.lab,'Actuator');
   inf0=strcmpi(Range.lab,'amp');inf0=strcmpi(Range.lab,'inf0');indt=strcmpi(Range.lab,'dt');
   inxd=strcmpi(Range.lab,'M_Def1_xc');inyd=strcmpi(Range.lab,'M_Def1_yc');
   inrd=strcmpi(Range.lab,'M_Def1_rc');incd=strcmpi(Range.lab,'M_Def1_Ecoef');
   for nb_sim = 1:size(Range.Res,1)
     disp(['Loading file ' num2str(nb_sim) ' out of ' num2str(size(Range.Res,1)) ' ...'])
     %check input format
     curr_param.input_freq   = Range.val(nb_sim,inf0);
     if isempty(curr_param.input_freq)  %no input parameters
       if any(inIN)
         input=Range.param.input.data{Range.val(nb_sim,inIN)};
         curr_param.input_freq =input.f0; curr_param.input_amp =input.A;
         try curr_param.nb_cycles = input.nbcycle; end     %check other signals
       else
         dbstack;keyboard;
       end
     end      
     curr_param.sampl_freq   = 1/Range.val(nb_sim,indt) ;
     curr_param.damage_loc_x = Range.val(nb_sim,inxd);
     curr_param.damage_loc_y = Range.val(nb_sim,inyd);
     curr_param.damage_size  = Range.val(nb_sim,incd);
     curr_param.damage_radius  = Range.val(nb_sim,inrd);    
     curr_param.nb_rep       = RO.nb_rep ;
     curr_param.n_act        = Range.val(nb_sim,inac);
     curr_param.file_path    =fullfile(wd,[Range.Res{nb_sim,2}.name '_RES.mat']) ;
     curr_param.jpar = nb_sim;
     
     for n_SNR = 1:length(RO.snr_ratio)
        curr_param.SNR  = RO.snr_ratio(n_SNR) ;       
       % Eventually create directory
       if isempty( curr_param.damage_radius) && (isempty( curr_param.damage_size)  || curr_param.damage_size==1)
         out_name = ['Healthy_SNR' num2str(curr_param.SNR )] ;
       else         
        out_name = ['Dam' num2str(curr_param.damage_radius) '_X' num2str(curr_param.damage_loc_x) '_Y' num2str(curr_param.damage_loc_y) '_SNR' num2str(curr_param.SNR )] ;
       end
       curr_param.dir_name = [RO.out_dir_name '/' out_name '/' num2str(curr_param.input_freq/1000) 'kHz_' num2str(curr_param.nb_cycles) 'cycles/Actionneur' num2str(curr_param.n_act) '/' ];
       if ~isdir(curr_param.dir_name); mkdir(curr_param.dir_name) ; end
       % Export simulated data to the SHM@PIMM format
       convert2SHM(Range,curr_param);
       % Update structure info lines
       %update_structure_info(curr_param,out_dir_name,out_name) ;
     end
   end
  
  %% clean end
elseif comstr(Cam,'@');out=eval(CAM);
elseif comstr(Cam,'cvs')
 out=sdtcheck('Revision');
else; error('%s unknown',CAM);
end
%% #End function
end


%% #DefaultRunFcn
function TimeRunFcn;
 dbstack; keyboard; 
end

function DefaultRunFcn;
  % Start transient using piezo actuation then continue with all sensors

 model=[];eval(iigui({'model','RO'},'MoveFromCaller'))
 % Build time opt:
 opt=stack_get(model,'','TimeOpt','get');
        
        % model=fe_case(model,'pcond','Piezo','d_piezo(''Pcond 1e12'')');
        % Do not constrain drilling DOF unless drilling stiffness set to 0
        % if all(model.Node(:,7)==0);model=fe_case(model,'fixdof','Drill',.06);end
if isfield(opt,'SaveTimes')
  dbstack;fprintf('ChangeName'); 
         % First transient with piezo as actuator
         tic;NStep=opt.Opt(5);opt.Opt(5)=opt.SaveTimes(2);
  opt.SaveFcn='q0=struct(''DOF'',[mdof;Case.DofIn],''def'',[u v a;0 0 0]);assignin(''caller'',''q0'',q0);';
         def=fe_time(stack_set(model,'info','TimeOpt',opt));
     fprintf('Done actuation transient in %.1f s\n',toc)
         
     %% Now do an observation transient
     opt.Opt([3 5])=[def.data(end) NStep-opt.Opt(5)];
     opt=rmfield(opt,{'SaveTimes','SaveFcn'}); model=fe_case(model,'Remove','V_IN');
     model=stack_set(model,'info','q0',q0);

     if isfield(opt,'InitAcceleration');opt.InitAcceleration=';';end %reuse existing

     def1=fe_time(stack_set(model,'info','TimeOpt',opt));
     fprintf('Done observation transient in %.1f s\n',toc)

     if norm(def.DOF-def1.DOF);
         in1=fe_c(def.DOF,def1.DOF,'ind');def.DOF=def.DOF(in1);def.def=def.def(in1,:);
     end % reordering DOF
     if abs(def.data(end)-def1.data(1))>=opt.Opt(4);dbstack;keyboard;end % need time check
     if norm(def.def(:,end)-def1.def(:,1))>eps*1e-8;dbstack;keyboard;end % need IC check
     def1.data(1)=[];def1.def(:,1)=[];Trigger=[zeros(length(def.data),1);ones(length(def1.data),1)];
     def=fe_def('appenddef',def,def1);toc 
     
else;tic;def=fe_time(stack_set(model,'info','TimeOpt',opt));toc %@EM
end
        
if isfield(model,'name'); def.name=model.name; end
sens=fe_case(model,'sens','V_OUT');
C1=fe_case('SensObserve',sens,def);
C1.name=['Elec_',RO.fname{:}];'xxx store range'

if exist('Trigger','var');C1.X{2}{end+1}='Trigger';C1.Y(:,end+1)=Trigger;end
      % find the "end" of the actuator signal
    r5=stack_get(model,'info','SimuInfo','get'); RO.Actuator=r5.Actuator;
    %  tol=1e-5;
    %  r5=abs(C1.Y(:,RO.Actuator)); t=C1.X{1};
    %  i6=fix(0.95*length(r5)); r5(i6:end)=[]; t(i6:end)=[];
    %  r5=flipud(r5); t=flipud(t);
    %  i5=(r5/max(r5)<tol); [~,i5]=max(find(cumprod(i5)>0)); RO.Tmax=t(i5);
    try;RO.Tmax=computeTmax(C1,RO);
    catch; warning('Tmax compute failed');
    end
    
   eval(iigui({'def','C1','model','RO'},'SetInCaller'))
end    


%% #flatParam : callback for range rebuilding
function [r1,RO]=flatParam(li,RO);

% Files, RO, Stat
%r5=who('-file',li.fname);
if isfield(li,'model'); r1=li;
else
  if isfield(RO,'wd'); li.wd=RO.wd; end
  r1=load(li.fname); %,'RO','model');
  if ~isfield(r1,'RO');
    fprintf('''%s'' with no RO.\n',li.fname);
    r1.RO=struct;
  end
end
if isfield(r1,'model') % get Mesh parameters
  % Step Mesh parameters :
  r2=stack_get(r1.model,'info','MeshInfo','getdata');
  % turn list as parameters
  if isfield(r2,'list')
    for j1=2:size(r2.list,1)
      st2=r2.list{j1,3}; if ~isempty(st2);st2=fieldnames(st2);end
      for j2=1:length(st2)
        r3=r2.list{j1,3};
        r1.RO.(sprintf('M_%s_%s',r2.list{j1,1},st2{j2}))=r3.(st2{j2});
      end
    end
    % compute distances between zp and defaults
    i1=strncmpi('pz',r2.list(:,1),2)|strncmpi('def',r2.list(:,1),2);
    st=r2.list(i1,1);
    st2={}; st2(2:1+length(st),1)=st; st2(1,end+1:end+length(st))=st;
    for j1=1:length(st)
      for j2=1:length(st)
        x1=r1.RO.(sprintf('M_%s_xc',st{j1})); y1=r1.RO.(sprintf('M_%s_yc',st{j1}));
        x2=r1.RO.(sprintf('M_%s_xc',st{j2})); y2=r1.RO.(sprintf('M_%s_yc',st{j2}));
        d=sqrt((x1-x2)^2+(y1-y2)^2);
        st2{j1+1,j2+1}=d;
      end
    end
    r1.RO.dist=st2;
  end
  % StepSimu : SimuInfo
  r2=stack_get(r1.model,'info','SimuInfo','getdata');
  if ~isempty(r2)
    st2=fieldnames(r2);
    for j2=1:length(st2); r1.RO.(st2{j2})=r2.(st2{j2}); end
  end
  if isfield(r2,'wdtest') % corresponding test data
    r1.RO.wdtest=fullfile(r2.wdtest,sprintf('Actionneur%i',r2.Actuator)); % This is specific to nida13
  end
  % StepRun : RunInfo
  r2=stack_get(r1.model,'info','RunInfo','getdata');
  if ~isempty(r2)
    st2=fieldnames(r2);
    for j2=1:length(st2); r1.RO.(st2{j2})=r2.(st2{j2}); end
  end
end

r1=r1.RO;
r1=feutil('rmfield',r1,{'Range','fname'});
%if isfield(li,'fname'); r1.fname=li.fname; end
end

%% #computeTmax : Compute Tmax
function [Tmax]=computeTmax(C1,RO);

if ~isfield(RO,'Actuator');RO=stack_get(RO.model,'','SimuInfo','get');end
tol=1e-5;
r5=abs(C1.Y(:,RO.Actuator)); t=C1.X{1};
i6=fix(0.95*length(r5)); r5(i6:end)=[]; t(i6:end)=[];
r5=flipud(r5); t=flipud(t);
i5=(r5/max(r5)<tol); [~,i5]=max(find(cumprod(i5)>0));
Tmax=t(i5);
end


%% #SHM@PIMM

%% #convert2SHM -2
function out=convert2SHM(Range,RO)
%feval(d_shm('@convert2SHM'),Range) or d_shm('convert2SHM',Range,RO)
if nargin==1;RO=[];end
if ~isfield(Range,'Res');out=[];return;end
if ~isfield(RO,'jpar');RO.jpar=1:size(Range.Res,1);end

out=cell(length(RO.jpar),2);
for j1=1:size(out,1)
  %clean
  C1=Range.Res{RO.jpar(j1),2};inPZ=1:size(C1.Y,2);inT=[];
  if strcmpi(C1.X{2}{end},'Trigger');inT=inPZ(end);inPZ(end)=[];end  
  st1=regexpi(C1.name,'actuator=','split'); if size(st1,2)~=2;dbstack;keyboard;end; 
  iACT=sscanf(st1{2},'%d',1);in1=[iACT setdiff(inPZ,iACT)]; 
  out{j1,1}=[C1.X{1} C1.Y(:,in1) ];
  if ~isempty(inT);out{j1,2}=C1.Y(:,inT);end
  %figure;plot(C1.X{1},C1.Y(:,iACT).*C1.Y(:,end))
  %figure;plot(C1.X{1},C1.Y(:,iACT).*(1-C1.Y(:,end)))
  
  if isfield(RO,'SNR')
   nb_act=length(inPZ);
   for n_rep = 1:RO.nb_rep
    % Add noise
    mean_total_energy = mean(sum(out{j1,1}(:,3:nb_act+1).^2));
    noise_generated = randn(size(out{j1,1}(:,3:nb_act+1))) ;
    mean_noise_energy = mean(sum(noise_generated.^2)) ;
    Time_Response = out{j1,1} ;
    Time_Response(:,3:nb_act+1) = Time_Response(:,3:nb_act+1) + 10^(-RO.SNR/20)*sqrt(mean_total_energy/mean_noise_energy)*noise_generated ;
    % Save data
    outname=fullfile(RO.dir_name,[ 'sim_data_rep_' num2str(n_rep) '.mat']);
    if ~isempty(inT);
     Trigger=C1.Y(:,inT);save(outname,'Time_Response', 'Trigger')
    else
     save(outname,'Time_Response');
    end
   end
  end
end

end


%% #check fields -2
function Ru=check_fields(Ru)

if isfield(Ru,'lc'); Ru.EltSize=Ru.lc; Ru=rmfield(Ru,'lc'); end
if isfield(Ru,'damage_loc_x'); Ru.DefXc=Ru.damage_loc_x;Ru=rmfield(Ru,'damage_loc_x'); end
if isfield(Ru,'damage_loc_y'); Ru.DefYc=Ru.damage_loc_y;Ru=rmfield(Ru,'damage_loc_y'); end
if isfield(Ru,'damage_radius'); Ru.DefDiam= Ru.damage_radius;Ru=rmfield(Ru,'damage_radius'); end
if isfield(Ru,'damage_size'); Ru.DefRed=Ru.damage_size;Ru=rmfield(Ru,'damage_size');end
if isfield(Ru,'actuators'); Ru.Actuator=Ru.actuators;Ru=rmfield(Ru,'actuators');end
if isfield(Ru,'duration'); Ru.tend=Ru.duration;Ru=rmfield(Ru,'duration');end
if isfield(Ru,'sampl_freq');Ru.dt=1/Ru.sampl_freq;Ru=rmfield(Ru,'sampl_freq');end
if isfield(Ru,'input_amp'); Ru.inA=Ru.input_amp;Ru=rmfield(Ru,'input_amp');end
if isfield(Ru,'input_freq'); Ru.inF=Ru.input_freq;Ru=rmfield(Ru,'input_freq');end

end


%% #genInput -2
function [allinput,R1]=genInput(R1)

RO=struct('type',[],'inA',[],'inF',[]);
if ~isfield(R1,'inA');R1.inA=1;end; RO.inA=R1.inA;
if ~isfield(R1,'inF');R1.inF=150e3;end; RO.inF=R1.inF;
if isstruct(R1.input);RO.type=1;tlab={R1.input.type,R1.input};
elseif iscell(R1.input);RO.type=1:length(R1.input);
  tlab=cellfun(@(x)x.type,R1.input,'uniformoutput',0);tlab(:,2)=R1.input;
else;dbstack;keyboard
end
Range=fe_range('grid',RO);

allinput=cell(size(Range.val,1),2);
for j1=1:size(Range.val,1)
  fcn=tlab{Range.val(j1,1)};amp=Range.val(j1,2);f0=Range.val(j1,3);
  name1=sprintf('%s_%iV_%ikHz',fcn,amp,f0/1e3);
  R2=tlab{Range.val(j1,1),2};if ~isfield(R2,'offset');R2.offset=0;end;
  if comstr(fcn,'sinus')
     tend=R2.nbcycle/f0;
    str1=sprintf('eval (t>%.15g & t<=%.15g).*%.15g.*sin(2*pi*%.15g*(t-%.15g))',...
     R2.offset,tend+R2.offset,amp,f0,R2.offset);
    cursig=struct('type',str1,'f0',f0,'nbcycle',R2.nbcycle,'A',amp,'name',sprintf('%s_%icycles',name1,R2.nbcycle));    
    if isfield(R2,'nbcycle');cursig.Tf=R2.nbcycle/f0;end;
  elseif comstr(fcn,'burst')    
    tend=R2.nbcycle/f0;
    str1=sprintf('eval (t>%.15g & t<=%.15g).*%.15g.*sin(2*pi*%.15g*(t-%.15g)).*sin(2*pi*%.15g*(t-%.15g)/(2*%i))',...
      R2.offset,tend+R2.offset,amp,f0,R2.offset,f0,R2.offset,R2.nbcycle);
    cursig=struct('type',str1,'f0',f0,'nbcycle',R2.nbcycle,'A',amp,'opt',[],'fcn',[],'name',sprintf('%s_%icycles',name1,R2.nbcycle));
  elseif comstr(fcn,'coshan')
    tend=R2.nbcycle/f0;
    str1=sprintf('eval (t>%.15g & t<=%.15g).*%.15g.*(1-cos((t-%.15g)*%.15g*2*pi/%.15g)).*cos((t-%.15g)*%.15g*2*pi)',...
      R2.offset,tend+R2.offset,amp/2.,R2.offset,f0,R2.nbcycle,R2.offset,f0);
    cursig=struct('type',str1,'f0',f0,'nbcycle',R2.nbcycle,'A',amp,'name',sprintf('%s_%icycles',name1,R2.nbcycle));
  elseif comstr(fcn,'sinhan')
    tend=R2.nbcycle/f0;
    str1=sprintf('eval (t>%.15g & t<=%.15g).*%.15g.*(1-cos((t-%.15g)*%.15g*2*pi/%.15g)).*sin((t-%.15g)*%.15g*2*pi)',...
      R2.offset,tend+R2.offset,amp/2.,R2.offset,f0,R2.nbcycle,R2.offset,f0);
    cursig=struct('type',str1,'f0',f0,'nbcycle',R2.nbcycle,'A',amp,'name',sprintf('%s_%icycles',name1,R2.nbcycle));
  elseif comstr(fcn,'sweep')
     if ~isfield(R2,'t1');R2.t1=R1.tend;end;
     if ~isfield(R2,'method');R2.method='linear';end;
     str1=sprintf('eval (t>%.15g & t<=%.15g).*%.15g.*chirp(t,%.15g,%.15g,%.15g,''%s'',-90)',...
      R2.offset,R2.t1+R2.offset,amp/2.,f0,R2.t1,R2.f1,R2.method);
    cursig=struct('type',str1,'f0',f0,'f1',R2.f1,'t1',R2.t1,'A',amp,'method',R2.method,'name',sprintf('%s_%ikHz_%s',name1,R2.f1/1e3,R2.method));
  elseif comstr(fcn,'tab')
    if size(R2.t,2)>1||size(R2.Y,2)>1;error('Signal with bad dimensions !');end
    cursig=struct('ID',1,'X',R2.t, 'Y',R2.Y,'data',[], 'xunit',[],'yunit',[],'unit',[],'name',R2.name);
  else;dbstack;keyboard
  end
  %save
  allinput(j1,:)={cursig.name,cursig};
end
R1=rmfield(R1,{'inA','inF'});
end
