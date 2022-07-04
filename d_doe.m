function [out,out1]=d_doe(varargin);

% D_DOE sample files for design of numerical experiments
%

%       Etienne Balmes, Guilherme Malacrida Alves
%       Copyright (c) 1990-2022 by SDTools, All Rights Reserved.
%       For revision information use d_tdoe('cvs')

if nargin==0
 fprintf('You should call a specific command, see:\n')
 sdtweb('_link','sdtweb(''_taglist'',''d_doe'')','See the list of commands')
 return
end

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*DEFNU,*NOSEM>
%persistent nmap;

[CAM,Cam]=comstr(varargin{1},1);carg=2;

%% #Script ------------------------------------------------------
if comstr(Cam,'script');[CAM,Cam]=comstr(CAM,7);

if comstr(Cam,'mount0d')
%% #Script0DMount : stepped sine on a spring model

sd _prafael;%sdtweb t_ident dec21_tdoe
% sdtweb _bp d_hbm FirstStab
RP=struct('MeshCfg','d_hbm(0D):DofSet:0dm1t','SimuCfg','SteppedSine{.5,1,10}:C0{0,15}:C1{2.5,10}', ...
         'NperPer',2e3,'Nper',1,'iteStab',20,'RunCfg','{run,dfr_ident@va}','PerFrac',.2);
d_tdoe('Solve',RP)

li={'MeshCfg{"d_hbm(0D):DofSet:0dm1t"}';';'
      'SimuCfg{RO{NperPer2e3,Nper1,iteStab20},"SteppedSine{5}:C0{0,15}:C1{2.5,10}"}';';'
      'RunCfg{Time,dfr_ident@va}'};
RT=struct('nmap',vhandle.nmap);
r2=sdtm.range(RT,li);%d2=mo2.nmap('CurTime');

%RP.Mesh='d_hbm(Mesh0D):t_vibrac(0Dm1tsvli)';dfr_ident('Load',RP); % Stresss rate relax + Dahl: OK

% RP=struct('MeshCfg','0D:cub','SimuCfg','DofLoad:Sine{1}:C0{-30,-15,0,15,30,45,60}:C1{2.5,10}');
% RP.RunCfg='{run,dfr_ident@va}'; 
% d_tdoe('Solve',RP);
 
 % {'ExitFcn'}    {'FirstStab'}    {1×1 struct}
 % {'info'   }    {'Range'    }    {1×1 struct}
 % {'ExitFcn'}    {'Default'  }    {1×1 struct}
 % sdtweb d_hbm SteppedSineDefFNL

elseif comstr(Cam,'duff2')
%% #ScriptDuff2DOF : stepped sine on a duffing

 li={'MeshCfg{d_hbm(Duffing2Dof),,CubFu}';';' 
     'SimuCfg{RO{NperPer2e3,Nper1},"SteppedSine{1}:C1(*.001=N){2.5,10}"}';';'
     'RunCfg{run,gui21@postproto}'};
 RT=struct('nmap',vhandle.nmap);
 RT.nmap('FirstStab')=struct('FinalCleanup','d_hbm@FirstStab','ite',10,'cond','d_hbm@CountIte');
 r2=sdtm.range(RT,li);%d2=mo2.nmap('CurTime');

li={'MeshCfg{d_hbm(Duffing2Dof),,CubFu}';';' 
     'SimuCfg{RO{NperPer2e3,Nper1},"SteppedSine{@ll(1,10,10)}:C1(*.001=N){2.5,10}"}';';'
     'RunCfg{run,gui21@postproto}'};
 RT=struct('nmap',vhandle.nmap);
 RT.nmap('FirstStab')=struct('FinalCleanup','d_hbm@FirstStab','ite',10,'cond','d_hbm@CountIte');
 r2=sdtm.range(RT,li);%d2=mo2.nmap('CurTime');



% sdtweb d_hbm firststab

 % RS.CbFcn=@gui21; 
 PA=sdtroot('paramVh'); PA.TDOE=struct('MeshCb','gui21');st=sdtroot('param.TDOE.MeshCb -safe');
 % sdtweb d_hbm FirstStab % xxx should work on harmonic convergence checking
 % sdtweb d_hbm SteppedSineDefFNL % xxx 
 % sdtweb cbi20b scriptproto


else; error('Script%s',CAM);
end
elseif comstr(Cam,'nmap'); [CAM,Cam]=comstr(CAM,5);
%% #nmap list of named experiments used for demos and non-regression tests ---

nmap=vhandle.nmap;
%% #SDT-contact_two_cube test cases -2

%% #CtcCube.A : load pressure exponential contact -3
li={'MeshCfg{"d_contact(cube)::n3e13{Kc1e12}"}',';', ...
     'SimuCfg{"Static{1e-8,chandle1}","Imp{100u,.1,chandle1,acall.}","EigOpt{5 5 0}"}', ...
     ';','RunCfg{nl_solve(Static),nlutil(HRbuild{q0m1}),run}'}';
nmap('CtcCube.A')=li; 


%% #CtcCube.B : load pressure and corner, exponential contact  -3
%    static followed by, hyperreduction
li={'MeshCfg{"d_contact(cube{loPC})::n3e13{Kc1e12,Lambda500}"}';
    ';';'SimuCfg{"Static{1e-8,chandle1}","Imp{50u,.1,chandle1,acall.}","EigOpt{5 5 0}"}';
    ';';'RunCfg{nl_solve(Static),nlutil(HRbuild{q0m1})}'};
nmap('CtcCube.B')=li; 

% CtcCube.C : exponential contact; static followed by, hyperreduction
% xxx add static load and point load 

%% #HE1 : hyperelastic test with one element -2
%  RO=struct('mat','simoA','Mesh','OneTrac','Case','DofSet:Sine{10}:C0{0}','NperPer',1e5,'Nper',3);RO.do='{run,va,pow}';dfr_ident('Load',RO);

li={'MeshCfg{"d_fetime(OneTrac{d2 2 2,MatSimoA}):Rivlin{-.2 -.2 -.4}"}';';' % RivlinCube experiment
      'SimuCfg{Imp{1m,3,chandle1,rt-1e-3},Sig{Tri(.1,/2,5)}}';';'
      'RunCfg{Time}'};
nmap('HE.1T')=li; 
li={'MeshCfg{"d_fetime(OneTrac):TopZa:SimoA"}';';' % RivlinCube experiment
      'SimuCfg{Exp{.2m,.2,chandle1},"SteppedSine{5}:C1(*100=%){60}"}';';'
      'RunCfg{Time}'};
nmap('HE.1Ta')=li; 

%% deal with outputs 
if comstr(Cam,'range')
  st=horzcat(li{:});
  fprintf('Running experiment\n %s\n',st)
  out=sdtm.range(struct,st);
elseif nargin==1; out=nmap;
else
  out=nmap(varargin{2});
end

elseif comstr(Cam,'range'); [CAM,Cam]=comstr(CAM,6);
%% #Range range_set_params ---------------------------------
  dbstack; keyboard;

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
  
  if isempty(Cam)&&isfield(val,'NeedInit')
    %% #NeedInit : actually fill -2
    [out,out1]=d_tdoe(val.type,Range,val.NeedInit);
    % out1 should contain the selected values
    return;
        
  elseif ~isempty(regexpi(Cam,'^m_','once'))
    %% #RangeM_ : MeshCfg parameters
    if iscell(val)
      out=struct('val',[],'lab',{{CAM}}, ...
        'param',struct(CAM,struct('type','pop','level',10)));
      [out.param.(CAM),out.val]=feval(fe_range('@popMerge'),out,CAM,val);
    else
      out=struct('val',val(:),'lab',{{CAM}}, ...
        'param',struct(CAM,struct('type','double','level',10)));
    end
  elseif ~isempty(regexpi(Cam,'^r_','once'))
    %% #RangeR_ : RunCfg parameters
    if iscell(val)
      out=struct('val',[],'lab',{{CAM}}, ...
        'param',struct(CAM,struct('type','pop','level',30)));
      [out.param.(CAM),out.val]=feval(fe_range('@popMerge'),out,CAM,val);
    else
      out=struct('val',val(:),'lab',{{CAM}}, ...
        'param',struct(CAM,struct('type','double','level',30)));
    end
    
  elseif ~isempty(regexp(Cam,'cfg$','once'))
    %% #RangeTypeCfg : support choices of
    if ischar(val); val={val}; end
    if iscell(val);
      [speed,r2]=ismember(lower(val),lower(Range.param.(CAM).choices));
      if ~all(speed); error('%s not found in %s',comstr(val(~speed),-30),CAM);end
    else; r2=val(:);
    end
    out=struct('val',r2(:),'lab',{{CAM}});
    
    
  elseif comstr(Cam,'check')||isempty(Cam)
    %% #Range or #RangeCheck : verify that all things needed are there
    if comstr(Cam,'check');[CAM,Cam]=comstr(CAM,6);end
    
    if carg<=nargin && isstruct(varargin{carg}); RO=varargin{carg};carg=carg+1;
    else; RO=struct();
    end
    if ~isfield(RO,'UsedField');RO.UsedField={};end
    st={'MeshCfg','SimuCfg','RunCfg'};
    if isempty(Range);
      Range=struct('param',[],'FileName',{{'ID'}},'val',[]);
      Range.param.MainFcn={'d_tdoe','range'}; % Allow callback for param setting
    end
    for j1=1:length(st) % Verify base params
      if isfield(Range.param,st{j1});
      elseif strcmpi(st{j1},'meshcfg');% Use reference implementation
        % d_mesh('MeshCfg',Range,'d_hbm(0D):DofSet:0dm1t');
        Range=d_mesh(st{j1},Range,val.(st{j1}));
      elseif strcmpi(st{j1},'simucfg');
       % Use reference implementation d_fetime
       % This updates the Range.param.MeshCfg.data{xxx} and Range.param.SimuCfg.data{1} 
        evt=sdth.sfield('addselected',struct('urn',val.(st{j1})), ...
            val,{'NperPer','iteStab','PerFrac','freq','tend','dt'});% xxxgetFieldList
        if isempty(evt.urn);
           fprintf('Skipping empty SimuCfg\n');
           if ~isfield(Range.param,'UsedField');Range.param.UsedField={};end
           Range.param.UsedField{end+1}='SimuCfg'; 
           continue;
        else
         Range=d_fetime(st{j1},Range,evt);
        end
      elseif isfield(val,st{j1});Range=d_tdoe(st{j1},Range,val.(st{j1}));
      else;Range=d_tdoe(st{j1},Range);
      end
      val.(st{j1})=Range.param.(st{j1}).choices;
    end
    if isstruct(val); % Allow automated fill from struct
      val.param=Range.param; val.FileName=Range.FileName;
      if isfield(Range.param,'UsedField')
       val=feutil('rmfield',val,Range.param.UsedField);
      end
      if isempty(CAM);
        Range=fe_range('BuildGrid -FileName',val,RO);
      else;%Range=d_tdoe('RangeCheckBuildSimple -FileName',RB);
        Range=fe_range(CAM,val,RO);
      end
    end
    Range.param=sdth.sfield('rmfield',Range.param,'UsedField');
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

  
elseif comstr(Cam,'solve'); [CAM,Cam]=comstr(CAM,6);
%% #Solve : parametric study  ---------------------------------------
if carg>nargin;RO=struct('info',1);else;RO=varargin{carg};carg=carg+1;end

if isempty(Cam)
 %  #Solve.Loop_
 % sdtweb d_mesh MeshCfg
 % sdtweb d_fetime SimuCfg
 % sdtweb d_tdoe RunCfg
 if isstruct(RO); 
  Range=d_tdoe('Range',RO); 
 end
 out=fe_range('Loop',Range,struct('ifFail','error'));
else
 
end
elseif nargin==3&&strcmpi(varargin{3},'steprun')
%% #stepRun  ---------------------------------------

Range=varargin{1};carg=2;evt=varargin{carg};carg=carg+1;
if isfield(Range,'Node')||isempty(Range)
 mo1=Range; Range=[]; RO=evt; % d_tdoe(mo1,RO,'steprun')
 S=RO.S;
 if ~isfield(evt,'RL');evt.RL.ifFail='error';end
else;
 RO=fe_range('ValEvtMerge',Range,evt); 
 mo1=Range.Res{1}; 
 if ~isfield(mo1,'nmap')&&isfield(Range,'nmap');mo1.nmap=Range.nmap;end
 S=sdth.findobj('_sub:',RO.urn);S=S.subs;if ischar(S);S={S};end
end
if ~isa(mo1.nmap,'vhandle.nmap'); error('Not an expected case');end
nmap=mo1.nmap; 
 %% #Solve.loop_level_30 do :  -3
 % earlier implementation see sdtweb d_tdoe LoadRun.Run
 for j1=1:length(S);js=j1; 
  if ischar(S{j1});[CAM,Cam]=comstr(S{j1},1);
  else;CAM=S{j1};Cam='cb';
  end
  switch Cam
  case {'time','run'}
 %% #stepRun.Time  -3
  if ~exist('mo1','var')&&isKey(nmap,'CurModel');mo1=nmap('CurModel');end
  op1=stack_get(mo1,'','TimeOpt','g');
  if isempty(op1)||strncmpi(Cam,'sta',3); op1=stack_get(mo1,'','TimeOptStat','g');end
  if ~isempty(stack_get(op1,'','Range'));op1.FinalCleanupFcn='';end
  % Actually run simulation
  [d1,mo1b]=fe_time(stack_set(mo1,'info','TimeOpt',op1));
  if iscell(d1); d1(cellfun(@isempty,d1(:,3)),:)=[];end
  if iscell(d1)&&size(d1,1)==1;d1=d1{1,3};
    d1=stack_set(d1,stack_get(op1,'','Range'));
  end
  i1=[]; if isfield(d1,'Y');i1=any(~isfinite(d1.Y),[1 3]); end
  if any(i1); 
    fprintf('Diverged for %s\n',comstr(d1.X{2}(i1)',-30));
    Range=stack_get(op1,'','Range','g');
    i1=find(squeeze(~any(~isfinite(d1.Y),[1 2])));Range.val(i1,:)
  end
  nmap('CurModel')=mo1b; nmap('CurTime')=d1;
  %eval(iigui({'d1','mo1b'},'SetInCallerC')) % set with comments
  otherwise
  if strncmpi(Cam,'dfrf',4)
 %% #stepRun.Dfrf  -3
   if ~isempty(stack_get(mo1,'','oProp'))
   elseif ~exist('mklserv_utils','file')|| ~exist('mklserv_client','file')
       warning('Consider installing mklserv_utils')
   else % Default assume complex sym
     mo1=stack_set(mo1,'info','oProp',mklserv_utils('opropCpxSym'));
   end
   d1=fe_simul(CAM,mo1);
   nmap('CurModel')=mo1; nmap('CurFreq')=d1;
  elseif strncmpi(Cam,'save',4)
   %% #stepRun.Save: intermediate save -3
   try
    RunRes=sdth.PARAM('RunRes');
    sdtm.save(RO.fStep{j1},'RunRes');
   catch err
    sdtw('_nb','Failed save step with %s',err.message);
   end
  else
   %% #stepRun.stepCb -3
   EndEval='';  stRes='RunRes';
   if isempty(Cam)
   elseif strcmpi(evt.RL.ifFail,'error')% Fail with error
     if ischar(CAM)&&~isempty(regexp(CAM,'[^\()]*=','once')); eval(CAM);
     else
      st=sdtm.urnCb(CAM); ans='';
      feval(st{:});  % Attempt to run a step d_shm@va/d_shm(va)
      if ~isempty(ans); sdth.PARAM(stRes,ans);
       mo1.nmap(stRes)=ans;
      end
     end
   else % Attempt try /catch
   try; 
     %% gui21@PostPro : allow callback to user function 
     if regexp(CAM,'[^\()]*='); eval(CAM);
     else
      st=sdtm.urnCb(Cam);
      feval(st{:})
     end
   catch e
     mo1.nmap('CurModel')=mo1; 
     error('%s not implemented',CAM);
   end
   end
   if ~isempty(EndEval);eval([EndEval ';']);end
  end
  end
 end


elseif comstr(Cam,'runcfg'); [CAM,Cam]=comstr(CAM,7);
%% #RunCfg : define runcfg  ---------------------------------------

Range=varargin{carg};carg=carg+1;
RO=varargin{carg};carg=carg+1;
if ischar(RO);RO=struct('urn',RO);end
if ~isfield(Range,'param');Range.param=struct;end
if ~isfield(Range.param,'RunCfg');
  Range.param.RunCfg=struct('type','pop','value',1,'level',30,...
      'choices',{{}},'data',{{}},'SetFcn',{{@d_tdoe,'stepRun'}},...
      'ShortFmt',1, ...
      'RepList', ... % 'TestEvtData',1,
      {{'R_(.*)_([^_]+)', ... % start with M_ and end with non_ parameter name
      struct('type',{'.','{}'},'subs',{'list',{'@token{1}',3}})
      }});
end
Range.param.RunCfg=feval(fe_range('@popMerge'),Range,'RunCfg',{RO.urn,RO});
if ~isfield(Range,'FileName')||isempty(Range.FileName);Range.FileName={'@RunCfg'};
elseif ~any(strcmpi(Range.FileName,'@RunCfg'));Range.FileName{end+1}='@RunCfg';
end


out=Range;
%% clean end
elseif comstr(Cam,'@');out=eval(CAM);
elseif comstr(Cam,'cvs'); out=sdtcheck('revision','$Revision: 7e3766c $  $Date: 2022-04-04 11:18:19 +0200 $');
else; error('%s unknown',CAM);
end 
end

function TgtMdl(varargin)
 %% #TgtMdl return tangent model as superelement (going through fe_time)
 mo1=evalin('caller','mo1');
 if evalin('caller','exist(''d1'',''var'')');d1=evalin('caller','d1');end
 
 if ~isfield(mo1,'NL')&&~isfield(mo1,'Case')
   opt=stack_get(mo1,'info','TimeOpt','g');
   if isempty(opt); opt=d_fetime('TimeOpt');end
   opt.Method='back'; 
   [mo1,C1,opt,z]=fe_time(opt,mo1);
   mo1.Case=C1;
 end
 [mo1,C1]=nl_spring('NLJacobianUpdate-tangentmdl',mo1,struct('u',0,'v',0));
 
 dbstack; keyboard; 
end
