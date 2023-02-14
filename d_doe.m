function [out,out1,out2]=d_doe(varargin);

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

li={'MeshCfg{d_hbm(0D):DofSet:0dm1t}';';'
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

 % xxx not yet functional see first example in cbi20b xxx 
li={'MeshCfg{d_hbm(Duffing2Dof),,CubFu}';';' 
     'SimuCfg{RO{NperPer2e3,Nper1,Methodnl_solve ModalNewmark},"SteppedSine{@ll(1,10,10)}:C1(*.001=N){2.5,10}"}';';'
     'RunCfg{Reduce,run,gui21@postproto}'};
 RT=struct('nmap',vhandle.nmap);
 RT.nmap('Reduce')='nl_solve(ReducFree 2 10 0 -float2 -SE)';
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

 key=''; if nargin>1; key=varargin{2};end
 if nargin>2;uo=varargin{3};carg=4; 
   if isfield(uo,'nmap');nmap=uo.nmap;else; nmap=vhandle.nmap;uo.nmap=nmap;end
   if isfield(uo,'Daq'); nmap('AcqTime')=uo.Daq.AcqTime;end
 else; nmap=vhandle.nmap;
 end
 if ~any(key=='{')
 elseif strncmp(key,'@',1);out=sdtm.urnCb(key);return; % @gui21{nmap,DaqShakerPres}
 else;[key,nmap]=sdtm.keyRep(nmap,key,'_ParSet');%sdtsys('nmap','VtDaq{dt,1e-6}')
 end
 if ~isKey(nmap,'AcqTime');nmap('AcqTime')=.3;end

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

%% #HE : hyperelastic testing -2
%% #HE.1T : hyperelastic test with one element -3
%  RO=struct('mat','simoA','Mesh','OneTrac','Case','DofSet:Sine{10}:C0{0}','NperPer',1e5,'Nper',3);RO.do='{run,va,pow}';dfr_ident('Load',RO);

li={'MeshCfg{"d_fetime(OneTrac{d2 2 2,MatSimoA}):Rivlin{-.2 -.2 -.4}"}';';' % RivlinCube experiment
      'SimuCfg{Imp{1m,3,chandle1,rt-1e-3},Sig{Tri(.1,/2,5)}}';';'
      'RunCfg{Time}'};
nmap('HE.1T')=li; 
li={'MeshCfg{"d_fetime(OneTrac):TopZa:SimoA"}';';' % RivlinCube experiment
      'SimuCfg{Exp{.2m,.2,chandle1},"SteppedSine{5}:C1(*100=%){60}"}';';'
      'RunCfg{Time}'};
nmap('HE.1Ta')=li; 

%% HE.relax : triangular with steps(Rep,Height,RaiseTime)
% sdtsys('urnsig','dt1m:Step{1,.1 .02,re5}:Step{-1,.1 .02,re5}');figure(1);plot(ans.Y)
% Think possibly Urnsig Table istair, iramp(v), 

%% #HBM : harmonic balance testing -2
nmap('Hbm.DoRange')={'nl_solve@DoTimeRangeb','ok';
     'nl_solve@chQueueSrc',struct;'nl_solve@chGetBuffer',struct
     'd_hbm@FirstStab','{outH 1:3,condd_hbm@ShootCheck}'};
nmap('Hbm.ExpList')={ ...
    ['SimuCfg{RO{NperPer200,Nper3,Methodnl_solve ModalNewmark},' ...
      '"SteppedSine{@ll(10,120,50)}:C1(*.001=N){1}"}'];
    ';'
    'RunCfg{Time,d_hbm@viewHarm,SetCI}'}; % Time{Profile}

%% #Hbm.OneDof nmap and list for reduced one DOF [RT,li]=d_doe('nmap','Hbm.OneDofRed'); -3
RT=struct('nmap',vhandle.nmap);
if ~isKey(nmap,'zeta'); nmap('zeta')=1e-2;end
RT.nmap('Reduce')='nl_solve(ReducFree 2 10 0 -float2 -SE)';
RT.nmap('SetCI')='ci=iiplot;cingui(''plotwd'',ci,''@OsDic(SDT Root)'',{''FnI'',''ImSw80'',''WrW49c''});;';
li={'MeshCfg{d_fetime(1DOF):MaxwellA{F2}}';';'
     'SimuCfg{ModalNewmark{1m,.1,fc,chandle1}}';';'
     'RunCfg{Reduce}'};
RT.nmap('CurExp')=li;nmap('Hbm.OneDofRed')=RT;

%% #Hbm.Gart : transient of Garteur testbed -3
RT=struct('nmap',vhandle.nmap);
if ~isKey(nmap,'NM'); nmap('NM')=10;end
if ~isKey(nmap,'dt'); nmap('dt')=1e-3;end
RT.nmap('NM')=nmap('NM'); RT.nmap('dt')=nmap('dt'); 
RT.nmap('Reduce')='nl_solve(ReducFree 2 $NM$ 1e3 -Float2 -SetDiag -SE)';
li={'MeshCfg{d_fetime(Gart):VtGart}';';' % Model Gart (2 DofLoad) sdtweb d_fetime MeshGart 
     % Case VtGart (defines Act:In1 and Act:In2 : 2 DofLoad combinaisons = 2 impacts) sdtweb d_fetime VtGart 
     'SimuCfg{ModalNewmark{$dt$,10,fc,chandle1}}';';'
     'RunCfg{Reduce}'};
RT.nmap('CurExp')=li;nmap('Hbm.Gart')=RT;
%% #Hbm.Duff : tabular with cubic non-linearity ready for VirtualTest -3
RT=struct('nmap',vhandle.nmap);
RB=struct('spec','BufTime 20 Overlap .75 fmin0 fmax60 -window hanning','ci',3);
RT.nmap('PostA')={'ExitFcn','Tip', ...
      struct('FinalCleanup',{{'nl_solve','PostCdof'}},'DOF',2.03,'DoFreq',RB)
      };
RT.nmap('ShowSpectro')='fe_simul(''fe_timeCleanup'')';
RT.nmap('Reduce')='nl_solve(ReducFree 2 10 0 -float2 -SE)';
RT.nmap('Transient')='nmap(''CurTime'')=fe_time(nmap(''CurModel''));';
RT.nmap('SetCI')='ci=iiplot;cingui(''plotwd'',ci,''@OsDic(SDT Root)'',{''FnI'',''ImSw80'',''WrW49c''});;';
li={'MeshCfg{"d_fetime(1DOF):MaxwellA{Z.01,Fds}:Duff"}',';', ...
      'SimuCfg{ModalNewmark{1m,$AcqTime$,sPostA},Sig{cnInput,Table(0 1n .1,1 0 0,istair)}}',';', ... % Sig=step
      'RunCfg{Reduce,Transient}'};
li=sdtm.keyRep(nmap,li);
RT.nmap('CurExp')=li;RT.tooltip='Single mass with cubic spring'; 
nmap('Hbm.Duff')=RT;


%% #Hbm.Fu : piecewise linear currently used in t_nlspring -3
RT=struct('nmap',vhandle.nmap);
RT.nmap('PostA')={'ExitFcn','Tip', ...
      struct('FinalCleanup',{{'nl_solve','PostCdof'}},'DOF',2.03,'DoFreq',RB)
      };
RT.nmap('ShowSpectro')='fe_simul(''fe_timeCleanup'')';
RT.nmap('Reduce')='nl_solve(ReducFree 2 10 0 -float2 -SE)';
RT.nmap('Transient')='nmap(''CurTime'')=fe_time(nmap(''CurModel''));';
RT.nmap('SetCI')='ci=iiplot;cingui(''plotwd'',ci,''@OsDic(SDT Root)'',{''FnI'',''ImSw80'',''WrW49c''});;';
li={'MeshCfg{"d_fetime(1DOF):MaxwellA{Z0}:FuA"}',';', ...
      'SimuCfg{ModalNewmark{1m,100,sPostA,fcShowSpectro}}',';', ...
      'RunCfg{Reduce,Transient,SetCI}'};
RT.nmap('CurExp')=li;
RT.tooltip='Single mass with gap contact'; 
nmap('Hbm.Fu')=RT;

%% #Hbm.SqBase : button groups for SqSig tab -3

 nmap('Hbm.SqBase')=struct('blocklist',{{ ...
   '@DefBut.Daq.Main';{'AcqSet',''}; '@DefBut.Daq.AcqSet';
   {'procSig','FRF Estimation'}; '@DefBut.Daq.procSig.FRFestimation';
   {'DoE','Not set'}}} ...
  );

%% #TV : time varying system tests -2

RT=struct('nmap',vhandle.nmap);
RT.nmap('PostA')='nl_solve@doFreq{spec{BufTime 2 Overlap .90 Fmax 500 -window hanning},ci3,nameAR}';
%RB=struct('spec','BufTime 2 Overlap .90 fmax500 -window hanning','ci',3);
%RT.nmap('PostA')={'FinalCleanupFcn','Tip', ...
%      struct('Cb',{{'nl_solve','PostCdof',struct('DOF',2.03,'DoFreq',RB,'name','AR')}})
%      };
li={'MeshCfg{"d_fetime(1DOF):ARTV{z.01}"}'; ';'
     'SimuCfg{ModalNewmark{.2m,50}}';';';'RunCfg{run,PostA}'};
%     'SimuCfg{ModalNewmark{.2m,50,sPostA}}';';';'RunCfg{run}'};
RT.nmap('CurExp')=li;
nmap('TV.AR')=RT;

%% #Tv.Hoffman with time varying stiffness -3
  RT=struct('nmap',vhandle.nmap);
  %RT.nmap('PostA')='nl_solve@doFreq{spec{BufTime 20 Overlap .90 Fmax 50 -window hanning},ci3}';
  RT.nmap('PostA')='d_squeal(ViewSpec{BufTime 50 Tmin 50 Overlap .90 Fmax 20 -window hanning},nameHoffKmuV)';
  RT.nmap('PostB')='d_contact@autoCycle{tclip50 20,dmBand1.2,ci3}';
  RT.nmap('PostInit')='d_squeal(LoadTime{ci[2 13]},$nmap)';

  li={'MeshCfg{d_contact(Hoffmann),TV,KmuV}';';'  % Mesh:Case:NL
   'SimuCfg{ModalNewmark{1m,400,uva111,rt-1e-4}SQ0{vq1Amp__5}}';';'
   'RunCfg{Time,PostB}'};
  RT.nmap('CurExp')=li;
nmap('TV.Hoff')=RT;


%% deal with outputs 
if comstr(Cam,'range')
   if nargin==2&&isKey(nmap,varargin{2}) % r1= cbi21('nmapRange','SoKBenchB_static')
    r1=nmap(varargin{2});li=r1.nmap('CurExp');
    fprintf('Running %s\n %s',varargin{2},sdtm.toString(li));
    out=sdtm.range(r1,li);
    if nargout>1; out1=out.nmap('CurModel');end
    if nargout>2; out2=out.nmap('CurTime');end
    
   else % mo2=daqsdt('nmaprange',{'SensA',';','RunVa'})
    out=sdtm.range(struct('nmap',nmap,'silent',-1),varargin{2:end});
   end

  %st=horzcat(li{:});
  %fprintf('Running experiment\n %s\n',st)
  %out=sdtm.range(struct,st);
elseif nargin==1; out=nmap;
elseif ~isempty(key)
  out=nmap(key);
  if nargout==2&&iscell(out)&&numel(out)==2;out1=out{2};out=out{1};end
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

error('Moved to sdtsys')


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
