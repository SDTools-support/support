function [out,out1,out2]=d_doe(varargin);

% D_DOE sample files for design of numerical experiments
%
%   Use d_doe('nmap') to see list  

%       Etienne Balmes, Guilherme Malacrida Alves
%       Copyright (c) 1990-2023 by SDTools, All Rights Reserved.
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
%RP=struct('MeshCfg','d_hbm(0D),DofSet,0dm1t','SimuCfg','SteppedSine{.5,1,10}:C0{0,15}:C1{2.5,10}', ...
%         'NperPer',2e3,'Nper',1,'iteStab',20,'RunCfg','{run,dfr_ident@va}','PerFrac',.2);
%d_tdoe('Solve',RP)

li={'MeshCfg{d_hbm(0D),DofSet,0dm1t}';
      'SimuCfg{RO{NperPer2e3,Nper1,iteStab20},"SteppedSine{5}:C0{0,15}:C1{2.5,10}"}';
      'RunCfg{Time,dfr_ident@va}'};
RT=struct('nmap',vhandle.nmap);RT.nmap('CurExp')=li;
r2=sdtm.range(RT);mo1b=r2.nmap('CurModel');%d2=mo2.nmap('CurTime');

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

 li={'MeshCfg{d_hbm(Duffing2Dof),,CubFu}';
     'SimuCfg{RO{NperPer2e3,Nper1},"SteppedSine{1}:C1(*.001=N){2.5,10}"}';
     'RunCfg{run,gui21@postproto}'};
 RT=struct('nmap',vhandle.nmap);
 RT.nmap('FirstStab')=struct('FinalCleanup','d_hbm@FirstStab','ite',10,'cond','d_hbm@CountIte');
 r2=sdtm.range(RT,li);%d2=mo2.nmap('CurTime');

 % xxx not yet functional see first example in cbi20b xxx 
li={'MeshCfg{d_hbm(Duffing2Dof),,CubFu}';
     'SimuCfg{RO{NperPer2e3,Nper1,Methodnl_solve ModalNewmark},"SteppedSine{@ll(1,10,10)}:C1(*.001=N){2.5,10}"}';
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
elseif comstr(Cam,'nmap');
%% #nmap list of named experiments used for demos and non-regression tests ---

[key,nmap,uo,carg]=sdtm.stdNmapArgs(varargin,CAM,carg);
%  key=''; if nargin>1; key=varargin{2};end
%  if nargin>2;uo=varargin{3};carg=4; 
%    if isfield(uo,'nmap');nmap=uo.nmap;else; nmap=vhandle.nmap;uo.nmap=nmap;end
%    if isfield(uo,'Daq'); nmap('AcqTime')=uo.Daq.AcqTime;end
%  else; nmap=vhandle.nmap;
%  end
%  if ~any(key=='{')
%  elseif strncmp(key,'@',1);out=sdtm.urnCb(key);return; % @gui21{nmap,DaqShakerPres}
%  else;[key,nmap]=sdtm.keyRep(nmap,key,'_ParSet');%sdtsys('nmap','VtDaq{dt,1e-6}')
%  end

 if ~isKey(nmap,'AcqTime');nmap('AcqTime')=.3;end

%% #SDT-contact CTC test cases -2

RT=struct('nmap',vhandle.nmap);
RT.nmap('SimBack')='SimuCfg{back{.2m,50,chandle1}}';

%% #CtcCube.A : load pressure exponential contact -3
li={'MeshCfg{d_contact(cube),,n3s13{Kc1e12}}', ...
     'SimuCfg{"Static{1e-8,chandle1}","Imp{100u,.1,chandle1,acall.}","EigOpt{5 5 0}"}', ...
     'RunCfg{SetM{q0/out,CurModel/out1},nl_solve(Static),SetM{CurModel/out},nlutil(HRbuild{q0m1}),run}'}';
RT.nmap('CubeA')=li; % CtcCube.A

%% #CtcCubeJ24 : developments for jacobian building tests -3
li={'MeshCfg{d_contact(cube),None{Traj{selSetNametop1,curveDownForward,o 1 0 0,storeprojM{q0/d2}}},n3s13{Kc1e12}}', ...
     'SimuCfg{"Imp{100u,.1,chandle1,acall.}","EigOpt{5 5 0}"}', ...
     'RunCfg{SetM{CurModel/model},nl_solve(fe_timeBack),SetM{CurModel/out},nlutil(HRbuild{q0m1}),run}'}';
t=(1:10)';
RT.nmap('DownForward')=struct('X',{{t,{'x';'z'}}},'Y',linspace(0,1,length(t))'*[1 -.1], ...
    'LabFcn','sprintf(''x%.4f z%.4f'',def.data(ch,2:3))','fun',[0 4]);

RT.nmap('CubeJ24')=li; % CtcCube.J24
% xxx Forced Motion 

%% #CtcCube.B : load pressure and corner, exponential contact  -3
%    static followed by, hyperreduction
li={'MeshCfg{d_contact(cube{loPC}),,n3s13{Kc1e12,Lambda500}}';
    'SimuCfg{"Static{1e-8,chandle1}","Imp{50u,.1,chandle1,acall.}","EigOpt{5 5 0}"}';
    'RunCfg{SetM{q0/out,CurModel/out1},nl_solve(Static),SetM{CurModel/out},nlutil(HRbuild{q0m1})}'};
RT.nmap('CubeB')=li; % CtcCube.B

  %RT.nmap('PostA')='nl_solve@doFreq{spec{BufTime 20 Overlap .90 Fmax 50 -window hanning},ci3}';
%  RT.nmap('PostA')='d_squeal(ViewSpec{BufTime 50 Tmin 50 Overlap .90 Fmax 20 -window hanning},nameHoffKmuV)';
%  RT.nmap('PostB')='d_contact@autoCycle{tclip50 20,dmBand1.2,ci3}';
%  RT.nmap('PostInit')='d_squeal(LoadTime{ci[2 13]},$nmap)';
% CtcCube.C : exponential contact; static followed by, hyperreduction
% xxx add static load and point load 

%% #CtcCub.Sclda : large motion checks -3
  % sdtweb d_contact meshScldCube
li={'MeshCfg{d_contact(ScldCube),None{TrajScld}}','SimBack','RunCfg{Time}'};
% RT.nmap('StickMesh')={};
RT.nmap('TrajScld')={ ... % sdtweb _bp d_contact caseTraj 
    'Traj{selSetNameWheel,curveDownForward,o 1 0 0,storemodel{DownForward/d2}}'
    'CbRefWheel';'CbStickWheel';'CbRefRail';'CbCtcGen'
    };
RT.nmap('TrajScldB')={ ...
    'Case{reset}'
    'Case{FixDof,Base,"inelt{innode{setname"Rail"}&selface&facing> .9 0 0 -1000}"}'
    'Case{DofSet,Top,"rb{inelt{proid8&selface&innode{z>.14}&facing> .9 0 0 1000},dir 1 3,curveDownForward,KeepDof}"}'
    'Case{Pcond,Scld,p_contact(''PcondScld'')}'
    'CbRefWheel';'CbRefRail';'CbStickWheel';'CbCtcGen'
    };
RT.nmap('TrajScldC')=struct('ToolTip','Cube with target Lc', ...
    'li',{{ ...
    'Case{reset}'
    'Case{FixDof,Base,"inelt{innode{setname"Rail"}&selface&facing> .9 0 0 -1000}"}'
    'Case{DofSet,Top,"rb{inelt{proid8&selface&innode{z>.14}&facing> .9 0 0 1000},dir 1 3,curveDownForward,KeepDof}"}'
    'Case{Pcond,Scld,p_contact(''PcondScld'')}'
    'InitQ0{elem0(vect{x0,y0,z-.01,"selwithnode{setNameWheel}"})}'
    'CbRefWheelLC';'CbRefRailLC';'CbStickWheel';'CbCtcGen'
    }});
RT.nmap('TrajScldD')=struct('ToolTip','Cube with target Lc and force', ...
    'li',{{ ...
    'Case{reset}'
    'Case{FixDof,Base,"inelt{innode{setname"Rail"}&selface&innode{z<.01}&facing> .9 0 0 -1000}"}'
    'Case{FixDof,Topd,"inelt{innode{setname"Wheel"}&selface&innode{z>.14}&facing> .9 0 0 1000} -DOF2"}'
    'Case{DofLoad,Top,"rb{inelt{proid8&selface&facing> .9 0 0 1000},dir 1 3,curveDownForward,KeepDof}"}'
    'Case{Pcond,Scld,p_contact(''PcondScld'')}'
    'CbRefWheelLC';'CbRefRailLC';'CbStickWheel';'CbCtcGen'
    }});
RT.nmap('ScLdA')=li; % ScLdA

l=0.0750; 
%% Wheel refinement selection
RT.nmap('CbRefWheel')={@fe_shapeoptim,'RefineHexaMesh','$projM', ...
 sprintf('ProNameWheel&innode {distFcn"{box{%.15g %.15g %.15g, %.15g %.15g %.15g ,1 0 0,0 1 0, 0 0 -1}}"}',...
 l*[1.5 0.5 1.15 .3 .3 .2]) };
RT.nmap('CbRefWheelLC')={@fe_shapeoptim,'RefineHexaMesh -lc2.5e-3 -lcmin2e-3','$projM', ...
 sprintf('ProNameWheel&innode {distFcn"{box{%.15g %.15g %.15g, %.15g %.15g %.15g ,1 0 0,0 1 0, 0 0 -1}}"}',...
 l*[1.5 0.5 1.15 .3 .3 .2]) };

 ...l*[1.5 0.5 1.1 .5 .5 .25]) };
   %sprintf('ProNameWheel&withnode {distfcn"{sphere{%.15g %.15g %.15g, .03}}"}',...
    % l*[1.5 .5 1])};
%    sprintf('withnode{setname Wheel}&withnode{z<%g}',1.2*l)};
% Wheel surface sticking 
RS=struct('sel',[ 'inElt{ProNameWheel & selface&innode{z<' sprintf('%.5g',1.02*l) '}&facing >.9 0 0 -5000}'],'distFcn',[]);
RS.distFcn=lsutil('gen',[],{struct('shape','sphere','rc',l*5,'xc',l*1.5,'yc',l*.5,'zc',(5+1)*l)});
RT.nmap('CbStickWheel')={@lsutil,'SurfStick','$projM',RS};
RB=struct(...'RefRail', ...
   ...sprintf('ProNameRail&innode {distfcn"{cyl{%.15g %.15g %.15g,r%.15g,n1 0 0,z%.15g %.15g}}"}', ...
  ... l*[0 .5 1,.1, .7 4]), ...
   ...l*[0 .5 1,.27, .7 4]), ...
   'CtcRail', ...
  ...  sprintf('ProNameRail&selface&facing >.9 0 0 1e4&innode {distfcn"{cyl{%.15g %.15g %.15g,r%.15g,n1 0 0,z%.15g %.15g}}"}', ...
  ...   l*[0 .5 1,.4, .7 4]), ...
    sprintf(['ProNameRail&selface&facing >.9 0 0 1e4&innode {distfcn"{'...
    'box{%.15g %.15g %.15g,%.15g %.15g %.15g,%.15g %.15g %.15g,%.15g %.15g %.15g,%.15g %.15g %.15g}}"}'], ...
     l*[0 .5 1,4 .4 .05],[ 1 0 0,0 1 0,0 0 1]), ...
   'CtcWheel',['ProNameWheel&selface&withnode', ...
      sprintf('{distfcn"{sphere{%.15g %.15g %.15g, .015}}"}',l*[1.5 .5 1])]);
RB.top=[.6 2.5]*l; 
%% rail Refinement
RB.RefRail=sprintf('ProNameRail & innode {distFcn"{box{%.15g %.15g %.15g, %.15g %.15g %.15g ,1 0 0,0 1 0, 0 0 -1}}"}',...
 [mean(RB.top) 0.5*l 1*l, diff(RB.top)/2 .3*l .3*l]) ;
RT.nmap('CbRefRail')={@fe_shapeoptim,'RefineHexaMesh','$projM',RB.RefRail}; 
RT.nmap('CbRefRailLC')={@fe_shapeoptim,'RefineHexaMesh -lc2.5e-3 -lcmin2e-3','$projM',RB.RefRail}; 
% -interMPC problem with slaves in slave surface
 % sprintf('withnode{setname Rail}&withnode{z<%g}',1.2*l)};
% fecom('shownodemark','distFcn{sphere{0 0 0, .05}}')

RT.nmap('CbCtcGen')={@ctc_utils,'generatecontactpair','$projM$', ...
     struct('slave',RB.CtcRail,'master',RB.CtcWheel,'ProId',201, ...
      'InitUnl0','d_contact@surfStick','StoreType',3)};

% xxx Missing InitUnl0 
RT.nmap('ScldCS')=struct('ToolTip','Sphere/cube contact', ...
    'li',{{['MeshCfg{d_contact(ScldCube{Kc1e9,Integ2}),' ...
      'None{TrajScld}}']
    'SimuCfg{"Imp{100u,.1,chandle1,BetaR7e-6}"}'
    'RunCfg{Time}'}});
RT.nmap('RangeLoopOpt')= struct('RangeLoopResKeys',{{'CurModel','CurTime'}}, ...
    'ifFail','error');

% select current experiment base on n field 
if isKey(nmap,'n')&&isKey(RT.nmap,nmap('n')); RT.nmap('CurExp')=RT.nmap(nmap('n')); end
nmap('Ctc')=RT;

%% #HE : hyperelastic testing -2
%% #HE.1T : hyperelastic test with one element -3
%  RO=struct('mat','simoA','Mesh','OneTrac','Case','DofSet:Sine{10}:C0{0}','NperPer',1e5,'Nper',3);RO.do='{run,va,pow}';dfr_ident('Load',RO);

li={'MeshCfg{d_fetime(OneTrac{d2 2 2,MatSimoA}),Rivlin{-.2 -.2 -.4}}'; % RivlinCube experiment
      'SimuCfg{Imp{1m,3,chandle1,rt-1e-3},Sig{Tri(.1,/2,5)}}'
      'RunCfg{Time}'};
nmap('HE.1T')=li; 
li={'MeshCfg{d_fetime(OneTrac),TopZa,SimoA}' % RivlinCube experiment
      'SimuCfg{Exp{.2m,.2,chandle1},"SteppedSine{5}:C1(*100=%){60}"}'
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
    'RunCfg{Time,d_hbm@viewHarm,SetCI}'}; % Time{Profile}

%% #Hbm.OneDof nmap and list for reduced one DOF [RT,li]=d_doe('nmap','Hbm.OneDofRed'); -3
RT=struct('nmap',vhandle.nmap);
if ~isKey(nmap,'zeta'); nmap('zeta')=1e-2;end
RT.nmap('Reduce')={'nl_solve','$RO',[],'ReducFree 2 10 0 -float2 -SE'};
RT.nmap('SetCI')='ci=iiplot;cingui(''plotwd'',ci,''@OsDic(SDT Root)'',{''FnI'',''ImSw80'',''WrW49c''});;';
li={'MeshCfg{d_fetime(1DOF),MaxwellA{F2}}';
     'SimuCfg{ModalNewmark{1m,.1,fc,chandle1}}';
     'RunCfg{Reduce}'};
RT.nmap('CurExp')=li;nmap('Hbm.OneDofRed')=RT;

%% #Hbm.Gart : transient of Garteur testbed -3
RT=struct('nmap',vhandle.nmap);
if ~isKey(nmap,'NM'); nmap('NM')=10;end
if ~isKey(nmap,'dt'); nmap('dt')=1e-3;end
RT.nmap('NM')=nmap('NM'); RT.nmap('dt')=nmap('dt'); 
RT.nmap('Reduce')={'nl_solve','$RO',[],'ReducFree 2 $NM$ 1e3 -Float2 -SetDiag -SE'};
%RT.nmap('Reduce')='nl_solve(ReducFree 2 $NM$ 1e3 -Float2 -SetDiag -SE)';
li={'MeshCfg{d_fetime(Gart),VtGart}'; % Model Gart (2 DofLoad) sdtweb d_fetime MeshGart 
     % Case VtGart (defines Act:In1 and Act:In2 : 2 DofLoad combinaisons = 2 impacts) sdtweb d_fetime VtGart 
     'SimuCfg{ModalNewmark{$dt$,10,fc,chandle1}}';
     'RunCfg{Reduce}'};
RT.nmap('CurExp')=li;li{3}='RunCfg{Reduce,Time}';RT.nmap('ExpRun')=li;
nmap('Hbm.Gart')=RT;
%% #Hbm.Duff : tabular with cubic non-linearity ready for VirtualTest -3
RT=struct('nmap',vhandle.nmap);
RB=struct('spec','BufTime 20 Overlap .75 fmin0 fmax60 -window hanning','ci',3);
RT.nmap('PostA')={'ExitFcn','Tip', ...
      struct('FinalCleanup',{{'nl_solve','PostCdof'}},'DOF',2.03,'DoFreq',RB)
      };
RT.nmap('ShowSpectro')='fe_simul(''fe_timeCleanup'')';
RT.nmap('Reduce')='nl_solve(ReducFree 2 10 0 -float2 -SE,$RO)';
RT.nmap('Transient')='nmap(''CurTime'')=fe_time(nmap(''CurModel''));';
RT.nmap('SetCI')='ci=iiplot;cingui(''plotwd'',ci,''@OsDic(SDT Root)'',{''FnI'',''ImSw80'',''WrW49c''});;';
li={'MeshCfg{d_fetime(1DOF),MaxwellA{Z.01,Fds},Duff}', ...
      'SimuCfg{ModalNewmark{1m,$AcqTime$,sPostA},Sig{cnInput,Table(0 1n .1,1 0 0,istair)}}', ... % Sig=step
      'RunCfg{Reduce,Transient}'};
li=sdtm.keyRep(nmap,li);
RT.nmap('CurExp')=li;RT.ToolTip='Single mass with cubic spring'; 
nmap('Hbm.Duff')=RT;


%% #Hbm.Fu : piecewise linear currently used in t_nlspring -3
RT=struct('nmap',vhandle.nmap);
RT.nmap('PostA')={'ExitFcn','Tip', ...
      struct('FinalCleanup',{{'nl_solve','PostCdof'}},'DOF',2.03,'DoFreq',RB)
      };
RT.nmap('ShowSpectro')='fe_simul(fe_timeCleanup)';
RT.nmap('Reduce')='nl_solve(ReducFree 2 10 0 -float2 -SE,$RO)';
%RT.nmap('Transient')='nmap(''CurTime'')=fe_time($RO);';
RT.nmap('Transient')={'fe_time','CurModel'};
RT.nmap('SetCI')='ci=iiplot;cingui(''plotwd'',ci,''@OsDic(SDT Root)'',{''FnI'',''ImSw80'',''WrW49c''});;';
li={'MeshCfg{d_fetime(1DOF),MaxwellA{Z0},FuA}'
      'SimuCfg{ModalNewmark{1m,100,sPostA,fcShowSpectro}}'
      'RunCfg{Step{CurTime},Reduce,change{nmap(LastContinue):LastContinue},Transient,SetCI}'};
RT.nmap('CurExp')=li;
RT.ToolTip='Single mass with gap contact'; 
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
li={'MeshCfg{d_fetime(1DOF),ARTV{z.01}}';
     'SimuCfg{ModalNewmark{.2m,50}}';'RunCfg{Step{CurTime},change{nmap(LastContinue):LastContinue},run,PostA}'};
%     'SimuCfg{ModalNewmark{.2m,50,sPostA}}';;'RunCfg{run}'};
RT.nmap('CurExp')=li;RT.ToolTip='Auto-resonance+time variation example'; 
nmap('TV.AR')=RT;

%% #Tv.Hoffman{n,Mmuv} with time varying stiffness -3
% see also sdtweb d_contact('HoffMann')
  RT=struct('nmap',vhandle.nmap);
  %RT.nmap('PostA')='nl_solve@doFreq{spec{BufTime 20 Overlap .90 Fmax 50 -window hanning},ci3}';
  RT.nmap('PostA')='d_squeal(ViewSpec{BufTime 50 Tmin 50 Overlap .90 Fmax 20 -window hanning},nameHoffKmuV)';
  RT.nmap('PostB')='d_contact@autoCycle{tclip50 20,dmBand1.2,ci3}';
  RT.nmap('PostC1')='d_squeal(ViewSpec{BufTime 2 Tmin 5 Overlap .90 Fmax 50 -window hanning},nameHoffman)';
  RT.nmap('PostC2')='d_contact@autoCycle{tclip50 20,dmBand20,ci3}';
  RT.nmap('PostD1')='d_squeal(ViewSpec{BufTime .4 fmin 2700 3000 Overlap .90 -window hanning},nameHoffman)';
  RT.nmap('PostD2')='d_contact@autoCycle{tclip.01 .02,ifBand50,dmBand500,aeBand100,ci3}';
  RT.nmap('PostInit')='d_squeal(LoadTime{ci[2 13]},$nmap)';

  li={'MeshCfg{d_contact(Hoffmann),TV,KmuV}';  % Mesh:Case:NL
   'SimuCfg{ModalNewmark{1m,400,uva111,rt-1e-4}SQ0{vq1Amp__5}}';
   'RunCfg{Time,PostB}'};
  RT.nmap('TVK.MN')=li;  % time varying stiffness
  % time and amplitude varying stiffness, MSSP case
  li={'MeshCfg{d_contact(Hoffmann{kx3.55e+08,kz3.55e+08,ks4.4e+04,mu0.55,cx56.5,cz56.5}),TV{KA}}';
   'SimuCfg{ModalNewmark{40u,4,uva111,rt-1e-4}SQ0{vq1Amp__.0005}}';
   'RunCfg{Time,PostInit,PostD1,PostD2}'};
  RT.nmap('TVKA.MN')=li;  
  %% 
  li={'MeshCfg{d_contact(Hoffmann),,KmuV}';  % Mesh:Case:NL
  'SimuCfg{ModalNewmark{10m,400,uva111,rt-1e-4}SQ0{vq1Amp__10}}';
  'RunCfg{Time,PostInit,PostA,PostB}'};
  RT.nmap('Kmuv.MN')=li; % non-linear damping (change viscous damping based on velocity)
  %% 
  li={'MeshCfg{d_contact(Hoffmann),Ctc,TExp}'; % Mesh:Case:NL
  'SimuCfg{Implicit{10m,400,uva111,rt-1e-4}SQ0{vq1Amp__1}}';
  ...'SimuCfg{ModalNewmark{10m,400,uva111,rt-1e-4}SQ0{vq1Amp__10}}';
  'RunCfg{Time,PostInit,PostA,PostB}'};
  RT.nmap('Texp.Imp')=li; % non-linear damping (change viscous damping based on velocity)

  %% sdtweb d_contact meshhoff
  li={'MeshCfg{d_contact(Hoffmann{cx.01,cz.01,kx100,kz110,ks32}),Ctc{},{1001,TExpCh{ke.2},1,KmuV}}';  % Mesh:Case:NL
  'SimuCfg{ModalNewmark{.5m,40,uva111,rt-1e-4}SQ0{vq1Amp__10}}';
  'RunCfg{Time,PostInit,PostC1,PostC2}'};
  RT.nmap('Texp.MN')=li; % non-linear damping (change viscous damping based on velocity)

  %%
  li={'MeshCfg{d_contact(Hoffmann),,}';  % Mesh:Case:NL
   ...'SimuCfg{Implicit{10m,400,uva111,rt-1e-4}SQ0{vq1Amp__5}}';
   'SimuCfg{ModalNewmark{1m,100,uva111,rt-1e-4}SQ0{vq1Amp__5}}';
   'RunCfg{Time,PostInit,PostA}'};
  RT.nmap('base.MN')=li; % non-linear damping (change viscous damping based on velocity)

%% 
li={'Steel',[1 fe_mat('m_elastic','SI',1)  210e9 .3 7800]};
r2=vhandle.nmap(li,'Map:MatDb');
NL=struct('type','nl_inout','lab','uMaxw','MatTyp',{{1}},'keepLin',1, ...
   'StoreType',3, ... % Just along z
  'adofi',[],'MexCb',{{nlutil('@uMaxw')}}, ...
  'Fu',struct('Einf',1,'cells',[4 .3 1 0;4 .4 5 0],'X',{{[-1;1],{'Coef';'dCoef'}}}, ...
  'Xlab',{{'Strain','Comp'}},'Y',[1+1e-5 0;2 0]));
r2('HyperLVol')=struct('NLdata',NL,'il','setpro 111 in-3');

nmap('MatDb')=r2;



if isKey(nmap,'n')&&isKey(RT.nmap,nmap('n')); 
    RT.nmap('CurExp')=RT.nmap(nmap('n')); 
    if isKey(nmap,'s')&&isequal(nmap('s'),'back') % d_doe('nmap','TV.Hoff{n,Texp.MN,s,back}')
     RT.nmap('CurExp')=strrep(RT.nmap('CurExp'),'ModalNewmark','Back');
    end
end

nmap('TV.Hoff')=RT; % d_doe('nmap','TV.Hoff{n,Texp.MN}')

% d_doe('nmap','Ex.Hoff.KmuV.Imp')
%  sdtm.toString(li(1:2:end)');RT.nmap('CurExp')=li;

 %RT=struct('nmap',nmap);d_doe('nmap','TV.Hoff');RT.nmap('CurExp')=li;
 % sdtm.stdNmapOut('call')
 if nargout==1;out=sdtm.stdNmapOut(nmap,key,nargout,CAM);
 elseif nargout>1;[out,out1,out2]=sdtm.stdNmapOut(nmap,key,nargout,CAM);
 else; sdtm.stdNmapOut(nmap,key,nargout,CAM);
 end
 %% #Nmap.End

elseif comstr(Cam,'solve'); [CAM,Cam]=comstr(CAM,6);
%% Solve : parametric study  ---------------------------------------
if carg>nargin;RO=struct('info',1);else;RO=varargin{carg};carg=carg+1;end
sdtw('_ewt','should use sdtsys')
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
%% stepRun  ---------------------------------------
error('Moved to sdtsys')


elseif comstr(Cam,'runcfg'); [CAM,Cam]=comstr(CAM,7);
%% #RunCfg : define runcfg  ---------------------------------------

sdtw('_ewt','should use sdtsys')

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
