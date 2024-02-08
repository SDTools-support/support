function [out,out1,out2]=d_squeal(varargin)
 
% D_SQUEAL Support for demonstrations related to squeal simulations
% Requires SDT/FEMLink, SDT/Contact, SDT/Param
%
% Starter scripts
% d_squeal('ScriptSimple')
% d_squeal('ScriptFromABQ')
%
% Functionalities
% d_squeal('SolveModes')
% d_squeal('ViewStab')
% d_squeal('ViewShapes')
%
% sdtweb _taglist d_squeal
%

%       Guillaume Vermot des Roches, Etienne Balmes
%       Copyright (c) 1990-2023 by SDTools, All Rights Reserved.
%       For revision information use d_squeal('cvs')

if nargin==0
 if sdtdef('isinteractive')
  fprintf('You should call a specific command, see:\n')
  sdtweb('_link','sdtweb(''_taglist'',''d_squeal'')','See the list of commands')
  
 else % allow auto testing
  d_squeal('ScriptSimple')
  d_squeal('ScriptFromABQ')
  % need MVR example
  
 end
 return
end

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*DEFNU,*NOSEM>

if ~ischar(varargin{1})||isempty(varargin{1})
  obj=varargin{1};evt=varargin{2};[CAM,Cam]=comstr(varargin{3},1);carg=4;
else; [CAM,Cam]=comstr(varargin{1},1);carg=2;obj=[]; evt=[];
end

if comstr(Cam,'script'); [CAM,Cam]=comstr(CAM,7);
 %% #Script/#Tutos ------------------------------------------------------------------
 if comstr(Cam,'simple')
  %% ScriptSimple: concepts illustrated on a simple SDT model -2
  
  %% #TutoSimple : Command demonstration on a simple case --------------------
  
  %% Step 1 : Complex Eigenvalue Analysis
  % Generate a sample brake squeal model, with resolved static
  cf=d_contact('ScriptBrake');
  
  % Generate a tangent model, based on the static contact states provided
  mo1=ctc_utils('tangentmdl inModel-evalFNL',cf.mdl.GetData,cf.Stack{'q0'});
  
  % Do CEA with a packaged call to nl_solve complex mode solver
  def=d_squeal('SolveModes',mo1,[5 25 1e3]);
  
  % Display a stability diagram
  d_squeal('ViewStab-gf1-reset',def)
  
  %% Step 2 : Property handling
  % Remove friction damping, and recompute modes
  ctc_utils('set euler=-1',cf.mdl)
  mo1=ctc_utils('tangentmdl inModel-evalFNL',cf.mdl.GetData,cf.Stack{'q0'});
  def2=d_squeal('SolveModes',mo1,[5 25 1e3]);
  
  d_squeal('ViewStab-gf1',def2)
  d_squeal('ViewShapes',cf,def)
  
  % Assign gyroscopic coupling and recompute modes
  data={'setname discSet',5,[0 0 1]}; % disc element selection, Amplitude, axis
  % declare rotation data in the current model
  mo1=fe_cyclic('OmegaGroupSetFirst',stack_set(mo1,'info','OmegaData',data));
  
  % resolve modes with gyro effects
  def7=d_squeal('SolveModes',mo1,[5 25 1e3]);
  d_squeal('ViewStab-gf1',def7)
  
  % Vary Mu
  ctc_utils('set Fv"Mu=0.6"',cf.mdl)
  mo1=ctc_utils('tangentmdl inModel-evalFNL',cf.mdl.GetData,cf.Stack{'q0'});
  
  def8=d_squeal('SolveModes',mo1,[5 25 1e3]);
  d_squeal('ViewStab-gf1',def8)
  
  d_squeal('ViewStab-Table',def8)
  
  % Use dMu/dV
  ctc_utils('setedit Fv"dMudV=-5"',cf.mdl)
  mo1=ctc_utils('tangentmdl inModel-evalFNL',cf.mdl.GetData,cf.Stack{'q0'});
  
  def9=d_squeal('SolveModes',mo1,[5 25 1e3]);
  def.label='dMu/dV=0';def9.label='dMu/dV=-5';
  d_squeal('ViewStab-gf1-reset',def)
  d_squeal('ViewStab-gf1',def9)
  d_squeal('ViewShapes',cf,def9)
  
  %% EndTuto
  
 elseif comstr(Cam,'fromabq')
  %% ScriptFromABQ: squeal from model and static state imports (ABAQUS case) -2
  
  %% #TutoFromAbaqus : CEA from an ABAQUS model and steady sliding solution --
  wd=d_squeal('WdAbqDemo'); % provides and/or load the directory containing files to be imported
  
  %% Step 1 : Model import, integrated
  % First define the ABAQUS computation file
  finp=fullfile(wd,'brake_squeal.inp');
  clear sdtroot
  % Read with FEMLinkTab
  sdtroot('SetFEMlink',struct('FileName',finp,...
   'ImportType','Custom',...
   'BuildListGen',{{strrep(finp,'.inp','.fil')}},...
   'BuildStepGen','steplast',...
   'BuildCb','Code(BuildModel-contact -module -tgStickNoMotion -useRes)',...  
   'Import','do'));
  
  %% Step 2 : Model import, low level
  % Read model, resolve implicit data and add static results from .fil
  mo1=abaqus('read-resolve',finp,'buildup',strrep(finp,'.inp','.fil'));
  mo1.unit='MM';
  % Setup contact from ABAQUS data and generate a load case and static results
  % associated to a given static step (here use last static step)
  % contact is here built using the contact module
  % tgStickNoMotion defined a sticking state for contact pairs with no applied motion
  % useRes sets contact states to reproduce the nodal resultant observed on the contact
  % surfaces from the static deformation provided by ABAQUS
  mo2=abaqus('BuildModel-steplast -contact -module -tgStickNoMotion -useRes',mo1);
  
  %% Step 3 : contact display
  % Display the model
  cf=feplot(mo2);
  
  % Display contact fields form displayed model
  ctc_utils('showpn',cf,cf.Stack{'q0'})
  
  %% Step 4 : Generate a tangent model, based on the static contact states provided
  % inModel sets a SE named tgtctc in the output model
  % evalFNL asks to evaluate the contact states from the static deformation
  %    this is necessary here due to the fact that the data imported from 
  %    ABAQUS do not provide the states in a format exploitable by SDT
  mo3=ctc_utils('tangentmdl inModel-evalFNL',mo2,stack_get(mo2,'curve','q0','get'));
  
  %% Step 5 : Do CEA with a packaged call to nl_solve complex mode solver
  def=d_squeal('SolveModes',mo3,[5 100 1e3]);
  def.label='base';
  d_squeal('ViewStab-unst-gf1-reset',def)
  d_squeal('ViewStab-Table',def)
  
  %% Step 6 : Variations
  % Deactivate friction damping: euler option
  %   euler option has the role of *MOTION, set to 1 for rotation
  %     friction damping is used by default
  proid=ctc_utils('get CProSlide',mo2); % recover identifiers of contact properties with applied motion
  % apply the new euler option
  mo2=ctc_utils(sprintf('set proid %s euler=-1',num2str(proid(:)')),mo2);
  % generate a new tangent model associated to the new options
  mo4=ctc_utils('tangentmdl inModel-evalFNL',mo2,stack_get(mo2,'curve','q0','get'));
  % New modes
  def4=d_squeal('SolveModes',mo4,[5 100 1e3]);
  def4.label='nofd';
  d_squeal('ViewStab-unst-gf1',def4)

  % Declare gyroscopic effects on the disc
  % it is now possible to move the rotation origin in SDT
  %mo4.Node(:,5:7)=mo4.Node(:,5:7)-ones(size(mo4.Node,1),1)*[3.7438 -7.6851 32.508];
  data={'matid4',5,[0 0 1],[3.7438 -7.6851 32.508]}; % disc element selection, Amplitude, axis, origin
  % declare rotation data in the current model
  mo4=fe_cyclic('OmegaGroupSetFirst',stack_set(mo4,'info','OmegaData',data));
  
  % resolve modes with gyro effects
  def7=d_squeal('SolveModes',mo4,[5 100 1e3]);
  def7.label='gyro, nofd';
  d_squeal('ViewStab-unst-gf1',def7)
  
  % define a stiffness as 1 times the representative stiffness
  mo2=ctc_utils(sprintf('set proid %s Fu="Kc=-1"',num2str(proid(:)')),mo2);
  % generate a new tangent model associated to the new options
  mo4=ctc_utils('tangentmdl inModel-evalFNL',mo2,stack_get(mo2,'curve','q0','get'));
  % New modes
  def8=d_squeal('SolveModes',mo4,[5 100 1e3]);
  def8.label='nofd, softer';
  d_squeal('ViewStab-unst-gf1',def8)
  
  %% Step 7 : reporting
  comgui('ImWrite',1)
  d_squeal('ViewShapes-unst',cf,def7);
  
  if 1==2 % command to capture animations
   d_squeal('ViewShapes-Write-movie-unst',cf,def7)
  end

  %% EndTuto 

 elseif comstr(Cam,'padmat'); [CAM,Cam]=comstr(CAM,7);
  %% #ScriptPadMat: sample material sensitivity study on the pad material
  
  
model=d_contact('tutoSimpleBrakeLS -s11 mo3'); % MeshSimplebrakeLS{xxx} xxxGV 
model=feutilb('submodel',model,'matid 2 3 & withnode{z>0}',model);
if comstr(Cam,'bend')
 % clamp braking_FSurf area
 % add support on plate
 model=feutilb('submodel',model,'matid 2 & withnode{z>0}',model);
 n1=feutil('getnode inelt{setname PressureArea}',model);
 n2=feutil('getnode z==min(z)',model);
 model=fe_case(model,'FixDof','Base',n2(:,1),'FixDof','Top',n1(:,1))
else
 model=fe_case(model,'reset');
 model.name='simple pad';
end
model.pl(1,:)=0;
pl1=[2 fe_mat('m_elastic','TM',6) 7600 7600 2800 .1 .036842 .23 1000 1000 3800 2.62e-9]
model.pl(1,1:length(pl1))=pl1;

% parametrize each value, couple in plane
RA=struct('type','ortho','Group',...
 {{'Ezz',3; % name is standardized for CMT
 'Exy',[1 2];
 'Nuzxy',[4 5];
 'Nu12',6;
 'Gzxy',[7 8];
 'Gxy',9}},...
 'ParTyp',1,'MatId',2);
mo1=fevisco('MatSplit-appendPar',model,RA);


mo1=stack_set(mo1,'info','EigOpt',[5 15 1e3]);
Ra0={'lab"Ezz" cur1 min.45 max1.45 nPoints5000 scale"lin"';
'lab"Exy" cur1 min.45 max1.45 nPoints5000 scale"lin"';
'lab"Nuzxy" cur1 min.45 max1.45 nPoints5000 scale"lin"';
'lab"Nu12" cur1 min.45 max1.45 nPoints5000 scale"lin"';
'lab"Gzxy" cur1 min.45 max1.45 nPoints5000 scale"lin"';
'lab"Gxy" cur1 min.45 max1.45 nPoints5000 scale"lin"';}
Ra0=cellfun(@(x)fe_range('buildvect',x),Ra0,'uni',0);
mo1=stack_set(mo1,'info','Range0',Ra0(:)')

RO=struct('name',[mo1.name '_mat'],'RangeType','type"randperm"5000',...
 'ListEig','eig 2005 15 1e3','matdes',[2 1],...
 'NormOpt',[-1 0 0 1e-12],'ListR','MinMaxm-flip',...
 'skipNorm',1,'NormAtEnd',1,'norm',[2],'KdNom',1); % xxx NormAtEnd and omattypev6

sdtdef('omattype','v6')
sdtdef('EigOOC',30)

MVR=d_squeal('SolveMVR',mo1,RO)

tic; [hist,dr]=fe2xf('PoleRangeReal-OutDef',MVR,struct('IndMode',7:15)); toc
%Elapsed time is 231.298584 seconds.


dm15('ViewStabHist',stack_set(hist,'curve','dr',dr))


[h1,d1]=fe2xf('poletrackbydist upSS1 sTol1e-3 dTol.1',hist,dr)
h1=stack_set(h1,'curve','dr',dr,'curve','dclust',d1)

dm15('ViewStabHist',h1)



  
 elseif comstr(Cam,'par')
  %% #ScriptPar: sample reduced parametered multi-model study

  d_contact('discBrake Fz1e3'); cf=feplot;
  % Define contact
  data=struct('name','squeal','ProId',1001,...
   'slave','matid 1 & selface & innode{z==.013 | z==-.013} & innode{y<0.58*x & y>-0.58*x} ',... % disc surface
   'master','matid 2 & selface & innode{z==.013 | z==-.013}',...
   'contact','linear','Fu','Kc 1e12 fixKcJac0','KcLin',1e12,... % pad surfaces
   'friction','coulomb arctan2','Fv','Mu 0.1 CtLin 50',...
   'Euler',1,'Vel',4,'TangStick',0,'Sel',1,'MAP','polar'); % pad surfaces

  ctc_utils('GenerateContactPair',cf.mdl,data); % generate contact model

  ctc_utils('ParInit',cf.mdl,...
   {'lab"knmu" min-6 max1 cur1 scale"log" nPoints250 info"Kn -proid 1001"'});


  RO=struct('name','knmu','NormOpt',[-1 0 0 1e-12],'ListR','MinMaxM','setPar',1,'SEPro',1);

  MVR=d_squeal('SolveMVR',cf.mdl.GetData,RO)

  [hist,dr]=fe2xf('PoleRange-def',MVR,struct('IndMode',1:40));

  fe2xf('PlotPoleSearch-reset',hist)


  % xxx also parcoexp
  % xxx sliding ods expansion sensty to formulation
  


 elseif comstr(Cam,'hoff')
  %% #ScriptHoffmann: Hoffmann squeal minimal model simulations
  % mode coupling phenomenon on a 2DOF system
  % calibrated limit cycle configuration


  %li{2}='SimuCfg{ModalNewmark{1m,400,uva111,rt-1e-4}SQ0{vq1Amp__5}}'
  % sdtweb d_fetime simuurn
  % sdtweb d_fetime sq0
  
  RT=d_doe('nmap','TV.Hoff{n,base.MN}') % exactly on limit cycle
  sdtm.range(RT);mo2=RT.nmap('CurModel');d2=RT.nmap('CurTime');

  % xxx this illustrates zero damping cycle, not influence of NL on level 
  d1=sdth.urn('dcpx',mo2);d1.data

  %% damping variation w/ stable|unstable transition from mu on velocity
  RT=d_doe('nmap','TV.Hoff{n,Kmuv.MN}');
 sdtm.range(RT);mo3=RT.nmap('CurModel');d3=RT.nmap('CurTime');

 mo4.NL{end}

'xxxEB clean figure 1 (labels) '


 %% Hoffmann with real exponential contact: use 100x kz target at unit pressure
 RT=d_doe('nmap','TV.Hoff{n,Texp.Imp}')
 sdtm.range(RT);mo3=RT.nmap('CurModel');d3=RT.nmap('CurTime');

 %% now using modal Newmark
 RT=d_doe('nmap','TV.Hoff{n,Texp.MN}');sdtm.range(RT);mo4=RT.nmap('CurModel');d4=RT.nmap('CurTime');nlutil('postv',mo4)


 %% xxxGV : example of call using MexIOa
 RT=d_doe('nmap','TV.Hoff{n,Texp.MN}');RT.nmap('CurExp')
 MexIOa=1;
 sdtm.range(RT);mo4=RT.nmap('CurModel');d4=RT.nmap('CurTime');



 'xxx need to verify residual'

  %% Now a time varying limit cycle sdtweb d_fetime hoff TVK/TVKA
  %% #Hoff.Mssp -3
  % see also sdtweb cbi20b Script23HoffTransient
  % sdtweb d_fetime TV.Hoffman for the PerK 
  % need version for exponential law 
  RS=d_doe('nmap','TV.Hoff{n,TVKA.MN}'); sdtm.range(RS);mo4=RS.nmap('CurModel');d4=RS.nmap('CurTime');nlutil('postv',mo4);
  figure(1);semilogy(d4.data,abs(d4.def))


  if 1==2
   mo1=mo4;
   t=linspace(0,.1,1000)'; tra=struct('def',[.1;1]*t'*.06,'DOF',mo1.DOF,'data',t,'v',[.1;1]*.06*ones(size(t')));
   NL=mo1.NL{1,3};NL.opt(15)=-.2;
   [C1,C2,r2]=vhandle.chandle.ioDoStep(NL,tra);figure(1);plot(t,C1.Y(:,2)./tra.def(2,:)')
   NL=mo1.NL{2,3};
   [C1,C2,r2]=vhandle.chandle.ioDoStep(NL,tra);figure(101);plot(t,C1.Y(:,2)./tra.v(2,:)')

   mo4.K=cellfun(@(x)x.GetData,mo4.K,'uni',0);
   r1=NL.Fv{1};figure(100);plot(r1(:,1),r1(:,2)+mo4.K{2}(1)*r1(:,1));xlim([-10 10])
   xlabel('Velocity');ylabel('Viscous force')
   cingui('plotwd',100,'@OsDic(SDT Root)',{'ImSw80','WrW49c'});
  end

  NL=mo4.NL{strcmpi(mo4.NL(:,2),'TabC'),3};
  r2=vhandle.chandle.ioDoStep(NL);
  r2.u=d4.def; r2.v=d4.v; r2.fc=zeros(size(r2.u));r2.out=zeros(14,size(d4.def,2));
  R1=sdth.urn('MeshEvt',mo4);Ku=[0 1;0 0]*R1.mu*R1.kz;Ks=mo4.K{3}-Ku;
  mkl_utils(r2);
  v=(NL.c*d4.v);
  p=[v.*r2.out(NL.iopt(1)+1,:)+(sum(d4.v.*(mo4.K{2}*d4.v),1)); % viscous power
     sum(d4.v.*(Ku*d4.def),1) % Friction power
     sum(d4.v.*(Ks*d4.def+mo4.K{1}*d4.a),1) % Mechanical power
     ];
  figure(1);clf;
  ga=subplot(211);plot(d4.data,d4.def);ylabel('Displacement')
  legend('x','z');grid on;set(gca,'xticklabel',[])
  ga(2)=subplot(212);plot(d4.data,[diag([1 1 -1])*p;sum(p)]);xlabel('Time [s]');ylabel('Power output')
  grid on;;legend({'Viscous','Friction','-Mechanical','Total'})
  cingui('objset',1,{'@OsDic',{'ImSw80'},'@axes',{'xlim',[98 100]},'@line',{'linewidth',2}});setlines
  set(ga(1),'position',[.1 .62 .8 .3])
  set(ga(2),'position',[.1 .15 .8 .45])
  cingui('plotwd',1,'@OsDic(SDT Root)',{'WrW49c'});

 else; error('Script %s unknown.',CAM);
 end
 
elseif comstr(Cam,'solve'); [CAM,Cam]=comstr(CAM,6);
 %% #Solve: demo for squeal related solvers ----------------------------------
 if comstr(Cam,'modes'); [CAM,Cam]=comstr(CAM,6);
  %% #SolveModes: CEA computation
  
  % recover a model
  model=varargin{carg}; carg=carg+1;
  [ob,u1,model]=sdth.GetData(model,'-mdl'); if isempty(obj); obj=ob; end
  % recover modal resolution options
  if carg<=nargin; eigopt=varargin{carg}; carg=carg+1;
  else; eigopt=fe_def('defeigopt',model);
  end
  if ~isstruct(eigopt); eigopt=struct('EigOpt',eigopt); end
  [eigopt,st,CAM]=cingui('paramedit -DoClean',[ ...
  '-noHyst(#3#"not to add hysteretic matrix type 4")' ...
  '-silent(#3#"to get a silent call")' ...
  ],{eigopt,CAM}); Cam=lower(CAM);
  
  % Initialize solver adapted to large models
  try; ofact('mklserv_utils -silent')
  catch
   sdtw('_nb','mklserv_utils not available')
   ofact('_sel','spfmex 32 .1');
  end
  
  MVR=stack_get(model,'SE','MVR','Get');
  
  if ~isempty(MVR)||(isempty(nl_spring('getpro',model))&&~isfield(eigopt,'AssembleCall'))
   %% #SolveModesLinear: no NL in model, compute complex modes -3
   assType=0;
   if ~isempty(MVR);
    if isa(MVR.Opt,'v_handle'); MVR=sdthdf('hdfreadref',MVR); end
    mo1=rmfield(MVR,'Case'); C1=MVR.Case;
    if isfield(MVR,'NL')&&~isempty(nl_spring('getpro',MVR))
     mo1=nl_solve('TgtMdlAssemble -evalFNL',mo1);
    end
   elseif (~isfield(model,'Node')||~isfield(model,'Elt'))&&isfield(model,'DOF')
    mo1=model; C1=fe_case(mo1,'getCase'); if ~isfield(C1,'mDOF'); C1.mDOF=mo1.DOF; end
    assType=1;
   else
    [mo1,C1]=fe_case('assemble NoPNoT -se -RSeA -matdes 2 3 1 9 7',model);
    C1.mDOF=mo1.DOF;
   end
   %RO.cpx=1; RO.EigOpt=eigopt;
   RO=struct('cpx',1); RO=sdth.sfield('addmissing',RO,eigopt);
   %i1=cellfun(@(x) isa(x,'v_handle'),mo1.K);
   %if any(i1); mo1.K(i1)=cellfun(@(x) x.GetData,mo1.K,'uni',0); end
   mo1.K=cellfun(@(x)sdth.GetData(x),mo1.K,'uni',0);
   if ~assType||~isequal(mo1.DOF,C1.DOF)
    mo1.K=feutilb('tkt',C1.T,mo1.K);
   end
   m=feutilb('sumkcoef',mo1.K,ismember(mo1.Opt(2,:),[2 20]));
   k=feutilb('sumkcoef',mo1.K,ismember(mo1.Opt(2,:),[1 5 8]));
   if sdtm.Contains(Cam,'-load'); nl_solve('modeENH 9 load');
   else;   nl_solve('modeENH 9'); % Compute enhanced real mode basis
   end
   if sdtm.Contains(Cam,'-gettr'); 
    TT=T.def; T.def=[]; st=feutilb('a*b',C1.T,'TT',struct('CallVar','TT'));
    T.def=TT; 
    if isfield(C1,'mDOF'); T.DOF=C1.mDOF; else; T.DOF=mo1.DOF; end
    clear TT
    out=T; return; 
   end
   [K,r1]=feval(nl_solve('@getTKTcpx'),mo1,T.def,[],[]); % sym call
   % compute reduced complex modes
   def=fe_ceig(horzcat(K,[1:size(K{1},1)+.99]'),1000,'lr normEmc'); % constant mode norm.
   if isfield(RO,'SubDef');
    if ~isempty(RO.SubDef); ins=eval(RO.SubDef);
    else; ins=[];
    end
   else
    if eigopt.silent; stc='ModeStabClean;'; else; stc='ModeStabClean'; end
    ins=nl_solve(stc,def,RO); if isempty(ins); sdtw('_nb','ModeStabClean found no modes'); end
   end
   if ~isempty(ins); def=fe_def('subdef',def,ins); end % remove non converged modes
   if sdtm.Contains(Cam,'-keept2'); 
    TT=T.def; T.def=[]; 
    feutilb('a*b',C1.T,'TT',struct('CallVar','TT'));
    T.def=TT; T.DOF=mo1.DOF; clear TT
    def.TR=T;
   else   % restit modes
    def.alpha=def.def; def.T=T; def.def=T.def*def.def;
    def.def=C1.T*def.def; def.DOF=C1.mDOF;
    if isfield(def,'defL');
     def.defL=T.def*def.defL;
     def.defL=C1.T*def.defL;
    end
   end
   out=def;

  else
   %% #SolveModesNL -3
   % define complex modes solver options
   if ~isstruct(eigopt); eigopt=struct('EigOpt',eigopt); end
   RA=struct('traj',1,'stat',1,'matdes',[2 3 1 9 7],'fulldof',1,'keepT',2,...
    'NormOpt',[-1 0 0 1e-12],'ext','nl_solve(''modeENHcpx9'')','evalFNL',1,...
    'SubDef','nl_solve(''ModeStabClean'',def,RO);');%,...
   if eigopt.silent; RA.NormOpt(5)=0; RA.SubDef='nl_solve(''ModeStabClean;'',def,RO);'; end
   if ~eigopt.noHyst; RA=checkEta(model,RA); end % add 4 if eta exists
   %'EigOpt',eigopt);
   RB=sdth.sfield('AddMissing',eigopt,RA);
   if isequal(RB.matdes,'auto'); RB.matdes=RA.matdes; end
   % Make sure that a static state is present, or do not initialize static states
   if isempty(stack_get(model,'curve','q0'))
    sdtw('_nb','no static state')
    RB.skip=1; RB.traj=0;
   end
   % Run solver
   def=nl_solve('ModeCPX',model,RB);
   if ~isfield(def,'SE'); out=def.Mode; else; out=def; end
  end
  
  if isfield(out,'def')
   out=d_squeal('ModeData',out); % add data entries
   out=d_squeal('ModeLab',out); % add specific LabFcn
  end

  if isa(obj,'nmap');sdtm.store(obj);end % mapped output options uses NextStore key
  
 elseif comstr(Cam,'time'); [CAM,Cam]=comstr(CAM,5);
  %% #SolveTime: procedures for transient simualtions
  if comstr(Cam,'red')
  %% #SolveTimeRed: Reduction phase for transient simulations
  % stra minrio
  % rfield TR or mf
  % xxx expect an SE, or a Reduc Call
  % xxx other fields ?
  if ~isempty(obj) % cb format: call with obj and evt
   % curModel 
   [RO,u2,model]=sdth.GetData(obj,'-mdl'); % xxx nmap=expM and SolveTimeRed option
   if isKey(RO.nmap,'SolveTimeRed');nmap=RO.nmap;RO=RO.nmap('SolveTimeRed');
   else; nmap=RO.nmap; 
   end
   if ~isfield(RO,'nmap');RO.nmap=nmap; end

   if ischar(evt); def=RO.nmap(evt); else; def=nmap('CurStep'); end 
   if isfield(def,'SE'); SE=def.SE; def=def.Mode; else; SE=model; end
  else % basic format with list of arguments
   SE=varargin{carg}; carg=carg+1;
   def=varargin{carg}; carg=carg+1;
   if carg<=nargin; RO=varargin{carg}; carg=carg+1; else; RO=struct; end
  end

  [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   '-minrio(#3#"strategy keep minimal damping and use real imaginary orth basis")' ...
   '-normE(#3#"renorm with elastic matrices")' ...
   '-TimeCont(#%s#"initialize for continuation")' ...
   '-Merge(#%g#"indices of NL to be merged")' ...
   ],{RO,CAM}); Cam=lower(CAM);

  if RO.minrio
   % reduction basis generation as post
   [u1,i1]=min(def.data(:,2));

   if isfield(def,'TR') % expressed in real basis, orth and restit
    TR=orth([real(def.def(:,i1)) imag(def.def(:,i1))]);
    TR=struct('DOF',def.TR.DOF,'def',feutilb('a*b',def.TR.def,TR),'data',def.data(i1,:));

   else % xxx need renorm
    def=fe_def('subdef',def,i1);
    TR=[real(def.def) imag(def.def)];
    TR=struct('DOF',def.DOF,'def',TR,'data',def.data);

   end
  end % stra

  if RO.normE % xxx post renorm with elastic matrices only
   % I guess this is needed for diagonal modal newmark
   TR=feutilb('placeindof',SE.DOF,TR);
   K=feutilb('tkt',TR.def,SE.K);
   d1=fe_eig({K{1},K{3},[]},2);
   TR.def=TR.def*d1.def;
  end

  if isfield(SE,'K') % EB version (no using RedOpt)
   TR.adof=(1:size(TR.def,2))'+.99; RO.K=1; % Reduce matrices
   SEf=SE;
   SE=nlutil('HrRedNL',SE,TR,RO);
   % else; % we would call fe_reduc for a redefined assembly
   while nnz(SE.K{end}==0);
       i1=length(SE.K);SE.K(i1)=[];SE.Klab(i1)=[];SE.Opt(:,i1)=[];
   end
  end
  if 1==2 % Recheck cpx mode / matrices
   mo1=SEf;  NL=mo1.NL{end};NL.b=sdth.GetData(NL.b);NL.c=sdth.GetData(NL.c);
   % xxx you won't get it with fe_ceig due to matrix9 that is not symmetric
   ic=ismember(mo1.Opt(2,:),[3 7 70 11]);mo1.Klab(ic)
   ik=ismember(mo1.Opt(2,:),[1 5 8 9]);mo1.Klab(ik)
   K={feutilb('sumkcoef',mo1.K,ismember(mo1.Opt(2,:),[2 10]))
    feutilb('sumkcoef',mo1.K,ic)
    feutilb('sumkcoef',mo1.K,ik)
    }'; def=fe_ceig(horzcat(K,SE.DOF),[1 5 20 1e3]);def.data
   full(NL.b(:,1:3:end)*speye(size(NL.b,2)/3)*NL.c(1:4:end,:))
  end

   if ~isempty(RO.TimeCont)
    %% #SolveRedOpt.TimeCont prepare time continue at TimeRed exit
    if comstr(RO.TimeCont,'{')
     SE=d_fetime('TimeOpt',struct('urn',RO.TimeCont(2:end-1)),SE);
    end
    if ~isfield(SE,'nmap');SE.nmap=vhandle.nmap;end
    SE.nmap('LastContinue')='do';
    % xxx matrices handling
    if length(SE.K)>3; in1=3+find(cellfun(@nnz,SE.K(4:end))==0);
     SE.K(in1)=[];SE.Klab(in1)=[];SE.Opt(:,in1)=[];
    end
    SE.NL{end}.iopt(1)=0;SE.NL{end}.FInd=0;SE.NL{end}.iotp(4)=0;'xxx iopt'
    i3=SE.NL{end}.iopt;
    i3=sum(i3(3:4)+i3(3))*i3(5); % xxx StoreType
    SE.FNL(i3,1)=0;SE.FNLlab(end+1:i3)={''};SE.FNLDOF(end+1:i3)=0;
    % now first run
    % figure(1);NL=SE.NL{end}; plot(sort(NL.vnl(2,:,2)))
    [d2,mo3]=fe_time(SE); %iicom('curveinit','def',d2)
    
    if norm(d2.def,'Inf')>1e100; error('Diverged');end

    co2=SE.nmap('LastContinue'); sdtm.store(RO)%iigui(RO,'setInNmap');
    clear out; return;
   end


   out=SE;

  elseif comstr(Cam,'continue')
   %% #SolveTimeContinue transient continuation and display
   if carg<=nargin;RT=varargin{carg};carg=carg+1;else;RT=evalin('caller','RO');end
   co2=RT.nmap('LastContinue'); %mo2=RT.nmap('CurModel'); co2=mo2.nmap('LastContinue');
   [~,R2]=sdtm.urnPar(CAM,'{}{}');
   if ~any(co2.u)&&~any(co2.v)||any(strcmpi(R2.Failed,'randv'));
       co2.v=rand(size(co2.v));
   end

   [r3,d1]=nl_solve('fe_timeChandleContinue',co2,[],struct);
   C3=feval(nlutil('@FNL2curve'),d1,struct('drop0',1));
   C3.Y(1,:)=C3.Y(2,:);'xxx'
   C3.X{1}(:,1)=C3.X{1}(:,1)-C3.X{1}(1);
   c3=iiplot(3,';');iicom('curveinit','Time',C3)

   %i3=~any(C3.Y);C3.Y(:,i3)=[];C3.X{2}(i3,:)=[];
   %C3=fe_def('def2curve',d1);C3.X{1}(:,1)=C3.X{1}(:,1)-C3.X{1}(1);C3.X{2}=fe_c(d1.DOF);
   %    C3.Xlab{1}={'Time','s';'Fc0',''};
   if isKey(RT.nmap,'ViewSpec')
    d_squeal(RT.nmap('ViewSpec'));
    %RT.nmap('ViewSpec')='ViewSpec(BufTime .4 Overlap .8 tmin .5 -window hanning fmin 1300 1700)';
   end
   if 1==2
    d_squeal('ViewInstFreq{dmBand100,harm1,do{ReEstY,f(t),a(t)},f 1500,ifBand100,aeBand100,clipBand500 3k,,tclip .01 .01}');
    c3=iiplot(3,';');nmap=c3.data.nmap.nmap;C4=nmap('SqInstFreq');
    c14=sdth.urn('iiplot(3).clone(14)');iicom(c14,'curveinit','A1',C4);
    iicom(c3,'polarx2');
    %vhandle.cdm.urnVec(C4,'{tlim1 100}{x,2}{y,l1,abs}{gf404,tight}');
    figure(403);clf;
    h=vhandle.cdm.urnVec(C4,'{tlim1 100}{x,2}{y,l1,abs}{x,1,col}{gf403,tight,os{imgrid,ImSw80},pd{WrW49c,FnY}}');
    set(h,'linewidth',2);ii_plp('ColorMapBand',parula(6));colorbar
   end

   %% EndSolveTime
  end

 elseif comstr(Cam,'mvr'); [CAM,Cam]=comstr(Cam,4);
  %% #SolveMVR: build reduced parametric studies with SDT/Param adapted to squeal setups
  
  model=sdth.GetData(varargin{carg},'-mdl'); carg=carg+1;
  if carg<=nargin; RO=varargin{carg}; carg=carg+1; else; RO=struct; end
  RB=struct('matdes',[2 3 1 9 7],'enh',9,'RangeType','Grid');
  RO=sdth.sfield('AddMissing',RO,RB);
  
  % Linearize coupling
  if ~isempty(nl_spring('getpro',model))
   model=nl_solve('TgtMdlBuild keepName -evalFNL',model,RO);
  end
  % Assemble with parameters
  if ~ischar(RO.matdes); RO.matdes=num2str(RO.matdes(:)'); end
  [model,Case]=fe_case(sprintf('assembleNoTNoP -matdes %s -1.1 -se',RO.matdes),model);
  model=fe_def('zCoef-default',model);
  model=stack_set(model,'case','Case 1',Case);
  ofact('mklserv_utils -silent');
  bopt={};

  % Phase 1: real mode basis
  if isfield(RO,'TR');
   TR=RO.TR; RO.TR=[]; RA=[];if isfield(TR,'TR'); TR=TR.TR; end
   if ~isstruct(TR); TR=struct('def',TR,'DOF',model.DOF); end % xxx recovery
  else % generate reduction basis
   RA=RO;
   if ~isfield(RA,'list')
    RC=struct('EigOpt',fe_def('defeigopt;',model),'NoT',1,'RAMoptim',0);
    RA=sdth.sfield('AddMissing',RO,RC);
    RA=fe2xf('BuildListRMinMax',model,RA);
   end
   if isempty(RO.enh)
    RA.NoT=0; MVR=fe2xf('Build',model,RA);
   else
   [TR,kd,K]=fe2xf('Build',model,RA); TR=TR.TR; bopt={kd,K}; clear K
   end
  end
  % Phase 2: enhance with non elastic matrices (here type 9)
  RB=struct('enhtyp',RO.enh); 
  if isfield(RO,'NormOpt'); RB.NormOpt=RO.NormOpt; end
  if isfield(RO,'norm'); RB.norm=RO.norm; end
  if isfield(RO,'Load'); RB.Load=1; end
  if isfield(RA,'NoT'); if RA.NoT; RB.KismDOF=1; else; RB.KiscDOF=1; end; end % add info if NoT is controlled
  if comstr(Cam,'tr') % and output TR
   fe2xf('BuildTrEnh-norm-NoT','TR',model,Case,RB,bopt{:});
   out=TR;
  else % and reduce
   if ~isempty(RO.enh); MVR=fe2xf('BuildTrEnh-norm-MVR','TR',model,Case,RB,bopt{:}); % NoT ??
   elseif exist('MVR','var') ;% ok, carry on
   else % directly use TR
    if ~isfield(TR,'adof'); TR.adof=[1:size(TR.def,2)]'+.99; end
    MVR=sdth.GetData(model);
    K=MVR.K; MVR.K=[]; MVR.TR=TR; MVR.DOF=TR.adof;
    MVR.K=feutilb('tkt',MVR.TR.def,'K');
   end % NoT ??
   % cleanup and set range
   MVR=feutil('rmfield',MVR,'Range');
   if ~isempty(RA); MVR.RBuild=RA; end
   Range=stack_get(MVR,'info','Range0','get');
   if ~isempty(Range)
    % check if some parameters are missing from zCoef
    lab=cellfun(@(x)x.lab,Range,'uni',0);
    i5=cellfun(@iscell,lab); if any(i5); lab=cat(1,lab{:});  end
    zC=stack_get(MVR,'info','zCoef','get');
    i1=cellfun(@(x)any(~cellfun(@isempty,x)),...
     cellfun(@(x)strfind(zC(2:end,4),x),lab,'uni',0),  'uni',1);
    if any(~i1) % some parameter has no effect on any matrix
     cellfun(@(x)sdtw('_nb','Parameter %s does not have any effect, skipped',x),lab(~i1),'uni',0);
     Range(~i1)=[];
    end
    if ~strcmpi(RO.RangeType,'custom')
     Range=fe_range(sprintf('Build%s',RO.RangeType),Range,struct('replace',40));
     MVR=stack_set(MVR,'info','Range',Range);
    end
    % cleanup Omega and non varying matrices
    try
     MVR=fe2xf('BuildCleanMVR',MVR);
    catch err, err.getReport
    end
    %    if ~ismember('Omega',Range.lab)
    %     [zCoef,i1]=stack_get(MVR,'info','zCoef','get');
    %     zCoef(:,4)=strrep(zCoef(:,4),'Omega','1');
    %     MVR.Stack{i1,3}=zCoef;
    %    end
   end
   out=MVR;
  end
  
 elseif comstr(Cam,'cmt')
 %% #SolveCMT: CMT reduction for squeal models
 
 cf=varargin{carg}; carg=carg+1;
 if carg<=nargin; RO=varargin{carg}; carg=carg+1; else; RO=struct; end
 [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
  '-write(#3#"to do captures")' ...
  '-linok(#3#"model already linearized")' ...
  'Init(#3#"only to init model")' ...
  '-PostInit(#3#"to skip init if already done")' ...
  '-saveInit(#3#"to save init model")' ...
  ],{RO,CAM}); Cam=lower(CAM);
 % Ensure pad selection is done
 stlab='';
 if ~RO.PostInit % do init
  [u1,cf.mdl.Elt]=feutil('eltidfix;',cf.mdl);
  if isempty(stack_get(cf.mdl,'info','CMTOpt')); error('Missing list'); end
  % Build tangent model (nl_springs and contact)
  % call LinMdl with usual arguments
  RC=struct('CaseC2SE',1,'TgtMdl',1,'KpAuto',1,'scaleKcDiag',1,'Kp',0,'KinC2SE',1,...
   'maxKCondenseSize',0,'KeepFixDof',1);
  RB=sdth.sfield('AddMissing',sdth.sfield('sub',RO,fieldnames(RC)),RC);
  if isfield(RO,'parList'); RB.parList=RO.parList;end
  if ~RO.linok; fe_cmt('LinMdl',cf.mdl,RB); end
  
  % other parts are automatically handled
  fe_cmt('ListCheck -noRem -keepSE',cf.mdl)
  
  % Critical reset of OmegaGroups ?
  feval(abaqus('@rotDynFSet'),cf.mdl);
  
  % Prepare reduction and handle custom file names
  if isempty(cf.Stack{'curve','AssembledModes'})&&~isempty(cf.Stack{'curve','dcpx'})
   cf.Stack{'curve','AssembledModes'}=cf.Stack{'dcpx'}.TR;
   if isfield(cf.Stack{'dcpx'},'label'); stlab=cf.Stack{'dcpx'}.label; end
  elseif ~isempty(cf.Stack{'curve','AssembledModes'})
   stlab=cf.Stack{'curve','AssembledModes'}.label;
  end
  if ~isfield(cf.mdl,'wd')||isempty(cf.mdl.wd); cf.mdl.wd=sdtdef('tempdir'); end
  if isempty(strfind(cf.mdl.file,'.mat')); cf.mdl.file=[cf.mdl.file '.mat']; end
  if isfield(RO,'KnKt')
   st1=comstr(sprintf('%s',RO.KnKt{:}),-36);
   cf.mdl.file=strrep(cf.mdl.file,'.mat',sprintf('_%s_.mat',st1));
  end
  cf.mdl.file=strrep(cf.mdl.file,'.mat',sprintf('%sCMT.mat',stlab));
  f1=fullfile(cf.mdl.wd,cf.mdl.file); f2=strrep(f1,'CMT.mat','db.mat')
  cf.mdl.fileSE=f2; % HDF
  if exist(f1,'file'); delete(f1); end
  if exist(f2,'file'); delete(f2); end
  
  ofact('mklserv_utils -silent')
  %4- Build CMT model and save
  % In all cases, recompute modes with penalized version
  fe_cmt('init -matdes 2 3 1 9 7 -FixDofAsSeC -setCon',cf)
  % xxx cases where interf has internal DOF SeRestit after init is incomplete
  % now should be ok
  
  if RO.saveInit
   model=cf.mdl.GetData; model=feutil('rmfield',model,'fileSE','file');
   sOpt=sdtdef('V6Flag');
   save(fullfile(cf.mdl.wd,strrep(cf.mdl.file,'.mat','_CMTInit.mat')),'model',sOpt{:});
  end
 else; f1=fullfile(cf.mdl.wd,cf.mdl.file);
 end
 if RO.Init; out=cf; return; end
 
 RA=cf.Stack{'info','EigOpt'};
 
 RA=struct('EigOpt',RA);
 RA.AssembleCall=nl_spring('AssembleCall -compose"NoP -matdes 2 3 1 9 7 -RSeA -load"');
 RA.matdes=[2 3 1 9 7];
 fe_case(cf.mdl,'reset');
 def=d_squeal('SolveModes',cf.mdl,RA); if ~isempty(stlab); def.TR.label=stlab; end
 cf.Stack{'curve','AssembledModes'}=def.TR;
 cf.Stack{'CMTOpt'}.Global='AssembledModes';
 % then do enh or not
 % then Build and save
 if isfield(RO,'KnKt') % advanced build
  %fe_cmt('Init -matdes 2 3 1 9 7',cf)
  %fe_cmt('ListConInit',cf.mdl)
  fe_cmt('ListConParSet',cf.mdl,RO.KnKt)
  fe_cmt('ListConParNTSplit',cf.mdl)
  fe_cmt('defred hdf -gENH 9',cf)
  fe_cmt('Build normENH hdf cutoff12e3 -sPost -matdes 2 3 1 9 7 -SepKi -forceStra1',cf)
 else
  fe_cmt('DefRed hdf cpfmax12e3',cf);
  fe_cmt('Build normENH hdf cutoff12e3 -sPost -matdes 2 3 1 9 7 -SepKi -forceStra1',cf);
  %fe_cmt('Init-Build normENH hdf cutoff7e3 -sPost -matdes 2 3 1 9 7 -SepKi -forceStra1',cf);
 end
 
 fe_cmt('save',cf,f1)
 
 out=cf;
 
 elseif comstr(Cam,'ods'); [CAM,Cam]=comstr(CAM,4);
  %% #SolveODS: ODS extraction functionalities -2
  if comstr(Cam,'tse')
   %% #SolveODSTSE: time series events detection -3

   C3=varargin{carg}; carg=carg+1; % sensors time respone
   if carg<=nargin; RC=varargin{carg}; carg=carg+1; else; RC=struct; end
   RO=struct('lab','','imw',0,'aThr',5,'aM',0.8,'iref',1,'delay',1000,'lpf',100,'tdrop',0,'addCol',{{}});
   RO=sdth.sfield('AddMissing',RC,RO);
   if ~isfield(RO,'aThrD'); RO.aThrD=RO.aThr; end
   
   Fs=1./unique(diff(C3.X{1})); Fs=mean(Fs); % 8kHz
   
   C4=C3; C5=C3; C4.name='Amp'; C5.name='IFreq';
   for j1=1:size(C4.Y,2)
    [A,F,phi]=feval(ii_hvd('@inst'),C3.Y(:,j1),Fs);
    C4.Y(:,j1)=A(:); C5.Y(:,j1)=F(:);
   end
   
   C4l=ii_mmif(sprintf('bandpass fmin0 fmax%g -struct',RO.lpf),C4);
   C5l=ii_mmif(sprintf('bandpass fmin0 fmax%g -struct',RO.lpf),C5);
   
   % events:  
   y4=max(C4l.Y,[],2);
   in1=find(y4>RO.aThr); % 40% of time
   in2=[1; find(diff(in1)>RO.delay)+1; length(in1)];
   
   t1=zeros(length(in2)-1,6);
   
   %j1=4; %ir=5
   for j1=1:size(t1,1)
    i2=find(y4(in1(in2(j1))+1:in1(in2(j1+1)))<=RO.aThrD,1,'first');
    if isempty(i2); i2=in2(j1+1)-in2(j1); end
    % in3=in1(max(1,in2(j1))-500):(in1(in2(j1))+i2); % 5000
    
    in33=in1(max(1,in2(j1))):(in1(in2(j1))+i2); % 5000
    fm=mean(C5l.Y(in33,RO.iref)); fsd=std(C5l.Y(in33,RO.iref));
    % we have the window: fmean,Amax, Amean  Fm(Amean)
    aa=abs(C4.Y(in33,:));
    aM=max(max(aa));
    a10=max(mean(aa(aa>aM*RO.aM),1));
    fm10=mean(C5l.Y(in33(max(aa,[],2)>aM*RO.aM),RO.iref));
    l1=length(in33); if l1/Fs<RO.tdrop; l1=-1; end
    
    t1(j1,1:6)=[fm fsd aM l1/Fs a10 fm10]; if l1==-1; continue; end
    
    ic=1;
    for j2=1:size(RO.addCol,1)
     C2=RO.addCol{j2,1}; %{C1,1,'p',{'mean','std','maM'}
     stl=RO.addCol{j2,3};
     if isequal(C2.X{1},C3.X{1}); in5=in33;
     else; % another time vector in context
      in5=C2.X{1}>=min(C4l.X{1}(in33))&C2.X{1}<=max(C4l.X{1}(in33));
      if isempty(in5)||~any(in5);in5=find(C2.X{1}>=min(C4l.X{1}(in33)),1,'first'); end
     end
     for j3=1:length(RO.addCol{j2,4})
      switch lower(RO.addCol{j2,4}{j3})
       case 'mean';  r5=mean(C2.Y(in5,RO.addCol{j2,2}));
       case 'std';   r5=std(C2.Y(in5,RO.addCol{j2,2}));
       case 'mam';   r5=mean(C2.Y(in5(max(aa,[],2)>aM*RO.aM),RO.addCol{j2,2})); % xxx only valid for in5=in33
       otherwise; sdtw('_nb','%s is not a known feature',RO.addCol{j2,4}{j3});
      end
      t1(j1,6+ic)=r5; ic=ic+1;
     end
    end
    
    if RO.imw
     figure(111); plot(C5l.X{1}(in33),C5l.Y(in33,:)); setlines
     legend(C3.X{2}(:)','location','eastoutside','interpreter','none')
     axis tight;xlabel('Time [s]'); ylabel(sprintf('LF%s inst. Freq',RO.lpf));
     title(sprintf('%s - Event %i',RO.lab,j1));
     cingui('PlotWd',111,'@OsDic(SDT Root)',{'ImSw80','Fnrtxy','WrW32c'});
     comgui('ImWrite',111)
     
     figure(112); plot(C4l.X{1}(in33),C4l.Y(in33,:)); setlines
     legend(C3.X{2}(:)','location','eastoutside','interpreter','none')
     axis tight;xlabel('Time [s]'); ylabel(sprintf('LF%s inst. Amp',RO.lpf))
     title(sprintf('%s - Event %i',RO.lab,j1));
     cingui('PlotWd',112,'@OsDic(SDT Root)',{'ImSw80','Fnrtxy','WrW32c'});
     comgui('ImWrite',112)
    end
    
   end
   
   in2(t1(:,4)<0,:)=[]; t1(t1(:,4)<0,:)=[];
   
   out={in1,in2}; % xxx
   
   % lets' do a table
   cmap=feval(cinguj('@rgb2gray'),[1 1 1; .18 .75 .1;1 .75 0;1 0 0;],'bt601','ft');
   crat=struct('clevel',[-1 0 10 40 2000 ],'cmap',cmap,...
    'Fcn',fegui('@CritFcn'));
   %      crph=struct('clevel',[-1 0 1 2 100 ],'cmap',cmap,...
   %       'Fcn',fegui('@CritFcn'));
   crt=struct('clevel',[-1 0 1 5 5000 ],'cmap',cmap,...
    'Fcn',fegui('@CritFcn'));
   
   lab={'Fmean','Amax','Av(Am20)','T(Am20)','Fmean(Am20)','Fstd'};
   if ~isempty(RO.addCol)
    lb1=RO.addCol(:,3);
    lb2=cellfun(@(x)x(:),RO.addCol(:,4),'uni',0);
    lb2=cellfun(@(x,y)cellfun(@(x1)[x x1],y,'uni',0),lb1,lb2,'uni',0);
    lb2=cat(1,lb2{:});
    lab=[lab lb2(:)'];
   end
   
   
   ua=struct('table',{num2cell([(1:size(t1,1))' t1(:,[1 3 5 4 6 2 7:end])])},...
    'name',RO.lab,'FigTag','SDT Root',...
    'ColumnName',{[['Event',lab];
    repmat({''},1,length(lab)+1);
    ['00',repmat({'0.00'},1,length(lab))];
    [{[] [] crat crat crt} repmat({[]},1,length(lab)-4)]]},...
    'RowHeight',18,'ColWidth',[ repmat(60,1,length(lab)+1)],'setSort',1);
   [un1,un2,u0]=comstr(ua.table,-17,'tab',ua);
   
   out1=ua;
   out2={C4l,C5l};

  elseif comstr(Cam,'ext')
   %% #SolveODSExt: extract ODS from time series events -3
   C3=varargin{carg}; carg=carg+1;
   RC=varargin{carg}; carg=carg+1;
   ua=varargin{carg}; carg=carg+1;
   ind=varargin{carg}; carg=carg+1; in1=ind{1}; in2=ind{2};
   XF0=varargin{carg}; carg=carg+1; C4l=XF0{1}; C5l=XF0{2};
   if ~isfield(RC,'fbw'); RC.fbw=100; end
   if ~isfield(RC,'aThr'); RC.aThr=5; end
   
   y4=max(C4l.Y,[],2);
   % first indicator: SVD on each window +MAC or SS
   %j1=4; ir=5;
   s1=zeros(size(ua.table,1),2); u1=struct('DOF',{C3.X{2}},'def',[],'data',[]);
   for j1=1:size(ua.table,1,1)
    i2=find(y4(in1(in2(j1))+1:in1(in2(j1+1)))<=RC.aThr,1,'first');
    if isempty(i2); i2=in2(j1+1)-in2(j1); end
    % in3=in1(max(1,in2(j1))-500):(in1(in2(j1))+i2); % 5000
    
    in33=in1(max(1,in2(j1))):(in1(in2(j1))+i2); % 5000
    C33=C3; C33.X{1}=C33.X{1}(in33); C33.Y=C33.Y(in33,:);
    %figure; plot(C33.X{1},C33.Y)
    %1/(C33.X{1}(end)-C33.X{1}(1)) % df 6Hz
    
    XF=ii_mmif('FFT -struct',C33);
    % recover peak around fmean
    fm=mean(C5l.Y(in33,RC.iref));
    
    fr=XF.X{1}>=fm-RC.fbw&XF.X{1}<=fm+RC.fbw;
    %fr=XF.X{1}>=fm-50&XF.X{1}<=fm+50;
    if nnz(fr)
     [uu,ss,vv]=svd(XF.Y(find(fr),:).',0);
     if length(ss)>1
      s1(j1,:)=[ss(1) ss(2,2)/ss(1)];
     else; s1(j1)=ss(1);
     end
     u1.def=[u1.def uu(:,1)];
    elseif isempty(u1.def); u1.def=zeros(size(XF.Y,2),1);
    else; u1.def(:,end+1)=0;
    end
    u1.data=[u1.data;fm];
   end
   
   out=u1; out1=s1;
     
  elseif comstr(Cam,'merge')
   %% #SolveODSMerge: merge extracted ODS batches exploiting references -3
   dbstack, keyboard
 
  else; error('SolveODS %s unknown.',CAM)
  end
 elseif comstr(Cam,'chandle')
  %% #SolveCHandle go to form with merged contact 
  RT=varargin{carg};carg=carg+1;
  mo1=RT.nmap('mod1');%(15)=.45;
  for j1=1:size(mo1.NL,1)
   NL=mo1.NL{j1,3};
   if ~isfield(NL,'CHandle');NL=feval(nl_contact('@legToChandle'),NL,mo1);end
   mo1.NL{1,3}=NL;
  end
  RT.nmap('mod1')=mo1;


 else; error('Solve %s unknown.',CAM)
 end
 %% #Mode: GUI cbks and posts
elseif comstr(Cam,'mode'); [CAM,Cam]=comstr(CAM,5);
 if comstr(Cam,'cpx'); [CAM,Cam]=comstr(CAM,4);
  %% #ModeCPX: call from GUI -2
 if comstr(Cam,'type'); out='nltgt';
 else; out=d_squeal('SolveModes',varargin{carg:end});
 end
  
 elseif comstr(Cam,'data')
  %% #ModeData: add indicators in data field
  def=varargin{carg}; carg=carg+1;
  if isfield(def,'fun')&&def.fun(2)==3 % cpx modes
   if isa(def.data,'v_handle'); def.data=def.data(:,:); end
   def.data(:,3)=1:size(def.data,1);
   if isfield(def,'Xlab')
    def.Xlab{end}{3}='Index';
    if ~ismember('SqC',def.Xlab{end})
     def.Xlab{end}{end+1}='SqC';  % Squeal Coefficient
     lj=feval(fe_ceig('@getLambda'),def.data(:,1:2));
     def.data(:,end+1)=100*2*real(lj)./abs(lj);
    end
   end
  end
  out=def;
  
   elseif comstr(Cam,'lab')
  %% #ModeLab: apply squeal labels to modes -2
  def=varargin{carg}; carg=carg+1;
  if isfield(def,'Xlab')&&iscell(def.Xlab{end})
   r1={...
    'Index','Mode %i',1;
    'Freq','f=%.4g Hz',1;
    'Damp',{'%s=%.2g%%','''\zeta'''},100;
    'SqC','SqC=%.2g%%',1;
    'VNS','V_s=%.3gmm^2/s^2',1;
    'NLA1','A_{NL1}=%.3g',1};
   r2=def.Xlab{end}; st=''; sta='';
   for j1=1:size(r1,1)
    i1=find(ismember(r2,r1{j1,1}));
    if ~isempty(i1)
     r3=r1{j1,2};
     if iscell(r3)
      r4=r3(2:end); r4=[r4(:)']; r4(2,:)={','}; r4=reshape(flipud(r4),[],1);
      sta=[sta r4{:}]; r3=r3{1};
     end
     st=[st r3 '\n']; sta=[sta ',' sprintf('%g*def.data(ch,%i)',r1{j1,3},i1)];
    end
   end
   def.LabFcn=['sprintf(''' st '''' sta ''');'];
  end
  out=def;
  
 else; error('Mode %s unknown.',CAM);
 end
 
 
elseif comstr(Cam,'parezz'); [CAM,Cam]=comstr(CAM,7);
 %% #ParEzz: lining transverse young modulus parametrization -1
 if comstr(Cam,'fromparam')
  %% #ParEzzFromParam: setup appropriate ParFcn in param fields if existing -2
  %   'ParFcn', {{@fe_caseg,'parMat',sprintf('E -matid %i',RunOpt.MatId)}}));
  % check if a parameter was defined in material
  model=varargin{carg}; carg=carg+1; par={};
  mat=fe_mat('getpl',model); st=fe_mat('typemstring',mat(:,2));
  mat=mat(ismember(st,'m_elastic.6'),1);
  if ~isempty(mat)
   for j1=mat(:)'
    r1=matgui('getstackpl',model,j1);
    if ~isempty(r1)&&isfield(r1{1,3},'param')
     st2=fieldnames(r1{1,3}.param);
     r1{1,3}.param.ParFcn={@d_squeal,sprintf('ParEzz -lab"%s" -sel"matid %i"',st2{1},j1),...
      r1{1,3}.param.(st2{1})};
     model=stack_set(model,r1);
    end
   end
  end
  out=model;
  
 else
  %% #ParEzzDefinition of parameter Ezz -2
  [CAM,Cam,RO.keepOther]=comstr('-keepOther',[-25 3],CAM,Cam);
  [CAM,Cam,RO.Eta]=comstr('-eta',[-25 2],CAM,Cam);
  [CAM,Cam,RO.GlobalE]=comstr('-globale',[-25 3],CAM,Cam);
  [CAM,Cam,RO.Hvar]=comstr('-hvar',[-25 3],CAM,Cam);
  [CAM,Cam,RO.sel]=comstr('-sel',[-25 4],CAM,Cam);
  [CAM,Cam,RO.lab]=comstr('-lab',[-25 4],CAM,Cam); if isempty(RO.lab); RO.lab='Ezz'; end
  RO.hasVisc=sdtcheck('fun','fevisco');
  model=sdth.GetData(varargin{carg}); carg=carg+1; par={};
  if carg<=nargin; par=varargin{carg}; carg=carg+1; end
  if ~iscell(par); par={par}; end
  if isempty(RO.sel)
   if carg<=nargin; RO.sel=varargin{carg}; carg=carg+1;
   else % get all elts
    mat=fe_mat('getpl',model); st=fe_mat('typemstring',mat(:,2));
    mat=mat(ismember(st,'m_elastic.6'),1);
    if ~isempty(mat);RO.sel=sprintf('matid %s',num2str(mat(:)'));
    else; sdtw('_nb no orthotropic material found'); out=model; return;
    end
   end
  end
  
  if (isempty(RO.GlobalE)||~RO.GlobalE)&&(isempty(RO.Hvar)||~RO.Hvar)&&~RO.hasVisc
   RO.GlobalE=1; sdtw('_nb','Can only use global variation of Young modulus')
  end
  
  % set a single ezz to all lining
  [i1,elt]=feutil(sprintf('findelt %s',RO.sel),model);
  if isempty(elt); sdtw('_nb no orthotropic material found'); out=model; return; end
  mpid=feutil('mpid',model.Elt); m1=unique(mpid(i1,1)); mpid(i1,1)=m1(1);
  model.Elt=feutil('mpid',model.Elt,mpid);
  if ~isempty(RO.Eta)&&RO.Eta % assign Eta to m1(1)
   model=fe_mat(sprintf('setmat %i Eta=%.15g',m1(1),RO.Eta),model);
  end
  
  if ~isempty(RO.GlobalE)&&RO.GlobalE
   % declare lining modulus variation
   mo1=fe_caseg('parMat',model,...
    sprintf('E1 E2 E3 nu23 nu31 nu12 G23 G31 G12 -matid %i',m1(1)), par{1});
   
  elseif ~isempty(RO.Hvar)&&RO.Hvar   % vary eta
   mo1=fe_caseg('parMat',model,sprintf('Eta -matid %i',m1(1)),par{1});
   
  else
   % matsplit for Ezz
   RA=struct('type','ortho','Group',{{RO.lab,3; % name is standardized for CMT
    'Other',[1 2 4 5 6 7 8 9]}},'ParTyp',1,'MatId',m1(1));
   mo1=fevisco('MatSplit-appendPar',model,RA);
   if RO.Eta % apply iso Eta on components
    in1=m1(1)*100+(1:2);
    for j1=in1(:)'; mo1=fe_mat(sprintf('setmat %i Eta=%.15g',j1,RO.Eta),mo1); end
    mo1=fe_caseg('parmat',mo1,sprintf('eta -matid%i',in1(1)),struct('lab',{{[RO.lab 'H']}},'zCoef',RO.lab));
   end
   if RO.keepOther
    if RO.Eta
     mo1=fe_caseg('parmat',mo1,sprintf('eta -matid%i',in1(2)),struct('lab',{{'OtherH'}},'zCoef','Other'));
    end
   else; mo1=fe_case(mo1,'remove','Other');
   end
  end
  % define Ezz parameter range0
  if ~isempty(par);
   i1=cellfun(@ischar,par);
   par(i1)=cellfun(@(x)fe_range('buildvect',x),par(i1),'uni',0);
   Ra0=stack_get(mo1,'info','Range0','get'); % Needed for fe_range build calls
   if isempty(Ra0); Ra0=par(:)'; 
   else; 
    lb1=cellfun(@(x)x.lab,par,'uni',0); i5=cellfun(@iscell,lb1); lb1(i5)=cellfun(@(x)x{:},lb1(i5),'uni',0);
    lb0=cellfun(@(x)x.lab,Ra0,'uni',0); i5=cellfun(@iscell,lb0); lb0(i5)=cellfun(@(x)x{:},lb0(i5),'uni',0);
    %[u1,i1]=setdiff(cellfun(@(x)x.lab{:},par,'uni',0),cellfun(@(x)x.lab{:},Ra0,'uni',0));
    [u1,i1]=setdiff(lb1,lb0);
    par=par(i1);
    Ra0=[Ra0; par(:)]; 
   end
   mo1=stack_set(mo1,'info','Range0',Ra0);
  end
  out=mo1;
 end
 
elseif comstr(Cam,'view'); [CAM,Cam]=comstr(CAM,5);
 %% #View: post-treatment displays associated to squeal
 if comstr(Cam,'stab')
  %% #ViewStab: display squeal stability -2
  def=varargin{carg}; carg=carg+1;
  if carg<=nargin; RO=varargin{carg}; carg=carg+1; else; RO=struct; end
  [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   '-write(#3#"to trigger capture")' ...
   '-OsDic("SDT Root"#%s#"figtag")' ...
   '-gf(#%i#"to provide a figure number")' ...
   '-ftarg("all"#%s#"provide a frequency target")' ...
   '-unst(#3#"to zoom in on unstable area")' ...
   '-reset(#3#to clear current figure")' ...
   '-table(#3#"to display pole table")'...
   '-lcmap("lines"#%s#"line colormap")' ...
   '-mk("dxos^<>v"#%s#"markers to be used")' ...
   '-linestyle(":"#%s#"linestyle")' ...
   '-mfcol(1#%i#"apply marker face color")' ...
   '-zmax(#%g#"damping saturation")' ...
   '-Hz(#3#"to apply hertz")' ...
   ],{RO,CAM}); Cam=lower(CAM);
  if isfield(def,'Mode'); def=def.Mode; end
  if isstruct(def); data=def.data; else; data=def; end
  RO.ylab='Damping [%]'; RO.yulim='@tight .001 vis max .1';
  if RO.Hz; RO.xcoef=1; RO.Hz='Hz'; else; RO.xcoef=1e-3; RO.Hz='kHz'; end
  r1=1e2*data(:,2);
  
  if RO.table
   %% #ViewStabTable -3
   table=struct('table',{num2cell([data(:,1) r1])},'name','stability',...
    'NeedClose',1,'ColumnName',{[{'Frequency [Hz]',RO.ylab};{'',''};{'0.00','0.000'}]},...
    'RowHeight',20,'FigTag','SDT Root');
   comstr(table.table,-17,'tab',table)
   
  else
   %% #ViewStabDiagram -3
   if isempty(RO.gf); RO.gf=double(figure);
   else; RO.gf=double(figure(RO.gf));
   end
   if isempty(findall(RO.gf,'type','line'))||RO.reset
    plot(RO.xcoef*data(:,1),r1,'linestyle',RO.linestyle,'marker','o','markersize',4);
    if isfield(def,'label'); legend(def.label,'location','southeast'); end
   else; % append plots
    delete(findall(RO.gf,'tag','now'));
    line(RO.xcoef*data(:,1),r1,'linestyle',RO.linestyle,'marker','o','markersize',4,'tag','poled');
    i1=length(findobj(RO.gf,'tag','poled'))+1;
    cingui('objset',RO.gf, ...       % Handle to the object to modify
     {'',sprintf('@out=setlines(%s(%i),'':'',''%s'',%i);',RO.lcmap,max(2,i1),RO.mk,RO.mfcol)}); % Line seq
    ob=findobj(RO.gf,'tag','legend');
    if ~isempty(ob)&&isfield(def,'label')
     st=get(ob,'String'); st=cat(2,st(1:end-1),def.label);
     legend(st,'location','southeast')
    else; legend('off')
    end
   end
   ii_plp(struct('po',[0 0],'marker','horizontal',...
    'LineProp',{{'color','r','linewidth',.5,'linestyle',':'}}),...
    get(RO.gf,'CurrentAxes'))
   xlabel(sprintf('Frequency [%s]',RO.Hz)); ylabel(RO.ylab);axis('tight');
   if isfield(def,'label');RO.tit=def.label; else; RO.tit=''; end
   % close-up
   if ~isempty(RO.ftarg)
    if comstr(RO.ftarg,'all') % do nothing
    elseif comstr(RO.ftarg,'[')
     RO.tit=[RO.tit sprintf(' %s %s',regexprep(RO.ftarg,'[\[\]]',''),RO.Hz)];
     RO.ftarg=eval([RO.ftarg ';']);
     cingui('objset',RO.gf,{'@axes',{'xlim',RO.ftarg}});
     cingui('objset',RO.gf,{'@axes',{'ylim','@tight -vis coef.001'}});
     
    else; error('ftarg %s unknown',RO.ftarg)
    end
   end
   if RO.unst;
    RO.tit=['Unstable ' RO.tit];
    cingui('objset',RO.gf,{'@axes',{'ylim',RO.yulim}});
    if ~isempty(get(get(RO.gf,'currentaxes'),'Legend')); legend('location','southeast'); end
   else; cingui('objset',RO.gf,{'@axes',{'ylim','@tight .001 vis'}});
    if ~isempty(get(get(RO.gf,'currentaxes'),'Legend')); legend('location','northeast'); end
   end
   if ~isempty(RO.zmax)
    cingui('objset',RO.gf,{'@axes',{'ylim',[NaN RO.zmax]}});
   end
   title(RO.tit,'interpreter','none');
   % PlotWd: printing options
   if ~isempty(RO.OsDic); stdic=sprintf('@OsDic(%s)',RO.OsDic); end
   cingui('PlotWd',RO.gf,stdic,{'ImToFigN','ImSw','WrW80c','Fnrtxy'})
   
   if RO.write
    comgui('ImWrite',RO.gf)
   end
   
  end
 elseif comstr(Cam,'shapes')
  %% #ViewShapes: display squeal shapes -2
  cf=varargin{carg}; carg=carg+1;
  if carg<=nargin; def=varargin{carg}; carg=carg+1; else; def=[]; end
  if carg<=nargin; RO=varargin{carg}; carg=carg+1; else; RO=struct; end
  [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   '-write(#3#"to trigger capture")' ...
   '-unst(#3#"to display unstable modes only")' ...
   '-fmin(-Inf#%g#"minimum frequency")'...
   '-fmax(Inf#%g#"maximum frequency")' ...
   '-movie(#3#"to animate modes and write gif instead of png")' ...
   '-lg("LgBRflip"#%s#"provide specific legend entry")'...
   '-aSet(#%s#"set to display as alpha")' ...
   '-doSel(#3#"to only do selection")' ...
   ],{RO,CAM}); Cam=lower(CAM);
  
  % reset display selection
  
  try; fecom(cf,'colorbaroff'); end; cf.def=[]; 
  clf(cf.opt(1));
  if isfield(RO,'cbSel'); feval(RO.cbSel{1},cf,RO,'ViewShapeSel',RO.cbSel{2:end})
  elseif ~isempty(cf.Stack{'fsel','CMTsel'}) % handle CMT selection
   fsel=cf.Stack{'fsel','CMTSel'};
   li=cf.Stack{'SeRestit'}.SEname;
   for j1=1:size(li,1); fsel.Edit.(li{j1}).value='groupall'; end
   fsel.Edit.Global.value='groupall';
   if isempty(RO.aSet)
    cf.Stack{'fsel','CMTSel'}=fsel;
    fe_cmt('Restit',cf); drawnow;
    fecom(cf,'view3')
    fecom(cf,'SetProp sel(1).fsProp','FaceAlpha',1,'EdgeAlpha',0.1);
   else
    fs1=fsel; fs2=fsel;
    li2={RO.aSet}; li1=setdiff(li(:,1),li2);
    if ~isequal(li1,li(:,1)) % we found a component to alphs
    for j1=1:length(li2); fs1.Edit.(li2{j1}).value=''; end
    for j1=1:length(li1); fs2.Edit.(li1{j1}).value=''; end
    else % no component run on global
     fs1.Edit.Global.value=sprintf('setname %s:exclude',RO.aSet);
     fs2.Edit.Global.value=sprintf('setname %s',RO.aSet);
    end
    cf.Stack{'fsel','SelAlpha'}=fs2;
    cf.Stack{'fsel','SelNoalpha'}=fs1;
    sel1=fe_cmt('RestitSelNoShow-NoPro-fsel"SelNoAlpha"',cf);
    sel2=fe_cmt('RestitSelNoShow-NoPro-fsel"SelAlpha"',cf);
    cf.SelF{1}=sel1;
    cf.SelF{2}=sel2;
    try; uaob=cf.ua.ob(1,:); catch; uaob=[0 1 0 0 1 0 1]; end
    cf.o(1)=sprintf('sel 1 def %i ch %i ty%i scc %g',uaob([5 4 2 7])); % mesh
    cf.o(2)=sprintf('sel 2 def %i ch %i ty1 scc %g',uaob([5 4  7])); % mesh
    fecom(cf,'scaleequal')
    fecom(cf,'setpropsel(1).fsProp','FaceAlpha',1,'edgealpha',.1)
    fecom(cf,'setpropsel(2).fsProp','FaceAlpha',.1,'edgealpha',.1)
    fecom(cf,'colordatasel1evala')
    fecom(cf,'colordatasel2evala')
   end
  elseif ~isempty(RO.aSet)
   cf.sel(1)=sprintf('setname %s:exclude -linface',RO.aSet);
   cf.sel(2)=sprintf('setname %s -linface',RO.aSet);
   %cf.o(1)='sel 1 def 1 ch 1 ty1 scc 32'; % mesh
   cf.o(2)=sprintf('sel 2 def %i ch %i ty%i scc %i',cf.ua.ob(1,[5 4 2 7]));
   drawnow
   iimouse('zoomreset',cf.ga)
   fecom(cf,'SetProp sel(1).fsProp','FaceAlpha',1,'EdgeAlpha',0.1);
   fecom(cf,'SetProp sel(2).fsProp','FaceAlpha',.1,'EdgeAlpha',0.1);
   fecom(cf,';colorsel1dataevala;colorsel2dataevala;scaleequal;');
  else
   cf.sel='eltname ~=SE -linface'; drawnow
   fecom(cf,'view3')
   fecom(cf,'SetProp sel(1).fsProp','FaceAlpha',1,'EdgeAlpha',0.1);
  end
  if ~isempty(cf.Stack{'cv3'})
   iimouse('view',cf.ga,cf.Stack{'cv3'}.cv)
  else;fecom(cf,'view3')
  end
  if RO.doSel; return; end
  % recover der
  if isempty(def); def=cf.def; end
  if isempty(def); error('no def to display'); end
  if isfield(def,'Mode'); def=def.Mode; end
  % unstable modes indices
  RO.iUnst=def.data(:,2)<-1e-5;
  
  data=def.data;
  RO.iDef=def.data(:,1)<RO.fmax&def.data(:,1)>RO.fmin;
  % recover channels with unstable modes
  if RO.unst;RO.iDef=find(RO.iUnst&RO.iDef);RO.ch='all';
  else; RO.iDef=find(RO.iDef); RO.ch=find(RO.iUnst(RO.iDef));
  end
  % assign def with clean labels
  if isfield(def,'label'); def.label=strrep(def.label,'_',' '); end
  if isempty(RO.iDef); sdtw('_nb','No mode to display'); return; end
  cf.def=fe_def('subdef',def,RO.iDef);
  cf.def.Legend=sdtroot('CbOsDic',[],{RO.lg});
  fecom(cf,'colordataevala');
  % Print options
  %cingui('objset',cf,{'@OsDic',{'CmJet'}});
  r1={}; if RO.movie; r1=[r1 'ImMovie']; end
  if ~isempty(r1);  cingui('PlotWd',cf,'@OsDic(SDT Root)',r1); end
 
  if RO.write
   fecom(cf,'ImWrite',struct('ch',RO.ch));
  end
  
 elseif comstr(Cam,'ods'); [CAM,Cam]=comstr(CAM,4);
  %% #ViewODS: post treatement of ods extraction on events -2
  if comstr(Cam,'events')
   %% #ViewODSEvents: some plots to illustrate events -3
   ua=varargin{carg}; carg=carg+1; % table
   XF=varargin{carg}; carg=carg+1; % Alp, Flp
   ind=varargin{carg}; carg=carg+1;
   RO=varargin{carg}; carg=carg+1;
   if ~isfield(RO,'imw'); RO.imw=0; end
   
   C4=XF{1}; C5=XF{2}; in1=ind{1}; in2=ind{2};
   
   UI=sdtroot('PARAMUI');
   cingui('PlotWd',UI.gf,'@OsDic(SDT Root)',{'WrW32c'})
   %comgui('ImWrite-JavaT',UI.gf,'BR206_800Hz_b1.png')
   [un1,un2,u0]=comstr(ua.table,-17,'tab',ua);
   if RO.imw; comgui('ImWrite-JavaT',UI.gf,comstr([RO.lab '.png'],-36)); end
   
   y4=max(C4.Y,[],2);
   gf=figure; plot(C4.X{1},y4);
   xlabel('Time [s]'); ylabel('Max Amp');
   axis tight
   ii_plp(struct('po',C4.X{1}(in1(in2))))
   title(sprintf('%s - Events instants',RO.lab))
   cingui('PlotWd',gf,'@OsDic(SDT Root)',{'ImSw60','Fnrtxy','WrW32c'});
   if RO.imw; comgui('ImWrite',gf); end
   
   
   gf=figure; plot(C5.X{1},C5.Y(:,RO.iref));
   xlabel('Time [s]'); ylabel(sprintf('LF%s inst. Freq [Hz]',RO.lpf));
   axis tight
   ii_plp(struct('po',C4.X{1}(in1(in2))))
   title(sprintf('%s - Events instants',RO.lab))
   cingui('PlotWd',gf,'@OsDic(SDT Root)',{'ImSw60','Fnrtxy','WrW32c'});
   if RO.imw; comgui('ImWrite',gf); end
   
   
  elseif comstr(Cam,'macauto')
   %% #ViewODSMACAuto: check ods consistency by events -3
  % dbstack,keyboard
   u1=varargin{carg}; carg=carg+1;
   s1=varargin{carg}; carg=carg+1;
   RO=varargin{carg}; carg=carg+1;
   if~isfield(RO,'imw'); RO.imw=0; end
        
     gf=figure; bar(s1(:,2)); axis tight
     ii_plp(struct('po',[0.1 0;0.05 0],'marker','horizontal'))
     xlabel('Event'); ylabel('SVD Rel2')
     title(sprintf('%s ODS Rel2',RO.lab))
     cingui('PlotWd',gf,'@OsDic(SDT Root)',{'ImSw60','Fnrtxy','WrW32c'});
     if RO.imw; comgui('ImWrite',gf); end
     
     
     
     gf=figure; ii_mac(gf,u1.def,'inda',find(s1(:,1)),'macautoa')
     title(sprintf('%s - MAC Auto ODS by Event',RO.lab))
     cingui('PlotWd',gf,'@OsDic(SDT Root)',{'ImLw50','Fnrtxy','WrW32c'});
     if RO.imw; comgui('ImWRite',gf); end

   
  else; error('ViewODS %s unknown.',CAM);
  end
  
elseif comstr(Cam,'time')
%% #ViewTime{par} : standardized viewing of time response 
  [st,RO]=sdtm.urnPar(CAM,'{}:{ci%g,tmin%g,tmax%g,ChSel%s,name%s}');  
  c2=sdth.urn('Dock.Id.ci');Time=stack_get(c2,'curve','Time','g');
  if isfield(RO,'tmin')
   Time=fe_def('subdef',Time,Time.X{1}(:,1)>=RO.tmin&Time.X{1}(:,1)<=RO.tmax);
  end
  stack_set(c2,'curve','Time',Time);


elseif comstr(Cam,'specevt')
%% #ViewSpecEvt : extract frequency from peaks
 c13=iiplot(13,';');
 ob=handle(c13.ua.ob(1));
 Z=ob.ZData; if contains(c13.ua.YFcn,'log');Z=10.^Z;end
 [r2,i2]=max(Z,[],2);
 RO.gf=102; 
 figure(RO.gf);ga(1)=subplot(122);plot(r2,ob.YData);
 axis tight; %xlim([0 9000])
 cingui('plotwd',101,'@OsDic(SDT Root)',{'ImSw80','WrW49c','imgrid'});

 [r3,i3]=max(Z,[],1);ga(2)=subplot(121);cla(ga(2))
 r3=r3'; r3(r3<.05*max(r3))=NaN;
 r3(:,2)=ob.XData';r3(:,3)=ob.YData(i3)';r3=r3(:,[2 3 1]);
 r3(end+1,:)=NaN;
 C3=struct('X',{{(1:size(r3,1))',{'Time';'Freq';'Amp'}}}, ...
     'Xlab',{{'Index','Comp'}},'Y',r3,'name','IFreq');
 try; Time=c13.Stack{'Time'};key=Time.info.key;
     C3.name=[ key '-iFreq'];
     PA=sdtroot('PARAMVH');
     if isfield(PA,'DirScan')&&~isempty(PA.DirScan); 
      if isfield(PA.DirScan,'PAR')&&isfield(PA.DirScan.PAR,key)
       C2=PA.DirScan.PAR.(key);
       r2=interp1(C2.X{1},C2.Y,r3(:,1),'linear');
       C3.Y=[C3.Y r2]; 
       C3.X{2}(end+(1:size(r2,2)),1:size(C2.X{2},2))=C2.X{2};
      end
      PA.DirScan.iFreq.(key)=C3; 
      PA.DirScan.iFreq.(key)=C3; 
     end
 end


 r4=r3(isfinite(r3(:,3)),2);
 RO.ylim=[min(r4) max(r4)]; RO.ylim=RO.ylim+[-1 1]*diff(RO.ylim)/10;

 p=patch(ga(2),'faces',(1:length(r3)),'Vertices',r3,'FaceVertexCData',r3(:,3),'edgecolor','flat','facecolor','none','linewidth',2);
 set(ga(2),'xlim',ob.XData([1 end]),'ylim',ob.YData([1 end]))
 ylabel(get(get(c13.ga,'ylabel'),'string'))
 title(get(sdth.urn('ii_legend',13),'string'));
 set(ga(1),'yticklabel',[],'xticklabel',[],'position',[.7 .11 .2 .81],'ylim',RO.ylim);
 set(ga(2),'position',[.15 .11 .53 .81],'ylim',RO.ylim);
 setappdata(RO.gf,'curve',C3)
 cingui('objset',RO.gf,{'@OsDic',{'ImGrid'}})

elseif comstr(Cam,'spec');[CAM,Cam]=comstr(CAM,5);
%% #ViewSpec _tro : standardized viewing of base spectrogram 
c2=sdth.urn('Dock.Id.ci');Time=stack_get(c2,'curve','Time','g');
if isempty(Cam)
  m1=c2.Stack{'curData'};projM=c2.data.nmap.nmap;
  st=m1.views;st=st{sdtm.regContains(st,'ViewSpec')};
  if strncmp(st,'d_',2);Cb=sdtm.urnCb(st,projM);
  else;Cb=sdtm.urnCb(['$' st],projM); % replace from projM entry
  end
  fprintf('Calling %s\n',comstr(Cb,-30))
  feval(Cb{:});
  return
else
  [st,RO]=sdtm.urnPar(CAM,'{Spec%s}:{ci%g,jframe%g,ChSel%s,name%s,jPar%g}');  
end
  if ~isfield(RO,'Failed');RO.Failed={};end
  i1=~cellfun(@isempty,regexpi(RO.Failed,'[ft](min|max)'));
  if any(i1)
    RO.Spec=horzcat(RO.Spec,RO.Failed{i1});RO.Failed(i1)=[];
  end
  if isfield(RO,'jPar')&&RO.jPar
    r2=stack_get(c2,'curve','Split','g');
    if ~isempty(r2);Time=r2(:,:,RO.jPar);end
  end
  if ~isfield(RO,'nmap')
   try;RO.nmap=c2.data.nmap.nmap;catch;RO.nmap=vhandle.nmap;c2.data.nmap.nmap=RO.nmap;end
  end
  if carg<=nargin&&isfield(varargin{carg},'Y')
      Time=varargin{carg};carg=carg+1;
      if ~isfield(Time,'evt');iicom(c2,'curveinit','Time',Time);end%Not click
  elseif isempty(Time)
   Time=evalin('base','Time');iicom(c2,'curveinit','Time',Time);
  end
%spec=ii_signal('spectro BufTime 1 Overlap .8 fmax 10e3 fmin 800 -window hanning',Time);C2=spec.GetData; 
if isa(Time,'curvemodel');
    if ~isfield(RO,'jPar');RO.jPar=1;end
    if RO.jPar==0;%d_squeal('ViewSpec{Overlap .9 fmin 3000 fmax 3250 -window hanning,jPar0},nameVA')
        Time=Time.Source; 
    else
     Time=Time(:,:,RO.jPar);
    end
end
Time.X{1}=Time.X{1}-Time.X{1}(1);
spec=ii_signal(['spectro',RO.Spec],Time,struct,RO);
%C2=spec.GetData; 
r3=fe_def('cleanentry',spec.Source.Edit);
spec.Source.X{1}=spec.Source.X{1}+r3.BufTime/2;
%spec.Source.PlotInfo=ii_plp('PlotInfo2D -type "surface"',spec.Source); % To allow horz shift
spec.Source.PlotInfo=ii_plp('PlotInfo2D -type "image"',spec.Source); % To allow horz shift
C2=spec.Source;
if isequal(C2.Xlab{2},'DOF'); RO.type='fft'; C2.PlotInfo={};C2.name='FFT'; 
else; RO.type='spec';
end
if isfield(RO,'name');C2.name=RO.name;
elseif isfield(Time,'info')&&isfield(Time.info,'key')
    C2.name=Time.info.key; 
end
spec.Source=C2;

if isfield(Time,'evt')&&isfield(Time.evt,'urn')
   %% save instant input freq with possible trigger shift 
   st=Time.evt.urn; 
   if strncmpi(st,'sw',2); 
       RO.dt=diff(Time.X{1}(1:2));
       [~,RG]=sdtm.urnPar(st,sdtsys('UrnSig{Sw}'));RG.T=Time.X{1};
       RG=feval(sdtsys('@Sw'),RO,RG,Time.X{1},[]);
       % use spectrogram to realign times 
       Time.evt.f=RG.f(:);
       %i1=20:10:size(spec.Source.X{1},1)-10;
       %r2=spec(i1,:,2);[~,i2]=max(abs(r2),[],2);r2=spec.Source.X{2}(i2);
   end
end

% line(r3.X{1},r3.Y(:,1),max(zlim(c13.ga))*ones(size(r3.X{1})),'linewidth',3,'color','r')

%C2=spec.GetData; C2.X{1}=C2.X{1}-C2.X{1}(1);iicom(c2,'curveinit','spec',C2)
%c2.ua.YFcn='r3=abs(r3);';iiplot

iicom(c2,'initSqSig{Hbm.SqBase}'); 
if ~isfield(RO,'ci');RO.ci=13;end
new=ishandle(RO.ci);
c13=sdth.urn(sprintf('Dock.Id.ci.Clone{%i}',RO.ci));
if RO.ci==13; setappdata(13,'SdtName','Spec');end
if new % Place in same tile as iiplot 
 cingui('objset',c13,{'@Dock',{'Name','Id','tile',c2.opt(1)}})
end
iicom(c13,'curveinit','spec',spec)
if strcmpi(RO.type,'fft') % #guess_cycle_freq -3
 r2=sum(abs(C2.Y),2);if C2.X{1}(1)==0; r2(1:2)=0;end
 [~,i2]=max(r2);RO.f=C2.X{1}(i2); out=RO; 
else;% If spectro
 if ~isfield(RO,'Failed')||~any(strcmpi(RO.Failed,'zlog'));
  c13.ua.YFcn='r3=abs(r3);';iiplot(c13);
 end
 if isfield(C2,'Y')
  r2=sum(sum(abs(C2.Y),3),1);if C2.X{2}(1)==0; r2(1:2)=0;end
 else; % Update needed
  r2=specMax;
 end
 [~,i2]=max(r2);
 RO.f=C2.X{2}(i2); out=RO; 
end 

c2.os_('p.','ImToFigN','ImSw80','WrW49c');c13.os_('p.','ImToFigN','ImSw80','WrW49c');

iimouse('interacturn',13,menu_generation('interact.surf3d'));
nmap=c2.data.nmap.nmap;
out.SpecEdit=spec.Source.Edit; % Save handle to spec parameters
out.Time=Time;
nmap('SqLastSpec')=out; % Save context of last spec
if nargout==0; clear out;end


% figure(1);plot(C2.X{2},r2,':+')
elseif comstr(Cam,'par')
%% #ViewPar{fs2,f(p),xxx,cuName}
c2=sdth.urn('Dock.Id.ci'); nmap=c2.data.nmap.nmap;
[~,RO]=sdtm.urnPar(CAM,'{}{fs%ug,u%s,cu%s,ci%i,it%g}');
if ~isfield(RO,'Failed');RO.Failed={};end
if ~isfield(RO,'cu');RO.cu='Time';end
if isKey(nmap,RO.cu);Time=nmap(RO.cu); else; Time=c2.Stack{RO.cu};end;
RO.back=0;
if isfield(RO,'it')% d_squeal('viewpar{a(f,p),cuSqInstFreq,it4 7}')
   Time=fe_def('subdef',Time,@(x)x(:,1)>RO.it(1)&x(:,1)<RO.it(2));
end
omethod=sdth.eMethods.omethod;
if any(sdtm.Contains(lower(RO.Failed),'f(p)'))
  %% #ViewPar.f(p) Display instant freq as function of pressure -3
  [r1,i2,st]=omethod('xvec',Time,1,{'Pres','Freq'});
  gf=102;
  figure(gf); plot(r1(:,1),r1(:,2))
  xlabel(st(1));ylabel(st(2));
  cingui('plotwd',gf,'@OsDic(SDT Root)',{'ImSw80','WrW49c'});
  RO.back=1;

end
if any(sdtm.Contains(lower(RO.Failed),'f(t)'))
  %% #ViewPar.f(t) Display instant freq as function of time attempt to show wheel pos -3
  [r1,i2,st]=omethod('xvec',Time,1,{'Time','Freq'});
  gf=sdth.urn('figure(103).os{@Dock,{name,SqSig},name,103 f(t),NumberTitle,off}');
  if length(Time.Xlab)==3&&isequal(Time.Xlab{3},'hdof')
   figure(gf);clf;ga=get(gf,'CurrentAxes'); if isempty(ga);ga=axes('parent',gf);end
   y=[sum(abs(Time.Y(:,:,1)).^2,2);NaN];
   st1={[r1(:,1);NaN],[r1(:,2);NaN],y,y,'edgecolor','interp','tag','iFreq', ...
       'facevertexalphadata',y/max(y),'edgealpha','interp','linewidth',2};
   patch(st1{:},'parent',ga);
   set(ga,'alim',[-.05 1]);
   iimouse('on');
  else;
   figure(gf);plot(r1(:,1),r1(:,2))
  end

  xlabel(st(1));ylabel(st(2));
  axis tight; wheelPosLines(c2);
  cingui('plotwd',gf,'@OsDic(SDT Root)',{'ImSw80','WrW49c'});
  RO.back=1;

end
if any(sdtm.Contains(lower(RO.Failed),'a(f)'))
  %% #ViewPar.a(f) xxx -3
  [r1,i2,st]=omethod('xvec',Time,1,{'Time','iFreq'});
  [r1,st]=getAmp(r1,Time,st);r1(:,4:end)=[];
  if strncmp(RO.Failed{2},'{',1)
    [~,RP]=sdtm.urnPar(RO.Failed{2},'{}{lp%g}');
    if ~isfield(RP,'dt');RP.dt=diff(r1(1:2),1);end
    RP.lp=sdtpy(sprintf('lowpass{8,%f %f,pa 5 .8,dososfiltfilt}',RP.lp,1/RP.dt));
    r1(:,4:5)=r1(:,2:3);
    r1(:,2)=RP.lp*r1(:,2);
    %r1(:,3)=(RP.lp*abs(r1(:,3)));
    r1(:,3)=10.^(RP.lp*log10(abs(r1(:,3))));
  end
  gf=sdth.urn('figure(104).os{@Dock,{name,SqSig},name,104 a(f),NumberTitle,off}');
  figure(gf);clf;ga=get(gf,'CurrentAxes'); if isempty(ga);ga=axes('parent',gf);end
   st1={[r1(:,2);NaN],[abs(r1(:,3));NaN],[r1(:,1);NaN],[r1(:,1);NaN],'edgecolor','interp','tag','iFreq', ...
       'linewidth',2};
   st2=st{3,1};if iscell(st2);st2=st2{1};end
   h(1)=patch(st1{:},'parent',ga,'DisplayName',st2);
   if size(r1,2)>3
    h(2)=line(r1(:,4),abs(r1(:,5)),'linestyle',':', ...
        'marker','none','color',[.5 .5 .5 .5],'markersize',1, ...
        'DisplayName',sprintf('Raw %s',st2));
    set(h(1),'linewidth',4)
   end
   set(ga,'alim',[-.05 1],'yscale','log');
   iimouse('on');
  xlabel(st{2,1});ylabel(st{3,1});h(3)=colorbar;
  h(3).Label.String=st{1};ii_plp('colormapband',parula(7));
  axis tight; 
  cingui('plotwd',gf,'@OsDic(SDT Root)',{'ImSw80','WrW49c'});
  RO.back=1;
end
if any(sdtm.Contains(lower(RO.Failed),'a(f,p)'))
  %% #ViewPar.a(f,p) Display instant freq as function of time attempt to show wheel pos -3
  [r1,i2,st]=omethod('xvec',Time,1,{'Pres','iFreq'});
  [r1,st]=getAmp(r1,Time,st);

  gf=104; 
  figure(gf);clf;ga=get(gf,'CurrentAxes'); if isempty(ga);ga=axes('parent',gf);end
   st1={[r1(:,2);NaN],[r1(:,3);NaN],[r1(:,1);NaN],[r1(:,1);NaN],'edgecolor','interp','tag','iFreq', ...
       'linewidth',2};
   patch(st1{:},'parent',ga);
   set(ga,'alim',[-.05 1]);
   iimouse('on');
  xlabel(st{2,1});ylabel(st{3,1});h=colorbar;
  h.Label.String=st{1};ii_plp('colormapband',parula(7));
  axis tight; wheelPosLines(c2);
  cingui('plotwd',gf,'@OsDic(SDT Root)',{'ImSw80','WrW49c'});
  RO.back=1;
end
if any(sdtm.Contains(lower(RO.Failed),'a(t)'))
  %% #ViewPar.a_t Display amplitude as function of time and instant freq -3
  [r1,i2,st]=omethod('xvec',Time,1,{'Time','iFreq'});
  [r1,st]=getAmp(r1,Time,st);
  gf=sdth.urn('figure(105).os{@Dock,{name,SqSig},name,105 a(t),NumberTitle,off}');

  figure(gf);clf;ga=get(gf,'CurrentAxes'); if isempty(ga);ga=axes('parent',gf);end
   st1={[r1(:,1);NaN],[abs(r1(:,3));NaN],[r1(:,2);NaN],[r1(:,2);NaN],'edgecolor','interp','tag','iFreq', ...
       'linewidth',2};
   patch(st1{:},'parent',ga);
   set(ga,'alim',[-.05 1]);
   iimouse('on');
  xlabel(st{1,1});ylabel(st{3,1});h=colorbar;
  h.Label.String=st{2};ii_plp('colormapband',parula(7));
  axis tight; wheelPosLines(c2);
  cingui('plotwd',gf,'@OsDic(SDT Root)',{'ImSw80','WrW49c'});
  if 1==2% size(r1,1)>3
   gf=sdth.urn('figure(106).os{@Dock,{name,SqSig},name,106 q2R(t),NumberTitle,off}');
   figure(gf);clf;ga=get(gf,'CurrentAxes'); if isempty(ga);ga=axes('parent',gf);end
   plot(r1(:,1),abs(r1(:,3:end)));axis tight
   set(ga,'ylim',max(get(ga,'ylim'))*[.01 1])
   xlabel(st{2,1});ylabel(strrep(st{3,1},'q1','qj'));
   legend('q_{2R}','q_{3R}','location','best')
   cingui('objset',gf,{'@OsDic',{'ImSw80{@line,""}'},'@PlotWd',{'@OsDic','WrW49c'}})
   wheelPosLines(c2);

  end
  RO.back=1;
end
if any(sdtm.Contains(lower(RO.Failed),'at'))
  %% #ViewPar.at : harmonic modulation d_squeal('viewpar{at,cuSqInstFreq}',C1) -3
  gf=202;figure(gf); 
  [r1,i2,st]=omethod('xvec',Time,1,{'Time','iFreq','TR'});
  % d_squeal('viewpar{at{5000},cuSqInstFreq}',C1)
  [~,r2]=sdtm.urnPar(RO.Failed{sdtm.Contains(lower(RO.Failed),'at')},'{N%g}{harm%g,dh%g,der%g}');
  if ~isfield(r2,'dh');r2.dh=1;end
  r1=interp1(r1(:,1),r1,linspace(r1(1),r1(end,1),r2.N),'linear','extrap');
  phi=r1(:,3)*pi/2;t=r1(:,1);
  RO.harm=cellfun(@(x)str2double(x(2:end)),Time.X{3});
  r3=interp1(Time.X{1}(:,1),permute(Time.Y,[1 3 2]),r1(:,1));
  if ~isfield(r2,'der')
  elseif r2.der==1 % v2a 
   for j1=1:length(RO.harm); r3(:,j1,:)=r3(:,j1,:)*RO.harm(j1);end
  else; 
      error('Not implemented')
  end
  if isfield(r2,'harm');RO.harm(~ismember(RO.harm,r2.harm))=0;end
  at=sum(r3.*exp(-1i*phi*(RO.harm(:)'-r2.dh)),2);
  C2=struct('X',{{r1(:,[3 2 1]),Time.X{2},{'at';'ct1'}}}, ...
      'Xlab',{{Time.Xlab{1}([3 2 1],:),Time.Xlab{2},'Harm'}},'Y',reshape(at,length(t),size(at,3)));
  C2.Y(:,:,2)=reshape(real(r3(:,RO.harm==1,:).*exp(-1i*phi)),length(t),size(at,3));
  C2.Ylab=202;C2.name='at';
  gf=202;c202=sdth.urn(sprintf('Dock.Id.ci.Clone{%i}',gf));
  iicom(c202,'curveinit','At',C2)
  %plot((phi),abs(at))
  %plot(at);xlabel('Real(at)');ylabel('Imag(at)');grid on
  cingui('plotwd',gf,'@OsDic(SDT Root)',{'ImToFigN','ImSw80','WrW49c'});

end

if isfield(RO,'ci')&&~isempty(RO.ci)
  %% ViewPar.ci20 display curve in a clone
  c12=sdth.urn(sprintf('Dock.Id.ci.Clone{%i}',RO.ci));
  if ~isfield(Time,'name');Time.name=RO.cu;end
  setappdata(c12.opt(1),'SdtName',RO.cu);
  c12.os_('p.','FnIy','ImToFigN','ImSw80','WrW49c');%
  iicom(c12,'curveInit',Time);

end
if RO.back; return;end

%[~,r2]=sdtm.urnPar('a{polar,fs2}','{}{fs%ug,u%s}');
C1=struct('X',{{Time.X{1}(:,1),Time.Xlab{1}(2:end,:)}},...
    'Xlab',{{Time.Xlab{1}(1,:),'Comp'}},'Y',Time.X{1}(:,2:end),'Ylab',2);
if isfield(Time,'name'); C1.name=Time.name;end
if isfield(Time,'info')&&isfield(Time.info,'key'); C1.name=Time.info.key;end
if isfield(RO,'fs') % Allow decimation
 C1=sdtpy.decimate(C1,round(1/diff(C1.X{1}(1:2))/RO.fs));
end
if isfield(Time,'ID');C1.ID=Time.ID;end
  % c11=iiplot(11);iicom(c11,'curveinit','Time',C1);
c12=sdth.urn('Dock.Id.ci.Clone{12}');setappdata(12,'SdtName','Par(t)');
c12.os_('p.','FnIy','ImToFigN','ImSw50','WrW49c');%
if ~isempty(C1.Y)
 iicom(c12,'curveinit','Parameters',C1);
 r1=[.13 .2 .8 .75];c12.ax(1,6:9)=r1;set(c12.ga,'position',r1);iiplot(c12);
end
if any(strcmpi(RO.Failed,'polar'))
  c12=get(12,'userdata');iicom(c12,'polar Comp(2)')
  C1=c12.Stack{'Parameters'};
  t=C1.X{1}(:,1);
  i1=sdtm.indNearest(t,(t(1):10:t(end))');
  go=handle(c12.ua.ob(1));r2=[go.XData(i1)',go.YData(i1)'];r2(:,3)=max(get(c12.ga,'zlim'));
  l=patch('parent',c12.ga,'vertices',r2,'faces',(1:length(i1))');
  set(l,'facecolor','none','marker','*','markersize',12, ...
      'FaceVertexCData',t(i1),'MarkerEdgeColor','flat','MarkerFaceColor','flat','linewidth',2)
  
end


 elseif comstr(Cam,'instfreq')
%% #ViewInstFreq : estimate instant frequency and modulation

if carg>nargin||comstr(Cam,'instfreq{') 
 c2=sdth.urn('Dock.Id.ci'); 
 Time=c2.Stack{'Time'};% (SqLastSpec).Time is preemptive 
 projM=c2.data.nmap.nmap; 
 if ~isKey(projM,'SqLastSpec');warning('Missing SqLastSpec entry');return;end
 RO=projM('SqLastSpec');
 [~,r2]=sdtm.urnPar(CAM,['{}{dmBand%g,aeBand%g,harm%g,do%s,tclip%g,' ...
     'clipBand%g,f%g,chRef%s,ifBand%g,ifSat%g,jPar%g,iu%g,hf%g,ci%i}']);
 RO=sdth.sfield('addmissing',r2,RO);
 if carg<=nargin; r2=varargin{carg};carg=carg+1;
     RO=sdth.sfield('addmissing',r2,RO);
 end
else
 Time=varargin{carg};carg=carg+1;projM=[];
 RO=varargin{carg};carg=carg+1;
end
hfs=ii_signal('@sqSig');
if ~isfield(RO,'Time'); RO.Time=Time;else; Time=RO.Time; end
if ~isfield(RO,'jPar')||RO.jPar==0; 
else 
    RO.Time=Time(:,:,RO.jPar);RO.Time.X{1}=RO.Time.X{1}-RO.Time.X{1}(1);
end
if ~isfield(RO,'harm');RO.harm=(0:5);end

if ~isfield(RO,'do');RO.do='';end
%try
    [out,RB]=hfs('obsPha',RO); % sdtweb ii_signal obspha
%catch err
%    sdtm.toString(err)
%    out=[];return
%end
RO=RB; if nargout>1; out1=RO;end

if isfield(RO,'ci')&&~isempty(RO.ci)
    %'xxx need to adjust how to display freq if not spectro'
    ci=iiplot(RO.ci,';');r2=struct('X',{{out.X{1}(:,1),{'Freq'}}},'Y',out.X{1}(:,2));
    ci.Stack{'spec'}.ID=struct('po',r2,'marker','xy'); 
    %if length(RO.f)==1;ci.ua.axProp={'ylim',RO.f*[.98 1.02]};
    %else; ci.ua.axProp={'ylim',[min(RO.f)*.98 max(RO.f)*1.02]};
    %end
    iiplot(ci);
end
if sdtm.Contains(lower(RO.do),'autoseg')
 r4=abs(out.Y(:,2));it=find(r4>max(r4)*RO.EvtTol);
 C3=feval(process_r('@SigEvt'),'init',Time,struct('it',it,'do','AutoSeg'));
 C3.DimPos=1:3; c2=sdth.urn('Dock.id');
 r4=Time.X{1}(C3.Source.Range.val(:,1:2));
 c2.Stack{'Time'}.ID=struct('po',r4,'marker','band');
 stack_set(c2,'curve','EvtTime',C3);
end
omethod=sdth.eMethods.omethod;
if ~isempty(projM) % Store in standard map 
 projM('SqLastSpec')=RO;
 projM('SqInstFreq')=out;% Store result 
end
st=regexprep(RO.do,'ReEstY[,]?',''); %'{ReEstY}'
if isequal(st,'{}');st='';end
if ~isempty(st);st=strrep(st,'}',',cuSqInstFreq}');end
if ~isempty(st); d_squeal(['viewpar' st]);end
%if sdtm.Contains(lower(RO.do),'f(t)'); d_squeal('viewpar{f(t),cuSqInstFreq}');end
 elseif comstr(Cam,'bandpass')
 %% #viewBandPass
 RO=varargin{carg};carg=carg+1;
 c2=sdth.urn('Dock.Id.ci');

 TA=fe_def('subchcurve',c2.Stack{'Time'},{'Channel',RO.ch}); 
 TA=fe_def('subdef',TA,@(x)x(:,1)>RO.tmin(1)&x(:,1)<RO.tmin(2));
 st={'Time'};r1=1/diff(TA.X{1}(1:2));RO.fmin(RO.fmin>r1)=r1;
 for j2=1:size(RO.fmin,1)
  TB=ii_mmif(sprintf('bandpass fmin %g %g -struct',RO.fmin(j2,:)),TA); 
  TA.Y(:,j2+1)=TB.Y(:,1);
  st1=sprintf('[%.0f - %.0f] Hz',RO.fmin(j2,:));
  st2=sprintf('[%.2f - %.2f] kHz',RO.fmin(j2,:)/1000);
  if length(st1)>length(st2);st1=st2;end
  st{j2+1}=st1;
 end
 TA.X{2}=repmat(TA.X{2}(1,:),size(TA.Y,2),1);
 %TA=fe_def('subdef',TA,@(x)x(:,1)>61.7&x(:,1)<64.3);TA.X{2}=repmat(TA.X{2}(1,:),3,1)

 gf=301;figure(gf);clf; cmap=lines(size(TA.Y,2));
 sdth.omethod('cleanvec',TA,'{x,1}{y,1}{gf301}');
 for j2=2:size(TA.Y,2)
  h=line(TA.X{1}(:,1),TA.Y(:,j2),'color',cmap(j2,:));
 end
 axis tight
 cingui('plotwd',gf,'@OsDic(SDT Root)',{'ImToFigN','ImSw80','WrW49c'});
 legend(st{:},'location','best')
 if isfield(RO,'ID')&&~isequal(RO.ID,0)
  ii_plp(c2.Stack{'Time'}.ID)
 end
 elseif comstr(Cam,'occ')
 %% #viewOcc : occurence tracking 
 c2=sdth.urn('Dock.Id.ci');projM=c2.data.nmap.nmap;
 C2=projM('SqInstFreq');RO=projM('SqLastSpec');
%r2=C2.Y;r2=r2./r2(:,1);r2=z*r2;r2=r2.*sum(abs(r2).^2,2).^(-.5);
[~,RO]=sdtm.urnPar(CAM,'{tref%ug}:{xlim%ug}');
RO.iref=sdtm.indNearest(C2.X{1}(:,1),RO.tref);
RO.lp=sdtpy('lowpass{8,10 17k,pa 5 .8,dososfiltfilt}');

r2=C2.Y;r2=r2./r2(:,1);r2=complex(RO.lp*real(r2),RO.lp*imag(r2));
r2=r2.*sum(abs(r2).^2,2).^(-.5);
r3=(abs(r2*r2(RO.iref,:)'));
gf=107;figure(gf); 
[r1,i2,st]=sdth.omethod('xvec',C2,1,{'Time','iFreq'});

cingui('plotwd',gf,'@OsDic(SDT Root)',{'ImToFigN','WrW49c'});
%gf=103:110; gf(~ishandle(i2))=[];
cingui('objset',gf,{'@dock',{'name','Squeal'},'@OsDic','ImSw80{@line,""}'})
h=plot(C2.X{1}(:,1),r3,'linewidth',2);xlabel(st{1});ylabel('MAC')
st1=cellfun(@(x)sprintf('%.0fs',x),num2cell(RO.tref),'uni',0);
ylim([.5 1.01]);grid on
cmap=lines(length(h));
for j1=1:length(RO.iref);
    line(RO.tref(j1),r3(RO.iref(j1),j1),'color',cmap(j1,:), ...
        'MarkerFaceColor',cmap(j1,:),'marker','d','markersize',12)
end
legend(h,st1{:},'location','best')

r4=r3;r4(r4<.9)=NaN;
c13=get(13,'userdata');ga=handle(c13.ga);
if ~isfinite(ga.YLim(1));iiplot(c13);end
clu=[];delete(findobj(ga,'type','line'))
for j1=1:size(r3,2)
 r5=r1(:,2);prop={};
 r5=ones(size(r5))*(ga.YLim*[.95;.05]);prop={'linewidth',12};
 r5(r3(:,j1)<.8|r3(:,j1)<max(r3,[],2))=NaN;
 clu(j1)=line(r1(:,1),r5,'color',cmap(j1,:),'parent',ga,prop{:});
 line(RO.tref(j1),r1(RO.iref(j1),2),'color',cmap(j1,:), ...
      'MarkerFaceColor',cmap(j1,:),'marker','d','markersize',12,'parent',ga)
 %,ones(size(r1,1),1)*max(ga.YLim)
end
cingui('objset',13,{'@OsDic','ImSw80{@line,""}'})

%ylim([0 1]);xlim([25 280]);grid on


 elseif comstr(Cam,'ci')
 %% #viewCi : dock init
 RO=varargin{2};Time=RO.Time; RO.out={};
  for j1=1:length(RO.ci)
   if nargout>0
   elseif j1==1;
       c2=sdth.urn('Dock.Id.ci');
       if ~isempty(c2.Stack{'Time'});stack_set(c2,'curve','Time',Time);iiplot(c2)
       else;iicom(c2,'curveinit','Time',Time);
       end
   else; ci=sdth.urn(sprintf('Dock.Id.ci.Clone{%i}',RO.ci(j1)));
   end
   ci=iiplot(RO.ci(j1),';');
   ci.os_('p.','ImToFigN','ImSw80','WrW49c');
   eval(sprintf('c%i=ci;',RO.ci(j1)));RO.out{end+1}=sprintf('c%i',RO.ci(j1));
  end  
 if nargout==0;  eval(iigui(RO.out,'SetInBaseC'));end

 else; error('View%s unknown.',CAM)
 end

 %% #Load: load and associated common model operations
elseif comstr(Cam,'load');[CAM,Cam]=comstr(CAM,5);
 if comstr(Cam,'ods'); [CAM,Cam]=comstr(CAM,4);
  %% #LoadODS: [data,d1,d2]=d_squeal('loadods',fsvd,RA); -2
  % RA with fields adatped to polytec read+ f/ch/split/basSet
 
  f1=varargin{carg}; carg=carg+1;
  if carg>nargin;RA=struct;  else; RA=varargin{carg}; carg=carg+1;  end
  if ~isfield(RA,'noRemerge')
   [CAM,Cam,RA.noRemerge]=comstr('-noremerge',[-25 3],CAM,Cam);
  end
  
  % Reload from svd file (or .mat)
  if strncmp(fliplr(f1),'ffu.',4)||strncmp(fliplr(f1),'vnu.',4) % skip non polytec files
   if isfield(RA,'fmat')&&exist(RA.fmat,'file')
    load(RA.fmat,'wire','d1')
   else
    UFS=ufread(f1);
    wire=UFS(1); d1=UFS(2);% xxx need refinement if more entries are present
    if isfield(d1,'dof')
     d1=struct('DOF',d1.dof(:,1),'def',d1.res.','data',d1.po,'fun',d1.fun,...
      'info',d1.header);
    end
    wire.tdof=fe_sens('tdof',feutil('getdof',wire.Node(:,1),[.01;.02;.03]));
    if isfield(RA,'fmat'); sdtm.save(RA.fmat,'wire','d1'); end
   end
   if isfield(RA,'rmPoints')&&~isempty(RA.rmPoints) % additionnal point removal
    wire=fe_sens('MeshSub',wire,RA.rmPoints(:));
    d1=feutilb('placeindof',wire.tdof(:,1),d1);
   end
   out=struct('TEST',wire); out1=d1; out2=[];
  else % polytec files
   RA.FastScan=sdtm.regContains(f1,'(TimeScan|FastScan)');
   [out,out1,out2]=polytec(['ODSLoad' CAM],f1,RA);
  end
  
  % post split, issue with merging
  if isfield(RA,'split')
   mo1=out.TEST; out.TEST=cell(1,length(RA.split));
   for j1=1:length(RA.split)
    out.TEST{j1}=fe_sens('MeshSub',mo1,RA.split{j1},struct('noKeep',1));
    if isfield(RA,'mirror')&&RA.mirror(j1)
     %out.TEST{j1}.tdof(:,1)=out.TEST{j1}.tdof(:,1)+.06;
     out.TEST{j1}.tdof(:,3:5)=-out.TEST{j1}.tdof(:,3:5)
    end
   end
  end
    
  % apply basSet if any
  if isfield(RA,'basSet')
   if iscell(out.TEST)
    for j1=1:length(out.TEST);
     out.TEST{j1}=fe_sens('MeshProject',out.TEST{j1},RA.basSet{j1});
    end
   else; 
    if iscell(RA.basSet); RA.basSet=RA.basSet{1}; end
    out.TEST=fe_sens('MeshProject',out.TEST,RA.basSet);
   end
  end
  % remerge if needed
  if ~RA.noRemerge&&iscell(out.TEST)
   TEST=out.TEST{1};
   for j1=2:length(out.TEST)
    TEST.Node=[TEST.Node;out.TEST{j1}.Node];
    TEST.Elt=[TEST.Elt; out.TEST{j1}.Elt];
    TEST.tdof=[TEST.tdof;out.TEST{j1}.tdof];
   end
   out.TEST=TEST;
  end
  
  if isfield(RA,'InFEM')
   if isstruct(out.TEST); out.TEST.InFEM=RA.InFEM;
   elseif isscell(out.TEST)
    for j1=1:length(out.TEST);
     if isstruct(out.TEST{j1}); out.TEST{j1}.InFEM=1; end
    end
   end
  end

  if ~iscell(out.TEST)&&~isempty(out1)
   out1=feutilb('placeindof',out.TEST.tdof,out1);
   if ~isfield(out1,'tdof')&&isfield(out1,'DOF')&&size(out1.DOF,2)==5
    out1.tdof=out1.DOF; out1=rmfield(out1,'DOF');
   end
   if ~isempty(out2); out2=feutilb('placeindof',out.TEST.tdof,out2);
    if ~isfield(out2,'tdof')&&isfield(out2,'DOF')&&size(out2.DOF,2)==5
     out2.tdof=out1.DOF; out1=rmfield(out2,'DOF');
   end
   end
  end
  
 elseif comstr(Cam,'def'); [CAM,Cam]=comstr(CAM,4);
  %% #LoadDef: complex mode loading common utility -2
  %d_squeal('LoadDef',cf,RO);
  cf=varargin{carg}; carg=carg+1;
  RO=varargin{carg}; carg=carg+1;
  if ~isfield(RO,'up'); RO.up=0; end
  if ~isfield(RO,'up1'); RO.up1=RO.up; end
  if isempty(cf); [wd,f1,ex1]=fileparts(RO.f1); f1=[f1 ex1];
  else; wd=cf.mdl.wd; f1=cf.mdl.file;
  end
  fdef=fullfile(wd,strrep(f1,'.mat','_defstab.mat'))
  if ~isempty(cf)&&~isempty(cf.Stack{'dcpx'})&&~RO.Reset&&~RO.ResetDef; def=cf.Stack{'dcpx'};
  elseif RO.Reset||RO.ResetDef||~exist(fdef,'file')
   if isempty(cf); error('file %s not generated, need model',fdef); end
   cf.Stack{'dcpx'}=[]; %cf.Stack{'info','EigOpt'}=[5 RO.nM 1e3];
   ofact('mklserv_utils -silent')
   def=d_squeal('SolveModes',cf.mdl);
   def.label=sprintf('stab_up%i',RO.up1);
   sOpt=sdtdef('V6Flag');
   save(fdef,'def',sOpt{:});
  elseif RO.hdf; def=sdthdf('hdfreadref-level0',fdef,'/def');
  else; load(fdef,'def');
  end
  if isempty(cf); out=def;
  else
   stack_set(cf.mdl,'curve','dcpx',def);
   %   cf.Stack{'curve','dcpx'}=def;
   if nargout>0; out=cf; out1=def; end
  end
  
 elseif comstr(Cam,'modif'); [CAM,Cam]=comstr(CAM,6);
  %% #LoadModif: apply saved modifications to model -2
  % at the moment, saved morphed components
  cf=varargin{carg}; carg=carg+1;
  RO=varargin{carg}; carg=carg+1;
  root=varargin{carg}; carg=carg+1;
  
  if iscell(root); im=[root(:); {''}];
  else
   im=strfind(RO.modif,'+');
   if isempty(im); im=[0 length(RO.modif)+1];
   else; im=[0 im length(RO.modif)+1];
   end
  end
  for jm=1:length(im)-1
   if iscell(im); fm=im{jm}
   else;   fm=fullfile(cf.mdl.wd,[root '_' RO.modif(im(jm)+1:im(jm+1)-1) '.mat'])
   end
   if ~exist(fm,'file');
    stm=RO.modif(im(jm)+1:im(jm+1)-1);
    if stm(1)=='@'; stm=stm(2:end); end
    stm=sdtm.urnCb(stm);
    cf=feval(stm{:});
   else
    %error('Modif unknown %s',RO.modif); end
    load(fm,'mo3')
    opt=stack_get(mo3,'info','optMesh','get');
    if isfield(opt,'opt'); opt=opt.opt; end
    if iscell(opt)
     for j1=1:size(opt,1)
      if ischar(opt{j1,1})
       if strncmpi(opt{j1,1},'rm',2)
        cf.mdl.Elt=feutil('removeelt withnode{nodeid}',cf.mdl,opt{j1,2}(:,1));
       elseif strncmpi(opt{j1,1},'add',3)
        st=comstr(opt{j1,1},4);
        switch lower(st)
         case '-combinemodel'
          cf.mdl=feutilb('CombineModel',cf.mdl,stack_rm(mo3,'info','optMesh'));
         case '-combine-autoconn'
          cf.mdl=feutilb('CombineModel',cf.mdl,stack_rm(mo3,'info','optMesh'));
          cf.mdl=fe_caseg('ConnectionSurfaceAuto',cf.mdl,opt{j1,2});
          %-aTol1 -matchS -dens -KpAuto -maxdist.001 -maxNdist.0001 -radius1',....
          % cf.mdl,'kbthread','proid14','proid13')
         case '-stack'
          cf.mdl=stack_set(cf.mdl,opt{j1,2});
         otherwise
        end
       else
        cf.mdl=fe_shapeoptim('AdjustLoop exact midInterp ',cf.mdl,opt(j1,2:3));
       end
      else
       cf.mdl=fe_shapeoptim('AdjustReplay',cf.mdl,opt(j1,:));
      end
     end
    else
     cf.mdl.Elt=feutil('removeeltwithnode{inelt{withnode{nodeid}}}',cf.mdl,opt);
    end
   end
  end
  if isfield(RO,'modif')&&~isempty(RO.modif)&&isfield(cf.mdl,'name')
   cf.mdl.name=[cf.mdl.name '_' RO.modif];
   cf.mdl.file=[cf.mdl.name '.mat']; % xxx not saved at the moment
  end
  
  ouf=cf; out1=RO;
  
 elseif comstr(Cam,'time'); [CAM,Cam]=comstr(CAM,5);
 %% #LoadTime : load a long time -2
 if comstr(CAM,'(')
   % d_squeal('loadTime(cbi20b,TimeScan_3100Hz_laser_1.mat,ci2)')
   if carg<=nargin;RO=varargin{carg};carg=carg+1;
   else
    [st,RO]=sdtm.urnPar(CAM,'{LoadFcn%s,fname%s}:{ci%g}');  
   end
 elseif comstr(Cam,'{')
   [st,RO]=sdtm.urnPar(CAM,'{}:{ci%g}');  
   nmap=varargin{carg}; carg=carg+1;
   RO.Time=fe_def('def2curve',nmap('CurTime'));
 else;RO=varargin{carg};carg=carg+1; 
 end
 if isfield(RO,'Time'); Time=RO.Time;wire=[];i1=-1;
 elseif ~isempty(RO.LoadFcn)&&~strcmp(RO.LoadFcn,'ufread')
  [FileName,i1]=feval(RO.LoadFcn,'pwd',RO.fname);[~,~,RO.ext]=fileparts(FileName);
 elseif exist(RO.fname,'file');FileName=RO.fname;i1=1; [~,~,RO.ext]=fileparts(FileName);
 else; error('File not found');
 end
 RO.out={'Time','wire'};
 if i1==-1 % Given in field
 elseif i1&&~strcmpi(RO.ext,'.svd') % File exist
     if isequal(RO.LoadFcn,'ufread');r1=ufread(FileName);
     elseif strcmpi(RO.ext,'.wav')
       [r2,fs]=audioread(FileName); t=1/fs:1/fs:1/fs*length(r2);
       r1=struct;r1.Time=struct('X',{{t(:),{'mic'}}},'Xlab',{{'Time','Sens'}},'Y',r2(:));
     elseif isfield(RO,'Failed') % xxx third arg to load a given variable only
      r1=load(FileName,RO.Failed{1}); 
     else;r1=load(FileName);
     end
     if isfield(r1,'Time');Time=r1.Time;
     elseif isfield(r1,'XF'); Time=r1.XF;if isfield(r1,'TEST');wire=r1.TEST; end
     else; st=fieldnames(r1);Time=r1.(st{1});
     end
     if isfield(r1,'wire');wire=r1.wire;else;wire=[];end
 else
     %% revise fasttime storage. 
     f2=strrep(RO.fname,'.mat','.svd');
     if exist(f2,'file'); 
     elseif ~isempty(RO.LoadFcn)
      f2=feval(RO.LoadFcn,'pwd',f2);%'17TestPolytec500\combi_time.svd');
     end
     if ~exist(f2,'file'); error('''%s'' Not found',f2);
     end
     r1=polytec('ReadList',f2);
     r2=struct('iSignal',1,'pointdomain','Time','FastScan',1);
     r2.iChannel=unique(vertcat(r1{2:end,2}));
     Time=polytec('ReadSignal',f2,r2);
     % Ref 3,4,5 are a triax with default and Ref 6 has a amplitude issue
     % xxx cleanup
     % -> only keep Ref 1 and 2
     %r1.X{2}(6:end,:)=[];r1.Y(:,6:end)=[];r1.Range.lab(6:9)=[];r1.Range.val(:,6:9)=[];
     r1=Time.Source;
     Time.Source=r1;    
     wire=polytec('ReadMesh',f2);
     save(strrep(FileName,'.svd','.mat'),'Time','wire');
 end
 if isfield(RO,'cf');cf=comgui('guifeplot',RO.cf);cf.mdl=wire; end
 if isfield(RO,'ci') % Initialize ci
  RO.Time=Time;d_squeal('ViewCi',RO);
 end
 if nargout==0;  eval(iigui(RO.out,'SetInBaseC'))
 elseif nargout==1;out=struct('Time',Time,'wire',wire);
 else; out=Time;out1=wire;
 end
 else; error('Load %s unknown.',CAM);
 end
 
elseif comstr(Cam,'fixdiscrot')
 %% #FixDiscRot: set a disc rotation penalized constraint to represent wheel momentum
 model=varargin{carg}; carg=carg+1;
 if carg<=nargin; RO=varargin{carg}; carg=carg+1; else; RO=struct; end
 [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
  'discsel("setname disc"#%s#"disc selection")' ...
  'dhubsel("setname disc"#%s#"disc + hub selection")' ...
  'rTol(1.1#%g#"bolt radius tolerance")' ...
  'nTol(1e-3#%g#"node detection tolerance")' ...
  'radVarTol(1e-3#%g#"bolt radius variation")' ...
  'GivenRad(#%g#"given bolt radius to look for")' ...
  'krot(1e8#%g#"rotation stiffness")' ...
  'name("discrot"#%s#"name root")' ...
  'keepDTop(#3#"to also keep top disc surface")' ...
  '-getH(#3#"just to get bolt holes")' ....
  ],{RO,CAM}); Cam=lower(CAM);
 
 if carg<=nargin
  r1=varargin{carg}; carg=carg+1; % axis from orig
  r3=varargin{carg}; carg=carg+1; % master nodeids
  if carg<=nargin; n1=varargin{carg}; carg=carg+1; else; n1=0; end % slave pos offset from orig
 else
  % recover rotation axis (assume Abaqus input for the moment)
  r1=stack_get(model,'info','Motion','get');
  %  r1=r1{2,3}.axis;
  try
   in2=cellfun(@ischar,r1(:,2)); in2(1)=false;
   r1(in2,2)=cellfun(@(x)feutil(sprintf('findnode setname "%s"',x),model),r1(in2,2),'uni',0);
  catch; if ~isfield(RO,'axis'); sdtw('_nb','motion handling problem'); else; r1=[]; end
  end
  % get disc nodes
  mo1=sdth.GetData(model);
  mo1.Elt=feutil(sprintf('selelt %s & selface',RO.discsel),mo1);
  mo1.Node=feutil('getnodegroupall',mo1);
  % get disc top surface based on motion normal
  n1=mo1.Node;
  if isempty(r1); r1=RO.axis;
  else % based on abaqus motion command
   r1=r1{1+find(cellfun(@(x)any(ismember(x,n1(:,1))),r1(2:end,2)),1,'first'),3}.axis;
  end
  n1=bsxfun(@minus,n1(:,5:7),r1(1:3));
  r2=basis(r1(4:6)-r1(1:3),[1 1 1]);
  n1=n1*r2; i1=n1(:,1)>=(1-RO.nTol)*max(n1(:,1)); % nodes in local rot basis in top surface
  mo1.Elt=feutil('selelt innode',mo1,mo1.Node(i1,1));
  mo1.Node=feutil('getnodegroupall',mo1);
  % detect bolt holes, get holes in surface and get ones with equal largest radius
  mo2=nl_mesh('holegroups',mo1);
  [EGroup,nGroup]=getegroup(mo2.Elt); r2=cell(1,nGroup);
  for jGroup=1:nGroup
   n2=feutil(sprintf('getnode group%i',jGroup),mo2);
   r2{jGroup}=feutilb('geofindcircle',mo2,struct('nodes',n2(:,1)'));
  end
  r2=cellfun(@(x)x{1},r2,'uni',0); r3=cellfun(@(x)x.radius,r2);
  [r4,i4]=sort(r3,'descend');
  if isempty(RO.GivenRad)
   r2=r2(abs(r3-max(r4(find(abs(diff(r4))<RO.radVarTol))))<=RO.radVarTol);
  else; r2=r2(abs(r3-RO.GivenRad)<RO.radVarTol);
  end
  % recover wheel center
  r3=basis(r1(4:6)-r1(1:3),[1 1 1]);
  n4=r1(1:3)'+max(n1(:,1))*r3(:,1); % center node
  if RO.getH; out=r2; out1=n4; return; end
  % select nodes around bolt holes
  r3=r2;
  for j2=1:length(r2)
   n2=feutil(sprintf(['findnode inelt{%s} & '...
    'cyl <= %.15g  o %.15g %.15g %.15g  %.15g %.15g %.15g -30 30'],...
    RO.dhubsel,r2{j2}.radius*1.1,r2{j2}.Origin,r2{j2}.normal),model);
   r3{j2}=n2;
  end
  r3=cat(1,r3{:});
 end
 % keep disc top surface, or not
 if RO.keepDTop;  r3=unique([r3;mo1.Node(:,1)]); end
 
 % add a slave node for RBE3 center
 %r2=basis(r1(4:6)-r1(1:3),[1 1 1]);
 %n4=r1(1:3)'+max(n1(:,1))*r2(:,1); % center node
 [model.Node,i4]=feutil('addnodeknownnew',model.Node,n4(:)'); %add
 n5=model.Node(i4,1); %added nodeid
 if isfield(RO,'whMass')&&~isempty(RO.whMass) % add wheel mass to slave
  %mass2 input: [n1 M I11 I21 I22 I31 I32 I33]
  model=feutil('AddElt-NewId',model,'mass2',[n5 RO.whMass]);
 end
    
 % Define a RBE3 constraint
 data=struct('SlaveSel',sprintf('nodeid %i',n5),...
  'MasterSel',sprintf('nodeid %s',num2str(r3(:)')),...
  'DOF', 123456,...
  'MasterDOF', 123);
 model=fe_case(model,'rbe3',RO.name,data);
  
 % Define rotation penalty with cbush connection to clamped master node
 [model.Node,i4]=feutil('addnodeknownnew',model.Node,n4(:)');
 n6=model.Node(i4,1); % new node to cbush connection

 % add clamped master with cbush connection on rotation
 pro=feutil('proidnew',model);
 el1=[n5 n6 0 pro 0 r1(4:6)-r1(1:3) -1]; % elt
 model=feutil('addelt-newId',model,'cbush',el1);
 stu='SI'; if isfield(model,'unit')&&~isempty(model.unit); stu=model.unit; end
 pro=[pro fe_mat('p_spring',stu,2) 0 0 0 RO.krot 0 0]; % property
 model.il(end+1,1:length(pro))=pro;
 model=fe_case(model,'fixdof',[RO.name 'm'],n6); % clamp
  
  out=model;

 
elseif comstr(Cam,'squeal_five_ecl');[CAM,Cam]=comstr(CAM,5);
  %% #Squeal_Five_ECL: load open access data from the squeal bench of ECL
  % Load the "Dataset of vibrational and acoustic measurements for squeal 
  % analysis from the laboratory brake setup Friction-Induced "
  % Define the first time the directory location if the extracted archive with
  % sdtdef('Squeal_FIVE_ECL-setpref','root_directory');
  
  try;wd=sdtdef('Squeal_FIVE_ECL');
  catch;error(['You must define the directory pointing to FIVE data, '...
    'sdtdef(''Squeal_FIVE_ECL-setpref'',''root_directory'');']);
  end
  r1=load(fullfile(wd,'data.mat'));
  r1=r1.data;
  r2=load(fullfile(wd,'time.mat'));
  r2=r2.time;
  r3={'TriAxe.1.X' 'TriAxe.1.Y' 'TriAxe.1.Z' 'TriAxe.2.X' 'TriAxe.2.Y'... 
   'TriAxe.2.Z' 'TriAxe.3.X' 'TriAxe.3.Y' 'TriAxe.3.Z' 'TriAxe.4.X'...
   'TriAxe.4.Y' 'TriAxe.4.Z' 'TriAxe.5.X' 'TriAxe.5.Y' 'TriAxe.5.Z'...
   'Micro.1' 'Micro.2' 'Micro.3' 'Micro.4' 'Micro.5' 'Micro.6' 'Micro.7'...
   'Micro.8' 'Micro.9' 'Micro.10' 'Micro.11' 'Micro.12' 'Micro.13'...
   'Micro.14' 'Micro.15' 'Proxy.1' 'Proxy.2' 'Proxy.3' 'Proxy.4'...
   'RotationSpeed' 'Torque' 'BrakePressure' 'BrakeTemperature'}';
  r4=struct('X',{{r2 r3}},'Y',r1,'Xlab',{{'time','channel'}});
  ci=iiplot(2); ci.Stack{'curve','Squeal_FIVE_ECL'}=r4;
  iicom('IIxOnly',{'Squeal_FIVE_ECL'}); 
  
elseif comstr(Cam,'wd'); [CAM,Cam]=comstr(CAM,3);
 %% #Wd: resolution of working directories to locate demo files --------------
 if comstr(Cam,'abqdemo'); 
  %% #WdAbqDemo: squeal abaqus demonstration -2
  f1=sdtcheck('PatchFile',struct('fname','brake_squeal.inp', ...
      'in','demo_squeal_abq.zip','back',1));
  out=fileparts(f1);
  
 else; error('Wd %s unknown.',CAM);
 end
elseif comstr(Cam,'cb'); [CAM,Cam]=comstr(CAM,3);
 %% #Cb
if comstr(Cam,'changech')
  ci=obj; i2=setdiff([2 13 6],ci.opt(1));
  for j2=1:length(i2)
      cj=get(i2(j2),'userdata');
      if isa(cj,'sdth');iicom(cj,'ch',ci.ua.ch);end
  end
  figure(ci.opt(1))
else; error('Cb%s',CAM)
end

%% #Admin
%% #Tuto: recover model from a specific tuto step -3
elseif comstr(Cam,'tuto'); 
 eval(sdtweb('_tuto',struct('file','d_squeal','CAM',CAM)));
 if nargout==0; clear out; end
elseif comstr(Cam,'cvs')
 out='$Revision: 1.165 $  $Date: 2024/01/16 15:47:25 $';
elseif comstr(Cam,'@');out=eval(CAM);else; error('%s unknown.',CAM)
end
end

%% #SubFunc - - --------------------------------------------------------------
function RA=checkEta(model,RA); % add 4 if eta exists
%% #checkEta: add matdes 4 if a loss factor is present - - ------------------2
pl=fe_mat('getpl',model); r1=0;
for j1=1:size(pl,1)
 try; r1=feutil(sprintf('getmat %i eta',pl(j1,1)),model); catch; r1=0; end
 if r1; break; end
end
if r1==0 % also check in SE
 eli=fesuper('s_',model);
 for j1=1:size(eli,1)
  try
  SE=stack_get(model,'SE',eli{j1,1},'get');
  if isfield(SE,'Opt')&&~isempty(SE.Opt)&&any(SE.Opt(2,:)==4); r1=1; break; end
  catch; r1=0; 
  end
 end
end
if r1==0 % check mdmp loss
 pro=nl_spring('getpro type"nl_modaldmp"',model);
 if ~isempty(pro)
  for j1=1:size(pro,1)
   try
    r2=nl_spring(sprintf('getproid%i loss',pro{j1,3}.il(1)),model); 
    if r2; r1=1; break; end
   catch; r1=0;
   end
  end
 end
end

if r1; RA.matdes=[RA.matdes 4]; RA.ext='nl_solve(''modeENHcpx9 4'')'; end

end
function [uo,st]=ParTick(obj,evt)
 %% #ParTick : parameter ticks
 ua=get(obj,'userdata');cf=get(ancestor(obj,'figure'),'userdata');
 C1=cf.Stack{'Time'};t=C1.X{1}(:,1);
 i1=evt.xtick(:); st=cell(size(i1));
 for j1=1:length(i1);
  i1(j1)=sdtm.indNearest(t,i1(j1));
  val=C1.X{1}(i1(j1),evt.TickInfo.ipar);
  st{j1}=sprintf(evt.TickInfo.valFmt,val);
 end
 ga=handle(obj);delete(findobj(ga,'tag','now'))
 drawnow;r3=get(ga,'ylim')'; if ~all(isfinite(r3));set(ga,'ylimmode','auto');r3=get(ga,'ylim')';end
 set(ga,'ylim',[1.01 -.01;-.01 1.01]*r3);
 if isfield(C1,'ID')
  r1=C1.ID.po(:,1);r2=get(obj,'xlim');
  if ~all(isfinite(r2));set(ga,'xlimmode','auto');r2=get(ga,'xlim')';end
  r1=r1(r1>r2(1)&r1<r2(2));
  li=line('parent',ga,'xdata',r1,'ydata',r3(2)*ones(size(r1)),'linestyle','none','marker','.','tag','now');
 end
 evt.ytick(:)=[1.05 -.05]*r3;
 evt.ztick=([.5 .5]*get(ga,'zlim')')*ones(size(evt.ytick));
 uo=evt;

end
function out=specMax
  %% #specMax
  c13=get(13,'userdata');  ob=handle(c13.ua.ob(1));
  if isprop(ob,'ZData');r2=ob.ZData; else;r2=ob.CData;end
  if ob.YData(1)==0;r2(1:2,:)=NaN;end
  r2=max(r2,[],2); 
  st=c13.ua.YFcn; 
  if sdtm.Contains(st,'log10(abs(r3))');r2=10.^r2;end
  gf=sdth.urn('figure(101).os{@Dock,{name,SqSig},name,101 SpecMax,NumberTitle,off}');
  figure(double(gf));
  plot(ob.YData,r2);xlabel('Frequency [Hz]');ylabel('Max_t(spectro)')
  if nargout>0; out=r2;end
end
function  [r1,st]=getAmp(r1,Time,st);
  % #getAmp
  if length(Time.Xlab)~=3||~isequal(Time.Xlab{3},'hdof')
   error('Not a valid case');
  elseif strcmpi(Time.X{2}{1},'q1R');
   r2=squeeze((Time.Y(:,:,1)));st{3,1}=sprintf('%s [%s]',Time.X{2}{1,1:2});
   r1(:,2+(1:size(r2,2)))=r2;
  else
   r1(:,3)=sum(abs(Time.Y(:,:,1)).^2,2);st{3,1}={'Amplitude'};
  end
end

function wheelPosLines(c2)
  r2=c2.Stack{'Time'};
  if isfield(r2,'ID')&&(~isfield(r2.ID,'marker')||~strcmp(r2.ID.marker,'band'));
   ii_plp(r2.ID);% Add wheel position lines
   go=findobj(103,'type','line','tag','now');
   if length(go.XData)>50;set(go,'visible','off');end
  end 
end