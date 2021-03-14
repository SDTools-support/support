function [out,out1,out2]=t_beam(varargin)

% Various tests beam stress/strain validation and 
%
% Contributed by E. Balmes

%#ok<*ASGLU>

if nargin==0; Cam='';
else; [CAM,Cam]=comstr(varargin{1},1);carg=2; %#ok<*NOSEM>
end

if isempty(Cam)  % These are tests that run automatically

ofact('lu;');
t_beam('analyticStrain rb')
t_beam('analyticStrain bend')
t_beam('analyticStrain disp')
t_beam('beam1t'); % tests of variable section elements
t_beam shape
t_beam('enersub')
%t_beam('sections') generation of section views
t_beam('dbval')
%t_beam('maps') %for now tested in beam1 example
t_beam('thermalload')
ofact('spfmex;');

%% #shape  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'shape') % shape verifications
    
    model=struct('Node',[1 0 0 0  0 0 0]);
    model.Elt=feutil('ObjectMass',1,[1.1 1.1 1.1]);
    model=feutil('extrude 0 1 0 0',model,[0 .2 1]);
    model.pl=m_elastic('dbval 1 Steel');
    model.il=p_beam('dbval 1 circle .01');
    model.Elt(1,1:7)=[Inf abs('beam1t')];
    model.Elt(2:end,7)=1; % orient along z

    [model,C1,def]=t_beam('analyticStrainBend-Back',model);
    r2=integrules('beam1',linspace(0,1,5)'*[1 0 0 0]);
    st={'w','N','Nr','Nrr','xi','Nw','NDN','jdet'};
    for j1=1:length(st);
        C1.GroupInfo{end}.(st{j1})=r2.(st{j1});
    end
    C1.GroupInfo{5}.Y(1,size(r2.w,1),1,1)=0;
    data=fe_case(model,'getdata','Force');
jGroup=1;
%rk=fe_caseg(sprintf('stressObserve 1 -jgroup %i -rule -1 10',jGroup),model,C1,def);
 rk=fe_caseg('StressObserve 1',model,C1,def);
 rm=fe_caseg('StressObserve 2',model,C1,def);

 figure(1);plot(rm.Node(:,1),[rm.tx*def.def rm.ty*def.def rm.tz*def.def ...
     rm.rz*def.def rm.ry*def.def],':x')
 figure(1);plot(rk.Node(:,1),[rk.k_y*def.def rk.phi_x*def.def],':x')


elseif comstr(Cam,'ng'); 
%% #Ng : verify global shape functions formulas for beam1t

EC=integrules('beam1',[-1 40]);
model=femesh('testbeam1t divide10');
L=.1; bas=[1 0 0;0 0 1;0 -1 0];
%  L=state(10); bas=reshape(state(1:9),3,3)/L^2; 
  k=zeros(12); J=L;
  Nt=zeros(3,12,size(EC.w,1));
  for jw=1:size(EC.w,1)
    r=EC.w(jw,1);N=EC.N(jw,:).*[1 L 1 L];
    % sdtweb t_stress('cutshapebeam');
    Nt(:,:,jw)=bas'*[ ...
        [1-r 0 0;0  N(1) 0  ;   0 0       N(1)  ]*bas ... % tA(local)=bas*tA
        [0 0 0  ;0  0    N(2);0 -N(2) 0]*bas ...  % rA
        [r 0 0;  0  N(3) 0;     0 0       N(3)  ]*bas ...    % tB
        [0 0 0;  0  0    N(4);0 -N(4) 0]*bas];    % rB
   %k=k+Nl'*(dd* J*EC.w(jw,4))*Nl;
  end

% theta_y=-dw/dx  
% theta_z= dv/dx
% w=-L*theta_y for rigid body mode
  
rb=feutilb('geomrb',model,[0 0 0]);
rb=feutilb('placeindof',feutil('getdof',model.Elt(2,1:2)',(1:6)'/100),rb);
rb.def

x=squeeze(Nt(1,:,:))';%x(:,~any(x))=[];
y=squeeze(Nt(2,:,:))';%y(:,~any(y))=[];
z=squeeze(Nt(3,:,:))';%y(:,~any(y))=[];
figure(1);plot(EC.w(:,1)*L,[x y])
figure(1);subplot(211);plot(EC.w(:,1)*L,y*rb.def);ylim([0 L])
figure(1);subplot(212);plot(EC.w(:,1)*L,z*rb.def);ylim([-L 0])
setlines([],{'-','--','-.'},'+ox*sdv^><ph')

y\(-EC.w(:,1)*L);

% Report from Sami Daouk page 30, formulas N,M, and page 32 for s_xx
% Cut in a beam
 
 
%% #Hexa8 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'hexa8')
    
model=femesh('testhexa8b divide 2 2');
data=struct('sel','groupall','dir',{{'x','y'}}, ...
    'DOF',[.01;.02]);  
model=fe_case(model,'dofset','Force',data,'fixdof','base','z==0');
def=fe_simul('static',model); %fe_c(q.DOF);ans(:,2)=num2cell(q.def)

C1=fe_stress('stress -gstate',model,def);model.DOF=def.DOF;
[opt,match]=elem0('GaussObserve',C1.GroupInfo{[8 3 4]},model,C1); 

% Mesh a cut-plane
r1=struct('orig',[0 1 .1],'normal',[0 1 1]);
r1.Elt=feutil('selelt seledgeAll',model);
r1.Elt(~isfinite(r1.Elt(:,1)),:)=[];
NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
i1=reshape(full(NNode(r1.Elt(:,1:2))),[],2);
r2=(model.Node(i1(:,1),5:7)-ones(size(i1,1),1)*r1.orig)*r1.normal(:);
r3=(model.Node(i1(:,2),5:7)-ones(size(i1,1),1)*r1.orig)*r1.normal(:);
i2=sign(r2).*sign(r3)>0;
i1(i2,:)=[];r2(i2,:)=[];r3(i2,:)=[];
r4=[-r3 r2]./((r2-r3)*[1 1]);
r2=model.Node(i1(:,1),5:7).*r4(:,[1 1 1])+ ...
    model.Node(i1(:,2),5:7).*r4(:,[2 2 2]);

r3=fe_fmesh('delaunay',r2*sp_util('basis',r1.normal,[0 0 0]));
r3.Node(:,5:7)=r2;cg=feplot(5);cg.model=r3;
r3=feutil('getpatch new',r3);
r3.ScaleColorMode='one';
r3=sdsetprop(r3,'fsProp','edgealpha',0);
cg.SelF{2}=r3;cg.o(3)='sel2 ch 0 ty1 '; % undeformed

%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #AnalyticStrain Stress computations on a beam
elseif comstr(Cam,'analyticstrain') %t_beam('analyticStrain')
[CAM,Cam]=comstr(CAM,15);

if carg>nargin
    model=struct('Node',[1 0 0 0  0 0 0]);
    model.Elt=feutil('ObjectMass',1,[1.1 1.1 1.1]);
    model=feutil('extrude 0 1 0 0',model,[0 .05 linspace(.1,1,10)]);
    model.pl=m_elastic('dbval 1 Steel');
    model.il=p_beam('dbval 1 circle .01');
    model.Elt(1,1:7)=[Inf abs('beam1t')];
    model.Elt(2:end,7)=1; % orient along z
else; model=varargin{carg};carg=carg+1;
end

%model.pl=[100 fe_mat('m_elastic','SI',5) ...
%      38.095e9  9.4e9 0.3   7.5e9 0    0     1640];
[CAM,Cam,RunOpt.Back]=comstr('-back',[-25 3],CAM,Cam);
data2=[];d2=[];
% analytic evaluation
if carg<=nargin&&iscell(varargin{carg})
  % All standard Defs
  list=varargin{carg};carg=carg+1;def=[]; %#ok<*NASGU>
  for j1=1:length(list)
   [mo1,C1,d1]=t_beam(sprintf('analyticStrain%s -back',list{j1}),model); 
   def=fe_def('appenddef',def,d1);
  end
  out=def;
  return;
elseif isempty(Cam)||comstr(Cam,'rb') % rigid body rotation 
 data=struct('sel','groupall','dir',{{2,'-2*x',1,'x',}}, ...
    'DOF',[.05;.03;.06;.02]);RunOpt.Type='rb';
elseif comstr(Cam,'bend')||comstr(Cam,'disp') 
 % Euler bernoully bending : theta = +- dw/dx
 data=struct('sel','groupall','dir',{{'x','-x.^2/2','2*x','x.^2'}}, ...
    'DOF',[.05;.03;.06;.02]);RunOpt.Type='bend';  
 % Constant shear
 data2=struct('sel','groupall','dir',{{'x.^2/2','-x.^3/6','x.^2','x.^3/3'}}, ...
    'DOF',[.05;.03;.06;.02]);RunOpt.Type='bend';  
end
model=fe_case(model,'dofset','Force',data);
model.name='';[i1,model.Elt]=feutil('eltidfix;',model);
def=fe_simul('static',model); %fe_c(q.DOF);ans(:,2)=num2cell(q.def)
if size(def.def,2)==1;def.lab={RunOpt.Type};end
if ~isempty(data2)
 model=fe_case(model,'dofset','Force',data2);
 d2=fe_simul('static',model);
end

[C1,model.DOF]=fe_mknl('init',model); 
def=feutilb('placeindof',model.DOF,def);
% eval(iigui({'def','model','C1'},'SetInBaseC'))
C1=fe_mknl('gstate-struct',model,C1,def);
k=fe_mknl('assemble',model,C1,fe_def('subdef',def,1),1);
%out1=fe_stress('stressgstate',model,def)
r1=C1.GroupInfo{5};
if RunOpt.Back; out=model;out1=C1;out2=def; return;end

dd=feval(beam1t('@buildDD'),C1.GroupInfo{4},0);EI=dd(5);EA=dd(1);
%constit=C1.GroupInfo{4};EI=constit(1)*constit(9);
%EA=model.pl(3)*model.il(6);

switch RunOpt.Type
case 'rb'
  if norm(r1.Y(:),'inf')>1e-4; error('RB should have no strain');end
case 'bend' % and axial strain
 i2=strcmpi(r1.X{1},'k_y');i3=strcmpi(r1.X{1},'k_z');
 r1.Y(i2,:,:)=r1.Y(i2,:,:)+EI;
 r1.Y(i3,:,:)=r1.Y(i3,:,:)+2*EI;
 if norm(r1.Y(:),'inf')>1e-4; error('Inconsistent bend');end
 % Now the same thing with StressCut
 mo1=fe_caseg('StressCut',struct('type','BeamGauss'),model);
 mo1.type='ObserveAtGauss';mo1.Rule=-3;
 cut=fe_caseg('StressCut-SelOut',mo1,model);
 cut.StressObs.CritFcn='';
 R1=fe_caseg('stressobserve ',cut,def);
 R1.Y(i2,:)=R1.Y(i2,:)+EI;
 R1.Y(i3,:,:)=R1.Y(i3,:,:)+2*EI;
 if norm(R1.Y(:),'inf')>1e-4; error('Inconsistent bend');end
 % verify the axial stress
 model=fe_case(model,'reset','fixdof','base','x==0','DofLoad','Tip',12.01);
 d1=fe_simul('static',model);
 R1=fe_caseg('stressobserve ',cut,d1);
 if norm(R1.Y(1,:)-1)>1e-10; error('problem with axial strain');end
 
 R2=fe_caseg('stressobserve ',cut,d2);
 i4=strcmpi(r1.X{1},'s_y');i5=strcmpi(r1.X{1},'s_z');
 R2.Y(i4,:)=R2.Y(i4,:)+dd(2); R2.Y(i5,:,:)=R2.Y(i5,:,:)+2*dd(2);
 if norm(R2.Y(i4|i5,:))>1e-4;error('on shear');end
 % verify torsion
 data2=struct('sel','groupall','dir',{{'x'}},'DOF',.04);
 model=fe_case(model,'reset','dofset','Force',data2);
 d2=fe_simul('static',model);
 R2=fe_caseg('stressobserve ',cut,d2);
 if norm(R2.Y(strcmpi(R2.X{1},'phi_x'),:)-dd(4))>1e-4;error('on torsion');end

 
case 'disp'
% Check the location of Gauss points - - - - - - - - - - - - - - - - - - -
r1=basis('rotate',[],'rz=10',1);r1=reshape(r1(7:15),3,3);
model.Node(:,5:7)=model.Node(:,5:7)*r1;
def=fe_simul('static',model); %fe_c(q.DOF);ans(:,2)=num2cell(q.def)

[C1,model.DOF]=fe_mknl('init',model); 
def=feutilb('placeindof',model.DOF,def);
C1=fe_mknl('gstate-struct',model,C1,def);
k=fe_mknl('assemble',model,C1,fe_def('subdef',def,1),1);

[opt,match]=elem0('GaussObserve',C1.GroupInfo{[8 3 4]},model,C1);
mo2=struct('Node',[(1:size(match.Node,1))'*[1 0 0 0],match.Node], ...
    'Elt',[]);
mo2.Elt=feutil('objectbeamline',1:size(mo2.Node,1));
mo2=feutil('addtest',model,mo2);
i2=stack_get(mo2,'info','OrigNumbering','getdata');
i2=double(i2(:,2)');i2=reshape([i2+.01;i2+.02;i2+.03],[],1);
d2=struct('def',[def.def;opt.CIN'*def.def],'DOF',[def.DOF;i2]);
feplot(mo2,d2)

otherwise

 rule=integrules('beam1',3);   
 jw=2;x=model.Node(1:end-1,5)*(1-rule.w(jw))+model.Node(2:end,5)*rule.w(jw);
 r2=squeeze(r1.Y(:,jw,:));figure(1);plot(x,r2,'x-')
       
end

if 1==2  % Plot shape functions in detail - - - - - -  - - - - 
r3=integrules('beam1',linspace(0,1,30)'*[1 0 0 0]);
figure(1);subplot(211);
x=r3.w(:,1);plot(x,[r3.N*[0 0 1/2 1]',x.^2/2])

    
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% #StressOutput/cut
elseif comstr(Cam,'stress'); [CAM,Cam]=comstr(CAM,9);

%% -> d_mesh ?
l=1;
mdl=struct('Node',[1 0 0 0  0  0 0;   
                    2 0 0 0  l  0 0]);
prop=[100 100 0 0 1.1]; % MatId ProId nx ny nz

mdl.Elt=feutil('ObjectBeamLine 1 2',prop);
mdl=feutil('Divide 100',mdl);
mdl.Elt=feutil('SetGroup 1 name beam1t',mdl);
mdl.Elt=feutil('Divide group 1 withnode{x<0.3}',mdl);
mdl.Elt=feutil('Divide group 2 withnode{x>0.7}',mdl);
mdl.Elt=feutil('Setgroup3 pro101',mdl);
%feutil('info',mm)
% Properties :
mdl.pl=m_elastic('dbval 100 steel'); % mat
mdl.il=[100 fe_mat('p_beam','SI',3) 0 comstr('I',-32) ...
                 0.1 0.05 0.05 0.01 0.015 0.015];
mdl.il=p_beam(mdl.il,'dbval 101 circle .05');
mdl=feutil('setpro 101',mdl,'StressOut', ...
   struct('Node',[1 0 0 0  0 0 0; 2 0 0 0  0 0.03 0],...
                       'lab',{{'Center';'vertpoint3cm'}}));
mdl=feutil('setpro 100',mdl,'StressOut','default');



%r1=p_beam('StressObserve',mdl,[100 100]); % low level call to observe

C1=fe_mknl('init',mdl);C1.GroupInfo{1,end}
cutg=fe_caseg('StressCut',struct('type','BeamGauss'),mdl);
cut=fe_caseg('stresscut -radius 10 -SelOut',cutg,mdl);

if 1==2 % Eb need to see why match is needed
  cut=fe_caseg('stresscut -radius 10 -SelOut',struct('type','gauss'),mdl);
end

% xxxJP : add a pure tension and pure torsion to analytic strain

d1=t_beam('analyticStrain',mdl,{'Rb','Bend'});
rb=feutil('geomrb',mdl);rb.lab={'x';'y';'z';'rx';'ry';'rz'};
d1=fe_def('appenddef',d1,feutilb('placeindof',d1.DOF,rb));
C1D=fe_caseg('stressobserve -crit""',cut,d1); % Observation as CURVE

squeeze(C1D.Y(:,1,:))

C1D % XXXeb mismatch between X{2} and dim2 of Y , only 1st group taken in

    
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% #3dstrain Droppers tension computation
elseif comstr(Cam,'3dstrain'); [CAM,Cam]=comstr(CAM,9);

error('moved to fe_caseg(''stressobserve'')');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% #StaticUpdatedLag
elseif comstr(Cam,'staticupdatedlag'); [CAM,Cam]=comstr(CAM,9);

model=beam1t('test');model.name='';  
model=fe_case(model,'FixDof','2D',[.02 .04 .06]);
model=feutil('setpro 112',model,'MAP', ...
 struct('dir',{{'1e4'}},'lab',{{'ten'}},'Formulation','of_mk','data0','Need'));
   
% check assembly and state updates
[Case,model.DOF]=fe_mknl('init',model);
m=fe_mknl('assembleNoT',model,Case,2);
g=-(m*sum(fe_c(model.DOF,.03))');
model=fe_case(model,'DofLoad','grav',struct('def',g,'DOF',model.DOF), ...
    'DofLoad','MidPoint',struct('def',500,'DOF',29.01));
q0=fe_simul('static',model);
feplot(model,q0);
r1=beam1t('viewTen',model,q0);
fe_mknl('assemble',model,Case,302,sum(q0.def,2)); %302 update rotation
IA=Case.GroupInfo{7};figure(1);plot(IA.data(11,:))

%fe_stress('stress-gstate',model,q0) : cannot work right now due to straindef

dc=q0;dc.def(:,2)=0;
k=fe_mknl('assemble',model,Case,dc,5);sparse(dc.def(:,2))

%% Test the updated lagrangian procedure
opt=struct('Opt',[],'Method','StaticNewton',...
   'Jacobian','ki=basic_jacobian(model,ki,0.,0.,opt.Opt);',...
   'NoT',1, ... % Don't eliminate constraints in model.K
   'AssembleCall','assemble -fetimeNoT -cfield1', ...
   'IterInit','opt=fe_simul(''IterInitLagUpdate'',model,Case,opt);');
model=fe_case(model,'DofLoad','grav',struct('def',g,'DOF',model.DOF), ...
    'DofLoad','MidPoint',struct('def',1,'DOF',29.03));
nstep=20;cf=feplot(model); ofact('lu');
%C1=fe_curve(sprintf('testramp NStep=%i t0=1 Yf=2000',nstep));C1.Y=C1.Y-C1.Y(1);
C1=struct('X',(0:nstep-1)','Y',2*ones(nstep,1));
C2=struct('X',(0:nstep-1)','Y',ones(nstep,1));C2.Y(2:end)=0;
model=fe_case(model,'setcurve','MidPoint',C1); % 20 steps gradual load
model=fe_case(model,'setcurve','grav',C2);
q0=fe_time(opt,model);
cf.def=q0;fecom ch1:20
n1=sortrows(model.Node,5);r2=fe_c(q0.DOF,n1(:,1)+.03)*q0.def;
figure(1);subplot(211);plot(r2)
r3=diff(fe_c(q0.DOF,29.03)*q0.def);subplot(212);plot(r3);axis tight

%% #Beam1t_compile stiffness

model=beam1t('test divide 1000');model.name='';  
model=fe_case(model,'FixDof','2D',[.02 .04 .06]);
model=feutil('setpro 112',model,'MAP', ...
 struct('dir',{{'1e4'}},'lab',{{'ten'}},'Formulation','of_mk','data0','Need'));
   
% check assembly and state updates
[Case,model.DOF]=fe_mknl('init',model);
k=fe_mknl('assembleNoT',model,Case,1);


%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #Beam1t Stress computations prestressed beam example 
elseif comstr(Cam,'beam1t') ;[CAM,Cam]=comstr(CAM,15);

model=femesh('testbeam1');model.name='';
model.Elt=feutil('set group 1 name beam1t',model);
model=fe_case(model,'fixdof','clamp',[1]);
[Case,model.DOF]=fe_mknl('init',model); %Case.GroupInfo{1,7}.data(11,:) 
model=fe_case(model,'dofLoad','tens',struct('DOF',2.01,'def',1e6));
def=fe_simul('Static',model);
fe_mknl('assemble',model,Case,300,def.def); % update state
if Case.GroupInfo{1,7}.data(11,:)~=1e6;error('Mismatch');end

k1=fe_mknl('assemble',model,Case,1); 

% Element computations for of_mk -> switch with CTable in variable section
model=feutil('setpro 112',model,'MAP', ...
   struct('dir',{{'1e6'}},'lab',{{'ten'}},'Formulation','of_mk'));
[C1,model.DOF]=fe_mknl('init',model); %Case.GroupInfo{1,7}.data(11,:) 
C1.GroupInfo{3}(3:4)=int32([12;6]);
k=fe_mknl('assemble',model,C1,1); % co
if normest(k-k1); error('mismatch');end


%% Variable section beam using interpolated constitutive fields

model=femesh('testbeam1 divide 30');model.il=[];
model.Elt=feutil('set group 1 name beam1t',model);
model=fe_case(model,'fixdof','clamp',[1]);

%xxx matgui('get_plstack') 
%matgui('get_stackpl')

pro=p_beam('database 112 circle .1');
pro.il(3:6)=-2; % -2 implies local interp, -1 interp during GetPl
pro.A=struct('X',[0;1],'Xlab',{{'x'}},'Y',.01*[1;2]);
pro.I1=struct('X',[0;1],'Xlab',{{'x'}},'Y',[1;2]);
pro.I2=struct('X',[0;1],'Xlab',{{'x'}},'Y',2*[1;2]);
pro.J=struct('X',[0;1],'Xlab',{{'x'}},'Y',3*[1;2]);

model=stack_set(model,'pro','CW',pro);
[Case,model.DOF]=fe_mknl('init',model);EC=Case.GroupInfo{1,8};
constit=Case.GroupInfo{1,4};
InfoAtNode=Case.GroupInfo{7};

NodePos=int32([1 4;1 2])';
for jElt=1:size(NodePos,2)
 [nodeE,nodeEt]=feval(elem0('@get_nodeE'),model.Node,NodePos(:,jElt),1,InfoAtNode,[],EC);
 ConstitInterp=elem0('@ConstitInterp');
 r1=feval(ConstitInterp,nodeE,[.5;.5],constit,EC.CTable);
 fprintf('J %.4g I1 %.4g I2 %.4g A %.4g\n',r1(8:11))
 of_time('cinterp',nodeE,[.5;.5],constit,EC.CTable,zeros(1));
 if any(constit-r1); error('Mismatch');end
end
%sdtweb of_mk.c#ctable
%sdtweb of_mk_pre.c#ctable
[Case,model.DOF]=fe_mknl('init',model);
k=fe_mknl('assembleNoT',model,Case,1);
%dk=diag(k);figure(1);plot(dk(1:6:end))


%% Variable section beam using MAP
model=femesh('testbeam1 divide 30');model.il=[];
model.Elt=feutil('set group 1 name beam1t',model);
model=fe_case(model,'fixdof','clamp',[1]);

pro=p_beam('database 112 circle .1');
model=stack_set(model,'pro','CW',pro);
[Case,model.DOF]=fe_mknl('init',model);
InfoAtNode=Case.GroupInfo{7};
pro.il(3:6)=-3; model=stack_set(model,'pro','CW',pro);% UsenodeE

InfoAtNode.data=(1:4)'*linspace(1,2,size(InfoAtNode.data,2));
InfoAtNode.lab={'J','I1','I2','A'};
pro.MAP=InfoAtNode;
model=stack_set(model,'pro','CW',pro);
[Case,model.DOF]=fe_mknl('init',model);InfoAtNode=Case.GroupInfo{7};
EC=Case.GroupInfo{1,8};
sdtkey('cvsnum>=1.036;','of_time') % Field interp needed
    
NodePos=int32([1 4;1 2])';
for jElt=1:size(NodePos,2)
 [nodeE,nodeEt]=feval(elem0('@get_nodeE'),model.Node,NodePos(:,jElt),1,InfoAtNode,[],EC);
 ConstitInterp=elem0('@ConstitInterp');constit=Case.GroupInfo{1,4};
 r1=feval(ConstitInterp,nodeE,[.5;.5],constit,EC.CTable);
 fprintf('J %.4g I1 %.4g I2 %.4g A %.4g\n',r1(8:11))
 of_time('cinterp',nodeE,[.5;.5],constit,EC.CTable,zeros(1));
 if any(constit-r1); error('Mismatch');end
end

% Need update
%disp({'p_solid','beam1','beam1t','elem0','matgui','of_time'});


% ------------------------------------------------------------
%% #EnerSub Energy fractions in pre-stressed beams
elseif comstr(Cam,'enersub') %t_beam('maps')

model=beam1t('test');
  model.il(1,4)=1e-5;
  [Case,model.DOF]=fe_mknl('init',model);
  T=1e8; Case.GroupInfo{1,7}.data(11,:)=T;
  
  k=fe_mknl('assemble',model,Case,1);
  m=fe_mknl('assemble',model,Case,2);

  g=(m*sum(fe_c(Case.DOF,.03))');
  kd=ofact(k); q=kd\g;ofact('clear',kd)
 def=struct('def',Case.T*q,'DOF',model.DOF);
eltid=feutil('eltidfix',model);
RO=struct('MatDes',1);
E=fe_caseg('enersub',model,Case,def,RO);
r1=fe_mknl('assemble',model,Case,251,def);
i2=find(r1);
if norm(r1(i2)./sum(E.Y(i2,:,:),3)-1)>1e-10
    error('Energy fraction mismatch');
end
data=fe_stress('stress-gstate',model,def);data=data.GroupInfo{5};
if ~isfield(data,'Y');error('stress computation failed');end
%squeeze(E.Y)

% ------------------------------------------------------------
%% #Sections Beam sections capture
elseif comstr(Cam,'sections') %t_beam('sections')
 
 wd=comgui('cd','o:\sdt.cur\tex\plots\');
 mm=p_beam('database');
 cf=feplot;
 for j1=1:length(mm)
  [r1,sectionmdl]=p_beam('ConvertTo1',mm(j1));
  if ~isempty(sectionmdl)  
   cf=feplot(sectionmdl);
   fecom colordatamat
   comgui('imwrite',fullfile(wd,[mm(j1).name '_section.png']))
  end
 end

% ------------------------------------------------------------
%% #dbval Database and subtype conversion
elseif comstr(Cam,'dbval') %t_beam('sections')
 
 pro=p_beam('database 100 rectangle .05 .01');
 il=p_beam('dbval 101 circle .05');
 if ~isequal(pro.il(6),0.05*0.01)||abs(il(6)-pi*0.05^2)>1e-5
  sdtw('_err','dbval or database lost input parameters')
 end
 il=p_beam(sprintf('dbval 1 BAR %.15g %.15g',.18,.3));
 [st,i1,i2]=fe_mat('type',il(2));
 if i2~=3; sdtw('_err','something change : for nastran section database used to return subtype 3'); end
 il2=p_beam('ConvertTo1',il);
 if abs(il2(6)-.18*.3)/.18*.3>1e-8
  sdtw('_err','dbval or database lost input parameters or error in section area computation')
 end
 
 % assembling test 
 mdl=femesh('testbeam1');mdl.name='';
 mdl.il(4,:)=[];
 r1=p_beam(sprintf('database 112 ROD %.15g',1)); 
 mdl.il(4,1:length(r1.il))=r1.il;
 mdl=fe_case('assemble-secdof-reset-matdes 2 3 1',mdl);
 %def=fe_eig(mdl,[5 1e3 20]);
 mdl.il(4,1:8)=p_beam('convertto1',r1.il); mdl.il(4,1)=112;
 mdl2=fe_case('assemble-secdof-reset-matdes 2 3 1',mdl);
 %def2=fe_eig(mdl,[5 1e3 20]);
 if norm(mdl2.K{3}-mdl.K{3},'inf')/norm(mdl.K{3},'inf')>1e-6
  sdtw('_err','nastran cross section beam assembling test failed')
 end
 
%% #ThermalLoad ------------------------------------------------------------
elseif comstr(Cam,'thermalload') %t_beam('thermalload')

model=femesh('testbeam1 divide 10');model.Elt(1,1:7)=[Inf abs('beam1t')];
model.Node(:,6:7)=model.Node(:,5)*[.1 .2];
model.pl=[100 fe_mat('m_elastic','SI',1) 210e9 .3 7800 0 0 1e-5 20];
model.il=feutil('getil',model);
 
% Uniform temperature of 30 degrees (nominal 20)
defT=struct('def',ones(size(model.Node,1),1)*30,'DOF',model.Node(:,1)+.20);
model=fe_case(model,'DofSet','ThermalState',defT);
[Case,model.DOF]=fe_mknl('init',model);

b=struct('def',fe_mknl('assemble NoT',model,Case,103),'DOF',model.DOF);

model=fe_case(model,'FixDof','base',1,'DofLoad','Thermal',b);

def=fe_simul('static',model);
fe_mknl('assemble',model,Case,300,def); % computation of tensions
r1=mean(Case.GroupInfo{7}.data(11,:))- ...
    model.pl(3)*model.il(6)*model.pl(8)*10; % E A alpha dt

if abs(r1)>1e-5;error('mismatch');end

%feplot(model,def);fecom showdefarrow



% ------------------------------------------------------------
%% #Maps Displays of orientation maps
elseif comstr(Cam,'maps') %t_beam('maps')
    
  cf=feplot(femesh('test2bay'));
  % Map is in very variable direction due to undefined nr
  % This is only ok for sections invariant by rotation
  beam1t('map',cf.mdl);fecom('view3'); 
  
  % Now define generator for bending plane 1
  i1=feutil('findelt eltname beam1',cf.mdl);
  cf.mdl.Elt(i1,5:7)=ones(size(i1))*[-.1 .9 0]; % vx vy vz
  beam1t('map',cf.mdl);fecom('view2'); 
%% #end
else; error('%s',CAM);
end
