function t_plate2(varargin)

% Various tests for composite shell development.
%
% Contributed by B. Desmorat and E. Balmes

%	Etienne Balmes
%   Copyright (c) 1990-2012 by SDTools, All Rights Reserved.
%   For revision see fe_norm('cvs')

if nargin==0; Cam='';
else; [CAM,Cam]=comstr(varargin{1},1);carg=2;
end
%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC>

if isempty(Cam)

t_plate2('patch4')
t_plate2('patch9')
model=evalin('base','model');%cf=feplot(model); %#ok<NODEF>
t_plate2('analyticStrain x',model)
t_plate2('analyticStrain xy',model)
t_plate2('analyticStrain s',model)
t_plate2('analyticStrain r',model)

t_plate2('jacobian')
t_plate2('cylinder')
t_plate2('sigma')
t_plate2('maps')
% sdtweb t_plate('ortho') % testing of composite

%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #Patch Computation with q4cs (irregular mesh)
elseif comstr(Cam,'patch')

model=ofdemos(CAM);
MAP=feutil('getnormalElt MAP -dir1',model);
MAP.normal(:,1)=0;MAP.normal(:,2)=1;MAP.normal(:,3)=0;
model=p_shell('setTheta -strategy 2',model,MAP);
MAP1=feutil('getnormalElt MAP -dir1',model);  % incorrect NASTRAN strategy
if norm(MAP1.normal-MAP.normal)>1e-10; 
    warning('feutil should be fixed, local basis mismatch');
end

% Get the orientation MAP that is truly used
MAP2=feutilb('shellmapnodepos',model.Node,model.Elt);
MAP2=struct('normal',MAP2.data(1:3,MAP2.NodePos(1,:))', ...
     'ID',MAP.ID,'vertex',MAP.vertex,'opt',1);

%if norm(MAP2.normal-MAP.normal)>1e-10; error('local basis mismatch');end;

Phi=0; H=0.005;E=200e+9; nu=0.3; E1=E; E2=E; G12=E/(2*(1+nu)); nu12=nu;
%E1  = 180e+9; E2  = 10.3e+9; G12 = 7.17e+9; nu12 = 0.28;
G13=G12; G23=G12;

model.pl=[100 fe_mat('m_elastic','SI',5) E1 E2 nu12 G12 G13 G23 0];
model.il=[110 fe_mat('p_shell',1,2) -H/2 0 0 0 0 0 0 100 H Phi 0];
r1=feutil('getnode groupall',model);
model = fe_case(model,'FixDof','CL encastrement ', ...
    sprintf('x==%.15g',min(r1(:,5))));
data =  struct('sel','groupall','dir',[1e+8 0 0]); 
data =  struct('sel','groupall','dir',{{1e+8, 0,0}},'DOF',[.01;.02;.03]); 

model =fe_case(model,'AddToCase 1','Fsurf','Charge Fsurf',data); 
def=fe_simul('static',model);

% mo2 : model sans prise en compte de setTheta
MAP.normal(:,1)=1;MAP.normal(:,2)=1;MAP.normal(:,3)=0;
mo2=p_shell('setTheta -strategy 2',model,MAP);
def2=fe_simul('static',mo2);
i1=fe_c(def.DOF,.06,'ind',2);
if norm(def.def(i1)-def2.def(i1))>1e-10;
  r1=[def.def def2.def];i1=find(any(r1,2));i1( abs(r1(i1,1)./r1(i1,2)-1)<1e-5)=[];
  fe_c(def.DOF(i1))
  cf=feplot(mo2,def);cf.def(2)=def2;fecom('show2def');
    error('displacements mismatch');
end

if 1==2
mo2=model;mo2.Elt=model.Elt(:,1:6);
def2=fe_simul('static',mo2);
if  norm(def.def(i1)-def2.def(i1))>1e-10
    error('displacements mismatch');
end
end

i1=fe_c(def.DOF,[3.01;4.01]);
if size(model.Elt,1)==10 && ...% 9 element patch with symmetric element
 abs(fe_c(def.DOF,[3.01;4.01],[1 -1])*def.def)>1e-10 % comparaison de ux pour les noeuds 3 et 4
 feplot(mo2,def);
 error('symmetry mismatch');
end
assignin('base','model',model);

%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #Cylinder Computation with q4cs 
elseif comstr(Cam,'cylinder')

model=feutil('object cylinder  0  0  0  0  0  1   0.2  10    20');
FEnode = model.Node;FEel0 = model.Elt;
femesh('set groupa1 name quad4 matid100 proid110');femesh('addsel');
model=femesh('model');  
% base locale initiale 
MAP1=feutil('getnormalElt MAP -dir1',model); % Xloc = etheta
MAP2=feutil('getnormalElt MAP -dir2',model); % Yloc = z
MAP3=feutil('getnormalElt MAP',model);       % Zloc = er
% nouvelle base locale par projection de Xloc sur -z
MAP=MAP1;MAP.normal(:,1)=0;MAP.normal(:,2)=0;MAP.normal(:,3)=-1;
model=p_shell('setTheta',model,MAP);

MAP1new=feutil('getnormalElt MAP -dir1',model); % Xlocnew = -z
MAP2new=feutil('getnormalElt MAP -dir2',model); % Ylocnew = etheta
MAP3new=feutil('getnormalElt MAP',model);       % Zlocnew = er
% verification
if norm(MAP1new.normal(1,:)-[0 0 -1])>1e-10; error('local basis mismatch');end;
if norm(MAP2new.normal(1,:)-MAP1.normal(1,:))>1e-10; error('local basis mismatch');end; 

%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #Cyl_irreg Computation with q4cs (cylinder with irregular mesh)
elseif comstr(Cam,'cyl_irreg')
    
model=ofdemos('patch9');
model.Node(:,5)=model.Node(:,5)-0.5;
model.Node(:,7)=sqrt(0.5^2-model.Node(:,5).^2);
% 
MAP=feutil('getnormalElt MAP -dir1',model); 
MAP.normal(:,1)=0;MAP.normal(:,2)=1;MAP.normal(:,3)=0;
model=p_shell('setTheta -strategy 2',model,MAP);
% Get the orientation MAP that is truly used
MAP2=feutilb('shellmapnodepos',model.Node,model.Elt);
MAP2=struct('normal',MAP2.data(1:3,MAP2.NodePos(1,:))', ...
     'ID',MAP.ID,'vertex',MAP.vertex,'opt',1);

 
 % Verifications des cartes imposees
 data=struct('dir',{{0,1,0}},'lab',{{'v1x','v1y','v1z'}});
 model=feutil('setpro 110',model,'MAP',data);
 Case=fe_mknl('init',model);
 MAP=Case.GroupInfo{1,7};
 
 
%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #Sigma Computation with q4cs (orthotropic and setTheta)
elseif comstr(Cam,'sigma')

femesh(';reset;testquad4'); % local basis (X,Y,Z)
femesh('divide 5 5');mo1=femesh('model0'); % mo1 local basis: (X  Y Z)
mo1.Node(:,6)=mo1.Node(:,6)*2;

E1  = 180e+9; E2  = 10.3e+9; G12 = 7.17e+9; nu12 = 0.28;G13=G12; G23=G12;

mo1.pl=[100 fe_mat('m_elastic','SI',5) E1 E2 nu12 G12 G13 G23 0];
mo1.il = [110   fe_mat('p_shell',1,2) -0.0025 0 0 0 0 0 0 100 0.005 0 0];

MAP=feutil('getnormalElt MAP -dir1',mo1);
MAP.normal(:,1)=1;MAP.normal(:,2)=0;MAP.normal(:,3)=0;
mo1=p_shell('setTheta',mo1,MAP);
mat_orient1=feutil('getnormalElt MAP -dir1',mo1);
mo1.Elt(2:end,8)=1e-10; % force non-zero

mo2=mo1;mo2.Elt(2:end,1:4)=mo2.Elt(2:end,[2 3 4 1]); % mo2 local basis: (Y -X Z)
MAP.normal(:,1)=-1;
mo2=p_shell('setTheta',mo2,MAP);
mat_orient2=feutil('getnormalElt MAP -dir1',mo2);
%mo2.il = [110   fe_mat('p_shell',1,2) -0.0025 0 0 0 0 0 0 100 0.005 90 0];
if norm(mat_orient1.normal+mat_orient2.normal)>1e-10; 
    error('Material orientation mismatch')
end

%C1=fe_mknl('init',mo1); C1.GroupInfo{7}.data(:,1)
%C2=fe_mknl('init',mo2); C2.GroupInfo{7}.data(:,1)

sig=cell(1,4); 
for j1=1:2
 model=eval(sprintf('mo%i',j1));
 model = fe_case(model,'FixDof','CL encastrement ','x==0 | x==1');
 data =  struct('sel','groupall','dir',[0 -1e+7 -1e+4]); 
 model =fe_case(model,'AddToCase 1','Fsurf','Charge Fsurf',data); 

if ~sp_util('issdt') % openfem
  Load  = fe_load(model);
  [Case,model.DOF]=fe_mknl('init',model);
  k=fe_mknl('assemble',model,Case,1);q=ofact(k,Load.def);
  def = struct('def',Case.T*q,'DOF',model.DOF);
else
  def=fe_simul('static',model);
end

% current call
r1=fe_stress('stress gstate-curve',model,def);sig{j1*2-1}=r1{3};
% New call through q4cs (account for orientation)
[Case,model.DOF]=fe_mknl('init',model);Case.GroupInfo{end}.defe=zeros(24,1);
Case=fe_mknl('gstate-struct',model,Case); % fill gstate structure
Case.GroupInfo{8}.defe=zeros(24,1);
Case.GroupInfo{7}.data
k=fe_mknl('assemble',model,Case,def,1);sig{j1*2}=Case.GroupInfo{5};
eval(sprintf('d%i=def;',j1)); % build d1 and d2 
eval(sprintf('C%i=Case;',j1)); % build d1 and d2 
end % Loop on cases

% Now compare results, displacements
i1=fe_c(d1.DOF,.06,'ind',2);i1=setdiff(i1,find(~any(d1.def,2)));
i1(abs([d1.def(i1)./d2.def(i1)-1])<1e-8)=[];
[d1.def(i1)./d2.def(i1)-1]

if norm(d1.def(i1,:)-d2.def(i1,:))>1e-10; 
 cf=feplot(mo1,d1);cf.def(2)=d2;fecom show2def
 error('displacements mismatch');
end;
% stress
r1=sig{2}.Y(:,[2 3 4 1]); r2=sig{4}.Y(:,1:4);
r2(7:8,:)=-r2(7:8,:); % shear has rotated : -> sign change
test=r1-r2; test(9,:)=0; % drilling shear is not expected to be invariant

if  norm(test,'inf')>2e-8*norm(r1,'inf'); error('stress Mismatch');end

%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #Q4CS_stress Sample stress computation with q4cs

femesh(';reset;testquad4'); % local basis (X,Y,Z)
%femesh('divide 30 30');mo1=femesh('model0');
femesh('divide 5 5');mo1=femesh('model0'); % mo1 local basis: (X  Y Z)

mo2=mo1;mo2.Elt(2:end,1:4)=mo2.Elt(2:end,[2 3 4 1]); % mo2 local basis: (Y -X Z)


E=200e+9; nu=0.3; E1=E; E2=E; G12=E/(2*(1+nu)); G13=G12; G23=G12; nu12=nu;
mo1.pl=[100 fe_mat('m_elastic','SI',5) E1 E2 nu12 G12 G13 G23 0];
mo1.il = [110   fe_mat('p_shell',1,2) -0.0025 0 0 0 0 0 0 100 0.005 0 0];
mo2.pl=mo1.pl; mo2.il=mo1.il;

sig=cell(1,4);
for j1=1:2
 model=eval(sprintf('mo%i',j1));
 model = fe_case(model,'FixDof','CL encastrement ','x==0 | x==1');
 data =  struct('sel','groupall','dir',[0 -1e+7 -1e+4]); 
 model =fe_case(model,'AddToCase 1','Fsurf','Charge Fsurf',data); 

 if ~sp_util('issdt') % openfem
  Load  = fe_load(model);
  [Case,model.DOF]=fe_mknl('init',model);
  k=fe_mknl('assemble',model,Case,1);q=ofact(k,Load.def);
  def = struct('def',Case.T*q,'DOF',model.DOF);
 else
  def=fe_simul('static',model);
 end

% current call : now implemented sdtweb elem0('stress_og')
r1=fe_stress('stress gstate-curve',model,def);sig{j1*2-1}=r1{3};
% New call through q4cs (account for orientation)
[Case,model.DOF]=fe_mknl('init',model);Case.GroupInfo{end}.defe=zeros(24,1);
Case=fe_mknl('gstate-struct',model,Case); % fill gstate structure
Case.GroupInfo{8}.defe=zeros(24,1);
k=fe_mknl('assemble',model,Case,def,1);sig{j1*2}=Case.GroupInfo{5};
sig{j1*2}.name='StressGstate';
eval(sprintf('d%i=def;',j1));
end % loop 1:2

% check multi-vector
d2=def;d2.def(:,2)=d2.def(:,1)*2;
r1=fe_stress('stress gstate-curve',model,d2);r1=r1{3};
if norm(sum(reshape(r1.Y(:,:,:,2),9,[])) ./ sum(reshape(r1.Y(:,:,:,1),9,[]))-2)>1e-10
    error('Mismatch on multi-stress computation');
end

z=mo1;z;Elt=[mo1.Elt(1:2,:);mo2.Elt(2,:)];

% Now compare results
% membrane
if norm(sig{1}.Y(1:3,1:4)./sig{2}.Y(1:3,1:4)-1)>1e-10
  error('Mismatch');
end
  % fe_stress call
[sig{2}.Y(1:3,1:4);sig{4}.Y(1:3,1:4)];% q4cs call

% bending
[sig{2}.Y(4:6,1:4);sig{4}.Y(4:6,1:4)];% q4cs call

% bending
[sig{2}.Y(7:8,1:4);sig{4}.Y(7:8,1:4)];% q4cs call

% mo1 local basis: (X  Y Z) ; mo2 local basis: (Y -X Z)
test=[(sig{2}.Y(1,[2 3 4 1])-sig{4}.Y(2,1:4))... % Nxx1 =  Nyy2
      (sig{2}.Y(2,[2 3 4 1])-sig{4}.Y(1,1:4))... % Nyy1 =  Nxx2
      (sig{2}.Y(3,[2 3 4 1])+sig{4}.Y(3,1:4))... % Nxy1 = -Nxy2
      (sig{2}.Y(4,[2 3 4 1])-sig{4}.Y(5,1:4))... % Mxx1 =  Myy2
      (sig{2}.Y(5,[2 3 4 1])-sig{4}.Y(4,1:4))... % Myy1 =  Mxx2
      (sig{2}.Y(6,[2 3 4 1])+sig{4}.Y(6,1:4))... % Mxy1 = -Mxy2
      (sig{2}.Y(7,[2 3 4 1])+sig{4}.Y(8,1:4))... % Qx1  = -Qy2
      (sig{2}.Y(8,[2 3 4 1])-sig{4}.Y(7,1:4))];  % Qy1  =  Qx2
% go for relative error (condensed from human readable entry that is kept)
test=sig{2}.Y(1:8,[2 3 4 1])+diag([-1 -1 1 -1 -1 1 1 -1])*sig{4}.Y([2 1 3 5 4 6 8 7],1:4);
r3=max(abs(sig{2}.Y(1:8,[2 3 4 1])),[],2);
test=diag(1./r3)*abs(test);
if max(max(abs(test))) > 2e-8; error('Mismatch');end

if 1==2 % for EB
  [C1,mo1.DOF]=fe_mknl('init',mo1);C1=fe_mknl('gstate-struct',mo1,C1); 
 C1.GroupInfo{8}.defe=zeros(24,1);k=fe_mknl('assemble',mo1,C1,def,1);z=C1.GroupInfo{5};
end

%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% #AnalyticStrain Stress computations on orthotropic plate example 
elseif comstr(Cam,'analyticstrain') %t_plate2('analyticStrain')
[CAM,Cam]=comstr(CAM,15);

if carg>nargin; model=ofdemos('composite');
else; model=varargin{carg};carg=carg+1; %#ok<*NOSEM>
end
%model.pl=[100 fe_mat('m_elastic','SI',5) ...
%      38.095e9  9.4e9 0.3   7.5e9 0    0     1640];

% analytic evaluation
if isempty(Cam) % z motion
   data=struct('sel','groupall','dir',{{0,0,'1*x.^3'}}, ...
    'DOF',[.04;.06;.03]);
  RunOpt.Type='z';
elseif comstr(Cam,'xy')
   data=struct('sel','groupall','dir',{{'1*x','1*y'}},'DOF',[.01;.02]);
  RunOpt.Type='xy';
elseif comstr(Cam,'x')
   data=struct('sel','groupall','dir',{{'1*x',0}},'DOF',[.01;.02]);
  RunOpt.Type='x';
elseif comstr(Cam,'s')
   data=struct('sel','groupall','dir',{{'1*y',0}},'DOF',[.01;.02]);
  RunOpt.Type='s';
elseif comstr(Cam,'r') % gamma_a=beta_a+w,a  and b1=theta2
  data=struct('sel','groupall','dir',{{'2*x+y',1,-2}},'DOF',[.03;.04;.05]);%no shear
  data=struct('sel','groupall','dir',{{'2*x+y',0,0}},'DOF',[.03;.04;.05]);%shear -1,2
  RunOpt.Type='r';
end
model=fe_case(model,'dofset','Force',data);
%model.DOF=feutil('getdof',model);
%r2=elem0('VectFromDirAtDof',model,data);
%feplot(model,r2);fecom('scc1')

def=fe_simul('static',model); 
%def.def(:,2:end)=[];%fe_c(q.DOF);ans(:,2)=num2cell(q.def)
%dbstop in fe_mknl at 560

mo1=model;mo1.pl=m_elastic('dbval 100 StrainShell');
mo1.il = p_shell('dbval 110 laminate 100 1 0'); % single ply
mo1=stack_set(mo1,'info','StressCritFcn', ...
    '[r1,dir]=p_shell(''stresscrit'',r1,[.01 1]);DIRS(cEGI(jElt),1:2,jDef)=dir(:,1);');
data={fe_stress('stress -atcenter',mo1,def)};
[C1,mo1.DOF]=fe_mknl('init',mo1);C1.GroupInfo{end}.defe=zeros(24,1);
C1=fe_mknl('gstate-struct',mo1,C1); % fill gstate structure
if 1==2 % Manual test of orient
   IA=C1.GroupInfo{7};
   IA.NodePos=ones(4,1)*(1:9);
   IA.data=IA.data(:,1)*ones(1,9);IA.data(4,:)=model.Elt(2:end,8)';
   IA.lab{4}='the';C1.GroupInfo{7}=IA;
   sp_util('diag',11)
   d1=def;k=fe_mknl('assemble',mo1,C1,d1,1);data{2}=C1.GroupInfo{5};   
end
d1=def;k=fe_mknl('assemble',mo1,C1,d1,1);data{2}=C1.GroupInfo{5};C0=C1;


%cf=feplot(mo1);cf.def=def;feplot('colordataelt',data{1});fecom('colorbar')
mo1.il = p_shell('dbval 110 laminate 100 1 90'); % single ply
%MAP=struct('dir',[0 -1 0],'DOF',[.01;.02;.03]);
MAP=struct('dir',{{0 -1 0}},'lab',{{'v1x','v1y','v1z'}});
mo1=feutil('setpro 110',mo1,'MAP',MAP);
%mo1=stack_set(mo1,'MAP','Group_1',MAP);
mo1=stack_set(mo1,'info','StressCritFcn', ...
    '[r1,dir]=p_shell(''stresscrit'',r1,[.01 1]);DIRS(cEGI(jElt),1:2,jDef)=dir(:,1);');
data{3}=fe_stress('stress -atcenter',mo1,def);
[C1,mo1.DOF]=fe_mknl('init',mo1);C1.GroupInfo{end}.defe=zeros(24,1);
C1=fe_mknl('gstate-struct',mo1,C1); % fill gstate structure
d1=def;k=fe_mknl('assemble',mo1,C1,d1,1);data{4}=C1.GroupInfo{5};


if strcmp(RunOpt.Type,'x') % constant membrane strain
    % 'no orient'
    r2=squeeze(mean(data{2}.Y,2));r2=r2(1:3,:);
    % 'orient'
    r2=squeeze(mean(data{4}.Y,2));r2=r2(1:3,:);
    IA=C1.GroupInfo{7};i1=find(IA.data(ismember(IA.lab,{'v1x','v1y','v1z'}),1));
    r2(i1,:)=r2(i1,:)-1;
    if norm(r2)>1e-10; error('Incorrect computation'); end

    return;%feplot(model,def)
elseif strcmp(RunOpt.Type,'xy') % constant membrane strain

    r1=reshape(data{2}.Y,9,[]);r1(1:2,:)=r1(1:2,:)-1;r1(9,:)=[];
    if norm(r1)>1e-10; error('Incorrect computation'); end
    r1=reshape(data{4}.Y,9,[]);r1(1:2,:)=r1(1:2,:)-1;r1(9,:)=[];
    if norm(r1)>1e-10; error('Incorrect computation'); end
    
    return;%feplot(model,def)
elseif strcmp(RunOpt.Type,'s') % constant membrane strain

    % shear in variable orient
    r1=reshape(data{2}.Y,9,[]);r1(3,:)=r1(3,:)-1;r1(9,:)=[];
    %if norm(r1)>1e-10; error('Incorrect computation'); end
    r1=reshape(data{4}.Y,9,[]);r1(3,:)=r1(3,:)+1;r1(9,:)=[];
    if norm(r1)>1e-10; error('Incorrect computation'); end
    
    return;%feplot(model,def)
elseif strcmp(RunOpt.Type,'r') % shear verification
    r1=reshape(data{2}.Y,9,[]);r1(7,:)=r1(7,:)-1;r1(8,:)=r1(8,:)+2;
    if norm(r1)>1e-10; error('Incorrect computation'); end
    IA=C1.GroupInfo{7};i1=sign(IA.data(ismember(IA.lab,{'v1y'}),1));
    r1=reshape(data{4}.Y,9,[]);r1(7,:)=r1(7,:)-1*i1;r1(8,:)=r1(8,:)+2*i1;
    if norm(r1)>1e-10; error('Incorrect computation'); end
    return

elseif norm((data{3}.data./data{1}.data)-1)>1e-10; 
    disp([data{3}.data data{1}.data])
    error('Mismatch in computed stress');
end

%[data{2}.Y(:,1:4) data{4}.Y(:,1:4)]

%  - - - - - - - - - -
  
% evaluate stress on ply
def=fe_eig(model);
mo1=model;mo1.pl=m_elastic('dbval 100 StrainShell');
mo1.il = p_shell('dbval 110 laminate 100 1 0'); % single ply
mo1=stack_set(mo1,'info','StressCritFcn','r1=p_shell(''stresscrit'',r1,.01);');
d1=fe_def('subdef',def,7);data=fe_stress('stress -atcenter',mo1,d1);
cf=feplot(model);cf.def=d1;feplot('colordataelt',data);

% q4cs stress evaluation with matrix assembly
[C1,mo1.DOF]=fe_mknl('init',mo1);C1.GroupInfo{end}.defe=zeros(24,1);
C1=fe_mknl('gstate-struct',mo1,C1); % fill gstate structure
d1=fe_def('subdef',def,7);k=fe_mknl('assemble',mo1,C1,d1,1);


% ------------------------------------------------------------
% Verification of oriented Jacobian - - - - - - - - - - - - - -
elseif comstr(Cam,'jacobian') %t_plate2('jacobian')

RO.opt=-2;RO.node=[0.5,0.25,0,9;0.75,0.5,0,16;0.5,0.75,0,11;0.35,0.5,0,10];
rule=integrules('quad4',RO.opt);rule.J=zeros(4,4);rule.bas=zeros(9,4);
nodeE=RO.node;
nodeEt= comstr({'x','y','z','Id'},-32);
rule.nodeE=nodeE;rule.nodeEt=int32(nodeEt);of_mk('buildndn',23,rule);
r2=integrules('quad4',RO.opt);r2.J=zeros(4,4);r2.bas=zeros(9,4);
nodeE=[RO.node rule.bas(1:3,:)'];
r2.nodeEt= comstr({'x','y','z','Id','v1x','v1y','v1z'},-32);
r2.nodeE=nodeE;of_mk('buildndn',23,r2);
if norm(rule.NDN-r2.NDN)+norm(rule.J-r2.J)>1e-10;error('mismatch');end
%   sd('compare',r2,rule)
 

% ------------------------------------------------------------
% Displays of orientation maps
elseif comstr(Cam,'maps') %t_plate2('maps')

model=femesh(strcat('testquad4 struct divide 2 10 back'));
model.Node(:,5:6)=model.Node(:,5:6)*diag([.15 .01]);
model.il = p_shell('dbval 110 laminate 100 1 0'); % single ply
Range=45;j1=1;
model=stack_set(model,'MAP','Group_1', ...
     struct('dir',[cos(Range(j1)*pi/180) sin(Range(j1)*pi/180) 0]));
cf=feplot(model);
feutilb('shellmapnodepos',cf)


% ------------------------------------------------------------
%Pour la definition des orientations par elements, voici les utilitaires qui
% me semblent utiles :

model=ofdemos('composite');
% Define material angle based on direction at element
MAP=feutil('getnormalElt MAP -dir1',model);
bas=basis('rotate',[],'rz=30;',1);
MAP.normal=MAP.normal*reshape(bas(7:15),3,3)';
model=p_shell('setTheta',model,MAP);
if norm(model.Elt(2:end,8)-30)>1e-10; error('Mismatch');end

% Obtain a MAP of material orientations
MAP=feutil('getnormalElt MAP -dir1',model);
feplot(model);fecom('showmap',MAP)

% Set elementwise material angles using directions given at nodes. 
% Here a global direction
MAP=struct('normal',ones(size(model.Node,1),1)*bas(7:9), ...
    'ID',model.Node(:,1),'opt',2);
model=p_shell('setTheta',model,MAP);
if norm(model.Elt(2:end,8)-30)>1e-10; error('Mismatch');end

% Using an analytic expression to define components of 
% material orientation vector at nodes
data=struct('sel','groupall','dir',{{'x-0','y+.01',0}},'DOF',[.01;.02;.03]);
model=p_shell('setTheta',model,data);
MAP=feutil('getnormalElt MAP -dir1',model);
feplot(model);fecom('showmap',MAP)

end
