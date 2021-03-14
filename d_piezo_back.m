function [out,out1,out2]=d_piezo(varargin); %#ok<*STOUT>

% D_PIEZO Support for demonstrations related to piezo-electricity
%
% Arnaud Deraemaeker, ULB and Etienne Balmes, SDTools


%       Copyright (c) 1990-2014 by SDTools, All Rights Reserved.
%       For revision information use d_piezo('cvs')

if nargin==0
 sdtweb _taglist d_piezo % see structure of d_piezo file    
 return
end

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*DEFNU,*NOSEM>
[CAM,Cam]=comstr(varargin{1},1);carg=2;
if carg>nargin||~isstruct(varargin{2});RO=struct('info',1);
else;RO=varargin{carg};carg=carg+1;
end

%% #Script -------------------------------------------------------------------
if comstr(Cam,'script');[CAM,Cam]=comstr(CAM,7);

if comstr(Cam,'patchshear')
%% #ScriptPatchShear : Demo script : see sdtweb pz_basics#pz_patch_num_shear 
%
%% BeginSource sdtweb('_example','pz_theory.tex#pz_patch_num_shear')

% See full example as MATLAB code in d_piezo('ScriptPatchShear')
% Build mesh - Meshing script can be viewed with sdtweb d_piezo('MeshPatch')
model=d_piezo('MeshPatch lx=1e-2 ly=1e-2 h=2e-3 nx=1 ny=1 nz=1');

% Define electrodes
model=p_piezo('ElectrodeMPC Top -ground',model,'z==2e-3');
model=p_piezo('ElectrodeMPC Bottom -Input "Free patch"',model,'z==0');

% Rotate basis to align poling direction with y (-90deg around x)
model.bas=basis('rotate',[],'rx=-90',1); %create local basis with id=1
model=feutil('setpro 1 COORDM=1',model); % assign basis with id=1 to pro=1
% Compute response due to applied difference of potential at very low freq
% to avoid rigid body mode
model=stack_set(model,'info','Freq',10);
def=fe_simul('dfrf',model); def.lab={'Free patch, shear'};
def.fun=[0 1]; def=rmfield(def,'data');

% Append mechanically constrained structure
%    can't call fe_simul because no free DOF
%    see code with sdtweb d_piezo('scriptFullConstrain')
def=d_piezo('scriptFullConstrain',model,def);
def.lab{2}='Constrained patch, shear';

% Deformed shape
cf=feplot(model,def); fecom('undef line');

% Electric field representation
p_piezo('viewElec EltSel "matid2" DefLen 20e-4 reset',cf);
% Decompose constitutive law
CC=p_piezo('viewdd -struct',cf); %

% Display and compute mean strains
a=p_piezo('viewstrain -curve -mean',cf); % Strain S
fprintf('Relation between mean strain on free structure and d_24\n');
E3=a.Y(9,1); disp({'E3 mean' a.Y(9,1)  1/2e-3 'E3 analytic'})

disp([a.X{1}(4) num2cell([a.Y(4,1)/E3 CC.d(2,4)']) ...
     {'d_24'}])

% Display and compute mean stresses
b=p_piezo('viewstress -curve -mean',cf); % Stress T
fprintf('Relation between mean stress on pure electric and e_24 \n');
disp([b.X{1}(4) num2cell([b.Y(4,2)/-E3 CC.e(2,4)']) ...
     {'e_24'}])
% Mean stress/strain
disp([b.X{1} num2cell(b.Y(:,2)) num2cell(a.Y(:,1)) a.X{1}])

% Theoretical values of Capacitance and charge density - free patch
CT=CC.epst_r(2,2)*8.854e-12*1e-2*1e-2/2e-3; %% Capacitance - free patch
CdensT=CC.epst_r(2,2)*8.854e-12/2e-3*1e12; %% charge density - free patch

% Theoretical values of Capacitance and charge density - constrained patch
CS=CC.epss_r(2,2)*8.854e-12*1e-2*1e-2/2e-3; %% Capacitance - constrained patch
CdensS=CC.epss_r(2,2)*8.854e-12/2e-3*1e12; %% charge density - constrained patch

% Represent charge density (C/S) value on the electrodes
%   - compare with analytical values
cut=p_piezo('electrodeviewcharge',cf,struct('EltSel','matid 1'));
b=fe_caseg('stressobserve',cut,cf.def);b=reshape(b.Y,[],2);
disp([{'','CdensT','CdensS'};{'Numeric';'Theoretical'} ...
  num2cell([mean(abs(b));CdensT CdensS])])

% Compute the value of the total charge (from reaction at electrical dof)
% Ccompare with analytical values
p_piezo('electrodeTotal',cf) %
disp('Theoretical values of capacitance')
disp([{'CT';'CS'} num2cell([CT;CS])])

%% EndSource

elseif comstr(Cam,'patch')
%% #ScriptPatch : Demo script for patch, see sdtweb pz_basics#pz_patch_num_ext 
%
%% BeginSource sdtweb('_example','pz_theory.tex#pz_patch_num_ext')

% See full example as MATLAB code in d_piezo('ScriptPatch')
% Build mesh - Meshing script can be viewed with sdtweb d_piezo('MeshPatch')
model=d_piezo('MeshPatch lx=1e-2 ly=1e-2 h=2e-3 nx=1 ny=1 nz=1');
% Define electrodes
model=p_piezo('ElectrodeMPC Top -ground',model,'z==2e-3');
model=p_piezo('ElectrodeMPC Bottom -Input "Free patch"',model,'z==0');
p_piezo('TabInfo',model)
model.pl=m_piezo('dbval 1 -elas 2 PIC_255');
p_piezo('TabDD',model) % creates the table with full set of matrices
% Compute response due to applied difference of potential at very low freq
% to avoid rigid body mode
model=stack_set(model,'info','Freq',10);
def=fe_simul('dfrf',model); def.lab={'Free patch, axial'};
def.fun=[0 1];def=feutil('rmfield',def,'LabFcn','data');

% Append mechanically constrained structure
%    can't call fe_simul because no free DOF
%    see code with sdtweb d_piezo('scriptFullConstrain')
def=d_piezo('scriptFullConstrain',model,def);
def.lab{2}='Constrained patch, axial';
% Deformed shape
cf=feplot(model,def);
% Electric field representation
p_piezo('viewElec EltSel "matid2" DefLen 20e-4 reset',cf);
fecom('colormap',[1 0 0]);fecom('undef line');iimouse('resetview');fecom('viewcv');
% Decompose constitutive law
CC=p_piezo('viewdd -struct',cf); %
% Display and compute mean strains
a=p_piezo('viewstrain -curve -mean',cf); % Strain S
fprintf('Relation between mean strain on free structure and d_3i\n');
E3=a.Y(9,1); disp({'E3 mean' a.Y(9,1)  1/2e-3 'E3 analytic'})

disp([a.X{1}(1:3) num2cell([a.Y(1:3,1)/E3 CC.d(3,1:3)']) ...
     {'d_31';'d_32';'d_33'}])
% Display and compute mean stresses
b=p_piezo('viewstress -curve -mean',cf); % Stress T
fprintf('Relation between mean stress on pure electric and e_3i\n');
disp([b.X{1}(1:3) num2cell([b.Y(1:3,2)/-E3 CC.e(3,1:3)']) ...
     {'e_31';'e_32';'e_33'}])

% Mean stress/strain
disp([b.X{1} num2cell(b.Y(:,2)) num2cell(a.Y(:,1)) a.X{1}])
% Theoretical values of Capacitance and charge density - free patch
CT=CC.epst_r(3,3)*8.854e-12*1e-2*1e-2/2e-3; %% Capacitance - free patch
CdensT=CC.epst_r(3,3)*8.854e-12/2e-3*1e12; %% charge density - free patch

% Theoretical values of Capacitance and charge density - constrained patch
CS=CC.epss_r(3,3)*8.854e-12*1e-2*1e-2/2e-3; %% Capacitance - free patch
CdensS=CC.epss_r(3,3)*8.854e-12/2e-3*1e12; %% charge density - free patch

% Represent charge density (C/S) value on the electrodes
%   - compare with analytical values
cut=p_piezo('electrodeviewcharge',cf,struct('EltSel','matid 1'));
b=fe_caseg('stressobserve',cut,cf.def);b=reshape(b.Y,[],2);
disp([{'','CdensT','CdensS'};{'Numeric';'Theoretical'} ...
  num2cell([mean(abs(b));CdensT CdensS])])

% Compute the value of the total charge (from reaction at electrical dof)
% Compare with analytical values
p_piezo('electrodeTotal',cf) %
disp('Theoretical values of capacitance')
disp([{'CT';'CS'} num2cell([CT;CS])])

%% EndSource

elseif comstr(Cam,'pz_patch_num_IDE');
%% #ScriptPatch_num_IDE : see sdtweb pz_tuto#pz_patch_num_IDE

%% BeginSource sdtweb('_example','pz_tuto.tex#pz_patch_num_IDE')

% See full example as MATLAB code in d_piezo('ScriptPz_Patch_Num_IDE')
% Meshing script can be viewed with sdtweb d_piezo('MeshIDEPatch') 
% Build mesh, electrodes and actuation
model=d_piezo('MeshIDEPatch nx=10 ny=5 nz=14 lx=400e-6 ly=300e-6 p=700e-6');

% Compute response due to applied difference of potential
% low freq response to avoid rigid body modes
model=stack_set(model,'info','Freq',10);
def=fe_simul('dfrf',model); def.lab={''};
def.fun=[0 1]; def=rmfield(def,'data');

% Plot deformed shape
cf=feplot(model,def); fecom('view3'); fecom('viewy-90'); fecom('viewz+90')
fecom('undef line'); fecom('triax')
% visualize electric field
cf.sel(1)={'groupall','colorface none -facealpha0 -edgealpha.1'};
p_piezo('viewElec EltSel "matid1" DefLen 50e-6 reset',cf);
fecom('scd 1e-10')
p_piezo('electrodeview -fw',cf); % to see the electrodes on the mesh
% Decompose constitutive law
CC=p_piezo('viewdd -struct',cf); %
% Compute mean value of fields and deduce equivalent d_ij
% Uniform field is assumed for analytical values
a=p_piezo('viewstrain -curve -mean',cf); % mean value of S1-6 and E1-3
fprintf('Relation between mean strain on free structure and d_3i\n');
E3=a.Y(9,1); disp({'E3 mean' a.Y(9,1) -1/700e-6 'E3 analytic'})
disp([a.X{1}(1:3) num2cell([a.Y(1:3,1)/E3 CC.d(3,1:3)']) ...
{'d_31';'d_32';'d_33'}])
% total charge on the electrodes
p_piezo('electrodeTotal',cf)
% charge density on the electrodes
cut=p_piezo('electrodeviewcharge',cf,struct('EltSel','matid 1'));
fecom('view3'); fecom('viewy-90'); fecom('viewz+90')
% Theoretical capacitance for uniform field
Ct=CC.epst_r(3,3)*8.854e-12*400e-6*300e-6/700e-6;
% total charge on the electrodes = capacitance (unit voltage)
C=p_piezo('electrodeTotal',cf);
% Differences are due to non-uniform field, this is to be expected
disp({'C_{IDE}' cell2mat(C(2,2)) Ct 'C analytic'})

%% EndSource

elseif comstr(Cam,'fullconstrain');
%% #ScriptFullConstrain : no mechanical displacement and enforced potential 
model=RO;

% Build analytic expression of displacement

r1=[min(model.Node(:,7)) max(model.Node(:,7))];
data=struct('sel','groupall','dir',{{'x*0','y*0','z*0' ...
    sprintf('(z-%.15g)/%.15g',r1(2),-r1(2)+r1(1))}}, ...
    'DOF',[.01;.02;.03;.21]);
def=elem0('VectFromDirAtDof',model,data,model.DOF);
def.name='Constrained Patch';

if carg<=nargin % Clean combine with earlier deformation
  d1=varargin{carg};carg=carg+1;
  if ~isfield(d1,'lab');d1.lab={d1.name}; end
  def.lab={def.name};def=feutil('rmfield',def,'sel','dir','name');
  def=fe_def('appenddef',d1,def);def.name='Reference solutions';
end

out=def;

elseif comstr(Cam,'pz_plate_4pzt_single');
%% #ScriptPlate_4pzt_single : see sdtweb pz_tuto#pz_plate_4pzt_single

%% BeginSource sdtweb('_example','pz_tuto.tex#pz_plate_4pzt_single')

% See full example as MATLAB code in d_piezo('ScriptPz_plate_4pzt_single')
model=d_piezo('MeshULBplate');  % creates the model
model=fe_case(model,'FixDof','Cantilever','x==0');
% Set modal default zeta based on loss factor of material 1
model=stack_set(model,'info','DefaultZeta',feutilb('getloss',1,model)/2);
cf=feplot(model); fecom('colordatagroup'); set(gca,'cameraupvector',[0 1 0])
p_piezo('TabDD',model)   % List piezo constitutive laws
p_piezo('TabInfo',model) % List piezo related properties
model=fe_case(model,'SensDof','Tip',1054.03); % Displ sensor
model=fe_case(model,'DofSet','V-Act',struct('def',1,'DOF',1682.21)); %Act
model=p_piezo('ElectrodeSensQ  1682 Q-Act',model); % Charge sensors
model=p_piezo('ElectrodeSensQ  1683 Q-S1',model);
model=p_piezo('ElectrodeSensQ  1684 Q-S2',model);
model=p_piezo('ElectrodeSensQ  1685 Q-S3',model);
% Fix dofs 1683-1685 to measure resultant (charge)
model=fe_case(model,'FixDof','SC*1683-1685',[1683:1685]+.21);
sens=fe_case(model,'sens');
% Compute static response due to actuator
d0=fe_simul('dfrf',stack_set(model,'info','Freq',0)); % direct refer frf at 0Hz
cf.def=d0; fecom(';view3;scd .1;colordatagroup;undefline')
% Compute frequency response function (full model)
if exist('mklserv_utils','file')&&sdtdef('isinteractive')
  ofact('method mklserv_utils -silent');
  f=linspace(1,100,400); % in Hz
else;
  f=linspace(1,100,50); % in Hz (just 50 points to make it fast)
end

d1=fe_simul('dfrf',stack_set(model,'info','Freq',f(:))); % direct refer frf
% Project response on sensors
C1=fe_case('SensObserve',sens,d1);C1.name='DFRF';C1.Ylab='V-Act';
C1.Xlab{1}={'Frequency','Hz'};
% state-space model
[s1,TR1]=fe2ss('free 5 10 0',model); %
C2=qbode(s1,f(:)*2*pi,'struct');C2.name='SS';

% Compare the two curves
C2.X{2}=sens.lab; C1.X{3}=nor2ss('lab_in',s1);C2.X{3}=nor2ss('lab_in',s1);%
C2=feutil('rmfield',C2,'lab'); C1.Ylab=C2.Ylab;
ci=iiplot; iicom(ci,'curveinit',{'curve',C1.name,C1;'curve',C2.name,C2});
iicom('submagpha');

%% EndSource

elseif comstr(Cam,'pz_plate_pzcomb_2');
%% #ScriptPlate_pzcomb_2 :  see sdtweb pz_tuto#pz_plate_pzcomb_2

%% BeginSource sdtweb('_example','pz_tuto.tex#pz_plate_pzcomb_2')

% See full example as MATLAB code in d_piezo('ScriptPz_plate_pzcomb_2')
model=d_piezo('MeshULBplate cantilever');  % creates the model
model=fe_case(model,'DofSet','V*1683-1682', ...
    struct('def',[1;-1],'DOF',[1682;1683]+.21));

% Compute static response
d0=fe_simul('dfrf',stack_set(model,'info','Freq',0)); % direct refer frf
cf=feplot(model); cf.def=d0;
fecom(';view3;scd .1;colordatagroup;undefline')
% Combined charge output (SC electrodes) % difference of charge 1684-1685
r1=struct('cta',[1 -1],'DOF',[1684;1685]+.21,'name','QS3+4');
model=p_piezo('ElectrodeSensQ',model,r1);
% Combined voltage output (OC electrodes) % difference of voltage 1684-1685
r1=struct('cta',[1 -1],'DOF',[1684;1685]+.21,'name','VS3+4');
model=fe_case(model,'SensDof',r1.name,r1);
[sys,TR]=fe2ss('free 5 10 0 -dterm',model);
 C1=qbode(sys,linspace(1,100,400)'*2*pi,'struct'); C1.name='OC';

% Now you need to SC 1684 and 1685 to measure charge resultant
model=fe_case(model,'FixDof','SC*1684-1685',[1684;1685]+.21);
[sys2,TR2]=fe2ss('free 5 10 0 -dterm',model);
C2=qbode(sys2,linspace(1,100,400)'*2*pi,'struct');C2.name='SC';

% invert channels and scale
C1.Y=fliplr(C1.Y); C1.X{2}= flipud(C1.X{2});
C2.Y(:,1)=C2.Y(:,1)*C1.Y(1,1)/C2.Y(1,1);
iicom('curvereset'),iicom('curveinit',{'curve',C1.name,C1;'curve',C2.name,C2 });
model=d_piezo('MeshULBplate -cantilever');
% Open circuit : do nothing on electrodes
d1=fe_eig(model,[5 20 1e3]);
% Short circuit : fix all electric DOFs
DOF=p_piezo('electrodeDOF.*',model);
d2=fe_eig(fe_case(model,'FixDof','SC',DOF),[5 20 1e3]);
r1=[d1.data(1:end)./d2.data(1:end)];
figure;plot(r1,'*','linewidth',2);axis tight
xlabel('Mode number');ylabel('f_{OC}/f_{SC}');

%% EndSource

elseif comstr(Cam,'pz_plate_mfc');
%% #ScriptPlate_mfc : see sdtweb pz_tuto#pz_plate_mfc

%% BeginSource sdtweb('_example','pz_tuto.tex#pz_plate_mfc')

% See full example as MATLAB code in d_piezo('ScriptPz_plate_MFC')
% Meshing script,open with sdtweb d_piezo('MeshMFCplate')
model=d_piezo('MeshMFCplate -cantilever')  % creates the model
cf=feplot(model); fecom('colordatagroup-EdgeAlpha.1');
data.def=[1 -1;1 1]'; % Define combinations for actuators
data.lab={'V-bend';'V-Tract'};
data.DOF=[1682.21; 1683.21];
model=fe_case(model,'DofSet','V_In',data);
% Add tip displacement sensors in z
model=fe_case(model,'SensDof','Tipt-z',156.03); % Z-disp
model=fe_case(model,'SensDof','Tip-x',156.01); % X-disp
model=fe_case(model,'SensDof','Tipb-z',2964.03); % Z-disp
% Static response
d0=fe_simul('dfrf',stack_set(model,'info','Freq',0)); %
feplot(model); cf=fecom
cf.def=d0; sens=fe_case(model,'sens'); C1=fe_case('SensObserve',sens,d0);
fecom(';view3;scd .1;colordataEvalA -edgealpha.1;undefline')
%% Rotate fibers
model.il(2,[20 44])=[45 -45];
d1=fe_simul('dfrf',stack_set(model,'info','Freq',0)); % static response
cf.def=d1; fecom('scd 1e-2'); C2=fe_case('SensObserve',sens,d1);

%% EndSource


elseif comstr(Cam,'pz_plate_triang');
%% #ScriptPlate_triang : see sdtweb pz_tuto#pz_plate_triang

%% BeginSource sdtweb('_example','pz_tuto.tex#pz_plate_triang')

% See full example as MATLAB code in d_piezo('ScriptPz_Plate_Triang')
% Meshing script can be viewed with sdtweb d_piezo('MeshTrianglePlate') 
% --- requires gmsh
model=d_piezo('MeshTrianglePlate'); 
cf=feplot(model); fecom('colordatapro'); fecom('view2')
%% Voltage actuation
model=fe_case(model,'SensDof','Tip',7.03); % Displ sensor
model=fe_case(model,'DofSet','V-Act',struct('def',1,'DOF',100001.21));

% static response to voltage actuation
d0=fe_simul('dfrf',stack_set(model,'info','Freq',0));
cf.def=d0; fecom('colordataz -alpha .4 -edgealpha .1')
fecom('scd -.03'); fecom('view3');
% Dynamic response (reduced modal model)
[sys,TR]=fe2ss('free 5 20 0 -dterm',model);
C1=qbode(sys,linspace(0,500,1000)'*2*pi,'struct'); C1.name='.';

%% Point load actuation
model=fe_case(model,'Remove','V-Act'); % remove piezo actuator
model=fe_case(model,'FixDof','Top piezo',100001); %SC piezo electrode

% Determine scaling factor, check b/l ratio and build point force
CC=p_piezo('viewdd -struct',model);
a=100; b=33.58;
zm=0.650e-3; V=1; e31=CC.e(1); A=-(e31*zm*V*b)/a;
bl= 2*sqrt(-CC.e(2)/CC.e(1));

data=struct('DOF',[7.03],'def',A); data.lab=fe_curve('datatype',13);
model=fe_case(model,'DofLoad','PointLoad',data);

% Static response to point load
d1=fe_simul('dfrf',stack_set(model,'info','Freq',0));
ind=fe_c(d1.DOF,7.03,'ind'); d1p=d1.def(ind);

% Dynamic response (reduced modal model)
[sys,TR]=fe2ss('free 5 20 0 -dterm',model);
C2=qbode(sys,linspace(0,500,1000)'*2*pi,'struct'); C2.name='-';

% Compare frequency responses
ci=iiplot;
iicom(ci,'curveinit',{'curve',C1.name,C1;'curve',C2.name,C2});
iicom('submagpha');

%% EndSource



elseif comstr(Cam,'pz_plate_4pzt_comb');
%% #ScriptPlate_4pzt_comb : DOC script plate with compination of electrodes
% see sdtweb pz_tuto#pz_plate_4pzt_comb

%% BeginSource sdtweb('_example','pz_tuto.tex#pz_plate_4pzt_comb')

% See full example as MATLAB code in d_piezo('ScriptPz_plate_4pzt_comb')
model=d_piezo('MeshULBplate -cantilever');  % creates the model

% combine electrodes to generate pure bending / pure traction
data.def=[1 -1 1 -1;1 1 1 1]'; % Define combinations for actuators
data.lab={'V-bend';'V-Tract'};
data.DOF=p_piezo('electrodeDOF.*',model);
model=fe_case(model,'DofSet','V_In',data);

% Compute static response
d0=fe_simul('dfrf',stack_set(model,'info','Freq',0)); % direct refer frf
cf=feplot(model); cf.def=d0;
fecom(';view3;scd .02;colordataEvalZ;undefline')
% Add tip displacement sensor in x and z
model=fe_case(model,'SensDof','Tip-z',1054.03); % Z-disp
model=fe_case(model,'SensDof','Tip-x',1054.01); % X-disp

% Make SS model and display FRF
[sys,TR]=fe2ss('free 5 30 0 -dterm',model);
C1=qbode(sys,linspace(1,100,400)'*2*pi,'struct');
C1.name='Bend-tract combination'; % Force name
C1.X{2}={'Tip-z';'Tip-x'}; % Force input labels
C1.X{3}={'V-bend';'V-tract'}; % Force output labels
ci=iiplot;iicom('CurveReset');iicom('curveinit',C1)

%% EndSource


elseif comstr(Cam,'pz_accA');
%% #pz_accA See sdtweb pz_tuto#pz_accA -2

%% BeginSource sdtweb('_example','pz_tuto.tex#pz_accA')

% See full example as MATLAB code in d_piezo('ScriptPz_accA')
% Meshing script can be viewed with sdtweb d_piezo('MeshBaseAccel') 
model=d_piezo('MeshBaseAccel');
cf=feplot(model); fecom('colordatagroup');
set(gca,'cameraposition',[-0.0604   -0.0787    0.0139])
% -MatID 2 requests a charge resultant sensor
% -vout requests a voltage sensor
model=p_piezo('ElectrodeMPC Top sensor -matid 2 -vout',model,'z==0.004');
% -ground generates a v=0 FixDof case entry
model=p_piezo('ElectrodeMPC Bottom sensor -ground',model,'z==0.003');
% Add a displacement sensor for the basis
model=fe_case(model,'SensDof','Base-displ',1.03);
% Remove the charge sensor (not needed)
model=fe_case(model,'remove','Q-Top sensor');
% Normal surface force (pressure) applied to bottom of wear plate for excitation:
data=struct('sel','z==0','eltsel','groupall','def',1e4,'DOF',.19);
model=fe_case(model,'Fsurf','Bottom excitation',data);

% Other parameters
  model=stack_set(model,'info','Freq',logspace(3,5.3,200)'); % frequencies for computation

 % compute response
 ofact('silent'); d1=fe_simul('dfrf',model);

 % Project on sensor
 sens=fe_case(model,'sens');

 % Build a clean "curve" for iiplot display
 C1=fe_case('SensObserve',sens,d1);C1.name='DFRF';C1.Ylab='Base-Exc';
 C1.X{2}={'Sensor output(V)';'Base Acc(m/s^2)';'Sensitivity (V/m/s^2)'};
 C1.Y(:,2)=C1.Y(:,2).*(-(C1.X{1}(:,1)*2*pi).^2); % Base acc =disp.*-w.^2
 C1.Y(:,3)=C1.Y(:,1)./C1.Y(:,2);                 % Sensitivity=V/acc
 C1=sdsetprop(C1,'PlotInfo','sub','1 1','show','abs','scale','xlog;ylog');
 C1.name='Free-Voltage';

 ci=iiplot; iicom(ci,'curveInit',C1.name,C1);iicom ch3;
% Remove pressure
 model=fe_case(model,'remove','Bottom excitation')

% Link dofs of base and impose unit vertical displacement
rb=feutilb('geomrb',feutil('getnode z==0',model),[0 0 0], ...
     feutil('getdof',model));
rb=fe_def('subdef',rb,3); % Keep vertical displacement
model=fe_case(model,'DofSet','Base',rb);

 % compute
 ofact('silent'); model.DOF=[]; d1=fe_simul('dfrf',model);

 % Project on sensor and create output
 sens=fe_case(model,'sens');
 C2=fe_case('SensObserve',sens,d1);C2.name='DFRF';C2.Ylab='Imp-displ';

 % Build a clean "curve" for iiplot display
 C2.X{2}={'Sensor output(V)';'Base Acc(m/s^2)';'Sensitivity (V/m/s^2)'};
 C2.XLab{1}='Freq [Hz]';
 C2.Y(:,2)=C2.Y(:,2).*(-(C2.X{1}*2*pi).^2); % Base acc
 C2=sdsetprop(C2,'PlotInfo','sub','1 1','show','abs','scale','xlog;ylog');
 C2.Y(:,3)=C2.Y(:,1)./C2.Y(:,2);% Sensitivity
 C2.name='Imp-Voltage';
 C2=feutil('rmfield',C2,'Ylab'); C1=feutil('rmfield',C1,'Ylab');
 ci=iiplot; iicom(ci,'curveinit',{'curve',C1.name,C1;'curve',C2.name,C2});
  % Meshing script,open with sdtweb d_piezo('MeshBaseAccel')
model=d_piezo('MeshBaseAccel');
model=fe_case(model,'remove','V-Top sensor');

% Short-circuit electrodes of accelerometer
model=fe_case(model,'FixDof','V=0 on Top Sensor', ...
    p_piezo('electrodedof Top sensor',model));

% Other parameters
  model=stack_set(model,'info','Freq',logspace(3,5.3,200)');

% Link dofs of base and impose unit vertical displacement
rb=feutilb('geomrb',feutil('getnode z==0',model),[0 0 0], ...
     feutil('getdof',model));
rb=fe_def('subdef',rb,3); % Keep vertical displacement
model=fe_case(model,'DofSet','Base',rb);

 % compute
 ofact('silent'); model.DOF=[]; d1=fe_simul('dfrf',model);

 % Project on sensor and create output
 sens=fe_case(model,'sens');
 C4=fe_case('SensObserve',sens,d1);C4.name='DFRF';C4.Ylab='Imp-displ';

 % Build a clean "curve" for iiplot display
 C4.X{2}={'Sensor output(C)';'Base Acc(m/s^2)';'Sensitivity (C/m/s^2)'};
 C4.XLab{1}='Freq [Hz]';
 C4.Y(:,2)=C4.Y(:,2).*(-(C4.X{1}*2*pi).^2); % Base acc
 C4=sdsetprop(C4,'PlotInfo','sub','1 1','show','abs','scale','xlog;ylog');
 C4.Y(:,3)=C4.Y(:,1)./C4.Y(:,2);% Sensitivity
 C4.name='Imp-Charge';

 % Normalize the sensitivities to plot on same graph
 C6=C2; % save C6 as non-normalized
 C2.Y(:,3)=C2.Y(:,3)./C2.Y(1,3);C4.Y(:,3)=C4.Y(:,3)./C4.Y(1,3);
 C2=feutil('rmfield',C2,'Ylab'); C4=feutil('rmfield',C4,'Ylab');
 iicom(ci,'curvereset');
 iicom(ci,'curveinit',{'curve',C2.name,C2;'curve',C4.name,C4});
 iicom('ch 3')

%% EndSource


elseif comstr(Cam,'pz_acc_shaker');
%% #pz_acc_shaker See sdtweb pz_tuto#pz_acc_shaker -2

%% BeginSource sdtweb('_example','pz_tuto.tex#pz_acc_shaker')

% See full example as MATLAB code in d_piezo('ScriptPz_acc_shaker')
% Meshing script,open with sdtweb d_piezo('MeshPiezoShaker')
model=d_piezo('MeshPiezoShaker');
cf=feplot(model); fecom('colordatapro');
set(gca,'cameraposition',[-0.0604   -0.0787    0.0139])
  % -input "In" says it will be used as a voltage actuator
model=p_piezo('ElectrodeMPC Top Actuator -input "Vin-Shaker"',model,'z==-0.01');
  % -ground generates a v=0 FixDof case entry
model=p_piezo('ElectrodeMPC Bottom Actuator -ground',model,'z==-0.012');

% Voltage sensor will be used - remove charge sensor
model=fe_case(model,'remove','Q-Top sensor');

% Frequencies for computation
model=stack_set(model,'info','Freq',logspace(3,5.3,200)');

 % compute - full model, voltage input on shaker
 ofact('silent'); model.DOF=[]; d1=fe_simul('dfrf',model);

 % Project on sensor and create output
 sens=fe_case(model,'sens');
 C5=fe_case('SensObserve',sens,d1);C5.name='DFRF';C5.Ylab='Shaker-Exc';

 % Build a clean "curve" for iiplot display
 C5.X{2}={'Sensor output(V)';'Base Acc(m/s^2)';'Sensitivity (V/m/s^2)'};
 C5.XLab{1}='Freq [Hz]';
 C5.Y(:,2)=C5.Y(:,2).*(-(C5.X{1}*2*pi).^2); % Base acc
 C5=sdsetprop(C5,'PlotInfo','sub','1 1','show','abs','scale','xlog;ylog');
 C5.Y(:,3)=C5.Y(:,1)./C5.Y(:,2);% Sensitivity
 C5.name='Shaker-Voltage'; ci=iiplot;
 C5=feutil('rmfield',C5,'Ylab');
 iicom('curveinit',C5); iicom('ch 3');
 iicom('curvereset'); C6=feutil('rmfield',C6,'Ylab');
iicom('curveinit',{'curve',C1.name,C1;'curve',C6.name,C6;'curve',C5.name,C5});
iicom('ch 3')

%% EndSource



%%
else; error('Script%s unknown',CAM);
    
end
%% #Mesh -------------------------------------------------------------------
elseif comstr(Cam,'mesh');[CAM,Cam]=comstr(CAM,5);

if comstr(Cam,'patch');[CAM,Cam]=comstr(CAM,6);
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
%% #MeshPatch : generic patch using volume elements
if carg<=nargin&&isstruct(varargin{carg});RO=varargin{carg};carg=carg+1;
else; RO=struct;
end

if isempty(CAM)&&isempty(fieldnames(RO))
%% Generic meshing of pre-defined patches 
% d_piezo('meshpatch','SmartM.MFC-P1.2814')
 RO=m_piezo('patch',varargin{carg});carg=carg+1;  
 if strcmpi(RO.shape,'rect')
   RO.lx=RO.dim(1); RO.ly=RO.dim(2); RO.h=RO.dim(3);
 end
 CAM='';  
end     
% Parameter handling
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   ' lx(10e-3#%g#"patch width/radius")'...
   ' ly(10e-3#%g#"patch length")'...
   ' h(2e-3#%g#"patch height")'...
   ' lcx(#%g#"characteristic length in plane")'...
   ' lcz(#%g#"characteristic length out of plane")'...
   ' nx(1#%g#"n_elt x direction/radial")'...
   ' ny(1#%g#"n_elt y direction")'...
   ' nz(1#%g#"n_elt z direction")'...
   ' Quad(#3#"quadratic elements if present")'...
   ' shape(rect#%s#"define shape rect/circle")'...
   ],{RO,CAM});

if ~isfield(RO,'il');RO.il=p_solid('dbval 2 d3 -3'); end % Volume by default

%% Mesh
if strcmpi(RO.shape,'rect')
 model=feutil('objectbeam 1 1',[0 0 0;RO.lx 0 0],RO.nx); % 2D
 model=feutil(sprintf('extrude %i 0 %.15g 0',RO.ny,RO.ly/RO.ny),model);
else
  error('Shape %s not yet supported',RO.shape);  
end

%% Define material properties
if isfield(RO,'pl'); model.pl=RO.pl;
else
 model.pl=m_piezo('dbval 1 -elas 2 SONOX_P502_iso');
end

% Now extrude based on il information
if strcmpi(fe_mat('typep',RO.il(1,2)),'p_solid') % Monolitic
  model=feutil(sprintf('extrude %i 0 0 %.15g',RO.nz,RO.h/RO.nz),model);
  RO.zelec=[0 RO.h];
elseif strcmpi(fe_mat('typep',RO.il(1,2)),'p_shell'); % Laminate
 [constit,integ,elmap,layer]=p_shell('buildconstit',[0 RO.il(1)],model.pl, ...
     RO.il,model);layer=layer.layer;
 RO.zelec=layer(layer(:,1)==101,2:3);
 
 layer=[layer(:,2);layer(end,3)];layer=feutil('refineline',layer);
 model=feutil('extrude 0 0 0 1',model,layer); 
else
  error('%s Not implemented',fe_mat('type',RO.il(1)));
end


% Set MatId and ProID to whatever material is a m_piezo
pl=feutil('getpl m_piezo',model); RO.MatId=pl(1);
model.Elt=feutil(sprintf('set group1 matid %i ProId %i', ...
    RO.MatId,RO.MatId),model);
% Integration rules for volumes
model=p_solid('default;',model);


if ~isfield(model,'unit'); model.unit='SI'; end
if RO.Quad;model=feutil('lin2quad',model);end

%% Build a MPC defining a single potential for the electrodes
% InputDOF=[];
% [model,InputDOF(end+1,1)]=p_piezo( ...
%     sprintf('ElectrodeMPC -MatId%i Top',RO.MatId),model, ...
%     sprintf('z==%.15g',RO.zelec(2)));
% [model,InputDOF(end+1,1)]=p_piezo('ElectrodeMPC Bottom',model, ...
%     sprintf('z==%.15g',RO.zelec(1)));
% 
% RO.InputDOF=InputDOF;
% model=stack_set(model,'info','MeshInfo',RO);

out=model;  % Send output to out variable
if nargout>1; out1=RO; end % Send MeshInfo back

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
elseif comstr(Cam,'plate')
%% #MeshPlate : generic integration of plate with geometric patches -2
% Requires fe_gmsh
% sdtweb t_piezo2('geo')


if 1==2 % Simple example
 RE=struct('Rect',[15 10 3 20 5 4 10;50 10 3  150 5 4  80]);
 RE.Circ=[100 100 10 4];
 d_piezo('meshplate',RE); 
end
% --------------------------------------------------------------------
%% Parameter handling
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   ' unit("SI"#%s#" unit system")'...
   ' Quad(#3#"quadratic elements if present")'...
   ' MatId(100#%g#"Starting matid for piezos")'...
   ],{RO,CAM(7:end)});
RunOpt.NeedOrient=1;
RO.PatchId=[];
if ~isfield(RO,'list');error('Struct with Circ fields no longer supported');
else
%% #MeshPlateList : Most general call based on list of features -3
 model=[]; 
 if ~isfield(RO,'MatId');RO.MatId=100;end % Offset for new matid
 if ~isfield(RO,'ProId');RO.ProId=100;end % Offset for new proid
 RO.Rect=[]; RO.Circ=[];  RO.Cam0=CAM; RO.Electrode=[];
 for j1=1:size(RO.list,1)
   %% Deal with lam(inate) column
   st=RO.list{j1,3}; if ischar(st);RB=struct;else; RB=st;st='';end 
   RB.CAM=st;

   if strncmpi(RO.list{j1,2},'lam',3)% Skip header 'lam'
   elseif isstruct(RO.list{j1,2})&& ...
           (isfield(RO.list{j1,2},'il')||isfield(RO.list{j1,2},'pl'))
     %% #global Predefined model with at least properties in pl and il -3
     model=RO.list{j1,2}; 
     [RB,st,CAM]=cingui('paramedit -DoClean',[ ...
      ' lx(200#%g#"plate length")'...
      ' ly(300#%g#"plate width")'...
      ' lc(5#%g#"Characteristic length for division")'...
      ' nx(#%g#"n_elt x direction if lc not given")'...
      ' ny(#%g#"n_elt y direction if lc not given")'...
      ' Sens(#3#"keep sensor mesh for main")'...
      ' unit("SI"#%s#" unit system")' ...
      ],{RB,st});
     if ~isfield(model,'Elt')||isempty(model.Elt); 
      %% Mesh base plate if base plate not already meshed
      if isempty(RB.nx);RB.nx=ceil(RB.lx/RB.lc);end
      mo1=feutil('objectbeam 1 1',[0 0 0;RB.lx 0 0],RB.nx);
      if isempty(RB.ny);RB.ny=ceil(RB.ly/RB.lc);end
      mo1=feutil(sprintf('extrude %i 0 %.15g 0',RB.ny,RB.ly/RB.ny),mo1);
      [un1,i2]=sortrows(round(mo1.Node/max(RB.lx,RB.ly)*1e9),[6 5]); % good grid
      mo1.Node=mo1.Node(i2,:); mo1=feutil('renumber -noOri;',mo1);
      model.Node=mo1.Node; model.Elt=mo1.Elt; 
      if RB.Sens
       mo1=feutil('rmfield',mo1,'Stack');mo1.tdof=mo1.Node(:,1)+.03;
       mo1.MeshInfo=RB;
       model=fe_case(model,'sensdof','Main',mo1);
      end
     elseif isfield(RO,'lc')&&isfield(RO,'type')&&strcmpi(RO.type,'global') 
      % Mesh plate from coarse
      mo1=fe_fmesh(sprintf('qmesh %g',RO.lc/2),model);
      model.Node=mo1.Node;model.Elt=mo1.Elt;model=stack_set(model,mo1.Stack);
     end
     if isfield(RB,'Ecoef')
      warning('Ecoef for plate must be checked and generalised XXXeb to discuss')
      [Typ,~,SubTyp]=fe_mat('type',RO.list{2,2}.pl(1,2));
      if comstr(Typ,'m_elastic')&&SubTyp==6;
       ind=[3 4 5 9 10 11];
       model.pl(1,ind)=model.pl(1,ind)*RB.Ecoef;
      else
       error('This case is not handled')
      end
     end
          
   elseif strncmpi(RO.list{j1,2},'mdl',3)
     %% Refine coarse model -3
     model=RB.mdl;
     mdl=fe_fmesh(sprintf('qmesh %g',RB.lc*2),model);
     model.Node=mdl.Node;model.Elt=mdl.Elt;
     
   elseif strncmpi(RO.list{j1,2},'baseid',5)
     %% Combination_of_base+patch information -3

     RD=model;RD.name=RO.list{j1,2};  % xxxeb need clarifyRB/
 
     RD.MatId=RO.MatId; RD.ProId=RO.ProId;
     [RD,RC]=m_piezo('patch',RD); % Mesh the patch, RC gives patch options
     % Added electrodes (should relocate)
     RO.MatId=RC.MatId; RO.ProId=RC.ProId;
     if isfield(RD,'rc')&&RD.rc==0;RD.rc=RB.rc;end % Not robut
     st=fieldnames(RD);for j2=1:length(st);RB.(st{j2})=RD.(st{j2});end     
     model.Node=RC.Node; 
     if isfield(RC,'Electrode')&&~isempty(RC.Electrode)
      RO.Electrode(end+(1:size(RC.Electrode)),1:2)=RC.Electrode;
     end
     
     switch lower(RB.shape)
     case 'circ'
      [model,RO,RB]=meshCirc(model,RO,RB);
     case 'rect'
      [model,RO,RB]=meshRect(model,RO,RB);
     otherwise; error('Shape%s unknown',RB.shape);
     end
         
     if isfield(RB,'pl');model=feutil('setmat',model,RB.pl);end% Set materials 
     if isfield(RB,'Ecoef') % Assume composite with single material
       % model = d_shm('cbMatDef MatID coef',model)
       RB.pl=feutil(sprintf('getpl %i',RB.il(10)),model);
       ind=10:4:size(RB.il,2);
       RB.il(ind(RB.il(ind)==RB.pl(1)))=RB.il(1); % Modify MatId 
       RB.pl(1)=RB.il(1); 
       model=feutil('setpro',model,RB.il); % Set properties
       model=feutil('setmat',model,RB.pl); % Set properties
       model=d_shm(sprintf('CbMatDef %i %.15g',RB.pl(1),RB.Ecoef),model);
     else;
      model=feutil('setpro',model,RB.il); % Set properties
     end
 
   elseif strncmpi(RO.list{j1,2},'meshdef',7) % mesh deformations
    model=fe_fmesh('meshdef',model,RO.list(j1,:));
   elseif strncmpi(RO.list{j1,1},'svs',3) 
    %% Extrude to nida SVS (requires visco license)
     % Assuming top layer to be laminate with piezo
     % extrude to [-z0 volume z0]
%       R2=struct('extrude','svs','lam','bmi_CoreA');
%       model=nida13('geoHGen',stack_set(model,'info','MeshInfo',R2));
%       R2=stack_get(model,'info','MeshInfo','get');R2.list=S1;
%       model=stack_set(model,'info','MeshInfo',R2);
%       RO=stack_get(mo1,'info','MeshInfo','get');
     il0=feutil('getil1',model); % assume proid 1 for shell
     if isempty(il0); error('Can''t find shell of proid 1'); end
     [st5,Unit,SubType]=fe_mat('typep',il0(2));
     if ~comstr(st5,'p_shell'); error('Expect p_shell for ProId1'); end
     if SubType==2; RB.ShellH=-feutil('getpro1 Z0',model)*2;
     elseif SubType==1; RB.ShellH=feutil('getpro1 h',model);
     end
     model.Elt=feutil('Orient 1:1e4 n 0 0 1e5',model);
     RunOpt.NeedOrient=0;
     if isfield(RB,'plcore'); ProIdcore=RB.plcore(1); else; ProIdcore=10; end
     sandCom=sprintf('makesandwich shell 0 0 %.15g volume %i %.15g shell 1 %.15g %.15g', ...
      -RB.ShellH/2,ProIdcore,-RB.height,RB.ShellH*[-.5 .5]);
     m1=fevisco(sandCom,model); 
     model.Node=m1.Node; model.Elt=m1.Elt; % XXX eb can't we propagate all model ?
     if ~isfield(model,'il'); model.il=[]; end
     model.il=p_solid(model.il,'dbval10 d3 -3'); % default for nida
     if ~isfield(RB,'plcore')
      sdtw('_nb','You should have given a .plcore containing core material property for SVS')
      model.pl=m_elastic(model.pl,'dbval10 steel');  % default for nida
     else
      model.pl(end+1,1:length(RB.plcore))=RB.plcore;
     end
     
%       mo1.Elt=mo2.Elt;mo1.Node=mo2.Node;
%       %% deal with mat/pro of SVS, 1 plate, 10 core
%        [pl,il]=nida13('matlam',struct('lam',RO.lam,'piezo',''));pl(2,1)=10;%core10
%        mo1=feutil('setmat',mo1,pl);mo1.pl=fe_mat('convertMM',mo1.pl); 
%        mo1=feutil('setpro',mo1,il);mo1.il=fe_mat('convertMM',mo1.il);
%        mo1.il=p_solid(mo1.il,'dbval 10 d3 -3');
 
   else;
     error('Need rewrite possibility of defect other');
     %dbstack; keyboard
   end
   if ~isfield(RO,'lc')&&isfield(RB,'lc');RO.lc=RB.lc;end
 end % Loop on list
 if isfield(model,'unit');
     model.pl=fe_mat(sprintf('convert%s',model.unit),model.pl);
     model.il=fe_mat(sprintf('convert%s',model.unit),model.il);
 end
 model=stack_set(model,'info','Electrodes',struct('data',RO.Electrode));
 if nargout==0;cf=feplot(model);fecom(cf,'colordatamat-edgealpha.1');
 else; out=model;
 end
 %% End of list building
end



if isfield(RO,'RefineToQuad');
   r2=feutil('refineToQuad',model);model.Node=r2.Node;model.Elt=r2.Elt;
   NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
   for j1=1:size(RO.Circ,1) % make true circles
    n1=feutil(sprintf('getnode matid1 & matid %i',RO.Circ(j1,5)),model);
    r2=RO.Circ(j1*ones(size(n1,1),1),1:2);
    r1=n1(:,5:6)-r2;r1=diag(sparse((RO.Circ(j1,3)./sqrt(r1.^2*[1;1]))))*r1+r2;
    model.Node(NNode(n1(:,1)),5:6)=r1;
   end
end

% Force vertical and x axis orientation for composite support
if RunOpt.NeedOrient
 model.Elt=feutil('orient 1:100 n 0 0 1e5;',model);
end
data=struct('sel','groupall','dir',{{1,0,0}},'DOF',[.01;.02;.03]);
%[Case,CaseName]=fe_case(model,'getcase');
model=p_shell('setTheta;',model,data);

if nargout==0; cf=feplot(model);fecom(cf,'colordatamat-edgealpha.1');
else; out=model;
end
out1=RO;

%% #MeshBaseAccel: basic acceleromter mesh
elseif comstr(Cam,'baseaccel')

% Mesh
model=feutil('object quad 1 1',[0 0 0;5 0 0;0 0 3],4,2);
model=feutil('object quad 2 2',model,[0 0 3;2.5 0 0;0 0 1],2,2);
model=feutil('object quad 4 4',model,[0 0 4;5 0 0;0 0 10],4,6);
model=feutil('Rev 20 o 0 0 0 360 0 0 1',model);
model.Node(:,5:7)=model.Node(:,5:7)/1000; % SI
model.unit='SI';

% @ad : need to put more details on properties
% Define properties
model.pl=[ ...
   % Elastic base - typical properties of Alumine
   1 fe_mat('m_elastic','SI',1) 4e11 0.22 3965 zeros(1,25);  
   % Piezo properties, see sdtweb('m_piezo#Database')
    m_piezo('dbval 2 -elas 3 SONOX_P502_iso');
   % Steel backmass
   m_elastic('dbval 4 steel') zeros(1,24)];
        
model=p_solid('default;',model);
model.name='Base Accel';

% Possibly use a quadratic mesh
if ~isempty(strfind(Cam,'quad'));model=feutil('lin2quad',model);end

%%% Electrodes and boundary conditions
% See sdtweb p_piezo#ElectrodeMPC
% Build MPCs defining a single potential for the electrodes 
  % -input "In" says this could be used as a voltage actuator
  % -MatID 2 requests a charge resultant sensor
  % -vout requests a voltage sensor
model=p_piezo('ElectrodeMPC Top sensor -matid 2 -vout',model,'z==0.004');
  % -ground generates a v=0 FixDof case entry
model=p_piezo('ElectrodeMPC Bottom sensor -ground',model,'z==0.003');

% Add a displacement sensor for the basis
model=fe_case(model,'SensDof','Base-displ',1.03); % Displ sensor


% Default frequency range for response computation
model=stack_set(model,'info','Freq',2*logspace(0,5,40)');

if nargout==0; feplot(model);fecom colordatamat
else; out=model;
end

%% #MeshPiezoShaker: calibration of a piezo with a shaker
elseif comstr(Cam,'piezoshaker')

% Mesh
model=feutil('object quad 4 4',      [0 0   0;10 0 0;0 0 -10],8,2);
model=feutil('object quad 5 5',model,[0 0 -10;10 0 0;0 0  -2],8,2);
model=feutil('object quad 4 4',model,[0 0 -12;10 0 0;0 0 -10],8,2);
model=feutil('Rev 20 o 0 0 0 360 0 0 1',model);
model.Node(:,5:7)=model.Node(:,5:7)/1000; % SI
model.unit='SI';
model=feutil('AddTest-NoOri-merge',d_piezo('MeshBaseAccel'),model);
%model=feutil('AddTest-NoOri',d_piezo('MeshBaseAccel'),model);
%if you do no wish to merge

model.DOF=[]; % Do not reuse DOF field of BaseAccel model
model.name='Shaker-accelero assembly';

% Piezo properties of shaker - SONOX P502_iso
model.pl=[ model.pl; 
   m_piezo('dbval 5 -elas 6 SONOX_P502_iso')];
model=p_solid('default;',model);

%%% Electrodes and boundary conditions (already exist for accelerometer)
% Remove input for accelerometer
model=fe_case(model,'remove','IN');

% Shaker
  % -input "In" says it will be used as a voltage actuator
model=p_piezo('ElectrodeMPC Top Actuator -input "Vin-Shaker"',model,'z==-0.01');
  % -ground generates a v=0 FixDof case entry
model=p_piezo('ElectrodeMPC Bottom Actuator -ground',model,'z==-0.012');
  % fix shaker at bottom
model=fe_case(model,'fixdof','Fixed-shaker','z==-0.022');

if nargout==0; feplot(model);fecom colordatamat
else; out=model;
end

%% #MeshLangevin: Mesh a Langevin type transducer 

elseif comstr(Cam,'langevin')
% Mesh
model=feutil('object quad 1 1',      [0 0   0;25 0 0;0 0 19],8,5);
model=feutil('object quad 2 2',model,[0 0   19;25 0 0;0 0  2],8,2);
model=feutil('object quad 1 1',model,[0 0   21;25 0 0;0 0 19],8,5);
model=feutil('Rev 20 o 0 0 0 360 0 0 1',model);
model.Node(:,5:7)=model.Node(:,5:7)/1000; % SI
model.unit='SI';
model.name='Langevin-transducer';

% Piezo properties of transducer - SONOX P502
model.pl=[ m_elastic('dbval 1 steel') zeros(1,24); 
            m_piezo('dbval 2 -elas 3 SONOX_P502')];
model=p_solid('default;',model);

%%% Electrodes and boundary conditions (already exist for accelerometer)

  % -input "In" says it will be used as a voltage actuator
model=p_piezo('ElectrodeMPC Top Actuator -input "Vin"',model,'z==0.019');
  % -ground generates a v=0 FixDof case entry
model=p_piezo('ElectrodeMPC Bottom Actuator -ground',model,'z==0.021');
  % add a charge sensor on the top electrode
model=p_piezo(['ElectrodeSensQ '  ...
    num2str(floor(p_piezo('electrodedof Top Actuator',model))) ' Q'],model);

if nargout==0; feplot(model);fecom colordatamat
else; out=model;
end  

%% #MeshLangConcrete: Langevin Transducer in concrete

elseif comstr(Cam,'langconcrete')
% Mesh
model=feutil('object quad 1 1',      [0 0   0;25 0 0;0 0 19],8,5);
model=feutil('object quad 2 2',model,[0 0   19;25 0 0;0 0  2],8,2);
model=feutil('object quad 1 1',model,[0 0   21;25 0 0;0 0 19],8,5);
model=feutil('object quad 4 4',model,[0 0   40;25 0 0;0 0 100],8,20);
model=feutil('object quad 4 4',model,[0 0   -100;25 0 0;0 0 100],8,20);
% initial node + 2 vectors for direction starting at 0

model=feutil('object quad 4 4',model,[25 0   0;100 0 0;0 0 19],32,5);
model=feutil('object quad 4 4',model,[25 0   19;100 0 0;0 0  2],32,2);
model=feutil('object quad 4 4',model,[25 0   21;100 0 0;0 0 19],32,5);
model=feutil('object quad 4 4',model,[25 0   40;100 0 0;0 0 100],32,20);
model=feutil('object quad 4 4',model,[25 0   -100;100 0 0;0 0 100],32,20);

 model=feutil('join group 4:10',model);

model=feutil('Rev 20 o 0 0 0 360 0 0 1',model);
model.Node(:,5:7)=model.Node(:,5:7)/1000; % SI
model.unit='SI';
model.name='Langevin-transducer-concrete';

% Piezo properties of transducer - SONOX P502
model.pl=[ m_elastic('dbval 1 steel') zeros(1,24); 
            m_piezo('dbval 2 -elas 3 SONOX_P502') ;
%            m_elastic('dbval 4 air') zeros(1,26) ];
           m_elastic('dbval 4 concrete') zeros(1,23) ];
model=p_solid('default;',model);

%%% Electrodes and boundary conditions (already exist for accelerometer)

  % -input "In" says it will be used as a voltage actuator
model=p_piezo('ElectrodeMPC Top Actuator -input "Vin"',model,'z==0.019');
  % -ground generates a v=0 FixDof case entry
model=p_piezo('ElectrodeMPC Bottom Actuator -ground',model,'z==0.021');
  % add a charge sensor on the top electrode
model=p_piezo(['ElectrodeSensQ '  ...
    num2str(floor(p_piezo('electrodedof Top Actuator',model))) ' Q'],model);

if nargout==0; feplot(model);fecom colordatamat
else; out=model;
end  
    

    
%% #MeshIDEPatch : patch with inter digitated electrodes 
elseif comstr(Cam,'ide')
% Parameter handling
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   ' lx(400e-6#%g#"patch width")'...
   ' ly(300e-6#%g#"patch height")'...
   ' p(700e-6#%g#"patch length")'...
   ' e(50e-6#%g#"IDE width")'...
   ' nx(10#%g#"n_elt x direction")'...
   ' ny(5#%g#"n_elt y direction")'...
   ' nz(14#%g#"n_elt z direction")'...
   ' Quad(#3#"quadratic elements if present")'...
   ],{RO,CAM});
    
% Mesh
model=feutil('objectbeam 1 1',[0 0 0;RO.lx 0 0],RO.nx);
model=feutil(['extrude ' num2str(RO.ny) ' 0 ' num2str(RO.ly/RO.ny) ' 0'],model);
model=feutil(['extrude ' num2str(RO.nz) ' 0 0 ' num2str(RO.p/RO.nz) ],model);

% Define material properties
model.pl=m_piezo('dbval 1 -elas 10 SONOX_P502_iso');   
% Integration rules for volumes
model=p_solid('default;',model);
model.unit='SI';

% Build a MPC and defining ground and actuation DOF
InputDOF=[];
[model,InputDOF(end+1,1)]=p_piezo('ElectrodeMPC IDE bottom -ground',model, ...
    ['y==0 | y==' num2str(RO.ly) '& z<' num2str(RO.e*1.01)]);
[model,InputDOF(end+1,1)]=p_piezo('ElectrodeMPC IDE top -input"Applied Voltage"' ...
    ,model,['y==0 | y==' num2str(RO.ly) '& z>' num2str(RO.p-RO.e*1.01)]);
RO.InputDOF=InputDOF;
model=stack_set(model,'info','MeshInfo',RO);

if ~isfield(model,'unit'); model.unit='SI'; end
if RO.Quad;model=feutil('lin2quad',model);end

if nargout==0;feplot(model); fecom colordatapro-edgealpha.1
else;out=model;
end  % Send output to out variable
if nargout>1; out1=RO; end % Send MeshInfo back    


elseif comstr(Cam,'ulbplate')
%% #MeshULBPlate : sample mesh of a plate with four patches - - -2
%  Note that d_piezo('MeshPlate') is much more general

%% generate the model - - - - - - - - - - - - - - - - - - - - - - - -
 model=struct('Node',[1 0 0 0 0 0 0],'Elt',[]);
 model=feutil('addelt',model,'mass1',1);
 % Note that the extrusion values are chosen to include the patch edges
 dx=[linspace(0,15,3) linspace(15,15+55,10) linspace(15+55,463-5,50) 463];
 model=feutil('extrude 0 1 0 0',model,.001*unique(dx));
 dy=[linspace(0,12,3) linspace(12,12+25,5) linspace(12+25,63,5) ...
    linspace(63,63+25,5) linspace(63+25,100,3)];
 model=feutil('extrude 0 0 1 0',model,.001*unique(dy));

 %Set patch areas as set different properties
 model.Elt=feutil('divide group 1 withnode{x>.015 & x<.07 & y>.013 & y<.037 }',model);
 model.Elt=feutil('divide group 2 withnode{x>.015 & x<.07 & y>.064 & y<.088 }',model);
 model.Elt=feutil('set group 1 proId3',model);
 model.Elt=feutil('set group 2 proId4',model);

%%%%%% Material Properties
model.pl=m_elastic('dbval 1 ', ...
     [1 fe_mat('m_elastic','SI',1) 42.5e9 .042 1490 3.35e9 .01], ... %composite layer property
     'dbval 2 ', ...
     [1 fe_mat('m_elastic','SI',1) 65e9 .3 7800 25e9 .01]); % PZT elastic properties 
model.pl=m_piezo(model.pl,'dbval 3 -elas2 Sample_ULB');
model.pl(model.pl(:,1)==3,4:5)=-model.pl(3,4:5); %polarization change 

%%%%% Laminate properties

model.il=p_shell('dbval 1 laminate 1 2.1666667e-4 0 1 2.1666667e-4 90 1 2.1666667e-4 0 1 2.1666667e-4 90 1 2.1666667e-4 0 1 2.1666667e-4 90', ...
'dbval 2 laminate 3 2.5e-4 0 1 2.1666667e-4 0 1 2.1666667e-4 90 1 2.1666667e-4 0 1 2.1666667e-4 90 1 2.1666667e-4 0 1 2.1666667e-4 90 3 2.5e-4 0');

%%% Piezo electrodes
%%%                                         NdNb LayerID   NdNb  LayerID
model.il=p_piezo(model.il,'dbval 3 shell 2 1682    1   0    1683  8 0');
model.il=p_piezo(model.il,'dbval 4 shell 2 1684    1   0    1685  8 0');

model.name='ULB_plate'; model.unit='SI';
if ~isempty(strfind(Cam,'cantilever'));
    model=fe_case(model,'FixDof','Cantilever','x==0 -DOF 1:6');
end

if ~isempty(strfind(Cam,'feplot'))||nargout==0; 
    cf=feplot(model); fecom colordatapro-edgealpha.1
    if nargout>0;out=cf;end
else;out=model;
end

elseif comstr(Cam,'mfcplate')
%% #MeshMFCPlate : sample mesh of a plate with two MFCs
%  Note that d_piezo('MeshPlate') is much more general

%% generate the model - - - - - - - - - - - - - - - - - - - - - - - -
mfcmdl=m_piezo('patchSmartM.MFC-P1.8528'); % 
%xxxad : is thickness always .2 ?

if 1==2 % Will need to go to new geometry formatting
 RG.list={'Name','Lam','shape'
   'Main_plate', mdl,'global lx=.4 ly=.3 lc=.02'
   'Act1','BaseId1+SmartM.MFC-P1.2814 -in','xc=.03 yc=.05 ang=30'}
end

% Create mesh
model=struct('Node',[1 0 0 0 0 0 0],'Elt',[]);
model=feutil('addelt',model,'mass1',1); 
% Note that the extrusion values are chosen to include the patch edges
dx=feutil('refineline 3',[0 15 15+mfcmdl.lx 463]);
dy=feutil('refineline 3',[0 11 11+mfcmdl.ly 50]);
model=feutil('extrude 0 1 0 0',model,dx);
model=feutil('extrude 0 0 1 0',model,dy);
model.Node(:,5:7)=model.Node(:,5:7)/1000; model.unit='SI'; % go to m

%Set patch area as set different properties
 model.Elt=feutil('divide group 1 withnode{x>.0151 & x<.0999 & y>.0111 & y<.0389 }',model);
 % xxxEB cleaner call to divide the groups ?
 model.Elt=feutil('set group 1 proId3',model);

 %%%%%% Material Properties for supporting plate
 model.pl=m_elastic('dbval 1 Aluminum'); % Aluminum
 model.il=p_shell('dbval 1 laminate 1 1e-3 0');
 
 % Merge properties with the supporting plate model
 model=feutil('setmat',model,mfcmdl.pl);
 % Add piezo laminate and then merge with base plate
 model=feutil('setpro',model, ...
   p_shell('dbval 2 Merge 1 +100 -100',feutil('setpro',model,mfcmdl.il)));
 % now define electrode
 model.il=p_piezo(model.il,'dbval 3 shell 2 1682    3   0    1683   9   0');
 p_piezo('tabpro',model);
 
if ~isempty(strfind(Cam,'cantilever'));
    model=fe_case(model,'FixDof','Cantilever','x==0 -DOF 1:6');
end

if ~isempty(strfind(Cam,'feplot'))||nargout==0; 
    cf=feplot(model); fecom colordatapro-edgealpha.1
    if nargout>0;out=cf;end
else;out=model;
end

elseif comstr(Cam,'triangleplate')
%% #MeshTrianglePlate : Mesh plate with triangular point load actuator

fname=fullfile(sdtdef('tempdir'),'demo/PzTrianglePlate.mat');
if exist(fname,'file')&&isempty(strfind(Cam,'reset'))
  fprintf('Loading %s\n',fname);
  load(fname); out=model;
else
%% Geometrical properties and reference nodes
L=414; w=314; a=100; b=33.58; c=(w-b)/2; 

% create reference nodes
mdl.Node = [1 0 0 0  0 0 0; 
          2 0 0 0  L 0 0; 
          3 0 0 0  L w 0;
          4 0 0 0  0 w 0 ;
          5 0 0 0  0 b+c 0;
          6 0 0 0  0 c  0;
          7 0 0 0  a w/2 0];
mdl.Elt=[];

%% Meshing - Need gmsh

% create reference lines and mesh plate
mdl=fe_gmsh('AddLine -loop1',mdl,[1 2; 2 3; 3 4; 4 5; 5 7; 7 6; 6 1]);
mdl=fe_gmsh('AddLine -loop2',mdl,[5 7; 7 6; 6 5]);
mdl.Stack{end}.PlaneSurface=[1; 2]; mdl.Stack{1,3}.LineLoop{2,2}=[5 6 10];
model=fe_gmsh('write tmp.geo -lc 20 -run -2 -v 0',mdl)
delete('tmp.msh');delete('tmp.geo');

% divide in groups and set ProIds
model.Elt=feutil('sel elt group 3',model);
model.Elt=feutil('divide group 1 Proid 14',model);
model.Elt=feutil('set group 1 proId 3',model);
model.Elt=feutil('set group 2 proId 1',model);
model.Elt=feutil('set group 1:2 matid 1',model);
model.Node(:,5:7)=model.Node(:,5:7)/1000;
model.Node=[model.Node; 100001 0 0 0 1 1 1]; % Add node for piezo electrode

%% Material and section properties

% Material properties
model.pl=[1 fe_mat('m_elastic','SI',1) 69e9 0.31 2700]; % Aluminum
model.pl=[model.pl zeros(1,25); 
m_piezo('dbval 101 -elas 102 SONOX_P502');
[103 fe_mat('m_elastic','SI',1) 2.6e9 .33 1500] zeros(1,25) ; % epoxy
[104 fe_mat('m_elastic','SI',1) 1 0 1] zeros(1,25)]; % dummy layer
 
 % ! Need to change orientation for SONOX_P502 : Poling in direction x
 % instead of z
 model.pl(2,[6 12])=model.pl(2,[12 6]); % d coefficients change
 model.pl(3,[3 5 8 ])=[model.pl(3,[5 3 ]) 0.4]; % E coeff change + adjust nu
 
% Section properties
model.il=p_shell('dbval 1 laminate 1 1e-3 0', ...
    ['dbval 2 laminate  104 3e-04 0 1 1e-3  0 103 6e-05 0 101 1.8e-4 0 103 6e-05 0']);

model.il=p_piezo(model.il,'dbval 3 shell 2 100001    4   0');
model.name='SSS_example'; model.unit='SI';

% orient all triangles
data=struct('sel','groupall','dir',{{'1','0',0}},'DOF',[.01;.02;.03]);
       model=p_shell('setTheta',model,data);

       % Boundary conditions       
model=fe_case(model,'FixDof','Cantilever1','x==0 | x==0.414 | y==0 | y==0.314  -DOF 1:6');

if ~isempty(strfind(Cam,'feplot'))||nargout==0
  cf=feplot(model); fecom colordatapro-edgealpha.1
  if nargout>0;out=cf;end
  if nargout==0
  fprintf('Saving %s\n',fname); sdtkey(mkdir(fileparts(fname)));
  save(fname,'model');
  end
else;out=model;
end
end

elseif comstr(Cam,'shunt')
%% #MeshShunt : Mesh plate for shunt example

% Create mesh
model=struct('Node',[1 0 0 0 0 0 0],'Elt',[]);
model=feutil('addelt',model,'mass1',1); 
% Note that the extrusion values are chosen to include the patch edges
dx=feutil('refineline 3',[0 50 100 300 350]);
dy=feutil('refineline 3',[0 25]);
model=feutil('extrude 0 1 0 0',model,dx);
model=feutil('extrude 0 0 1 0',model,dy);
model.Node(:,5:7)=model.Node(:,5:7)/1000; model.unit='SI'; % go to m

% Divide in groups
  model.Elt=feutil('divide group 1 withnode{x<.0499 }',model);
  model.Elt=feutil('divide group 2 withnode{x<0.0999 }',model);
  model.Elt=feutil('set group 1 proid 101',model);
  model.Elt=feutil('set group 2 proid 201',model);
  model.Elt=feutil('set group 3 proid 1',model);

% Define material properties
model.pl=m_elastic('dbval 1 steel');
model.pl=m_piezo(model.pl,'dbval 2 -elas 3 PIC_255');

% Define section properties
model.il=p_shell('dbval 1 laminate 1 2e-3 0');
model.il=p_shell(model.il,'dbval 2 laminate 2 5e-4 0 1 2e-3 0 2 5e-4 0');
model.il=p_shell(model.il,'dbval 3 laminate 2 5e-4 0 1 2e-3 0 2 5e-4 0');
model.il=p_piezo(model.il,'dbval 101 shell 2 10001 1 0 10002 3 0');% Electrodes of actuators
model.il=p_piezo(model.il,'dbval 201 shell 3 20001 1 0 20002 3 0');% Electrodes of sensors

% BC
 model=fe_case(model,'FixDof','Cantilever','x==0 -DOF 1:6')
 
 if ~isempty(strfind(Cam,'feplot'))||nargout==0; 
    cf=feplot(model); fecom colordatapro-edgealpha.1
    if nargout>0;out=cf;end
else;out=model;
end
 

elseif comstr(Cam,'simplecomp')
%% #MeshSimpleComp : simple volume representative of a MFC layout

%%% define geometry
opt=struct( ...
  'w',  100e-6, ...    % width of cell
  'h',  100e-6, ...    % height of cell
  'd',  100e-6, ...    % depth of cell
  'n',  2, ...         % cell subdivision along x
  'ncy',5, ...        % number of unit cells along y
  'ncx',4, ...;        % number of unit cells along x
  'nz',4);             % model subdivision along z


%%%% build generic cell
model=feutil('objectbeam 1 1',[0 0 0;opt.w*opt.ncx 0 0],2*opt.ncx);
model.Elt([3:4:end 4:4:end],3:4)=2; % layered properties 1 and 2
model=feutil(sprintf('divide %i',opt.n),model);
model=feutil(sprintf('extrude %i 0 %.15g 0',opt.ncy,opt.d),model);
model=feutil(sprintf('extrude %i 0 0 %.15g',opt.nz,opt.h/opt.nz),model);

% Use separate MatId for top and bottom piezo to avoid internal electrode
model.Elt(feutil('findelt matid 2 & innode {z<=}',model,opt.h/2),9:10)=3;
model.Elt=feutilb('SeparatebyProp',model.Elt);

%%%% Define material properties
model.pl=[100 fe_mat('m_elastic','SI',1) 66e9 .3 7800 0 ;  % Elastic for PZT
          1 fe_mat('m_elastic','SI',1) 2.6e9 .3 1500 0]; % Epoxy
% Piezo properties for the fibers
%        [ProId Type ElasMatId d31 d32 d33 epsi1 epsi2 epsi3  EDType] 
eps0=8.854e-12;
pl=[2 fe_mat('m_piezo','SI',1) 100 ...
  -185e-12 -185e-12 440e-12 1850*eps0 1850*eps0 1850*eps0 0 ];
model=feutil('setmat',model,pl); % Piezo top
pl(1)=3; model=feutil('setmat',model,pl); % Piezo top

% Integration rules for volumes
model=p_solid('default;',model);
model.name='MFC representative volume';model.unit='SI';

% The are 3 electrodes : bottom, middle and top of piezo volumes
% which are polarized along z

InputDOF=[];
% Build a MPC defining a single potential for the electrodes
[model,InputDOF(end+1,1)]=p_piezo('ElectrodeMPC Middle',model,'z==5e-5');
[model,InputDOF(end+1,1)]=p_piezo('ElectrodeMPC Top',model,'z==1e-4');
[model,InputDOF(end+1,1)]=p_piezo('ElectrodeMPC Bottom',model,'z==0');

% Clamp one edge
model=fe_case(model,'FixDof','Clamp','y==0 -DOF 1 2 3');
% set middle electrode to potential V=0
model=fe_case(model,'FixDof','V=0 on middle',InputDOF(1));

if ~isempty(strfind(Cam,'feplot'))||nargout==0; out=feplot(model);
    fecom colordatapro-edgealpha.1
else;out=model;
end

elseif comstr(Cam,'sampleacoustic')
%% #MeshSampleAcoustic : simple example of an acoustic generator

mdl=struct('Node',[1 0 0 0  0 0 0],'Elt',feutil('addelt','mass1',1));

lz={linspace(0,2,3) % backing damper, 2mm, 3 mesh points
    linspace(0,1,4) % Piezo, 1 mm 4 mesh points
    linspace(0,.1,4) % Front plate, .1mm, 4 mesh ponts
    linspace(0,5,10) % Water
    };

% Mesh layers
r1=lz{1}; RO.zBottom=r1(end);
for j1=2:length(lz); % Build z values
    r1=[r1 lz{j1}+r1(end)]; %#ok<AGROW> 
    if j1==2;RO.zTop=r1(end);end
end 
mdl=feutil('extrude 0 0 0 1',mdl,unique(r1));
% set material properties
i0=1;
for j1=1:length(lz); 
 i0=i0(end)+(1:length(lz{j1})-1);
 mdl.Elt(i0,3:4)=j1; % set MatId ProId
end

% Generate a sector
mdl=feutil('extrude 1 1 0 0',mdl);
mdl=feutil('rev 5 o 0 0 0 90 0 0 1',mdl);
mdl.Node(:,5:7)=mdl.Node(:,5:7)/1000; mdl.unit='SI';
RO.zBottom=RO.zBottom/1000;RO.zTop=RO.zTop/1000;

% display geometry
% feplot(mdl);fecom('colordatamat-edgealpha.1')

% Define material properties
mdl.pl=m_elastic('dbval 1 Epoxy_Resin','dbval 3 steel','dbval 4 Water');
mdl=feutil('setmat 1 eta=.1',mdl); % Set loss factor

mdl.pl=m_piezo(mdl.pl,'dbval 2 -elas 102 PZT_Ansys');
mdl=p_solid('default;',mdl);
mdl.Elt=feutilb('SeparatebyProp',mdl.Elt);

% Define piezo electrodes
mdl=p_piezo('ElectrodeMPC Bottom-ground',mdl,sprintf('z==%.15g',RO.zBottom));
mdl=p_piezo('ElectrodeMPC Top-Input"In"',mdl,sprintf('z==%.15g',RO.zTop));

% Symmetric structural response
mdl=fe_case(mdl,'FixDof','EdgeX','y==0 -DOF 2');
mdl=fe_case(mdl,'FixDof','EdgeY','x==0 -DOF 1');

mdl=stack_set(mdl,'info','Freq',logspace(log10(20e3),log10(30e3),10)');
def=fe_simul('dfrf',mdl);

cf=feplot(mdl,def); 
cf.sel={'matId~=4','colordataEvalZ -edgealpha.1'}; % Mechanical response

elseif comstr(Cam,'gammas')
%% #MeshGammaS : build a weigting for surface control
% RO.Rect={[xc lx nx, yc ly ny alpha],'gauss'} 
if isfield(RO,'Elt');model=RO;RO=varargin{carg};carg=carg+1;
else; model=varargin{carg};carg=carg+1;
end
RO.lc=0;
if carg<=nargin; def=varargin{carg};carg=carg+1;else; def=[]; end

for j1=1:size(RO.Rect,1); % For now single rectangle
 r2=RO.Rect{j1,1}; r3=r2(7)*pi/180;
 RO.x=[cos(r3) sin(r3) 0]'; RO.y=[-sin(r3) cos(r3) 0]'; % Local patch orient
 n1=[r2([1 4]) 0;r2(2)*RO.x(:)';r2(5)*RO.y(:)'];
 r2=model.Node(:,5:7)-n1(ones(size(model.Node,1),1),:);
 r3=r2*RO.y; r2=r2*RO.x; 
 i1=(r2>=-RO.lc & r2<=RO.Rect{j1,1}(2)+RO.lc ...
     & r3>=-RO.lc & r3 <=RO.Rect{j1,1}(j1,5)+RO.lc); % Selection rectangle
 mo2=model;mo2.Elt=feutil('SelElt withnode',model,model.Node(i1,1));
 mo2.Elt=feutil('setgroupall matid 1 proid 1',mo2);
 mo2.pl=[1 fe_mat('m_elastic','SI',2) 1 1 1];
 mo2.il=[1 fe_mat('p_solid','SI',2) 0 0 -3];
 mo2.Stack={};
 %mo2=fe_case(mo2,'assemble -SE -matDes 2 -NoT'); % uniform
 cut=fe_caseg('StressCut-SelOut',struct('type','Gauss'),mo2);
 r2=cut.StressObs;
 %r2=fe_homo('AssGaussObs',mo2);
 n2=r2.Node(:,5:7)*[RO.x(:) RO.y(:)];n2=n2-n1(ones(size(n2,1),1),1:2);
 n2=n2*diag(1./[n1(2,1);n1(3,2)])*2-1; % Normalized -1:1 rectangle
 x=n2(:,1);y=n2(:,2);
 if strcmpi(RO.Rect{j1,2},'flat'); w=ones(size(r2.wjdet));
 elseif strcmpi(RO.Rect{j1,2},'gauss');
   w=exp(-(x.^2+y.^2));
 end 
 %
 m=r2.trans.cta'*diag(sparse(r2.wjdet.*w))*r2.trans.cta;
 out=struct('K',{{m}},'Klab',{{'Weight'}},'DOF',fix(r2.DOF)+.03); 
 break; % single patch
end
if ~isempty(def)
  def=feutilb('placeindof',out.DOF,def);
  out.K=feutilb('tkt',def.def,out.K); out.DOF=(1:size(def.def,2))'+.99;
  out.Klab{1}='GGt';
end
if nargout==0 % Display for verification
    cf=feplot;cf.def={diag(out.K{1}),out.DOF};fecom colordataz
end

elseif comstr(Cam,'homomfcp2')
%% #MeshHomoMfcP2 : mesh a representative volume for P2-type MFC

[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   ' lx(300e-3#%g#"patch length")'...
   ' ly(300e-3#%g#"patch width")'...
   ' lz(180e-3#%g#"patch height")'...
   ' dd(4e-2#%g#"patch height")'...
   ' rho(#%g#"Volume fraction of piezo")'...
   ],{RO,CAM});

% Create mesh
model=struct('Node',[1 0 0 0 0 0 0],'Elt',[]);
model=feutil('addelt',model,'mass1',1);

model.dx=feutil(['refineline ' num2str(RO.dd)],[0 RO.lx/2*(1-RO.rho) RO.lx/2*(1+RO.rho) RO.lx]);
model.dy=feutil(['refineline ' num2str(RO.dd)],[0 RO.ly]);
model.dz=feutil(['refineline ' num2str(RO.dd)],[0 RO.lz]);
model=feutil('extrude 0 1 0 0',model,model.dx);
model=feutil('extrude 0 0 1 0',model,model.dy);
model=feutil('extrude 0 0 0 1',model,model.dz);
model.unit='SI'; % go to m

% Divide in groups
  model.Elt=feutil(['divide group 1 withnode{x<' num2str(RO.lx/2*(1-RO.rho)-1e-4) '}'],model);
  model.Elt=feutil(['divide group 2 withnode{x<' num2str(RO.lx/2*(1+RO.rho)-1e-4) '}'],model);
  model.Elt=feutil('set group 1 matid 1',model);
  model.Elt=feutil('set group 2 matid 2',model);
  model.Elt=feutil('set group 3 matid 1',model);

% Define material properties

if 1==1
% Piezo material with no coupling (dij=0) and epoxy-resin properties
 model.pl=m_piezo('dbval 1 -elas 4 SONOX_P502_iso');
 model.pl(2,3:5)=[2.9e9 0.3 1200];
 model.pl(1,[6 9 12 14 16])=0; model.pl(1,[22 26 30])=4.25*8.854e-12;
else
 model.pl=m_elastic('dbval 1 Epoxy_Resin');
end

model.pl=m_piezo(model.pl,'dbval 2 -elas3 SONOX_P502');
% correct for values in IJMSS paper
model.pl(4,6:11)=[0.44 0.393 0.41 19.48e9 19.48e9 19.14e9];
model.il=p_solid('dbval 1 d3 -3');

% Define electrodes
  % -input "In" says it will be used as a voltage actuator
model=p_piezo('ElectrodeMPC Top Actuator -input "Vin"',model,['z==' num2str(RO.lz)]);
  % -ground generates a v=0 FixDof case entry
model=p_piezo('ElectrodeMPC Bottom Actuator -ground',model,'z==0');
out=model;

elseif comstr(Cam,'homomfcp1')
%% #MeshHomoMfcP1 : mesh a representative volume for P1-type MFC

[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   ' lx(180e-3#%g#"patch length")'...
   ' ly(180e-3#%g#"patch width")'...
   ' lz(1080e-3#%g#"patch height")'...
   ' e(90e-3#%g#"ide size")'...
   ' dd(4e-2#%g#"patch height")'...
   ' rho(#%g#"Volume fraction of piezo")'...
   ],{RO,CAM});

    % Create mesh
model=struct('Node',[1 0 0 0 0 0 0],'Elt',[]);
model=feutil('addelt',model,'mass1',1);

model.dx=feutil(['refineline ' num2str(RO.dd)],[0 RO.lx/2*(1-RO.rho) RO.lx/2*(1+RO.rho) RO.lx]);
model.dy=feutil(['refineline ' num2str(RO.dd)],[0 RO.e RO.ly-RO.e RO.ly]);
model.dz=feutil(['refineline ' num2str(RO.dd)],[0 RO.lz]);
model=feutil('extrude 0 1 0 0',model,model.dx);
model=feutil('extrude 0 0 1 0',model,model.dy);
model=feutil('extrude 0 0 0 1',model,model.dz);
model.unit='SI'; % go to m

% Divide in groups
  model.Elt=feutil(['divide group 1 withnode{x<' num2str(RO.lx/2*(1-RO.rho)-1e-5) '}'],model);
  model.Elt=feutil(['divide group 2 withnode{x<' num2str(RO.lx/2*(1+RO.rho)-1e-5) '}'],model);
  model.Elt=feutil('set group 1 matid 1',model);
  model.Elt=feutil('set group 2 matid 2',model);
  model.Elt=feutil('set group 2 Proid 2',model);
  model.Elt=feutil('set group 3 matid 1',model);

% Define material properties

% Piezo material with no coupling (dij=0) and epoxy-resin properties
 model.pl=m_piezo('dbval 1 -elas 4 SONOX_P502_iso');
 model.pl(2,3:5)=[2.9e9 0.3 1200];
 model.pl(1,[6 9 12 14 16])=0; model.pl(1,[22 26 30])=4.25*8.854e-12;

 model.pl=m_piezo(model.pl,'dbval 2 -elas3 SONOX_P502');
 % correct for values in IJMSS paper
 model.pl(4,6:11)=[0.44 0.393 0.41 19.48e9 19.48e9 19.14e9];
 model.il=p_solid('dbval 1 d3 -3','dbval 2 d3-3');


%Define electrodes
model=p_piezo('ElectrodeMPC Bottom Actuator -ground',model, ...
    ['y==0 | y==' num2str(RO.ly) '& z<' num2str(RO.e*1.01)]);
model=p_piezo('ElectrodeMPC Top Actuator -input"Vin"' ...
    ,model,['y==0 | y==' num2str(RO.ly) '& z>' num2str(RO.lz-RO.e*1.01)]);

out=model;

elseif comstr(Cam,'infsphere')
%% #MeshInfSphere : mesh a spherical piezo with infinite sphere around it
Nb = 4;
f1 = 1e3;
fend = 5e5;
freq = linspace(f1,fend,Nb);

%Geometry
Rin = 0.005;            %Internal radius
r4 = Rin/sqrt(3);       %Node position to create volley ball 
r0 = 4.1*Rin;             %External Cube (end PML)
r1 = 4*r4;             %Internal Cube (end of computation region)
r3 = 2.5*r4;            %Radius Layer 1

%Transfinite Mesh properties
% n0 = 3;                 %Slice PML
% n1 = 11;                 %Slice External region 
% n2 = 13;                 %Layer 1 to cube
% n3 = 13;                 %Interal Radius to Layer 1 
% progr = 1;
% 
n0 = 2;                 %Slice PML
n1 = 13;                 %Slice External region 
n2 = 13;                 %Layer 1 to cube
n3 = 13;                 %Interal Radius to Layer 1 
progr = 1;

fileID = fopen('temp.geo','w');
fprintf(fileID,'r0 = %f ; \r\n r1 = %f ; \r\n r3 = %f ; \r\n r4 = %f ; \r\n',[r0 r1 r3 r4]);
fprintf(fileID,'n0 = %f ; \r\n n1 = %f ; \r\n n2 = %f ; \r\n n3 = %f ; \r\n ',[n0 n1 n2 n3]);
fclose(fileID);


%%%%%%%%%%
%% MESH %%
%%%%%%%%%%

%model=fe_gmsh('write -run -3 -v 0','spherevdb.geo');
%model=fe_gmsh('write -run ','spherevdb.geo');
model=fe_gmsh('write spherevdb.msh -run','spherevdb.geo')

[model.Elt,RemovedElt]=feutil('RemoveElt EltName mass2',model);
[model.Elt,RemovedElt]=feutil('RemoveElt EltName beam1',model);
[model.Elt,RemovedElt]=feutil('RemoveElt EltName quad4',model);

cf=feplot(model);
com1 =['groupall & innode{x>=0 & y>=0}'];
cf.sel(1) = com1;
% fecom('colordatamat')
% % % pause
% comgui('closeall');


% model=feutil('Lin2Quad',model);
namenum = model.Elt(1,2:end);
namenum = namenum(namenum~=0);
posmat = eval([ char(namenum) '(''prop'') ']);
model.Elt(2:end,posmat(3))=model.Elt(2:end,posmat(1));
model.Elt(2:end,posmat(1))=1;
model.Elt(2:end,posmat(2))=1;


%% MATERIAL DEFINITION
Eext = 6e9;
nuext = 0.2;
rhoext = 900; 

Vs = sqrt(Eext/(2*(1+nuext)*rhoext)); 
Vp = sqrt((1-nuext)*Eext/((1-2*nuext)*(1+nuext)*rhoext)); 
c1 = rhoext*Vp; 
c2 = rhoext*Vs; 

model.pl=[];
model.pl(end+1,1:9)= [1  fe_mat('m_elastic','SI',1) Eext nuext rhoext 0 0.01 0 0];


%% PROPERTIES
model=p_solid('default',model);
model.unit='SI';



dbstack; keyboard

%% #MeshEnd
else;error('Mesh%s unknown',CAM);
end
%% #Solve -------------------------------------------------------------------
elseif comstr(Cam,'solve');[CAM,Cam]=comstr(CAM,6);

if comstr(Cam,'reducebase+f1')
%% #Base+f1 Enforce base acceleration with correction for first frequency -2

eval(iigui({'model','Case','Load','out'},'MoveFromCaller'));
eigopt=fe_def('DefEigOpt',model);

% robustness should have a get_k, get_f
if ~isequal(model.Opt(2,:),[2 1])
 error('Implementation assumes M and K')
end
def=fe_eig({model.K{1},model.K{2},Case.T,model.DOF},eigopt);

freq=stack_get(model,'info','Freq','get'); 
if isempty(freq);freq=def.data(1)/10;end
Z=feutilb('sumkcoef',model.K,[-(freq(1)*2*pi)^2 1]);
% bset=feutilb('placeindof',model.DOF,Load.bset); not correct
bset=struct('def',Case.TIn,'DOF',model.DOF);
Zin=-(Z*bset.def);
% Remove inertial load of moving frame
%Zin=Zin-model.K{1}*def.def*((Zin.'*def.def).');
u=ofact({Z,Case.T},Zin); % +bset.def; % Reference response in global frame
u=u-def.def*(u.'*model.K{1}*def.def)'; % orthogonalize with respect to modes

[T,fr]=fe_norm([def.def u],model.K{1:2});
out.KIn=cell(1,length(model.K));
for j1=1:length(out.KIn)
 out.KIn{j1}=((model.K{j1}*bset.def)'*T)'; 
end
    
out.K=feutilb('tkt',T,model.K); out.TR.def=T;out.TR.data=fr/2/pi;
out.DOF=(1:size(T,2))'+.99;

assignin('caller','T',out); clear out;

elseif comstr(Cam,'patchcapacity')
%% #SolvePatchCapacity : reduction for capacity resolution

model=[];Sens=[];eval(iigui({'model','Case','Load','out','Sens','RunOpt'},'MoveFromCaller'));
m=feval(fe_reduc('@get_mat'),model.K,2,model,RunOpt);
k=feval(fe_reduc('@get_mat'),model.K,1,model,RunOpt);

Freq=stack_get(model,'info','Freq','get');
if ~isempty(Freq)
  %% #FrfSnapShot % Called from fe_reduc
  MVR=struct('T',Case.T,'TIn',Case.TIn,'full',1, ...
    'K',{feutilb('tkt',Case.T,model.K)},'Klab',{model.Klab}, ...
   'BIN',[],'w',Freq(:)*2*pi,'DOF',model.DOF);
  MVR.BIN=cell(size(model.K));
  for j1=1:length(model.K);MVR.BIN{j1}=(-(model.K{j1}*MVR.TIn).'*MVR.T)';end
  
  d1=feval(fe_simul('@DfrfBIn'),MVR);
  % Case.TIn'*model.K{model.Opt(2,:)==1}*d1.def
  C1=struct('X',{{d1.data,Sens.lab}},'Xlab',{{'Freq','Sens'}}, ...
      'Y',(Sens.cta*d1.def).');
  
 % Estimate of grounded modes
 T=d1.def(:,2:end)-repmat(d1.def(:,1),1,size(d1.def,2)-1);
 % cf.def={T,model.DOF};fecom scalecolorone
 [T,fr]=fe_norm(T,model.K{1:2});fr=fr/2/pi;
 
 % Recompute at resonances
 if isfield(model,'Iterate')
  MVR.w=fr(fr>min(Freq)&fr<max(Freq))*2*pi; d3=feval(fe_simul('@DfrfBIn'),MVR);
  T=[d3.def(:,2:end)-repmat(d3.def(:,1),1,size(d3.def,2)-1) T];
  % cf.def={T,model.DOF};fecom scalecolorone
  [T,fr]=fe_norm(T,model.K{1:2});fr=fr/2/pi;
  'xxx missing .data(:,2) for capa influence'
 else; d3=[]; 
 end
 
 d2=struct('def',T,'DOF',model.DOF,'data',fr);
 if 1==2
  fecom colordata21-colorbartitle"V";fecom scalecolorone
  cf.SelF{1}.vert0= cf.SelF{1}.vert0*diag([1 1 10]);cla;feplot %#ok<NODEF>
 end
 
  MVR=struct('T',[],'TIn',d1.def(:,1),'full',0, ...
    'K',{feutilb('tkt',d2.def,model.K)}, ...
   'BIN',[],'w',Freq(:)*2*pi,'d3',d3, ...
   'TR',struct('def',[d2.def d1.def(:,1)],'DOF',model.DOF, ...
    'adof',(1:size(d2.def,2)+1)'+.99,'data',[d2.data(:,1);Inf]));
  MVR.BIN=cell(size(model.K));
  for j1=1:length(model.K);MVR.BIN{j1}=(-(model.K{j1}*MVR.TIn).'*d2.def)';end
  Sens.cta=Sens.cta*MVR.TR.def;Sens.DOF=MVR.TR.adof;
  MVR.Sens=Sens;
  MVR.C1=C1; % Full observation 
  
  clear out; assignin('caller','T',MVR); % fe_reduc output is T
  return
  %d1r=feval(fe_simul('@DfrfBIn'),MVR);

 
end 

eigopt=fe_def('DefEigOpt',model);
%def=fe_eig({m,k,Case.T,model.DOF},eigopt); % Free modes, grounded electrodes


bsr=struct('def',d1.def(:,1),'DOF',model.DOF);
bsr.K={}; for j1=1:length(model.K);bsr.K{j1}=(bsr.def'*model.K{j1}*d2.def).';end

% Reduced model solve with BIn
MVR=model; MVR.K=feutilb('tkt',d2.def,MVR.K);
fr=linspace(Freq(1),Freq(end),1000)';
for j1=1:length(fr)
  coef=[-(fr(j1)*2*pi)^2 1 0];
  Z=feutilb('sumkcoef',MVR.K,coef);
  Zin=feutilb('sumkcoef',bsr.K,coef); 
  u=d2.def*(Z\Zin)+bsr.def;
  dbstack; keyboard
end

 

T=[d1.def(:,1) d2.def];


% do a lanczos to get accurate estimate of in plane resonance
T=[];
freq=def.data(:,1);freq=mean(freq(1:min(7,length(freq))));
ite=1;
while ite<4
 Z=feutilb('sumkcoef',model.K,[-(freq(1)*2*pi)^2 1]);
 bset=struct('def',Case.TIn,'DOF',model.DOF);
 Zin=-(Z*bset.def); u=ofact({Z,Case.T},Zin);
 uf=u+bset.def; % Full u with input
 [T,fr]=fe_norm([uf T],model.K{1:2});fr=fr/2/pi;
 if abs(freq/fr(1)-1)<.01; break; else; freq=fr(1);end 
end

din=struct('def',T,'DOF',model.DOF,'data',fr);
dbstack; keyboard
% cf=feplot;cf.def={uf,model.DOF}; fecom colordata21


[T,fr]=fe_norm([uf def.def],model.K{1:2});
out.KIn=cell(1,length(model.K));
for j1=1:length(out.KIn)
 out.KIn{j1}=((model.K{j1}*bset.def)'*T)'; 
end
    
out.K=feutilb('tkt',T,model.K); out.TR.def=T;out.TR.data=fr/2/pi;
out.DOF=(1:size(T,2))'+.99;

assignin('caller','T',out); clear out;
assignin('caller','PostFcn','assignin(''caller'',''Load'',Load)'); 

elseif comstr(Cam,'patchdfrf')
%% #SolvePatchDfrf : capacity of an isolated patch

mo1=RO;carg=3;

% Compute in-plane modes to get frequency range
Freq=stack_get(mo1,'info','Freq','get');def=[];
if isempty(Freq);
 [SE,CE]=fe_case(fe_case(mo1,'fixdof','plane',.03), ...
     'assemble -matdes 2 1 -SE NoT');
 ie=fe_c(SE.DOF,.21,'ind'); is=setdiff(1:length(SE.DOF),ie);
 r1=diag(SE.K{2});r1(ie)=sqrt(mean(r1(is))/mean(abs(r1(ie))));r1(is)=1;
 r1=diag(sparse(r1));
 SE.K=feutilb('tkt',r1,SE.K);
 def=fe_eig({SE.K{1:2},CE.T,SE.DOF},[5 20 1e3]);
 def=fe_def('subdef',def,def.data(:,1)>1);
end
d1=[];
%feplot(mo1,def);fecom colordataevalA

% DFRF Compute forced frequency response 
% Define frequencies in Hz
 mo1=fe_case('sensmatch',mo1);
 if ~isempty(strfind(Cam,'snap'))
  % New strategy based on enforced potential
  % TZT qr + T Z TIn u = 0, in SS on must handle sens.cta*(T qr + TIn u)
  st='Call d_piezo(''SolvePatchCapacity'')-se-needsens-matdes2 1 3';
  Freq=stack_get(mo1,'info','RedFreq','get');
  if ~isempty(Freq)
   SE=stack_set(mo1,'info','Freq',Freq);
   if ~isempty(strfind(Cam,'iterate'));SE.Iterate=1;end
   SE=fe_reduc(st,SE);
  else;   SE=fe_reduc(st,mo1);
  end
  SE.T=eye(size(SE.K{1},1)+1,size(SE.K{1},2));
  SE.TIn=[zeros(size(SE.T,2),1);1];
  SE.DOF=(1:size(SE.T,1))'+.99;
  Freq=stack_get(mo1,'info','Freq','get');
  if isempty(Freq);Freq=stack_get(mo1,'info','RedFreq','get');
    Freq=logspace(log10(SE.w(1)),log10(SE.w(end)),1024)';
  end
  SE.w=Freq(:)*2*pi;
  SE.K{2}=SE.K{2}*(1+1i*fe_def('defeta',mo1)); % Loss factor
  d1=feval(fe_simul('@DfrfBIn'),SE);
  d1.TR=SE.TR; % Restitution of potential does not work
  
  C1=struct('X',{{d1.data,SE.Sens.lab}}, ...
      'Xlab',{{'Freq','Sens'}},'Y',(SE.Sens.cta*d1.def).');
  C1.ID=struct('po',SE.TR.data(:,1));
  C1.name='FrfSnapShot';
  out2=SE;
  % fe2ss('calld_piezo(''SolvePatchCapacity'')',mo1);
 elseif ~isempty(strfind(Cam,'ssold')) 
  %% Old strategy for ss building based on electrodes
  r1=stack_get(mo1,'info','Mesh','Get');
  r1=struct('data',[fix(r1.InputDOF) 0]);
  mo1=stack_set(mo1,'info','Electrodes',r1);
  mo1=fe_case(mo1,'remove','VIn');
  if isempty(stack_get(mo1,'info','EigOpt'))
   error('An ''info'',''EigOpt'' stack entry is needed');
  end
  [sys,TR]=fe2ss('free -fe_norm_tol=1e100',mo1);d1=TR;
  a=sys.a; b=sys.b; c=sys.c; d=sys.d; 
  % Remove rigid modes
  iz=find(TR.data(:,1)>20);iz=[iz;iz+size(TR.data,1)];
  ic=1; 
  s2=struct('a',a(iz,iz),'b',b(iz,:),'c',c(ic,iz),'d',d(ic,:));
  % damping of residual at 20 %
  % s2.a(end,end)=-2*0.20*sqrt(-s2.a(end,size(s2.a,2)/2));
  out2=s2; 
  
  if isempty(Freq);Freq=logspace(-3,log10(1.5),3000)'*TR.data(end,1);end
  C1=qbode(s2,Freq*2*pi,'struct');
  C1.ID=struct('po',TR.data(:,1)*[1 0], ...
        'LineProp',{{'LineStyle',':','color','k'}});
  C1.name='SSold';
 elseif ~isempty(strfind(Cam,'ss')) 
  %% SDT 6.6 strategy for ss building based on electrodes
  r1=stack_get(mo1,'info','Mesh','Get');
  [sys,TR]=fe2ss('free',mo1); % 

  d1=TR;
  a=sys.a; b=sys.b; c=sys.c; d=sys.d; 
  % Remove rigid modes
  iz=find(TR.data(:,1)>20);iz=[iz;iz+size(TR.data,1)];
  ic=1; 
  s2=struct('a',a(iz,iz),'b',b(iz,:),'c',c(ic,iz),'d',d(ic,:));
  % damping of residual at 20 %
  % s2.a(end,end)=-2*0.20*sqrt(-s2.a(end,size(s2.a,2)/2));
  out2=s2; 
  
  if isempty(Freq);Freq=logspace(-3,log10(1.5),3000)'*TR.data(end,1);end
  C1=qbode(s2,Freq*2*pi,'struct');
  C1.ID=struct('po',TR.data(:,1)*[1 0], ...
        'LineProp',{{'LineStyle',':','color','k'}});
  C1.name='SS';
  
 else
  if isempty(Freq);
   Freq=logspace(-3,log10(2),40)'*def.data(1);
   model=stack_set(mo1,'info','Freq',Freq);
  else; model=mo1; 
  end
  ofact('silent'); d1=fe_simul('dfrf',model);
  % Compute capacity (assumed to have Top input and Q-Top output)
  C1=fe_case('Sensobserve',mo1,'Q-Top',d1);
  if ~isempty(def);
   C1.ID=struct('po',def.data(:,1)*[1 0], ...
      'LineProp',{{'LineStyle',':','color','r'}});
  end
 end
 C1=sdsetprop(C1,'PlotInfo','sub','1 1','show','abs','scale','xlog;ylin');

if nargout==0;   iicom('curveinit',C1.name,C1);
else;out1=d1; out=C1;
end
%% #SolveEnd -2
else;error('Solve%s unknown',CAM);
end
elseif comstr(Cam,'pcond')
%% #Pcond coef : rescales electrical DOFs for better conditionning
%  default coef is 1e8
[CAM,Cam,r1]=comstr('cond',[-25 2],CAM,Cam); if isempty(r1);r1=1e8;end
 Case=evalin('caller','Case');model=evalin('caller','model');
 T=evalin('caller','T');
 pc=ones(size(Case.T,1),1);
 pc(fe_c(model.DOF,.21,'ind'))=r1;pc=diag(sparse(pc));
 T=pc*T;
 if isfield(Case,'TIn')
     Case.TIn=pc*Case.TIn; assignin('caller','Case',Case);
 end
 assignin('caller','T',T);
 %z=diag(k);mean(z(fe_c(model.DOF,.21,'ind')))/mean(z(fe_c(model.DOF,.21,'ind',2)))

 
%% #im : figure formatting ---------------------------------------------------
elseif comstr(Cam,'im')

% see also sdtweb('nida13','imw'), sdtweb('oscarui','imw')
 if nargin==2&&~isstruct(varargin{2}) % generate the calling string
   %comp12t('imw','Dispersion_.png')
   pw0=pwd;
   if ~ischar(varargin{2}); % Apply reshaping to figure
     gf=varargin{2};if isa(gf,'sdth');cf=gf; gf=gf.opt(1);end
     if ~ishandle(gf);figure(gf);plot([0 1]);end
     val=d_piezo(CAM);
     [st1,st2]=fileparts(pwd);
     if strcmpi(st2,'tex'); % force plotwd here
      val=[val {'@PlotWd',fullfile(pwd,'plots')}];       
     end
     if strcmpi(get(gf,'tag'),'iiplot');
        ci=get(gf,'userdata'); ci.ua.axProp=val;iiplot(ci);
     else; cingui('objset',gf,val)
     end
   elseif strcmpi(varargin{2},'.') 
    %% oscarui('imfw','.')
    st=sprintf('imwrite-objSet"@oscarui(''%s'')"-ftitle',varargin{1});
    comgui(st);
   else
       dbstack;keyboard
    cd(nida13('wd','plots'));
    st=sprintf('imwrite-objSet"@nida13(''%s'')"-ftitle%s',varargin{1:2});
    comgui(st);
    cd(pw0);
   end
 % z=sdt_table_generation('Rep');comstr(z{strcmpi(z(:,1),'LargeWide'),3},-30)
 elseif comstr(Cam,'imxfa') 
 % #ImXFA : basic display of transfers
   out={'@figure',{'position',[NaN,NaN,800,481]},'@exclude',{'legend.*'}, ...
       '@text',{'FontSize',14},'@axes',{'FontSize',14,'box','on'}, ...
       '@xlabel',{'FontSize',14},'@ylabel',{'FontSize',14}, ...
       '@zlabel',{'FontSize',14},'@title',{'FontSize',14}};
 elseif comstr(Cam,'imfa') 
 % #ImFa feplot for Cassem examples
   out={'@figure',{'position',[NaN,NaN,600,450]},'@exclude',{'legend.*'}, ...
       '@headsub',{'FontSize',14,'interpreter','none'}};

 elseif comstr(Cam,'imcheck') 
%% #Imcheck options
RO=varargin{carg};carg=carg+1;
if ~isfield(RO,'PlotWd');
  try; RO.PlotWd=oscarui('wd@PlotWd');
  catch; error('You must provide a PlotWd');
  end
end
if ~isfield(RO,'HtmWidth'); RO.HtmWidth=[];end % Default width
if ~isfield(RO,'RelPath'); RO.RelPath=2;end % Path for file history

out=RO;

%% 
 else; error('%s unknown',CAM);
 end
 
%% clean end
elseif comstr(Cam,'@');out=eval(CAM);
elseif comstr(Cam,'cvs')
 out='$Revision: 313 $  $Date: 2016-12-20 13:10:32 +0100 (Tue, 20 Dec 2016) $';
else; error('%s unknown',CAM);
end 
%% #End function
end


function [model,RO,RB]=meshRect(model,RO,RB);
%% #MexhRect Insert a rectangle
%if ~isempty(i2);RB.MatId=i2;end
%      RO.Rect(end+1,1:7)=[RB.xc RB.lx 0 RB.yc RB.ly 0 RB.ang]; 
%      if ~isempty(i2); RO.Rect(end,8:9)=i2; end % Mat/Pro of shape

%% Add rectangular patches
% RO.Rect=[xc lx nx, yc ly ny ang (MatId ProId)] 
  [RB,st,CAM]=cingui('paramedit -DoClean',[ ...
      ' nx(0#%g#" x divsion")'...
      ' ny(0#%g#" y division")'...
      ' xc(#%g#" position of corner")'...
      ' yc(#%g#" position of corner")'...
      ' ang(0#%g#" rotation of patch")'...
        ],{RB,RB.CAM});

 r3=[RO.MatId RO.ProId];r3(r3==0)=[];RO.PatchId(end+1,1)=r3(1);
 
 r3=RB.ang*pi/180;
 RB.x=[cos(r3) sin(r3) 0]'; RB.y=[-sin(r3) cos(r3) 0]'; % Local patch orient
 n1=[RB.xc RB.yc 0;RB.lx*RB.x(:)';RB.ly*RB.y(:)'];
 if RB.nx==0;RB.nx=ceil(RB.lx/RO.lc); end % Accept Lc rather than nx
 if RB.ny==0;RB.ny=ceil(RB.ly/RO.lc); end % Accept Lc rather than ny
 
 mo1=feutil(sprintf('objectquad %i %i',RO.MatId,RO.ProId),n1,RB.nx,RB.ny);
 r2=model.Node(:,5:7)-n1(ones(size(model.Node,1),1),:);
 r3=r2*RB.y; r2=r2*RB.x;
 i1=(r2>=-RO.lc & r2<=RB.lx+RO.lc & r3>=-RO.lc & r3 <=RB.ly+RO.lc);
 %feplot(feutil('addtest',model,mo2));fecom view3; fecom('shownodemark',model.Node(i1,1))
 mo2=model; % remove elements of the plate material only 
 [model.Elt,mo2.Elt]=feutil('removeElt withnode & proid ~=', ...
  model,model.Node(i1,1),RO.PatchId);
 mo2.Elt=feutil('selelt seledge',mo2);
 mo2=feutil('addtestmerge-NoOri;',mo2,mo1);
 mo2.Node=feutil('optimmodel',mo2);
 mo2=fe_gmsh('addline',mo2,'group1','seed',2);
 mo2=fe_gmsh('addline',mo2,'group2 & seledge','seed',2);
 GM=stack_get(mo2,'geom','GMSH','get');
 GM.PlaneSurface=[1 -2];
 mo2=stack_set(mo2,'geom','GMSH',GM);
 mo2=fe_gmsh(sprintf('write del.geo -lc %.15g -run -2 -v 0',RO.lc),mo2);
 if isempty(mo2.Elt); 
     error('Lc %g too small to fit, patch is empty',RO.lc);
 end
 mo2.Elt=feutil('selelt eltname tria',mo2);
 mo2.Elt=feutil('set group1 matid 1 proid 1',mo2);
 mo2=feutil('addtest-merge-NoOri;',mo2,mo1);
 model=feutil('addtest-merge-NoOri;',model,mo2);

end

function [model,RO,RB]=meshCirc(model,RO,RB);
%% #MeshCirc Add circular patches
% RO global options, RB current options of patch
% RO.Circ=[xc yc rc lc (MatId ProId)] 

[RB,st,CAM]=cingui('paramedit -DoClean',[ ...
      ' xc(#%g#" center [length unit]")'...
      ' yc(#%g#" center [length unit]")'...
      ' ang(0#%g#" rotation of patch")'...
      ' MatId(#%g#" target MatId")'...
      ' ProId(#%g#" target ProId")'...
        ],{RB,RB.CAM});
 if ~isfield(RB,'lc')||isempty(RB.lc);RB.lc=RO.lc;end
 
 r3=[RO.MatId RO.ProId];r3(r3==0)=[];RO.PatchId(end+1,1)=r3(1);
 mo1=feutil(sprintf('objectcircle %.15g %.15g 0  %.15g 0 0 1 %i', ...
     RB.xc,RB.yc,RB.rc,2*pi*RB.rc/RB.lc));
 n1=[RB.xc,RB.yc 0];
 r3=model.Node(:,5:7)-n1(ones(size(model.Node,1),1),:);
 i1=sqrt(sum(r3.^2,2))<=RB.rc+RB.lc; 
 if nnz(i1)==0; warning('No node for Circle(MatId%i), skipping',RO.MatId);
     return;
 end
 % feplot(feutil('addtest',model,mo2));fecom view3; fecom('shownodemark',model.Node(i1,1))
 mo2=model;[model.Elt,mo2.Elt]=feutil('removeElt withnode & proid ~=', ...
     model,model.Node(i1,1),RO.PatchId);
 mo2.Elt=feutil('selelt seledge',mo2); 
 mo2.Node=feutil('optimmodel',mo2);
 mo2=fe_gmsh('addline',mo2,'group1','seed',2);
 mo2=fe_gmsh('AddFullCircle -loop2',mo2,[n1;n1+[RB.rc 0 0]; 0 0 1]);  
 GM=stack_get(mo2,'geom','GMSH','get');
 GM.PlaneSurface=[2 0;1 -2];
 %GM.Post={'Mesh.SubdivisionAlgorithm=1;'}; % GMSH implement 2quad
 RO.RefineToQuad=1;
 mo2=stack_set(mo2,'geom','GMSH',GM);
 if isempty(mo2.Elt); 
     error('Lc %g too small to fit, patch is empty',RO.lc);
 end

 mo3=fe_gmsh(sprintf('write del.geo -lc %.15g -run -2 -v 0',RO.lc),mo2);
 if isempty(mo3.Elt)
   feplot(feutil('addtest',mo2,mo1));
 end
 mo2=mo3; mo2.Elt=feutil('selelt eltname quad| eltname tria',mo2);
 mo2.Elt=feutilb('separatebyprop',mo2);
 if isempty(mo2.Elt); error('Lc too small to fit patch');end
 mo2.Elt=feutil('set group2 matid 1 proid 1',mo2); % patch
 mo2.Elt=feutil(sprintf('set group1 matid %i proid %i',RO.MatId,RO.ProId),mo2);
 model=feutil('addtest-merge-NoOri;',model,mo2);

end
