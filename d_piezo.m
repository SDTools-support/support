function [out,out1,out2]=d_piezo(varargin); %#ok<*STOUT>

% D_PIEZO Support for demonstrations related to piezo-electricity
%
% See <a href="matlab: sdtweb _taglist d_piezo">TagList</a>
%
% Arnaud Deraemaeker, ULB and Etienne Balmes, SDTools


%       Copyright (c) 1990-2025 by SDTools, All Rights Reserved.
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

if comstr(Cam,'buildc1');
  r1=varargin{2};
  C1=struct('X',[],'Xlab',[],'Y',[]);
  w=varargin{3};
  if isa(r1,'ss')
    C1.Y=freqresp(r1,w); C1.Y=permute(C1.Y,[3 1 2]); 
    C1.X{1}=w(:)/2/pi;C1.Xlab{1}={'Frequency','Hz',[18 0 0 0 -1]};
    C1.X{2}=r1.OutputName;C1.Xlab{2}='Out';
    C1.X{3}=r1.InputName;C1.Xlab{3}='In';
    C1.Ylab=[2 3];
  elseif isfield(r1,'a')
    C1=qbode(r1,w,'{stra2,otStruct}');
    if isfield(r1,'InputName');C1.X{3}=r1.InputName; end
    if isfield(r1,'OutputName');C1.X{2}=r1.OutputName; end
    
  else
     C1.X{1}=varargin{2}; 
     C1.X{2}={varargin{4}}; C1.X{3}={varargin{5}}; % Freq - Output Label - Input Label
     C1.Xlab{1}={'Frequency','Hz',[18 0 0 0 -1]};
     C1.Xlab{2}='Outputs'; C1.Xlab{3}='In';
     C1.Y=varargin{3}; C1.Ylab=[2 3]; C1.name='DFRF';
     C1=sdsetprop(C1,'PlotInfo','sub','magpha','scale','xlin;ylog');
  end
out=C1;


%% #Script -------------------------------------------------------------------
elseif comstr(Cam,'script');[CAM,Cam]=comstr(CAM,7);

if comstr(Cam,'tutopzpatchext')
%% #TutoPzPatchExt : Piezoelectric extension patch - Statics -2
% see sdtweb pz_fe#tutopzpatchext 

%% BeginSource sdtweb('_example','pz_fe.tex#tutopzpatchext')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzPatchExt')
d_piezo('DefineStyles');

%% Step 1 Build mesh - Define electrodes 
% Meshing script can be viewed with sdtweb d_piezo('MeshPatch')
model=d_piezo('MeshPatch lx=1e-2 ly=1e-2 h=2e-3 nx=1 ny=1 nz=1');
% Define electrodes
model=p_piezo('ElectrodeMPC Top -ground',model,'z==2e-3');
model=p_piezo('ElectrodeMPC Bottom -Input "Free patch"',model,'z==0');
p_piezo('TabInfo',model)
%% Step 2 Define material properties
model.pl=m_piezo('dbval 1 -elas 2 PIC_255');
p_piezo('TabDD',model) % creates the table with full set of matrices
%% Step 3 Compute static response
% to avoid rigid body mode
model=stack_set(model,'info','Freq',10);
def=fe_simul('dfrf',model); def.lab={'Free patch, axial'};
def.fun=[0 1]; def=feutil('rmfield',def,'data','LabFcn');

% Append mechanically constrained structure
%    can't call fe_simul because no free DOF
%    see code with sdtweb d_piezo('scriptFullConstrain')
def=d_piezo('scriptFullConstrain',model,def);
def.lab{2}='Constrained patch, axial';
%% Step 4 Visualize deformed shape
cf=feplot(model,def);
% Electric field representation
p_piezo('viewElec EltSel "matid1" DefLen 20e-4 reset',cf);
fecom('colormap',[1 0 0]);fecom('undef line');iimouse('resetview');
cf.mdl.name='E-field'; % Sets figure name
d_piezo('SetStyle',cf); feplot(cf);
%% Step 5 : check constitutive law
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
%% Step 6 Check capacitance values
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
  iimouse('zoom reset');

% Compute the value of the total charge (from reaction at electrical dof)
% Compare with analytical values
p_piezo('electrodeTotal',cf) %
disp('Theoretical values of capacitance')
disp([{'CT';'CS'} num2cell([CT;CS])])
cf.mdl.name='Charge'; % Sets figure name

d_piezo('SetStyle',cf); feplot(cf);
fecom('ch 1')


%% EndSource EndTuto

elseif comstr(Cam,'tutopzpatchshear')
%% #TutoPzPatchShear : Piezoelectric shear patch - statics -2

% see sdtweb pz_fe#tutopzpatchshear 
%
%% BeginSource sdtweb('_example','pz_fe.tex#tutopzpatchshear')

% Init working directory for figure generation
d_piezo('SetPlotwd');
d_piezo('DefineStyles');

% See full example as MATLAB code in d_piezo('ScriptPzPatchShear')
d_piezo('DefineStyles');

%% Step 1 Build mesh and define electrodes
%Meshing script can be viewed with sdtweb d_piezo('MeshPatch')
model=d_piezo('MeshPatch lx=1e-2 ly=1e-2 h=2e-3 nx=1 ny=1 nz=1');

% Define electrodes
model=p_piezo('ElectrodeMPC Top -ground',model,'z==2e-3');
model=p_piezo('ElectrodeMPC Bottom -Input "Free patch"',model,'z==0');

% Rotate basis to align poling direction with y (-90 deg around x)
model.bas=basis('rotate',[],'rx=-90',1); %create local basis with id=1
model=feutil('setpro 1 COORDM=1',model); % assign basis with id=1 to pro=1
%% Step 2 Compute static response
% to avoid rigid body mode
model=stack_set(model,'info','Freq',10);
def=fe_simul('dfrf',model); def.lab={'Free patch, shear'};
def.fun=[0 1]; def=feutil('rmfield',def,'data','LabFcn');

% Append mechanically constrained structure
%    can't call fe_simul because no free DOF
%    see code with sdtweb d_piezo('scriptFullConstrain')
def=d_piezo('scriptFullConstrain',model,def);
def.lab{2}='Constrained patch, shear';
%% Step 3 Vizualise deformed shape
cf=feplot(model,def); fecom('undef line');

% Electric field representation
p_piezo('viewElec EltSel "matid1" DefLen 20e-4 reset',cf);
cf.mdl.name='E-field'; % Sets figure name
d_piezo('SetStyle',cf); feplot(cf);
%% Step 4 : Check constitutive law
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
    iimouse('zoom reset');

%% Step 5 Check capacitance
% Compute the value of the total charge (from reaction at electrical dof)
% Ccompare with analytical values
p_piezo('electrodeTotal',cf) %
disp('Theoretical values of capacitance')
disp([{'CT';'CS'} num2cell([CT;CS])])

%% EndSource EndTuto

elseif comstr(Cam,'tutopzdiskimpedance')
%% #TutoPzDiskImpedance : Piezoelectric disk impedance -2

% see sdtweb pz_fe#tutopzdiskimpedance 
%
%% BeginSource sdtweb('_example','pz_fe.tex#tutopzdiskimpedance')

% Init working directory for figure generation
d_piezo('SetPlotWd')
% See full example as Matlab code in d_piezo('ScriptTutoPzDiskImpedance')
d_piezo('Definestyles');

%% Step 1 Build and represent mesh and electrodes
model=d_piezo('MeshPIC181disk th=2e-3 r=8e-3 ner=10 nez=4 nrev=16');
feplot(model); cf=fecom; cf.mdl.name='PIC 181 piezo disk mesh'; iimouse('resetview')
d_piezo('setstyle',cf)
% Visualize electrodes
fecom('curtabCase',{'Top Actuator';'Bottom Actuator'}) % 
fecom(';showline;proviewon;triax') % 
cf.mdl.name='PIC 181 piezo disk electrodes'
d_piezo('setstyle',cf)
%% Step 2 : Define range of frequencies and compute dynamic response
frq=linspace(20e3,200e3,256);
def=fe_simul('dfrf',stack_set(model,'info','Freq',frq)); 


% visualize potential
feplot(model,def); cf=fecom; 
fecom(';showpatch;colordata21;'); cf.mdl.name='PIC 181 piezo disk voltage'
d_piezo('setstyle',cf) ; 
cf.osd_('cbtr{string,Voltage(V)}')
fecom('colorscaleone') %To have the correct scale

% View electric field
fecom(';showline;scd 1e-4')
p_piezo('viewElec EltSel "matid1" DefLen 1e-4',cf); 
cf.mdl.name='PIC 181 piezo disk E-field'
% To have a single color change clim (must be done with axProp to bypass normal)
st=cf.ua.axProp; st(3:4)={'@axes',{'clim',[480 510]}};cf.ua.axProp=st;
d_piezo('setstyle',cf)
cf.osd_('cbtr{string,E(V/m)}')
%% Step 3: Compute q/V as a function of the frequency
sens=fe_case(model,'sens');
C1=fe_case('SensObserve -DimPos 2 3 1',sens,def);
C1=sdsetprop(C1,'PlotInfo','sub','magpha','scale','xlin;ylog');
ci=iiplot; 
iicom(ci,'curveInit',C1.name,C1); iicom('submagpha');
d_piezo('setstyle',ci);
%% Step 4: Compute and plot electric impedance
% extract impedance
C2=C1; C2.Y=1./(2*pi*1i*C2.X{1}.*C2.Y); C2.X{2}={'Imp(Ohm)'};
iicom(ci,'curveInit',C2.name,C2); iicom('submagpha');
d_piezo('setstyle',ci);

%% EndSource EndTuto

elseif comstr(Cam,'tutopzbeamcol')
%% #TutoPzBeamCol : 3D beam with collocated sensors and actuators -2

% see sdtweb pz_fe#tutopzbeamcol 
%
%% BeginSource sdtweb('_example','io_theory.tex#tutopzbeamcol')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzBeamCol')
d_piezo('DefineStyles');

%% Step 1 : meshing, boundary conditions and point sensor definition
model = femesh('test ubeam');
% BC : fix top
model=fe_case(model,'FixDOF','Clamp','z==0');
% Introduce a point displacement sensor and visualize
% sdtweb sensor#slab % URN based definition of sensors
model = fe_case(model,'SensDOF','Point Sensors',{'104:x'});
cf=feplot(model); iimouse('resetview');

fecom('showFiAlpha') % Make mesh transparent 

% Visualize sensor, specify arrow length and width
sdth.urn('Tab(Cases,Point Sensors){Proview,on,deflen,.25}',cf)
sdth.urn('Tab(Cases,Point Sensors){arProp,"linewidth,2"}',cf)

cf.mdl.name='Ubeam PS1'; % Model name for title
d_piezo('SetStyle',cf); feplot(cf);

% Insert the number for the sensor :
fecom('textnode',104,'fontsize',14)
%% Step 2 : Introduce collocated force actuator and compute response
model=fe_case(model,'DofLoad SensDof','Collocated Force','Point Sensors:1') 
% 1 for first sensor if there are multiple
% compute static response and visualize
model=stack_set(model,'info','oProp',mklserv_utils('oprop','CpxSym'));
d0=fe_simul('dfrf',stack_set(model,'info','Freq',0)); % Static response
feplot(model,d0); fecom(';scd .3;undef line'); 

% Visualize sensors
fecom curtabcases 'Point Sensors'
sdth.urn('Tab(Cases,Point Sensors){Proview,on,deflen,.25}',cf)
sdth.urn('Tab(Cases,Point Sensors){arProp,"linewidth,2"}',cf)
fecom('textnode',104,'fontsize',14); 

% Title
cf.mdl.name='Ubeam PS1 Static';
d_piezo('SetStyle',cf); feplot
%% Step 3 :  Compute dynamic response in freq band of first 5 modes and plot
def=fe_eig(model,[5 10 0]);
frq=linspace(0,def.data(6)-(def.data(6)-def.data(5))/2,300);
d1=fe_simul('dfrf',stack_set(model,'info','Freq',frq)); % Dynamic response

% Construct projection matrix in sens.cta and plot collocated FRF
sens=fe_case(model,'sens');
% Plot the 4 FRFs, FRFs 1 and 4 are collocated
C1=fe_case('SensObserve -DimPos 2 3 1',sens,d1);
C1=sdsetprop(C1,'PlotInfo','sub','magpha','scale','xlin;ylog');
ci=iiplot; 
iicom(ci,'curveInit',C1.name,C1); iicom('submagpha');
d_piezo('setstyle',ci)
%% Step 4 : multiple collocated sensors and actuators
% Introduce two sensors and visualize
model = fe_case(model,'SensDOF','Point sensors',{'104:x';'344:y'});
cf=feplot(model);
% Visualize sensors :
fecom('showfialpha') %
fecom proviewon
fecom curtabcases 'Point Sensors' % Shows the case 'Point Sensors'

% Improve figure
sdth.urn('Tab(Cases,Point Sensors){Proview,on,deflen,.25}',cf)
sdth.urn('Tab(Cases,Point Sensors){arProp,"linewidth,2"}',cf)

% Title
cf.mdl.name='Ubeam MS1';

% Insert the number for the sensor :
fecom('textnode',[104 344],'FontSize',14)
d_piezo('SetStyle',cf)
%% Step 5 : Introduce collocated force actuators
model=fe_case(model,'DofLoad SensDof','Collocated Force','Point sensors:1:2')

% compute static response and visualize (two static responses)
d0=fe_simul('dfrf',stack_set(model,'info','Freq',0)); % Static response
feplot(model,d0); fecom(';scd .3; undef line')
fecom('textnode',[104 344],'FontSize',14)

% Style
d_piezo('SetStyle',cf);  cf.os_('LgMl-FontSize14');% Keep both mdl.name and title
%% Step 6 :  Compute dynamic response in freq band of first 5 modes and plot
frq=linspace(0,def.data(6)-(def.data(6)-def.data(5))/2,300);
d1=fe_simul('dfrf',stack_set(model,'info','Freq',frq)); % Dynamic response

% Construct projection matrix in sens.cta and project resp on sensor
sens=fe_case(model,'sens');

% Plot the 4 FRFs, FRFs 1 and 4 are collocated
C1=fe_case('SensObserve -DimPos 2 3 1',sens,d1);
C1=sdsetprop(C1,'PlotInfo','sub','magpha','scale','xlin;ylog');
ci=iiplot; 
iicom(ci,'curveInit',C1.name,C1); iicom('submagpha');
d_piezo('setstyle',ci);
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopzbeamncol')
%% #TutoPzBeamNCol : 3D beam with non collocated sensors and actuators -2

% see sdtweb io_theory#tutopzbeamncol 
%
%% BeginSource sdtweb('_example','io_theory.tex#tutopzbeamncol')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example in d_piezo('ScriptTutoPzBeamNCol')
d_piezo('DefineStyles'); % Init styles for figures

%% Step 1 : meshing, BC and sensors def
model = femesh('test ubeam');
% BC : fix top
model=fe_case(model,'FixDOF','Clamp','z==0');
% Define sensors
model = fe_case(model,'SensDOF','Point Sensors',{'104:x';'207:y'});
cf=feplot(model); iimouse('resetview');

% Make mesh transparent :
fecom('showfialpha') %

% Visualize sensor
 fecom proviewon
 fecom curtabcases 'Point Sensors' % Shows the case 'Point Sensors'

% Improve figure
% Arrow length and thickness
sdth.urn('Tab(Cases,Point Sensors){Proview,on,deflen,.25}',cf)
sdth.urn('Tab(Cases,Point Sensors){arProp,"linewidth,2"}',cf)

cf.mdl.name='Ubeam MSNC'; % Model name for title
d_piezo('SetStyle',cf); feplot(cf);

% Insert the number for the sensor :
fecom('textnode',[104 207],'fontsize',14)
%% Step 2 : Define point actuators
% relative force between DOFs 207x and 241x and one point loads at DOFs 207y
data  = struct('DOF',[207.01;241.01;207.02],'def',[1 0;-1 0;0 1]);
model=fe_case(model,'DofLoad','Actuators',data); %
cf=feplot(model); fecom('showline')
 
fecom curtabcases 'Actuators' % Shows the case 'Actuators'

% Improve figure
% Arrow length and thickness
sdth.urn('Tab(Cases,Actuators){Proview,on,deflen,.25}',cf)
sdth.urn('Tab(Cases,Actuators){arProp,"linewidth,2"}',cf)

cf.mdl.name='Ubeam MANC 1'; % Model name for title
d_piezo('SetStyle',cf); feplot(cf);

% Insert the number for the sensor :
fecom('textnode',[241 207],'fontsize',14)
% Visualize second combination
cf.CStack{'Actuators'}.Sel.ch=2;sdth.urn('Tab(Cases,Actuators)',cf) % second
cf.mdl.name='Ubeam MANC 2'; % Model name for title
d_piezo('SetStyle',cf); feplot(cf);
%% Step 3 : compute static response and visualize
model=stack_set(model,'info','oProp',mklserv_utils('oprop','CpxSym'));
d0=fe_simul('dfrf',stack_set(model,'info','Freq',0)); % Static response
feplot(model,d0); fecom(';scd .1; undef line;')


% Title
d_piezo('SetStyle',cf); cf.os_('LgMl-FontSize14');% Keep both mdl.name and title
fecom('textnode',[241 207],'FontSize',14)
%% Step 4 :  Compute dynamic response in freq band of first 5 modes and plot
% Compute modes and frequencies
def=fe_eig(model,[5 10 0]);
% compute response
frq=linspace(0,def.data(6)-(def.data(6)-def.data(5))/2,300);
d1=fe_simul('dfrf',stack_set(model,'info','Freq',frq)); % Dynamic response

% Construct projection matric in sens.cta and project resp on sensor
sens=fe_case(model,'sens');

% Plot FRFs
C1=fe_case('SensObserve -DimPos 2 3 1',sens,d1);
C1=sdsetprop(C1,'PlotInfo','sub','magpha','scale','xlin;ylog');
ci=iiplot;
iicom(ci,'curveInit',C1.name,C1); iicom('submagpha');
d_piezo('setstyle',ci);
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopzbeamsurfvol')
%% #TutoPzBeamSurfVol : 3D beam with Surface and Volume forces -2

% see sdtweb io_theory#tutopzbeamsurfvol 
%
%% BeginSource sdtweb('_example','io_theory.tex#tutopzbeamsurfvol')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example in d_piezo('ScriptTutoPzBeamSurfVol')
 d_piezo('DefineStyles'); % Define styles for figures
 %% Step 1 Apply  a volumic load and represent
model = femesh('testubeam');
data=struct('sel','groupall','dir',[0 32 0]);
data2=struct('sel','groupall','dir',{{0,0,'(z-1).^3.*x'}});
model=fe_case(model,'FVol','Constant',data, ...
                     'FVol','Variable',data2);

% Visualize loads
cf=feplot(model); iimouse('resetview');

% Make mesh transparent :
fecom('showfialpha') %

% Visualize Load
fecom proviewon

% Improve figure
sdth.urn('Tab(Cases,Constant){deflen,.5,arProp,"linewidth,2"}',cf)
fecom curtabcases 'Constant' % Shows the case 'Constant'
% Set style and print
cf.mdl.name='Ubeam VLoad-Cst'; % Model name for title
d_piezo('SetStyle',cf); feplot(cf);
% Visualize variable Load and print
  sdth.urn('Tab(Cases,Variable){deflen,.5,arProp,"linewidth,2"}',cf)
 fecom curtabcases 'Variable' % Shows the case 'Variable'
 cf.mdl.name='Ubeam VLoad-Var'; % Model name for title
 d_piezo('SetStyle',cf); feplot(cf);
% Visualize Constant and Variable loads with colors
Load = fe_load(model); cf=feplot(model,Load); cf.mdl.name='Ubeam VLoad Color'; 
% display as color-code to see change of vol force with z and x
fecom(';showpatch;colordataz;scd .0001;'); 
d_piezo('SetStyle',cf); feplot(cf); fecom('colorbar on')
%% Step 2 :  Apply a surface load case in a model using selectors
data=struct('sel','x==-.5', ...
             'eltsel','withnode {z>1.25}','def',1,'DOF',.19);
model=fe_case(model,'Fsurf','Surface load',data); cf=feplot(model);
% Visualize Load
fecom proviewon
fecom curtabcases 'Surface Load' % Shows the case 'Constant'
fecom showline

sdth.urn('Tab(Cases,Surface load){deflen,.5,arProp,"linewidth,2"}',cf)
cf.mdl.name='Ubeam SLoad'; % Model name for title
d_piezo('SetStyle',cf); feplot(cf);
%% Step 3 : Applying a surfacing load case in a model using node lists
data=struct('eltsel','withnode {z>1.25}','def',1,'DOF',.19);
NodeList=feutil('findnode x==-.5',model);
data.sel={'','NodeId','==',NodeList};
model=fe_case(model,'Fsurf','Surface load 2',data); cf=feplot(model);
fecom proviewon
fecom curtabcases 'Surface load 2' %
fecom showline

sdth.urn('Tab(Cases,Surface load){deflen,.5,arProp,"linewidth,2"}',cf)
cf.mdl.name='Ubeam SLoad2'; % Model name for title
d_piezo('SetStyle',cf); feplot(cf);
%% Step 4 : Applying a surfacing load case in a model using sets

% Define a face set
[eltid,model.Elt]=feutil('eltidfix;',model);
i1=feutil('findelt withnode {x==-.5 & y<0}',model);i1=eltid(i1);
i1(:,2)=2; % fourth face is loaded
data=struct('ID',1,'data',i1,'type','FaceId');
model=stack_set(model,'set','Face 1',data);

% define a load on face 1
data=struct('set','Face 1','def',1,'DOF',.19);
model=fe_case(model,'Fsurf','Surface load 3',data); cf=feplot(model);
sdth.urn('Tab(Cases,Surface load 3){deflen,.5,arProp,"linewidth,2"}',cf)
fecom proviewon
fecom curtabcases 'Surface load 3' %

cf.mdl.name='Ubeam SLoad3'; % Model name for title
d_piezo('SetStyle',cf); feplot(cf);
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopzbeamuimp')
%% #TutoPzBeamUImp : 3D beam with imposed displacement -2

% see sdtweb io_theory#tutopzbeamuimp 
%
%% BeginSource sdtweb('_example','io_theory.tex#tutopzbeamuimp')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example in d_piezo('ScriptTutoPzBeamUimp')
d_piezo('DefineStyles');
%% Step 1 Meshing and BC
model = femesh('test ubeam');
% BC : Impose displacement -  Fix all other dofs for Base
model=fe_case(model,'FixDof','Clamping','z==0 -DOF 1 3');
%% Step 2 Apply base displacement in y-direction
% find node z==0
nd=feutil('find node z==0',model); 
data.DOF=nd+.02; data.def=ones(length(nd),1);
model=fe_case(model,'DofSet','Uimp',data);
cf=feplot(model); iimouse('resetview')

% Visualize
fecom proview on
fecom curtabcases 'Uimp' % Shows the case 'Constant'
% Set style and print
cf.mdl.name='Ubeam Uimp'; % Model name for title
d_piezo('SetStyle',cf); feplot(cf);
%% Step 3 : Introduce a point displacement sensor and visualize
model = fe_case(model,'SensDOF','Point Sensors',{'104:y'});
cf=feplot(model);

% Make mesh transparent :
fecom('showfialpha') %

% Visualize sensor and actuator
fecom proviewon
sdth.urn('Tab(Cases,Point Sensors){arProp,"linewidth,2"}',cf)
sdth.urn('Tab(Cases,Uimp){arProp,"linewidth,2"}',cf)
fecom('curtabCase','#(Uimp|Point)')

cf.mdl.name='Ubeam Uimp Sens Act'; % Model name for title
d_piezo('SetStyle',cf); feplot(cf);

% Insert the number for the sensor :
fecom('textnode',[104],'fontsize',14)
%% Step 4 : Compute modes and transfer functions
def=fe_eig(model,[5 10 0]);
feplot(model,def)

% Compute dynamic response in freq band of first 5 modes and plot
frq=linspace(0,def.data(6)-(def.data(6)-def.data(5))/2,300);
d1=fe_simul('dfrf',stack_set(model,'info','Freq',frq)); % Dynamic response

% Construct projection matrix in sens.cta and plot collocated FRF
sens=fe_case(model,'sens');
% Plot the 4 FRFs, FRFs 1 and 4 are collocated
C1=fe_case('SensObserve -DimPos 2 3 1',sens,d1);
C1=sdsetprop(C1,'PlotInfo','sub','magpha','scale','xlin;ylog');
ci=iiplot; 
iicom(ci,'curveInit',C1.name,C1); iicom('submagpha');
d_piezo('setstyle',ci);
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopzbeamdispvelacc')
%% #TutoPzBeamDispVelAcc : 3D beam with displ,vel and acc sensor -2

% see sdtweb io_theory#tutopzbeamdispvelacc
%
%% BeginSource sdtweb('_example','io_theory.tex#tutopzbeamdispvelacc')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example in d_piezo('ScriptTutoPzBeamDispVelAcc')
d_piezo('DefineStyles');

%% Step 1 : meshing, BC and compute modeshapes
model = femesh('test ubeam');
model=fe_case(model,'FixDOF','Clamp','z==0');
def=fe_eig(model,[5 10 0]);
%% Step 2 : Introduce a point displ/vel/acc sensor and collocated force
model = fe_case(model,'SensDOF','Sensors',{'104:x';'104:vx';'104:ax'});
model=fe_case(model,'DofLoad SensDof','Collocated Force','Sensors:1'); 
% 1 for first sensor if there are multiple
%% Step 3 : Compute response and plot FRFs
frq=linspace(10,def.data(6)-(def.data(6)-def.data(5))/2,300);

d1=fe_simul('dfrf',stack_set(model,'info','Freq',frq)); % Dynamic response

% Construct projection matrix in sens.cta and plot collocated FRF
sens=fe_case(model,'sens');

% Plot the 3 FRFs
C1=fe_case('SensObserve -DimPos 2 3 1',sens,d1);
C1=sdsetprop(C1,'PlotInfo','sub','magpha','scale','xlin;ylog');
ci=iiplot; iicom(ci,'curveInit',C1.name,C1); iicom('submagpha');
d_piezo('setstyle',ci);
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopzplate4pzt')
%% #TutoPzPlate4pzt : Cantilever plate with 4 piezo patches -2

% see sdtweb pzio_theory#tutopzplate4pzt
%
%% BeginSource sdtweb('_example','pzio_theory.tex#tutopzplate4pzt')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example in d_piezo('ScriptTutoPzPlate4Pzt')
d_piezo('DefineStyles');

%% Step 1 - Build model and visualize
model=d_piezo('MeshULBplate');  % creates the model
model=fe_case(model,'FixDof','Cantilever','x==0'); % Clamp plate
% Set modal default zeta = 0.01
model=stack_set(model,'info','DefaultZeta',0.01);
cf=feplot(model); fecom('colordatagroup'); set(gca,'cameraupvector',[0 1 0])
cf.mdl.name='Plate_4pzt';
d_piezo('SetStyle',cf); feplot(cf);
p_piezo('TabDD',model);      % List piezo constitutive laws
r1=p_piezo('TabInfo',model); % List piezo related properties
%% Step 2 - Define actuators and sensors and visualize
nd=feutil('find node x==463 & y==100',model);
model=fe_case(model,'SensDof','Tip',{[num2str(nd) ':z']}); % Displ sensor
i1=p_piezo('TabInfo',model);i1=i1.Electrodes(:,1);
model=fe_case(model,'DofSet','V-Act',struct('def',1,'DOF',i1(1)+.21, ...%Act
    'Elt',feutil('selelt proid 104',model))); % Elt defined for display
model=p_piezo(sprintf('ElectrodeSensQ  %i Q-Act',i1(1)),model); % Charge sensors
model=p_piezo(sprintf('ElectrodeSensQ  %i Q-S1',i1(2)),model);
model=p_piezo(sprintf('ElectrodeSensQ  %i Q-S2',i1(3)),model);
model=p_piezo(sprintf('ElectrodeSensQ  %i Q-S3',i1(4)),model);
% Fix ElectrodeSensQ dofs to measure resultant (charge)
model=fe_case(model,'FixDof','SC*S1-S3',i1(2:end)+.21);
cf=feplot(model); fecom('view3')
cf.mdl.name='Plate_4pzt'; d_piezo('SetStyle',cf); feplot(cf);
fecom('showfialpha')
fecom('proviewon')

% Arrow length and thickness
sdth.urn('Tab(Cases,Tip){Proview,on,deflen,20}',cf)
sdth.urn('Tab(Cases,Tip){arProp,"linewidth,2"}',cf)
fecom('curtabCase',{'Tip';'V-Act'})
d_piezo('SetStyle',cf); feplot(cf);

% Insert the number for the sensor :
fecom('textnode',[nd],'fontsize',14)
cf.mdl.name='Plate_4pzt_Vact-Tip'; d_piezo('SetStyle',cf); feplot(cf);
fecom('curtabCase',{'Q-S1';'Q-S2'})
cf.mdl.name='Plate_4pzt_QS1-2';d_piezo('SetStyle',cf); feplot(cf);
%% Step 3 Compute static and dynamic response
model=stack_set(model,'info','oProp',mklserv_utils('oprop','CpxSym'));
d0=fe_simul('dfrf',stack_set(model,'info','Freq',0)); % direct refer frf at 0Hz
cf=feplot(model,d0); fecom(';view3;scd 20;colordatagroup;undefline')
cf.mdl.name='Plate_4pzt'; d_piezo('setstyles',cf);
f=linspace(1,100,400); % in Hz
d1=fe_simul('dfrf',stack_set(model,'info','Freq',f(:))); % direct refer frf
sens=fe_case(model,'sens');

% Plot FRFs
C1=fe_case('SensObserve -DimPos 2 3 1',sens,d1);
C1=sdsetprop(C1,'PlotInfo','sub','magpha','scale','xlin;ylog');
ci=iiplot; 
iicom(ci,'curveInit',C1.name,C1); iicom('submagpha');
d_piezo('setstyle',ci);
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopzaccshaker')
%% #TutoPzAccShaker : Piezoelectric shaker and accelero -2

% see sdtweb pzio_theory#tutopzaccshaker
%
%% BeginSource sdtweb('_example','pzio_theory.tex#tutopzaccshaker')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzAccShaker')
d_piezo('DefineStyles');

%% Step 1 - Build mesh and visualize
% Meshing script,open with sdtweb d_piezo('MeshPiezoShaker')
model=d_piezo('MeshPiezoShaker');
cf=feplot(model); fecom('colordatapro');
%% Step 2 - Define actuators and sensors
  % -input "In" says it will be used as a voltage actuator
model=p_piezo('ElectrodeMPC Top Actuator -input "Vin-Shaker"',model,'z==-0.01');
  % -ground generates a v=0 FixDof case entry
model=p_piezo('ElectrodeMPC Bottom Actuator -ground',model,'z==-0.012');
% Voltage sensor will be used - remove charge sensor
model=fe_case(model,'remove','Q-Top sensor');
% Visualize electrodes
fecom(';showline;proviewon')
fecom('curtabCase',{'Top sensor';'Bottom sensor';'Top Actuator';'Bottom Actuator'}) % 
% Visualize Vin electrode
fecom curtabcases Vin-Shaker %
% Visualize VSens electrode
fecom curtabcases 'V-Top sensor'
r1=p_piezo('TabInfo',model); % List piezo related properties
%% Step 3: compute and visualize dynamic response
ofact('silent'); f=logspace(3,5.3,400)';
model=stack_set(model,'info','oProp',mklserv_utils('oprop','CpxSym'));
d1=fe_simul('dfrf',stack_set(model,'info','Freq',f(:))); % direct refer frf

% Project on sensor and create output
sens=fe_case(model,'sens');
C1=fe_case('SensObserve -DimPos 2 3 1',sens,d1);
C1.X{2}{1}='V-sensor(V)'
ci=iiplot; 
iicom(ci,'curveInit',C1.name,C1); iicom('submagpha');
d_piezo('setstyle',ci);
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopzmeshingbasics')
%% #TutoPzMeshingBasics : Plate with 4 pzt patches:manual meshing -2

% see sdtweb pzplatemeshing#tutopzmeshingbasics
%
%% BeginSource sdtweb('_example','pzplatemeshing.tex#tutopzmeshingbasics')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example in d_piezo('ScriptTutoPzMeshingBasics')
d_piezo('DefineStyles');

%% Step 1 - Mesh the plate
 model=struct('Node',[1 0 0 0 0 0 0],'Elt',[]);
 model=feutil('addelt',model,'mass1',1);
 % Note that the extrusion values are chosen to include the patch edges
 dx=[linspace(0,15,3) linspace(15,15+55,10) linspace(15+55,463-5,50) 463];
 model=feutil('extrude 0 1 0 0',model,.001*unique(dx));
 dy=[linspace(0,12,3) linspace(12,12+25,5) linspace(12+25,63,5) ...
    linspace(63,63+25,5) linspace(63+25,100,3)];
 model=feutil('extrude 0 0 1 0',model,.001*unique(dy));
 %% Step 2 - Set patch areas and set different properties
 model.Elt=feutil('divide group 1 withnode{x>.015 & x<.07 & y>.013 & y<.037 }',model);
 model.Elt=feutil('divide group 2 withnode{x>.015 & x<.07 & y>.064 & y<.088 }',model);
 model.Elt=feutil('set group 1 proId3',model);
 model.Elt=feutil('set group 2 proId4',model);
%% Step 3 - Material Properties
model.pl=m_elastic('dbval 1 Aluminum');
model.pl=m_piezo(model.pl,'dbval 3 -elas2 SONOX_P502_iso');
% To avoid warning due to the use of simplified piezo properties.
model=p_piezo('DToSimple',model);
  model.pl=m_elastic('dbval 1 Aluminum');
  d=zeros(1,18); d([11 13])=560e-12; d([3 6])=-185e-12; d(9)=440e-12;
  eps=zeros(1,9); eps([1 5 9])= 8.854e-12*1850;
%                                      elasid |dij coeff  | dielectric coeffs
  model.pl=[ model.pl zeros(1,24);
            3 fe_mat('m_piezo','SI',2) 2        d           eps;
            2 fe_mat('m_elastic','SI',1) 54e9 0.41 7740 zeros(1,25)];
%% Step 4 - Laminate properties and piezo electrodes
model.il=p_shell('dbval 1 laminate 1 1.2e-3 0', ...
  'dbval 2 laminate 3 2.5e-4 0 1 1.2e-3 0 3 2.5e-4 0');
%%% Piezo electrodes
%%%                                         NdNb LayerID   NdNb  LayerID
model.il=p_piezo(model.il,'dbval 3 shell 2 1682    1   0    1683  3 0');
model.il=p_piezo(model.il,'dbval 4 shell 2 1684    1   0    1685  3 0');
cf=feplot(model); fecom('showmap'); fecom('view3');
% scale properly
fecom('scalecoeff 1e-10'); fecom('showmap')
d =  2.5244e-06
%% Step 5 - Compute and display response to static imposed voltage
model=fe_case(model,'FixDof','Cantilever','x==0 -DOF 1:6');
model=fe_case(model,'DofSet','V-Act',struct('def',1,'DOF',1682.21)); %Act
model=p_piezo('ElectrodeSensQ  1683 Q-S1',model);
model=p_piezo('ElectrodeSensQ  1684 Q-S2',model);
model=p_piezo('ElectrodeSensQ  1685 Q-S3',model);
model=fe_case(model,'SensDof','Tip',1054.03); % Displ sensor top right corner
sens=fe_case(model,'sens');
model=fe_case(model,'FixDof','SC*1683-1685',[1682:1685]+.21);
d0=fe_simul('dfrf',stack_set(model,'info','Freq',0)); % direct refer frf at 0Hz
cf=feplot(model,d0); fecom(';view3;colordatagroup;undefline');
cf.mdl.name='Plate_4pzt'; d_piezo('SetStyle',cf); 
d=sens.cta(4,:)*d0.def % Tip displacement is positive.
%% Step 6 - Piezoelectric patches on the bottom only
 model=struct('Node',[1 0 0 0 0 0 0],'Elt',[]);
 model=feutil('addelt',model,'mass1',1);
 % Note that the extrusion values are chosen to include the patch edges
 dx=[linspace(0,15,3) linspace(15,15+55,10) linspace(15+55,463-5,50) 463];
 model=feutil('extrude 0 1 0 0',model,.001*unique(dx));
 dy=[linspace(0,12,3) linspace(12,12+25,5) linspace(12+25,63,5) ...
    linspace(63,63+25,5) linspace(63+25,100,3)];
 model=feutil('extrude 0 0 1 0',model,.001*unique(dy));

% Set patch areas as set different properties
 model.Elt=feutil('divide group 1 withnode{x>.015 & x<.07 & y>.013 & y<.037 }',model);
 model.Elt=feutil('divide group 2 withnode{x>.015 & x<.07 & y>.064 & y<.088 }',model);
 model.Elt=feutil('set group 1 proId3',model);
 model.Elt=feutil('set group 2 proId4',model);

%%  Material Properties
model.pl=m_elastic('dbval 1 Aluminum');
model.pl=m_piezo(model.pl,'dbval 3 -elas2 SONOX_P502_iso');
% To avoid warning due to the use of simplified piezo properties.
model=p_piezo('DToSimple',model);
%%  Laminate properties and piezo electrodes
% This is where an offset must be specified
%(else z0=-7.25e-4 (half of total thickness which is not correct).
model.il=p_shell('dbval 1 laminate 1 1.2e-3 0', ...
   'dbval 2 laminate z0=-8.5e-4 3 2.5e-4 0 1 1.2e-3 0 ');
%%% Piezo electrodes
%%%                                         NdNb LayerID
model.il=p_piezo(model.il,'dbval 3 shell 2 1682    1   0    ');
model.il=p_piezo(model.il,'dbval 4 shell 2 1684    1   0    ');

%%  Static response computation
model=fe_case(model,'FixDof','Cantilever','x==0 -DOF 1:6');
model=fe_case(model,'DofSet','V-Act',struct('def',1,'DOF',1682.21)); %Act
model=p_piezo('ElectrodeSensQ  1684 Q-S2',model);
model=fe_case(model,'SensDof','Tip',1054.03); % Displ sensor top right corner
sens=fe_case(model,'sens');
model=fe_case(model,'FixDof','SC*1682-1684',[1682 1684]+.21);
d0=fe_simul('dfrf',stack_set(model,'info','Freq',0)); % direct refer frf at 0Hz
feplot(model,d0); fecom(';view3;colordatagroup;undefline');
d2=sens.cta(2,:)*d0.def % Tip displacement is positive.
disp(['difference of static response ' num2str((d2-d)/d*100) '%'])
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopzmeshingauto')
%% #TutoPzMeshingAuto : Plate with 4 pzt patches:auto meshing -2

% see sdtweb pzplatemeshing#tutopzmeshingauto
%
%% BeginSource sdtweb('_example','pzplatemeshing.tex#tutopzmeshingauto')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example in d_piezo('ScriptTutoPzMeshingAuto')
d_piezo('DefineStyles');

%% Step 1 : model of host plate -
 model=struct('Node',[1 0 0 0 0 0 0],'Elt',[]);
 model=feutil('addelt',model,'mass1',1);
 % Note that the extrusion values are chosen to include the patch edges
 dx=[linspace(0,15,3) linspace(15,15+55,10) linspace(15+55,463-5,50) 463];
 model=feutil('extrude 0 1 0 0',model,unique(dx));
 dy=[linspace(0,12,3) linspace(12,12+25,5) linspace(12+25,63,5) ...
    linspace(63,63+25,5) linspace(63+25,100,3)];
 model=feutil('extrude 0 0 1 0',model,unique(dy));
model.unit='mm';

% Material Properties
model.pl=m_elastic('dbval 1 Aluminum');
% Laminate properties
model.il=p_shell('dbval 1 -punit mm laminate 1 1.2 0') % this is to specify in mm
model=fe_case(model,'FixDof','Cantilever','x==0 -DOF 1:6');
%%% Step 2: Add patches
RG.list={'Name','Lam','shape'
   'Main_plate', model,''  % Base structure
   'Act1', ... % name of patch
   'BaseId1 -Rect.Sonox_P502_iso.5525TH0_25 +Rect.Sonox_P502_iso.5525TH0_25', ...
     struct('shape','LsRect', ... % Remeshing strategy (lsutil rect here)
       'xc',15+55/2,'yc',12+25/2,'alpha',0,'tolE',.1)
    'Act2', ... % name of patch
   'BaseId1 -Rect.Sonox_P502_iso.5525TH0_25 +Rect.Sonox_P502_iso.5525TH0_25', ...
     struct('shape','LsRect', ... % Remeshing strategy (lsutil rect here)
       'xc',15+55/2,'yc',63+25/2,'alpha',0,'tolE',.1)
       };
mo1=d_piezo('MeshPlate',RG);
mo1=stack_rm(mo1,'info','Electrodes'); % Obsolete stack field to be removed
% To avoid warning due to the use of simplified piezo properties.
mo1=p_piezo('DToSimple',mo1);

%% Step 3 : compute response
nd=feutil('find node x==463 & y==100',model);
elnd=floor(p_piezo('electrodedof.*',mo1)); % Nodes associated to electrodes
mo1=fe_case(mo1,'SensDof','Tip',nd+.03); % Displ sensor
mo1=fe_case(mo1,'DofSet','V-Act',struct('def',1,'DOF',elnd(1)+.21)); %Act
mo1=p_piezo(['ElectrodeSensQ  ' num2str(elnd(1)) ' Q-Act'],mo1); % Charge sensors
mo1=p_piezo(['ElectrodeSensQ ' num2str(elnd(2)) ' Q-S1'],mo1);
mo1=p_piezo(['ElectrodeSensQ ' num2str(elnd(3)) ' Q-S2'],mo1);
mo1=p_piezo(['ElectrodeSensQ ' num2str(elnd(4)) ' Q-S3'],mo1);

% Fix last 3 elec dofs to measure resultant (charge)
mo1=fe_case(mo1,'FixDof', ...
['SC*' num2str(elnd(2)) '/' num2str(elnd(3)) '/' num2str(elnd(4)) ], ...
elnd([2:4])+.21);
sens=fe_case(mo1,'sens');

d1=fe_simul('dfrf',stack_set(mo1,'info','Freq',0)); % direct refer frf at 0Hz
d1t=sens.cta(1,:)*d1.def % Extract tip displ
feplot(mo1,d1);
fecom('colordatapro'); fecom('view3');
%% Step 4 : use local remeshing element size is 7 mm
 model=feutil('objectquad 1 1',[0 0 0;1 0 0;0 1 0], ...
    feutil('refineline 7',[0 463]), ...
    feutil('refineline 7',[0 100]));
     model.unit='mm';
         
% Material Properties
model.pl=m_elastic('dbval 1 Aluminum');
% Laminate properties
model.il=p_shell('dbval 1 -punit mm laminate 1 1.2 0') % this is to specify in mm
model=fe_case(model,'FixDof','Cantilever','x==0 -DOF 1:6');

%% Add patches
RG.list={'Name','Lam','shape'
   'Main_plate', model,''  % Base structure
   'Act1', ... % name of patch
   'BaseId1 -Rect.Sonox_P502_iso.5525TH0_25 +Rect.Sonox_P502_iso.5525TH0_25', ...
     struct('shape','LsRect', ... % Remeshing strategy (lsutil rect here)
       'xc',15+55/2,'yc',12+25/2,'alpha',0,'tolE',.1)
    'Act2', ... % name of patch
   'BaseId1 -Rect.Sonox_P502_iso.5525TH0_25 +Rect.Sonox_P502_iso.5525TH0_25', ...
     struct('shape','LsRect', ... % Remeshing strategy (lsutil rect here)
       'xc',15+55/2,'yc',63+25/2,'alpha',0,'tolE',.1)
       };
mo2=d_piezo('MeshPlate',RG);
mo2=stack_rm(mo2,'info','Electrodes'); % Obsolete stack field to be removed
% To avoid warning due to the use of simplified piezo properties.
mo2=p_piezo('DToSimple',mo2);

nd=feutil('find node x==463 & y==100',model);
elnd=floor(p_piezo('electrodedof.*',mo2)); % Nodes associated to electrodes
mo2=fe_case(mo2,'SensDof','Tip',nd+.03); % Displ sensor
mo2=fe_case(mo2,'DofSet','V-Act',struct('def',1,'DOF',elnd(1)+.21)); %Act
mo2=p_piezo(['ElectrodeSensQ  ' num2str(elnd(1)) ' Q-Act'],mo2); % Charge sensors
mo2=p_piezo(['ElectrodeSensQ ' num2str(elnd(2)) ' Q-S1'],mo2);
mo2=p_piezo(['ElectrodeSensQ ' num2str(elnd(3)) ' Q-S2'],mo2);
mo2=p_piezo(['ElectrodeSensQ ' num2str(elnd(4)) ' Q-S3'],mo2);

% Fix last 3 elec dofs to measure resultant (charge)
mo2=fe_case(mo2,'FixDof', ...
['SC*' num2str(elnd(2)) '/' num2str(elnd(3)) '/' num2str(elnd(4)) ], ...
elnd([2:4])+.21);
sens=fe_case(mo2,'sens');

d2=fe_simul('dfrf',stack_set(mo2,'info','Freq',0)); % direct refer frf at 0Hz
d2t=sens.cta(1,:)*d2.def;

cf=feplot(mo2,d2);
fecom(';colordatapro;view3;undef line'); 
d_piezo('SetStyle',cf); feplot(cf);
[d1t d2t]
%% Step 5 : use a finer mesh to check convergence
ref=[5 3 2];
dt=[d2t];

for ij=1:length(ref)

 model=feutil('objectquad 1 1',[0 0 0;1 0 0;0 1 0], ...
    feutil(['refineline ' num2str(ref(ij))],[0 463]), ...
    feutil(['refineline ' num2str(ref(ij))],[0 100]));
     model.unit='mm';
% Material Properties
model.pl=m_elastic('dbval 1 Aluminum');
% Laminate properties
model.il=p_shell('dbval 1 -punit mm laminate 1 1.2 0') % this is to specify in mm
model=fe_case(model,'FixDof','Cantilever','x==0 -DOF 1:6');

%% Add patches
RG.list={'Name','Lam','shape'
   'Main_plate', model,''  % Base structure
   'Act1', ... % name of patch
   'BaseId1 -Rect.Sonox_P502_iso.5525TH0_25 +Rect.Sonox_P502_iso.5525TH0_25', ...
     struct('shape','LsRect', ... % Remeshing strategy (lsutil rect here)
       'xc',15+55/2,'yc',12+25/2,'alpha',0,'tolE',.1)
    'Act2', ... % name of patch
   'BaseId1 -Rect.Sonox_P502_iso.5525TH0_25 +Rect.Sonox_P502_iso.5525TH0_25', ...
     struct('shape','LsRect', ... % Remeshing strategy (lsutil rect here)
       'xc',15+55/2,'yc',63+25/2,'alpha',0,'tolE',.1)
       };
mo2=d_piezo('MeshPlate',RG);
mo2=stack_rm(mo2,'info','Electrodes'); % Obsolete stack field to be removed
% To avoid warning due to the use of simplified piezo properties.
mo2=p_piezo('DToSimple',mo2);

nd=feutil('find node x==463 & y==100',model);
elnd=floor(p_piezo('electrodedof.*',mo2)); % Nodes associated to electrodes
mo2=fe_case(mo2,'SensDof','Tip',nd+.03); % Displ sensor
mo2=fe_case(mo2,'DofSet','V-Act',struct('def',1,'DOF',elnd(1)+.21)); %Act
mo2=p_piezo(['ElectrodeSensQ  ' num2str(elnd(1)) ' Q-Act'],mo2); % Charge sensors
mo2=p_piezo(['ElectrodeSensQ ' num2str(elnd(2)) ' Q-S1'],mo2);
mo2=p_piezo(['ElectrodeSensQ ' num2str(elnd(3)) ' Q-S2'],mo2);
mo2=p_piezo(['ElectrodeSensQ ' num2str(elnd(4)) ' Q-S3'],mo2);

% Fix last 3 elec dofs to measure resultant (charge)
mo2=fe_case(mo2,'FixDof', ...
['SC*' num2str(elnd(2)) '/' num2str(elnd(3)) '/' num2str(elnd(4)) ], ...
elnd([2:4])+.21);
sens=fe_case(mo2,'sens');

d2=fe_simul('dfrf',stack_set(mo2,'info','Freq',0)); % direct refer frf at 0Hz
d2t=sens.cta(1,:)*d2.def;

dt=[dt; d2t];
end

gf=figure; plot([7 ref],1e6*dt,'linewidth',2); set(gca, 'XDir','reverse');
set(gca,'Fontsize',15); v=get(gca,'XLim'); hold on;
plot(v,1e6*[d1t d1t],'r','linewidth',2);
legend('Local remeshing','Conforming mesh')
xlabel('mesh size (mm)'); ylabel('tip displacement (mm/V)')
%% Step 6 : with circular patches
  model=feutil('objectquad 1 1',[0 0 0;1 0 0;0 1 0], ...
    feutil('refineline 7',[0 463]), ...
    feutil('refineline 7',[0 100]));
     model.unit='mm';

%%%%%% Material Properties
model.pl=m_elastic('dbval 1 Aluminum');
%%%%% Laminate properties
model.il=p_shell('dbval 1 -punit mm laminate 1 1.2 0') % this is to specify in mm

model=fe_case(model,'FixDof','Cantilever','x==0 -DOF 1:6');
RG.list={'Name','Lam','shape'
   'Main_plate', model,''  % Base structure
   'Act1', ... % name of patch
   'BaseId1 -Disk.Sonox_P502_iso.RC10TH0_25 +Disk.Sonox_P502_iso.RC10TH0_25', ...
   struct('shape','lscirc','xc',15+55/2,'yc',12+25/2),
   'Act2', ... % name of patch
   'BaseId1 -Disk.Sonox_P502_iso.RC10TH0_25 +Disk.Sonox_P502_iso.RC10TH0_25', ...
   struct('shape','lscirc','xc',15+55/2,'yc',63+25/2)};
%
mo3=d_piezo('MeshPlate',RG);
mo3=stack_rm(mo3,'info','Electrodes'); 
% Old Stack not necessary or should be set to 0
mo3.pl([3 5 7 9],7)=0; 
% Set damping to zero in Noliac otherwise complex static response
% To avoid warning due to the use of simplified piezo properties.
mo3=p_piezo('DToSimple',mo3);

nd=feutil('find node x==463 & y==100',model);
elnd=floor(p_piezo('electrodedof.*',mo3)); % Nodes associated to electrodes
mo3=fe_case(mo3,'SensDof','Tip',nd+.03); % Displ sensor
mo3=fe_case(mo3,'DofSet','V-Act',struct('def',1,'DOF',elnd(1)+.21)); %Act
mo3=p_piezo(['ElectrodeSensQ  ' num2str(elnd(1)) ' Q-Act'],mo3); % Charge sensors
mo3=p_piezo(['ElectrodeSensQ ' num2str(elnd(2)) ' Q-S1'],mo3);
mo3=p_piezo(['ElectrodeSensQ ' num2str(elnd(3)) ' Q-S2'],mo3);
mo3=p_piezo(['ElectrodeSensQ ' num2str(elnd(4)) ' Q-S3'],mo3);

% Fix last 3 elec dofs to measure resultant (charge)
mo3=fe_case(mo3,'FixDof', ...
['SC*' num2str(elnd(2)) '/' num2str(elnd(3)) '/' num2str(elnd(4)) ], ...
elnd([2:4])+.21);
sens=fe_case(mo3,'sens');

d3=fe_simul('dfrf',stack_set(mo3,'info','Freq',0)); % direct refer frf at 0Hz
d3t=sens.cta(1,:)*d3.def; % First electrode is on top now ?

feplot(mo3,d3)
fecom(';colordatapro;view3;undef line')
d_piezo('setstyle',cf);
iimouse('view',gca,[ -1499 -1955 1462 -248.6 -325.7 ...
 275.9 0.30 0.40 0.87 2.20]); % obtained with iimouse('cv')
%% Step 7 : Check convergence

ref=[5 3 2];
dt=[d3t];

for ij=1:length(ref)

 model=feutil('objectquad 1 1',[0 0 0;1 0 0;0 1 0], ...
    feutil(['refineline ' num2str(ref(ij))],[0 463]), ...
    feutil(['refineline ' num2str(ref(ij))],[0 100]));
     model.unit='mm';
% Material Properties
model.pl=m_elastic('dbval 1 Aluminum');
% Laminate properties
model.il=p_shell('dbval 1 -punit mm laminate 1 1.2 0') % this is to specify in mm
model=fe_case(model,'FixDof','Cantilever','x==0 -DOF 1:6');

RG.list={'Name','Lam','shape'
   'Main_plate', model,''  % Base structure
   'Act1', ... % name of patch
   'BaseId1 -Disk.Sonox_P502_iso.RC10TH0_25 +Disk.Sonox_P502_iso.RC10TH0_25', ...
   struct('shape','lscirc','xc',15+55/2,'yc',12+25/2),
   'Act2', ... % name of patch
   'BaseId1 -Disk.Sonox_P502_iso.RC10TH0_25 +Disk.Sonox_P502_iso.RC10TH0_25', ...
   struct('shape','lscirc','xc',15+55/2,'yc',63+25/2)};

%
mo3=d_piezo('MeshPlate',RG);
mo3=stack_rm(mo3,'info','Electrodes');
% Old Stack not necessary or should be set to 0
mo3.pl([3 5 7 9],7)=0; 
% Set damping to zero in Noliac otherwise complex static response
% To avoid warning due to the use of simplified piezo properties.
mo3=p_piezo('DToSimple',mo3);

nd=feutil('find node x==463 & y==100',model);
elnd=floor(p_piezo('electrodedof.*',mo3)); % Nodes associated to electrodes
mo3=fe_case(mo3,'SensDof','Tip',nd+.03); % Displ sensor
mo3=fe_case(mo3,'DofSet','V-Act',struct('def',1,'DOF',elnd(1)+.21)); %Act
mo3=p_piezo(['ElectrodeSensQ  ' num2str(elnd(1)) ' Q-Act'],mo3); % Charge sensors
mo3=p_piezo(['ElectrodeSensQ ' num2str(elnd(2)) ' Q-S1'],mo3);
mo3=p_piezo(['ElectrodeSensQ ' num2str(elnd(3)) ' Q-S2'],mo3);
mo3=p_piezo(['ElectrodeSensQ ' num2str(elnd(4)) ' Q-S3'],mo3);

% Fix last 3 elec dofs to measure resultant (charge)
mo3=fe_case(mo3,'FixDof', ...
['SC*' num2str(elnd(2)) '/' num2str(elnd(3)) '/' num2str(elnd(4)) ], ...
elnd([2:4])+.21);
sens=fe_case(mo3,'sens');

d3=fe_simul('dfrf',stack_set(mo3,'info','Freq',0)); % direct refer frf at 0Hz
d3t=sens.cta(1,:)*d3.def; % First electrode is on top now ?

dt=[dt; d3t];
end

gf=figure; plot([7 ref],1e6*(dt),'linewidth',2); set(gca, 'XDir','reverse');
set(gca,'Fontsize',15)
xlabel('mesh size (mm)'); ylabel('tip displacement (mm/V)')
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopzmeshingmfc')
%% #TutoPzMeshingMFC : Plate with MFCs : meshing -2

% see sdtweb pzplatemeshing#tutopzmeshingmfc
%
%% BeginSource sdtweb('_example','pzplatemeshing.tex#tutopzmeshingmfc')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example in d_piezo('ScriptTutoPzMeshingMFC')
d_piezo('DefineStyles');

%% Step 1 : Create mesh. Geometric properties in the manual
RO=struct('L',463,'w',50,'a',85,'b',28,'c',15,'d',11);

% create a rectangle with targetl = 3 mm
 model=feutil('objectquad 1 1',[0 0 0;1 0 0;0 1 0], ...
    feutil('refineline 5',[0 RO.c+[0 RO.a] RO.L]), ...
    feutil('refineline 5',[0 RO.d+[0 RO.b] RO.w]));
%%%%% Material Properties for supporting plate
 model.pl=m_elastic('dbval 1 -unit MM Aluminum'); % Aluminum
 model.il=p_shell('dbval 1 -punit MM laminate 1 1 0');
 model.unit='MM';

RG.list={'Name','Lam','shape'
   'Main_plate', model,''  % Base structure
   'Act1', ... % name of patch
   'BaseId1 +SmartM.MFC-P1.8528 -SmartM.MFC-P1.8528', ... % Layout definition
     struct('shape','LsRect', ... % Remeshing strategy (lsutil rect here)
       'xc',RO.c+RO.a/2,'yc',RO.d+RO.b/2,'alpha',0,'tolE',.1)
       };
mo1=d_piezo('MeshPlate',RG); 
cf=feplot(mo1); fecom(';colordatapro;view3');
cf.mdl.name='MFC plate mesh'; % Model name for title
d_piezo('setstyle',cf)

%% EndSource EndTuto

elseif comstr(Cam,'tutopztowerred')
%% #TutoPzTowerRed : Concrete tower : reduced models -2

% see sdtweb reduction#tutopztowerred
%
%% BeginSource sdtweb('_example','reduction.tex#tutopztowerred')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzTowerRed')
d_piezo('DefineStyles');

%% Step 1 : Build the model and define actuator and sensor
model=d_piezo('MeshTower');
%% Step 2 : Build Reduced basis (3 modes) and compute TF

% M and K matrices have been built with fe_mknl in d_piezo('MeshTower')
M=model.K{1}; K=model.K{2};
% Excitation
Load=fe_load(model); b=Load.def;

% Projection matrix for output
sens=fe_case(model,'sens');
[Case,model.DOF]=fe_mknl('init',model); % Build Case.T for active dofs
cta=sens.cta*Case.T;

% Compute Modeshapes
def=fe_eig(model);
% Rearrange with active dofs only
def.def=Case.T'*def.def; def.DOF=Case.T'*def.DOF;

% Build reduced basis on three modes
T=def.def(:,1:3);

% Reduce matrices and compute response - explicit computation
Kr=T'*K*T; Mr=T'*M*T; br=T'*b;

% Define frequency vector for computations
w=linspace(0,30*2*pi,512);% Extended frequency range

% Solution with full and reduced matrices - explicit computation
% Use loss factor =0.02 = default for SDT (1% modal damping)
for i1=1:length(w)
    U(:,i1)=(K*(1+0.02*1i)-w(i1)^2*M)\b; 
    %Default loss factor is 0.02 in SDT
    Ur(:,i1)=(Kr*(1+0.02*1i)-w(i1)^2*Mr)\br; %Reduced basis
end

% Extract response on sensor and visualize in iicom
out1=cta*U; out2=cta*(T*Ur);

% Change output format to be compatible with iicom
C0m=d_piezo('BuildC1',w'/(2*pi),out1.','d-top','F-top'); C0m.name='Full';
C1m=d_piezo('BuildC1',w'/(2*pi),out2.','d-top','F-top'); C1m.name='3md';
ci=iiplot;iicom(ci,'curveinit',{'curve',C0m.name,C0m;'curve',C1m.name,C1m}); 
iicom('submagpha')
d_piezo('setstyle',ci); set(gca,'XLim',[0 10])
%% Step 3 - Compute with fe_simul Full/3md

% Full model response
C0=fe_simul('dfrf -sens',stack_set(model,'info','Freq',w/(2*pi)));
% -sens option projects result on sensors only
C0.name='Full'; C0.X{2}={'d-top'}; C1.X{3}={'F-top'};

% With reduced basis 3 modes
model = stack_set(model,'info','EigOpt',[5 3 0]); % To keep 3 modes
% Build super-element with 3 modes
SE1=fe_reduc('call d_piezo@modal -matdes 2 1 3 4',model); 

% Make model with a single super-element
SE0 = struct('Node',[],'Elt',[]);
mo1 = fesuper('SEAdd 1 -1 -unique -initcoef -newID se1',SE0,SE1) ;

% Define input/output
mo1=fe_case(mo1,'SensDOF','Output',21.01);
mo1=fe_case(mo1,'DofLoad','Input',21+.01,1);

% Compute response with fe_simul 
C1=fe_simul('dfrf -sens',stack_set(mo1,'info','Freq',w/(2*pi))); 
C1.name='3md+SE'; C1.X{2}={'d-top'};
ci=iiplot;iicom(ci,'curveinit',{'curve',C0.name,C0;'curve',C1.name,C1}); 
iicom('submagpha')
d_piezo('setstyle',ci);
set(gca,'XLim',[0 10])
%% Step 4 : Reduced basis 2 : 3 modes + static corr - explicit 

Tstat=K\b;
T=fe_norm([def.def(:,1:3) Tstat],M,K); % Orthonormalize vectors

% Reduced matrices
Kr=T'*K*T; Mr=T'*M*T; br=T'*b;

% Ref solution
for i1=1:length(w)
     Ur2(:,i1)=(Kr*(1+0.02*1i)-w(i1)^2*Mr)\br;
end

out3=cta*(T*Ur2);

C2m=d_piezo('BuildC1',w'/(2*pi),out3.','d-top','F-top'); C2m.name='3md+Stat';
iicom(ci,'curveinit',{'curve',C0m.name,C0m;'curve',C2m.name,C2m}); 
iicom('submagpha'); d_piezo('setstyle',ci);
%% Step 5 with fe_simul
model=d_piezo('meshtower');
model = stack_set(model,'info','EigOpt',[5 3 0]); % To keep 3 modes

% Reduce model with 3 modes and static correction
mo1=fe_reduc('free -matdes -SE 2 1 3 4 -bset',model); 
% Compute response with fe_simul and represent
C2=fe_simul('dfrf -sens',stack_set(mo1,'info','Freq',w/(2*pi))); % 
C2.name='3md+Stat-SE'; C2.X{2}={'d-top'};
ci=iiplot;iicom(ci,'curveinit',{'curve',C0.name,C0;'curve',C2.name,C2}); 
iicom('submagpha'); d_piezo('setstyle',ci);
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopztowerss')
%% #TutoPzTowerSS : Concrete tower : reduced state-space models -2

% see sdtweb state_space#tutopztowerss
%
%% BeginSource sdtweb('_example','state_space.tex#tutopztowerss')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzTowerSS')
d_piezo('DefineStyles');

%% Step 1 : Build the model and define actuator and sensor
model=d_piezo('MeshTower');
% Step 2 : State-space models
[sys,TR] = fe2ss('free 5 3 0 -dterm',model);
[sys2,TR2] = fe2ss('free 5 3 0 ',model);

w=linspace(0,30*2*pi,512);% Frequency range
% Convert to curve object and relabel X (fe2ss uses dofs and not sens/act names)
C1=qbode(sys,w,'struct-lab');C1.name='SS-dterm'; C1.X{2}={'d-top'}; 
C2=qbode(sys2,w,'struct-lab');C2.name='SS-mode'; C2.X{2}={'d-top'}; 
ci=iiplot;
iicom(ci,'curveinit',{'curve',C1.name,C1;'curve',C2.name,C2}); iicom('submagpha')
d_piezo('setstyle',ci);
%% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopztowerssuimp')
%% #TutoPzTowerSSUimp : Concrete tower : reduced ss models Uimp -2

% see sdtweb state_space#tutopztowerssuimp
%
%% BeginSource sdtweb('_example','state_space.tex#tutopztowerssuimp')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzTowerSSUimp')
d_piezo('DefineStyles');

%% Step 1 : Build the model and define actuator and sensor
model=d_piezo('MeshTower');
model=d_piezo('meshtower');
model=fe_case(model,'FixDof','Clamped',[1.06]); % Leave x free for imposed displ
model=fe_case(model,'Remove','F-top'); % Remove point force
model=fe_case(model,'DOFSet','Uimp',[1.01]); % Imposed horizontal displ
%% Step 2 : Reference method - exact solution + Inertial term neglected - full model
% Build matrices
[model,Case] = fe_case('assemble NoT -matdes 2 1 Case -SE',model) ;

% Full model
K0 = feutilb('tkt',Case.T,model.K); % Assemble matrices taking into account BCs
F1 = -Case.T'*model.K{2}*Case.TIn; % Loading due to imposed displacement
F2 = -Case.T'*model.K{1}*Case.TIn; % Inertial term
%
% compute response in freq domain
w=linspace(0,100,512);
for i=1:length(w)
    U0(:,i)=Case.T*((K0{2}*(1+0.02*1i)-w(i)^2*K0{1})\(F1-w(i)^2*F2)); 
    % Take into account mass term
    U1(:,i)=Case.T*((K0{2}*(1+0.02*1i)-w(i)^2*K0{1})\F1); % Stiffness term only
end
CTA = fe_c(model.DOF,21.01); u0 = CTA*U0; u1 = CTA*U1;

% Change output format to be compatible with iicom
C1=d_piezo('BuildC1',w'/(2*pi),u0.','d-top','Uimp'); C1.name='M and K';
C2=d_piezo('BuildC1',w'/(2*pi),u1.','d-top','Uimp'); C2.name='K only';
ci=iiplot;iicom(ci,'curveinit',{'curve',C1.name,C1;'curve',C2.name,C2}); 
iicom('submagpha')
d_piezo('setstyle',ci);
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopztowerssaimp')
%% #TutoPzTowerSSAimp : Concrete tower : reduced ss models Aimp -2

% see sdtweb state_space#tutopztowerssaimp
%
%% BeginSource sdtweb('_example','state_space.tex#tutopztowerssaimp')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzTowerSSAimp')
d_piezo('DefineStyles');

%% Step 1 : Build the model and define actuator and sensor
model=d_piezo('MeshTower');

model=fe_case(model,'FixDof','Clamped',[1.06]); % Leave x free for imposed displ
model=fe_case(model,'Remove','F-top'); % Remove point force
model=fe_case(model,'DOFSet','UImp',[1.01]); % Set an imposed displacement
%% Step 2 : Regular method with RHS M*Tin (relative displacement)
% --------- full model
[model,Case] = fe_case('assemble NoT -matdes 2 1 Case -SE',model) ;

K0 = feutilb('tkt',Case.T,model.K); % Assemble matrices taking into account BCs
% Compute TIn as static response to imposed displacement
TIn=fe_simul('static',model); TIn=TIn.def; 
% Loading due to imposed acceleration
F = -Case.T'*model.K{1}*TIn; 

% compute response in freq domain (relative displacement)
w=linspace(1,100,512);
for i=1:length(w)
    U0r(:,i)=Case.T*((K0{2}*(1+0.02*1i)-w(i)^2*K0{1})\F);
end
Uimp=1./(-w.^2); % Imposed displacement for a unit imposed acceleration

% Absolute displacement
for i=1:length(w)
U0(:,i)=U0r(:,i)+Uimp(i)*TIn;
end

CTA = fe_c(model.DOF,21.01); u0 = CTA*U0; u0r = CTA*U0r;

% Change output format to be compatible with iicom
C0=d_piezo('BuildC1',w'/(2*pi),u0.','d-top','Aimp'); C0.name='Full';
C1=d_piezo('BuildC1',w'/(2*pi),u0r.','dr-top','Aimp'); C1.name='Full';
%% Step 3 - State-space with SDT(relative) using a DofLoad
model2=d_piezo('Meshtower');
model2=fe_case(model2,'FixDof','Clamped',[1.01 1.06]); % Block all interface dofs
model2=fe_case(model2,'Remove','F-top'); % Remove point force
[model2,Case2] = fe_case('assemble NoT -matdes 2 1 Case -SE',model2) ;

% Load when using relative displacements
SET.DOF=model2.DOF; SET.def=Case2.T*F; 
model2=fe_case(model2,'DofLoad','Aimp',SET); % 

% state-space model
sysr=fe2ss('free 5 5 0 -dterm',model2);
C2=qbode(sysr,w,'struct-lab');C2.name='fe2ss-5md+st'; C2.X{2}={'dr-top'}; 
ci=iiplot;iicom(ci,'curveinit',{'curve',C1.name,C1;'curve',C2.name,C2}); 
iicom('submagpha')
d_piezo('setstyle',ci);
%% Step 4 : state-space model for absolute displacements - Aimp

% This is a CB basis which is renormalized (so free BCs and rigid body mode)
% TR2.data is needed for nor2ss hence the normalization.
TR2 = fe2ss('craigbampton 5 5 -basis',model); 
sysu= nor2ss(TR2,model) ;
C3=qbode(sysu,w,'struct-lab');C3.name='fereduc+nor2ss'; 
C3.X{2}={'d-top'}; C3.X{3}={'Aimp'}
ci=iiplot;iicom(ci,'curveinit',{'curve',C0.name,C0;'curve',C3.name,C3}); 
iicom('submagpha');d_piezo('setstyle',ci); 
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopzplate4pztss')
%% #TutoPzPlate4pztSS : Plate with 4 pzt patches: ss model -2

% see sdtweb state_space#tutopzplate4pztss
%
%% BeginSource sdtweb('_example','state_space.tex#tutopzplate4pztss')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example in d_piezo('ScriptTutoPzPlate4Pztss')
d_piezo('DefineStyles');

%% Step 1 - Build model and visualize
model=d_piezo('MeshULBplate');  % creates the model
model=fe_case(model,'FixDof','Cantilever','x==0'); % Clamp plate
% Set modal default zeta = 0.01
model=stack_set(model,'info','DefaultZeta',0.01);
%% Step 2 - Define actuators and sensors and visualize
nd=feutil('find node x==463 & y==100',model);
model=fe_case(model,'SensDof','Tip',{[num2str(nd) ':z']}); % Displ sensor
i1=p_piezo('TabInfo',model);i1=i1.Electrodes(:,1);
model=fe_case(model,'DofSet','V-Act',struct('def',1,'DOF',i1(1)+.21, ...%Act
    'Elt',feutil('selelt proid 104',model))); % Elt defined for display
model=p_piezo(sprintf('ElectrodeSensQ  %i Q-Act',i1(1)),model); % Charge sensors
model=p_piezo(sprintf('ElectrodeSensQ  %i Q-S1',i1(2)),model);
model=p_piezo(sprintf('ElectrodeSensQ  %i Q-S2',i1(3)),model);
model=p_piezo(sprintf('ElectrodeSensQ  %i Q-S3',i1(4)),model);
% Fix ElectrodeSensQ dofs to measure resultant (charge)
model=fe_case(model,'FixDof','SC*S1-S3',i1(2:end)+.21);
%% Step 2 Compute dynamic response full/state-space and compare
model=stack_set(model,'info','oProp',mklserv_utils('oprop','CpxSym'));
f=linspace(1,100,400); % in Hz

% Full model
C1=fe_simul('dfrf -sens',stack_set(model,'info','Freq',f(:))); % direct refer frf
C1.X{2}(1)={'Tip'};C1.name='Full';

% state-space model
[sys,TR1]=fe2ss('free 5 10 0 -dterm',model); %
C2=qbode(sys,f(:)*2*pi,'struct-lab');C2.name='SS 10 modes+static';
C2.X{2}(1)={'Tip'};

% Compare the two curves
ci=iiplot;
iicom(ci,'curveinit',{'curve',C1.name,C1;'curve',C2.name,C2});
iicom('submagpha'); d_piezo('setstyles',ci)
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopzplate4pztsscomb')
%% #TutoPzPlate4pztSS : Plate with 4 pzt patches: ss model a&s combi -2

% see sdtweb state_space#tutopzplate4pztsscomb
%
%% BeginSource sdtweb('_example','state_space.tex#tutopzplate4pztsscomb')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example in d_piezo('ScriptTutoPzPlate4PztSSComb')
d_piezo('DefineStyles');

%% Step 1 - Build model and define actuator combinations + sttatic response
model=d_piezo('MeshULBplate cantilever');  % creates the model
model=stack_set(model,'info','DefaultZeta',0.01); % Set modal damping zeta = 0.01
edofs=p_piezo('electrodeDOF.*',model);
model=fe_case(model,'DofSet','V1+2', ...
    struct('def',[1;-1],'DOF',edofs(1:2)));
% Static response
d0=fe_simul('dfrf',stack_set(model,'info','Freq',0)); % direct refer frf
cf=feplot(model); cf.def=d0;
fecom(';view3;scd 20;colordatagroup;undefline')
cf.mdl.name='Plate_4pzt_Comb_static'; d_piezo('SetStyle',cf); feplot(cf);
%% Step 2 - Define sensor combinations
% Combined charge output (SC electrodes) % difference of charge 1684-1685
r1=struct('cta',[1 -1],'DOF',edofs(3:4),'name','QS3+4');
model=p_piezo('ElectrodeSensQ',model,r1);
% Combined voltage output (OC electrodes) % difference of voltage 1684-1685
r1=struct('cta',[1 -1],'DOF',edofs(3:4),'name','VS3+4');
model=fe_case(model,'SensDof',r1.name,r1);
model=fe_case(model,'pcond','Piezo','d_piezo(''Pcond 1e8'')');
%% Step 3 - Compute dynamic response with state-space model
[sys,TR]=fe2ss('free 5 10 0 -dterm',model);
w=linspace(1,100,400)*2*pi;
C1=qbode(sys,w,'struct-lab'); C1.name='OC';

% Now you need to SC 1057 and 1058 to measure charge resultant
model=fe_case(model,'FixDof','SC*3-4',edofs(3:4));
[sys2,TR2]=fe2ss('free 5 10 0 -dterm',model);
C2=qbode(sys2,w,'struct-lab');C2.name='SC';

% Flip channels for C1
C1.Y=fliplr(C1.Y); C1.X{2}(1)={'VS3+4'}; C1.X{2}(2)={'QS3+4'}

% Scale to compare
C2.Y(:,1)=C2.Y(:,1)*C1.Y(1,1)/C2.Y(1,1);
ci=iiplot;iicom('curvereset');
iicom('curveinit',{'curve',C1.name,C1;'curve',C2.name,C2 });
iicom('submagpha'); d_piezo('setstyles',ci)
%% Step 5 - Compute OC and SC frequencies
model=d_piezo('MeshULBplate -cantilever');
% Open circuit : do nothing on electrodes
d1=fe_eig(model,[5 20 1e3]);
% Short circuit : fix all electric DOFs
DOF=p_piezo('electrodeDOF.*',model);
d2=fe_eig(fe_case(model,'FixDof','SC',DOF),[5 20 1e3]);
r1=[d1.data(1:end)./d2.data(1:end)];
gf=figure;plot(r1,'*','linewidth',2);axis tight; set(gca,'Fontsize',15)
xlabel('Mode number');ylabel('f_{OC}/f_{SC}');

% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopztowerssuimpcb')
%% #TutoPzTowerSSUimpCB : Concrete tower : CB reduced ss models Uimp -2

% see sdtweb state_space_CB#tutopztowerssuimpcb
%
%% BeginSource sdtweb('_example','state_space_CB.tex#tutopztowerssuimpcb')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzTowerSSUimpCB')
d_piezo('DefineStyles');

%% Step 1 : Build Model, CB reduction, explicit response
model=d_piezo('meshtower');
model=fe_case(model,'FixDof','Clamped',[1.06]); % Leave x free for imposed displ
model=fe_case(model,'Remove','F-top'); % Remove point force
model=fe_case(model,'DOFSet','UImp',[1.01]); % Set imposed displ

%Full model matrices
[model,Case] = fe_case('assemble NoT -matdes 2 1 Case -SE',model) ;
K0 = feutilb('tkt',Case.T,model.K); % Full-model

% CB matrices and T matrix
CB    = fe_reduc('craigbampton 5 5',model); TR = CB.TR; KCB=CB.K;

% Build loads full and reduced models
F = -Case.T'*model.K{2}*Case.TIn;
F1=-CB.K{2}(2:end,1);
F2=-CB.K{1}(2:end,1);

w=linspace(1,100,512);

for i=1:length(w)
 U0(:,i)=Case.T*((K0{2}*(1+0.02*1i)-w(i)^2*K0{1})\F); % full model
 U1r(:,i)=((KCB{2}(2:end,2:end)*(1+0.02*1i)-w(i)^2*KCB{1}(2:end,2:end)) ...
 	\(F1-w(i)^2*F2)); % CB stiffness and inertia
 U2r(:,i)=((KCB{2}(2:end,2:end)*(1+0.02*1i)-w(i)^2*KCB{1}(2:end,2:end)) ... 
    \(F1)); % CB stiffness
 U3r(:,i)=((KCB{2}(2:end,2:end)*(1+0.02*1i)-w(i)^2*KCB{1}(2:end,2:end))... 
    \(-w(i)^2*F2)); %CB inertia
end

CTA = fe_c(model.DOF,21.01); u0 = CTA*U0;
U1=TR.def*[ones(1,length(w)); U1r]; u1 = CTA*U1;
U2=TR.def*[ones(1,length(w)); U2r]; u2 = CTA*U2;
U3=TR.def*[ones(1,length(w)); U3r]; u3 = CTA*U3;

% Change output format to be compatible with iicom
C1=d_piezo('BuildC1',w'/(2*pi),u0.','d-top','Uimp'); C1.name='full model';
C2=d_piezo('BuildC1',w'/(2*pi),u1.','d-top','Uimp'); C2.name='CB K and M';
C3=d_piezo('BuildC1',w'/(2*pi),u2.','d-top','Uimp'); C3.name='CB K';
C4=d_piezo('BuildC1',w'/(2*pi),u3.','d-top','Uimp'); C4.name='CB M';
ci=iiplot;iicom(ci,'curveinit',{'curve',C1.name,C1;'curve',C2.name,C2; ...
	'curve',C3.name,C3;'curve',C4.name,C4}); 
iicom('submagpha')
d_piezo('setstyle',ci);
%% Step 2: keep bottom and top translation in the CB basis

% Ss from full model - 5 modes + static
sys0=fe2ss('free 5 5 0 -dterm',model);
C0=qbode(sys0,w,'struct-lab');C0.name='fe2ss'; C0.X{2}={'d-top'}

% Top and bottom DOF to be kept in CB reduction
SET.DOF=[1.01; 21.01]; SET.def=eye(2); 
model=fe_case(model,'DOFSet','UImp',SET); 

% Build CB reduced model
model = stack_set(model,'info','EigOpt',[5 5 0]);
SE1= fe_reduc('CraigBampton -SE -matdes 2 1 3 4 ',model); % 
% To keep only the matrices
SE1 = rmfield(SE1,{'il','pl','Stack','TR','mdof'}) ; 
SE1.Node=feutil('getnode',SE1,unique(fix(SE1.DOF)));SE1.Elt=[];

% Define input/output
SE1=fe_case(SE1,'SensDOF','d-top',21.01);
SE1=fe_case(SE1,'DOFSet','UImp',[1.01]); % 

% Build state-space based on CB reduced model
sys=fe2ss('free 5 5 0 -dterm',SE1);
C5=qbode(sys,w,'struct-lab');C5.name='fe2ss CB'; 
C5.X{2}={'d-top'}; C5.X{3}={'Uimp'};
ci=iiplot;iicom(ci,'curveinit',{'curve',C0.name,C0;'curve',C5.name,C5}); 
iicom('submagpha')
d_piezo('setstyle',ci);
%% Step 3 : Apply Raze transform before making the state-space model

% Transform the initial CB matrices
KCB=SE1.K;
N=size(KCB{1},1);
T1=[ eye(2) zeros(2,N-2) ; - KCB{1}(3:end,3:end)\KCB{1}(3:end,1:2) eye(N-2)];

Mr=T1'*(KCB{1})*T1;
Kr=T1'*(KCB{2})*T1;

% Replace matrices in the CB model
SE2=SE1;
SE2.K{1}=Mr;
SE2.K{2}=Kr;

sys2=fe2ss('free 5 5 0 -dterm',SE2); 

% Make CB matrices with Raze's transform directly -noMCI
SE3= fe_reduc('CraigBampton -SE -matdes 2 1 3 4 -noMCI',model);

% To keep only the matrices
SE3 = rmfield(SE3,{'il','pl','Stack','TR','mdof'}) ; 
SE3.Node=feutil('getnode',SE3,unique(fix(SE3.DOF)));SE3.Elt=[];

% Define input/output
SE3=fe_case(SE3,'SensDOF','d-top',21.01);
SE3=fe_case(SE3,'DOFSet','UImp',[1.01]); % 

% Build state-space based on CB reduced model
sys3=fe2ss('free 5 5 0 -dterm',SE3);
C6=qbode(sys2,w,'struct-lab');C6.name='fe2ss CB+Raze'; 
C6.X{2}={'d-top'}; C6.X{3}={'Uimp'};
C7=qbode(sys3,w,'struct-lab');C7.name='fe2ss CB+Raze-direct'; 
C7.X{2}={'d-top'}; C7.X{3}={'Uimp'};

ci=iiplot;iicom(ci,'curveinit',{'curve',C0.name,C0;'curve',C6.name,C6;'curve',C7.name,C7}); 

iicom('submagpha')
d_piezo('setstyle',ci);
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopztowerssaimpcb')
%% #TutoPzTowerSSAimpCB : Concrete tower : CB reduced ss models Aimp -2

% see sdtweb state_space_CB#tutopztowerssaimpcb
%
%% BeginSource sdtweb('_example','state_space_CB.tex#tutopztowerssaimpcb')

% Init working directory for figure generation
d_piezo('SetPlotwd');
%% Step 1: reference solution -direct computation

% Build model
model=d_piezo('meshtower');
model=fe_case(model,'FixDof','Clamped',[1.06]); 
% Leave x free for imposed displ
model=fe_case(model,'Remove','F-top'); % Remove point force

% Initial matrices without imposed displacements
[model0,Case0] = fe_case('assemble NoT -matdes 2 1 Case -SE',model) ;
K0 = feutilb('tkt',Case0.T,model0.K); % Full-model 41x41

% Matrices with Block displ at basis
model=fe_case(model,'DOFSet','UImp',[1.01]); % 

% Full model matrices
[model,Case] = fe_case('assemble NoT -matdes 2 1 Case -SE',model) ;
Kcc = feutilb('tkt',Case.T,model.K); % Blocked at basis 40x40

% Put manually unitary values to be exact
TIn=zeros(40,1); TIn(1:2:end)=1;

% The load is built after applying the boundary condition to the matrices.
F = -(Kcc{1}*TIn);

w=linspace(0,100,512);
eta=0.02;

% Reference solution - full model
for i=1:length(w)
    U0(:,i)=Case.T*((Kcc{2}*(1+eta*1i)-w(i)^2*Kcc{1})\F);
end

CTA = fe_c(model.DOF,21.01); u0 = CTA*U0;
% Change output format to be compatible with iicom
C0=d_piezo('BuildC1',w'/(2*pi),u0.','dr-top','Aimp'); C0.name='Full';
%% Step 2 : state-space model based on CB reduction

% CB reduction on two DOFs
SET.DOF=[1.01; 21.01]; SET.def=eye(2); 
model=fe_case(model,'DOFSet','UImp',SET); 

% Build CB reduced model
model = stack_set(model,'info','EigOpt',[5 5 0]); % No mass shift!
SE1= fe_reduc('CraigBampton -SE -matdes 2 1 3 4 -bset ',model); % 

SE0=SE1; % Save superelement model

% To keep only the matrices
SE1 = rmfield(SE1,{'il','pl','Stack','TR','mdof'}) ; 
SE1.Node=feutil('getnode',SE1,unique(fix(SE1.DOF)));SE1.Elt=[];

% Force vector using CB matrices
KCB=SE1.K;
TInCB= [1; -KCB{2}(2:end,2:end)\KCB{2}(2:end,1)];
FCB= -KCB{1}*TInCB; % From CB matrices only


% Block 1.01 and define a DofLoad and a new observation matrix instead
SE1=fe_case(SE1,'FixDOF','BC',1.01); % BC for relative displacement
SET.DOF=SE1.DOF(2:end); SET.def= FCB(2:end);
SE1=fe_case(SE1,'DOFLoad','UImp',SET); % Equivalent load
SE1=fe_case(SE1,'SensDOF','Output',SET.DOF(1));

sys=fe2ss('free 5 5 0 -dterm',SE1);
C1=qbode(sys,w,'struct-lab');C1.name='fe2ss 5md+st';
C1.X{2}={'dr-top'}; C1.X{3}={'Aimp'};
ci=iiplot;iicom(ci,'curveinit',{'curve',C0.name,C0;'curve',C1.name,C1}); 
iicom('submagpha')
d_piezo('setstyle',ci); 
%% Step 3: State-space model with absolute displacement

SE2=SE0; % Initial model without BCs
TR2=fe_eig(SE2,[5 7 0]); % Compute all7  modes of reduced model

% Imposed acceleration at base
SET.DOF=[1.01]; SET.def=eye(1); %
SE2=fe_case(SE2,'DOFSet','UImp',SET); % Impose acc at bottom

% build state-space model with nor2ss
sys2= nor2ss(TR2,SE2) ;
C2=qbode(sys2,w,'struct-lab');C2.name='nor2ss';
C2.X{2}={'d-top'};  C2.X{3}={'Aimp'}; 

% convert reference solution to absolute displacement
u0a=u0-1./w.^2;
C0a=d_piezo('BuildC1',w'/(2*pi),u0a.','d-top','Aimp'); C0a.name='Full';
ci=iiplot;iicom(ci,'curveinit',{'curve',C0a.name,C0a;'curve',C2.name,C2}); 
iicom('submagpha')
d_piezo('setstyle',ci); 

%% EndSource EndTuto

elseif comstr(Cam,'tutopzplate4pztsscb')
%% #TutoPzPlate4pztSSCB : Plate with 4 pzt patches: CB ss model -2

% see sdtweb state_space_CB_piezo#tutopzplate4pztsscb
%
%% BeginSource sdtweb('_example','state_space_CB_piezo.tex#tutopzplate4pztsscb')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzPlate4pztSSCB')
d_piezo('DefineStyles');

%% Step 1 - Reference solution with fe2ss
model=d_piezo('MeshULBplate');  % creates the model
model=fe_case(model,'FixDof','Cantilever','x==0'); % Clamp plate
% Set modal default zeta = 0.01
model=stack_set(model,'info','DefaultZeta',0.01);

% Define sensors and actuators
nd=feutil('find node x==463 & y==100',model); % Tip node 
i1=p_piezo('TabInfo',model);i1=i1.Electrodes(:,1); % Electrode info
scell={'Tip:z{unit,*1e3=\mu m}'; ...
    'Act:q{unit,*1e12=pC}';'S1:q{unit,*1e12=pC}';'S2:q{unit,*1e12=pC}';'S3:q{unit,*1e12=pC}'};
dIn=struct('def',[1;0;0;0],'DOF',i1+.21,'lab',{{'Act:v'}}, ...%Act
    'Elt',feutil('selelt proid 104',model));

% Define Node names 
model=sdth.urn('nmap.Node.set',model, ...
    {'Tip',nd;'Act',i1(1);'S1',i1(2);'S2',i1(3);'S3',i1(4)});
% Define sensors 
model=fe_case(model,'SensDof','Sensors',scell);
% Define input (and 0 tension to measure charge)
model=fe_case(model,'DofSet','In',dIn); % Elt defined for display

% Ref solution with state-space 10 modes+static
f=linspace(1,100,400); % in Hz
s0=fe2ss('free 5 10 0 -dterm',model); %
C0=qbode(s0,f(:)*2*pi,'struct-lab');C0.name='Free+Static';
%% Step 2 - CB reduction
model=fe_case(model,'DOFSet','In',[nd+.03; i1+.21]); % 

% Build CB matrices
model = stack_set(model,'info','EigOpt',[5 10 0]);
SE1= fe_reduc('CraigBampton -SE -matdes 2 1 3 4 -bset',model); % 
% To keep only the matrices
SE1 = rmfield(SE1,{'il','pl','Stack','TR'}) ; 
SE1.Node=feutil('getnode',SE1,unique(fix(SE1.DOF)));SE1.Elt=[];

% Define node names
SE1=sdth.urn('nmap.Node.set',SE1, ...
    {'Tip',nd;'Act',i1(1);'S1',i1(2);'S2',i1(3);'S3',i1(4)});

% Define input/output
SE1=fe_case(SE1,'SensDof','Sensors',scell);
SE1=fe_case(SE1,'DofSet','In',rmfield(dIn,'Elt')); % Elt defined for display
SE1=stack_set(SE1,'info','DefaultZeta',1e-2);

% Build state-space model based on reduced CB matrices
s1=fe2ss('free 5 10 0 -dterm',SE1);
%% Step 3: Now with Raze's transform, using -noMCI option in fe_reduc
SE2= fe_reduc('CraigBampton -SE -matdes 2 1 3 4 -noMCI',model); % Removes M-coupling 
% To keep only the matrices
SE2 = rmfield(SE2,{'il','pl','Stack','TR'}) ; 
SE2.Node=feutil('getnode',SE2,unique(fix(SE2.DOF)));SE2.Elt=[];

%Define node names
SE2=sdth.urn('nmap.Node.set',SE2, ...
    {'Tip',nd;'Act',i1(1);'S1',i1(2);'S2',i1(3);'S3',i1(4)});
% Define input/output
SE2=fe_case(SE2,'SensDof','Sensors',scell);
SE2=fe_case(SE2,'DofSet','In',rmfield(dIn,'Elt')); % Elt defined for display
SE2=stack_set(SE2,'info','DefaultZeta',1e-2);

% Build state-space model based on reduced CB matrices
s2=fe2ss('free 5 10 0 -dterm',SE2);
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopzpatchnumide')
%% #TutoPzPatchNumIDE : Piezo patch with IDE -2

% see sdtweb pz_composite#tutopzpatchnumide
%
%% BeginSource sdtweb('_example','pz_composite.tex#tutopzpatchnumide')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzPatchNumIDE')
d_piezo('DefineStyles');

%% Step 1 - Build mesh and compute static response
% Meshing script can be viewed with sdtweb d_piezo('MeshIDEPatch')
% Build mesh, electrodes and actuation
model=d_piezo(['MeshIDEPatch nx=10 ny=5 nz=14 lx=400e-6' ...
 'ly=300e-6 p0=700e-6 e0=50e-6']);
% Transform in mm
model.Node(:,5:7)=model.Node(:,5:7)*1000; % From m to mm
model.unit='mm';
% Convert material properties
model.pl = fe_mat('convert SI mm',model.pl);
% Compute response due to V and visualize
model=fe_case(model,'pcond','Piezo','d_piezo(''Pcond'')');
% low freq response to avoid rigid body modes
model=stack_set(model,'info','Freq',10);
def=fe_simul('dfrf',model);
% Plot deformed shape
cf=feplot(model,def); fecom('view3'); fecom('viewy-90'); fecom('viewz+90')
fecom('undef line'); fecom('triax') ; iimouse('zoom reset')
cf.mdl.name='patch_IDE_deformed';d_piezo('SetStyle',cf); feplot(cf);
%% Step 2 - visualize electric field
cf.sel(1)={'groupall','colorface none -facealpha0 -edgealpha.1'};
p_piezo('viewElec EltSel "matid1" DefLen 50e-3 reset',cf);
fecom('scd 1e-10')
p_piezo('electrodeview -fw',cf); % to see the electrodes on the mesh
iimouse('zoom reset')
cf.mdl.name='patch_IDE_EField';d_piezo('SetStyle',cf); feplot(cf);
%% Step 3 - Compare effective values of constitutive law
% Decompose constitutive law
CC=p_piezo('viewdd -struct',cf); %
% Compute mean value of fields and deduce equivalent d_ij
% Uniform field is assumed for analytical values
a=p_piezo('viewstrain -curve -mean',cf); % mean value of S1-6 and E1-3
fprintf('Relation between mean strain on free structure and d_3i\n');
E3=a.Y(9,1); disp({'E3 mean' a.Y(9,1) -1/700e-3 'E3 analytic'})
disp([{'Sx/E3';'Sy/E3';'Sz/E3'} num2cell([a.Y(1:3,1)/E3 CC.d(3,1:3)']) ...
{'d_31(mm/muV)';'d_32(mm/muV)';'d_33(mm/muV)'}])
%% Step 4 - Charge visualisation and capacitance
p_piezo('electrodeTotal',cf);
% charge density on the electrodes
feplot(model,def);
cut=p_piezo('electrodeviewcharge',cf,struct('EltSel','matid 1'));
fecom('view3'); fecom('viewy-90'); fecom('viewz+90'); iimouse('zoom reset');
%iimouse('trans2d 0 0 0 1.6 1.6 1.6')
set(cf.ga,'climmode','auto');
cf.mdl.name='patch_IDE_charge_elec';d_piezo('SetStyle',cf); feplot(cf);
% - Theoretical capacitance for uniform field
Ct=model.pl(1,22)*400e-3*300e-3/700e-3;
% total charge on the electrodes = capacitance (1 muV actuation)
C=p_piezo('electrodeTotal',cf);
% Differences are due to non-uniform field, this is to be expected
disp({'C_{IDE}' cell2mat(C(2,2)) Ct 'C analytic'})
%% Step 5 - Stress and strain  visualisation
% Stress field using fe_stress
c1=fe_stress('stressAtInteg -gstate',model,def);
cf.sel='reset';cf.def=fe_stress('expand',model,c1);
cf.def.lab={'T11';'T22';'T33';'T23';'T13';'T12';'D1';'D2';'D3'}; %
fecom('colordata 99 -edgealpha.1');
fecom('colorbar',d_imw('get','CbTR','String','Stress/Voltage [kPa/muV]'));
iimouse('trans2d 0 0 0 1.6 1.6 1.6')
%- Strain  visualisation
% Replace with 'PiezoStrain' material
mo2=model; mo2.pl=m_piezo('dbval 1 PiezoStrain')
% Now represent strain fields using fe_stress
c1=fe_stress('stressAtInteg -gstate',mo2,def);
cf.sel='reset';cf.def=fe_stress('expand',mo2,c1);
cf.def.lab={'S11';'S22';'S33';'S23';'S13';'S12';'E1';'E2';'E3'};
fecom('colordata 99 -edgealpha.1');
fecom('colorbar',d_imw('get','CbTR','String','Strain/Voltage [1/muV]'));
iimouse('trans2d 0 0 0 1.6 1.6 1.6')

%% EndSource EndTuto

elseif comstr(Cam,'tutopzMFCP2homo')
%% #TutoPzMFCP2homo : P2-type MFC homogenization -2

% see sdtweb pz_composite#tutopzMFCP2homo
%
%% BeginSource sdtweb('_example','pz_composite.tex#tutopzMFCP2homo')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzMFCP2Homo')
d_piezo('DefineStyles');

%% Step 1 - Meshing of RVE
% Meshing script can be viewed with sdtweb d_piezo('MeshHomoMFCP2')
Range=fe_range('grid',struct('rho',[0.001 linspace(0.1,0.9,9) .999], ...
  'lx',.300,'ly',.300,'lz',.180,'dd',0.04));
%% Step 2 - Loop on volume fraction and compute homogenize properties
for jPar=1:size(Range.val,1)

RO=fe_range('valCell',Range,jPar,struct('Table',2));% Current experiment

% Create mesh
model= ...
d_piezo(sprintf(['meshhomomfcp2 rho=%0.5g lx=%0.5g ly=%0.5g' ...
  'lz=%0.5g dd=%0.5g'],[RO.rho RO.lx RO.ly RO.lz RO.dd]));


% Define the six local problems
RB=struct('CellDir',[max(model.dx) max(model.dy) max(model.dz)],'Load', ...
       {{'e11','e22','e12','e23','e13','vIn'}});
 % Periodicity on u,v,w on x and y face
 % periodicity on u,v only on the z face
  RB.DirDofInd={[1:3 0],[1:3 0],[1 2 0 0]};
 % Voltage DOFs are always eliminated from periodic conditions

 % Compute the deformation for the six local problems
 def=fe_homo('RveSimpleLoad',model,RB);
 % Represent the deformation of the RVE
 cf=comgui('guifeplot-reset',2);cf=feplot(model,def); fecom('colordatamat')
% Compute stresses, strains and electric field
a1=p_piezo('viewstrain -curve -mean -EltSel MatId1 reset',cf); % Strain epoxy
a2=p_piezo('viewstrain -curve -mean -EltSel MatId2 reset',cf); % Strain piezo
b1=p_piezo('viewstress -curve -mean- EltSel MatId1 reset',cf); % Stress epoxy
b2=p_piezo('viewstress -curve -mean- EltSel MatId2 reset',cf); % Stress piezo

% Compute charge on electrodes
 mo1=cf.mdl.GetData;
  ind1=fe_case(mo1,'getdata','Top Actuator');ind1=fix(ind1.InputDOF);
  mo1=p_piezo('electrodesensq TopQ2',mo1,struct('MatId',2,'InNode',ind1));
  mo1=p_piezo('electrodesensq TopQ1',mo1,struct('MatId',1,'InNode',ind1));
  c1=fe_case('sensobserve',mo1,'TopQ1',cf.def); q1=c1.Y;
  c2=fe_case('sensobserve',mo1,'TopQ2',cf.def); q2=c2.Y;


% Compute average values:
a0=a1.Y(1:6,:)*(1-RO.rho)+a2.Y(1:6,:)*RO.rho;
b0=b1.Y(1:6,:)*(1-RO.rho)+b2.Y(1:6,:)*RO.rho;
q0=q1+q2; % Total charge is the sum of charges on both parts of electrode


% Compute C matrix
C11=b0(1,1)/a0(1,1); C12=b0(1,2)/a0(2,2); C22=b0(2,2)/a0(2,2);
C44=b0(4,4)/a0(4,4); C55=b0(5,5)/a0(5,5); C66=b0(6,3)/a0(6,3);
sE=inv([C11 C12; C12 C22]);

% Extract mechanical engineering constants
E1(jPar)=1/sE(1,1); E2(jPar)=1/sE(2,2); nu12(jPar)=-sE(1,2)*E1(jPar);
nu21(jPar)=-sE(1,2)*E2(jPar);G12(jPar)=C66; G23(jPar)=C44; G13(jPar)=C55;

% Extract piezoelectric properties
e31(jPar)=b0(1,6)*RO.lz; e32(jPar)=b0(2,6)*RO.lz;
d=[e31(jPar) e32(jPar)]*sE; d31(jPar)=d(1); d32(jPar)=d(2);

% Extract dielectric properties
eps33(jPar)=-q0(6)*RO.lz/(RO.lx*RO.ly);
eps33t(jPar)=eps33(jPar)+ [d31(jPar) d32(jPar)]*[e31(jPar); e32(jPar)];

end % Loop on rho0 values
%% Step 3 - Homogeneous properties as a function of volume fraction
rho0=Range.val(:,strcmpi(Range.lab,'rho'));

out=struct('X',{{rho0,{'E_T','E_L','nu_{LT}','G_{LT}','G_{Tz}','G_{Lz}', ...
    'e_{31}','e_{32}','d_{31}','d_{32}','epsilon_{33}^T'}'}},'Xlab',...
    {{'\rho','Component'}},'Y',[E1' E2' nu21' G12' G13' G23' e32' e31' ...
    d32' d31' eps33t'/8.854e-12]);

ci=iiplot; iicom('CurveReset');
iicom(ci,'CurveInit','P2-MFC homogenization',out);
d_piezo('setstyle',ci);

%% EndSource EndTuto

elseif comstr(Cam,'tutopzMFCP1homo')
%% #TutoPzMFCP1homo : P1-type MFC homogenization -2

% see sdtweb pz_composite#tutopzMFCP1homo
%
%% BeginSource sdtweb('_example','pz_composite.tex#tutopzMFCP1homo')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzMFCP1homo')
d_piezo('DefineStyles');

%% Step 1 - Meshing or RVE and definition of volume fractions
% Meshing script can be viewed with sdtweb d_piezo('MeshHomoMFCP1')
Range=fe_range('grid',struct('rho',[0.001 linspace(0.1,0.9,9) .999], ...
  'lx',.18,'ly',.18,'lz',1.080,'e',0.09,'dd',0.04));
%% Step 2 - Loop on volume fractions and computation of homogenized properties
for jPar=1:size(Range.val,1)

RO=fe_range('valCell',Range,jPar,struct('Table',2));% Current experiment

% Create mesh
model=...
d_piezo(sprintf(['meshhomomfcp1 rho=%0.5g lx=%0.5g ly=%0.5g lz=%0.5g' ...
   ' e=%0.5g dd=%0.5g'],[RO.rho RO.lx RO.ly RO.lz RO.e RO.dd]));

%%
RB=struct('CellDir',[max(model.dx) max(model.dy) max(model.dz)],'Load', ...
           {{'e33','e11','e23','e12','e13','vIn'}});
 % It seems fe_homo reorders the strains
 % Periodicity on u,v,w on x and z face
 % periodicity on u,w only on the y face
 % Voltage DOFs are always eliminated from periodic conditions
 RB.DirDofInd={[1:3 0],[1 0 3 0],[1:3 0]};
 def=fe_homo('RveSimpleLoad',model,RB);

cf=comgui('guifeplot-reset',2);
cf=feplot(model,def); fecom('colordatamat'); fecom('triax')
% Electric field for Vin
p_piezo('electrodeview -fw',cf); % to see the electrodes on the mesh
cf.sel(1)={'groupall','colorface none -facealpha0 -edgealpha.1'};
p_piezo('viewElec EltSel "matid1:2" DefLen 0.07 reset',cf);
fecom('scd 1e-10')
% Compute stresses, strains and electric field
a1=p_piezo('viewstrain -curve -mean -EltSel MatId1 reset',cf); % Strain S epoxy
a2=p_piezo('viewstrain -curve -mean -EltSel MatId2 reset',cf); % Strain S piezo
b1=p_piezo('viewstress -curve -mean- EltSel MatId1 reset',cf); % Stress T
b2=p_piezo('viewstress -curve -mean- EltSel MatId2 reset',cf); % Stress T

% Compute charge
 mo1=cf.mdl.GetData;
  ind1=fe_case(mo1,'getdata','Top Actuator');ind1=fix(ind1.InputDOF);
  mo1=p_piezo('electrodesensq TopQ2',mo1,struct('MatId',2,'InNode',ind1));
  mo1=p_piezo('electrodesensq TopQ1',mo1,struct('MatId',1,'InNode',ind1));
  c1=fe_case('sensobserve',mo1,'TopQ1',cf.def); q1=c1.Y;
  c2=fe_case('sensobserve',mo1,'TopQ2',cf.def); q2=c2.Y;

% Compute average values:
a0=a1.Y(1:6,:)*(1-RO.rho)+a2.Y(1:6,:)*RO.rho;
b0=b1.Y(1:6,:)*(1-RO.rho)+b2.Y(1:6,:)*RO.rho;
q0=q1+q2; %total charge is the sum of charges

% Compute C matrix
C11=b0(3,2)/a0(3,2); C12=b0(3,1)/a0(1,1); C22=b0(1,1)/a0(1,1);
C44=b0(5,5)/a0(5,5); C55=b0(4,4)/a0(4,4); C66=b0(6,3)/a0(6,3);
sE=inv([C11 C12; C12 C22]);
E1(jPar)=1/sE(1,1); E2(jPar)=1/sE(2,2); nu12(jPar)=-sE(1,2)*E1(jPar);
nu21(jPar)=-sE(1,2)*E2(jPar);
e33(jPar)=b0(3,6)*RO.lz;
e31(jPar)=b0(1,6)*RO.lz;
eps33(jPar)=-q0(6)*RO.lz/(RO.ly*RO.lx);
d=[e33(jPar) e31(jPar)]*sE; d33(jPar)=d(1); d31(jPar)=d(2);
eps33t(jPar)=eps33(jPar)+ [d33(jPar) d31(jPar)]*[e33(jPar); e31(jPar)];
G12(jPar)=C44;
G23(jPar)=C66;
G13(jPar)=C55;

end % loop on Rho values
%% Step 3 - Plot homogeneous properties as a function of volume fraction
rho0=Range.val(:,strcmpi(Range.lab,'rho'));

out=struct('X',{{rho0,{'E_L','E_T','nu_{LT}','G_{LT}','G_{Lz}','G_{Tz}', ...
    'e_{31}','e_{33}','d_{31}','d_{33}','epsilon_{33}^T'}'}},'Xlab',...
    {{'\rho','Component'}},'Y',[E1' E2' nu12' G12' G13' G23' e31' e33' ...
    d31' d33' eps33t'/8.854e-12]);
ci=iiplot;
iicom('CurveReset');
iicom(ci,'CurveInit','P1-MFC homogenization',out);
d_piezo('setstyle',ci)

%% EndSource EndTuto

elseif comstr(Cam,'tutopzMFCPlate')
%% #TutoPzMFCP1homo : Cantilever Plate with MFC transducers -2

% see sdtweb pz_composite#tutopzMFCPlate
%
%% BeginSource sdtweb('_example','pz_composite.tex#tutopzMFCPlate')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzMFCPlate')
d_piezo('DefineStyles');

%% Step 1 - Build mesh and visualize
% Meshing script,open with sdtweb d_piezo('MeshMFCplate')
model=d_piezo('MeshMFCplate -cantilever');  % creates the model
cf=feplot(model); fecom('colordatagroup-EdgeAlpha.1');
%- Define actuators and sensors
r1=p_piezo('electrodedof',model);
data.def=[1 -1;1 1]'; % Define combinations for actuators
data.lab={'V-bend';'V-Tract'};
data.DOF=vertcat(r1{:,2})+.21;
model=fe_case(model,'DofSet','V_{In}',data);
% Force use of DofSet (compatibility mode) by setting .ver=1
model.Stack{strcmpi(model.Stack(:,2),'Electrodes'),3}.ver=1;

% Add tip displacement sensors in z
nd1=feutil('find node x==463 & y==50',model);
nd2=feutil('find node x==463 & y==0',model);
model=fe_case(model,'SensDof','Tipt-z',nd1+.03); % Z-disp
model=fe_case(model,'SensDof','Tip-x',nd1+.01); % X-disp
model=fe_case(model,'SensDof','Tipb-z',nd2+.03); % Z-disp
%- Compute static response
d0=fe_simul('dfrf',stack_set(model,'info','Freq',0)); %
cf=feplot(model,d0); sens=fe_case(model,'sens');
C1=fe_case('SensObserve -dim 2 3 1',sens,d0);
fecom(';view3;scd 20;colordataEvalA;undefline')
%% Step 2 - Rotate fibers
model.il(2,[20 44])=[45 -45];
d1=fe_simul('dfrf',stack_set(model,'info','Freq',0)); % static response
cf.def=d1; fecom('scd 10');  C2=fe_case('SensObserve - dim 2 3 1',sens,d1);

%% EndSource EndTuto

elseif comstr(Cam,'tutopztriangle')
%% #TutoPzTriangle : Triangular point load actuator -2

% see sdtweb pz_composite#tutopztriangle
%
%% BeginSource sdtweb('_example','pz_composite.tex#tutopztriangle')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzTriangle')
d_piezo('DefineStyles');

%% Step 1 - Build Mesh using gmsh and visualize
% Meshing script can be viewed with sdtweb d_piezo('MeshTrianglePlate')
% --- requires gmsh
model=d_piezo('MeshTrianglePlate');
cf=feplot(model); fecom('colordatapro'); fecom('view2')
%% Step 2 - Define actuators and sensors - static response
model=fe_case(model,'SensDof','Tip',7.03); % Displ sensor
model=fe_case(model,'DofSet','V-Act', ...
  struct('def',[-1; 1],'DOF',[100001; 100002]+.21));
% - Compute static response to voltage actuation
d0=fe_simul('dfrf',stack_set(model,'info','Freq',0));
cf.def=d0; fecom('colordataz -alpha .8 -edgealpha .1')
fecom('scd -.03'); fecom('view3');
%% Step 3 - Compute dynamic response with state-space model
[sys,TR]=fe2ss('free 5 20 0 -dterm',model);
C1=qbode(sys,linspace(0,500,500)'*2*pi,'struct-lab'); C1.name='SS20modes';

%% Point load actuation
model=fe_case(model,'Remove','V-Act'); % remove piezo actuator
model=fe_case(model,'FixDof','Piezos',[100001;100002]); 
%SC piezo electrodes

% Determine scaling factor, check b/l ratio and build point force
CC=p_piezo('viewdd -struct',model);
a=100; b=33.58;
zm=0.650e-3; V=1; e31=CC.e(1); A=-(e31*zm*V*b)/a; A=A*2; 

% Two triangles
bl= 2*sqrt(-CC.e(2)/CC.e(1));

data=struct('DOF',[7.03],'def',A); data.lab=fe_curve('datatype',13);
model=fe_case(model,'DofLoad','PointLoad',data);

% Static response to point load
d1=fe_simul('dfrf',stack_set(model,'info','Freq',0));
ind=fe_c(d1.DOF,7.03,'ind'); d1p=d1.def(ind);

% Dynamic response (reduced modal model)
[sys,TR]=fe2ss('free 5 20 0 -dterm',model);
C2=qbode(sys,linspace(0,500,1000)'*2*pi,'struct-lab'); C2.name='SS20modes';
% Compare frequency responses
ci=iiplot;
C1.X{2}={'dtip'}; C2.X{2}={'dtip'}; 
iicom(ci,'curveinit',{'curve',C1.name,C1;'curve',C2.name,C2});
iicom('submagpha'); d_piezo('setstyle',ci)
% End of script

%% EndSource EndTuto

elseif comstr(Cam,'tutopzaccel')
%% #TutoPzAccel : Piezo Accelerometer: sensitivity -2

% see sdtweb pz_applications#tutopzaccel
%
%% BeginSource sdtweb('_example','pz_applications.tex#tutopzaccel')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzAccel')
d_piezo('DefineStyles');

%% Step 1 - Build Mesh and visualize
% Meshing script can be viewed with sdtweb d_piezo('MeshBaseAccel')
model=d_piezo('MeshBaseAccel');
cf=feplot(model); fecom('colordatagroup');
set(gca,'cameraposition',[-0.0604   -0.0787    0.0139])
cf.mdl.name='accelero_mesh'; d_piezo('SetStyle',cf); feplot(cf);
%% Step 2 - Define sensors and actuators
% -MatID 2 requests a charge resultant sensor
% -vout requests a voltage sensor
model=p_piezo('ElectrodeMPC Top sensor -matid 2 -vout',model,'z==0.004');
% -ground generates a v=0 FixDof case entry
model=p_piezo('ElectrodeMPC Bottom sensor -ground',model,'z==0.003');
% Add a displacement sensor for the basis
model=fe_case(model,'SensDof','Base-displ',1.03);
%% Step 3 - Response with imposed displacement
% Remove the charge sensor (not needed)
model=fe_case(model,'remove','Q-Top sensor');

% Link dofs of base and impose unit vertical displacement
n1=feutil('getnode z==0',model);
rb=feutilb('geomrb',n1,[0 0 0],fe_c(feutil('getdof',model),n1(:,1),'dof'));
rb=fe_def('subdef',rb,3); % Keep vertical displacement
model=fe_case(model,'DofSet','Base',rb);

% Other parameters
f=linspace(1e3,2e5,200)';
model=stack_set(model,'info','Freq',f); % freq. for computation

 % Reduced ss-model
[sys,TR1]=fe2ss('free 5 10 0 -dterm',model,5e-3); %

C1=qbode(sys,f(:)*2*pi,'struct-lab');C1.name='SS-voltage';
C1.X{2}={'Sensor output(V)';'Base Acc(m/s^2)';'Sensitivity (V/m/s^2)'}; %outputs

% C1 compute accel and sensitivity
 C1.Y(:,2)=C1.Y(:,2).*(-(C1.X{1}*2*pi).^2); % Base acc
 C1.Y(:,3)=C1.Y(:,1)./C1.Y(:,2);% Sensitivity
 
 ci=iiplot; iicom(ci,'curveinit',{'curve',C1.name,C1});  iicom('ch3');
 d_piezo('setstyle',ci)

 %% Step 4 - Compare charge and voltage mode for sensing
 
model=d_piezo('MeshBaseAccel');
model=fe_case(model,'remove','V-Top sensor');

% Short-circuit electrodes of accelerometer
model=fe_case(model,'FixDof','V=0 on Top Sensor', ...
    p_piezo('electrodedof Top sensor',model));

% Other parameters
  model=stack_set(model,'info','Freq',f);

% Link dofs of base and impose unit vertical displacement
n1=feutil('getnode z==0',model);
rb=feutilb('geomrb',n1,[0 0 0],fe_c(feutil('getdof',model),n1(:,1),'dof'));
rb=fe_def('subdef',rb,3); % Keep vertical displacement
model=fe_case(model,'DofSet','Base',rb);

% Reduced model
[sys2,TR1]=fe2ss('free 5 10 0 -dterm',model,5e-3); %
C2=qbode(sys2,f(:)*2*pi,'struct');C2.name='SS-charge';

C2.X{2}={'Sensor output(V)';'Base Acc(m/s^2)';'Sensitivity (normalized)'}; %outputs
C1.X{2}={'Sensor output(q)';'Base Acc(m/s^2)';'Sensitivity (normalized)'};

% C2 compute accel and sensitivity
 C2.Y(:,2)=C2.Y(:,2).*(-(C2.X{1}*2*pi).^2); % Base acc
 C2.Y(:,3)=C2.Y(:,1)./C2.Y(:,2);% Sensitivity

 % Normalize sensitivity to first freq
 C1.Y(:,3)=C1.Y(:,3)/abs(C1.Y(1,3));
 C2.Y(:,3)=C2.Y(:,3)/abs(C2.Y(1,3))

 ci=iiplot; iicom(ci,'curveinit',{'curve',C1.name,C1;'curve',C2.name,C2}); ...
 iicom('ch3');  d_piezo('setstyle',ci)

%% EndSource EndTuto

elseif comstr(Cam,'tutopzaccshakersens')
%% #TutoPzAccShakerSens : Piezo Accelerometer: sensitivity wt Pz shaker -2

% see sdtweb pz_applications#tutopzaccshakersens
%
%% BeginSource sdtweb('_example','pz_applications.tex#tutopzaccshakersens')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzAccShakerSens')
d_piezo('DefineStyles');

%% Step 1 - Build mesh and visualize
% Meshing script,open with sdtweb d_piezo('MeshPiezoShaker')
model=d_piezo('MeshPiezoShaker');
cf=feplot(model); fecom('colordatapro');
set(gca,'cameraposition',[-0.0604   -0.0787    0.0139])
iimouse('resetview'); 
cf.mdl.name='Acc_Shaker_Mesh'; % Model name for title
d_piezo('SetStyle',cf); feplot(cf);
%% Step 2 - Define actuators and sensors
  % -input "In" says it will be used as a voltage actuator
model=p_piezo('ElectrodeMPC Top Actuator -input "Vin-Shaker"',model,'z==-0.01');
  % -ground generates a v=0 FixDof case entry
model=p_piezo('ElectrodeMPC Bottom Actuator -ground',model,'z==-0.012');
% Voltage sensor will be used - remove charge sensor
model=fe_case(model,'remove','Q-Top sensor');

% Replace by an acceleration sensor for the basis
model=fe_case(model,'remove','Base-displ');
model = fe_case(model,'SensDOF','Accel',{'1:az'});
%% Step 3 - Compute response, voltage input on shaker
f=linspace(1e3,2e5,200)';
model=stack_set(model,'info','Freq',f);

% Reduced ss-model
model=fe_case(model,'pcond','Piezo','d_piezo(''Pcond'')');
[sys,TR1]=fe2ss('free 5 45 0 -dterm ',model,1e-3); %

C1=qbode(sys,f(:)*2*pi,'struct-lab');C1.name='SS-voltage';
C1.X{2}={'Sensor output(V)';'Base Acc(m/s^2)';'Sensitivity (V/m/s^2)'}; %outputs

% Compute sensitivity
C1.Y(:,3)=C1.Y(:,1)./C1.Y(:,2);% Sensitivity
 
ci=iiplot; iicom(ci,'curveinit',{'curve',C1.name,C1});  iicom('ch3');
d_piezo('setstyle',ci)

%% EndSource EndTuto

elseif comstr(Cam,'tutopzshunt')
%% #TutoPzShunt : Cantilever plate with patches : RL shunt damping -2

% see sdtweb pz_applications#tutopzshunt
%
%% BeginSource sdtweb('_example','pz_applications.tex#tutopzshunt')

% Init working directory for figure generation
d_piezo('SetPlotwd');
% See full example as MATLAB code in d_piezo('ScriptTutoPzShunt')
d_piezo('DefineStyles');

%% Step 1 - Build mesh and visualize
% Meshing script can be viewed with sdtweb d_piezo('MeshShunt')
model=d_piezo('meshshunt');
model=stack_set(model,'info','DefaultZeta',1e-4);
feplot(model); cf=fecom; fecom('colordatapro')
cf.mdl.name='Shunt_mesh'; d_piezo('setstyle',cf); feplot(cf);
%% Step 2 - Define actuators and sensors
% Actuators
data.def=[1 -1 0 0; 0 0 1 -1]'; % Define combinations for actuators
data.DOF=p_piezo('electrodedof.*',model); edof=data.DOF;
model=fe_case(model,'DofSet','Vin',data);
% Sensors
r1=struct('cta',[1 -1],'DOF',edof(1:2),'name','Qs');
model=p_piezo('ElectrodeSensQ',model,r1);
nd=feutil('find node x==350 & y==25',model);
model=fe_case(model,'SensDof','Tip',nd+.03);
sens=fe_case(model,'sens');
 %% Step 3: compute response
 w=linspace(0,1e3,1e4)'*2*pi;
 [sys,TR]=fe2ss('free 5 30 0 -dterm',model);
 C1=qbode(sys,w,'struct'); C1.name='no shunt';
 C1.X{2}={'Qs';'tip-displ'}; %outputs 
 C1.X{3}={'Vin-u1';'Vin-u2'}; %inputs 
 ci=iiplot;
 iicom('CurveReset');iicom('curveinit',C1); iicom(ci,'xlim[0 30]')
 d_piezo('setstyle',ci)
%% Step 4 - Determine parameters for shunt tuning
 % Extract w1 and W1 and compute alpha_1
C=C1.Y(:,1);
% Find poles and zeros of impedance (1/jwC)
if ~exist('findpeaks','file'); warning('Skipping step, signal toolbox');
return;end
[pksPoles,locsPoles]=findpeaks(abs(1./C)); Wi=w(locsPoles);
[pksZeros,locsZeros]=findpeaks(abs(C)); wi=w(locsZeros);
% concentrate on mode of interest (mode 1)
 W1=w(locsPoles(1)); w1=w(locsZeros(1));
% Compute alpha for mode of interest
a1=sqrt((W1^2-w1^2)/W1^2);
% Compute Cs2 for mode of interest
ind1=1; i2=locsZeros(1); i3=locsZeros(2);
dw2=w(i3)-w(i2); wCs2=w(i2)+dw2/2;
[y,i]=min(abs(w-wCs2)); Cs2=abs(C(i));
%% Determine shunt parameters (R and L) and apply it to damp 1st mode
% Tuning using Yamada's rule
d=1; r=sqrt((3*a1^2)/(2-a1^2));
L_Yam=1/d^2/Cs2/W1^2; R_Yam=r/Cs2/W1;
%% Step 5 - Compute dynamic response with optimal shunt
w=linspace(0,40,1e3)*2*pi;
C1=qbode(sys,w,'struct'); C1.name='no shunt';
C1.X{2}={'Qs';'tip-displ'}; %outputs 
C1.X{3}={'Vin-u1';'Vin-u2'}; %inputs 

% Implement shunt using feeback - requires control toolbox - compute FRF
if ~exist('feedback','file'); warning('Skipping step, control toolbox');
  return;
end
A=tf([L_Yam R_Yam 0],1); % RL shunt in tf form
sys2=feedback(sys,A,1,1,1);
% qbode does not work with feedback so use freqresp from control toolbox
 C=freqresp(sys2,w); a=C(:); C2=C1; C2.Y=reshape(a,4,1000)';
 C2.name='RL shunt'; C1.X=C2.X;

% Plot and compare curves
iicom(ci,'curveinit',{'curve',C1.name,C1;'curve',C2.name,C2});
iicom(ci,'ch 4'); d_piezo('setstyle',ci);
comgui('imwrite',ci)

%% EndSource EndTuto

elseif comstr(Cam,'fullconstrain');
%% #ScriptFullConstrain : no mechanical displacement and zero potential -2 
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




else; error('Script%s unknown',CAM);
    
end

elseif comstr(Cam,'mesh');[CAM,Cam]=comstr(CAM,5);
%% #Mesh -------------------------------------------------------------------

if comstr(Cam,'patch');[CAM,Cam]=comstr(CAM,6);
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

elseif comstr(Cam,'pic181disk');[CAM,Cam]=comstr(CAM,6);
%% #MeshPIC181disk : mesh a disk made of PIC181 bulk piezo

% Parameter handling
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
    ' th(2e-3#%g#"disk thickness")'...
    ' r(8e-3#%g#"disk radius")'...
    ' ner(10#%g#"nb elts along radius")'...
    ' nez(4#%g#"nb elts along thickness")'...
    ' nrev(16#%g#"nb elts along circunf")'...
   ],{RO,CAM});


%% Mesh
%% Make mesh
nd=[0 0 0; RO.r 0 0; RO.r 0 RO.th; 0 0 RO.th];
model=feutil('object quad 1 1',nd,RO.ner,RO.nez); % Piezo
model=feutil(['Rev ' num2str(RO.nrev) ' o 0 0 0 360 0 0 1'],model);
model.unit='SI'; model.name='PIC 181 disk';

%% Define material properties
pl=m_piezo('dbval 1 -elas 2 PIC_181');

%Damping
pl(2,7)=0.02; % 1% damping in piezo
model.pl=pl; model=p_solid('default;',model);

%% Define electrodes, actuator and sensor
  % -input "In" says it will be used as a voltage actuator
model=p_piezo('ElectrodeMPC Top Actuator -input "Vin"',model,['z==' num2str(RO.th)]);
  % -ground generates a v=0 FixDof case entry
model=p_piezo('ElectrodeMPC Bottom Actuator -ground',model,'z==0');
  % add a charge sensor on the top electrode
model=p_piezo(['ElectrodeSensQ '  ...
    num2str(floor(p_piezo('electrodedof Top Actuator',model))) ' Q'],model);

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
model=[]; 
if isfield(RO,'Elt');model=RO;RO=varargin{carg};carg=carg+1;end
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   ' unit("SI"#%s#" unit system")'...
   ' Quad(#3#"quadratic elements if present")'...
   ' MatId(100#%g#"Starting matid for piezos")'...
   ' Orient([]#%g#"possibly give orientation")'...
   ],{RO,CAM(7:end)});
RunOpt.NeedOrient=~isempty(RO.Orient);
RO.PatchId=[];
if ~isfield(RO,'list');error('.list expected for meshing');
else
%% #MeshPlateList : Most general call based on list of features -3
 if ~isfield(RO,'MatId');RO.MatId=100;end % Offset for new matid
 if ~isfield(RO,'ProId');RO.ProId=100;end % Offset for new proid
 RO.Rect=[]; RO.Circ=[];  RO.Cam0=CAM; RO.Electrode=[];
 for j1=1:size(RO.list,1)
   %% Deal with lam(inate) column
   st=RO.list{j1,3}; if ischar(st);RB=struct;else; RB=st;st='';end 
   RB.CAM=st;
   if ischar(RO.list{j1,2})&&strncmpi(RO.list{j1,2},'selelt',6)
     [RO.EltOri,model.Elt]=feutil(['RemoveElt',RO.list{j1,2}(7:end)],model);
     RO.list{j1,2}=model;
   end
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
     if ~isfield(model,'unit');error('Expecting model.unit to be defined');
     else; RB.unit=model.unit;
     end
     
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
     if ~isfield(RB,'Ecoef')
     elseif length(RB.Ecoef)>2
       i1=find(model.pl(:,1)==RB.Ecoef(1)); 
       if isempty(i1);error('Material %i not found',RB.ECoef(1));end
       ind=find(RB.Ecoef~=1);
       model.pl(i1,ind)=model.pl(i1,ind).*RB.Ecoef(ind);
     else
      model=d_shm(sprintf('cbMatDef %i %.15g', ...
          RO.list{2,2}.pl(1),RB.Ecoef),model);
     end
     model=fe_mat('convert',model);
     if isfield(model,'pl')&&~isempty(model.pl)&&any(model.pl(:,1)==RO.MatId)
       RO.MatId=max(model.pl(:,1))+1;
     end
     if isfield(model,'il')&&~isempty(model.il)&&any(model.il(:,1)==RO.ProId)
       RO.ProId=max(model.il(:,1))+1;
     end
   elseif strncmpi(RO.list{j1,2},'mdl',3)
     %% Refine coarse model -3
     model=RB.mdl;
     mdl=fe_fmesh(sprintf('qmesh %g -noTest',RB.lc*2),model);
     model.Node=mdl.Node;model.Elt=mdl.Elt;
     
   elseif strncmpi(RO.list{j1,2},'baseid',5)
     %% Combination_of_base+patch information -3

     % Mesh the patch, RD gives patch options
     RB.MatId=RO.MatId;RB.ProId=RO.ProId;RB.name=RO.list{j1,2};RC=RB;
     [model,RB]=m_piezo('patchaddpro',model,RB);RB.CAM=RC.CAM;
     % Added electrodes (should relocate)
     RO.MatId=RB.MatId; RO.ProId=RB.ProId;
     if size(RO.list,2)>2&&isstruct(RO.list{j1,3})
      RB0=RB;% should be removed
      RB=sdth.sfield('MergeI',RO.list{j1,3},RB, ...
          {'MatId','ProId','Electrode','CAM','il','shape','lx','ly', ...
          'rc','idc','MakeSandwich'});
     end
     if isfield(RB,'Electrode')&&~isempty(RB.Electrode)
      RO.Electrode(end+(1:size(RB.Electrode,1)),1:2)=RB.Electrode;
     end
     
     switch lower(RB.shape)
     case {'circ','annulus'}
      [model,RO,RB]=meshCirc(model,RO,RB);
     case 'rect'
      [model,RO,RB]=meshRect(model,RO,RB);
     otherwise; 
       if strncmpi(RB.shape,'ls',2)% LsSphere,LsCyl, ...
         RB.shape=strrep(lower(RB.shape),'ls',''); 
         %lsutil('viewls',model,{RB})
         RC=stack_get(model,'info','LsutilOpt','get');
         if isfield(RB,'tolE');RC.tolE=RB.tolE;end
         if isfield(RB,'doCut'); model=feval(RB.doCut{:},model,{RB},RC);
         elseif isfield(RB,'NoSplit')
          def=lsutil('gen',model,{RB}); 
          i2=intersect(fix(def.DOF(def.def>0)),feutil(RB.onSurf,model));
          if isempty(i2);error('Problem');end
          [model.Node,model.Elt]=fevisco(['MakeSandwich -node ' RB.MakeSandwich],model, ...
           sprintf('selface & innode %s',sprintf('%i ',i2)));
         else
             model=lsutil('cut',model,{RB},RC);
         end
         [model.Elt,elt]=feutil('removeelt proid',model,RB.ProId);
         model.Elt=feutil('addelt',model.Elt,elt);
       else
         error('Shape%s unknown',RB.shape);
       end
     end
         
     if isfield(RB,'pl')&&~isempty(RB.pl);% Set materials
         model=feutil('setmat',model,RB.pl);
     end
     if isfield(RB,'delamply')
      %% #DelamPly extrude for delaminated ply and remove unused rigids -3
      RO.Delam=struct('PatchId',RO.PatchId,'MatId',RO.MatId,'delamply',RB.delamply);
      
     elseif isfield(RB,'Ecoef') % Assume composite with single material
       % model = d_shm('cbMatDef MatID coef',model)
       RB.pl=feutil(sprintf('getpl %i',RB.il(10)),model);
       ind=10:4:size(RB.il,2);
       RB.il(ind(RB.il(ind)==RB.pl(1)))=RB.il(1); % Modify MatId 
       RB.pl(1)=RB.il(1); 
       model=feutil('setpro',model,RB.il); % Set properties
       model=feutil('setmat',model,RB.pl); % Set properties
       model=d_shm(sprintf('CbMatDef %i %.15g',RB.pl(1),RB.Ecoef),model);
     elseif isfield(RB,'il');
      model=feutil('setpro',model,RB.il); % Set properties
     end
     RO.MatId=RO.MatId+1;RO.ProId=RO.ProId+1;
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
     [model,RO]=RefineToQuad(model,RO);
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
   if isfield(RB,'PostCb');model=feval(RB.PostCb{:},model,RB);end
 end % Loop on list
 if isfield(RO,'EltOri')
   model.Elt=feutil('addelt',RO.EltOri,model.Elt);
 end
 if isfield(model,'unit');
     model.pl=fe_mat(sprintf('convert%s;',model.unit),model.pl);
     model.il=fe_mat(sprintf('convert%s;',model.unit),model.il);
 end
 
 model=stack_set(model,'info','Electrodes',struct('data',RO.Electrode,'ver',1));
 if nargout==0;cf=feplot(model);fecom(cf,'colordatamat-edgealpha.1');
 else; out=model;
 end
 %% End of list building
end
[model,RO]=RefineToQuad(model,RO);
if isfield(RO,'Delam')
      model.Elt=feutil('Orient 1:1e4 n 0 0 1e5',model);
      RB.z=cumsum([0;RB.plies(:,2)])+RB.z0; 
      RB.newz=RB.z(1:RB.delamply+1);
      RB.oldz=RB.z(RB.delamply+1:end);
      error('Not yet implemented need shifting');
      sandCom=sprintf('makesandwich shell 0 0 %.15g shell %i %.15g %.15g', ...
      -(max(RB.oldz)-min(RB.oldz))/2+(max(RB.newz)-min(RB.newz))/2,RB.il(1)+1,0,0);
      m1=fevisco(sandCom,model); 
      model.Node=m1.Node; model.Elt=m1.Elt; 
      model.Elt=feutil('removeelt eltname rigid & withnode {proid}',model,RB.il(1));
    
      for j2=1:size(model.il,1)
       if ~strcmpi(fe_mat('typepstring',model.il(j2,2)),'p_shell.2');
           continue;
       end
       il=model.il(j2,:);r2=reshape(RB.plies(1:RB.delamply,:)',1,[]);
       i3=9+(1:length(r2));
       if norm(il(i3)-r2)<1e-10; il(i3)=[];il(1,size(model.il,2))=0;end % remove plies
       model.il(j2,:)=il; % Set back
      end
      il=[RB.il(1:9) r2];il(1)=il(1)+1;
      model=feutil('setpro',model,il);
 end

% Force vertical and x axis orientation for composite support
if RunOpt.NeedOrient
 sdtw('_ewt','obsolete needs check');
 % model.Elt=feutil('orient 1:100 n 0 0 1e5;',model);
 data=struct('sel','groupall','dir',{{1,0,0}},'DOF',[.01;.02;.03]);
 %[Case,CaseName]=fe_case(model,'getcase');
 model=p_shell('setTheta;',model,data);
end

if nargout==0; cf=feplot(model);fecom(cf,'colordatamat-edgealpha.1');
else; out=model;
end
out1=RO;


elseif comstr(Cam,'baseaccel')
%% #MeshBaseAccel: basic acceleromter mesh

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


elseif comstr(Cam,'piezoshaker')
%% #MeshPiezoShaker: calibration of a piezo with a shaker

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

elseif comstr(Cam,'tower')
%% #MeshTower: High-rise concrete building (central core)
% Data for the tower
L = 140 ;   % Height Tower
W = 20 ;    % Diameter Tower
th = 0.3 ;  % Thickness wall Tower

A = (W^2-(W-2*th)^2) ; MTOT = A*2200*L ;
I1 = W^4/12-(W-2*th)^4/12; I2 = I1 ;

%Create mesh
model=struct('Node',[1 0 0 0 0 0 0],'Elt',[]);
model=feutil('addelt',model,'mass1',1);
dz=feutil(sprintf('refineline %f',L/20),[0 L]);
model=feutil('extrude 0 0 1 0',model,dz);
model.unit='SI'; % go to m

%% Add properties
% Thick beam so k1 and k2 are important (shear correction factors). EB hyp
model.il=p_beam(sprintf('dbval 1 box %f %f %f %f',W,W,th,th));
model.pl=[1 fe_mat('m_elastic','SI',1) 30e9 0.2 2200 0 0.02 ];

% Keep only 2D vrtical DOFs (y, thetaxy) , clamp node 1,
model = fe_case(model,'FixDOF','Clamped',[1.01 1.02 1.06],'FixDOF','2-D motion',[.02 .03 .04 .05]);

% Assemble M and K matrices in model.K{1} and model.K{2}
model = fe_mk(model);

% Build Case to get Case.T...
[Case,model.DOF]=fe_mknl('init',model);

% Define excitation
nd=feutil('find node y==140',model);
model=fe_case(model,'DofLoad','F-top',nd+.01,1);

% Define sensor
model=fe_case(model,'SensDOF','d-top',{[num2str(nd) ':x']});
model.mdof=Case.DOF;
out=model;

elseif comstr(Cam,'langevin')
%% #MeshLangevin: Mesh a Langevin type transducer 

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



elseif comstr(Cam,'langconcrete')
%% #MeshLangConcrete: Langevin Transducer in concrete

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
    
elseif comstr(Cam,'ide')
%% #MeshIDEPatch : patch with inter digitated electrodes     

% Parameter handling
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   ' lx(400e-6#%g#"patch width")'...
   ' ly(300e-6#%g#"patch height")'...
   ' p0(700e-6#%g#"patch length")'...
   ' e0(50e-6#%g#"IDE width")'...
   ' nx(10#%g#"n_elt x direction")'...
   ' ny(5#%g#"n_elt y direction")'...
   ' nz(14#%g#"n_elt z direction")'...
   ' Quad(#3#"quadratic elements if present")'...
   ],{RO,CAM});

% Mesh
sp_util('epsl',1e-6); % Very small mesh in SI units check tolerance
model=feutil('objectbeam 1 1',[0 0 0;RO.lx 0 0],RO.nx);
model=feutil(['extrude ' num2str(RO.ny) ' 0 ' num2str(RO.ly/RO.ny) ' 0'],model);
model=feutil(['extrude ' num2str(RO.nz) ' 0 0 ' num2str(RO.p0/RO.nz) ],model);

% Define material properties
model.pl=m_piezo('dbval 1 -elas 10 SONOX_P502_iso');   
% Integration rules for volumes
model=p_solid('default;',model);
model.unit='SI';

% Build a MPC and defining ground and actuation DOF
InputDOF=[];
[model,InputDOF(end+1,1)]=p_piezo('ElectrodeMPC IDE bottom -ground',model, ...
   ['y==0 | y==' num2str(RO.ly) '& z<' num2str(RO.e0*1.01)]);
[model,InputDOF(end+1,1)]=p_piezo('ElectrodeMPC IDE top -input"Applied Voltage"' ...
   ,model,['y==0 | y==' num2str(RO.ly) '& z>' num2str(RO.p0-RO.e0*1.01)]);
%  [model,InputDOF(end+1,1)]=p_piezo('ElectrodeMPC IDE bottom -ground',model, ...
%      ['z<' num2str(RO.e0*1.01) '| z>' num2str(RO.p0-RO.e0*1.01) ' & y==0 ']);
%  [model,InputDOF(end+1,1)]=p_piezo('ElectrodeMPC IDE bottom top -input"Applied Voltage"' ...
%      ,model,['z<' num2str(RO.e0*1.01) '| z>' num2str(RO.p0-RO.e0*1.01) ' & y==' num2str(RO.ly)]);
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
 model=feutil('extrude 0 1 0 0',model,unique(dx));
 dy=[linspace(0,12,3) linspace(12,12+25,5) linspace(12+25,63,5) ...
    linspace(63,63+25,5) linspace(63+25,100,3)];
 model=feutil('extrude 0 0 1 0',model,unique(dy));
model.unit='mm';

 %%%%%% Material Properties
model.pl=m_elastic('dbval 1 Aluminum');
%%%%% Laminate properties
model.il=p_shell('dbval 1 -punitmm laminate 1 1.2 0');% this is to specify in mm
model=fe_case(model,'FixDof','Cantilever','x==0 -DOF 1:6');

% Add patches

RG.list={'Name','Lam','shape'
   'Main_plate', model,''  % Base structure
   'Act1', ... % name of patch
   'BaseId1 +Rect.Sonox_P502_iso.5525TH0_25 -Rect.Sonox_P502_iso.5525TH0_25', ... % Layout definition
     struct('shape','LsRect', ... % Remeshing strategy (lsutil rect here)
       'xc',15+55/2,'yc',12+25/2,'alpha',0,'tolE',.1)
    'Act2', ... % name of patch
   'BaseId1 +Rect.Sonox_P502_iso.5525TH0_25 -Rect.Sonox_P502_iso.5525TH0_25', ... % Layout definition
     struct('shape','LsRect', ... % Remeshing strategy (lsutil rect here)
       'xc',15+55/2,'yc',63+25/2,'alpha',0,'tolE',.1)
       };
RG.unit=model.unit;
mo1=d_piezo('MeshPlate',RG); 
mo1=p_piezo('DToSimple',mo1);

mo1.name='ULB_plate'; 
if ~isempty(strfind(Cam,'cantilever'));
    mo1=fe_case(mo1,'FixDof','Cantilever','x==0 -DOF 1:6');
end

if ~isempty(strfind(Cam,'feplot'))||nargout==0; 
    cf=feplot(mo1); fecom colordatapro-edgealpha.1
    if nargout>0;out=cf;end
else;out=mo1;
end

elseif comstr(Cam,'mfcplate')
%% #MeshMFCPlate : sample mesh of a plate with two MFCs
%  Note that d_piezo('MeshPlate') is much more general

%% generate the model - - - - - - - - - - - - - - - - - - - - - - - -

% Create mesh. Geometric properties in the manual
RO=struct('L',463,'w',50,'a',85,'b',28,'c',15,'d',11);

% create a rectangle with targetl = 3 mm
 model=feutil('objectquad 1 1',[0 0 0;1 0 0;0 1 0], ...
    feutil('refineline 5',[0 RO.c+[0 RO.a] RO.L]), ...
    feutil('refineline 5',[0 RO.d+[0 RO.b] RO.w]));
%%%%% Material Properties for supporting plate
 model.pl=m_elastic('dbval 1 -unit MM Aluminum'); % Aluminum
 model.il=p_shell('dbval 1 -punit MM laminate 1 1 0');
 model.unit='MM';
 
 
% ADD two 85x28mm MFC - P1 patches
RG.list={'Name','Lam','shape'
   'Main_plate', model,''  % Base structure
   'Act1', ... % name of patch
   'BaseId1 +SmartM.MFC-P1.8528 -SmartM.MFC-P1.8528', ... % Layout definition
     struct('shape','LsRect', ... % Remeshing strategy (lsutil rect here)
       'xc',RO.c+RO.a/2,'yc',RO.d+RO.b/2,'alpha',0,'tolE',.1)
       };
mo1=d_piezo('MeshPlate',RG);
%mo1=stack_rm(mo1,'info','Electrodes'); % xxxeb I am not sure why this is present ?

%p_piezo('tabpro',mo1);
 
if ~isempty(strfind(Cam,'cantilever'));
    mo1=fe_case(mo1,'FixDof','Cantilever','x==0 -DOF 1:6');
end

if ~isempty(strfind(Cam,'feplot'))||nargout==0; 
    cf=feplot(mo1); fecom colordatapro-edgealpha.1
    if nargout>0;out=cf;end
else;out=mo1;
end

elseif comstr(Cam,'triangleplate')
%% #MeshTrianglePlate : Mesh plate with triangular point load actuator

fname=sdtcheck('PatchFile',struct('fname','PzTrianglePlate.mat','in','piezo.zip','back',1));
% fname=fullfile(sdtdef('tempdir'),'sdtdemos/PzTrianglePlate.mat');
if exist(fname,'file')&&isempty(strfind(Cam,'reset'))
  fprintf('Loading %s\n',fname);
  load(fname); out=model;
else
%% Geometrical properties and reference nodes
L=414; w=314; a=100; b=33.58; c=(w-b)/2+50; 

% create reference nodes
mdl.Node = [1 0 0 0  0 0 0; 
          2 0 0 0  L 0 0; 
          3 0 0 0  L w 0;
          4 0 0 0  0 w 0 ;
          5 0 0 0  0 b+c 0;
          6 0 0 0  0 c  0;
          7 0 0 0  a c+b/2 0];
mdl.Elt=[];

mdl.Node(5:7,4)=10; % to remesh locally around the triangle

%% Meshing - Need gmsh

% create reference lines and mesh plate
mdl=fe_gmsh('AddLine -loop1',mdl,[1 2; 2 3; 3 4; 4 5; 5 7; 7 6; 6 1]);
mdl=fe_gmsh('AddLine -loop2',mdl,[5 7; 7 6; 6 5]);
mdl.Stack{end}.PlaneSurface=[1; 2]; mdl.Stack{1,3}.LineLoop{2,2}=[5 6 10];
model=fe_gmsh('write tmp.geo -lc 15 -run -2 -v 0',mdl);
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
model.il=p_shell('dbval 1 kirchhoff 1e-3 0', ...
    ['dbval 2 laminate  103 6e-05 0 101 1.8e-4 0 103 6e-05 0 1 1e-3  0 103 6e-05 0 101 1.8e-4 0 103 6e-05 0 ']);

model.il=p_piezo(model.il,'dbval 3 shell 2 100001    2   0 100002 6 0 ');


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

 %% New meshing script
 model=struct('Node',[1 0 0 0 0 0 0],'Elt',[]);
model=feutil('addelt',model,'mass1',1); 
% Note that the extrusion values are chosen to include the patch edges
dx=feutil('refineline 3',[0 50 100 300 350]);
dy=feutil('refineline 3',[0 25]);
model=feutil('extrude 0 1 0 0',model,dx);
model=feutil('extrude 0 0 1 0',model,dy);
model.unit='mm';

%%%%%% Material Properties
model.pl=m_elastic('dbval 1 steel');
%%%%% Laminate properties
model.il=p_shell('dbval 1 laminate 1 2e-3 0'); % Thickness should be specified in m

% Add patches
RG.list={'Name','Lam','shape'
   'Main_plate', model,''  % Base structure
   'Act1', ... % name of patch
   'BaseId1 +Rect.PIC_255.5025TH0_5 -Rect.PIC_255.5025TH0_5', ... % Layout definition
     struct('shape','LsRect', ... % Remeshing strategy (lsutil rect here)
       'xc',25,'yc',12.5,'alpha',0,'tolE',.1)
    'Act2', ... % name of patch
   'BaseId1 +Rect.PIC_255.5025TH0_5 -Rect.PIC_255.5025TH0_5', ... % Layout definition
     struct('shape','LsRect', ... % Remeshing strategy (lsutil rect here)
       'xc',75,'yc',12.5,'alpha',0,'tolE',.1)
       };
mo1=d_piezo('MeshPlate',RG); 
mo1=stack_rm(mo1,'info','Electrodes'); % old version should remove
mo1=fe_case(mo1,'FixDof','Cantilever','x==0 -DOF 1:6');
mo1=p_piezo('DToSimple',mo1);

mo1.name='Pz_shunt'; 

if ~isempty(strfind(Cam,'feplot'))||nargout==0; 
    cf=feplot(model); fecom colordatapro-edgealpha.1
    if nargout>0;out=cf;end
else;out=mo1;
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

elseif comstr(Cam,'simplebeam'); [CAM,Cam]=comstr(CAM,11);
    %% #MeshSimpleBeam : mesh a simple beam for sensor placement illustration
    % Create mesh
model=struct('Node',[1 0 0 0 0 0 0],'Elt',[]);      % Create a model and call the first node which is placed at (0,0,0) position
model=feutil('addelt',model,'mass1',1); 
dx=[0:1:350]';
dy=[0:2.5:25]';

model=feutil('extrude 0 1 0 0',model,dx);
model=feutil('extrude 0 0 1 0',model,dy);

%Switch to SI units
model.Node(:,5:7)=model.Node(:,5:7)/1000; model.unit='SI';

% Define material properties
model.pl=m_elastic('dbval 1 steel');

% Define section properties (including the thickness)
model.il=p_shell('dbval 1 laminate 1 2e-3 0');

% Boundar Conditions
model=fe_case(model,'FixDof','Cantilever','x==0 -DOF 1:6');

model=stack_set(model,'info','DefaultZeta',0.004);
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

elseif comstr(Cam,'baffledpiston'); [CAM,Cam]=comstr(CAM,14);
    %% #MeshBaffledPiston : mesh a circular piston -2
    


% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
if carg<=nargin&&isstruct(varargin{carg});RO=varargin{carg};carg=carg+1;
else; RO=struct;
end

if isempty(CAM)&&isempty(fieldnames(RO))
%% Generic meshing of pre-defined patches 
RO.Rev=20;
RO.nelt=10;
CAM='';  
end
 
% Parameter handling
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   ' Rev(20#%g#"Number of revolutions")'...
   ' selt(10e#%g#"Size of element along radius in mm")'...
   ],{RO,CAM});

model=struct('Node',[1 0 0 0 0 0 0],'Elt',[]);
model=feutil('addelt',model,'mass1',1); 
% 
%dx=feutil('refineline 10',[0 100]);
dx=feutil(sprintf('refineline %i',RO.selt),[0 100]);
model=feutil('extrude 0 1 0 0',model,dx);
%model=feutil('Rev 20 o 0 0 0 360 0 0 1',model);
model=feutil(sprintf('Rev %i o 0 0 0 360 0 0 1',RO.Rev),model);
model.Node(:,5:7)=model.Node(:,5:7)/1000; model.unit='SI'; % go to m
% 
%  %%%%%% Material Properties for supporting plate
  model.pl=m_elastic('dbval 1 Aluminum'); % Aluminum
  model.il=p_shell('dbval 1 laminate 1 1e-3 0');
  out=model;



%% #MeshEnd 
else;error('Mesh%s unknown',CAM);
end

%% #Mat -------------------------------------------------------------------
elseif comstr(Cam,'mat');[CAM,Cam]=comstr(CAM,4);

if comstr(Cam,'pic255a')
%% #MatPIC255a : compute full material properties from datasheet (2013) for PIC255
% PIC 255 data 2013
in.eps0=8.854e-12;
in.d33=400e-12; in.d31=-180e-12; in.d15=550e-12;
in.g31=-11.3e-3; in.g33=25e-3;
in.e33t=1750*in.eps0; in.e11t=1650*in.eps0;
in.s11E=16.1e-12; in.s33E=20.7e-12; in.c33D=11e10;
in.k33=0.69; in.k15=0.66; in.k31=0.35; in.kp=0.62; in.kt=0.47;
in.rho=7800;

%Data computed from datasheet info
ou.Ep=1/in.s11E; ou.Ez=1/in.s33E;
tp.s12E=-in.s11E + 2*in.d31^2/(in.kp^2*in.e33t);
ou.vp=-ou.Ep*tp.s12E;
ou.Gp=ou.Ep/(2*(1+ou.vp));
tp.s55E=in.d15^2/(in.e11t*in.k15^2); ou.Gzp=1/tp.s55E;

ou.vzp=0.3; % (hopothesis)
ou.vpz=ou.Ep/ou.Ez*ou.vzp;

% %Checking
% k31b=sqrt(d31^2/(e33t*s11E));
% k33b=sqrt(d33^2/(e33t*s33E));
% g31b=d31/e33t;
% g33b=d33/e33t;

% in = input values, out=output values, tp=temp values
out.in=in; out.out=ou; out.tp=tp;

elseif comstr(Cam,'pic181a')
%% #MatPIC181a : compute full material properties from datasheet (2013) for PIC181
% PIC 181 data 2013
in.eps0=8.854e-12;
in.d33=265e-12; in.d31=-120e-12; in.d15=475e-12;
in.g31=-11.2e-3; in.g33=25e-3;
in.e33t=1250*in.eps0; in.e11t=1100*in.eps0; %Datasheet value is 1500 but should be lower than e33, so corrected (how ?)
in.s11E=11.8e-12; in.s33E=14.2e-12; in.c33D=16.6e10;
in.k33=0.66; in.k15=0.63; in.k31=0.32; in.kp=0.56; in.kt=0.46;
in.rho=7800;

%Data computed from datasheet info
ou.Ep=1/in.s11E; ou.Ez=1/in.s33E;
tp.s12E=-in.s11E + 2*in.d31^2/(in.kp^2*in.e33t);
ou.vp=-ou.Ep*tp.s12E;
ou.Gp=ou.Ep/(2*(1+ou.vp));
tp.s55E=in.d15^2/(in.e11t*in.k15^2); ou.Gzp=1/tp.s55E;

ou.vzp=0.3; % (hopothesis)
ou.vpz=ou.Ep/ou.Ez*ou.vzp;

% %Checking
% k31b=sqrt(d31^2/(e33t*s11E));
% k33b=sqrt(d33^2/(e33t*s33E));
% g31b=d31/e33t;
% g33b=d33/e33t;

% in = input values, out=output values, tp=temp values
out.in=in; out.out=ou; out.tp=tp;

elseif comstr(Cam,'pic181b')
%% #MatPIC181b : compute full material properties from datasheet (2013) for PIC181
% PIC 181 data 2023
in.eps0=8.854e-12;
in.d33=265e-12; in.d31=-120e-12; in.d15=475e-12;
in.g31=-10.6e-3; in.g33=23e-3;
in.e33t=1250*in.eps0; in.e11t=1200*in.eps0; %Datasheet value is 1500 but should be lower than e33, so corrected (how ?)
in.s11E=11.8e-12; in.s33E=13.3e-12; in.c33D=17e10;
in.k33=0.66; in.k15=0.63; in.k31=0.32; in.kp=0.56; in.kt=0.46;
in.rho=7850;

%Data computed from datasheet info
ou.Ep=1/in.s11E; ou.Ez=1/in.s33E;
tp.s12E=-in.s11E + 2*in.d31^2/(in.kp^2*in.e33t);
ou.vp=-ou.Ep*tp.s12E;
ou.Gp=ou.Ep/(2*(1+ou.vp));
tp.s55E=in.d15^2/(in.e11t*in.k15^2); ou.Gzp=1/tp.s55E;

ou.vzp=0.3; % (hopothesis)
ou.vpz=ou.Ep/ou.Ez*ou.vzp;

% %Checking
% k31b=sqrt(d31^2/(e33t*s11E));
% k33b=sqrt(d33^2/(e33t*s33E));
% g31b=d31/e33t;
% g33b=d33/e33t;

% in = input values, out=output values, tp=temp values
out.in=in; out.out=ou; out.tp=tp;

elseif comstr(Cam,'pic255b')
%% #MatPIC255b : compute full material properties from datasheet (2023) for PIC255
% PIC 255 data 2023
in.eps0=8.854e-12;
in.d33=400e-12; in.d31=-180e-12; in.d15=550e-12;
in.g31=-11.8e-3; in.g33=25e-3;
in.e33t=1800*in.eps0; in.e11t=1750*in.eps0;
in.s11E=16e-12; in.s33E=19e-12; in.c33D=15.4e10;
in.k33=0.69; in.k15=0.65; in.k31=0.35; in.kp=0.62; in.kt=0.47;
in.rho=7800;

%Data computed from datasheet info
ou.Ep=1/in.s11E; ou.Ez=1/in.s33E;
tp.s12E=-in.s11E + 2*in.d31^2/(in.kp^2*in.e33t);
ou.vp=-ou.Ep*tp.s12E;
ou.Gp=ou.Ep/(2*(1+ou.vp));
tp.s55E=in.d15^2/(in.e11t*in.k15^2); ou.Gzp=1/tp.s55E;

ou.vzp=0.3; % (hopothesis)
ou.vpz=ou.Ep/ou.Ez*ou.vzp;

% %Checking
% k31b=sqrt(d31^2/(e33t*s11E));
% k33b=sqrt(d33^2/(e33t*s33E));
% g31b=d31/e33t;
% g33b=d33/e33t;

% in = input values, out=output values, tp=temp values
out.in=in; out.out=ou; out.tp=tp;

elseif comstr(Cam,'ferroperm')
%% #MatFerroperm : compute full material properties from xls sheet
% file='Ferroperm MatData.xls', 2018
% file='Ferroperm_materials_data_sep2023.xlsx', 2023
file=varargin{2}; [a,b,c]=xlsread(file);
ind=3:size(c,2)-2;
names=c(3,ind);
ind=ind-2;

s11E=a(36,ind);
nu12=a(34,ind);
nu13=a(35,ind);
s13E=a(38,ind);
s33E=a(39,ind);
s44E=a(40,ind);
s55E=s44E;
rho0=a(33,ind);
eps1t=a(1,ind);
eps3t=a(2,ind);
d310=a(12,ind);
d330=a(13,ind);
d150=a(14,ind);

for ij=1:length(names)

% Mechanical properties
Ep=1/s11E(ij); % Corresponds to Y11E
Ez=1/s33E(ij); % Corresponds to Y33E - Ep is usually larger than Ez, this is not coherent ... OK seems to depend on the material ...
nup=nu12(ij); % corresponds to nu12
nuzp=nu13(ij); % It seems that nu13 in datasheet is in fact nu31 (from the values of s11E s33E and s13E)
nupz=nuzp*Ep/Ez;
Gp=Ep/(2*(1+nup));
Gzp=1/s44E(ij);
Gpz=1/s55E(ij);
rho=rho0(ij);
% Permittivity
e1t=eps1t(ij);
e3t=eps3t(ij);
% Piezoelectric properties
d31=d310(ij);
d33=d330(ij);
d15=d150(ij);

d=zeros(1,18); d([11 13])=d15; d([3 6])=d31; d(9)=d33;
  eps=zeros(1,9); eps([1 5])= 8.854e-12*e1t; eps(9)=8.854e-12*e3t;
out(ij).pl=[1 fe_mat('m_piezo','SI',2) 2 d  eps;
    2  fe_mat('m_elastic','SI',6) Ep Ep Ez nupz nuzp nup Gpz Gzp Gp rho zeros(1,18)]; 
out(ij).name=names{ij};
out(ij).type='m_piezo';
out(ij).unit='SI';

end


%% #MatEnd
else;error('Mat%s unknown',CAM);
end


%% #Figstyles --------------------------------------------------------------

elseif comstr(Cam,'setplotwd');[CAM,Cam]=comstr(CAM,7);
%% #SetPlotWD : defines default dir for figures -2
 wd=fullfile(fileparts(which('d_piezo')),'../piezo/figures');
 if ~exist(wd,'dir');wd=fullfile(sdtdef('tempdir'),'plots');sdtkey('mkdir',wd);end
 sdtroot('setproject',struct('PlotWd',wd,'Report','')); % non empty report would open word

elseif comstr(Cam,'definestyles');[CAM,Cam]=comstr(CAM,11);
    %% #DefineStyles: defines default styles for feplot and iiplot figures -2
    
sdtroot('SetOsDic',{
    'cms_wide',{'Position',[NaN NaN 800 600]} % ImSw100 ImSw80{@line,""} : removes the @line part of the style ImSw80
    'cms_tall',{'Position',[NaN NaN 600 800]}
    'cms_tw',{'Position',[NaN NaN 800 800]}
    'cms_f14' ,{'@axes',{'fontsize',14},'@title',{'fontsize',14},'@ii_legend',{'fontsize',14}}
    'cms_grid',{'@axes',{'xgrid','on','ygrid','on','zgrid','on'}}});
li=sdtroot('cbosdicget',[],'cms_grid');
disp(sdtm.toString(sdtroot('cbosdicget',[],'ImSw80')));
% d_piezo('reset') ProjectWd='@tempdir/sdtdemos/piezo' 
% To change transparency
sdtroot('setosdic',{'FiAlpha',{'@PlotInfo',{'sel-linface','','showpatch','','coloredgealpha.5','','colorfacealpha.1',''}}})

elseif comstr(Cam,'setstyle');[CAM,Cam]=comstr(CAM,8);
  %% #SetStyle: apply style to a feplot or iiplot -2
    
 curfh=varargin{2};
 d_piezo('DefineStyles'); % Not very smart but robust

    %Apply style to figure for output % xxxAD needs to be replaced by
    %cf.osd calls
    if curfh.opt(1,3)==3 % setslyle.feplot This is a feplot
        % xxxAd
     cf=curfh;
     cf.osd_('cms_wide','LgMl-FontSize14','FnI') % size/ratio + Legend type (LgMl=model+defo title) 14pt + Output file name based on ii_legend
     % LgMl{Fontsize,14,FontFace,bold} % xxx a faire 
     feplot(cf);
    
     % for feplot add triax with 14 pt text
     fecom('triaxon','@text',{'Fontsize',14});

     
    else % setstyle.iiplot
     ci=curfh;
     ci.osd_('cms_wide','cms_grid')
     ci.ua.axProp={'@OsDic',{'cms_f14'}};% Needs to be done for each channel
     
     ci.os_({'@axes',{'linewidth',1},'@line',{'linewidth',2}})
     comgui('imftitle',ci); iiplot % for refresh
     ci.os_('@OsDic(SDT Root)',{'FnIy2'})

    end    
    
    %xxxAD : this should be predefined using styles 
    if size(varargin)>2
     fname=varargin{3};
         cingui('plotwd',curfh,'FileName',{'@Plotwd',fname,'.jpeg'},...
     'printOpt',{'-r600'}) ;
    end
    
    % No change of style for figure output
    curfh.os_('FgWysiwyg');
    %% #FigstylesEnd -2

%% #Oldstyles (obsolete) ---------------------------------------------------
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
 elseif comstr(Cam,'imxfa1') 
 % #ImXFA1 : basic display of transfers -2
   out={'@figure',{'position',[NaN,NaN,800,481]},'@exclude',{'legend.*'}, ...
       '@text',{'FontSize',14},'@axes',{'FontSize',14,'box','on','position',[.13 .15 .775 .815]}, ...
       '@xlabel',{'FontSize',14},'@ylabel',{'FontSize',14}, ...
       '@zlabel',{'FontSize',14},'@title',{'FontSize',14}};
   %  @line,{'linewidth',2}
 elseif comstr(Cam,'imxfa') 
 % #ImXFA : basic display of transfers -2
   out={'@figure',{'position',[NaN,NaN,800,481]},'@exclude',{'legend.*'}, ...
       '@text',{'FontSize',14},'@axes',{'FontSize',14,'box','on'}, ...
       '@xlabel',{'FontSize',14},'@ylabel',{'FontSize',14}, ...
       '@zlabel',{'FontSize',14},'@title',{'FontSize',14}};
 elseif comstr(Cam,'imraise')
 % #Imraise -2
    out={'@EndFcn','ga=findall(gf,''type'',''axes'');for j1=1:length(ga);ga(j1).Position=ga(j1).Position+[0 .05 0 0];end'};
 elseif comstr(Cam,'imfa') 
 % #ImFa feplot for Cassem examples -2
   out={'@figure',{'position',[NaN,NaN,600,450]},'@exclude',{'legend.*'}, ...
       '@ii_legend',{'FontSize',14,'interpreter','none'}};

 elseif comstr(Cam,'imcheck') 
%% #Imcheck options -2
RO=varargin{carg};carg=carg+1;
if ~isfield(RO,'PlotWd');
  try; RO.PlotWd=oscarui('wd@PlotWd');
  catch; error('You must provide a PlotWd');
  end
end
if ~isfield(RO,'HtmWidth'); RO.HtmWidth=[];end % Default width
if ~isfield(RO,'RelPath'); RO.RelPath=2;end % Path for file history

out=RO;

%% #OldstylesEnd
 else; error('%s unknown',CAM);
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
%% #SolvePatchCapacity : reduction for capacity resolution -2

model=[];Sens=[];eval(iigui({'model','Case','Load','out','Sens','RunOpt'},'MoveFromCaller'));
m=feval(fe_reduc('@get_mat'),model.K,2,model,RunOpt);
k=feval(fe_reduc('@get_mat'),model.K,1,model,RunOpt);

Freq=stack_get(model,'info','Freq','get');
if ~isempty(Freq)
  %% FrfSnapShot % Called from fe_reduc
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
%% #SolvePatchDfrf : capacity of an isolated patch - 2

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
%% #Pcond: rescales electrical DOFs for better conditionning
% model=fe_case(model,'pcond','Piezo','d_piezo(''Pcond 1e8'')');
%  default coef is 1e8
[CAM,Cam,r1]=comstr('cond',[-25 2],CAM,Cam); if isempty(r1);r1=1e8;end
 Case=evalin('caller','Case');model=evalin('caller','model');
 T=evalin('caller','T');
 pc=ones(size(Case.T,1),1);
 i1=fe_c(model.DOF,.21,'ind'); 
% i2=fe_c(model.DOF,.19,'ind'); if ~isempty(i2);pc(i2)=1e-3;disp('xxxd_piezo');end
 if isfield(Case,'TIn')
     i1(any(Case.TIn(i1,:),2))=[];pc(i1)=r1;
     pc=diag(sparse(pc)); T=pc*T; Case.pc=pc;
     Case.TIn=pc*Case.TIn; assignin('caller','Case',Case);
 else;
     pc(i1)=r1; pc=diag(sparse(pc)); T=pc*T;
 end
 assignin('caller','T',T);
 %z=diag(k);mean(z(fe_c(model.DOF,.21,'ind')))/mean(z(fe_c(model.DOF,.21,'ind',2)))

 

 
%% #TutoDo: recover model from a specific tuto step
elseif comstr(Cam,'tuto'); 
 eval(sdtweb('_tuto',struct('file','d_piezo','CAM',CAM)));
 if nargout==0; clear out; end
elseif comstr(Cam,'@');out=eval(CAM);
elseif comstr(Cam,'cvs')
 out=sdtcheck('revision','$Revision: c3f9e0b $  $Date: 2021-03-16 18:05:59 +0100 $ ');
else; error('%s unknown',CAM);
end 
 
end


function [model,RO,RB]=meshRect(model,RO,RB);
%% #fct-MeshRect: Insert a rectangle
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
if RB.lx==0||RB.ly==0; RB=struct; return; end
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
%% #fct-MeshCirc: Add circular patches
% RO global options, RB current options of patch
% RO.Circ=[xc yc rc lc (MatId ProId)] 
if RB.rc==0; RB=struct; RO.RefineToQuad=1; return; end
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
 if strcmpi(RB.shape,'annulus')
  RO.PatchId=RO.ProId;
  for j1=3:2:length(RB.OD)
   mo1=feutil(sprintf('objectcircle %.15g %.15g 0  %.15g 0 0 1 %i', ...
     RB.xc,RB.yc,RB.OD(j1)/2,2*pi*RB.rc/RB.lc),mo1);
  end
 end
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
 if strcmpi(RB.shape,'annulus')
  r3=[2 0;1 -2];
  %mo2=fe_gmsh('AddFullCircle -loop3',mo2,[n1;n1+[RB.ID/2 0 0]; 0 0 1]);  
  %GM=stack_get(mo2,'geom','GMSH','get');GM.PlaneSurface=[3 0;2 -3;1 -2];
  for j1=3:2:length(RB.OD)
   j2=(j1+1)/2+1;
   mo2=fe_gmsh(sprintf('AddFullCircle -loop%i',j2), ...
       mo2,[n1;n1+[RB.OD(j1)/2 0 0]; 0 0 1]);  
   r3=[j2 0;r3(1) -j2;r3(2:end,:)];
  end
  GM=stack_get(mo2,'geom','GMSH','get');GM.PlaneSurface=r3;
 else
  GM=stack_get(mo2,'geom','GMSH','get');GM.PlaneSurface=[2 0;1 -2];
 end
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
 if strcmpi(RB.shape,'annulus')
  [EGroup,nGroup]=getegroup(mo2.Elt);
  mo2.Elt=feutil(sprintf('set group %i matid 1 proid 1',nGroup),mo2); % patch
  i3=RB.OD(2:2:end); RO.iCyc={};
  for j1=1:length(i3);
   if i3(j1)==0; 
    mo2.Elt=feutil(sprintf('set group %i matid 1 proid 1',nGroup-j1),mo2); % patch
    RO.Cyc{j1,2}=1;
   else
    RO.Cyc{j1,2}=RO.ProId-max(i3)+i3(j1)+1;   
    mo2.Elt=feutil(sprintf('set group %i matid %i proid %i',nGroup-j1, ...
        RO.Cyc{j1,2}*[1 1]),mo2); % patch 
   end
   RO.Cyc{j1,3}=RB.OD(j1*2-1)/2;
   RO.Cyc{j1,1}=feutil('findnode group & group ',mo2,nGroup-j1+1,nGroup-j1);
  end
 else
  mo2.Elt=feutil('set group2 matid 1 proid 1',mo2); % patch
  mo2.Elt=feutil(sprintf('set group1 matid %i proid %i',RO.MatId,RO.ProId),mo2);
 end
 model=feutil('addtest-merge;',model,mo2);
 if isfield(RO,'Cyc');
   i2=double(stack_get(model,'info','OrigNumbering','get'));
   i2=sparse(i2(:,1),1,i2(:,2));
   RO.Cyc(:,1)=cellfun(@(x)full(i2(x)),RO.Cyc(:,1),'uni',0);
 end
 model=stack_rm(model,'info','OrigNumbering');
 if ~isfield(RO,'Circ');RO.Circ=[];end
 RO.Circ(end+1,1:6)=[RB.xc RB.yc RB.rc RB.lc RO.MatId RO.ProId];
end
%% #fct-RefineToQuad: Refine to second degree elements
function [model,RO]=RefineToQuad(model,RO);

if isfield(RO,'RefineToQuad');
   r2=feutil('refineToQuad;',model);model.Node=r2.Node;model.Elt=r2.Elt;
   NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
   for j1=1:size(RO.Circ,1) % make true circles
    if size(RO.Circ,2)<5;break;end
    if isfield(RO,'Cyc')&&~isempty(RO.Cyc)
     for j2=1:size(RO.Cyc,1)
       n1=feutil('selelt proid & seledge',model,RO.Cyc{j2,2});
       n1=feutil(sprintf('getnode inelt{ withnode %s}', ...
          sprintf('%i ',RO.Cyc{j2,1})),model.Node,n1);
       r2=RO.Circ(j1*ones(size(n1,1),1),1:2);r1=n1(:,5:6)-r2;   
       r1=diag(sparse((RO.Cyc{j2,3}./sqrt(r1.^2*[1;1]))))*r1;
       model.Node(NNode(n1(:,1)),5:6)=r1+r2;
     end
    else
     n1=feutil(sprintf('getnode matid1 & matid %i& proid %i', ...
        RO.Circ(j1,5:6)),model);
     r2=RO.Circ(j1*ones(size(n1,1),1),1:2);r1=n1(:,5:6)-r2;
     if size(RO.Circ,2)>6&&RO.Circ(j1,7) % annulus
      i2=sqrt(sum(r1.^2,2))>.5*RO.Circ(j1,3)+.25*RO.Circ(j1,7);
      r1(i2,:)=diag(sparse((RO.Circ(j1,3)./sqrt(r1(i2,:).^2*[1;1]))))*r1(i2,:);
      i2=~i2;
      r1(i2,:)=diag(sparse((RO.Circ(j1,7)/2./sqrt(r1(i2,:).^2*[1;1]))))*r1(i2,:);
     else
      r1=diag(sparse((RO.Circ(j1,3)./sqrt(r1.^2*[1;1]))))*r1;
     end
    end
    model.Node(NNode(n1(:,1)),5:6)=r1+r2;
   end
   RO=feutil('rmfield',RO,'RefineToQuad');
end
end




function SE=modal()
%% #modal : just modal 

% Transfers model to function
eval(iigui({'model','Case'},'MoveFromCaller'))
eval(iigui({'RunOpt'},'GetInCaller'))

% clean up model to remove info about matid/proid
SE=feutil('rmfield',model,'pl','il','unit','Stack');
%dbstack; keyboard

% Compute modeshapes
def=fe_eig({model.K{1},model.K{2},Case.T,Case.mDOF},stack_get(model,'info','Eigopt','g'));

% reduce matrics based on N modes
SE.K = feutilb('tkt',def.def,model.K);
SE.DOF=[1:size(SE.K{1},1)]'+.99; % .99 DOFS
SE.Klab={'m' 'k' 'c' '4'};
SE.TR=def; SE.TR.adof=SE.DOF;

% assign to output
assignin('caller','T',SE)
assignin('caller','Case','');


end