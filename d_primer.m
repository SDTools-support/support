function [out,out1]=d_primer(varargin); %#ok<*NOSEM,*STOUT>

% d_primer sample meshes used for support
%
% Use 
%     sdtweb('_taglist','d_primer') to view current contents
%     d_primer('tuto')  % to see integrated tutorials

%       Copyright (c) 1990-2021 by SDTools, All Rights Reserved.
%       For revision information use d_primer('cvs')

if nargin==0; d_primer('tuto'); return; end

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*DEFNU>

%% Argument parsing : generates CAM (command) and RO
RO=struct('info',1);
if ~ischar(varargin{1});obj=varargin{1};evt=varargin{2};
    [CAM,Cam]=comstr(varargin{3},1);carg=4;
    RO=fe_range('ToValue',obj,evt);
else; obj=[];evt=[];[CAM,Cam]=comstr(varargin{1},1);carg=2;
end
if carg>nargin||~isstruct(varargin{carg});
else;RO=varargin{carg};carg=carg+1;
end


%% #Truss : tutorials associated with truss mesh -----------------------------
if comstr(Cam,'mesh');[CAM,Cam]=comstr(CAM,6);

if comstr(Cam,'feutil2')
%% #TutoFeutil2 feutil example 2 - - - - - - - - - - - - - - - - - -2
% source primer.tex

%% Step1 Direct finite element model declaration and basic uses

% This is a direct declaration of a 2 bay truss model, defining nodes and
% elements.

model=struct('Node',[],'Elt',[]); % initialize empty model
%          NodeID   unused   x y z, see sdtweb('node')
model.Node=[ 1      0 0 0    0 1 0;
             2      0 0 0    0 0 0;
             3      0 0 0    1 1 0;
             4      0 0 0    1 0 0;
             5      0 0 0    2 0 0;
             6      0 0 0    2 1 0;
             7      0 0 0    1 1 1]; % reference node

model.Elt=[ % see sdtweb('elt')
 % declaration of element group for longerons
    Inf     abs('beam1') % see sdtweb('beam1') 
 %node1 node2   MatID ProID nodeR ,  zeros to fill the matrix 
    1       3      1    1     7       0
    3       6      1    1     7       0
    2       4      1    1     7       0
    4       5      1    1     7       0
 % declaration of element group for diagonals
    Inf     abs('beam1')
    2       3      1    2     7       0
    4       6      1    2     7       0
 % declaration of element group for battens
    Inf     abs('beam1')
    3       4      1    3     7       0
    5       6      1    3     7       0];

%% Now display the result
cf=comgui('guifeplot -project "SDT Root"',3); % Robust open in figure(3)
cf.model=model;       % create feplot axes
fecom(';view2;textnode;triax;'); % manipulate axes 

%% Step2 FeutilMesh ----------------------------------------------------------
% In practice you will generally prefer to generate the model using FEUTIL
% in this case this can be done as follows

 model=struct('Node',[],'Elt',[]); % initialize empty model

 % Build the first bay
 model.Node=[1 0 0 0  0 0 0;2 0 0 0    0 1 0;
             3 0 0 0  1 0 0;4 0 0 0    1 1 0];
 model.Elt=feutil('ObjectBeamLine 1 3 0 2 4 0 3 4 0 1 4');
 
 % Put it in the main model, translate, add the second bay, see result
 mo1=model; % Keep initial bay
 mo1=feutil('TransSel 1 0 0',mo1); % translate the 1st bay
 model=feutil('AddTestCombine',model,mo1);
 
 % Display
 feutil('info',model) % information on model
 cf=feplot(model); % display in feplot

%% Step3 : use femesh to create a more complex structure
% The following creates a 3-D truss based on the same basic cell.

 model.Elt=feutil('RemoveElt group2',model);
 model.Elt=feutil('Divide group 1 InNode 1 4',model);
 model.Elt=feutil('SetGroup1 name bar1',model);
 
 mo1=model;
 mo1=feutil('RepeatSel 10 1 0 0',mo1); % repeat and translate the 1st bay
 model=feutil('AddTestCombine',model,mo1);
 
 mo1=feutil('Rotatesel 1 60 1 0 0',mo1);
 model=feutil('AddTestCombine',model,mo1);

 mo1=model; mo1.Elt=feutil('SelElt group3:4',mo1);
 mo1=feutil('RotateSel 2 -60 1 0 0',mo1);
 model=feutil('AddTestCombine',model,mo1);
 
 
%% Step4 : discuss beam orientation
 model.Elt=feutil('SelElt group3:8',model);
 i1=feutil('FindElt eltname beam1',model);
 model.Elt(i1,5)=0; % vx
 model.Elt(i1,6)=0; % vy
 model.Elt(i1,7)=1; % vz
 
 % Define properties:
 model.il=p_beam('dbval 1 bar .02 .05');
 model.pl=m_elastic('dbval 1 Steel');
 model.unit='SI'; % unit of the mesh
 model.name='Truss';

 % Display result in feplot:
 cf=feplot(model);
 fecom(';triaxon;view3;view y+180');
 
 % now verify beam sections using a 3D view
 p_beam('Show3D',cf); feplot
 iimouse('view',gca,[ -22.47 -34.28 25.82 4.255 0.5565 0.4713 0.30 0.40 0.87 1.09]);

% EndTuto

elseif comstr(Cam,'mesh')
%% #TutoMesh2 : plate example in SDT primer  - - - - - - - - - - -1
% TutoSource primer.tex#TutoMesh2

%% Step1 init geometry
md0=struct('Node',[],'Elt',[],'name','primer','unit','SI'); % init model
mdl0.Node=[1 0 0 0 0 0 0 ;
           2 0 0 0 0 1 0 ];
mdl0.Elt=feutil('ObjectBeamLine 1 2'); % create beam between node 1 and 2
%% Step2 generate a plate in the y direction
mdl0=feutil('Extrude 6 2 0 0',mdl0); % extrude the beam along x as 6 quads
mdl0=feutil('RepeatSel 2 0 1 0',mdl0); % repeat mdl0
model=mdl0; % first step of final model
%% Step 3 add a second plate oriented vertically
mdl0.Elt=feutil('ObjectBeamLine', feutil('FindNode x>=0 & y==0',model)); % edge
mdl0=feutil('Extrude 2 0 0 1',mdl0); % extrude the edge
model=feutil('AddTestMerge',model,mdl0); % add mdl0 to model
cf=feplot(model);fecom('ShowPatch'); % plot model

%% EndTuto 

elseif comstr(Cam,'beambar')
%% #TutoBeamBar : illustrate mixed beam/bar using FEMESH - - - - - - - - - -1

%% Step1 define the mesh

% Define nodes sdtweb('node')
FEnode = [1  0 0 0  0 0 0;2  0 0 0  1 0 0;
          3  0 0 0  0 1 0;4  0 0 0  1 1 0];
% make an outer frame with beam elements
% see sdtweb('femesh')   sdtweb('elt')
femesh(';object beamline 1 2 0 4 3 1;addsel;info FEelt');

% make diagonals with bar elements
femesh(';object beamline 1 4 0 2 3;addsel;set group2 name bar1');

% repeat the cell 10 times and keep the result
femesh(';selgroup1:2;repeatsel 10 1 0 0;');
FEelt=FEel0;

% add a beam at the end of the last bay
femesh('objectbeamline',femesh('findnode x==10'));
femesh addsel;


%% Step2 Define the properties of both groups - - - - - - - - - - - - - - -

% Set the MatID and SecID values for elements
femesh(';setgroup1 3 mat1sec1;setgroup2 mat2 sec2');
model=femesh;  % Get a copy of the model

% Display the model, the model is accessible by cf.mdl

% Set the materials
model.pl=[m_elastic('dbval 1 aluminum');
           m_elastic('dbval 2 steel')];
       

cf=feplot(model); 
fecom('colordatamat') % Display using material color

% Set the beam properties          
cf.mdl.il=[p_beam('dbval 1 rectangle .1 .1')
           p_beam('dbval 2 circle .03')];
fecom('colordatapro') % Display properties

% Open model properties figure
fecom('promodelinit')

%% Step3 compute modes
% Set eigen value solver options, method 6, 20 modes, shift 1e3 
cf.mdl=stack_set(cf.mdl,'info','EigOpt',[6 20 1e3]);

% compute modes and display the result
cf.def= fe_eig(cf.mdl);
fecom(';scaledef1;view2;ch4');

% EndTuto

else;error('Truss%s unknown',CAM);
end

elseif comstr(Cam,'pretest')
%% #TutoPreTest : illustrate geometry definition
% TutoSource primer.tex#TutoPreTest

%% Step 1. Define nodes
cf=comgui('guifeplot -reset',2);

node=[1001 0 0 0   0 0 0  ;  1003 0 0 0 0.2 0 0
      1007 0 0 0 0.6 0 0  ;  1009 0 0 0 0.8 0 0
      1013 0 0 0 1.2 0 0  ;  1015 0 0 0 0 0.2 0
      1016 0 0 0 0.2 0.2 0;  1018 0 0 0 0.6 0.2 0
      1019 0 0 0 0.8 0.2 0;  1021 0 0 0 1.2 0.2 0
      1029 0 0 0 0 0 0.2  ;  1030 0 0 0 0.2 0 0.2
      1032 0 0 0 0.6 0 0.2;  1033 0 0 0 0.8 0 0.2
      1035 0 0 0 1.2 0 0.2];
fecom('AddNode',node)

%% Step 2. Define connectivity
% define straight edges
L=[1001 1003 1007 1009 1013]; fecom('AddLine',L);
L=[1015 1016 1018 1019 1021]; fecom('AddLine',L);
L=[1029 1030 1032 1033 1035]; fecom('AddLine',L);
% Define 5 L shaped edges as single 4th group
L=[1015 1001 1029 0 1016 1003 1030 0 1018 1007 1032 0 ...
   1019 1009 1033 0 1021 1013 1035 0]; fecom('AddLine',L);

%% Step 3. Define and show sensors
 sdof=[35.02 ; 21.03 ; 18.03 ; 32.02 ; 19.03 ; 33.02 ;
       16.03 ; 30.02 ; 35.03 ; 21.02]+1000;
cf.mdl=fe_case(cf.mdl,'SensDof','Test',sdof);
fecom('CurtabCases -ViewOn','Test');
fecom('TriaxOn');fecom('TextNode','GroupAll','FontSize',12)

%% EndTuto

elseif comstr(Cam,'kart')
%% #TutoKart : kart identification example
% TutoSource primer.tex#TutoKart

% The data can be downloaded with
demosdt('download http://www.sdtools.com/contrib/kart_example.mat')
% www.sdtools.com/contrib/kart_example.mat
%----------

% WIRE FRAME DEFINITION (present in WIREFRAME) built with
%----------------------

%% Step 1 Define nodes
WIREFRAME.Node=[ ...
  1 0 0 0 0.613 0 0      ;   2 0 0 0 1.133 -0.198 0
  3 0 0 0 1.548 -0.208 0 ;   4 0 0 0 1.538 -0.487 0
  5 0 0 0 1.133 -0.51 0  ;   6 0 0 0 0.613 -0.742 0
  7 0 0 0 0.32 0 0       ;   8 0 0 0 0.772 -0.287 0
  9 0 0 0 0.772 -0.521 0 ;  10 0 0 0 0.33 -0.7483 0
  11 0 0 0 0.41 -0.218 0 ;  12 0 0 0 0.38 -0.545 0
  13 0 0 0 0.1 -0.22066 0;  14 0 0 0 0.1 -0.55167 0
  15 0 0 0 0 -0.2177 0  ; 16 0 0 0 0 -0.544285 0

];

%% Step 2 Define connectivity
WIREFRAME.Elt=[];feplot(WIREFRAME)
fecom('TextNode');

% You can use the contextmenu (right click) cursor ... -> 3dline pick

feplot('addline',[4 3 2 1 7 15 16 10 6 5 4]) % add a first line
feplot('addline',[8 11 13 0 9 12 14 13 0 9 8 0 12 11]) % other line

%% Step 3 Load data for pole estimation
%-----------------

% download and load test data :
demosdt('download http://www.sdtools.com/contrib/kart_example.mat')

%% Step 4 : open interface and visualize test data

ci=idcom; % open interface and get pointers ci
iicom('curveinit','Test',TESTDATA);
 % .w Frequencies, .xf responses, .dof actuator/sensor info

iicom('SubMagPha');  % frequency response data plotted

%% Step 5 : Identify poles
% poles must now be identified one by one

idcom('e .01 44.6');
iicom('CurTabIdent')

% pole estimated at freq of 44.6Hz with possible range of .01 percent
% around that frequency.
% estimate checked on freq plot, if correct it is added

idcom('ea');

%% Step 6 : Add additional poles

idcom('e .01 48.8'); idcom('ea');
idcom('e .01 95.3'); idcom('ea');
idcom('e .01 125.3');idcom('ea');
idcom('e .01 141');  idcom('ea');

% five poles estimated in total and are stored in ci.Stack{'IdMain'}.po.
% modeshapes are based on narrowband estimate.

%% Step 7 : synthesize FRF based on estimated poles

idcom('est'); % build a broad-band identification

%% Step 8 : visualize mode shapes 
%---------------

cf=feplot;  cf.model=WIREFRAME; % plot WIREFRAME test model
cf.def=ci.Stack{'IdMain'}; fecom('view3') % display identified shapes

%% EndTuto

%% #Plate -------------------------------------------------------------------
elseif comstr(Cam,'plate');[CAM,Cam]=comstr(CAM,6);

if comstr(Cam,'prestress')
%% #PlatePresstress : sample plate with a pre-stressed beam

% A rectangular plate
model=feutil('objectquad 1 1 ',[0 0 0; 0 .400 0; 3 0 0],2,15);

% Add beams on center line and define pre-tension
elt=feutil('objectbeamline',feutil('findnode y==.2',model), ...
    [2 2 0 0 1]);  % MatId ProId vx,vy,vz orientation
model=feutil('addelt',model,'beam1t',elt(2:end,:));


E=2e10;
model.pl=[1 fe_mat('m_elastic','SI',1) E .2 2000 % Plate
          2 fe_mat('m_elastic','SI',1) 210e9 .3 7800]; % Beam
model.il=p_shell('dbval 1 Mindlin .1'); % Plate
model.il=p_beam(model.il,'dbval 2 I .03 .02 0.02  .005 .005 .005'); % I beam

model=feutil('setpro 2',model,'MAP', ...
   struct('dir',{{'1.5e6'}},'lab',{{'ten'}}));

feplot(model);fecom('curtabElProp','p_beam_2');fecom('view3');
fprintf('Click on Show3d or Showsection to view beams')
% a view to actually see something of the 3D section
iimouse('view',gca,[-7.93 -10.2 7.566 0.04831 0.2 -0.00040 0 0 1 1.03]);


% Define the boundary conditions
model=fe_case(model,'FixDof','left support','x==0 -DOF 1 2 3',...
    'FixDof','right support','x==3 -DOF 1 2 3');
 
def=fe_eig(model,[5 10 1e2 11]); 
 
% plotting 
cf=feplot(model);cf.def=def;
cf.sel(1)={'groupall','colordata z'};fecom view3


elseif comstr(Cam,'spring')
%% #PlateSpring : two plates connected by springs

% Two rectangular plates with duplicated nodes at the interface
model=feutil('objectquad 1 1 ',[0 0 0; 0 .400 0; 1 0 0],2,5);
model=feutil('objectquad 1 1 ',model,[1 0 0;0 .400 0; 1 0 0],2,5);
model=feutil('unjoin 1 2',model);

% Add springs, first ,odes at the connections
elt=[feutil('findnode group1 & x==1',model) ...
    feutil('findnode group2 & x==1',model)];

elt(:,3)=-123456; % Connect 6DOF
elt(:,5)=2; % ProID
model=feutil('addelt',model,'celas',elt);

% Define properties
model.pl=m_elastic('dbval 1 steel');
model.il=p_shell('dbval 1 Mindlin .02'); % Plate
model.il=p_spring(model.il,'dbval 2  1e6 '); % springs


% Define the boundary conditions
model=fe_case(model,'FixDof','left support','x==0 -DOF 1 2 3 4 5 6');

def=fe_eig(model,[5 20 0]);

cf=feplot(model,def);
cf.sel(1)={'groupall','colordata z'};fecom view3

elseif comstr(Cam,'maias')
%% #PlateMaias : bulk file of the MAIAS testbed

cd(sdtdef('tempdir')); 
if ~exist('Maias_Testbed.bdf','file')
 sdtcheck('PatchGet',struct('url','http://www.sdtools.com/contrib', ...
  'fname','Maias_Testbed.zip','back',1));
 unzip Maias_Testbed.zip
end
fname=which('Maias_Testbed.bdf');
model=nasread(fname);

%% #PlateEnd
else;error('Plate%s unknown',CAM);
end
%% #Beam -------------------------------------------------------------------
elseif comstr(Cam,'beam');[CAM,Cam]=comstr(CAM,5);

if comstr(Cam,'gravity')
%% #BeamGravity : check of stresses under gravity loading

model=struct('Node',[1 0 0 0  0 0 0;2 0 0 0  0 0 1; 3 0 0 0 1 0 1;
           4 0 0 0 2 0 1; 5 0 0 0 2 0 0], ...     %;3 0 0 0  1 0 1;4 0 0 0  1 0 0
     'Elt',feutil('addelt','beam1t', ...
          [1:2  1 1;2:3  1 1;3:4  1 1;4:5  1 1;1,3 1 1;3,5 1 1]), ...    
     'pl',[1 fe_mat('m_elastic','SI',1)   210e9   0.3   7850], ...   
     'il',[1 fe_mat('p_beam','MM',3)  0  comstr('I',-32') 127 76 76 4 7.6 7.6]); 
     %UB 127x76x13 PFC -b-h-tw-tf
     
model.il=fe_mat('convertSI',model.il); % make units SI
il=p_beam('convertTo1',model.il);% see 
         % p_beam('show3D',cf)
         % cf=feplot(model);beam1t('map',cf.mdl);
% boundary conditions

topload=struct('def',ones(3,1)*-1e7,'DOF',(2:4)'+.03);
model=fe_case(model,'fixdof','Clamped',[1.01 1.02 1.03 5.01 5.02 5.03]', ...
    'grav','Gravity',struct('dir',9.81*[0;0;-1]), ...
    'DofLoad','TopLoad',topload);
out=model;

%% #BeamEnd
else;error('Beam%s unknown',CAM);
end


%% clean end
elseif comstr(Cam,'@');out=eval(CAM);
%% #TutoDo: recover model from a specific tuto step -3
elseif comstr(Cam,'tuto'); 
 eval(sdtweb('_tuto',struct('file','d_primer','CAM',CAM)));
 if nargout==0; clear out; end

elseif comstr(Cam,'cvs')
 out=sdtcheck('revision','$Revision: 6e8b3bf $  $Date: 2021-03-16 18:06:23 +0100 $ ');
else; error('%s unknown',CAM);
end 
%% #End function
end
