function [out,out1,out2]=d_mesh(varargin); %#ok<*NOSEM,*STOUT>

% D_MESH sample meshes used for support
%
% Use 
%     sdtweb('_taglist','d_mesh') to view current contents
%     d_mesh('tuto')  % to see integrated tutorials

%       Copyright (c) 1990-2021 by SDTools, All Rights Reserved.
%       For revision information use d_mesh('cvs')

if nargin==0; d_mesh('tuto'); return; end

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
if comstr(Cam,'truss');[CAM,Cam]=comstr(CAM,6);

if comstr(Cam,'bmesh')
%% #TutoBmesh direct FEM declaration - - - - - - - - - - - - - - - - - -2

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

elseif comstr(Cam,'beambar')
%% #TutoBeamBar : illustrate mixed beam/bar using FEMESH - - - - - - - - - -2

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

elseif comstr(Cam,'fesh')
%% #TutoFesh building models with femesh - - - - - - - - - - - - - - - - - -2

%% Step1 : Initialize model
model=struct('Node',[1 0 0 0  0   0   0;     2 0 0 0 0   0   0.15;
                     3 0 0 0  0.4 1.0 0.176; 4 0 0 0 0.4 0.9 0.176],...
             'Elt',[],'unit','SI','name','GARTEUR');
%% Step2 Fuselage
model.Elt=feutil('ObjectBeamLine 1 2',model);
model=feutil('Extrude 0  1.0 0.0 0.0',model,...
             [linspace(0,.55,5) linspace(.65,1.4,6) 1.5]);
%% Step3 vertical tail
n1=feutil('FindNode z==.15 & x>=1.4',model);
mo0=model; mo0.Elt=feutil('ObjectBeamLine',n1);
mo0=feutil('Extrude 3 0 0 .1',mo0);
model=feutil('AddTestCombine-noori',model,mo0);
%% Step4 Vertical horizontal tail
n1=feutil('FindNode z==.45',model)
mo0=model; mo0.Elt=feutil('ObjectBeamLine',n1);
mo0=feutil('Extrude 0  0.0 0.2 0.0',mo0,[-1 -.5 0 .5 1]);
model=feutil('AddTestCombine;-noori',model,mo0);

%% right drum
mo0=model; mo0.Elt=feutil('ObjectBeamLine 3 4');
mo0=feutil('Extrude 1 .4 0 0',mo0);
mo0=feutil('Divide',mo0,[0 2/40 15/40 25/40 1],[0 .7 1]);
model=feutil('AddTestCombine;-noori',model,mo0);

%% left drum
mo0=feutil('SymSel 1 0 1 0',mo0);
model=feutil('AddTestCombine;-noori',model,mo0);

%% wing
n1=feutil('FindNode y==1 & x>=.55 & x<=.65',model);
mo0=model; mo0.Elt=feutil('ObjectBeamLine',n1);
mo0=feutil('Divide',mo0,[0 1-.762 1]);
mo0=feutil('Extrude 0  0.0 -1.0 0.0',mo0,[0 0.1 linspace(.15,.965,9) ...
                                          linspace(1.035,1.85,9) 1.9 2.0]);
model=feutil('AddTestCombine;-noori',model,mo0);

%% Connection plate
n1=feutil('FindNode y==0.035 | y==-0.035 & x==.55',model)
mo0=model; mo0.Elt=feutil('ObjectBeamLine',n1);
mo0=feutil('Divide 2',mo0);
mo0=feutil('TransSel -.02 0 0',mo0);
mo0=feutil('Extrude 0 1 0 0',mo0,[0 .02 .12 .14]);
i1=intersect(feutil('FindNode group6',model),feutil('FindNode group1',mo0));
mo0=feutil('TransSel 0.0 0.0 -0.026',mo0);
model=feutil('AddTestCombine;-noori',model,mo0);
%% Step5 Stiff links for the connection
mo0=model; mo0.Elt=feutil('Object mass',i1);
mo0=feutil('Extrude 1 0 0 -.026',mo0);
mo0.Elt=feutil('set group1 name celas',mo0);
%% Step6 set connected DOFs and spring value
mo0.Elt(2:end,3)=12345; % master dof
mo0.Elt(2:end,4)=0; % same dof as master
mo0.Elt(2:end,7)=1e12; % stiffness
model=feutil('AddTestCombine;-noori',model,mo0); % add springs to main model
%% Step7 Make a group of the part covered by the constraining layer
model.Elt=feutil('Divide group 6 InNode {x>.55 & y<=.85 & y>=-.85}',model);
%% Step8 Tip masses
i1=feutil('FindNode y==0.93 | y==-0.93 & x==0.42',model)
mo0=model; mo0.Elt=feutil('Object mass',i1,[0.2 0.2 0.2]); %200g
model=feutil('AddTestCombine;-noori',model,mo0);
i1=feutil('FindNode z==.45 & y==0',model)
mo0=model; mo0.Elt=feutil('Object mass',i1,[0.5 0.5 0.5]); %500g
model=feutil('AddTestCombine;-noori',model,mo0);
model=feutil('Join mass1',model); % all mass in the same group
%% Step9 Orient plates that will need an off-set
model.Elt=feutil('Orient 4:8 n 0 0 3',model);
i1=feutil('FindElt group4:5',model);
model.Elt(i1,9)=0.005; % drums (positive off-set)
i1=feutil('FindElt group6:7',model);
model.Elt(i1,9)=-0.005; % wing
i1=feutil('FindElt group8',model);
model.Elt(i1,9)=0.008; % wing
%% Step10 Deal with material and element properties identifier:
model.Elt=feutil('Set group1 mat1 pro3',model);
model.Elt=feutil('Set group2:7 mat1 pro1',model);
model.Elt=feutil('Set group8 mat2 pro2',model);
model.Elt=feutil('Set group6 pro4',model);
%% Step11 Define associated properties:
model.pl=[m_elastic('dbval 1 aluminum');
          m_elastic('dbval 2 steel')];
model.il = [1 fe_mat('p_shell','SI',1)  2 1 0     .01
            2 fe_mat('p_shell','SI',1)  2 1 0     .016
            3 fe_mat('p_shell','SI',1)  2 1 0     .05
            4 fe_mat('p_shell','SI',1)  2 1 0     .011];
%% Step12 Display in feplot
 cf=comgui('guifeplot -project "SDT Root"',3); % Robust open in figure(3)
 cf.model=model;       % display model
 fecom(';sub 1 1;view3; colordatamat-edgealpha.1'); % 1 subplot, specify view, color,

%EndTuto

else;error('Truss%s unknown',CAM);
end

%% #Plate -------------------------------------------------------------------
elseif comstr(Cam,'plate');[CAM,Cam]=comstr(CAM,6);

if comstr(Cam,'feut')
%% #TutoFeut direct FEM declaration - - - - - - - - - - - - - - - - - -2

%% Step1 : Initialize model
model=struct('Node',[1 0 0 0  0   0   0;     2 0 0 0 0   0   0.15;
                     3 0 0 0  0.4 1.0 0.176; 4 0 0 0 0.4 0.9 0.176],...
             'Elt',[],'unit','SI','name','GARTEUR');
%% Step2 Fuselage
model.Elt=feutil('ObjectBeamLine 1 2',model);
model=feutil('Extrude 0  1.0 0.0 0.0',model,...
             [linspace(0,.55,5) linspace(.65,1.4,6) 1.5]);
%% Step3 vertical tail
n1=feutil('FindNode z==.15 & x>=1.4',model);
mo0=model; mo0.Elt=feutil('ObjectBeamLine',n1);
mo0=feutil('Extrude 3 0 0 .1',mo0);
model=feutil('AddTestCombine-noori',model,mo0);
%% Step4 Vertical horizontal tail
n1=feutil('FindNode z==.45',model)
mo0=model; mo0.Elt=feutil('ObjectBeamLine',n1);
mo0=feutil('Extrude 0  0.0 0.2 0.0',mo0,[-1 -.5 0 .5 1]);
model=feutil('AddTestCombine;-noori',model,mo0);

%% right drum
mo0=model; mo0.Elt=feutil('ObjectBeamLine 3 4');
mo0=feutil('Extrude 1 .4 0 0',mo0);
mo0=feutil('Divide',mo0,[0 2/40 15/40 25/40 1],[0 .7 1]);
model=feutil('AddTestCombine;-noori',model,mo0);

%% left drum
mo0=feutil('SymSel 1 0 1 0',mo0);
model=feutil('AddTestCombine;-noori',model,mo0);

%% wing
n1=feutil('FindNode y==1 & x>=.55 & x<=.65',model);
mo0=model; mo0.Elt=feutil('ObjectBeamLine',n1);
mo0=feutil('Divide',mo0,[0 1-.762 1]);
mo0=feutil('Extrude 0  0.0 -1.0 0.0',mo0,[0 0.1 linspace(.15,.965,9) ...
                                          linspace(1.035,1.85,9) 1.9 2.0]);
model=feutil('AddTestCombine;-noori',model,mo0);

%% Connection plate
n1=feutil('FindNode y==0.035 | y==-0.035 & x==.55',model)
mo0=model; mo0.Elt=feutil('ObjectBeamLine',n1);
mo0=feutil('Divide 2',mo0);
mo0=feutil('TransSel -.02 0 0',mo0);
mo0=feutil('Extrude 0 1 0 0',mo0,[0 .02 .12 .14]);
i1=intersect(feutil('FindNode group6',model),feutil('FindNode group1',mo0));
mo0=feutil('TransSel 0.0 0.0 -0.026',mo0);
model=feutil('AddTestCombine;-noori',model,mo0);
%% Step5 Stiff links for the connection
mo0=model; mo0.Elt=feutil('Object mass',i1);
mo0=feutil('Extrude 1 0 0 -.026',mo0);
mo0.Elt=feutil('set group1 name celas',mo0);
%% Step6 set connected DOFs and spring value
mo0.Elt(2:end,3)=12345; % master dof
mo0.Elt(2:end,4)=0; % same dof as master
mo0.Elt(2:end,7)=1e12; % stiffness
model=feutil('AddTestCombine;-noori',model,mo0); % add springs to main model
%% Step7 Make a group of the part covered by the constraining layer
model.Elt=feutil('Divide group 6 InNode {x>.55 & y<=.85 & y>=-.85}',model);
%% Step8 Tip masses
i1=feutil('FindNode y==0.93 | y==-0.93 & x==0.42',model)
mo0=model; mo0.Elt=feutil('Object mass',i1,[0.2 0.2 0.2]); %200g
model=feutil('AddTestCombine;-noori',model,mo0);
i1=feutil('FindNode z==.45 & y==0',model)
mo0=model; mo0.Elt=feutil('Object mass',i1,[0.5 0.5 0.5]); %500g
model=feutil('AddTestCombine;-noori',model,mo0);
model=feutil('Join mass1',model); % all mass in the same group
%% Step9 Orient plates that will need an off-set
model.Elt=feutil('Orient 4:8 n 0 0 3',model);
i1=feutil('FindElt group4:5',model);
model.Elt(i1,9)=0.005; % drums (positive off-set)
i1=feutil('FindElt group6:7',model);
model.Elt(i1,9)=-0.005; % wing
i1=feutil('FindElt group8',model);
model.Elt(i1,9)=0.008; % wing
%% Step10 Deal with material and element properties identifier:
model.Elt=feutil('Set group1 mat1 pro3',model);
model.Elt=feutil('Set group2:7 mat1 pro1',model);
model.Elt=feutil('Set group8 mat2 pro2',model);
model.Elt=feutil('Set group6 pro4',model);
%% Step11 Define associated properties:
model.pl=[m_elastic('dbval 1 aluminum');
          m_elastic('dbval 2 steel')];
model.il = [1 fe_mat('p_shell','SI',1)  2 1 0     .01
            2 fe_mat('p_shell','SI',1)  2 1 0     .016
            3 fe_mat('p_shell','SI',1)  2 1 0     .05
            4 fe_mat('p_shell','SI',1)  2 1 0     .011];
%% Step12 Display in feplot
 cf=comgui('guifeplot -project "SDT Root"',3); % Robust open in figure(3)
 cf.model=model;       % display model
 fecom(';sub 1 1;view3; colordatamat-edgealpha.1'); % 1 subplot, specify view, color,

% EndTuto 


elseif comstr(Cam,'plate')
%% #Tutoplate : old d_plate example  - - - - - - - - - - -2

% Illustration of FEMESH capabilities and standard finite element computations
% for a few plate/shell examples
%
% Use sdtweb('fem') to access the introduction to FEMESH. The associate
% basic example is treated in the D_TRUSS demo.

% For a first case, one will build a square plate by defining two nodes, a beam
% connecting these nodes, extruding it to form a quadrilateral and dividing the
% element in a regular 4 by 4 mesh

%% Step 1 : Create mesh 
node=[0 0 0;0 1 0;1 1 0;1 0 0];
model=feutil('ObjectQuad 10 11',node,4,4);

%% Step 2 Display the model, the model is accessible by cf.mdl
cf=feplot(model); fecom('view3');

%% Step 3 : Material properties
% To compute the modes of this model one defines material and section properties

% data structure format of FE model : see doc fem
cf.mdl.pl=m_elastic('dbval 10 Aluminum');
cf.mdl.il=p_shell('dbval 11 kirchoff 5e-2');

%% Step 4 : Boundary conditions
% defines boundary conditions here clamped on left edge (x==0)
% see sdtweb('fe_case#FixDof')
cf.mdl=fe_case(cf.mdl,'fixdof','Clamped edge','x==0');

%% Step 5 : Compute and display mode shapes
% Set eigen value solver options, method 6, 20 modes, shift 1e3 
cf.mdl=stack_set(cf.mdl,'info','EigOpt',[6 20 1e3]);
% data are accessible by cf.Stack

% computes and displays modes
cf.def=fe_eig(cf.mdl);

%% Step 6 : Plot strain energy per element
% Display strain energy at elements, see sdtweb fe_stress#feplot
Ek=fe_stress('ener -MatDes 1 -curve',cf.mdl,cf.def);
feplot('ColorDataElt',Ek); fecom('view3')

%% Step 7 Quarter plate with distributed static load
% We will now consider a typical problem used to evaluate the quality of
% finite element model predictions.

%   Build model

cf.mdl=fe_case(cf.mdl,'reset',... % 1st reset the existing case 
 'FixDOF','simple support','y==1 | x==1 -DOF 3',...
 'FixDOF','Bending',[.01;.02;.06], ... % fix all x y and thetaz DOF
 'FixDOF','symmetry_y','x==0 -DOF 5',... % fix thetay for all nodes located at x=0
 'FixDOF','symmetry_x','y==0 -DOF 4' ); % fix thetax for all nodes located at y=0

%% Step 8  Now define a uniform volume load
data=struct('sel','groupall','dir',[0 0 9.81]);
cf.mdl=fe_case(cf.mdl,'FVol','Gravity',data);

fe_case(cf.mdl,'info') 
% load/case data are accessible by cf.CStack. On can also visualises Case
% entries in feplot : fecom(cf,'CurTabCases'); fecom promodelviewon

% One can build the load at nodes
Load=fe_load(cf.mdl); % Build the defined load

%% Step 9 : The response is then easily computed and displayed using
cf.def=fe_simul('Static',cf.mdl);
fecom(';sub 1 1 1 3 2;colordataevala');

%%  3D model of a short tube 
%Step 10 - meshing
model=struct('Node',[1  0 0 0  0.0 0.0 0.0;
                     2  0 0 0  0.0 1.0 0.0],'Elt',[]);
model.Elt=feutil('ObjectBeamLine 1 2',model);
model=feutil('rev 10 o 0 0 .1 360 0 1 0',model);
model=feutil('Divide 1 3',model);
cf=feplot(model); fecom('showline'); % display in feplot

model=fe_case(model,'reset'); % no boundary conditions

%% Step 11 Define material and section properties
model.pl=m_elastic('dbval 1 Aluminum'); % (Matid 1) material 
model.il=p_shell('dbval 1 kirchoff 5e-2'); % (ProId 1) Kirchoff plate, thickness 5e-2
model.unit='SI'; % mesh unit
model.name='tube'; % model name

%% Step 12 Mode shape computation and visualisation 
%Set eigen value solver options, method 5, 20 modes, shift 1e3 
model=stack_set(model,'info','EigOpt',[5 20 1e3]);

% computes  modes
def1=fe_eig(model);

% displays model and modes
cf=feplot(model,def1);
fecom(';ch7;scd.05');

% Here one computes the modes of a structure showing a symmetry along a plan
% The example is axisymmetric, but only the symmetry along plan xz is
% looked at
% The structure can be simplified by applying Symmetric and Antisymmetric
% boundary conditions to the pipe cut in xz at half length

%% Step 13 : Axisymmetric model -  build model
model=struct('Node',[1  0 0 0  0.0 0.0 0.0;
                     2  0 0 0  0.0 2.0 0.0],'Elt',[]);
model.Elt=feutil('ObjectBeamLine 1 2',model);
model=feutil('rev 10 o 0 0 .1 360 0 1 0',model);
model=feutil('Divide 2 40',model);

model=fe_case(model,'reset'); % no boundary conditions

%% Step 14 Define material and section properties
model.pl=m_elastic('dbval 1 Aluminum');
model.il=p_shell('dbval 1 kirchoff 5e-2');
model.unit='SI'; % mesh unit
model.name='tube'; % model name

% Set eigen value solver options, method 5, 25 modes, shift 1e3 
model=stack_set(model,'info','EigOpt',[5 25 1e3]);

%% Step 15 : Compute and display modes with specific symmetries
% Set half length pipe
mo1=model;
mo1.Elt=feutil('SelElt InNode{y<=1}',model);
 % Compute modes with the symmetry condition
 mo1=fe_case(mo1,'fixdof','Symm','y==1 -DOF 2 4 6'); % y, \theta_x, \theta_z
 defS=fe_eig(mo1);
 % Compute modes with the antisymmetry condition
 mo1=fe_case(mo1,'reset');
 mo1=fe_case(mo1,'fixdof','AntiSymm','y==1 -DOF 1 3 5'); % x, z, \theta_y
 defA=fe_eig(mo1);
 % Concatenation of all deformations
 r1=[defS.data ;defA.data];
 [r1,i1]=sort(r1,'ascend');
 defAS=defS;
 defAS.def=[defS.def defA.def];
 defAS.def=defAS.def(:,i1); defAS.data=r1;
  
% Modes of the full pipe
model=stack_set(model,'info','EigOpt',[5 50 1e3]);
def0=fe_eig(model);
[i1,r1]=feutil('findnode y>=1',model);
def1=feutilb('placeindof',defS.DOF,def0);

%% Step 16 : compare modes with full pipe
% The results obtained for the half length pipe are identical to the full
% pipe modes
figure(4); ii_mac(def1,defAS,'macploterror')
figure(5); semilogy(100*abs([def1.data(7:end)./defAS.data(7:end)-1]));

% EndTuto 
 
elseif comstr(Cam,'prestress')
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

elseif comstr(Cam,'table')
%% #PlateTable : Vibrating table
if 1==2
 fname='O:\sdtdata\ref_test\ansys\cdb\table.cdb';
 model=ans2sdt('read',fname);
 pl=feutil('getmat',model,101); pl(1)=1;
 il=feutil('getpro',model,1);
end

% Parameters : 
% Table Geometry :
l1=1000; % table half width
h1=400; % table height
angle=20; % angle to fix the actuators
% Satelite geometry
b1=750; %  Base radius bottom
b2=600; % Base radiu top
b3=300; % Base height
cw=1800; % Cube width
ch=2500; % Cube height
% Force Cell's height:
f1=100; % Height of cell 
fn=12; % Number of cells
% Meshinf refinement
lref=2; % Mesh refinment in plane
href=3; % Mesh refinment in height
rref=2; % Sat base refinement

% Material parameters
RO.pl=m_elastic('dbval 1 -unit TM Aluminum', ...% Aluminum plates
     'dbval 3 -unit TM steel', ...% Steel plates for satelite's base
     'dbval 4 -unit TM steel'); % Steel Beams  satelite
RO.il=p_shell('dbval 1 -punit TM Mindlin 30', ... % Thickness 30mm
      'dbval 3 -punit TM Mindlin 20'); % Base thickness 20mm
RO.il=p_beam(RO.il,'dbval 4 -punit TM box 100 100 3 3'); % Rectangular sections 200mm


il2=[2 fe_mat('p_spring;','TM',2) repmat(1e2,1,6)]; % Actuator cbush
il5=[5 fe_mat('p_spring;','TM',2) repmat(1e10,1,6)]; % force cell cbush
proil={'pro','Act',struct('il',il2,'type','p_spring','NLdata', ...
      struct('type','nl_inout','Fu',@abs,'isens',1));
       'pro','FCell',struct('il',il5,'type','p_spring')};
% Intermediate variables
gap=2*l1/6; % gap between each stiffener

%% Table
mo1=feutil('objectquad 1 1 ',[-l1 -l1 0;l1 -l1 0;l1 l1 0;-l1 l1 0],6*lref,6*lref);
mo2=mo1;mo2.Elt=feutil('GetEdgeLine',mo1); % Get edge for further extrusion
mo1=feutil(sprintf('RepeatSel 2 0 0 -%.16f',h1),mo1);
% First direction extrusion + rep
mo3=mo2;
mo3.Elt=feutil(sprintf('selelt innode{x==-%.16f}',l1),mo3);
i1=h1/href; % Height of a single mesh
mo3=feutil(sprintf('Extrude %.16f 0 0 %.16f',href,-i1),mo3);
mo3=feutil(sprintf('RepeatSel 7 %.16f 0 0',gap),mo3);
mo1=feutil('AddtestMerge-noori;',mo1,mo3);
% First direction extrusion + rep
mo3=mo2;
mo3.Elt=feutil(sprintf('selelt innode{y==-%.16f}',l1),mo3);
mo3=feutil(sprintf('Extrude %.16f 0 0 %.16f',href,-i1),mo3);
mo3=feutil(sprintf('RepeatSel 7 0 %.16f 0',gap),mo3);
mo1=feutil('AddtestMerge-noori;',mo1,mo3);

%% Accroches verins
%xxx : merge fixdofs and connection surfaces
mo2=mo1;
mo2.Elt=feutil(sprintf('selelt innode{x==-%.16f & y<=%.16f & y>=0 & z>=-%.16f}',...
 l1,l1-gap,h1/2),mo2);
mo3=feutil(sprintf('RotateNode o -%.16f %.16f 0 -%.16f 0 0 1',...
 l1,l1-gap,angle),mo2);
mo2.Elt=feutil('GetEdgeLine',mo2);
mo2=feutil(sprintf('Rev-OptimDegen 1 o -%.16f %.16f 0 -%.16f 0 0 1',...
 l1,l1-gap,angle),mo2);
mo2=feutil('AddtestMerge-noori;',mo2,mo3);

% Get the 4 lists of two nodes by symmetry
mo2=feutil('AddtestMerge-noori;',mo2,feutil('node mir "eq -1 1 0 0";',mo2));
mo2=feutil('AddtestMerge-noori;',mo2,feutil('node mir "eq 1 1 0 0";',mo2));
mo1=feutil('AddtestMerge-noori;',mo1,mo2);
mo1.Elt=feutil('orient',mo1);
%% Horizontal actuators
[EltId,mo1.Elt]=feutil('EltIdFix;',mo1.Elt);
mo2=mo1;
% Get face to build the two points
mo2.Elt=feutil(sprintf('selelt withnode{y==0} & innode{x<=-%.16f}',l1),mo2);
EltId=feutil('EltId',mo2);
r1=feutil('getnormal',mo2);
EltId=EltId(r1(:,2)>0.8);
mo2.Elt=feutil('selelt eltid',mo2,EltId);
mo2.Node=feutil('optimmodel',mo2);
n1=mean(mo2.Node(:,5:7)); 
dir=(n1-[-l1 l1-gap n1(3)]);dir=dir/norm(dir);
n2=n1+dir*l1;
% Get the 4 lists of two nodes by symmetry
n1n2=[[n1;n2];feutil('node mir "eq -1 1 0 0";',[n1;n2])];
n1n2=[n1n2;feutil('node mir "eq 1 1 0 0";',n1n2)];
n1n2=mat2cell(n1n2,2*ones(1,4),3);
for j1=1:length(n1n2)
 [mo1.Node,i1]=feutil('AddNode',mo1.Node,n1n2{j1});
 nid=mo1.Node(i1,1);
 % Get face to build the two points
 mo1=fe_caseg('Connection surface 123456 -MaxDist1;',mo1,sprintf('surf_act%i',j1), ...
    sprintf('NodeId %i',nid(1)), ...    % Selection of nodes to connect
    sprintf('withnode{x<-%.16f | x>%.16f | y<-%.16f | y>%.16f}',l1,l1,l1,l1));
 mo1=fe_case(mo1,'FixDof',sprintf('fix_act%i',j1),sprintf('NodeId %i',nid(2)));
 mo1.Elt=feutil('addelt',mo1.Elt,'cbush',[nid' 2 2 0 0 0 0 -1]); % EDID=-1 local coordinates
end

%% Vertical actuators
n1n2=[-2*l1/3 -2*l1/3 -h1;
      -2*l1/3 -2*l1/3 -h1-l1;
       2*l1/3 -2*l1/3 -h1;
       2*l1/3 -2*l1/3 -h1-l1;
      -2*l1/3  2*l1/3 -h1;
      -2*l1/3  2*l1/3 -h1-l1;
       2*l1/3  2*l1/3 -h1;
       2*l1/3  2*l1/3 -h1-l1];
n1n2=mat2cell(n1n2,2*ones(1,4),3);
for j1=1:length(n1n2)
 [mo1.Node,i1]=feutil('AddNode',mo1.Node,n1n2{j1});
 nid=mo1.Node(i1,1);
 % Get face to build the two points
 mo1=fe_caseg('Connection surface 123456 -MaxDist1;',mo1,sprintf('surf_act%i',j1+4), ...
    sprintf('NodeId %i',nid(1)), ...    % Selection of nodes to connect
    sprintf('innode{z==-%.16f}',h1));
 mo1=fe_case(mo1,'FixDof',sprintf('fix_act%i',j1+4),sprintf('NodeId %i',nid(2)));
 mo1.Elt=feutil('addelt',mo1.Elt,'cbush',[nid' 2 2 0 0 0 0 -1]); % EDID=-1 local coordinates
end

%% Base+top_plate
mo2=struct('Node',[],'Elt',[]);
[mo2.Node,i1]=feutil('addnode',mo2.Node,[0 0 b3;cw/2 0 b3;0 cw/2 b3]);
r1=[1 2 3 b2 b2 rref rref 4];
mo2=feutil(...
 sprintf('ObjectHoleInPlate %i %i %i %.16f %.16f %i %i %i',r1),mo2);

mo3=struct('Node',[],'Elt',[]);
[mo3.Node,i1]=feutil('addnode',mo3.Node,[b2 0 b3;b2 0 0;b1 0 0]);
mo3.Elt=feutil('addelt',mo3.Elt,'beam1',i1([1 2;2 3]));
mo3=feutil(sprintf('Rev %i o 0 0 0 360 0 0 1',8*rref),mo3);
mo2=feutil('AddTestMerge-noori;',mo2,mo3);

mo3=mo2;
mo3.Elt=feutil(sprintf('selelt innode{z==%.16f}',b3),mo3);
mo3.Elt=feutil('GetEdgeLine',mo3);
mo3.Elt=feutil(sprintf('selelt innode{x==%.16f}',-cw/2),mo3);
mo3=feutil(sprintf('Extrude %i %.16f 0 0',rref*2,cw/2/rref),mo3);
mo3=feutil('node',mo3,struct('trans',[0 0 ch]));
mo2=feutil('AddTestMerge-noori;',mo2,mo3);

mo2.Elt=feutil('SetGroupall Matid3 Proid3',mo2);

%%  Satelite
mo3=mo2;
mo3.Elt=feutil(sprintf('selelt innode{z==%.16f}',b3),mo3);
mo3.Elt=feutil('GetEdgeLine',mo3);
mo3.Elt=feutil(sprintf('selelt withnode{cyl >=%.16f o 0 0 0 0 0 1}',b2+1),mo3);
mo3=feutil(sprintf('RepeatSel 3 0 0 %.16f',ch/2),mo3);
mo4=struct('Node',[],'Elt',[]);
mo4.Node=feutil('addnode',mo4.Node,[cw/2 -cw/2 b3;cw/2 cw/2 b3]);
mo4.Elt=feutil('AddElt','mass1',[1;2]);
mo4=feutil(sprintf('RepeatSel 3 %.16f 0 0',-cw/2),mo4);
mo4=feutil(sprintf('Extrude %i 0 0 %.16f',rref*2,ch/2/rref),mo4);
mo3=feutil('AddTestMerge-noori;',mo3,mo4);

mo4=struct('Node',[],'Elt',[]);
mo4.Node=feutil('addnode',mo4.Node,...
 [-cw/2 -cw/2 b3;0 -cw/2 b3+ch/2;cw/2 -cw/2 b3]);
mo4.Elt=feutil('AddElt','beam1',[1 2;2 3]);
mo4=feutil(sprintf('Divide %i',rref),mo4);
mo4=feutil(sprintf('RepeatSel 2 0 0 %.16f',ch/2),mo4);
mo5=feutil('RotateNode o 0 0 0 180 0 0 1',mo4);
mo4=feutil('AddTestMerge-noori;',mo4,mo5);
mo3=feutil('AddTestMerge-noori;',mo3,mo4);

mo4=struct('Node',[],'Elt',[]);
mo4.Node=feutil('addnode',mo4.Node,...
 [cw/2 -cw/2 b3;cw/2 cw/2 b3+ch/2]);
mo4.Elt=feutil('AddElt','beam1',[1 2]);
mo4=feutil(sprintf('Divide %i',2*rref),mo4);
mo4=feutil(sprintf('RepeatSel 2 0 0 %.16f',ch/2),mo4);
mo5=feutil('RotateNode o 0 0 0 180 0 0 1',mo4);
mo4=feutil('AddTestMerge-noori;',mo4,mo5);
mo3=feutil('AddTestMerge-noori;',mo3,mo4);

mo3.Elt=feutil('SetGroupall Matid4 Proid4',mo3);

%% Satelite + Base + height shift
mo2=feutil('AddTestMerge-noori;',mo2,mo3);
mo2=feutil('node',mo2,struct('trans',[0 0 f1]));
mo1=feutilb('combinemodel CompatMatPro',mo1,mo2);

%% Force cells (nodes + connection surfaces + cbush)
% Nodes on the table surface
mo2=feutil(sprintf('object circle 0 0 0 %.16f 0 0 1 %i',(b1+b2)/2,fn));
[mo1.Node,i1]=feutil('AddNode',mo1.Node,mo2.Node(:,5:7)); 
nid1=mo1.Node(i1,1);
% Node on the base surface
mo2=feutil('node',mo2,struct('trans',[0 0 f1]));
[mo1.Node,i2]=feutil('AddNode',mo1.Node,mo2.Node(:,5:7)); 
nid2=mo1.Node(i2,1);
for j1=1:length(nid1)
 mo1=fe_caseg('Connection surface 123456 -MaxDist1',mo1,sprintf('force_cell_table%i',j1), ...
    sprintf('NodeId %i',nid1(j1)), ...    % Selection of nodes to connect
    'innode{z==0}');
 mo1=fe_caseg('Connection surface 123456 -MaxDist1',mo1,sprintf('force_cell_base%i',j1), ...
    sprintf('NodeId %i',nid2(j1)), ...    % Selection of nodes to connect
    sprintf('innode{z==%.16f}',f1));
 mo1.Elt=feutil('addelt',mo1.Elt,'cbush',[nid1(j1) nid2(j1) 5 5 0 0 0 0 -1]); % EDID=-1 local coordinates
end

%% Mat/pro prop
mo1.pl=RO.pl;mo1.il=RO.il; mo1=stack_set(mo1,proil);

%% Clean model
mo1=feutil('joinall',mo1); 
mo1=fe_case('merge',mo1,'#force_cell.*','force_cell');
mo1=fe_case('merge',mo1,'#surf_act.*','surf_act');
mo1=fe_case('merge',mo1,'#fix_act.*','fix_act');
out=mo1;

elseif comstr(Cam,'gen'); [CAM,Cam]=comstr(CAM,4);
%% #PlateGen_eric : variants of plates for parametric studies-------------2
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
 '-BC("free"#%s#"Clamping option (default : free plate)")' ...
 '-div(10#%i#"Division in X and Y")' ...
 '-PatchWidth(0.5#%g#"Relative patch width")' ...
 '-Sens("patch"#%s#"Small patch of sensors in the middle of the plate")' ...
 '-Load("botleft"#%s#Unit input at bottom left of the plate")'
 ],{RO,CAM}); Cam=lower(CAM);
if comstr(Cam,'mesh'); [CAM,Cam]=comstr(CAM,5);
 %% #PlateGenMesh : Mesh with varying number of elts + clamp strat---------3
 node=[0 0 0;0 1 0;1 1 0;1 0 0];
 model=feutil('ObjectQuad 1 1',node,RO.div,RO.div);
 
 dof=[]; % List of fixed dofs
 if strcmpi(RO.BC,'free'); dof=[]; % Free plate
 else; error('BC %s unknown');
 end
 %Fix dofs
 if ~isempty(dof)
  model=fe_case(model,'fixdof','fixed',dof);
 end
 % Default mat + damp
 model.pl=m_elastic('dbval 1 Aluminum');
 model.il=p_shell('dbval 1 kirchoff 5e-2');

 % Add loss factor (hysteretic damping)
 model=feutil('setmat 1 eta=0.02',model);

elseif comstr(Cam,'sens'); [CAM,Cam]=comstr(CAM,5);
 %% #PlateGenSens : Place sensors (right side, ...)------------------------3
 if ~isfield(RO,'model'); error('Provide RO.model'); end
 model=RO.model;
 if strcmpi(RO.Sens,'patch'); % patch of sensors in the middle of the plate
  i1=[.5-RO.PatchWidth/2 .5+RO.PatchWidth/2];
  n1=feutil('findnode',model,sprintf('x>=%g&x<=%g&y>=%g&y<=%g',i1,i1));
  tdof=feutil('getdof',n1,.03); % Sensors only in z direction
 else; error('unknown Sens command');
 end
 % Place sensors over the model
 model=fe_case(model,'SensDof append trans','Test',tdof);
elseif comstr(Cam,'load'); [CAM,Cam]=comstr(CAM,5);
 %% #PlateGenLoad : Place Load (first quarter of beam,...)-----------------3
 if ~isfield(RO,'model'); error('Provide RO.model'); end
 model=RO.model;
 if strcmpi(RO.Load,'botleft');
  % Add Load at the middle left of the beam
  i1=feutilb('AddNode -nearest',[0.1 0.1 0],model.Node(:,5:7));
  i1=model.Node(find(i1),1);
  model=fe_case(model,'DofLoad','Load',i1+.03);
 else; error('unknown Load command');
 end
else; error('Command d_mesh BeamGen%s unknown',CAM);
end
out=model;
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
%% #BeamNL ->sdtweb d_fetime testBeamNL
elseif comstr(Cam,'gen'); [CAM,Cam]=comstr(CAM,4);
%% #BeamGen_eric : variants of single beam for parametric studies---------2
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
 '-BC("clampl"#%s#"Clamping option (default : clamp beam left side)")' ...
 '-div(20#%i#"Number of elements")' ...
 '-Sens("right"#%s#"Where to place sensors")' ...
 '-Load("midleft"#%s#"Unit input at beam middle left")' ...
 '-eta(.02#%g#"Default loss factor")'
 ],{RO,CAM}); Cam=lower(CAM);
if comstr(Cam,'mesh'); [CAM,Cam]=comstr(CAM,5);
 %% #BeamGenMesh : Mesh with varying number of elts + clamp strat---------3
 if RO.div<5; sdtw('_nb','At least 5 elements for the beam'); RO.div=5; end
 model=femesh(sprintf('testbeam1 -divide %i',RO.div));
 dof=[]; % List of fixed dofs
 if strcmpi(RO.BC,'clampl'); dof=[dof;1]; % fix left side of beam
 elseif strcmpi(RO.BC,'clampr'); dof=[dof;2]; % fix right side of beam
 elseif strcmpi(RO.BC,'clamplr'); dof=[dof;1;2]; % fix both sides of beam
 else; error('BC %s unknown');
 end
 % Keep only z displacement and y rotation
 dof=[dof;.01;.02;.04;.06];
 %Fix dofs
 model=fe_case(model,'fixdof','fixed',dof);
 % Add loss factor (hysteretic damping)
 model=feutil(sprintf('setmat 100 eta=%.g',RO.eta),model);

elseif comstr(Cam,'sens'); [CAM,Cam]=comstr(CAM,5);
 %% #BeamGenSens : Place sensors (right side, ...)------------------------3
 if ~isfield(RO,'model'); error('Provide RO.model'); end
 model=RO.model;
 if strcmpi(RO.Sens,'right'); 
  n1=feutil('findnode',model,'x>0.5');
  tdof=feutil('getdof',n1,.03); % Sensors only in z direction
 else; error('unknown Sens command');
 end
 % Place sensors over the model
 model=fe_case(model,'SensDof append trans','Test',tdof);
elseif comstr(Cam,'load'); [CAM,Cam]=comstr(CAM,5);
 %% #BeamGenMesh : Place Load (first quarter of beam,...)-----------------3
 if ~isfield(RO,'model'); error('Provide RO.model'); end
 model=RO.model;
 if strcmpi(RO.Load,'midleft');
  % Add Load at the middle left of the beam
  i1=feutilb('AddNode -nearest',[0.25 0 0],model.Node(:,5:7));
  i1=model.Node(find(i1),1);
  model=fe_case(model,'DofLoad','Load',i1+.03);
 else; error('unknown Load command');
 end
else; error('Command d_mesh BeamGen%s unknown',CAM);
end
out=model;
%% #BeamEnd
else;error('Beam%s unknown',CAM);
end
%% #Volumes: pointers to volume test cases -----------------
elseif 1==2
%% #VolExpTime ->sdtweb d_fetime Expvol % explicit time integration -2
%% #VolSleeper ->sdtweb demosdt sleeper % explicit time integration -2

elseif comstr(Cam,'dd'); [CAM,Cam]=comstr(CAM,3);
%% #DD: sample mesh for domain decomposition functionalities -----------------
if carg<=nargin; RO=varargin{carg}; carg=carg+1; else; RO=struct; end
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
 'etree(2#%i#"level of etree decomposition")' ...
 '-noSplit(#3#"no to split model before output")'...
 ],{RO,CAM}); Cam=lower(CAM);
 if comstr(Cam,'plate')  
  %% #DDPlate -2
  % plate model with basic elimination tree
  [CAM,Cam,RO.div]=comstr('-div',[-25 1],CAM,Cam);
  if isempty(RO.div);  RO.div=[10 10];
  elseif length(RO.div)==1;  RO.div=repmat(RO.div,1,2);
  end
  r1=[0 0 1; 2 0 1;2 2 1;0 2 1];
  model=femesh('objectquad 1 1',r1,RO.div(1),RO.div(2));
  model=fe_mat('defaultpl',model); il=p_shell('default'); model.il=il.il;
  dof=feutil('getdof',model);
  if RO.etree
   switch RO.etree
    case 2
     RO.partSel={'x<1','x>1','x==1'};
     RO.etree=[3 3 0];
     K=cell(3); M=cell(3); ind=cell(3,1); r1=cell(3,1);
    case 4
     RO.partSel={'x<1 & y<1','x<1 & y>1','y==1 & x<1','x>1 & y<1','x>1 & y>1','y==1 & x>1','x==1'};
     RO.etree=[3 3 7 6 6 7 0];
     K=cell(7); M=cell(7); ind=cell(7,1); r1=cell(7,1);
    case 8
     RO.partSel={
      sprintf('x<1 & y>1 & y>-x+%g',...
      2+2/RO.div(1)+1/RO.div(1)/10)                  %  1
      sprintf('x<1 & y>1 & y<-x+%g',2-1/RO.div(1))   %  2
      sprintf('x<1 & y>1 & y>-x+%g & y<-x+%g',...
      2-1/RO.div(1),2+2/RO.div(1)+1/RO.div(1)/10)    %  3
      sprintf('x<1 & y<1 & y<x-%g',1/RO.div(1))      %  4
      sprintf('x<1 & y<1 & y>x+%g',...
      2/RO.div(1)+1/RO.div(1)/10)                    %  5
      sprintf('x<1 & y<1 & y>=x-%g & y<x+%g',...
      1/RO.div(1),2/RO.div(1)+1/RO.div(1)/10)        %  6
      'x<1 & y==1'                                   %  7
      sprintf('x>1 & y<1 & y<-x+%g',2-1/RO.div(1))   %  8
      sprintf('x>1 & y<1 & y>-x+%g',...
      2+2/RO.div(1)+1/RO.div(1)/10)                  %  9
      sprintf('x>1 & y<1 & y>-x+%g & y<-x+%g',...
      2-1/RO.div(1),2+2/RO.div(1)+1/RO.div(1)/10)    % 10
      'x>1 & y>1 & y<x'                              % 11
      sprintf('x>1 & y>1 & y>x+%g',2/RO.div(1))      % 12
      sprintf('x>1 & y>1 & y>=x-%g & y<x+%g',...
      1/RO.div(1),2/RO.div(1)+1/RO.div(1)/10)        % 13
      'x>1 & y==1'                                   % 14
      'x==1'          };                             % 15
     RO.etree=[3 3 7 6 6 7 15 10 10 14 13 13 14 15 0];
     K=cell(15); M=cell(15); ind=cell(15,1); r1=cell(15,1);
     
    otherwise; error('testplate with %i leaf nodes not pre set',RO.etree);
   end
   model=fe_case(model,'fixdof','zrot',.06);
   RO.DOF=[.01;.02;.03;.04;.05];
  end
  
 elseif comstr(Cam,'cube') 
  %% #DDCube -2
  [CAM,Cam,RO.div]=comstr('-div',[-25 1 3],CAM,Cam);
  if isempty(RO.div); RO.div=[10 10 10]; end
  model=femesh(sprintf('testhexa8b -divide %i %i %i',RO.div));
  dof=feutil('getdof',model);
  if RO.etree
   switch RO.etree
    case 2
     RO.partSel={sprintf('x<%.15g',.5-.1/RO.div(1))
      sprintf('x>%.15g',.5+.1/RO.div(1))
      'x==.5'   };
     RO.etree=[ 3 3 0];
     K=cell(3); M=cell(3); ind=cell(3,1); r1=cell(3,1);
    case 4
     RO.partSel={sprintf('x<%.15g & y<%.15g',.5-.1./RO.div([1 2]))
      sprintf('x<%.15g & y>%.15g',.5-.1/RO.div(1),.5+.1/RO.div(2))
      sprintf('x<%.15g & y==%.15g',.5-.1/RO.div(1),.5)
      sprintf('x>%.15g & y<%.15g',.5+.1/RO.div(2),.5-.1/RO.div(2))
      sprintf('x>%.15g & y>%.15g',.5+.1./RO.div([1 2]))
      sprintf('x>%.15g & y==%.15g',.5+.1/RO.div(1),.5)
      'x==.5' };
     RO.etree=[3 3 7 6 6 7 0];
     K=cell(7); M=cell(7); ind=cell(7,1); r1=cell(7,1);
     
    otherwise; error('testcube with %i leaf nodes not pre set',RO.etree);
   end
  end
  RO.DOF=[.01;.02;.03];
 else; error('Unknown %s model',Cam);
 end % test
 % Direct split preparation if asked
 if ~RO.noSplit
  [m,k]=fe_caseg('assemble -matdes 2 1',model);
  for j1=1:length(RO.partSel)
   r2=feutil(sprintf('findnode %s',RO.partSel{j1}),model);
   model=feutil('addsetnodeid',model,sprintf('Node_%i',j1),r2(:,1));
   r1{j1}=reshape(repmat(r2(:,1),1,length(RO.DOF))',1,[])'+repmat(RO.DOF,size(r2,1),1);
   dof=fe_case('gettdof',model);
   ind{j1}=fe_c(dof,r1{j1},'ind',1);
  end
  model=fe_dd('split',model,RO.etree,ind);
 end
 if nargout>0; out=model; else; feplot(model); end

%% #Rve : classical representative volume element meshes ---------------------
elseif comstr(Cam,'rve');[CAM,Cam]=comstr(CAM,4);
if comstr(Cam,'1fiber')
%% #Rve1Fiber : 3d mesh of a cube with an internal fiber -3

 [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   'vf(0.5#%g#"Fiber volume fraction") '...
   'r(#%g#"Fiber radius") '...
   'mat(BerthCarDx#%s#"Material command") '...
   'x(5#%g#"Cell x, mu m") y(2#%g#"Cell y")  z(5.1#%g#"Cell z") '...
   ' ny(1#%g#"n Cell y") scale(1e-3#%g#"space scale")'],{RO,CAM});

if isempty(RO.r) % Enforce volume fraction
 r1=(RO.x*RO.z)*RO.vf; % Fiber surface
 RO.r=sqrt(r1/pi); % Fiber radius
end

% Mesh cell: - - 
mo1=feutil(sprintf(...
    'ObjectHoleInBlock 0 0 0  1 0 0  0 1 0  %g %g %g  %g  8 8 %i 2',...
    RO.x,RO.z,RO.y,RO.r,RO.ny)); 
% mdl1=feutil('ObjectCylinder x1 y1 z1 x2 y2 z2 r divT divZ',model)
% mdl1=feutil('ObjectCylinder 0 0 0  0 1 0 0.7 32 3')
mo2=feutil(sprintf('ObjectDisk-nodeg 0 %g 0  %g  0 1 0  32  3',-RO.y/2,RO.r));%x y z r nx ny nz Nseg NsegR
mo2=feutil(sprintf('Extrude %i 0 %.15g 0',RO.ny,RO.y/RO.ny),mo2); 
mo2=feutil('Proid',mo2,[1 2]);mo2=feutil('MatId',mo2,[1 2]);
mo1=feutil('AddTestMerge epsl.001 -NoOrig;',mo1,mo2);
if isfield(RO,'quad')&&RO.quad;mo1=feutil('lin2quad',mo1);end
mo1=p_solid('default;',stack_rm(mo1,'','OrigNumbering'));
if isfield(RO,'pl'); mo1.pl=RO.pl; % Material is given MatId 1 & 2
else; mo1.pl=d_mesh(['matrve' RO.mat]);
end

% Build cyclic symetry: - -
mo1=fe_homo('DftBuild',mo1,[RO.x RO.y RO.z]);

mo1=fe_case(mo1,'park','Matrix','matid 1');
if RO.scale==1e-3
 mo1.Node(:,5:7)=mo1.Node(:,5:7)*RO.scale; % Original values in mu m
 mo1.unit='TM';mo1.pl=fe_mat('convertTM',mo1.pl);
else; mo1.unit='US';
end
r0=fe_case(mo1,'getdata','Symmetry');

pl=sortrows(feutil('getpl',mo1),1); % 1 matrix, 2 fiber
em=pl(1,3); ef=pl(2,3); hm=pl(1,7); hf=pl(2,7);vf=RO.vf;
et=1/(vf/ef+(1-vf)/em);
mo1.ViscInfo.eta= ...
  {'axial' hf*(vf/(vf+em/(ef*(1-vf))))+hm*(1-vf)/(1-vf + ef/em*vf);
   'transverse',hf*(et/ef)*vf+hm*(et/em)*(1-vf);
  };
mo1.name='1Fiber';

elseif comstr(Cam,'layered'); [CAM,Cam]=comstr(CAM,8);
 %% #RVELayered - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
  'mat("@d_mesh(''matRVElayers'')"#%s#"Material command") '...
  ...'-damp(#31#"to apply damping to i material")' ...
  'x(5#%g#"Cell x, mu m") y(5#%g#"Cell y")  z(5#%g#"Cell z") '...
  'nx(1#%g#"n Cell x") ny(1#%g#"n Cell y") nz(2#%g#"n Cell z")'],{RO,CAM});
 % Mesh Cell
 node=[0 0 0;RO.x 0 0;RO.x RO.y 0;0 RO.y 0]; 
 if length(RO.z)>1; RO.zHeights=RO.z(:);RO.nz=length(RO.z)-1;
 else; RO.zHeights=linspace(0,RO.z,RO.nz+1);
 end
 mo1=feutil('objectquad 1 1',node,RO.nx,RO.ny);
 mo1=feutil('extrude 0 0 0 1',mo1,RO.zHeights);
 mo1.unit='SI'; mo1.name='layered_cube';
 % Set materials
 mo1=p_solid('default;',mo1);
 if comstr(RO.mat,'@');  mo1.pl=eval(RO.mat(2:end)); % mo1=feutil('setmat',mo1,RO.mat)
 else; mo1.pl=d_mesh(['matrve' RO.mat]);
 end
 % Assign materials per z layer (sequence on existing materials)
 mpid=feutil('mpid',mo1);
 for j1=1:length(RO.zHeights)-1
  i1=feutil(sprintf('findelt innode{z>=%.15g & z<=%.15g}',RO.zHeights(j1:j1+1)),mo1);
  mpid(i1,1)=remi(j1,size(mo1.pl,1));
 end
 mo1.Elt=feutil('mpid',mo1.Elt,mpid);
 
 % Apply 1D conditions
 mo1=fe_case(mo1,'fixdof','1D',(1:2)'/100);
 % Default load
 mo1=fe_case(mo1,'DofLoad','Test',struct('def',1,'DOF',mo1.Node(1)+.01));
  
elseif comstr(Cam,'constitinterp'); [CAM,Cam]=comstr(CAM,8);
% This is tested in t_thermal

%% #RVEConstitInterp : interpolation of constitutive law - - - - - - - - - - -
mo1=femesh('testhexa8 divide 1 1 10');
r1=mo1.Node(:,7).^2; % externally defined 
data=struct('dir',{{'7800', ...
    @(x)r1, ... % Values at nodes 
    '210e9*(z.^1)'}}, ... % Analytic formula using x,y,z,r
    'lab',{{'rho','eta','D33'}});
% The labels must come from the list ConstitLab for example
% feval(p_solid('@ConstitLab'),'m_elastic.1')

pro=struct('il',111,'type','p_solid','MAP',data);
mo1=stack_set(mo1,'pro','WithMap',pro);
 
%% EndScript

elseif comstr(Cam,'celas'); [CAM,Cam]=comstr(CAM,8);
%% #RveCelas : two spring for homogeneization testing

[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
         '-Nc(50#%g#number of cells")' ...
         '-E1(10#%g#"first spring") -E2(20#%g#"second spring") ' ...
         '-loss1(0#%g#"first spring") -loss2(.1#%g#"second spring") ' ...
         ' -xfrac(.2#%g#"position of center point within dx")' ...
         ' -dx(20#%g#"step length ")' ...
      ],{RO,CAM}); Cam=lower(CAM); 
RO.k1=RO.E1/(RO.dx*RO.xfrac); % K=E A /L
RO.k2=RO.E2/(RO.dx*(1-RO.xfrac)); % K=E A /L

r1=(0:RO.Nc)*RO.dx;r1=[r1;r1+RO.dx*RO.xfrac];r1=r1(:);r1(end)=[];
mo1=struct('Node',[1 0 0 0  0 0 0],'Elt',feutil('addelt','mass1',1));
mo1=feutil('extrude 0  1 0 0',mo1,r1(:));
mo1.Elt=feutil('set group1 name celas proid 1 matid 1',mo1);
% Alternate spring types. Use orientation 3 for easier display
mo1.Elt(2:2:end,3:10)=repmat([3 3 1 0 RO.k1 0 RO.k1*RO.loss1 0],RO.Nc,1);
mo1.Elt(3:2:end,3:10)=repmat([3 3 2 0 RO.k2 0 RO.k2*RO.loss2 0],RO.Nc,1);
d1=struct('def',cos(mo1.Node(:,5)/(RO.Nc*RO.dx)*2*pi),'DOF',mo1.Node(:,1)+.03);

% Homogeneous spring
mo1=fe_case(mo1,'fixdof','1D',[1 2 4 5 6]'/100);
mo1.name='';mo2=mo1; 
mo2.Elt(2:2:end,3:10)=repmat([3 3 1 0 1/(RO.dx*RO.xfrac)     0 0 0],RO.Nc,1);
mo2.Elt(3:2:end,3:10)=repmat([3 3 2 0 1/(RO.dx*(1-RO.xfrac)) 0 0 0],RO.Nc,1);
[mo2,C2]=fe_case(mo2,'Assemble -MatDes 1 -SE');
% Heterogenous spring
mo1=fe_case(mo1,'Assemble -MatDes 1 3 -SE');mo1.Klab{2}='B';
mo1.K{3}=mo2.K{1};mo1.Klab{3}='homog';mo1.name='RveCelas';mo1.DOF=C2.DOF;

if nargout==0;feplot(mo1,d1);fecom(';colordataevalz;undefline');
else; out=mo1;out1=d1;out2=RO;
end
return

elseif comstr(Cam,'nida'); [CAM,Cam]=comstr(CAM,8);
%% #RveNida : honeycomb cell, 
% This example was originally contributed by Mikhail Guskov
% honeycomb core general parameters, notations ref [th. Florens, p23]

[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   ' dim([5, 5, 20.4, .1, .2]#%g#"cell dimensions [a b thick t tp]")'...
   ' ang([30 0]#%g#"deg, cell angle, orientation")'...
   ' h(1.12#%g#"skin thickness")'...
   ' lc(3#%g#"element target size")'...
   ' Quad(0#3#"quadratic elements if present")'...
   ],{RO,CAM(7:end)});

tmp.a  = RO.dim(2);
tmp.aa = RO.ang(1);
tmp.t  = RO.dim(4); 
tmp.t1 = RO.dim(5); 
tmp.th = atand((tmp.t - tmp.t1*sind(tmp.aa))/(tmp.t1*cosd(tmp.aa))); 
tmp.d1 = tmp.t*cosd(tmp.aa)/2;
tmp.d2 = tmp.t1*tand(tmp.th)/2;

% keypoint definition (notations see notes mg2012-10-08)
% kp{1} : OABCDEF
kp{1}(1,1:3) = [0,0,0]; % O
kp{1}(2,1:3) = kp{1}(1,1:3)+ tmp.a*[sind(tmp.aa) cosd(tmp.aa),0]; % A
kp{1}(3,1:3) = kp{1}(2,1:3) + [tmp.a,0,0]; % B
kp{1}(4,1:3) = kp{1}(3,1:3)...
    +[tmp.a*sind(tmp.aa),-tmp.a*cosd(tmp.aa),0]; % C
kp{1}(6,1:3) = kp{1}(1,1:3)...
    + [tmp.a*sind(tmp.aa),-tmp.a*cosd(tmp.aa),0]; % D
kp{1}(5,1:3) = kp{1}(6,1:3) + [tmp.a,0,0]; % E
kp{1}(7,1:3) = kp{1}(4,1:3) + [tmp.a,0,0]; % F

% contour 1 : cell section wireframe
ic = 1; ib = 0;
ib = ib+1; cont(ic).br(ib).p = [1:6 1];
ib = ib+1; cont(ic).br(ib).p = [4 7];

% kp{2} : detailed shape keypoints
ip = 0;
% O123  1+
ip=ip+1; kp{2}(ip,1:3) = kp{1}(1,:) + [tmp.d1,0,0]; % O1
ip=ip+1; kp{2}(ip,1:3) = kp{1}(1,:) + [-tmp.d2, tmp.t1/2,0]; % O2
ip=ip+1; kp{2}(ip,1:3) = kp{1}(1,:) + [-tmp.d2,-tmp.t1/2,0]; % O3
% A123  4+
ip=ip+1; kp{2}(ip,1:3) = kp{1}(2,:) + [tmp.d2,-tmp.t1/2,0]; % A1
ip=ip+1; kp{2}(ip,1:3) = kp{1}(2,:) + [tmp.d2, tmp.t1/2,0]; % A2
ip=ip+1; kp{2}(ip,1:3) = kp{1}(2,:) + [-tmp.d1,0,0]; % A3
% B123  7+
ip=ip+1; kp{2}(ip,1:3) = kp{1}(3,:) + [-tmp.d2,-tmp.t1/2,0]; % B1
ip=ip+1; kp{2}(ip,1:3) = kp{1}(3,:) + [tmp.d1,0,0]; % B2
ip=ip+1; kp{2}(ip,1:3) = kp{1}(3,:) + [-tmp.d2, tmp.t1/2,0]; % B3
% C123  10+
ip=ip+1; kp{2}(ip,1:3) = kp{1}(4,:) + [-tmp.d1,0,0]; % C1
ip=ip+1; kp{2}(ip,1:3) = kp{1}(4,:) + [tmp.d2,-tmp.t1/2,0]; % C2
ip=ip+1; kp{2}(ip,1:3) = kp{1}(4,:) + [tmp.d2, tmp.t1/2,0]; % C3
% D123  13+
ip=ip+1; kp{2}(ip,1:3) = kp{1}(5,:) + [-tmp.d2, tmp.t1/2,0]; % D1
ip=ip+1; kp{2}(ip,1:3) = kp{1}(5,:) + [-tmp.d2,-tmp.t1/2,0]; % D2
ip=ip+1; kp{2}(ip,1:3) = kp{1}(5,:) + [tmp.d1,0,0]; % D3
% E123  16+
ip=ip+1; kp{2}(ip,1:3) = kp{1}(6,:) + [tmp.d2, tmp.t1/2,0]; % E1
ip=ip+1; kp{2}(ip,1:3) = kp{1}(6,:) + [-tmp.d1,0,0]; % E2
ip=ip+1; kp{2}(ip,1:3) = kp{1}(6,:) + [tmp.d2,-tmp.t1/2,0]; % E3
% F123  19+
ip=ip+1; kp{2}(ip,1:3) = kp{1}(7,:) + [tmp.d1,0,0]; % F1
ip=ip+1; kp{2}(ip,1:3) = kp{1}(7,:) + [-tmp.d2, tmp.t1/2,0]; % F2
ip=ip+1; kp{2}(ip,1:3) = kp{1}(7,:) + [-tmp.d2,-tmp.t1/2,0]; % F3

% contour 2: detailed shape
%E1-O1-A1-B1-C1-D1-D3-C2-F3-F2-C3-B2-B3-A2-A3-O2-O3-E2
ic = 2; ib=1;
cont(ic).br(ib).p = [16 1 4 7 10 13 15 11 21 20 12 8 9 5 6 2 3 17 16];

r1=kp{2};
r1(end+1,:)=r1(8,:)-diff(r1([10;1],:));
r1(end+1,:)=r1(15,:)-diff(r1([10;1],:));
mo1=struct('Node',[(1:size(r1,1))'*[1 0 0 0] r1],'Elt', ...
    feutil('addelt',feutil('addelt','quad4',[16 17 3 1   1 1; ...
    1 2 6 4   1 1;4 5 9 7   1 1;7 8 12 10   1 1;10 11 15 13   1 1;
    11 12 20 21   1 1;
    16 13 10 1 2 2;1 10 7 4 2 2
    8 22 20 12 2 2;15 11 21 23 2 2]), ...
    'tria3',[1 2 3   1 1;4 5 6   1 1;7 8 9   1 1;10 11 12   1 1]), ...
    'Stack',{{}},'pl',[],'il',[],'unit','MM');
RO.dir=[norm(diff(r1([2 20],:))) norm(diff(r1([5 16],:)))];

% Refine the cell edges
i1=ceil(RO.dim(2)/(RO.lc*2));
[elt,mo1.Elt]=feutil('removeelt eltind 2:7',mo1);
mo1=feutil(sprintf('divide %i 1',i1),mo1);
el1=mo1.Elt;mo1.Elt=feutil('selelt eltname quad4',mo1.Node,elt);
mo1=feutil(sprintf('divide %i %i',i1,i1),mo1);
mo1=feutil('addelt',mo1,el1);
mo1=feutil('addelt',mo1,feutil('selelt eltname tria',mo1.Node,elt)); % Add triangles

if isfield(RO,'extrude')&&strcmpi(RO.extrude,'svs')
else
 % Now extrude a volume
 r1=RO.h(1)+linspace(0,RO.dim(3),ceil(RO.dim(3)/RO.lc));
 r1=[0 r1 r1(end)+RO.h(1)];
 mo1=feutil('extrude 0 0 0 1',mo1,r1);
 mo1.Elt=feutil('removeelt matid 2 & withoutnode {z==0 | z==}',mo1,r1(end));
 mo1.Node=feutil('getnodegroupall',mo1);
 % Set skin property to something uniform
 mpid=feutil('mpid',mo1);
 mpid(feutil('findelt withnode {z==0 | z==}',mo1,r1(end)),1:2)=2;
 mo1.Elt=feutil('mpid',mo1.Elt,mpid);
 % Define properties
 %[pl,il]=nida13('matbottom');
 pl=[2 fe_mat('m_elastic','SI',6) 69e9 69e9 8.1e9 0.03 0.03 0.03 4.8e9 4.8e9 4.8e9 1553.571];
 pl(1)=2; mo1=feutil('setmat',mo1,fe_mat(['convert' mo1.unit],pl));
 mo1.pl=m_elastic(mo1.pl,'dbval 1 -unit MM Aluminum');
 mo1=p_solid('default;',mo1);
 if RO.Quad;mo1=feutil('lin2quad',mo1);end
 mo1=fe_homo('DftBuild',mo1,RO.dir);
end
 mo1=stack_set(mo1,'info','MeshInfo',RO);


out=mo1; out1=RO;

%% #RveEnd
else;error('RVE%s unknown',CAM);
end
if isfield(RO,'Per') % Per provides the number of periodic directions
 %% RvePer : add periodicity conditions to the RVE
 r2=[RO.x RO.y RO.z];
 if isfield(RO,'epsl');st1=sprintf('epsl%.15g',RO.epsl);else;st1='';end
 mo1=fe_homo(['DftBuild' st1 ';'],mo1,r2(1:RO.Per));
end
% clean output of RVE models
if nargout>0
 out=mo1;
else
 cf=feplot(mo1);  %fecom(cf,'colordatamat')
 %mo2=mo1; mo2.Elt=[]; mo2=feutil('AddElt',mo2,'beam1',r0.IntYNodes); cf=feplot(mo2);% check
 mo1=fe_mat('defaultil',mo1);  mo1=fe_mat('defaultpl',mo1);
 dd=fe_eig(mo1,[5 30 1e3]);
 cf=feplot(mo1,dd); fecom(cf,'colordatamat')
end

%% #PreTest -------------------------------------------------------------------
elseif comstr(Cam,'pretest');[CAM,Cam]=comstr(CAM,8);


if comstr(Cam,'accbeam')
%% #PreTestAccbeam : beam on a shake table
 model=feutil('Objecthexa 1 1',[-.5 -.5 0;1 0 0;0 1 0;0 0 .1],3,3,1);
 model=feutil('objectbeam 2 2',model,[0 0 .1;0 0 1],10);
 model.pl=m_elastic('dbval 1 aluminum','dbval 2 Steel');
 model=feutil('setmat 1 rho = 100',model); % Lower density for holes in table
 
 model.il= p_beam('dbval 2 rectangle .05 .1');
 model=p_solid('default;',model);
 
 % Fix beam root
 [model.Node,i2]=feutil('addnode',model.Node,[0 0 0;0 0 .1]);
 model=fe_caseg('Connection surface 123456 -MaxDist0.5',model,'Beam', ...
   'matid2 & z<=.1', ...                          % Selection of nodes to connect
   'matid1'); % Selection of elements for matching

 % Define loads and sensors on table and beam
 [model.Node,i2]=feutil('addnode',model.Node,[-.5 -.5 0;-.5 .5 0;.5 -.5 0;0 0 1]);
 i2=model.Node(i2,1);
 sdof=[i2(1)+(1:3)'/100;i2(2)+[.01;.03];i2(3)+.03;i2(4)+(1:3)'/100];
 
 model=fe_case(model,'SensDof','Sensors',sdof, ...
     'DofLoad','In',struct('def',eye(7),'DOF',sdof(1:7)));
 %fecom('shownodemark',model.Node(i2,1))
 

 if nargout==0;def=fe_eig(model);feplot(model,def)
 else; out=model;
 end
 
 
else;error('PreTest%s unknown',CAM);
end
    
%% #Periodic -------------------------------------------------------------------
elseif comstr(Cam,'periodic');[CAM,Cam]=comstr(CAM,6);

if comstr(Cam,'plate')
%% #PeriodicPlate example of computation of periodic modes for a simple plate
%
% WARNING : PERIODIC COMPUTATIONS IS A NON SUPPORTED CAPABILITY
%  if you are interested in mounting a related project, please contact SDTools

model=quad4('test'); model=fe_case(model,'reset');
model=stack_set(model,'info','EigOpt',[5 10 0]);

% Define periodicity condition by giving step of repetition
%  display corresponding nodes given in 'cyclic','Symmetry' case entry.
model=fe_cyclic('build -1 1.0 0.0 0.0 ',model);
feplot(model);fecom('curtabCase','Symmetry');fecom ProViewOn

% Use fe_cyclic eig to compute modes with a periodicity of 1,3 and 3.2 cells.
d1=fe_def('appenddef',  ...
  fe_cyclic('eig 1',model), ...
  fe_cyclic('eig 3',model), ...
  fe_cyclic('eig 3.2',model));
d1.LabFcn='sprintf(''%i@ %.4g Hz (nCell=%i)'',ch,def.data(ch,1:2))';

cf=feplot; cf.def=d1;fecom(';view3;ColorDataEvalZ');

fecom('ch12');  % Visualize the travelling wave with 3 cell periodicity

% Note how frequencies come in pairs and depend on the wavelength
disp([reshape(d1.data(11:30,1),2,10)' reshape(d1.data(31:50,1),2,10)'])


%% #PeriodicEnd
else;error('Periodic%s unknown',CAM);
end

%% #ShapedBar : shaped bar with hole for stress study
elseif comstr(Cam,'shapedbar');[CAM,Cam]=comstr(CAM,10);

    
 [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
         'L1(300#%g#1st segment")' ...
         'L2(300#%g#2nd segment")' ...
         'L3(250#%g#3rd segment")' ...
         'L1R(3#%g#fillet start")' ...
         'L3R(3#%g#fillet start")' ...
         'R3(12.5#%g# hole radius")' ...
         'H1(25#%g#1st height")' ...
         'H2(50#%g#2nd height")' ...
         'lc(15#%g# characteristic mesh size")' ...
         'lcf1(0#%g# characteristic mesh size arounf filet 1")' ...
         'lcf2(0#%g# characteristic mesh size arounf hole")' ...
         'th(5#%g#thickness")' ...
         'nz(1#%i#number of elements in z direction")' ...
         'quad(#3#"quadratic elements")' ...
         'solve(1#%g#"0 mesh, 1 static, 2 eig")' ...
         'cf(#%g#"feplot figure number")' ...
      ],{RO,CAM}); Cam=lower(CAM); 

 if 1==2
  a=sdt_locale('defbutInMFile.ShapedBar','d_mesh');
  %fe_range('FromRO',a)
  Range=fe_range('grid',struct('L1R',[20 60],'quad',[0 1]));
  Range=fe_range('AddCompute',Range,struct('SetFcn',{{'d_mesh','ShapedBar'}}));
  Res=fe_range('Loop',Range);Res=Res.Res;
 end

 mdl=[];
 if RO.lcf1==0;RO.lcf1=min(RO.lc,RO.L1R/3);end % default 1/3 of radius 
 if RO.lcf1==0;RO.lcf1=RO.lc;end 
 if RO.lcf2==0;RO.lcf2=RO.lc;end 

 mdl.Node=[1  0 0 0          0                     RO.H1/2  0;
           2  0 0 RO.lcf1    RO.L1-RO.L1R          RO.H1/2  0;
           3  0 0 0          RO.L1                 RO.H2/2  0;
           4  0 0 0          RO.L1+RO.L2           RO.H2/2  0;
           5  0 0 0          RO.L1+RO.L2+RO.L3R    RO.H1/2  0;
           6  0 0 0          RO.L1+RO.L2+RO.L3     RO.H1/2  0;
           7  0 0 0          RO.L1+RO.L2+RO.L3     0        0;
           8  0 0 0          RO.L1+RO.L2/2+RO.R3   0        0;
           9  0 0 0          RO.L1+RO.L2/2-RO.R3   0        0;
           10 0 0 0          0                     0        0
           11  0 0 0         RO.L1+RO.L2/2         0        0;
           ];
 mdl.Elt=[];
 mdl=fe_gmsh('addline-loop1',mdl,[10 1;1 2]);
 %NNode=sparse(mdl.Node(:,1),1,1:size(mdl.Node,1));
 
 %% Possibly add first filet
 if RO.L1R==0 % No Filet
    mdl=fe_gmsh('addline-loop1',mdl,[2 3;3 4]);
 else
  if RO.L1R>=abs(RO.H1-RO.H2)/2;i2=3; % large radius
  else;
   [mdl.Node,i2]=feutil('addnode',mdl.Node,[RO.L1 RO.H1/2+RO.L1R 0]);
   mdl.Node(i2,4)=RO.lcf1; % characteristic length
  end
  mdl=fe_gmsh('addCircleArc-loop1-tangent',mdl,...
    [1 0 0; mdl.Node(mdl.Node(:,1)==2,5:7);mdl.Node(mdl.Node(:,1)==i2,5:7)]);
  if i2==3;i2=[3 4];else;i2=[i2 3;3 4];end;mdl=fe_gmsh('addline-loop1',mdl,i2);
 end
 
 %% Add second filet
 if RO.L3R>=abs(RO.H1-RO.H2)/2;i2=4; % large radius
 else;[mdl.Node,i2]=feutil('addnode',mdl.Node,[RO.L1+RO.L2 RO.H1/2+RO.L3R 0]);
     mdl=fe_gmsh('addline-loop1',mdl,[4 i2]);
 end
 mdl=fe_gmsh('addCircleArc-loop1-tangent2',mdl,...
    [1 0 0; mdl.Node(mdl.Node(:,1)==i2,5:7);mdl.Node(mdl.Node(:,1)==5,5:7)]);
 mdl=fe_gmsh('addline-loop1',mdl,[5 6;  6 7; 7 8]);
 %% Add half circle
 n1=feutil('getnode nodeid 11 8 9',mdl);
 % add mid point
 [mdl.Node,i2]=feutil('addnode',mdl.Node,[RO.L1+RO.L2/2 RO.R3 0]);
 n1(end+1,:)=mdl.Node(i2,:);mdl.Node(i2,4)=RO.lcf2;

 mdl=fe_gmsh('addCircleArc-loop1',mdl,n1([1 2 4],5:7));
 mdl=fe_gmsh('addCircleArc-loop1',mdl,n1([1 4 3],5:7));
 %mdl=fe_gmsh('addCircleArc-loop1',mdl,n1(1:3,5:7));
 mdl=fe_gmsh('addline-loop1',mdl,[9 10]);
 mdl.Stack{end}.PlaneSurface=1;
 
 % Give GMSH options and run
 RB=struct('lc',RO.lc,'Run','-2 ','sel','eltname tria'); 
 if RO.quad; RB.Run=[RB.Run ' -order 2'];end % Use second order
 mo1=fe_gmsh('write',mdl,RB);

 % Make symmetric
 mo1b=mo1; mo1b.Node(:,6)=-mo1b.Node(:,6);
 mo1=feutil('AddTestMerge',mo1,mo1b);
 mo1=feutil(sprintf('Extrude %i 0 0 %g',RO.nz,RO.th),mo1);

 % proid matid  
 mpid=feutil('mpid',mo1);
 i1=feutil('FindElt groupall',mo1);
 mpid(i1,1:2)=1; % proid=1 matid=1 
 mo1.Elt=feutil('mpid',mo1,mpid);
 mo1=feutil('joinall',mo1);[eltid,mo1.Elt]=feutil('eltidfix;',mo1);
 
 mo1.unit='MM';
 
 %% Boundary conditions
 n1=feutil('getnode x==',mo1,max(mo1.Node(:,5)));
 %% EdgeLoad 
 n1=feutil('getnode x==',mo1,max(mo1.Node(:,5)));
 d1=struct('def',ones(size(n1,1),1),'DOF',n1(:,1)+.01);
 d1.def=d1.def/norm(d1.def); % Unit load distributed on edge nodes
 mo1=fe_case(mo1,'FixDOF','clamped','x==0', ...
     'FixDof','Guide',sprintf('x==%.15g -DOF 2 3',max(mo1.Node(:,5))), ...
     'DofLoad','EdgeLoad',d1);
 %% material
 mo1.pl=m_elastic('dbval 1 steel -unit MM');
 mo1.il=p_solid('dbval 1 d3 -3');
 % Modes
 
 switch RO.solve
     case 0 % mesh
         def=[];
     case 1 % Static
         def=fe_simul('static',mo1);
     case 2 % Eigenvalue
         def=fe_eig(mo1,[5 20]);
     otherwise
         error('solve=%i unknown',RO.solve);
 end
 
 if ~isempty(RO.cf)
  cf=comgui('guifeplot',RO.cf); cf=feplot(mo1);
  if ~isempty(def);cf.def=def;end
  fecom('colordatastress MisesAtNode -edgealpha.1');fecom('colormapjet(7)');
  %fecom('renderer zbuffer');
 end
 
 if 1==2
  % z=fe_stress('stressMises AtNode',mo1,dd);
  fecom('colordataevalx -edgealpha0')
  % fecom('colordatastress MisesAtNode')
  fecom(cf,'SetProp sel(1).fsProp','EdgeColor','none');
  fecom colorbar
 end
 
 if ~isempty(evt) %From range loop
  def.name=strrep(strrep(RO.FileName{1},'root_',''),'_Compute=1','');
  assignin('caller','Res',{mo1,def});
 else;
  out=mo1;
 end
 
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------
%% #TestBlade an isolated blade example
elseif comstr(Cam,'testblade'); [CAM,Cam]=comstr(CAM,10);

[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
 'd1(12#%i#"divisions along length")' ...
 'd2(1#%i#"divisions along tickness")' ...
 'd3(5#%i#"divisions along width")' ...
 'mat(iso#%s#"material variants")' ...
 ],{RO,CAM}); Cam=lower(CAM);

node = [1  0  0;20  0  0;0 .2 0; 0  0 5];
model=feutil('Objecthexa 1 1',node,RO.d1,RO.d2,RO.d3); % creates model 
r1=basis('rect2cyl',model.Node);r1(:,6)=r1(:,6)-1*(r1(:,5)).*r1(:,7)/10;
model.Node=basis('cyl2rect',r1);
model.unit='CM'; model=m_elastic('default',model);
model=p_solid('default',model);
[eltid,model.Elt]=feutil('eltidfix;',model);

if strncmpi(RO.mat,'ortho',5)
%% Define a case with multiple orthotropic materials and orientation
 cEGI=find(eltid~=0);
 model.pl=repmat(fe_mat('convert SICM',[1 fe_mat('m_elastic','SI',6), ...
       80e9 50e9 8e9, 0.4 0.06  0.027, 4.7e9  8.6e9  7.5e9, 1550]),3,1);
 model.pl(2:3,1)=2:3; model.pl(2,3:5)=model.pl(2,3:5)*1.1;
 model.pl(3,3:5)=model.pl(2,3:5)*1.1;
 model.Elt(feutil('findelt innode{x<5}',model),9)=3;
 model.Elt(feutil('findelt innode{x>14}',model),9)=2;
 data=struct('EltId',eltid(cEGI),'bas',eltid(cEGI));
 NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
 for jElt=1:length(cEGI)
  n1=model.Node(NNode(model.Elt(cEGI(jElt),1:3)),:);
  p=diff(n1(:,5:7));p=sp_util('basis',p(1,:),p(2,:));
  data.bas(jElt,7:15)=p(:)';
 end
 model=stack_set(model,'info','EltOrient',data);
end

model=fe_case(model,'fixdof','base','r<2'); model.name='blade';

if nargout==0; feplot(model);fecom showpatch
else; out=model;
end 
 
 
%% #Mat
elseif comstr(Cam,'mat');[CAM,Cam]=comstr(CAM,4);
  
if comstr(Cam,'rve');[CAM,Cam]=comstr(CAM,4);% #MatRve -2
 if comstr(Cam,'berth'); 
  % #MatRveBerth comp12('matrveBerth car dx') -3
  % Berthelot2007423, table 2 page 424
typ=fe_mat('m_elastic','SI',1);
pl=[ ...
 2 typ 345e9  0.33 1200  0 0 % Carbon Fiber table 1, xxx Et=15 GPa
 2 typ 73e9 0.33 1200  0 0 % Glass fibre table 1
 2 typ 135e9 0.37 1200  12e9 0 % Kevlar fibre table 1
 1 typ 3.21e9 .338 1200 1.2e9 0.014  % Matrix DX210, x eta_lt=1.06%
 1 typ 2.8e9 .33 1200 1.08e8 0.0157  % Matrix SR1500
 ];
st={'car','gla','kev'};pl(:,1)=0;
for j1=1:3; if ~isempty(strfind(Cam,st{j1}));pl(j1,1)=2;end;end
if ~any(pl(1:3,1));pl(1,1)=2; end
        
st={'dx','sr'};
for j1=1:2; if ~isempty(strfind(Cam,st{j1}));pl(3+j1,1)=1;end;end
if ~any(pl(4:5,1));pl(4,1)=1; end
out=pl(pl(:,1)~=0,:);

elseif comstr(Cam,'1fiber') 
%% #MatRve1Fiber orthotropic material for a single fiber -3
if carg<=nargin;RO=varargin{carg};carg=carg+1;
else; RO=struct('vf',.5);
end
if carg<=nargin; pl=varargin{carg};carg=carg+1;
else; pl=comp12('matRveBerthCarDx');
end
vf=RO.vf; Ef=pl(1,3);nuf=pl(1,4); Gf=Ef/2/(1+nuf);
Em=pl(2,3); num=pl(2,4); Gm=Em/2/(1+num);

rho=vf*pl(1,5)+(1-vf)*pl(2,5); % Simple mass rule
E1=[vf*Ef+(1-vf)*Em];E2=(Ef*Em)/(Ef*(1-vf)+Em*vf); E3=1.02*E2;
Nu12=[vf*nuf+(1-vf)*num]; % Jones 3.19 page 94
Nu31=Nu12*E2/E1; Nu23=Nu12/2;
G12=(Gf*Gm)/(Gf*(1-vf)+Gm*vf); % Jone 3.26
{'Direct inverse G', G12 Gf*vf+Gm*(1-vf)} 
G2=G12*1.01; G3=G12*1.02; 
%[MatId typ Ex Ey Ez NuYZ NuZX NuXY GYZ GZX GXY rho a1 a2 a3 T0 eta]
Ex=E3; Ey=E1; Ez=E2; NuYZ=Nu12; NuZX=Nu23; NuXY=Nu31;
pl=[1 fe_mat('m_elastic','MM',6) Ex Ey Ez NuYZ NuZX NuXY G12 G2 G3 rho];
[a,b,c,dd]=p_solid('buildconstit',[1;1],pl,p_solid('dbval 1 d3 -3'));
disp(dd.dd)
out=pl([1;1],:);out(:,1)=[2;1];

 elseif comstr(Cam,'layers');
 %% #MatRveLayers
  typ=fe_mat('m_elastic','SI',1);
  out=[ ...
   1 typ 345e9   0.33 1200  0 .2 % Carbon Fiber table 1, xxx Et=15 GPa
   2 typ 86.25e9 0.33  600  0  0
   ];

%% #MatRveEnd -3
else;error('MatRve%s',CAM);
 end
else; error('Mat%s unknown',CAM);
end

elseif comstr(Cam,'stepmesh'); [CAM,Cam]=comstr(CAM,8);
%% #stepMesh : default meshing step  ---------------------------------------

Range=varargin{2};carg=2;
evt=varargin{carg};carg=carg+1;
model=evt.data;

assignin('caller','Res',model);

elseif comstr(Cam,'meshcfg'); [CAM,Cam]=comstr(CAM,8);
%% #MeshCfg : obtain/generate mesh  ---------------------------------------
% Default MeshCfg callback for simulation experiments
% see also sdtweb d_fetime SimuCfg
% sdth.findobj('_sub:~','d_hbm(Mesh):(Case):(NL)');comstr(ans,-30)
Range=RO;%varargin{carg};carg=carg+1;
RO=varargin{carg};carg=carg+1;
if ischar(RO);RO=struct('urn',RO);end
st=sdth.findobj('_sub:~',RO.urn);js=1; 
% 1 : create/load base mesh
if length(st)>1&&exist(st(1).subs,'file')% d_hbm(Mesh0D:cub)
  mo1=feval(st(js).subs,['Mesh' st(js+1).subs],Range,RO); 
  RO.MeshCb=st(js); RO.name=st(js+1).subs; js=js+2;
else
  st1=st;
  st=sdtroot('param.Project.MeshCb -safe');
  if isempty(st); 
      error('Project.MeshCb not set and ''%s'' not found',st1(1).subs);
  end
end
%% 2: deal with case building (possibly Mesh::NL for empty)
  if strcmpi(st(js).type,'()'); 
    RO.MeshCb=st(js);js=js+1; % :CaseFun(Case)
  end
  if js>length(st); RO.Case='';else; RO.Case=st(js).subs; js=js+1; end
  if ~isempty(RO.Case)
   RO.name=sprintf('%s:%s',RO.name,RO.Case);
   if js<length(st)&&strcmp(st(js).type,'{}') %% Mesh:V{data}:NL
    RO.CaseVal=st(js).subs; js=js+1;  
   else; RO.CaseVal={};
   end
   mo1=feval(RO.MeshCb.subs,'Case',mo1,RO); RO.MeshCb=st(1); 
  end

%% 3; now add the NLdata 
  if js<=length(st)&&strcmpi(st(js).type,'()'); 
    RO.MeshCb=st(js);js=js+1; % :MatFun(Mat)
  end
NLdata=[];
if js<=length(st); %  'd_hbm(Mesh0D):d_hbm(NL0Dm1t)' % sdtweb d_hbm NL 
  il=feutil('getil',mo1);
  if js+1<=length(st)&&strcmpi(st(js+1).type,'()');
    RO.MeshCb=st(js);js=js+1;
  end
  RO.NL=st(js).subs; 
  if ~iscell(RO.NL);RO.NL={il(1) RO.NL};end 
  for j2=1:2:length(RO.NL)
    if isempty(RO.NL{j2+1}); continue;end
    RO.name=sprintf('%s:%s',RO.name,RO.NL{j2+1});% Mesh:Case:NL name convention
    RN=RO;RN.NL=RO.NL{j2+1};NLdata=feval(RO.MeshCb.subs,'NL',mo1,RN);
    if isempty(NLdata); error('Expecting non empty NLdata');end
    if ~isfield(NLdata,'type');error('missing .type,  ''nl_inout'' is usual');end
    i1=RO.NL{j2}; if ischar(i1);i1=str2double(i1);end
    mo1=feutil(sprintf('setpro %i',i1),mo1,'NLdata',NLdata,'name',RN.NL);
  end
end

if ~isfield(mo1,'name');mo1.name=RO.name;end
if ~isfield(Range,'param');Range.param=struct;end
if ~isfield(Range,'val')%mo1=d_mesh('MeshCfg','cbi21(CoupStiff):StaticA:');
  out=mo1; 
else
 r1=struct('type','pop','value',1,'level',10,...
      'choices',{{}},'data',{{}},'SetFcn',{{@d_mesh,'stepMesh'}},...
      'ShortFmt',1, ...
      'RepList', ... % 'TestEvtData',1,
      {{'M_(.*)_([^_]+)', ... % start with M_ and end with non_ parameter name
      struct('type',{'.','{}'},'subs',{'list',{'@token{1}',3}})
      }});
 Range=sdtm.range('SafeParamInit-name',Range,'MeshCfg',r1);
 Range=sdtm.range('popMerge',Range,'MeshCfg',{mo1.name,mo1});
 if isfield(Range,'MeshCfg'); Range.MeshCfg=Range.param.MeshCfg.choices;end
 out=Range;
end

%% clean end
elseif comstr(Cam,'@');out=eval(CAM);
%% #Tuto: recover model from a specific tuto step -3
elseif comstr(Cam,'tuto'); 
 eval(sdtweb('_tuto',struct('file','d_mesh','CAM',CAM)));
 if nargout==0; clear out; end

elseif comstr(Cam,'cvs')
 out=sdtcheck('revision');
 %out='$Revision: 520 $  $Date: 2020-07-24 18:27:50 +0200 (Fri, 24 Jul 2020) $';
else; error('%s unknown',CAM);
end 
%% #End function
end
