function [out,out1,out2]=d_mesh(varargin); %#ok<*NOSEM,*STOUT>

% D_MESH sample meshes used for support
%
% Use 
%     sdtweb('_taglist','d_mesh') to view current contents
%     d_mesh('tuto')  % to see integrated tutorials
% See <a href="matlab: sdtweb _taglist d_mesh">TagList</a>


%       Copyright (c) 1990-2025 by SDTools, All Rights Reserved.
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
anglea=20; % angle to fix the actuators
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
 l1,l1-gap,anglea),mo2);
mo2.Elt=feutil('GetEdgeLine',mo2);
mo2=feutil(sprintf('Rev-OptimDegen 1 o -%.16f %.16f 0 -%.16f 0 0 1',...
 l1,l1-gap,anglea),mo2);
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
if any(Cam=='{')
 [st,RO]=sdtm.urnPar(CAM,'{div%g}:{o%s,i%s}');
else
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
 '-BC("free"#%s#"Clamping option (default : free plate)")' ...
 '-div(10#%i#"Division in X and Y")' ...
 '-PatchWidth(0.5#%g#"Relative patch width")' ...
 '-Sens("patch"#%s#"Small patch of sensors in the middle of the plate")' ...
 '-Load("botleft"#%s#Unit input at bottom left of the plate")'
 ],{RO,CAM}); Cam=lower(CAM);
end
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
 sdtw('_ewt','report EB obsolete')
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
 sdtw('_ewt','report EB obsolete')
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
elseif comstr(Cam,'track'); [CAM,Cam]=comstr(CAM,5);
 %% #BeamTrack : standard track model -----------------3
 
 [~,RO]=sdtm.urnPar(CAM,['ms(130#%g#"sleeper mass")' ...
   'kb(70e3#%g#"ballast stiffness")' ...
   'ncell(15#%g#"number of sleepers")' ...
   'unit(SI#%s#"unit system")' ...
   'lc(.15#%g#"refine length")' ...
   'pk(75e3#%g#"pad stiffness")']);
 if length(RO.ncell)<2;RO.ncell(2)=0;end
 if ~isfield(RO,'xsens'); RO.xsens=[];end

 model=struct('Node',[1  0 0 0   -.3 0 0; % first rail node
                     2  0 0 0     0 0 0;  % center rail node
                     %3  0 0 0    .3 0 0; % end rail node (not used)
                     4  0 0 0     0 0 -.07;
                     5  0 0 0     0 0 -.07-.22], ...
    'Elt',[%Inf abs('beam1');1 2 100 100 0 0;2 3 100 100 0 0;% Rail
           Inf abs('celas') 0;2 4 -3 0 200 0 0;4 5 -3 0 300 0 0;
           Inf abs('mass1') 0; 4 RO.ms*[1 1 1]  0 0 0]);
model=feutil(sprintf('RepeatSel %i .6 0 0',RO.ncell(1)),model);

% Now add rail
x=[model.Node(:,5); RO.xsens(:)]; if isfield(RO,'addMass');x=[x;vertcat(RO.addMass{:,1})];end
x=sort(unique(round(x*1000)))/1000;
x=feutil(sprintf('refineline%.15g',RO.lc),[x(1)-.3;x;x(end)+.3]);
[model.Node,i1]=feutil('addnode',model.Node,x*[1 0 0]);i1=model.Node(i1,1);
model=feutil('addelt',model,'beam1',[i1(1:end-1) i1(2:end) ones(length(i1)-1,2)*100]);
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

    error('Moved to meshblade')

 
 
elseif comstr(Cam,'mat');[CAM,Cam]=comstr(CAM,4);
%% #Mat erial database callback : mdl=d_mesh('mat',mdl,struct('mat','val','nmap',nmap));

if isempty(Cam)
 if nargin==2; model=struct; if carg==2;RO=varargin{2};carg=3;end
 else; model=RO;RO=varargin{carg};carg=carg+1;
 end
 if ischar(RO);RO=struct('mat',RO);end
 st=regexp(RO.mat,'(^[\d]+)','tokens'); 
 % 10Simo % affect MatId10
 if isfield(RO,'mpid');elseif isfield(model,'Elt')
     i1=feutil('mpid',model);RO.mpid=unique(i1,'rows');RO.mpid(RO.mpid(:,1)==0,:)=[];
 else; RO.mpid=[1 1];
 end
 if isfield(RO,'pl');
 elseif isempty(st);RO.pl=RO.mpid(1);
 else; st=st{1}{1};RO.pl=str2double(st);RO.mat=comstr(RO.mat(length(st)+1:end),1);
 end

 switch regexprep(RO.mat,'([^{,]*).*','$1') % Command before {}
 case 'SimoA'
  % #MatSimoA : sample Mooney Rivlin in large def rewritten from sdtweb dfr_ident matsimo
  r1=m_hyper('urn','SimoA{1,1,0,30,3,f5 20,g .33 .33,rho2.33n}');
  r1.il='set pro 1 isop100';
  model=sdtm.setMatPro(model,r1);
  out=model;return;
 case 'HyOpenFEM' % Hyperelastic used by OpenFEM RivlinCube
  r2=m_hyper('urn','Ref{.3,.2,.3,.3,.1,rho1u,unSI,isop100}');
 case 'HyOpenFEMd' % hyperelastic with damped cells
  r2=m_hyper('urn','Ref{.3,.2,.3,.3,.1,rho1u,unSI,g .1 .1,f 2 20,isop100}');r2.isop=100;
  % RT.nmap('MatCur')='m_hyper(urnRef{.3,.2,.3, 3k,1u,rho1u,unSI,g .1 .1,f 1 100,isop100})';

 case 'HyUP' 
  %% #MatHyUP Hyperelastic to test UP formulation
  [st,r2]=sdtm.urnPar(RO.mat,'{}{pro%s,mat%s,kappa%ug,NL%ug}');
  RO.k1=4200; RO.k2=33;RO.kappa=930000; RO.rho=1000e-6;% default values 
  RO=sdth.sfield('addmissing',r2,RO);
  RO.mu=2*(RO.k1+RO.k2); RO.lambda=RO.kappa-4/3*(RO.k1+RO.k2); 
  if abs((RO.lambda+2/3*RO.mu)/RO.kappa-1)>1e-4;error('Mismatch on Kappa');end
  RO.E=RO.mu*(3*RO.lambda+2*RO.mu)/(RO.lambda+RO.mu);
  RO.nu=RO.lambda/2/(RO.lambda+RO.mu); 
  r2.pl=[100,fe_mat('m_elastic','MM',1),RO.E,RO.nu, RO.rho, RO.E/2/(1+RO.nu)];
  model.info=RO;model.unit='MM';
  model=fe_case(model,'pcond','UP','m_hyper(''pcond 1e0'')');
  if isfield(r2,'NL')
   NLdata=struct('type','nl_inout','opt',[0,0,0], ...
     'MexCb',{{m_hyper('@hyper_g'),struct}},'adofi',-.99*ones(29,1));
   if isKey(RO.nmap,'g1')
     r1=RO.nmap('g1');
     NLdata.g=r1(:,1); NLdata.f=r1(:,2); NLdata
   end
   model=feutil('setpro 1 isop100',model,'NLdata',NLdata);
  end
  
 case 'DamA'
  %% #MatDamA : test case for damage testing
  r2=struct; r2.NLdata=struct('type','nl_inout','opt',zeros(1,3), ...
     'iopt',int32([0 0 6 2]),'adofi',[-.96;-.97], ...
     'MexCb',{{nlutil('@dama_g'),[]}},'wy',1e3,'gamma',1,'snl',[],'StoreType',3);
  r2.pl=m_elastic('dbval 1 steel -unit TM');r2.unit='TM';
 case 'DamB'
  %% #MatDamB : test case for surface damage testing
  r2=struct; r2.NLdata=struct('type','nl_inout','opt',zeros(1,3), ...
     'iopt',int32([0 0 6 2]),'adofi',[-.96;-.97], ...
     'MexCb',{{nlutil('@dama_g'),[]}},'wn',1e3,'gn',1,'wt',1e3,'gt',1,'snl',[],'StoreType',3);
  % Stiff 3rd direction enforces sliding. 1 & 2 give sliding stiffness
  r2.pl=[1 fe_mat('m_elastic','TM',4) 2e3 0 3e3 0 0 1e4 1e-10];r2.unit='TM';

 case 'PadA'
  %% #MatPadA R. Zhuravlev, table 2.3 PhD Thesis, ENSAM, 2017 https://pastel.archives-ouvertes.fr/tel-01744302
  r2=m_hyper('urn','PadA{2.5264,-0.9177,0.4711,1200,3,f .35,g .5688,rho1n,tyYeoh,unTM}');
  r2.isop=100; 
 otherwise; 
   if isfield(RO,'nmap')&&isKey(RO.nmap,RO.mat)
     r2=RO.nmap(RO.mat);
     if ischar(r2);sdtm.addLog(RO.nmap,r2); 
         r2=sdtm.urnCb(r2);
         r2=feval(r2{:});
     end
   else; r2=RO.mat;
   end
   if ischar(r2)
    %% attempt to use URN such as m_hyper(urnSimoA{1,1,0,30 ...
    r2=sdtm.urnCb(r2);r2=feval(r2{:}); % Attempt at using URN callback
    %catch
    % error('Mat%s not implemented',RO.mat)
    %end
   end
 end
 r2.pl(1)=RO.pl(1); 
 if isfield(model,'Elt')
  %% standard affect to model 
  r2.mpid=RO.mpid; out=sdtm.setMatPro(model,r2);
 else;out=r2; %r2=d_mesh('mat','PadA') % simply return structure
 end 
elseif comstr(Cam,'rve');[CAM,Cam]=comstr(CAM,4);
 % #MatRve -2
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
elseif comstr(Cam,'case'); [CAM,Cam]=comstr(CAM,5);
 %% #Case : step case (called by sdtsys StepMesh)
 sdtu.logger.entry(struct('message',sprintf('d_mesh(''Case%s'')',CAM),'type','diag'));

 if isfield(RO,'Elt');model=RO;end
 RO=varargin{carg};carg=carg+1;

 [CAM,Cam]=comstr(RO.Case,1);
 if comstr(Cam,'compa')
  %% #CaseCompA
  % Use lsutil to get material orientation map
  %RO.nmap('InitMatOrient')
  % leading / trailing edge
  
 RL=struct('OrientLine',struct('starts',1,'dir',[0 0 1]), ...
      'SurfSel','SelFace & facing >.7 0 -1e5 0','surforient',10, ...
      'Face1',[1 2 5]);
 eltid=feutil('eltidfix;',model);
 r1=RO.CaseVal(strncmpi(RO.CaseVal,'theta',5));
 if ~isempty(r1);r1=str2double(r1{1}(6:end));else; r1=30;end
 RL.EltIdAng=eltid(eltid~=0);RL.EltIdAng(:,2)=r1; % as two columns
 out=fe_shapeoptim('MapOrient',model,RL);


 elseif comstr(Cam,'parvecut')
  %% #CaseParVeCut (naca with matched layers)
  % add parametrization of visc layer with cut

  % xxx should be done elsewhere
  % cleanup materials and set a visc property
  
  %pl=double(model.nmap.('Map:MatName').('visc'));
  %pl=[pl fe_mat('type','m_elastic','MM',1) 1e3 .45 1.190e-6 0];
  %model.pl(model.pl(:,1)==pl(1),1:length(pl))=pl;
  %model.pl(:,6)=0; % xxx G

  % add bc
  model=fe_case(model,'fixdof','footclamp','z==0');

  % prepare matcut: do a selection of elements with order to follow
  % need to match cst and visc layers to process
  n1=feutil('findnode inelt{matname cply} & inelt{matname visc}',model);
  el1=feutil('selelt matname cply & selface & innode{nodeid}',model,n1);
  el2=feutil('selelt matname visc & selface & innode{nodeid}',model,n1);
  % quad nodes: match
  in1=find(isfinite(el1(:,1)));
  [i1,i2]=ismember(sort(el2(isfinite(el2(:,1)),1:4),2),...
   sort(el1(isfinite(el1(:,1)),1:4),2),'rows');
  % matched elements
  r1=[el2(in1(i1),7) el1(in1(i2(i1)),7)];
  % add visc w/o cply
  el1=feutil('selelt matname visc',model); eid=feutil('eltid',el1); eid(eid==0)=[];
  eid=setdiff(eid,r1(:,1)); eid(:,2)=NaN;
  r1=[r1;eid];
  % there should be not cply left
  el1=feutil('selelt matname cply',model); eid=feutil('eltid',el1); eid(eid==0)=[];
  if ~isempty(setdiff(eid,r1(:,2))); sdtw('_nb','some constraind layer elements missed'); end

  % define cuti into matrix
  model=fe_caseg('ParMatCut',model,'matname visc | matname cply',...
   struct('lab','CstVisc','EltId',r1,'val',1));

  out=model;

 elseif comstr(Cam,'empty'); out=model;return
 else; error('Case%s unknown',CAM)
 end


elseif comstr(Cam,'mesh'); [CAM,Cam]=comstr(CAM,5);
 %% if not found edit MeshCmd to Cmd
%% #Mesh : MeshCfg mesh command implementations (called by sdtsys StepMesh)
sdtu.logger.entry(struct('message',sprintf('d_mesh(''Mesh%s'')',CAM),'type','diag'));

if comstr(Cam,'plate')
%% #MeshPlate {div,o%s,i%s}
 [st,RO]=sdtm.urnPar(CAM,'{div%g}:{o%s,i%s}');
% [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
%  '-BC("free"#%s#"Clamping option (default : free plate)")' ...
%  '-div(10#%i#"Division in X and Y")' ...
%  '-PatchWidth(0.5#%g#"Relative patch width")' ...
%  '-Sens("patch"#%s#"Small patch of sensors in the middle of the plate")' ...
%  '-Load("botleft"#%s#Unit input at bottom left of the plate")'
%  ],{RO,CAM}); Cam=lower(CAM);
 %% #PlateGenMesh : Mesh with varying number of elts + clamp strat---------3
 node=[0 0 0;0 1 0;1 1 0;1 0 0];
 model=feutil('ObjectQuad 1 1',node,RO.div,RO.div);
 
 dof=[]; % List of fixed dofs
 if ~isfield(RO,'BC')||strcmpi(RO.BC,'free'); dof=[]; % Free plate
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

%elseif comstr(Cam,'sens'); [CAM,Cam]=comstr(CAM,5);
 %% Place sensors (right side, ...)-
 if strncmpi(RO.o,'patch',5); % patch of sensors in the middle of the plate
  RO.PatchWidth=str2double(RO.o(6:end));
  i1=[.5-RO.PatchWidth/2 .5+RO.PatchWidth/2];
  n1=feutil('findnode',model,sprintf('x>=%g&x<=%g&y>=%g&y<=%g',i1,i1));
  tdof=feutil('getdof',n1,.03); % Sensors only in z direction
 else; error('unknown Sens command');
 end
 % Place sensors over the model
 model=fe_case(model,'SensDof append trans','Test',tdof);

 %% Load : Place Load (first quarter of beam,...)-----------------3
 if ~isfield(RO,'i')||strcmpi(RO.i,'bl');
  % Add Load at the middle left of the beam
  i1=feutilb('AddNode -nearest',[0.1 0.1 0],model.Node(:,5:7));
  i1=model.Node(find(i1),1);
  model=fe_case(model,'DofLoad','In',struct('def',1,'DOF',i1+.03,'curve',{{'Input'}}));
 else; error('unknown Load command');
 end
out=model;
elseif comstr(Cam,'blade')
%% #MeshBlade

if any(Cam=='{');
  [st,RO]=sdtm.urnPar(CAM,'{}:{div%g,mat%s,quad%31}');
 if isempty(RO.div);RO.div=[12 1 5];end
 RO.d1=RO.div(1);RO.d2=RO.div(2);RO.d3=RO.div(3);
else
 [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
  'd1(12#%i#"divisions along length")' ...
  'd2(1#%i#"divisions along tickness")' ...
  'd3(5#%i#"divisions along width")' ...
  'mat(iso#%s#"material variants")' ...
  ],{RO,CAM}); Cam=lower(CAM);
end
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
 %% orient illustrated in d_mesh CaseCompA

 % data=struct('EltId',eltid(cEGI),'bas',eltid(cEGI));
 % NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
 % for jElt=1:length(cEGI)
 %  n1=model.Node(NNode(model.Elt(cEGI(jElt),1:3)),:);
 %  p=diff(n1(:,5:7));p=sp_util('basis',p(1,:),p(2,:));
 %  data.bas(jElt,7:15)=p(:)';
 % end
 % model=stack_set(model,'info','EltOrient',data);
end

model=fe_case(model,'fixdof','base','r<2'); model.name='blade';

if nargout==0; feplot(model);fecom showpatch
else; out=model;
end 

elseif comstr(Cam,'nacaskin')
 %% #MeshNacaSkin : non generic example

if ~isfield(RO,'projM');RO.projM=RO.nmap;end
model=RO.projM('CurModel'); % EdgeN contains Thickness


%% mesh area defining ply extent
% show level cuts for target thickness
cf=feplot(2,';');
fevisco('CompThick{PlyLevels 1 4.4 5}',model)
% xxx possible extension mesh using contours
go=findobj(2,'tag','iso');
[i1,i2]=unique([go.Faces(:,[1 4]);go.Faces(:,[2 3])],'rows');
n1=(go.Vertices(i1(:,1),:)+go.Vertices(i1(:,2),:))/2;
i1(:,3)=round(sqrt(sum((go.Vertices(i1(:,1),:)-go.Vertices(i1(:,2),:)).^2,2)));
mo1=struct('Node',[i1(:,1)*[1 0 0 0] n1],'Elt',feutil('addelt','beam1',go.Faces(:,[1 2])));
mo1=sdtu.fe.selCoarsenF2(mo1,struct('ToFace',[1 15 .7]));mo1.Node=feutil('getnodegroupall',mo1);

i1=fe_gmsh('lineloops',mo1.Elt);NNode=sdtu.fe.NNode(mo1);
for j1=1:length(i1)
 r2=diff(mo1.Node(NNode(i1{j1}([1 end])),5:7));
 if r2(3)<-10;  i1{j1}=fliplr(i1{j1});r2=-r2;end
 if r2(3)<1&&r2(1)<0; i1{j1}=fliplr(i1{j1});r2=-r2;end % 
 r3(j1)=mo1.Node(NNode(i1{j1}(1)),5);
end
[~,i2]=sort(r3);i1=i1(i2);

% create reference lines and mesh plate
i2=fliplr([i1{1} fliplr(i1{4}) fliplr(i1{2})])';mo1=fe_gmsh('AddLine -loop1',mo1,[i2 i2([2:end 1])]);
i2=fliplr([i1{2} fliplr(i1{3})])';mo1=fe_gmsh('AddLine -loop2',mo1,[i2 i2([2:end 1])]);
i2=[fliplr(i1{3})]';mo1=fe_gmsh('AddLine -loop3',mo1,[i2 i2([2:end 1])]);
mo1.Stack{end}.PlaneSurface=[1;2;3]; 
% xxx GV 
mo1.Stack{end}.Mesh=struct('Algorithm',5, ...
    'RecombinationAlgorithm',3,'RecombineAll',1,'MeshSizeMin',2,'MeshSizeMax',7);
S=scatteredInterpolant(mo1.Node(:,5),mo1.Node(:,7),mo1.Node(:,6));S.Method='natural';
mo1.Node(:,6)=0;
area=fe_gmsh('write @tempdir/tmp.geo -lc 5 -run -2 -v 0 -vol2',mo1);
area.Node(:,6)=S(area.Node(:,5),area.Node(:,7));

mpid=feutil('mpid',area);
mpid(mpid(:,2)==138,1:2)=1; % 1 layer
mpid(mpid(:,2)==139,1:2)=2; % 2 layer
mpid(mpid(:,2)==140,1:2)=3; % 3 layer
area.Elt=feutil('mpid',area,mpid);
RO.nmap('area')=area;

%% Mesh LongSkin 

skin=RO.projM('CurModel');
skin.Elt=feutil('selelt selface &facing <.8 0 0 1000',skin);skin=feutil('divide 3 3',skin);

% round off the skin tip
n1=feutil('getnode z==',skin,max(skin.Node(:,7)));
n2=feutil('getnode inelt{InNode&seledge}',skin,n1(:,1));
n1(ismember(n1(:,1),n2(:,1)),:)=[];
n1(:,7)=n1(:,7)+3;fecom('shownodemark',n1(:,5:7))
skin.Node(sdtu.fe.NNode(skin,n1(:,1)),:)=n1;

%cf.sel(1)='reset';
%cf.sel(2)='urn.LineTopo{starts 1 1520,cos.6}';cf.o(2)='sel2 ty1';
%cf.sel(2)='urn.LineTopo{starts 2 1503,cos.6}';cf.o(2)='sel2 ty1';

i1=feutil('geolinetopo',skin,struct('starts',[1 1520;2 1503],'cos',.6));

mo5=skin; mo5.Elt=feutil('selelt seledge & innode',skin,i1{2});
mo5=feutil('rev 5 o 0 10 0  90 1 0 0',mo5);
skin.Node=mo5.Node;skin.Elt=feutil('addelt',skin.Elt,mo5.Elt);

mo5=skin; mo5.Elt=feutil('selelt seledge & innode',skin,i1{1});
mo5=feutil('rev 5 o 0 -10 0  90 -1 0 0',mo5);
skin.Node=mo5.Node;skin.Elt=feutil('addelt',skin.Elt,mo5.Elt);
%feplot(skin,'showfipro')
skin.name='Skin with foot transition';

RO.projM('skin')=skin;


elseif comstr(Cam,'naca')
 %% #MeshNaca : sample blade with NACA profile and orthotropic material

 if carg>nargin
  RO=struct('nmap',vhandle.nmap);
 else
  RO=varargin{carg}; carg=carg+1;
  if isfield(RO,'urn');CAM=RO.name;end
 end
 if comstr(Cam,'naca{')
  r1=sdth.findobj('_sub{}',CAM);
  r2=sdth.findobj('_sub:',r1(2).subs{1});
  %'xxx interpret further parameters'
  list=RO.nmap.('Map:Bprofiles').(r2(1).subs);
  if length(r2)>1;RM=RO.nmap.('Map:Bplies').(r2(2).subs);
   %RO.plyList=RM.plyList; RO.OrientLine=RM.OrientLine;
   RO=sdth.sfield('AddMissing',RO,RM);
  else;RO.ver='box';
  end
 else % xxx
  if ~isKey(RO.nmap,'Map:profiles') % Default blade
   RO.nmap=d_mesh('nmap');
   RO.nmap('Profile')=RO.nmap.('Map:Bprofiles').('HyFoilA');
  end
  list=RO.nmap('Profile'); % recover profile
 end
 [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
  'xn(20#%i#"refinement in long side")' ...
  'yn(2#%i#"refinement in thickness")' ...
  'zs(.1#%g#"refinement step length in transverse")' ... % xxx ug
  'sRef(400#%i#"profile curve refinement factor")' ...
  'extr(5#%i#extrude in length")' ...
  'unit(mm#%s#"unit system")' ...
  'interp(#3#"interpolate sections along radius")' ...
  ],{RO,CAM}); Cam=lower(CAM);

 if isfield(RO,'ver')&&strcmpi(RO.ver,'box') 
 %% #MeshNaca.box create blade box volume from section list
 RO.iSec=4;

 for jl=2:size(list,1) % loop on profile along radius
  if ischar(list{jl,3}); % resolve section here
   RO.curProf=RO.nmap.('Map:Bsections').(list{jl,3});
   RO=buildContour(RO);% xxx should call splineR  
  else; RO.contour=list{jl,3}; % profile directly provided as nodal contour
  end
  mo1=sdtu.fe.boxInFlatContour(RO); % 2D profile

  % skewness  x,y,ry,scale sdtu.fe.splineR
  %% offset
  mo1.Node(:,5)=mo1.Node(:,5)+list{jl,2};
  mo1.Node(:,7)=mo1.Node(:,7)+list{jl,1};
  if isfield(RO.curProf,'skew')
   if RO.ilim>4; no=list{jl,6}; an=list{jl,7};
   else;   no=RO.curProf.skew.Orig; an=RO.curProf.skew.angle;
   end
   n2=mo1.Node(:,5:7) -no; %   n2=n2-no;
   n2=n2*[cos(an*pi/180) sin(an*pi/180) 0;    -sin(an*pi/180) cos(an*pi/180) 0;    0 0 1] +no;%   n2=n2+no;
   mo1.Node(:,5:7)=n2;
  end
  list{jl,RO.iSec}=mo1; 

 end % loop on section list
 model=sdtu.fe.boxCutsToVol(list(2:end,RO.iSec),RO);
 model=stack_set(model,'info','MeshOpt',RO);
 model.unit=RO.unit(1:2);
 RO.nmap('CurModel')=model; 
 d_mesh('MeshNacaSkinArea',RO)
 out=model;
 return

 elseif RO.interp; %isfield(RO,'interp')
  % #MeshNaca.Interp goal is to prepare appealing demonstration blade profiles
  % extract sections, interpolate extents, refine section positions

  for jl=1:size(list,1) % extract section info: create refined section profile for each entry
   RO.curProf=RO.nmap.('Map:Bsections').(list{jl,3});
   %% Generate profile curve (angle, xpos ypos), provide
   r1 = spline(RO.curProf.section(:,1)',...
    [RO.curProf.EndSlope;RO.curProf.section(:,2:3);RO.curProf.EndSlope]');
   s1=linspace(0,2*pi,RO.sRef*2*round(RO.xn/2)); % angle refinement
   r2=ppval(r1,s1); % interpolated curve
   %figure(11); plot(r2(1,:),r2(2,:),'-b',RO.curProf.section(:,2),RO.curProf.section(:,3),'or'),axis equal
   n0=[(1:size(r2,2)-1)' zeros(size(r2,2)-1,3) r2(:,1:end-1)' zeros(size(r2,2)-1,1)];
   list{jl,3}=n0;
  end

  % recover reasonable radius div
  n0=cell2mat(list(:,3));
  rd=round(max(abs(diff(sort(cell2mat(list(:,1))))))/...
   (max(n0(:,5))-min(n0(:,5)))*RO.xn);
  % create a new interpolated list based on provided sections
  rdv=ceil((list{end,1}-list{1})/rd);

  %% interpolate xz contour along R w/ cubic
  % prepare xz contour
  r1=list(:,1:3); r1(:,3)=cellfun(@(x)max(x(:,5)),r1(:,3),'uni',0);
  r1=cell2mat(r1);
  r2=[[r1(:,1);flipud(r1(1:end-1,1))] [r1(:,2);flipud(r1(1:end-1,3))]];
  % figure(3);  plot(r2(:,2),r2(:,1))

  % now prepare angular interpolation of xz contour
  r3=mean(r2);
  r22=r2-r3;  r23=sqrt(sum(r22.^2,2));  r24=angle(r22(:,1)+1i*r22(:,2));

  r10=linspace(r24(1),r24(end),200*rdv);
  %pp=spline(r24,r23);  %r11=ppval(pp,r10)
  r11=interp1(r24,r23,r10,'cubic');
  r4=[real(r11.*exp(1i*r10))' imag(r11.*exp(1i*r10))']+r3;
  %clf(33);figure(33); plot(r2(:,2),r2(:,1),'o--'); hold on;figure(33);  plot(r4(:,2),r4(:,1))

  %% create new section list from interpolated contour
  % setup section positions
  rlist=linspace(list{1},list{end,1},ceil((list{end,1}-list{1})/rd))'; rlist(end,3)=0;
  [u1,i1]=max(r4(:,1));
  r41=r4(1:i1,:); r42=r4(i1+1:end,:);
  for j1=1:size(rlist,1)
   [r5,i5]=min(abs(r41(:,1)-rlist(j1))); rlist(j1,2)=r41(i5,2);
   [r5,i5]=min(abs(r42(:,1)-rlist(j1))); rlist(j1,3)=r42(i5,2);
  end
  %clf(33);figure(33); plot(r2(:,2),r2(:,1),'o--'); hold on;figure(33);
  %plot([rlist(:,2);flipud(rlist(1:end-1,3))], [rlist(:,1);flipud(rlist(1:end-1,1))]);

  %% interpolate skew angle XXX origin
  rlist(:,7)=interp1(cell2mat(list(:,1)),cell2mat(list(:,7)),rlist(:,1),'makima');
  ori=cell2mat(list(:,6));
  ori=interp1(cell2mat(list(:,1)),ori,rlist(:,1));

  %% interpolate thickness
  rlist(:,5)=interp1(cell2mat(list(:,1)),cell2mat(list(:,5)),rlist(:,1),'makima');

  %% fill the new rlist
  % propagate section length
  rlist(:,4)=rlist(:,3);
  rlist=num2cell(rlist);

  for j1=1:size(rlist,1)
   rlist{j1,6}=ori(j1,:); % place interpolated origin
   % find closest section profile and scale it to relatice extent and thickness
   [u1,i1]=min(abs(cell2mat(list(:,1))-rlist{j1,1}));
   n0=list{i1,3}; % closest profile
   n0(:,5)=n0(:,5)*(rlist{j1,4}-rlist{j1,2})/(list{i1,4}-list{i1,2}); % scale xextent
   n0(:,6)=n0(:,6)*rlist{j1,5}/list{i1,5}; % scale thickness (yextent)
   rlist{j1,3}=n0; % store
  end
  list=rlist;
 end % end intepr

 %% create blade volume from section list
 RO.EdgeN=[]; RO.TipN=[]; RO.ExtN=[]; % nodes on extra/intra with vertical lines for thickness
 RO.ilim=size(list,2)+1; RO.ilie3=size(list,2)+2;
 for jl=1:size(list,1) % loop on profile along radius
  if ischar(list{jl,3}); % resolve section here
   RO.curProf=RO.nmap.('Map:Bsections').(list{jl,3});
   % Generate profile curve (angle, xpos ypos), provide
   r1 = spline(RO.curProf.section(:,1)',...
    [RO.curProf.EndSlope;RO.curProf.section(:,2:3);RO.curProf.EndSlope]');
   s1=linspace(0,2*pi,RO.sRef*2*round(RO.xn/2)); % angle refinement
   r2=ppval(r1,s1); % interpolated curve
   %figure(11); plot(r2(1,:),r2(2,:),'-b',RO.curProf.section(:,2),RO.curProf.section(:,3),'or'),axis equal

   % prepare contour model, with closed line
   n0=[(1:size(r2,2)-1)' zeros(size(r2,2)-1,3) r2(:,1:end-1)' zeros(size(r2,2)-1,1)];

  else; n0=list{jl,3}; % profile directly provided as nodal contour
  end

  % prepare a grid to morph
  %x1=linspace(min(RO.curProf.section(:,2)),max(RO.curProf.section(:,2)),RO.xn); x1=x1(:); % xxx should be max of profile
  x1=linspace(min(n0(:,5)),max(n0(:,5)),RO.xn); x1=x1(:); % xxx should be max of profile
  y1=[min(r2(2,:)) max(r2(2,:))];
  node=[x1(1) y1(1) 0;x1(end) y1(1) 0;x1(end) y1(2) 0;x1(1) y1(2) 0];
  mo1=feutil('objectquad 1 1',node,RO.xn,1);

  if jl>1
   mo1=feutil('renumber-noori',mo1,max(list{jl-1,RO.ilim}.Node(:,1)));
  end

  NNode=sparse(mo1.Node(:,1),1,1:size(mo1.Node,1));
  % prepare vertical top to bottom edges info
  el1=feutil('selelt seledge',mo1); r1=sort(el1(isfinite(el1(:,1)),1:2),2);
  el2=feutil('selelt seledgeAll',mo1); r2=sort(el2(isfinite(el2(:,1)),1:2),2);
  el2=setdiff(r2,r1,'rows');
  % then sort top to bottom
  [~,i2]=sort(reshape(mo1.Node(NNode(el2),6),[],2),2,'descend'); i2=i2(:,1)~=1;
  el2(i2,:)=el2(i2,[2 1]);

  % select edge, identify closest node on curve and move
  % Surfstick xxxeb
  mo2=mo1; mo2.Elt=el1; mo2.Node=feutil('getnodegroupall',mo2);
  %[n2,i2]=feutil('addnode-nearest epsl10',n0,mo2.Node(:,5:7));
  % Keep verticals
  % look for closest abscissa and choose point with minimal height correction
  n2=mo2.Node;
  for j1=1:size(n2,1)
   r2=abs(n2(j1,5:6)-n0(:,5:6)); [r3,in2]=sort(abs(r2(:,1)),'ascend');
   [r4,in3]=min(r2(in2(1:10),2));
   n2(j1,5:6)=n0(in2(in3),5:6);
  end

  NNode=sparse(mo1.Node(:,1),1,1:size(mo1.Node,1));
  %mo1.Node(NNode(mo2.Node(:,1)),5:7)=n2(i2,5:7);
  mo1.Node(NNode(mo2.Node(:,1)),5:7)=n2(:,5:7);

  % cleanup: coalesce exit tip
  %[r1,i1]=sort(abs(mo1.Node(:,5)-max(RO.curProf.section(:,2))),'ascend');
  [r1,i1]=sort(abs(mo1.Node(:,5)-max(n0(:,5))),'ascend');
  n1=mo1.Node(:,[1 1]); n1(i1(2),2)=n1(i1(1),1);
  RO.ExtN=[RO.ExtN;n1(i1(1),1)];
  % check height of front tip
  %[r1,i1]=sort(abs(mo1.Node(:,5)-min(RO.curProf.section(:,2))),'ascend');
  [r1,i1]=sort(abs(mo1.Node(:,5)-min(n0(:,5))),'ascend');
  %if abs(diff(mo1.Node(i1(1:2),6)))<abs(diff(y1))/20 % xxx should always be the case
  n1(i1(2),2)=n1(i1(1),1);
  RO.TipN=[RO.TipN;n1(i1(1),1)];
  %end
  mo1=feutil('renumber-noori',mo1,n1);
  n1=sparse(n1(:,1),1,n1(:,2)); el2=full(n1(el2));

  % EdgeN [bottom up] for verticals
  RO.EdgeN=[RO.EdgeN;el2];

  % we should divide after extrusion to allow profile variations
  % the grid should remain the same, extrusion should connect parts
  % use divide to clean quality if needed

  % offset
  mo1.Node(:,5)=mo1.Node(:,5)+list{jl,2};
  mo1.Node(:,7)=mo1.Node(:,7)+list{jl,1};

  % skewness
  if isfield(RO.curProf,'skew')||RO.ilim>4
   if RO.ilim>4; no=list{jl,6}; an=list{jl,7};
   else;   no=RO.curProf.skew.Orig; an=RO.curProf.skew.angle;
   end
   n2=mo1.Node(:,5:7) -no; %   n2=n2-no;
   n2=n2*[cos(an*pi/180) sin(an*pi/180) 0;    -sin(an*pi/180) cos(an*pi/180) 0;    0 0 1] +no;%   n2=n2+no;
   mo1.Node(:,5:7)=n2;
  end

  % also create a mid plane section
  mo2=mo1;
  n2=feutil('getnodegroupall',mo2); mo2.Node=n2;
  mo2=feutil('divideelt 2 1',mo2);
  n3=feutil('getnodegroupall',mo2);
  [~,i3]=setdiff(n3(:,1),n2(:,1)); NN=sparse(n3(:,1),1,1:size(n3,1));
  n3=n3(i3,:); [~,i3]=sort(n3(:,1),'ascend');
  r3=[RO.TipN(end);n3(i3,1);RO.ExtN(end)];
  n3=feutil('getnode',mo2,r3);
  el3=[r3(1:end-1,1) r3(2:end,1)]; % middle section

  % generate base section and mid plane section
  if jl==1
   model=struct('Node',mo1.Node,'Elt',[]);
   mo3=struct('Node',n3,'Elt',[]);
  else;
   %mo1=feutil('renumber-noori',mo1,max(list{jl-1,4}.Node(:,1))); % earlier
   model.Node=[model.Node;mo1.Node];
   el1=[list{jl-1,RO.ilim}.Elt(:,1:4) mo1.Elt(:,:)];
   el1(~isfinite(el1(:,1)),:)=[];
   model.Elt=feutil('addelt',model.Elt,'hexa8',el1);

   [mo3.Node,i3]=feutil('AddNodeKnownNew',mo3.Node,n3(:,5:7));
   NN=sparse(n3(:,1),1,mo3.Node(i3,1));
   el3=full(NN(el3));
   el4=[list{jl-1,RO.ilie3} fliplr(el3)];
   mo3.Elt=feutil('addelt',mo3.Elt,'quad4',el4);

  end% base section addition

  list{jl,RO.ilim}=mo1; list{jl,RO.ilie3}=el3;

 end % loop on section list

 eltip=feutil('selelt seledgeAll & innode{nodeid}',model,RO.TipN);
 elted=feutil('selelt seledgeAll & innode{nodeid}',model,RO.ExtN);

 model.Elt=feutil('JoinAll',model.Elt); % cleanup

 if RO.extr % link profile sequence with an adequate mesh length
  rd=round(max(abs(diff(sort(cell2mat(list(:,1))))))/...
   (max(model.Node(:,5))-min(model.Node(:,5)))*RO.xn);
  mo1=model;
  if ~RO.interp
   model=feutil(sprintf('divideelt 1 1 %i',rd),model);
  end
  %mo3=feutil(sprintf('divideelt %i 1',rd),mo3);
  % now identify intra/extra by identifying equivalent refinement on faces from EdgeN
  %RO.midPlane=mo3;

  % store tip line
  mo2=struct('Node',mo1.Node,'Elt',eltip); mo2.Node=feutil('getnodegroupall',mo2);
  if ~RO.interp
  mo2=feutil(sprintf('divideelt %i',rd),mo2);mo2.Node=feutil('getnodegroupall',mo2);
  end
  [n2,i2]=feutil('AddNode-nearest',model.Node,mo2.Node(:,5:7));
  model=feutil('AddSetNodeId',model,'TipLine',n2(i2,1));
  % store edge line
  mo2=struct('Node',mo1.Node,'Elt',elted); mo2.Node=feutil('getnodegroupall',mo2);
  if ~RO.interp
  mo2=feutil(sprintf('divideelt %i',rd),mo2);mo2.Node=feutil('getnodegroupall',mo2);
  end
  [n2,i2]=feutil('AddNode-nearest',model.Node,mo2.Node(:,5:7));
  model=feutil('AddSetNodeId',model,'EdgeLine',n2(i2,1));

  %% measure thickness and store
  in1=isfinite(model.Elt(:,1));
  r5=[reshape(model.Elt(in1,[3 4 7 8]),[],1) reshape(model.Elt(in1,[2 1 6 5]),[],1)];
  r5=unique(r5,'rows');
  RO.EdgeNN=r5;
  % generate the thickness map from EdgeNN and add centers here
  NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
  n61=model.Node(NNode(RO.EdgeNN(:,1)),5:7);
  n62=model.Node(NNode(RO.EdgeNN(:,2)),5:7);
  r6=.5*sqrt(sum((n62-n61).^2,2));
  RO.EdgeNN(:,3)=r6; % thickness from "mid" as half total
  RO.EdgeNN(:,4)=n61(:,3); % rad, could be mean of n61 and n62

  model=feutil('divideelt 1 2',model); % create midPlane in mesh
  [~,model.Elt]=feutil('eltidfix;',model);

  [n5,i5]=feutil('addnode-nearest',model.Node,.5*(n61+n62));
  def=struct('DOF',[r5(:,1)+.98;r5(:,2)+.98;n5(i5,1)+.98],'def',[0*r6;0*r6;r6]);
  [~,idof]=unique(round(def.DOF*100));
  def.DOF=def.DOF(idof); def.def=def.def(idof,:);

  RO.PThick=def;
  model=feutil('addsetnodeid',model,'ExtraDos_n',unique(RO.EdgeNN(:,1)));
  model=feutil('addsetnodeid',model,'IntraDos_n',unique(RO.EdgeNN(:,2)));
  model=feutil('addsetnodeid',model,'midPlane_n',unique(n5(i5,1)));

  model=feutil('addsetfaceid',model,'extrados_face','selface & innode{setname ExtraDos_n}');
  model=feutil('addsetfaceid',model,'intrados_face','selface & innode{setname IntraDos_n}');

  % midPlane: do face selection from extrados and midPlane nodes, remove degen faces
  mo1=model; mo1.Elt=feutil('selelt setname extrados_face:underlying & selface',mo1);
  % to remove degen faces, use temporatry unique eltid
  el0=mo1.Elt; elid0=feutil('eltid',el0);
  [eltid,mo1.Elt]=feutil('eltidfix;',mo1);
  data=fe_quality('MeasDegen-silent2',mo1);
  for j1=1:size(data.data,2); data.EltId{j1}=data.EltId{j1}(logical(data.data{j1})); end
  data=cat(1,data.EltId{:}); % degen elts on temporaty eltid
  EEid=sparse(eltid+1,1,1:length(eltid)); % convert to eltind
  mo1.Elt=el0; % get back to initial elts as represent faceselections
  mo1.Elt=feutil('removeelt eltind',mo1,full(EEid(data+1))); % remove found eltind
  mo1.Elt=feutil('selelt innode{nodeid}',mo1,n5(i5,1)); % keep midPlane
  model=feutil('addsetfaceid',model,'midPlane',mo1.Elt); % add face from selface

  %model=feutil('addsetfaceid',model,'midPlane',...
  % sprintf('setname extrados_face:underlying & selface & innode{nodeid %s}',num2str(n5(i5,1)')));

 end

 % refine in thickness xxx
 if RO.yn>1; mo1=feutil(sprintf('divideelt %i 1 ',RO.yn),mo1); end

 if isfield(RO,'foot')&&RO.foot
  % hub diameter and thickness
  % bladeori

  % generate cylinder mesh, match surf base foot
  % extrude foot and apply sitcknode shape

  mo1=model;
  mo1.Elt=feutil('selelt selface & innode{z==min(z)}',model);

  % extrude over transition
  mo1.Node=feutil('getnodegroupall',mo1); n1=mo1.Node;
  mo1=feutil('extrude 1 0 0 -3',mo1); % xxx foot transition height

  % place new nodes on cylinder face
  n2=mo1.Node(~ismember(mo1.Node(:,1),n1(:,1)),:); n2i=n2(:,1);
  %n2(:,6)=1.5*n2(:,6);
  n2(:,6)=n2(:,6)+2*sign(n2(:,6));


 % mof=feutil('objectannulus 50 0 -200 192 197  0.707 0.707 0 36 1');
  %mof=feutil('extrude 0 0.707 0.707 0',mof,[-60:10:60]);

 % match=struct('Node',n2(:,5:7));
  %match=feutilb('matchsurf-radiusinf',mof,match,'selface');
  NNode=sparse(mo1.Node(:,1),1,1:size(mo1.Node,1));
  mo1.Node(NNode(n2i),5:7)=n2(:,5:7); %match.StickNode;

%  feplot(feutilb('combinemodel',mo1,mof))

% prepare arc
% extrude and match foot nodes + extend wit "fillet" 
% remove nodes inside foot contour, do addnode -noearest to make compatible


na1=197*[0 sin(20*pi/180) cos(20*pi/180);0 sin(-20*pi/180) cos(-20*pi/180)]
na1=(na1+[50 0 -200])

b1=basis('rotate',[1 1 0  0 0 0  1 0 0  0 1 0  0 0 1],'rz=-45',1);
na1=(na1-[50 0 0])*reshape(b1(7:15),3,3)+[50 0 0];

moa=feutil(sprintf('objectarc 50 0 -200 %.15g %.15g %.15g %.15g %.15g %.15g 36 1',reshape(na1',1,[])))
moa=feutil('extrude 0 1 0 0',moa,[-100:10:120]);


  feplot(feutilb('combinemodel',mo1,moa))

% then remove nodes of moa inside contour
% then do addnode-nearest to close surface
% then cut edges and place edges on yz plane









 end


 model=stack_set(model,'info','MeshOpt',RO);
 model.unit='MM';% xxx
 if isfield(RO,'plyList'); 
  mo1=model;
  model=fevisco('MeshPlies',model,RO);
  model=stack_set(model,'info','OrigVolume',mo1);
 end
 out=model;



elseif comstr(Cam,'cfg'); [CAM,Cam]=comstr(CAM,4);
%% #MeshCfg : obtain/generate mesh  ---------------------------------------
error('Moved to sdtsys StepMesh')

else; % Redirect to Cmd
        eval(iigui({'d_mesh(CAM,varargin{2:end})',nargout},'OutReDir'));
end

elseif comstr(Cam,'nmap');
%% #nmap named maps ---

[key,nmap,uo,carg]=sdtm.stdNmapArgs(varargin,CAM,carg);RO.inKeys=nmap.keys;
projM=nmap;

%% #Cin: #Map:Cin parameter formating and tooltip for meshing ----2

cinM=projM('Map:Cin');% Create a default CinCell structures to be used for vhandle.uo 
cinM.add={
'gr:RefineHexaMesh','Hexa mesh refinement in fe_shapeotim', ...
   {[ ... %% RefineHexaMesh
  '-divlc(#%g#"to give a target lc for uniform division and mpc coupling at once")' ...
  '-InterMPC(#3#"to add a MPC at interfaces to old mesh")' ...
  '-lc(#%g#"ask for chararecteristic length on target")' ...
  '-lcmin(#%g#"do not iterate if estimated refinement would be smaller than that")' ...
  '-keepi(#3#"keep initial selection volume with lc iterations")' ....
  ]}};

cinM.add={
'gr:Cube','2cube example', ...
   {'ALoad','cbuild','cqoff','drill','exp','flipSlaveMaster','largeBase',...
     'mixtq','nCube','offset','OnPlate','penta','quad','refine','single',...
     'tab','tetra'}
      };
cinM.add={
'gr:ScldCube','Large displacement cube', ...
  {['Kc(1e12#%g#"NomStiff") ' ...
   'stepx(0#%g#"step x")', ...
   'Integ(-2#%g#"contact integration strategy")', ...
   'lxyz(.075 .075 .075#%g#"cube dim")' ...
   'DownForward(#%s#"Possibly trajectory as parameter")' ...
   'InitQ0(#%s#"command to init q0")' ...
   ],'largeBase','ALoad'} ... % do not redefine
  };
% vhandle.uo('',C3.info,rail19('nmap.Map:Cin'))

%% #Bsections #Bprofiles: geometry data for mesh naca -2

RA=struct('nmap',vhandle.nmap,'ToolTip','Composite blade example');
p0=struct('section',... % (angle, xpos ypos)
 [[0; acos(.9); pi/2+asin(.3); pi; 3*pi/2; pi+acos(-.9); 2*pi] ...
 [1 0;.9 .003; .3 .05;0 0;.5 -.02;.9 -.001;1 0]],...
 'EndSlope',[0 0]);
p1=p0; p2=p0;
p1.section(:,2)=120*p1.section(:,2); p2.section(:,2)=60*p2.section(:,2);
p1.section(:,3)=1e2*p1.section(:,3); p2.section(:,3)=5e1*p2.section(:,3);
p3=p2; p3.skew=struct('Orig',[40 0 0],'angle',-45); % skewness
secM=vhandle.nmap;
secM('naca66_120')=p1; secM('naca66_60')=p2; secM('naca66_60_rn45')=p3; 

p1=p0;
p1.section(:,2)=250*p1.section(:,2);p1.section(:,3)=1e2*p1.section(:,3);
secM('naca66_250e7')=p1;

p1=p0; % xxx pre-list may be better 'rad', 'xoff', 'sec', 'len', 'hei', 'skewO', 'skewA'
p1.section(:,2)=300*p1.section(:,2);p1.section(:,3)=1e2*p1.section(:,3);
secM('naca66_300e7')=p1;
p1=p0;
p1.section(:,2)=350*p1.section(:,2);p1.section(:,3)=1e2*p1.section(:,3);
secM('naca66_350e7')=p1;
p1=p0;
p1.section(:,2)=520*p1.section(:,2);p1.section(:,3)=8e1*p1.section(:,3);
secM('naca66_520e5p6')=p1;
p1=p0;
p1.section(:,2)=555*p1.section(:,2);p1.section(:,3)=6e1*p1.section(:,3);
secM('naca66_555e4p2')=p1;
p1=p0;
p1.section(:,2)=450*p1.section(:,2);p1.section(:,3)=5e1*p1.section(:,3);
secM('naca66_450e3p5')=p1;
p1=p0;
p1.section(:,2)=350*p1.section(:,2);p1.section(:,3)=2e1*p1.section(:,3);
secM('naca66_350e1p4')=p1;
p1=p0;
p1.section(:,2)=200*p1.section(:,2);p1.section(:,3)=5e-1*p1.section(:,3);
secM('naca66_200e0p35')=p1;


RA.nmap('Map:Bsections')=secM;

profM=vhandle.nmap;
%'rad', 'xoff', 'sec' -> x,y,ry,sec,scale
profM('HyFoilA')={'rad','xoff','sec';0,0,'naca66_120';180,30,'naca66_60'};
profM('HyFoilB')={'rad','xoff','sec';0,0,'naca66_120';180,30,'naca66_60_rn45'};
%
profM('HyFoilC')={'rad', 'xoff', 'sec', 'len', 'hei', 'skewO', 'skewA';
 0,0,'naca66_250e7',250,7,[175 0 0],0
 137,-2 ,'naca66_300e7',298,7,[175 0 0],15
 192,-7 ,'naca66_350e7',343,7,[175 0 0],30
 356,-20,'naca66_520e5p6',500,5.6,[175 0 0],45
 438,-15,'naca66_555e4p2',540,4.2,[175 0 0],60
 493,0,'naca66_450e3p5',450,3.5,[175 0 0],70
 530,75,'naca66_350e1p4',425,1.4,[175 0 0],77
 548,200,'naca66_200e0p35',400,.035,[175 0 0],80};

RA.nmap('Map:Bprofiles')=profM;

plyM=vhandle.nmap;
%% #plyA base plies test -3
RP=struct('plyList',{{ ... 
 'name','thick','matid','theta','rstop'
 'ply1' .15  1   0  Inf
 'ply2' .15  2  45  150
 'ply3' .15  3  45   75
 'ply4' .15  4  45   50 
 'core' Inf  5   0  Inf % Symmetry 
  }},...
  'OrientLine',struct('starts',10,'dir',[0 0 1]));
plyM('plyA')=RP;
%% #plyB various thickness and split core to control mesh size -3
RP=struct('plyList',{{ ... 
 'name','thick','matid','theta','rstop'
 'ply1' .15   1 0    Inf
 'ply2' .25   2 45   150
 'ply3'  .05  3 45   75
 'ply4'  .15  4 45   50 
 'core'  1  5  0   Inf
  'core'  1  5  0   Inf
  'core' Inf  5  0   Inf
  }},...
  'OrientLine',struct('starts',10,'dir',[0 0 1]));
plyM('plyB')=RP;
% unsymm case with various plies
RP=struct('plyList',{{ ... 
 'name','thick','matid','theta','rstop'
 'ply1' .15  1   0  Inf
 'ply2' .15  2  45  150
 'ply3' .15  3  45   75
 'ply4' .15  4  45   50 
 'core' Inf  5   0  Inf
 'ply5' .2  6  45   50 
 'ply1' .15  1  0  Inf
  }},...
  'OrientLine',struct('starts',10,'dir',[0 0 1]));
plyM('plyC')=RP;
% case with visco on extrados and rstop as contour
moC=struct('Node',[1 0 0 0 10 0 10;2 0 0  0 40 0 200 ;3 0 0 0  70 0 200;4 0 0 0 90 0 10],...
 'isClosed',1);
RC=feval(lsutil('@dToPoly'),'init',moC);
RP=struct('plyList',{{ ... 
 'name','thick','matid','theta','rstop','hforced','htrans'
 'cply' .3       101      0       RC      0         0 
 'visc' .2        201      0       RC      1         1
 'ply1' .15         1      0      Inf      1         1
 'ply2' .15         3     45       75      0         1
 'core' Inf         5      0      Inf      0         1 % Symmetry 
 'ply2' .15         3     45       75      0         1
 'ply1' .15         1      0      Inf      1         1
  }},...
  'notsym',1,  'OrientLine',struct('starts',10,'dir',[0 0 1]));
plyM('plyVe')=RP;

%% #plyD base plies test -3
RP=struct('plyList',{{ ...
 'name','thick','matid','theta','rstop'
 'ply1' .15  1   0  Inf
 'ply2' .15  2  45  500
 'ply3' .15  3  45  350
 'ply4' .15  4  45  150
 'core' Inf  5   0  540 % Symmetry
 }},...
 'OrientLine',struct('starts',10,'dir',[0 0 1]));
plyM('plyD')=RP;

RA.nmap('Map:Bplies')=plyM;
%% #Naca_Map:MatName definition -3
matM={'cply','m_elastic(dbval101 UD)'
 'ply1','m_elastic(dbval1 ortho1 -unitMM)' % To check unit conversion
 'ply2','m_elastic(dbval2 ortho1)' 
 'ply3','m_elastic(dbval3 ortho1)' 
 'ply4','m_elastic(dbval4 ortho1)'
 'ply5','m_elastic(dbval5 ortho1)' 
 'visc','m_visco(database201 Smactane 50_G)'
 };
matM=vhandle.nmap(matM,[],'Map:MatDB');
RA.nmap('Map:MatDB')=matM;


RA.nmap('NacaAA')={'MeshCfg{d_mesh(Naca{HyFoilA:plyA,zs.1,yn1,unitmm}):empty}','RunCfg{feplot}'};
RA.nmap('NacaBA')={'MeshCfg{d_mesh(Naca{HyFoilB:plyA,zs.1,yn1,unitmm}):empty}','RunCfg{feplot}'};
RA.nmap('NacaAB')={'MeshCfg{d_mesh(Naca{HyFoilA:plyB,zs.1,yn1,unitmm}):empty}','RunCfg{feplot}'};
RA.nmap('NacaAC')={'MeshCfg{d_mesh(Naca{HyFoilA:plyC,zs.1,yn1,unitmm}):empty}','RunCfg{feplot}'};
RA.nmap('NacaAVe')={'MeshCfg{d_mesh(Naca{HyFoilA:plyVe,zs.1,yn1,unitmm}):ParVeCut}','RunCfg{feplot}'};
RA.nmap('NacaCD')={'MeshCfg{d_mesh(Naca{HyFoilC:plyD,interp,zs.1,yn1,unitmm}):empty}','RunCfg{feplot}'};
projM('Naca')=RA;

 % sdtm.stdNmapOut('call')
 if nargout==1;out=sdtm.stdNmapOut(nmap,key,nargout,CAM);
 elseif nargout>1;[out,out1,out2]=sdtm.stdNmapOut(nmap,key,nargout,CAM);
 else; sdtm.stdNmapOut(nmap,key,nargout,CAM);
 end

%% clean end
elseif comstr(Cam,'@');out=eval(CAM);
%% #Tuto: recover model from a specific tuto step -3
elseif comstr(Cam,'tuto'); 
 eval(sdtweb('_tuto',struct('file','d_mesh','CAM',CAM)));
 if nargout==0; clear out; end

elseif comstr(Cam,'cvs')
 out=sdtcheck('revision','$Revision: eeb3f39 $  $Date: 2022-03-31 12:12:08 +0200 $ ');
else; error('%s unknown',CAM);
end 
%% #End function
end
%% #SubFunc ------------------------------------------------------------------
%% dThick: thickness measure for meshNACA models - - -------------------------
function out=dThick(xyz,R1,mo1)

if isempty(xyz); out=[]; return; end
if ischar(xyz)
 if strcmpi(xyz,'init')
  % recover midplane and thickness values at nodes
  mo5=mo1; mo5.Elt=feutil('selelt setname midPlane',mo5);
  mo5=fe_quality('cleanDegenRecast',mo5); % we have degen elts, need to treat quad4
  mo5.Node=feutil('getnodeGroupall',mo5);
  mo5=feutil('quad2tria',mo5); % quicker match on tria
  R1.moPlane=mo5;
  opt=stack_get(mo1,'info','MeshOpt','getdata'); R1.PThick=opt.PThick; % thickess value at mid
  out=R1; 
  % xxx check thickness value
 elseif strcmpi(xyz,'ro');out=R1;
 else; error('%s unknown',xyz);
 end
 return
end
% compute thickness at any node: project on midplane interp thickness value from sticknode
% remove distance to sticknde along normal, remove offset
mo5=R1.moPlane;
n0=struct('Node',xyz,'InterpNormal',1);
r5=feutilb('matchsurf-radiusinf',mo5,n0,'groupall');
% then you want value at StickNode and distance from sticknode along normal
r6=feutilb('mpcfrommatch -entry',mo5,r5); % gives interp from sticknode to matched elt
r6=feutilb('placeindof',fe_c(r6.DOF,.01,'dof',1),r6);
r6.c=r6.c(1:6:end,:);
%
d1=R1.PThick;
r6.DOF=r6.DOF+.97;
r6=feutilb('placeindof',d1.DOF,r6);
r7=-r6.c*d1.def; % thickness value at sticknode location (interp on midplane)

out=r7-sqrt(sum(((r5.Node-r5.StickNode).*r5.InterpNormal).^2,2)); % substract distance to sticknode

if isfield(R1,'thickness'); % offset
 out=out-R1.thickness;
end

end



function   RO=buildContour(RO)
   %sdtw('_ewt','sdtu.fe.splineR document calls t_fmesh')
   % Generate profile curve (angle, xpos ypos), provide
   r1 = spline(RO.curProf.section(:,1)',...
    [RO.curProf.EndSlope;RO.curProf.section(:,2:3);RO.curProf.EndSlope]');
   s1=linspace(0,2*pi,RO.sRef*2*round(RO.xn/2)); % angle refinement
   r2=ppval(r1,s1); % interpolated curve
   %figure(11); plot(r2(1,:),r2(2,:),'-b',RO.curProf.section(:,2),RO.curProf.section(:,3),'or'),axis equal

   % prepare contour model, with closed line
   n0=[(1:size(r2,2)-1)' zeros(size(r2,2)-1,3) r2(:,1:end-1)' zeros(size(r2,2)-1,1)];
   RO.contour=n0; RO.contourXY=r2;
end

