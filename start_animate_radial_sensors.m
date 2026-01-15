
% Start the example by building a tube FEM model
femesh('reset');
FEnode=[1  0 0 0  0.0 0.0 0.0; 2  0 0 0  0.0 1.0 0.0];
femesh(';objectbeamline 1 2;rev 10 o 0 0 .1 360 0 1 0;divide 1 3;addsel');
model=femesh;
model=fe_case(model,'reset','fixdof','base','y==0', ...
    'DofLoad','in',2.03);
model.pl=m_elastic('dbval 1 Aluminum');
model.il=p_shell('dbval 1 kirchoff 5e-2');

% For the sake of this example find normals to be used as sensor directions
MAP=feutil('getnormal map Node',model);
i1=~ismember(MAP.ID,feutil('findnode y==1',model));
MAP.ID(i1)=[];MAP.normal(i1,:)=[];

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Define radial sensors. This is done by the experimentalist
%     SensId      Node   dir_x dir_y dir_z
tdof=[MAP.ID+.01  MAP.ID MAP.normal        ];

model=fe_case(model,'SensDof trans','test',tdof); 
cf=feplot(model);

% Compute response
def=fe_eig(model,[5 20]); 
def.data(:,2)=.01; % damping
XF=nor2xf(def,model,linspace(0,2000,1024)','struct');
% go to the .w .xf format and display in IDCOM
XF1=struct('w',XF.X{1},'xf',reshape(XF.Y,[],10),'dof',XF.dof);

% Identify a modeshapes and display it
ci=idcom; ci.Stack{'curve','Test'}=XF1;iicom('submagpha')
idcom(sprintf(';e .01 %.2f;ea',def.data(1))) % estimate a pole and accept it

% Note how 'Test' is specified to force use of sensor orientations in 'Test'
% 
cf.def={ci.Stack{'IdMain'},'Test'}; 
fecom ShowLine




