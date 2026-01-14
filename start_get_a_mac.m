
%%Getting ready ----------------------------------------------------------
% Dowload test files can take minutes, you can do it by hand
cd(sdtdef('tempdir'))
demosdt('download-back http://www.sdtools.com/contrib/cc_mode.op2')
demosdt('download-back http://www.sdtools.com/contrib/cc_test.uff')
demosdt('download-back http://www.sdtools.com/contrib/cc_test_wire.unv')
ls('cc_*')



%% Read test node positions from UNV ------------------------------------
test=ufread('cc_test_wire.unv');
test=feutil('rmfield',test,'pl','il','header','Stack');% remove unused fields

% Read responses which will give sensor direction information
UFS=ufread('cc_test.uff');

test.tdof=UFS(1).dof(:,1); % define sensor dirs
% Here there is a reorientation needed. See details in 
% http://www.sdtools.com/helpcur/ccdemo.html#cccor
test.tdof=fe_sens('tdof',test.tdof);
test.tdof(:,3:5)=test.tdof(:,3:5)*[0,0,-1;-1,0,0;0,1,0];


%% Read mesh and modes from NASTRAN typically generated with -----------
% PARAM,POST,-2 or -4
cf=feplot(2);nasread('cc_mode.op2');
cf.def=cf.Stack{'BOPHIG.1'};

% Perform topology correlation

cf.mdl=fe_case(cf.mdl,'sensdof','Test',test);
fe_case(cf,'sensmatch radius 30','Test')
fecom(cf,'curtabCases','Test');fecom(cf,'Proviewon');
sens=fe_case(cf.mdl,'sens');

% View FEM modes in test wire frame
def=cf.Stack{'BOPHIG.1'};def.name='FEM';
def=fe_def('subdef',def,7:size(def.data,1));
def=feutilb('placeindof',sens.DOF,def);
dr=struct('def',sens.cta*def.def,'DOF',sens.tdof, ...
    'data',def.data(:,1));
cg=feplot(5);feplot(cg,test,dr);fecom(cg,'colordata evaly')

%% Now compute full and partial MAC
%idenfity modes ----------------------------------------------
ci=idcom;ci.Stack{'Test'}=UFS(1);iicom('submagpha')
idcom(';e 1415 .01;ea;e .01 3609 ;ea')

% View MAC
figure(1);ii_mac(ci.Stack{'IdMain'},def,'sens',sens,'macplot')

% -----------------------------------------------------------
% Correlate on a subset of DOFs
%  - first strategy : define less sensors (not illustrated)
%  - second strategy : do the observation by hand

ID=fe_def('subdof',ci.Stack{'IdMain'},fe_c(test.tdof,.03,'ind'));
dr=struct('def',sens.cta*def.def,'DOF',sens.tdof, ...
    'data',def.data(:,1));
dr=fe_def('subdof',dr,fe_c(test.tdof,.03,'ind'));
figure(1);ii_mac(ID,dr,'macplot')


