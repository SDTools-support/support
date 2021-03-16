% This demo answers the question
%
%  generate transfer functions from modes read in NASTRAN
%  generate State-space model from external FEM solution

%%Getting ready ----------------------------------------------------------
% Dowload test files can take minutes, you can do it by hand
cd(sdtdef('tempdir'))
demosdt('download-back http://www.sdtools.com/contrib/cc_mode.op2')
ls('cc_*')

%% Read mesh and modes from NASTRAN typically generated with -----------
% PARAM,POST,-2 or -4
model=nasread('cc_mode.op2');
%  For the structure of model, see sdtweb('model')
% FEM results are stored in model.Stack, see sdtweb('stack')
% Access to stack can be named based, see sdtweb('stack_get')
cf=feplot(2);cf.model=model;
cf.def=cf.Stack{'OUG(1)'}; fecom(';colordata evaly;ch7')

% Define inputs and outputs (here on specific DOF)
%  for details on DOF numbering, see sdtweb('mdof')
%  for more general sensors, see sdtweb('sensor')
cf.mdl=fe_case(cf.mdl,'DofLoad','IN',425.02);
cf.mdl=fe_case(cf.mdl,'SensDof','OUT',[425.02;911.02]);

def=cf.Stack{'OUG(1)'}; % Get Modes
def.data=def.data(:,1); % if multiple columns, damping in column 2
Freq=linspace(1000,3500,1024)';cf.Stack{'info','Freq'}=Freq;

% Compute for 2 damping levels
ci=iiplot(3);
nor2xf(def,.0,cf.mdl,'iiplot "FRF(0)"'); % No damping
nor2xf(def,.02,cf.mdl,'iiplot "FRF(.02)"'); % 2% damping
iicom('submagpha')

% The same using a state-space model
ss=nor2ss(def,.02,cf.mdl);
if exist('bode','file'); 
    figure(10);bode(ss,cf.Stack{'Freq'}*2*pi); % Using BODE in control toolbox
else qbode(ss,cf.Stack{'Freq'}*2*pi,'iiplot "SS"');    % Using SDT QBODE
end
