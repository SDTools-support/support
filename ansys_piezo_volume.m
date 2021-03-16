% Script to illustrate PIEZO analysis using ANSYS element formulations
%
% This scrip is very close to the MFC_actuator demo, it just uses an ANSYS
% model. 
% This requires ADDITIONAL ANSYS .RST and .EMAT files (.CDB is optional)

% Check that recent enough versions of functions are available
sdtkey('cvsnum >1.075','p_piezo')
sdtkey('cvsnum >1.073','fe_reduc')


% Read ANSYS model into UPCOM element - - - - - - - - - - - - - - - - 
cf=feplot; ans2sdt('BuildUp File')
fe_case(cf.mdl,'reset') % Reset boundary conditions to write in SDT
RO.Case=2; % Piezo is thick or thin

% Assemble model so that SDT will not try to do it (we want to use ANSYS
% element formulations)
[model,Case]=fe_case('assemble -SE -matdes 2 1 3 NoT',cf.mdl.GetData);
model=feutil('rmfield',model,'mind','bas','wd','file');

% The are 2 electrodes : bottom and top of piezo volumes
% which are polarized along z

InputDOF=[];
% Build a MPC defining a single potential for the electrodes
[model,InputDOF(end+1,1)]=p_piezo('ElectrodeMPC Bottom',model,'z==0');
if RO.Case==2
    [model,InputDOF(end+1,1)]=p_piezo('ElectrodeMPC Top',model,'z==-2e-5');
else
    [model,InputDOF(end+1,1)]=p_piezo('ElectrodeMPC Top',model,'z==2e-6');
end
% Clamp one edge
model=fe_case(model,'FixDof','Clamp','x==0 -DOF 1 2 3');
% set bottom electrode to potential V=0
model=fe_case(model,'FixDof','V=0 on bottom',InputDOF(1));
% Fix electric DOFs of non piezo volumes (ANSYS keeps electric elements 
% on all)
model=fe_case(model,'FixDof', ...
    'NonPiezo',model.DOF(find(diag(model.K{2}==0))));

% Define the electrodes
InputDOF(:,2)=0; % closed circuit electrodes (actuators)
data=struct('data',InputDOF(2,:), ... % each potential and type (act/sen)
 'def',1,'DOF',InputDOF(2,1), ... % combine electrodes into inputs
 'lab',{{'Act'}});


model=stack_set(model,'info','Electrodes',data);
model=stack_set(model,'info','DefaultZeta',.01);

% Force setting of case to avoid reassembly
[Case,CaseName]=fe_case(model,'getcase');
Case=fe_case(model,'gett');model=stack_set(model,'case',CaseName,Case);
[sys,TR]=fe2ss('free 5 10 0 -dterm',model);

qbode(sys,linspace(1e5,3e6,512)'*2*pi,'iiplot "Test" -po');

% See the static response to the actuators
cf.def=fe_simul('static',model);


% - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Transient response analysis with voltage excitation (at node 106,
% node 1 is fixed 0 volt as ground point).

t=linspace(0,4e-4,1e4);
%c1=fe_curve('test step 1e-4',t);
%c1=fe_curve('test CosHan 5e4 10 1',t);
c1=fe_curve('test ramp 0 1e-4 1',t);

x=lsim(sys,c1.Y,t); % Just charge
figure(1);plot(t,x)

% Extract all states to allow animation, do not use d term there
[sys2,TR]=fe2ss('free 5 10 0',model);
N=size(TR.def,2);
sys2=ss(sys2.a,sys2.b,eye(N,2*N),zeros(N,1));
x=lsim(sys2,c1.Y,t); % all modal amplitudes
def=struct('def',x','DOF',(1:N)'+.99,'TR',TR, ...
    'data',t(:),'LabFcn','sprintf(''%.1f \\mus'',def.data(ch)*1e6)', ...
    'fun',[0 4]);
cf.def=def;
cf.sel={'matid1','ColorDataEvalX'};
fecom ch100;fecom colorbar % Then clik on the animation




