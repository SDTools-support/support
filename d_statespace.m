function d_statespace(varargin)

% Examples of SDT based state-space usage
%
% ResultantDofSet : case with enforced acceleration and resultant
% ResultantDofLoad : case with applied load an resultants at two locations

if nargin==0; Cam='';
else; [CAM,Cam]=comstr(varargin{1},1);carg=2;
end
%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC>
if isempty(Cam)

%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Problem with resultant response
elseif comstr(Cam,'resultant');[CAM,Cam]=comstr(CAM,10); 

if comstr(Cam,'dofset')
%% Our customer ask me about the sensor definition as nodal reaction force.
%When we want to get the state space model with input as the imposed DOF
%displacement(DofSet) excitation 
%and output as the reaction force on the DOF,  how can we define the input
%and output setting?".


model=femesh('testbeam1 divide 10');
model=fe_case(model,'FixDof','2D',[.01;.03;.04;.05], ...
    'DofSet','IN',struct('DOF',[1.02;1.06],'def',[1;0]), ...
    'SensDof','Tip',2.02);
% Load input
if 1==1
model=fe_case(model,'remove','Tip', ...
    'DofLoad','IN',struct('DOF',[1.02;1.06],'def',[1;0]));
end
% Define resultant sensor
Sensor=struct('EltSel','withnode 1');
Sensor.SurfSel='nodeid 1';
Sensor.dir=[0.0 1.0 0.0];
Sensor.type='resultant';
model=fe_case(model,'SensDof append resultant','Reaction',Sensor);
model=fe_case(model,'sensmatch','Reaction');


def=fe_eig(model);
cf=feplot(model,def);fecom view2

% This does not work (not a supported case)
%sys=fe2ss('free 5 6',model);

model=stack_set(model,'info','EigOpt',[5 6 0]);
SE=fe_reduc('craigbampton',model);
dr=fe_eig({SE.K{:},[]});dr.def=SE.TR.def*dr.def;dr.DOF=SE.TR.DOF;
sys=nor2ss(dr,model); % Acc In to disp out
qbode(sys,2*pi*linspace(0,2e3,1000)','iiplot "Test" -po');

elseif comstr(Cam,'dofload')
%------------------------------------------------------------------
%% Now a case with force applied at free tip (node 2) 
% and reaction at cantilever point (node 1)
%

model=femesh('testbeam1');
if ~isempty(strfind(Cam,'very'))
 model=feutil('divide',model,[0:.1:.8 .9:.01:.98 linspace(.99,1,5)]);
 RO.typ='VeryFine';
elseif ~isempty(strfind(Cam,'fine'))
 model=feutil('divide 10',model,[0:.1:.8 .9:.01:1]);
 RO.typ='Fine';
else % coarse
 RO.typ='Coarse'; model=feutil('divide 10',model);
end
model=fe_case(model,'FixDof','2D',[1;.01;.03;.04;.05], ...
    'DofLoad','IN',struct('DOF',[2.02;2.06],'def',[1;0]));
% Define resultant sensor
Sensor=struct('EltSel','withnode 1');
Sensor.SurfSel='nodeid 1';
Sensor.dir=[0.0 1.0 0.0];
Sensor.type='resultant';
model=fe_case(model,'SensDof append resultant','ReactionBase',Sensor);
model=fe_case(model,'sensmatch','ReactionBase');
Sensor.EltSel='withnode 2';Sensor.SurfSel='nodeid 2';
model=fe_case(model,'SensDof append resultant','ReactionTip',Sensor);
model=fe_case(model,'sensmatch','ReactionTip');

model=stack_set(model,'info','EigOpt',[5 6 1e3],'info','DefaultZeta',.01);
sys=fe2ss('free -dterm',model);
qbode(sys,2*pi*linspace(1,2e3,1e4)',sprintf('iiplot "%s" -po',RO.typ));
ci=iiplot;ci.Stack{RO.typ}.X{2}={'Base';'Tip'};
iicom(';sub 1 1 1 1 1;ch1:2;xlog')

%% ---------------------------------------------------------------------
else;
    error('Resultant%s unknown',CAM);
end

%% ---------------------------------------------------------------------
elseif comstr(Cam,'cvs')
 out=sdtcheck('revision');
    %out='$Revision: 134 $  $Date: 2014-02-26 15:27:26 +0100 (Wed, 26 Feb 2014) $';
else;error('%s unknown',CAM);
end
