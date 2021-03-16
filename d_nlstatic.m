function [out,out1]=d_nlstatic(varargin)

%D_NLSTATIC: Example database for non linear static state resolution
%            and associated pre-stressed analyses
%
%
%       See also : fe_time, nl_solve, fe_case, mattyp

%	Etienne Balmes, Guillaume Vermot des Roches
%       Copyright (c) 1990-2021 by SDTools, All Rights Reserved.


if nargin==0; return; end % auto test skip

[CAM,Cam]=comstr(varargin{1},1);carg=2; 
%#ok<*ASGLU,*NASGU,*NOSEM>

if comstr(Cam,'nlstress')
%% #NLStress: basic computation of NL stresses on a sample element
%mo1=femesh('testhexa20b')
mo1=femesh('testhexa8b');
n1=basis('gnode','ry=5',mo1.Node);
dStat=struct('data',0,'DOF',feutil('getdof',(1:3)'/100,mo1.Node(:,1)),...
 'def',reshape(n1(:,5:7)-mo1.Node(:,5:7),[],1));

mo1=stack_set(mo1,'curve','StaticState',dStat);
sdtw('_ewt','@EB review proper outputs for Large Transformation elastic stresses')

Ek=fe_stress('Ener -matdes 5 -curve',mo1,dStat); % Need discussion @gv

Ek.Y
%cg=feplot(55); cg.model=mo1;cg.def=d1;fecom('scc1');

% Next step, problem with "pres-stressed stiffness"
mo2=mo1;
mo2=stack_set(mo2,'curve','StaticState',dStat);
[mo2,C2]=fe_case('assemble NoT -se -matdes 2 1 5',mo2);
% Compute the Green Lagrange strain tensor 
C2=fe_stress('stress-gstate -matdes 5',mo2,dStat);C2=C2.GroupInfo{5};

d1=fe_eig([mo2.K([1 2]) mo2.DOF],[5 10 1e3]);
d5=fe_eig([mo2.K([1 3]) mo2.DOF],[5 10 1e3]);
% A rigid body rotation does not affect the system modes:
[d1.data d5.data]
 
elseif comstr(Cam,'large')
%% #Large: sample resolution of a static with nlgeom, with advanced controls
model=demosdt('largeTransform');
model=stack_set(model,'info','TimeOpt', ...
   struct('Opt',[],'Method','nl_solve NewtonStatic',...
   'Jacobian','ki=basic_jacobian(model,ki,0.,0.,opt.Opt);',...
   'NoT',1, ... % Don't eliminate constraints in model.K
   'AssembleCall','assemble -fetimeNoT -cfield1', ...
   'IterInit','opt=fe_simul(''IterInitNLStatic'',model,Case,opt);',...
   'MaxIter',200,'RelTol',1e-6,'rIncTol',1,'MaxSDI',0,'CutB',2.5));

model=fe_case(model,'setcurve','PointLoad', ...
    fe_curve('testramp NStep=20 Yf=2e-6')); % 20 steps gradual load
d1=fe_time(model);
ci=iiplot(d1);iicom('ch',{'DOF',288.03}) % View response

elseif comstr(Cam,'modal')
%% #Modal: sample prestressed modal analysis


%% -------------------------------------------------------
elseif comstr(Cam,'@');out=eval(CAM);
elseif comstr(Cam,'cvs')
 out=sdtcheck('revision');
  %out='$Revision: 206 $  $Date: 2015-03-16 14:24:10 +0100 (Mon, 16 Mar 2015) $';
else; error('%s unknown',CAM);    
end
