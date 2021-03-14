% Modal transient example on a structure actuated by piezoelectrics

cf=ULB_plate_pzt('model feplot');  % creates the model
% define the electrode state
ULB_plate_pzt('cantileverclosed',cf);
% combine electrodes 2 by 2
cf.Stack{'Electrodes'}.def=[1 -1 0 0;0 0 1 -1]';
cf.Stack{'Electrodes'}.lab={'Ua';'Ub'};
cf.mdl=stack_set(cf.mdl,'info','DefaultZeta',.01);

% Generate modal model with static correction - - - - - - - - - - - - - 
if sdtkey('cvsnum','fe2ss')<1.033; error('FE2SS too old');end

cf.mdl=fe2ss('free 5 10 0 -se2',cf.mdl);
% Define transient loading (cleaning the set call would be needed)
cf.Stack{'MVR'}.Case.Stack{1,3}.curve={ ...
    fe_curve('TestRicker .1 2') fe_curve('TestRicker .1 2')};

TimeOpt=fe_time('TimeOpt Newmark .25 .5 0 1e-3 1000');
TimeOpt.NeedUVA=[1 1 0];
%TimeOpt.AssembleCall='assemble -fetime';
def=fe_time(TimeOpt,cf.mdl); 

cf.def={def,cf.Stack{'MVR'}.TR}; % Animate with restitution
fecom('view3'); fecom('animtime3')

