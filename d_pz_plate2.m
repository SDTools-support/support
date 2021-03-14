
%% Example answering the test case of Gerard Coffignal
%       Copyright (c) 2014-2015 by SDTools, All Rights Reserved.


%% Generic meshing : see documentation 
sdtweb d_piezo('MeshPlate')

mdl=struct('Node',[],'Elt',[], ... % empty model
    'pl', ... % plate
      [1 fe_mat('m_elastic','SI',1) 70e9 .3 2700], ... 
     'il', ... % laminate definition (1 3mm layer at 0)
     p_shell('dbval 1 -punit mm -unit mm laminate    1 3.0 0  '), ...
     'unit','mm'); 
mdl.il=fe_mat('convertMM',mdl.il);
mdl.pl=fe_mat('convertMM',mdl.pl);
feutilb('_write',mdl)


 RN.list={'Name','Lam','shape'
   'Main_plate', mdl,struct('type','global','lx',500,'ly',50,'lc',8)
   'A1','BaseId1 +SmartM.MFC-P1.2814 -SmartM.MFC-P1.1020 ', ...
     struct('shape','rect','xc',25,'yc',28,'alpha',0)
   'A2','BaseId1 +SmartM.MFC-P1.2814 -SmartM.MFC-P1.1020 ', ...
     struct('shape','rect','xc',25,'yc',2,'alpha',0)
   'S1','BaseId1 +SmartM.MFC-P1.2814 -SmartM.MFC-P1.1020 ', ...
     struct('shape','rect','xc',100,'yc',28,'alpha',0)
   'S2','BaseId1 +SmartM.MFC-P1.2814 -SmartM.MFC-P1.1020 ', ...
     struct('shape','rect','xc',100,'yc',2,'alpha',0)
   };
model=d_piezo('MeshPlate',RN);model.name='Plate with piezo';
model=fe_case(model,'fixdof','Base','x==0');

cf=feplot(model);
p_piezo('electrodeinfo',model);
matgui('jil',cf);matgui('jpl',cf); % Di

%% Open documentation
sdtweb pz_shell#combined_electrodes

% combine electrodes to generate pure bending / pure traction
data.def=[1 -1 0 0  0 0 0 0 % A1
          0 0 1 -1  0 0 0 0 % A2
          0 0 0 0   1 -1 0 0 % S1
          0 0 0 0   0 0 1 -1 % S2          
]'; % Define combinations for actuators
data.lab={'A1-Bend';'A1-Bend';'S1-Bend';'S2-Bend'};
data.DOF=p_piezo('electrodeDOF.*',model);
model=fe_case(model,'DofSet','V_In',data);
model=stack_rm(model,'info','Electrodes');

def=fe_eig(model,[5 20 0]); % Compute modes

% Forced response at 95% of first frequency
d0=fe_simul('dfrf',stack_set(model,'info','Freq',def.data(1)*.95)); 

% Do the combination A1=100;A2=100;S1=0;S2=0;
q0=d0; q0.def=q0.def*[100 100 0 0]';q0.data=q0.data(1);
q0.LabFcn='sprintf(''Combined response at %g Hz'',def.data(ch,1))';
cf.def=q0; fecom('colordataEvalZ-edgeAlpha.1')

