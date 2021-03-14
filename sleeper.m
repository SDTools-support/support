% this example illustrates the use of the FEMESH preprocessor to build a
% solid model of a U-beam, compute the associated modes, and display strain
% energy levels
%
% See also demos gartfe, beambar, demo_fe, d_plate
%          doc   fem, feco, dfeplot

femesh('reset');
% Base nodes
FEnode=[1 0 0 0  0 0 0;
        3 0 0 0  0 0 0
        4 0 0 0  0 0.175 0;
        5 0 0 0  0.285 0 0;
        6 0 0 0 0.285 0.175 0;];

% define U shaped base surface
FEelt=[Inf abs('quad4');4 6 5 3 1 1];
FEel0=[Inf abs('quad4');4 6 5 3 1 1];

%FEel0=[Inf abs('quad4');1 2 4 3 1 1;11 12 10 9 1 1];
femesh(';divide -new 1 8;addsel;');
%FEel0=[Inf abs('quad4');4 6 7 9 1 1];
femesh(';divide -new 1 10;addsel;');

femesh('join group 1:3');

% extrude the base to form the beam model
%femesh(';selgroup2;extrude 200 0 0 .0125;orientel0');
femesh(';selgroup2;extrude 100 0 0 .025;orientel0');

model=femesh('model0');
model.pl = m_elastic('dbval 1 concrete'); % isotropic elastic (steel)
model.il = p_solid('dbval 1 d3 2');    % integration rule for 3D element

% defining boundary conditions : fix base
%model=fe_case(model,'FixDof','Clamped end','z==0')

% Eigen solver 6, 10 modes, no shift, info level 11
% xxx shift should not be 0 in Free/Free
model=stack_set(model,'info','EigOpt',[6 20 1e3 11]);

% Compute the modes
def = fe_eig(model);

% Finally, one can display the model and strain energy distribtion
cf=feplot(model,def); 

fecom(';view3;view s-90;colorfacew;showpatch');
fecom('colordata ener k');

