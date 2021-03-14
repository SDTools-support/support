% ----------------Description :---------------------
% This function models the GARTEUR SM-AG-19 Testbed
%       using BEAMS,One Celas and Tip Masses
%
% Contributed by Linda Mthembu, University of Johannesburg
% and edited by SDTools
% --------------------------------------------------

%------------Variables to be changed later--------------
v=0.3;p=2700;G=2.7692e10; uE1=7.2e10;

% ------------------------------------------------------
%                   Meshing
femesh('reset')
%[NodeId PID DID GID x   y     z]
FEnode= [01 0 0 0 0.000 0.00 0.075;02 0 0 0 0.120 0.00 0.075;
         03 0 0 0 0.600 0.00 0.075;04 0 0 0 0.770 0.00 0.075;
         05 0 0 0 1.450 0.00 0.075;06 0 0 0 1.450 0.00 0.170;
         07 0 0 0 1.450 0.20 0.455;08 0 0 0 1.450 0.10 0.455;
         09 0 0 0 0.600 1.00 0.171;10 0 0 0 0.400 1.00 0.171;
         11 0 0 0 0.500 1.00 0.171;12 0 0 0 0.600 1.00 0.171;
         13 0 0 0 0.600 0.90 0.171;14 0 0 0 0.600 0.00 0.171];
 
%-----1st segment of fuselage Group 1-----
l= 0.60; h= 0.150; w= 0.050;
pro=[1 1 0 0 1]; %ProId and direction
femesh('objectbeamline 1 2',pro);femesh('repeatsel 5 0.12 0 0');
femesh('addsel;');

%----2nd segment of fuselage Group 2----
l= 0.85; h= 0.150; w= 0.050;
femesh('objectbeamline 3 4',[2 2 0 0 1]);femesh('repeatsel 5 0.17 0');
femesh('addsel;');

% ------Vertical Tail Plane Group 3----
l= 0.300; h= 0.100; w= 0.010;
femesh('objectbeamline 5 6',[3 3 0 0 1]);femesh('repeatsel 4 0 0 0.095')
femesh('addsel;');

%---- Horizontal Tail Plane Group4-----
l= 0.400; h= 0.100; w= 0.010;
femesh('objectbeamline 7 8',[4 4 0 0 1]);femesh('repeatsel 4 0 -0.100 0');
femesh('addsel;');

%-----Right Drum Group 5--------
l= 0.400; h= 0.100; w= 0.010;
femesh('objectbeamline 10 11',[5 5 0 0 1]);femesh('repeatsel 4 0.1 0 0');
femesh('addsel;');

%----Left Drum    Group 6------
femesh(';symsel 1 0 1 0;addsel;');

%------Main WING Span Group 7------
l= 2.00; h= 0.100; h= 0.010;
femesh('objectbeamline 12 13',[7 7 0 0 1]);femesh('repeatsel 20 0 -0.100 0')
femesh('addsel;');

%------Wing to Fuselage connection Group 8----
il=[14];
femesh('object mass',il);
femesh('extrude 1 0 0 -0.096');
femesh('set groupa 1 name celas');
% set connected DOFs and spring value
FEel0(2:end,3)=123456;FEel0(2:end,4)=0;FEel0(2:end,7)=1e10;
femesh('addsel;');

%------Tip masses Group 9------
femesh('object mass',femesh('findnode y==1.00 | y==-1.00 & x==0.40'));
FEel0(2:end,2:4)=0.18; % drum masses
femesh('addsel');

femesh('object mass',femesh('findnode z==0.455 & y==0'));
FEel0(end,2:4)=.5; % tail mass
femesh('addsel');femesh(';join mass1;');

%                Assemble Model 
% --------------------------------------------------------
model = femesh('model');

% ---------Can also orient the beam elements using:-------
i1=feutil('findelt eltname beam1',model); % Indices of elements rows
%  Z-Axis as 1st Bending Plane(Up and Down of whole structure)
model.Elt(i1,5:7)=ones(size(i1))*[0 0 1]; % vx vy vz

model.pl= ...
    [1 fe_mat('m_elastic','SI',1) uE1 v p G;  % Mat Prop ID 1 = Aluminum
     2 fe_mat('m_elastic','SI',1) uE1 v p G % fuselage segment2
     3 fe_mat('m_elastic','SI',1) uE1 v p G % vertical tail plane
     4 fe_mat('m_elastic','SI',1) uE1 v p G % horizontal ail plane
     5 fe_mat('m_elastic','SI',1) uE1 v p G % Drum
     7 fe_mat('m_elastic','SI',1) uE1 v p G % Wings
     ];
 
% was initially 
il1=p_beam('dbval 1 rectangle 0.150 0.050'); % Sectn Prop
il2=p_beam('dbval 2 rectangle 0.150 0.050');
il3=p_beam('dbval 3 rectangle 0.1 0.01');
il4=p_beam('dbval 4 rectangle 0.01 0.1');
il5=p_beam('dbval 5 rectangle 0.01 0.4');
il7=p_beam('dbval 7 rectangle 0.01 0.1'); 
model.il = [il1;il2;il3;il4;il5;il7]; % Different Sectn Prop

% Easier reading with
model.il=p_beam('dbval 1 bar .05 .15', ... % fuselage segment1
          'dbval 2 bar .05 .15', ... % fuselage segment2
          'dbval 3 bar .01 .1', ... % vertical tail plane
          'dbval 4 bar .1 .01', ... % horizontal tail plane
          'dbval 5 bar .1 .01', ... % drum
          'dbval 7 bar .1 .01');  % wings

%---------------------------------------------------

%              Set eigen value solver options, 
%      method 5, NModes, Shift(1/10 1st flexible mode), Print, Thres 
model=stack_set(model,'info','EigOpt',[5 20 1e3]); 
def=fe_eig(model); % Modes
cf=feplot(model,def);fecom ch7

