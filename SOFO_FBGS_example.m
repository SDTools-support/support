%%% This is an example of defining fiber optic strain sensors in a concreate
%%% cantilever beam. Two types of strain sensors are considered : (i) short
%%% gauge sensors (point sensors) such as FBGS, and (ii) long gauge sensors,
%%% such as the SOFO system.
%%%
%%% Arnaud Deraemaeker, ULB, april 2010

% creating the model
femesh('reset');
FEnode=zeros(8,7);
FEnode(:,1)=linspace(1,8,8)';
FEnode(:,5:7)=[0 0 0;1 0 0;1 .1 0;0 .1 0;0 0 .1 ;1 0 .1;1 .1 .1;0 .1 .1];
FEel0=[Inf abs('hexa8') 0 0 0 0 0; 1 2 3 4 5 6 7 8 1 1 1];
femesh(';divide 50 5 5;addsel'); femesh(';set group 1 matid 1;set group 1 proid 1');
model=femesh('model'); feplot(model);
model.pl=m_elastic('dbval 1 Concrete');
model.il= p_solid('dbval 1 d3 10002'); %%% for better 3D elements (non compatible) 
[Eltid,model.Elt]=feutil('eltidfix',model.Elt);
model=fe_case(model,'fixdof','Encastrement','x==0'); % cantilever boundary condition
% computing the modeshapes,
def=fe_eig(model,[5 20 500]);
cf=fecom; cf.def=def; 

%%% Define local strain sensors (FBGS)
%%% Assume straight fiber parallel to beam axis

y=0.09; % y-position of the optical fibre
z=0.09; % z-position of the optical fibre
nsens=50; % number of FBGS (point sensors)
L=max(model.Node(:,5)); % length of the beam
d0=L/100; % distance from edge on both sides
pas=(L-2*d0)/(nsens-1); % distance between the sensors
Sensors=[[1:nsens]' [d0:pas:L-d0]' y*ones(nsens,1) z*ones(nsens,1) ones(nsens,1) zeros(nsens,1) zeros(nsens,1)];
model=fe_case(model,'SensDof append strain','FBGS',Sensors); % Define the sensors
model=fe_case(model,'sensmatch','FBGS'); % match with FE model
sens=fe_case(model,'sens','FBGS'); % create projection matrix in sens.cta

%%% look at projected modes
def=feutilb('placeindof',sens.DOF,def);
md=sens.cta*def.def;

x=Sensors(:,2); % x position of the FBGS

figure(10);
i1=[1 3 10 20];plot(x,md(:,i1));legend(cellstr(num2str(i1(:)))); hold on

%--------------------------------------------------------------------------
%%%%%%% long gauge sensors (SOFO) defined by two nodes
%%%% idea is to define to translation sensors and make the difference
%%%% (relative displacement sensor)

nsens=10; %% number of SOFO sensors (long gauge)
d0=L/100; %% distance from the edge on both sides
d1=L/50; %% distance between two sensors
SensA=[[1:nsens]' [d0:L/nsens:L-L/nsens+d0]' y*ones(nsens,1) z*ones(nsens,1) -ones(nsens,1) zeros(nsens,1) zeros(nsens,1)];
SensB=[[1:nsens]' [d0+L/nsens-d1:L/nsens:L-d0]' y*ones(nsens,1) z*ones(nsens,1) ones(nsens,1) zeros(nsens,1) zeros(nsens,1)];

model=fe_case(model,'SensDof append trans','SensA',SensA); % left point translation sensor in x
model=fe_case(model,'SensDof append trans','SensB',SensB); % right point translation sensor in x
model=fe_case(model,'sensmatch','SensA'); % match with FE model
model=fe_case(model,'sensmatch','SensB'); % match with FE model
sensA=fe_case(model,'sens','SensA'); % create projection matrix
sensB=fe_case(model,'sens','SensB'); % create projection matrix

sensA.cta=sensA.cta+sensB.cta; % sum the two projection matrices

md2=sensA.cta*def.def; % project the modeshapes on SOFO sensors

xlg=(SensA(:,2)+SensB(:,2))/2;%%% define mid position of long gauge sensors
bl=SensB(:,2)-SensA(:,2);%% divide SOFO reading by base length for comparison

figure(10); plot(xlg,md2(:,1)./bl,'r*'); plot(xlg,md2(:,3)./bl,'r*');
plot(xlg,md2(:,10)./bl,'r*'); plot(xlg,md2(:,20)./bl,'r*'); 
xlabel('x'); ylabel('\epsilon_x')
legend([cellstr(num2str(i1(:)));'SOFO'])


