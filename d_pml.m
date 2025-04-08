function [out,out1,out2]=d_pml(varargin); %#ok<*STOUT>

% D_PML Support for demonstrations related to Perfectly Matched Layers
% (PML)
%
% See <a href="matlab: sdtweb _taglist d_pml">TagList</a>
%
% Etienne Balmes, SDTools, Arnaud Deraemaeker, ULB


%       Copyright (c) 1990-2025 by SDTools, All Rights Reserved.
%       For revision information use d_pml('cvs')

if nargin==0
 sdtweb _taglist d_pml % see structure of d_pml file    
 return
end

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*DEFNU,*NOSEM>
[CAM,Cam]=comstr(varargin{1},1);carg=2;
if carg>nargin||~isstruct(varargin{2});RO=struct('info',1);
else;RO=varargin{carg};carg=carg+1;
end

%% #Script -------------------------------------------------------------------
if comstr(Cam,'script');[CAM,Cam]=comstr(CAM,7);

if comstr(Cam,'pml_1d')
%% #1D_PML : Demo script : see sdtweb ... 

elseif comstr(Cam,'ulb1d')
%% #Script_ULB1D : analytic verification on column
% Derived from PhD of Cedric Dumoulin at ULB

%% #Shear/pressure frequency domain -2

 % Pressure wave
 li={'MeshCfg{"d_pml(Ulb1d{v1,quad,Lc2})"}',';','RunCfg{dfrf{10 100},d_pml(View1DPS -cf3)}'};
 mo2=sdtm.range(struct,li);

 % Shear wave 
 li={'MeshCfg{"d_pml(Ulb1dSW{v1,quad,Lc2})"}',';','RunCfg{dfrf{10 50},d_pml(View1DPS -cf3)}'};
 mo2=sdtm.range(struct,horzcat(li{:}));

%comgui('imwrite',2);

%% #Pressure_wave_time domain formulation -2

% Not clean
mo1=d_pml('MeshUlb1D',struct('subtype',1,'Lc',2,'quad',1,'EdgeMpc',0,'pow',2,'Freq',[10 100],'Np',20));feplot(mo1);
d1=fe_simul('dfrf',stack_set(mo1,'info','Freq',[10;100]));d_pml('View1Dps -cf3',mo1,d1);feutilb('_write',mo1)

% Clean regular mesh
mo1=d_pml('MeshUlb1D',struct('subtype',1,'Lc',2,'quad',1,'EdgeMpc',0,'pow',2,'Lp',30));feplot(mo1);
d1=fe_simul('dfrf',stack_set(mo1,'info','Freq',[10;100]));d_pml('View1Dps -cf3',mo1,d1);feutilb('_write',mo1)


%% #Shear_wave_time domain formulation -2

mo1=d_pml('MeshUlb1D SW',struct('subtype',1,'Lc',2,'quad',1,'EdgeMpc',0));
% -3000 smooth, -2000 oscilations
mo1=feutil('setpro 100 a0 2000 pow2 form -3000',mo1);d1=fe_simul('dfrf',stack_set(mo1,'info','Freq',[10;100]));d_pml('View1Dps -cf3',mo1,d1);feutilb('_write',mo1)
mo1=feutil('setpro 100 a0 2000 pow2 form -2000',mo1);d1=fe_simul('dfrf',stack_set(mo1,'info','Freq',[10;100]));d_pml('View1Dps -cf3',mo1,d1);feutilb('_write',mo1)

if 1==2
 [SE,CE]=fe_case(mo1,'assemble -matdes 2 3 1 -SE -NoT',mo1);
 IA=CE.GroupInfo{1,7};IA.lab(~any(IA.data,2))
 p_pml('viewTopo',SE,CE);set(findall(1,'type','line','marker','.'),'color','k')
 s=10*2i*pi;Z=feutilb('tkt',CE.T,feutilb('sumkcoef',SE.K,[s^2 s 1]));
 [L,U]=lu(real(Z));i1=CE.DOF(abs(diag(U))>1e6);
 [fe_c(CE.DOF) num2cell(feutilb('dtkt',[],feutilb('tkt',CE.T,SE.K)))]
 [fe_c(CE.DOF) num2cell(diag(U))]
 
 [u,z]=eig(full(Z));u=Z\speye(size(Z,1),1);z=1;
 cf.def=struct('def',CE.T*u,'DOF',SE.DOF,'data',abs(diag(z)))
 
  iz=fe_c(CE.DOF,n1(end+[-5:0],4),'ind');
  [fe_c(CE.DOF(iz)) num2cell(full(Z(iz,iz))\speye(length(iz),1))]

 [u,s]=svd(full(Z));s=diag(s);u(abs(u)<1e-6)=0;
 u(abs(u)<1e-5)=0;%figure(1);semilogy(sort(abs(diag(Z))))
 i1=CE.DOF(any(u,2));fecom('shownodemark',i1)
 disp(feval(p_pml('@DofLabT'),unique(round(rem(i1,1)*100))/100))
 %feval(p_pml('@DofLabT'));
end

%% #ShearTransient

mo1=d_pml('MeshUlb1D SW',struct('subtype',1,'Lc',2,'quad',1,'Transient',1));
mo1=feutil('setpro 100 a0 2000 pow2 form -3000',mo1);
cf=feplot(mo1);[CE,un1,mo1.DOF]=fe_case(mo1,'gett');

dt=.1e-5; opt=struct('Method','Newmark','Opt',[.25 .5 0 dt .07/dt],'NeedUVA',[1 0 0]);
ofact spfmex;%ofact umfpack
d1=fe_time(opt,mo1);d1.TR=struct('def',CE.T,'DOF',mo1.DOF);
d_pml('View1Dtime',mo1,d1);

cf.def=fe_def('subdef',d1,1:100:size(d1.def,2));fecom(';showfimdef;ch5;scc5e10;view1');

% #ShearTransient_Unl  xxxEB : disp/stress formulation
mo1=d_pml('MeshUlb1D SW',struct('subtype',1,'Lc',2,'quad',1,'Transient',1));
mo1=feutil('setpro 100 a0 2000 pow2 form -3000',mo1);
dt=1e-5;opt=struct('Method','back','Opt',[.25 .5 0 dt .07/dt],'AssembleCall', ...
    'assemble -matdes 2 3 1 -load -fetime -exitfcn "vout=p_pml(''''assembleexit'''',model,Case,Load,RunOpt);"');
opt.Residual=nl_spring('residualcall');
%opt.Method='back';[SE,CE,op2,o]=fe_time(opt,mo1);
opt.RelTol=-1e-6;
opt.Method='NLNewmark';d1=fe_time(opt,mo1);
d_pml('View1Dtime',mo1,d1);


%%
elseif comstr(Cam,'v1')
%% #v1 : very basic functional script

RO.fmax=100; RO.Lp=10; RO.Lc=[5 1 1];RO.len=.5;[model,RO]=d_pml('MeshHBeam',RO);
RO.fmax=10; RO.Lp=3; RO.Lc=[RO.Lp 1 1];RO.len=1;RO.quad=0;[model,RO]=d_pml('MeshHBeam',RO);
%feplot(model);fecom colordatapro-edgealpha.1
freq=logspace(-4,2,100)';
[d1,SE]=p_pml('SolveDfrf',stack_set(model,'info','Freq',freq));feplot(SE,d1);fecom('colordataEvalA');iicom('curveinit','Test',fe_def('subDof',d1,RO.In));iicom(';submagpha;xlog');

%% Very small example impact of tuning alpha x
RO.fmax=10; RO.Lp=3; RO.Lc=[RO.Lp 1 1];RO.len=1;RO.quad=0;
RO.fmax=[logspace(1,3,3)';5000]/2/pi;RO.jdisp=3;
RO.freq=logspace(0,2,300)';RO.volRule=-2;
[model,RO]=d_pml('MeshHBeam',RO);

model.il(model.il(:,1)==101,3)=31;d_pml('SolveFmaxStudy',model,RO);% unidir qxyz

d_pml('SolveFmaxStudy',model,RO);% Nominal
model.il(model.il(:,1)==101,3)=11;d_pml('SolveFmaxStudy',model,RO);% fix edge

z=load('D:\balmes\Dropbox\sdtdata\collab\hadrien\Avancement\PML_test.mat');
C2=fe_def('subchCurve',z.iiplot_data.Stack{3,3},{'Capt',9});
%C2.Y=repmat(C2.Y,1,size(C1.Y,2));C2.X{2}=(1:size(C2.Y,2))';
C2.Y(:,2)=z.iiplot_data.Stack{2,3}.Y(:,9);
C2.Y(:,3)=z.iiplot_data.Stack{4,3}.Y(:,9);
C2.Y(:,4)=z.iiplot_data.Stack{5,3}.Y(:,9);
cj=iiplot(14);iicom(cj,'curveinit','Hadrien',C2);iicom(cj,';submagpha;xlim2 80;chall;xlog')

 [d1,SE]=p_pml('SolveDfrf',stack_set(model,'info','Freq',freq));r2=fe_def('subDof',d1,RO.In+.01);
 C1.X{2}=RO.fmax(1:jpar);C1.Y(:,jpar)=r2.def(:);C1.Ylab='Resp';
iicom('curveinit','Test',C1);iicom(';ch7;xlim5 30');get(ci.ax(1,4),'ylim')


else; error('Script%s unknown',CAM);
    
end
%% #Mesh -------------------------------------------------------------------
elseif comstr(Cam,'mesh');[CAM,Cam]=comstr(CAM,5);


if comstr(Cam,'hbeam');[CAM,Cam]=comstr(CAM,6);
%% #MeshHbeam : horizontal beam test

if ~isfield(RO,'len');RO.len=1;end
if ~isfield(RO,'fmax');RO.fmax=100;end
if ~isfield(RO,'Form');RO.Form=0;end
if ~isfield(RO,'pow');RO.pow=1;end
if ~isfield(RO,'Lp')||length(RO.Lp)==1;
 %% Manual meshing
 model=feutil('objecthexa 1 1',[0 0 0;eye(3)], ...
    feutil(sprintf('refineline %.15g',RO.len),[0 RO.Lc(1) RO.Lc(1)+RO.Lp]), ...
    feutil(sprintf('refineline %.15g',RO.len),[0 RO.Lc(2)]), ...
    feutil(sprintf('refineline %.15g',RO.len),[0 RO.Lc(3)]));
 model.Elt(feutil('findelt innode{x>=}',model,RO.Lc(1)),10)=101;%PML
 model.pl=[1 fe_mat('m_elastic','SI',1) 6.8777e+07 4.0476e-01 1700 0];
 model=p_solid('default;',model);
 model=feutil('setpro',model,[101 fe_mat('p_pml','SI',1) RO.Form RO.pow RO.fmax(1)*2*pi RO.Lc(1) RO.Lp]); %qx
else
 %% Auto cube meshing RO.Lp=[3 1 0  0 0 0];
 if isnumeric(RO.Lc); 
   RO.Lc={feutil(sprintf('refineline %.15g',RO.len),[0 RO.Lc(1)])
         feutil(sprintf('refineline %.15g',RO.len),[0 RO.Lc(2)])
         feutil(sprintf('refineline %.15g',RO.len),[0 RO.Lc(2)])};
 end
 model=feutil('objecthexa 1 1',[0 0 0;eye(3)],RO.Lc{:});
 if ~isfield(RO,'pl')
  model.pl=[1 fe_mat('m_elastic','SI',1) 6.8777e+07 4.0476e-01 1700 0];
 else
  model.pl=[];
  for j1=1:size(RO.pl);
   model.pl(end+1,1:length(RO.pl{j1,2}))=RO.pl{j1,2};
   model.Elt(feutil(['findelt' RO.pl{j1,1}],model),9:10)=model.pl(end,1);
  end
 end
 model.Elt=feutil('orient;',model); 
 model=p_solid('default;',model); model.unit='SI';
 RA=struct('Form',RO.Form,'Lp',RO.Lp, ...
     'pow',RO.pow,'a0',RO.fmax(1)*2*pi,'Lc',RO.len);
 RA=sdth.sfield('AddSelected',RA,RO,{'CycBuild'});
 if isfield(RO,'PML')&&~strcmpi(RO.PML,'none');
  RA.MeshType=RO.PML; model=p_pml('addPML',model,RA);
 elseif ~isfield(RO,'PML');model=p_pml('addPML',model,RA);
 elseif isfield(RO,'CycBuild'); model=fe_cyclic(RO.CycBuild,model);
 end
end
if ~isfield(RO,'quad')||RO.quad==0;RO.quad=0;
elseif RO.quad==5;    model=q16p('h8Toh125',model);
else; model=feutil('lin2quad',model);
end

% Integration rule at node
if isfield(RO,'volRule');model.il(model.il(:,1)==1,4)=RO.volRule;end
model.Elt=feutilb('SeparatebyProp',model);

d1=feutil('findnode x==0',model);RO.In=d1(1);
Load=[];
if RO.quad~=5
 data=struct('sel','x==0','def',1,'DOF',.19);
 Load=fe_load(fe_case(model,'Fsurf','Left',data));
 if ~any(Load.def); Load=[];end
end
if isempty(Load);
 Load=struct('def',ones(size(d1)),'DOF',d1+.01);
end
model=fe_case(model,'DofLoad','Left',Load);
if isfield(RO,'freq');model=stack_set(model,'info','Freq',RO.freq);end

if isfield(RO,'PML')&&~strcmpi(RO.PML,'none'); % New element implementation
 i1=feutil('getdof',feutil('findnode x==',model,max(model.Node(:,5))),(13:21)'/100); 
 model=fe_case(model,'FixDof','EndPml',sprintf('x==%g -DOF13',max(model.Node(:,5))));
 model=p_pml('DispMpc',model);
 model=fe_case(model,'pcond','PML','p_pml(''Pcond 1e5'')');
 model=stack_set(model,'info','oProp',{'method','umfpack'});
end
model.name='HBeam';
out=model;out1=RO;

elseif comstr(Cam,'ulb1d');[CAM,Cam]=comstr(CAM,6);
%% #MeshUlb1D : sample ULB test for S and P waves
if any(Cam=='{') % struct('v',1,'quad',1,'Lc',2)
  %[CAM,RO]=sdtm.urnPar(CAM,struct('cst',{{'v','%g';'quad','%3';'Lc','%g'}}),RO);
  [CAM,RO]=sdtm.urnPar(CAM,'{}{v%ug,quad%3,Lc%ug}',RO);
end

[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
 'dim([2 2 40]#%g#"Lx Ly Lz of standard material")' ...
 'pow(2#%g#"defaut attenuation")' ...
 'L(20#%g#"PML length")' ...
 'Lc(1#%g#" default mesh size")' ...
 'quad(0#3#" if non zero use lin2quad")' ...
 'SW(0#3#" shear wave case otherwise PW")' ...
 'subtype(2#%g#"1 time domain, 2 frequency domain")' ...
 'EdgeMpc(#3#" constrain equal motion on corresponding edges")' ...
 'Transient(#3#" initialize for transient")' ...
 ],{RO,CAM}); Cam=lower(CAM);

RO.lx=feutil(sprintf('refineline %g',RO.Lc),[0 RO.dim(1)]);
RO.ly=feutil(sprintf('refineline %g',RO.Lc),[0 RO.dim(2)]);
model=feutil('objecthexa 4 4',[-RO.dim(1:2)/2 0;eye(3)], ...
    RO.lx,RO.ly, ...
    feutil(sprintf('refineline %g',RO.Lc),[0 RO.dim(3)]));
model.pl=[4  fe_mat('m_elastic','SI',1) 30e9 0.2 2200 0 0.01 1e-5 20];
model=p_solid('default;',model);
model.unit='SI';
RP=struct('Lp',[0 0 0    0 0 RO.L],'pow',RO.pow,'a0',-1,'Lc',RO.Lc, ...
    'ProId',100-6,'subtype',RO.subtype); 
if isfield(RO,'quad')&&RO.quad;RP.quad=RO.quad;end
if isfield(RO,'Freq');RP.Freq=RO.Freq;end
if isfield(RO,'Np');RP.Np=RO.Np;end


if RO.subtype==2 % Frequency domain
 RP.a=struct('X',{{[10;115],{'fe';'fi'}}},'Y',[5 0;15 20]');
 %pro=struct('name','pml','il',100,'unit','SI','type','p_pml','a', ...
 %   );
 %mo1=stack_set(mo1,'pro','PML',pro);
end
mo1=p_pml('addPML',model,RP);

% Reorder nodes to ease analysis
[n1,i1]=sortrows(mo1.Node,[5 6 7]);
mo1=feutil('renumber',mo1,[n1(:,1) (1:size(mo1.Node,1))' ]);


%% Now define the frequency evolution for fe and fi (frequency domain)

n1 = feutil('findnode z==0',mo1); % Bottom nodes
% Nodes on central line
if rem(length(RO.lx),2)&&rem(length(RO.ly),2)
 n2=feutil('getnode y==0 & x==0',mo1);n2=sortrows(n2,7);
else % Nodes on edge line
 n2=feutil('getnode y== & x==',mo1,min(mo1.Node(:,5)), ...
     min(mo1.Node(:,6)));
 n2=sortrows(n2,7);
end
mo1=fe_caseg('connectionequaldof -DOF 1 2 3 -safe',mo1,'Base', ...
    repmat(n1(1),size(n1,1)-1,1),n1(2:end,1));

if RO.SW 
%% #Shear_wave boundary conditions and load -3

 if RO.subtype==1% Time
  mo1 = fe_case(mo1,'FixDof','TopEdge',[sprintf('z==%f',max(n2(:,7))) '-DOF 1 2 3']);
  mo1=fe_case(mo1,'pcond','PMLp','p_pml(''Pcond 1e8'')');
  mo1=stack_set(mo1,'info','oProp',{'method','umfpack'});
  mo1=fe_case(mo1,'FixDof','EdgePml',['inelt {proid100 &selface}  -DOF' ...
    ' 13 16 ' ... % qx qy =0 
    ' 22 23 24 25 26' ... % sigma.._x
    ' 27 28 29 30 31 ' ... % sigma.._y
    ' 32 33 34 35' % keep Szx_z
    ]);

 else % Frequency
  mo1 = fe_case(mo1,'FixDof','TopEdge',sprintf('z==%f',max(n2(:,7))));
 end
 
 mo1 = fe_case(mo1,'FixDof','Shear','z>0 -DOF 2 3');
 %mo1 = fe_case(mo1,'FixDof','NoY',sprintf('y==%f | y==%f -DOF 2', ...
 %    min(mo1.Node(:,6)),max(mo1.Node(:,6))));

 mo1=fe_case(mo1,'DOFSet','ShearDisp',struct('def',1,'DOF',n1(1)+.01));
 % Put sensors on central line
 mo1=fe_case(mo1,'SensDof','Out',n2(:,1)+.01);
 
 RO.EdgeS={'connectionequaldof -safe-DOF  19 36';
           'connectionequaldof -safe-DOF  1'};
   
else;
%% #PWave_wave boundary conditions and load -3
 if RO.subtype==1% Time
  mo1 = fe_case(mo1,'FixDof','TopEdge',[sprintf('z==%f',max(n2(:,7))) '-DOF 1 2 3']);
  mo1=fe_case(mo1,'pcond','PMLp','p_pml(''Pcond 1e8'')');
  mo1=stack_set(mo1,'info','oProp',{'method','umfpack'});
  %feval(p_pml('@DofLabT'));
  mo1=fe_case(mo1,'FixDof','EdgePml',['inelt {proid100 &selface}  -DOF' ...
    ' 22 23 24 25 26' ... % sigma.._x
    ' 27 28 29 30 31 ' ... % sigma.._y
    ' 32 33 35 36' % keep Szz_z
    ]);
  mo1=fe_case(mo1,'FixDof','Compression','z>0 -DOF 1 2 13 14 15 16 17 18 19 20');%u qx_z
 else % Frequency
  mo1 = fe_case(mo1,'FixDof','TopEdge',sprintf('z==%f',max(n2(:,7))));
  mo1=fe_case(mo1,'FixDof','Compression','z>0 -DOF 1 2 ');%u qx_z
 end
  
 mo1=fe_case(mo1,'DOFSet','PDisp',struct('def',1,'DOF',n1(1)+.03));
 % Put sensors on central line
 mo1=fe_case(mo1,'SensDof','Out',n2(:,1)+.03);
 RO.EdgeS={'connectionequaldof -safe -DOF  21 34';
           'connectionequaldof -safe -DOF  3'};
end
if RO.Transient
 %% #Ulb1D_transient check time integration -3
 [C1,name]=fe_case(mo1,'getcase');
 i1=strcmpi(C1.Stack(:,1),'dofset');
 C1.Stack{i1,1}='DofLoad';
 C1.Stack{i1,3}.curve=fe_curve('test ricker dt=.01 A=1'); 
 % C1=ii_mmif('fft -struct -fmax 300',fe_curve('test ricker dt=.01 A=1',linspace(0,1,1e5)));iicom('curveinit','Ricker',C1)
 mo1=stack_set(mo1,'case',name,C1);
end
if isfield(RO,'EdgeMpc')&&RO.EdgeMpc 
 %% #EdgeMPC to have really 1D -3
 n1=feutil('getnode z>=40',mo1);n1=sortrows(n1,[7 5 6]);n1=reshape(n1(:,1),4,[])';
 mo1=fe_caseg(RO.EdgeS{1},mo1,'eqdof' ,[n1(:,1);n1(:,1);n1(:,1)],[n1(:,2);n1(:,3);n1(:,4)]);
 %mo1=fe_case(mo1,'fixdof','trans',n1(:)+.03);

 n1=feutil('getnode z>0&z<40',mo1);n1=sortrows(n1,[7 5 6]);n1=reshape(n1(:,1),4,[])';
 mo1=fe_caseg(RO.EdgeS{2},mo1,'eqdofb' ,[n1(:,1);n1(:,1);n1(:,1)],[n1(:,2);n1(:,3);n1(:,4)]);
end

out=mo1;
elseif comstr(Cam,'pml_1d');[CAM,Cam]=comstr(CAM,6);
%% #MeshPML_1D : 1D wave propagation with PML
sdtw('_ewt','obsolete');
RO=struct('v',1,'n',4);

%% Create mesh
model=feutil('Objectquad 1 1',2*[-1 -1 0;1 -1 0;1 1 0;-1 1 0]*1e-4,2,2);
st=sprintf('Extrude 0 0 0 1'); model=feutil(st,model,linspace(0,RO.n*20,RO.n*10+1)*1e-4);
model= feutil('lin2quad',model);

%% PML Properties

%MelId = 5; 

posx = [-10 -5 5 10]*1e-4; % out of range, not used here
posy = [-10 -5 5 10]*1e-4; % out of range, not used here
posz = [-2 -1 RO.n*20-10 RO.n*20]*1e-4; % first two values out of range, not used here

nf = 5; freq = linspace(1,100,nf)*3e3;
fxe = linspace(0,0,nf); % x-dir
fye = linspace(0,0,nf); % y-dir
fze = linspace(0,0,nf); % z-dir
fxp = linspace(20,20,nf); % x-dir
fyp = linspace(20,20,nf); % y-dir
fzp = linspace(20,20,nf); % z-dir
m = 2; 

%% divide in 2 groups and build pl and il
model.Elt=feutil(sprintf('DivideGroup 1 withnode{z<=%f}',posz(3)-1e-6),model);
model.Elt=feutil('SetGroup 2 Mat 1 Pro 2',model);  

tan = 0.01; % damping
model.pl(1,:)=[1  fe_mat('m_elastic','SI',1) 30e9 0.2 2200 0 tan 1e-5 20];
model.il=p_solid('dbval d3 -3');

% xxxEB : need to add the declaration of PML using p_pml
% parameters are :
% 1) limits of the PML domain (posx posy posz when PML is a rectang domain)
% 2) fxe,fye,fze, fxp, fyp, fzp as a function of frequency (constant here)
% 3) m = order of the fe and fp functions
%

out=model;
    
%% MeshEnd
else;error('Mesh%s unknown',CAM);
end
%% #Solve -------------------------------------------------------------------
elseif comstr(Cam,'solve');[CAM,Cam]=comstr(CAM,6);

if comstr(Cam,'reducebase+f1')
%% #Base+f1 Enforce base acceleration with correction for first frequency -2

elseif comstr(Cam,'fmaxstudy')
%% #SolveFmaxStudy : adjust properties -2

model=RO; RO=varargin{carg};carg=carg+1;

for jpar=1:length(RO.fmax)
 model.il(2,4:5)=[2 RO.fmax(jpar)*2*pi]; %pow,a0
 if any(remi(model.il(2,3),[],1)==[0 4]); d1=fe_simul('dfrf',model);
 else;[d1,SE]=p_pml('SolveDfrf',model);% obsolete
 end
 r2=fe_def('subDof',d1,RO.In+.01);
 if jpar==1;C1=struct('X',{{r2.data}},'Xlab',{{'Freq','fmax'}},'Y',r2.def(:));end
 C1.X{2}=RO.fmax(1:jpar)*2*pi;
 %C1.Y(:,jpar)=r2.def(:).*(r2.data(:,1)*2*pi).^2;C1.Ylab='Acc';
 C1.Y(:,jpar)=r2.def(:);C1.Ylab='Disp';
 if isfield(RO,'jdisp')&&RO.jdisp==jpar
  cf=comgui('guifeplot -reset',2);cf.model=model;cf.def=d1;
  fecom('colordataEvalX');
 end
end
ci=comgui('guiiiplot -reset',2);
iicom(ci,'curveinit','Test',C1);iicom(ci,';submagpha;chall;xlogall');

elseif comstr(Cam,'perred')
%% #SolvePerRed : periodic reduction -2

%  fid=fopen('xxx.txt','w'); fprintf(fid,'%c',comstr(RO,-30,struct('NoClip',1)));fclose(fid);edit xxx.txt

SE=RO;if carg<=nargin;RO=varargin{carg};carg=carg+1;else;RO=struct;end

RR=struct('BuildList',{{'feval',{'p_pml','SolveDfrf', ...
 '@stack_set(SE3,{''info'',''Freq'',r1{4}.Freq;''info'',''oProp'',{''method'',''umfpack''}})', ...
struct('Range',struct('ncx',[5]),'Freq',[0.1])}}}, ...
'matdes',[2 1 3 4],'fe_coor','lrisvd', ...
'Phase2',{{'fe_homo','DftRedp2LU'}}, ...
'Assemble','assemble -matdes 2 1 3 4  -se NoT-cfield=1 Load -exitFcn"vout=dyn_solve(''''ReduceHetAssembly'''',model,Case,vout,out);"', ...
'P2Sets',{{ ...
 'solid',(1:3)'/100,struct('type','svd','tol',1e-8)
 %'pmld', (13:21)'/100,struct('type','svd','tol',1e-8)
 'pmlq', (13:21)'/100,struct('type','CaseT')
 ...%'pmls', (22:46)'/100,struct('type','svd','tol',1e-4,'coef',1e-6,'doZero',1e-5)
 'pmls', (22:46)'/100,struct('type','CaseT')
 }});
 [SE,dref]=fe_homo('dftRedList',SE,RR);
 SE=stack_set(SE,'curve','Learn',dref);
 data=fe_case(SE,'getdata','Symmetry');
 [i1,SE.Elt]=feutil('eltidfix;',SE);
 RO.NodeId0=1;RB.bas=1;
 mt=fesuper(sprintf('SEadd -trans %i %g 0 0 %i %i %i %s -bas %i',...
    RO.n_slice,data.trans(1),size(data.IntNodes,1),RO.NodeId0,100000,'s1', ...
    RB.bas),...
    struct('Node',[],'Elt',[]),SE,(1000)*[1 1]); 
 
 Load=fe_load(SE);mt=fe_case(mt,'DofLoad','In',Load);
 out=SE; out1=mt; 
 

else;error('Solve%s unknown',CAM);
end
%% #View -------------------------------------------------------------------
elseif comstr(Cam,'view');[CAM,Cam]=comstr(CAM,5);

if comstr(Cam,'1dps')
%% #view1DPS : P/S wave comparison with analytic result
if carg>nargin; eval(iigui({'mo1','d1'},'GetInCaller'))
else
 mo1=RO; d1=varargin{carg};carg=carg+1;
end
mat=feutil('getmat 4 -struct',mo1,4);

if ~isempty(fe_case(mo1,'getdata','Shear'))
 vw = sqrt(mat.E/2/(1+mat.Nu)/mat.Rho);st='Analytic Shear';% VW=v shear
else % Compression wave
 K=mat.E*(1-mat.Nu)/(1+mat.Nu)/(1-2*mat.Nu);
 vw=sqrt(K/mat.Rho);st='Analytic Pressure';
end
%[fn,Nodfn1] = feutil('findnode x==0 & y==0',mo1);  

sens=fe_case(mo1,'sens');
C1=fe_case('sensobserve',sens,d1);
if size(C1.Y,1)==1&&ischar(C1.Xlab{1})&&strcmpi(C1.Xlab{1},'In')
 i1=size(C1.Y);C1.Y=reshape(C1.Y,i1(2:end));
 C1.X(1)=[];C1.Xlab(1)=[];
end
NNode=sparse(mo1.Node(:,1),1,1:size(mo1.Node,1));
C1.X{2}=mo1.Node(NNode(sens.tdof(:,2)),7);
C1.Xlab{2}={'Position','m',[],''};
C1.Xlab{1}={'Freq','Hz',[],'%.0fHz'};
C1.Ylab={'U','m',[]};
ks = 2*pi/vw*d1.data(:,1)';
C1.Xlab{3}='Mdl';C1.X{3}={'PML';st};
C1.Y(:,:,2)=exp(-1i*ks(:)*C1.X{2}');C1.DimPos=[2 1 3];

ci=comgui('gui iiplot -reset;',2);
cingui('plotwd',ci,'@OsDic',{'ImToFigN','ImSw60','WrW49c','FNN'});
iicom(ci,'curveinit','C1',C1);iicom(ci,'showPiSubRi');
iicom('ch',{'Freq',1:size(C1.Y,1);'Mdl',1:2})

try;d_imw('LiL1M0S2',ci);end
% Possibly display shapes
[CAM,Cam,cf]=comstr('cf',[-25 31],CAM,Cam);
if ~isempty(cf)&&cf
    cf=comgui('guifeplot -reset;',cf);
    feplot(mo1,d1);fecom(cf,'ShowFiMdef');
end

elseif comstr(Cam,'1dtime')
%% #view1Dtime : validation of time response to Ricker
mo1=RO;
d1=varargin{carg};carg=carg+1;

sens=fe_case(mo1,'sens');
NNode=sparse(mo1.Node(:,1),1,1:size(mo1.Node,1));
step=find(d1.data(:,1)-d1.data(1)>.005,1,'first');
C1=fe_case('sensobserve',sens,fe_def('subdef',d1,1:step:size(d1.def,2)));
C1.X{2}=mo1.Node(NNode(sens.tdof(:,2)),7);
C1.PlotInfo=ii_plp('PlotInfo2D -type surface',C1);
C1=sdsetprop(C1,'PlotInfo','ua.YFcn','r3=r3;');
ci=comgui('guiiiplot -reset',2);
iicom(ci,'curveinit',C1);


else;error('View%s unknown',CAM);
end
%% End function
elseif comstr(Cam,'cvs');
 out=sdtcheck('revision');
 %out='$Revision: 529 $  $Date: 2020-11-02 14:25:17 +0100 (Mon, 02 Nov 2020) $'; return;
end




