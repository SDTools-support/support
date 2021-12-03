function [out,out1,out2]=t_stress(varargin) %#ok<*STOUT>

% Various stress/strain validation illustrations and tests 
% see details in automated tests

% Contributed by E. Balmes


if nargin==0; Cam='';
else; [CAM,Cam]=comstr(varargin{1},1);carg=2;
end
%#ok<*ASGLU,*STRNU,*NOSEM>

if isempty(Cam)  % These are tests that run automatically

t_stress('basic')
t_stress('expand')
% t_stress('map')  illustrate display of a principal stress MAP
% sdtweb('_tracker',224) resultant
t_stress('cuthexa') % handling of stress cuts
t_stress('cutgauss') % observation at Gauss points
t_stress('cutbeam') % handling of stress cuts
t_stress('cutUbeam') % various tests with ubeam model
t_stress('cutRivlin') % various tests with ubeam model

%% #basic - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'basic') 

demosdt nostop; 
%ofutil of_mk
model = femesh('test ubeam plot');model.Elt=feutil('orient',model);
model=fe_case(model,'fixdof','clamped edge','z==0');

def=fe_eig(model,[105 10 1e3]);%cf.def=def; 

model.Elt(1,1:8)=[Inf abs('hexa8') 0 0];
tic;a=fe_stress('stress',model,def);toc
model.Elt(1,1:8)=[Inf abs('hexa8b') 0];
b=fe_stress('stress atnode',model,def);toc

if norm(a.data-b.data,'inf')/norm(a.data,'inf')>1e-10; error('Not consistent');end

% base functions of fe_stress
fe_stress tensortopology
[a,b]=fe_stress('tensortopologyMeca3D');
 
% simple tests on volume stress - - - - - - - - - - - - - - - - - - - -

C1=fe_stress('stress-gstate',model,def);C1=C1.GroupInfo{5};
r2=feval(fe_stress('@Principal'),41,C1.Y,[],a);
r3=feval(fe_stress('@Principal'),4,C1.Y,[],a);
% r3=sqrt(.5*((r1(1,:)-r1(2,:)).^2+(r1(2,:)-r1(3,:)).^2+(r1(3,:)-r1(1,:)).^2+6*sum(r1(4:end,:).^2,1)));
if norm(r2(:)./r3(:)-1)>1e-10; error('Mismatch');end
if 1==2
   r1=rand(6,1e6);
   tic;r2=feval(fe_stress('@Principal'),4,r1,[],a);toc
   tic;r2=feval(fe_stress('@Principal'),41,r1,[],a);toc
end
%% #expand Tests of stress reexpansion to nodes - - - - - - - - - - - - - - - -
elseif comstr(Cam,'expand') 

d_ubeam
model.Elt(1,1:8)=[Inf abs('hexa8b') 0];

Case=fe_stress('stress atinteg',model,cf.def);

en=fe_stress('expand',model,Case);
fecom('colordata node',en)



%% --------------------------------------------------------------------
% display of principal stress map
elseif comstr(Cam,'testhexa') 

if carg<=nargin&&isstruct(varargin{carg});RO=varargin{carg};carg=carg+1; %#ok<*NASGU>
else; RO=struct;end
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   'scale([1 1 1]#%g#"scale") ' 'quad(0#3#"quadratic elements") ' ],{RO,CAM});

model=femesh('testhexa8b divide 2 2');
data=struct('sel','groupall','dir',{{'x','y'}}, ...
    'DOF',[.01;.02]);  
model=fe_case(model,'dofset','Force',data,'fixdof','base','z==0');

model.Node(:,5:7)=model.Node(:,5:7)*diag(RO.scale);
if RO.quad;model=feutil('lin2quad',model);end


def=fe_simul('static',model); %fe_c(q.DOF);ans(:,2)=num2cell(q.def)
out=model; out1=def;

% --------------------------------------------------------------------
%% #ShowMap display of principal stress map
elseif comstr(Cam,'showmap') 


[model,def]=t_stress('testhexa');

% A MAP with only the principal stress component
feplot(model);
MAP=fe_stress('stress -gstate -post "fe_caseg(''stressprincipalmap'')"',model,def);model.DOF=def.DOF;
MAP.arProp={'edgecolor','flat'};
MAP.color=sqrt(sum(MAP.normal.^2,2));
fecom('showmap',MAP)

% A MAP with two components at each integration point
MAP2=struct('vertex',[MAP.vertex;MAP.vertex], ...
    'normal',[MAP.p2;MAP.p3],'color',[sqrt(sum(MAP.p2.^2,2));
    sqrt(sum(MAP.p3.^2,2))],'arProp',{{'edgecolor','flat'}});
MAP2.normal=MAP2.normal/max(abs(MAP2.color)); % renorm
fecom('showmap',MAP2)

% --------------------------------------------------------------------
%% #cut display of stress cut plane
elseif comstr(Cam,'cut');[CAM,Cam]=comstr(CAM,4); 

if comstr(Cam,'hexa') % #CutHexa
%% test with hexa element   
[model,def]=t_stress('testhexa'); cf=feplot(model);
C1=fe_stress('stress-gstate',model,def); C1=C1.GroupInfo{5};

% display as a separate mesh
RO=struct('Origin',[0 1 .1],'axis',[0 1 1]);
mo3=fe_caseg('stresscut',RO,model); % cut mesh
%mo3=model;mo3.Elt=model.Elt(1:3,:);mo3.Node=feutil('getnodegroupall',mo3);
cut=fe_caseg('stresscut -selout',mo3,cf);
cut.StressObs.CritFcn='';
R1=fe_caseg('stressobserve',cut,def);
r2=[mean(R1.Y,2) mean(reshape(C1.Y,6,[]),2)];
if norm(r2*[-1;1])/norm(r2)>1e-10;error('Mismatch');end

r2=struct('StressObs',cut.StressObs.trans,'fvcs',0);
r2.StressObs.Node=cut.StressObs.Node;r2.StressObs.CritFcn='';
rb=feutil('placeindof',r2.StressObs.DOF,feutilb('geomrb',model,[0 0 0]));
R1=fe_caseg('stressobserve',r2,rb);R1=reshape(permute(R1.Y,[2 1 3]),[],6);
cf.def=rb;

rb3=feutilb('placeindof',feutil('getdof',[.01;.02;.03],mo3.Node(:,1)), ...
 feutilb('geomrb',mo3,[0 0 0]));
rb3.def-R1;
iz=reshape(1:45,3,[])';
reshape(cf.SelF{2}.cna{1}(iz,:)*rb.def,[],6)

cf.SelF{2}.cna{1}=cut.StressObs.trans.cta(iz,:);feplot

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Test with a detailed beam
if 1==1 % degre 1
  model=femesh('testhexa8b divide 10  1 1');
elseif 1==2 % degre2 integrules('hexa20')
 model=femesh('testhexa20 divide 10  1 1');
 model.il=[];model=p_solid('default',model);
end

model.Node(:,6:7)=model.Node(:,6:7)/10;
model.Elt=model.Elt([1 5:end 2:4],:);
model=fe_case(model,'fixdof','base','x==0','fixdof','2D',.02);
def=fe_eig(model,[5 20 0]);

cf=feplot(model);cf.def=def;

n1={[1 0 0 0 0 .05 .05;2 0 0 0 1 .05 .05];
    [1 0 0 0 0 .05 .1;2 0 0 0 1 .05 .1];
    [1 0 0 0 0 .05 .0;2 0 0 0 1 .05 .0]};
mo2=[];RO.div=100;
for j1=1:length(n1)
 mo1=struct('Node',n1{j1});
 mo1.Elt=feutil(['objectbeamline ' sprintf(' %i',n1{j1}(:,1))]);
 mo1=feutil(sprintf('divide %i',RO.div),mo1);
 if isempty(mo2);mo2=mo1;
 else;mo2=feutil('addtest',feutil('rmfield',mo2,'Stack'),mo1);
 end
end
mo1=mo2;

fe_caseg('stresscut',mo1,cf);cla;feplot % generation observation
fecom('fps10')

% Test with cuts in variables
cut=fe_caseg('stresscut -SelOut',mo1,cf);
cutd=cut;st={'cta','DOF','X','Xlab'};
for j1=1:length(st)
 cutd.StressObs.(st{j1})=cutd.StressObs.trans.(st{j1});
end
cutd.StressObs.CritFcn=[];
R1=fe_caseg('stressobserve',cutd,def);

fecom scc.1;fecom('scalecolorone');
cf.SelF{2}.StressObs.CritFcn= ... % 11 garde sigma_xx
     '[r1,dir,lab]=feval(fe_stress(''@Principal''),11,r1,[],''Mecha3D'');'; 
r1=fe_caseg('stressobserve',cf.sel(2),cf.def); % pas echelle
chc=cf.ua.ob(1,4); % channel affiche dans feplot
scale=cf.ua.ob(1,7);
r2=get(cf.o(3),'FaceVertexCData');% couleur feplot
% xxx bug inversion echelle de couleur 
if sdtdef('isinteractive')
 figure(1);
 h=plot(r1.X{2}(:,5),[squeeze(r1.Y(1,:,chc))' r2],'+');
 set(h(2),'marker','o')
else
 r1=norm([squeeze(r1.Y(1,:,chc))'-r2]); if norm(r1)>eps;error('mismatch');end %#ok<*NBRAK>
end
if 1==2 % When color was interpolated
 cg=feplot(10);
 mo2=cf.sel(2).StressObs.IntMesh;[eltid,mo2.Elt]=feutil('eltidfix;',mo2);
 a=fe_quality('meas jacobian 2',mo2);
 st=['eltid ' sprintf('%i ',a.EltId{1}(a.data{1}<1e-15))];
  
 mo2.Node(:,6:7)=mo2.Node(:,6:7)*100;cg.mdl=mo2;
 cg.sel={st,'colordatamat -alpha 1 -edgealpha .1'};fecom colorfacew
 %iimouse cv
 iimouse('view',gca,[ -2.414 -3.748 1.605 0.5 0.05 0.05 0.19 0.25 0.95 9.41]);
 % comgui imwriteplots/2D_interpmesh.png
end

elseif comstr(Cam,'gauss') 
%% #CutGauss test with hexa element   
[model,def]=t_stress('testhexa -scale 1 1 .3 -quad'); 
cf=feplot(model,def);

cut=fe_caseg('StressCut -selout',struct('type','Gauss'),model);
fe_caseg('stresscut',cut,cf)  % displays field at nodes

if 1==2  % Display principal stress
 model.il=p_solid('dbval 111 d3 -1');
 cut=fe_caseg('StressCut -selout',struct('type','Gauss'),model);
 cut.vert0=repmat(cut.vert0,3,1);
 cut=fe_sens('SelArrowArrow',cut);
 cut.dir=fegui('@PrincipalStress');
 cut.fvcs=[];cut.f1=[]; cut.DefLen=.3;cut.opt(1,6)=2;
 fe_caseg('stresscut',cut,cf)
 fecom(cf,'SetProp sel(1).fsProp','FaceAlpha',0.1,'EdgeAlpha',0.1);

 sdtweb feplot('ResultantDir')
   
end

% verify ability to rebuild stiffness matrix ____disassembly____
mo1=cf.mdl.GetData;
r1=cut.StressObs;
[mo1,C1]=fe_case('assemble -matdes 2 1 -se NoT',mo1);
if norm(feval(elem0('@get_lambda'),C1,1,1)-r1.Lambda{1})>1e-10
    error('Mismatch');
end
EC=C1.GroupInfo{8};EC.w(:,4);
% should use ddg=vhandle.matrix.stressCutDDG(obs,RO);
% mkl_utils(struct('jac',xxx))
i1=size(r1.Lambda{1},1);
RM.ddg=vhandle.matrix.stressCutDDG(struct('alloc',[i1 size(r1.X{2},1)],'ddg',ones(i1)));

lambda=inv(r1.Lambda{1});i2=size(lambda,1);i1=size(r1.cta,1)/i2;
in1=repmat(reshape(1:size(r1.cta,1),i2,[]),i2,1);
in2=repmat(1:size(r1.cta,1),i2,1);
lambda=sparse(in1(:),in2(:),repmat(lambda(:),1,i1)*diag(r1.wjdet));
k1=r1.cta'*(lambda)*r1.cta; % cta is stress obs here 
if normest(k1-mo1.K{2}) /normest(k1)>1e-10; 
  full([diag(k1)./diag(mo1.K{2})]);figure(1);plot(ans)
  error('Mismatch');
end
%% Check NL implementation in uMax and jacobian calls
% sdtweb _textag 'For the representation of bushings'
% sdtweb nlutil umaxlo
%ja=mkl_utils(struct('jac',1,'simo',NL));dd=reshape(ja(1:81,1),9,9);

%% surface stress at nodes
mo1=fe_caseg('StressCut -selout',struct('type','conform','sel','selface & innode{z==0}'),model);
cut=fe_caseg('StressCut -selout',mo1,model);cut=cut.StressObs;
fecom('shownodemark',cut.Node(:,5:7))




elseif comstr(Cam,'ubeam'); 
% #CutUbeam #StressCutUbeam case with analytic expression of strain for verification
model=demosdt('demoubeam -noplot');
model.pl=m_elastic('dbval 1 Strain');
model.DOF=feutil('getdof',model);

data=struct('sel','groupall','dir',{{'1*x','y.^2','z.^2'}}, ...
    'DOF',[.01;.02;.03]);
def=elem0('VectFromDirAtDof',model,data,model.DOF);
if ~isempty(strfind(Cam,'back'));out=model;out1=def; return; end
    
RO=struct('Origin',[0 0 .5],'axis',[0 0 1]);
%sdtweb VectFromDir
[mo3,d3]=fe_caseg('StressCut',RO,model,def);
cg=feplot(10);feplot(cg,mo3,d3);fecom colordata19

cf=feplot(2);feplot(model,def);
fe_caseg('StressCut',RO,cf.mdl,fe_def('subdef',def,1));

% Stress on a line given by two points and steps
% vertical line
RO=struct('Origin',[-.45 -.45 0;-.45 -.45 1],'steps',linspace(0,1,20));
[mo3,o3]=fe_caseg('StressCut',RO,cf.mdl,fe_def('subdef',def,1));
figure(1);d1=feutilb('placeindof',o3.DOF,def);
r1=reshape(o3.cta*d1.def,6,[])';plot(r1(:,1:3));legend(o3.X{1}{1:3})

% horizontal line
RO=struct('Origin',[-.45 -.5 .125;-.45 .5 .125],'steps',linspace(0,1,20));
[mo3,o3]=fe_caseg('StressCut',RO,cf.mdl,fe_def('subdef',def,1));
figure(1);d1=feutilb('placeindof',o3.DOF,def);
r1=reshape(o3.cta*d1.def,6,[])';plot(r1(:,1:3));legend(o3.X{1}{1:3})
%figure(1);plot(d3.def(:,1:3));legend(d3.lab{1:3})

cut=fe_caseg('StressCut-selout',RO,cf.mdl);
data=fe_caseg('StressObserve',cut,def);


%% Generate a complex model with mixed elements
demosdt('demoubeam'); cf=feplot;
[cf.mdl.Elt,elt]=feutil('removeelt innode{x<=0}',cf.mdl);
mo1=struct('Node',cf.mdl.Node,'Elt',elt);
mo1=feutil('hexa2penta',mo1);
cf.mdl=feutil('addtestmerge',cf.mdl,mo1);
mo1=cf.mdl.GetData;
mo1=fe_case('reset',mo1);
mo1=feutil('renumber-noori',mo1,max(mo1.Node(:,1))+1);
mo1.Node=basis('gnode','tz=2.5;',mo1.Node);
r1=feutil(sprintf('getnode z==%.15g',max(cf.mdl.Node(:,7))),cf.mdl);
[n1,i1]=feutil('addnode;',mo1.Node,r1);
mo1.Node=basis('gnode','tz=.5;',mo1.Node);
cf.mdl=feutil('combinemodel -compatnodeelt-compatmatpro',cf.mdl,mo1);
r1=[r1(:,1) mo1.Node(i1,1) 123*ones(length(i1),1)];
cf.mdl=feutil('addelt',cf.mdl,'rigid',r1);
%r2=r1(:,
cf.sel='reset';

% attempt to use stresscut
fecom('showpatch')
%cf.mdl=feutil('lin2quad',cf.mdl); % xxx invalid with rigid better stress interpolation
cf.mdl=fe_case(cf.mdl,'fixdof','clamprot',[4:6]'/100);
def=fe_eig(cf.mdl,[5 10 1e3]);
cf.def=def;

r1=struct('EltSel','withnode {z==3} & eltname ~=rigid', ...
    'SurfSel','inelt{innode{z==3.25}}', ...
    'type','Resultant');
sel=fe_caseg('stresscut -selout',r1,cf.mdl.GetData);
sel.StressObs.CritFcn='r1=r1(3,:);'; % Observe Z
cf.Stack{'sel','Cut1'}=sel;

fe_caseg('stresscut',sel,cf) % Display

% adapt transparencies
fecom(cf,'SetProp sel(1).fsProp','FaceAlpha',0.01,'EdgeAlpha',0.2);
fecom('colorbar', ...
    'units','normalized','position',[.88 .5 .04 .4], ...
    'YAxisLocation','left','FontSize',12, ...
    '@title',{'String','unit','FontSize',14})


elseif comstr(Cam,'rivlin'); 
%% #CutRivlin verify units of resultant cut sensor

model=femesh('testhexa8 divide 5 5 2');
%model=femesh('testhexa8 divide 1 1 1');
model=fe_mknl(model,'NoT');
[u1,u2,u3,dd]=p_solid('buildconstit 3 1',[100;111],model.pl,model.il);

epsi=dd.dd\[1;2;3;0;0;0]; % Epsilon for target sigma (must be diagonal)
% when applying a load
data=struct('sel','groupall','dir',{{sprintf('%.15g*x',epsi(1)), ...
    sprintf('%.15g*y',epsi(2)),sprintf('%.15g*z',epsi(3))}}, ...
    'DOF',[.01;.02;.03]);
rb=elem0('VectFromDirAtDof',model,data,model.DOF);
cf=feplot(model,rb);
cut=fevisco('feplotStress v=3 MatId 100 DefLen .2',model)
cut=fevisco('feplotStress v=3 MatId 100 DefLen .2 -rule -1',model)
C1=fe_caseg('stressobserve -crit""',cut,rb);C1.Y

r1=struct('EltSel','groupall', 'SurfSel','inelt{innode{z==1}}', ...
    'type','Resultant');
cut=fe_caseg('stresscut -selout',r1,model);

C1=fe_caseg('stressobserve -crit""',cut,rb);C1.Y

mean(abs(C1.Y),2)
reshape(cut.StressObs.cta{2}*rb.def,3,[])

% Visualize : resultant
d1=feutilb('placeindof',cut.StressObs.aDOF, ...
    struct('def',model.K{2}*rb.def,'DOF',model.DOF));
i1=feutil('findnode inelt {seledge}',cut.StressObs);
d1.def(fe_c(d1.DOF,[i1;.03],'ind'),:)=0;% noz, or edge : response is zero
d1.def(:,2)=cut.StressObs.cta{2}*rb.def;cf.def=d1

d1.def=cut.StressObs.cta{2}*rb.def;
i1=feutil('findnode inelt {seledge}',cut.StressObs);
d1.def(fe_c(d1.DOF,i1,'ind',2),:)=0;% on edge
d1.def(fe_c(d1.DOF,[.02;.03],'ind'),:)=0;% on edge
d1.def=cut.StressObs.cta{1}*d1.def;cf.def=d1;
cf.sel={'selface & innode {z==1}','colordataEvalX'};

% Visualize : resultant density
cf.def={cut.StressObs.cta{1}*cut.StressObs.cta{2}*rb.def,cut.StressObs.aDOF}
fecom colordataevaly-alpha1
 
% Display resultant field : illustrates that the surface load generating
% identical resultants does not correspond exactly to the stress field
%  for the in-plane components one sees a propagation effect 
%  with some oscillations. 
r1.DefLen=.3; r1.Field=1;cf.def=rb;cut=fevisco('feplotResultant',cf,r1);

C1=fe_caseg('stressobserve -crit""',cut,rb);
%cut.StressObs.CritFcn='r1=r1(1,:);'; % Observe X
%fe_caseg('stresscut',cut,cf)


elseif comstr(Cam,'shell'); 
%% #CutShell #StressCutShell : shell element stress/resultant

model=femesh('testquad4 divide 2 2');
%feplot(model);
Test=feutil('objectbeam 1 1',[.1 .1 0;.9 .9 0],5);
% sdtweb elem0 gaussobserve
cut=fe_caseg('StressCut -selout -radius 3',Test,model);


elseif comstr(Cam,'beam'); 
%% #Cutbeam #StressCutBeam : beam element stress/resultant
% Report from Sami Daouk page 30, formulas N,M, and page 32 for s_xx
% Cut in a beam

model=femesh('testbeam1 -divide 10');model.Node(:,5)=model.Node(:,5).^2;
model.Elt(1,1:7)=[Inf abs('beam1t')];
model.Elt(2:end,6)=1; % Y first bending plane
model.Node=feutil('getnode groupall',model);
model=fe_case(model,'fixdof','base','x==0','fixdof','2d',[.02;.06]);
def=fe_eig(model,[5 10 1e3]);
cf=feplot(model,def);
q0=fe_def('subdef',def,1);def=fe_def('subdef',def,1:3);

% low level
[Case,model.DOF]=fe_mknl('init',model);
re=fe_caseg('stressObserve 1 -jgroup 1 -rule -1 5',model,Case,q0);
figure(1);subplot(211);plot(re.Node(:,1),re.s_z*def.def)
figure(1);subplot(212);plot(re.Node(:,1),re.k_z*def.def)

% newer integration 
mo1=feutil('refinebeam .02',model);
cut=fe_caseg('StressCut-SelOut',mo1,cf.mdl);cut.StressObs.CritFcn='';
R1=fe_caseg('stressobserve ',cut,def);
figure(11);plot(R1.X{2}(:,5),squeeze(R1.Y(6,:,:)),'+')

% Generate a mesh at beam Gauss points
mo1=fe_caseg('StressCut',struct('type','BeamGauss'),cf.mdl);
mo1.type='ObserveAtGauss';mo1.Rule=-3;
cut=fe_caseg('StressCut-SelOut',mo1,cf.mdl);
cut.StressObs.CritFcn='';
R1=fe_caseg('stressobserve ',cut,def);
figure(11);plot(R1.X{2}(:,5),squeeze(R1.Y(6,:,:)),':+')

% Verify invariance for axial load
model=fe_case(model,'DofLoad','Tip',struct('def',10,'DOF',2.01));
q0=fe_simul('static',model);
EA=model.pl(model.pl(:,1)==100,3)*model.il(model.il(:,1)==112,6);
re=fe_caseg('stressObserve 1 -jgroup 1 -rule -1 5',model,Case,q0);
figure(1);clf;plot(re.Node(:,1),EA*(re.e_xx*q0.def))
R1=fe_caseg('stressobserve ',cut,q0);
figure(11);plot(R1.X{2}(:,5),EA*squeeze(R1.Y(1,:,:)),':+')


elseif comstr(Cam,'subcut')
%% #CutSubcut
%  unclean attempt at subselection
  r1=cuts.StressObs; r1.Elt=feutil('divideingroups',r1);
  r1.Elt=feutil('selelt group1',r1);
  i1=ismember(r1.Node(:,1),feutil('findnode groupall',r1));
  i2=reshape(1:size(r1.Node,1)*size(r1.X{1},1),[],size(r1.Node,1));
  i2=i2(:,i1);i2=i2(:);
  r1.cta{1}=r1.cta{1}(i2,i2);  r1.cta{2}=r1.cta{2}(i2,:);
  i2=reshape(1:size(r1.Node,1)*3,[],size(r1.Node,1));
  i2=i2(:,i1);i2=i2(:);
  r1.trans.cta=r1.trans.cta(i2,:);
  r1.Node=r1.Node(i1,:);
  cut2=cuts;cut2.StressObs=r1;
  cut2.vert0=cut2.vert0(i1,:);
  cf=feplot(r1,d1); fe_caseg('StressCut',cut2,cf);

else

 % A general interface for FEM post-processing
 % - generate the view geometry using section and
 %   line commands
 % - interpolate based wire frame element groups (multiple interp, not yet
 %   implemented)
 % - support observation generation and post-processing of specific curve

 cf=feplot(2);model=demosdt('demoubeam -noplot');cf.mdl=model;
 def=fe_eig(demosdt('demoubeam -noplot'),[5 20 1e3]);cf.def=def;
 
 RO=struct('Origin',[0 0 .05],'axis',[0 0 1],'planes',[0;.5;1]);
 mo3=fe_caseg('StressCut',RO,model); % generate wire-frame mesh
 
 [mo3,r3]=fe_caseg('stresscut',mo3,cf.mdl); % get output info
 
 % should be associated with specific .def handle
 
 fe_caseg('stresscut',mo3,cf.mdl); % display as color
 %sdtweb fe_caseg('stresscutcolor')
 cf.SelF{2}.StressObs.CritFcn= ...
     '[r1,dir,lab]=feval(fe_stress(''@Principal''),4,r1,[],''Mecha3D'');'; 
 fecom(cf,'colorbar','units','normalized','position',[.88 .5 .04 .4], ...
    'YAxisLocation','left','FontSize',12, ...
    '@title',{'String','Pa','FontSize',14})
 
 cg=feplot(5);feplot(cg,model,def);fecom colordatastressmises
 fecom animcoloroff;
 fecom(cg,'colorbar','units','normalized','position',[.88 .5 .04 .4], ...
    'YAxisLocation','left','FontSize',12, ...
    '@title',{'String','Pa','FontSize',14})
    
 % observe response at specific DOF
 sel=cf.SelF{2}; r1=fe_caseg('stressobserve',sel,def)
 
 figure(1);plot(cg.sel.fvcs(:,1))
 line(1:size(r1.Y,2),squeeze(r1.Y(1,:,1)),'color','r')
 
end

%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Stress computations in volume 
elseif comstr(Cam,'analyticstrain') %t_plate2('analyticStrain')
[CAM,Cam]=comstr(CAM,15);

if carg>nargin; model=femesh('testubeam');
else; model=varargin{carg};carg=carg+1;
end
model=feutil('lin2quad',model);
% analytic evaluation
if isempty(Cam) % z motion
   data=struct('sel','groupall','dir',{{0,0,'1*z.^2'}}, ...
    'DOF',[.01;.02;.03]);
  RunOpt.Type='z';
end

model.DOF=feutil('getdof',model);
def=elem0('VectFromDirAtDof',model,data,model.DOF);

mo1=model;
C1=fe_stress('stress -gstate',mo1,def);C1.jGroup=1;
[opt,match]=elem0('GaussObserve',C1.GroupInfo{C1.jGroup,[8 3 4]},mo1,C1);
r1=C1.GroupInfo{C1.jGroup,5};

r2=[match.Node reshape(r1.Y,6,[])'];
figure(1);plot(r2(1:4:end,3),r2(1:4:end,4:end),'+')

%% #Homo : test of elements using homogeneization
elseif comstr(Cam,'homo') 

st={'hexa8b','pyra5','pyra13'};
for j1=1:length(st)
 model=femesh(['teststruct' st{j1}]);model.Elt=feutil('set group1 matid 100',model.Elt);
 model.name='';
 [a,b]=fe_homo('rvekubc',fe_case(model,'reset')); 
 dd=feutil('getdd',[100 111 3 1],model);
 if norm(b.dd-dd.dd,'inf')/norm(dd.dd,'inf')>1e-9; error('mismatch');end
end




%% #Aniso test of anisotropic materials
% sdtweb feform#feelas3d
% sdtweb t_constit('elastic') % Basic test of anisotropic properties
elseif comstr(Cam,'aniso') 


model=femesh('testhexa8');
model=fe_case(model,'fixdof','Base',2:9);
%feplot(mo1);

C11=168.4e9; C12=121.4e9; C44=75.4e9; % GPa
C=[C11 C12 C12 0 0 0;C12 C11 C12 0 0 0;C12 C12 C11 0 0 0;
   0 0 0 C44 0 0;    0 0 0 0 C44 0;    0 0 0 0 0 C44]; 
bas=basis('bunge',[5.175 1.3071 4.2012]);
pl1=m_elastic('formulaPlAniso 100',C,bas');
pl2=m_elastic('formulaPlAniso 100',C,eye(3));

%pl1=m_elastic('dbval 100 steel');pl2=pl1;

% Assemble rotated material law
model.pl=pl1;SE=fe_case(model,'assemble -matdes 2 1 -SE'); 

% Rotate mesh thne assemble
mo1=model;mo1.Node(:,2)=1; mo1.bas=[1 1 0   0 0 0 bas(:)'];
[mo1.Node,bas1]=basis('nodebas',mo1);mo1.bas=[];
mo1.pl=pl2;
mo1=fe_case(mo1,'assemble -matdes 2 1 -SE');

k0=full(SE.K{2});k1=full(mo1.K{2});
r1=[k0./feutilb('tkt',bas,k1)]-1; if norm(r1)>1e-10; error('Mismatch');end

% Orient using a vector MAP (implemented for Matrix 5, VectFromDir)
p=bas';
data=struct('dir',{num2cell([p(:,1);p(:,2)]')}, ...
    'lab',{{'v1x','v1y','v1z','v2x','v2y','v2z'}});
 pro=struct('il',111,'type','p_solid','MAP',data);
 mo2=stack_set(model,'pro','WithMap',pro);mo2.pl=pl2;
 
 
 C2=fe_mknl('init',mo2);InfoAtNode=C2.GroupInfo{7};
 if ~isstruct(InfoAtNode);error('struct expected');end
 mo2.K={};mo2=fe_case(mo2,'assemble -matdes 2 1 5 -SE');
 k2=full(mo2.K{3});
 r1=[k0./k2]-1; if norm(r1)>1e-10; error('Mismatch');end

% MAP at elements see: sdtweb t_cyclic('calme12ViscMap')

% Orient using a material basis in p_solid (as done in NASTRAN)
% This is only valid for a fixed basis
mo2=model; mo2.pl=pl2; % non rotated material
mo2.bas=[100 1 0  0 0 0  reshape(bas',1,[])];
mo2=feutil('setpro111 coordm=100',mo2);
 
[co,in,ID,r1]=p_solid('buildconstit 3 1',[100 111],pl1,model.il);
[co,in,ID,r2]=p_solid('buildconstit 3 1',[100 111],pl2,mo2.il,mo2);
r3=sd('compare',r1,r2,struct('RelTol',1e-12));
if ~isempty(r3); disp(r3); error('Mismatch');end

 mo2.K={};mo2=fe_case(mo2,'assemble -matdes 2 1 5 -SE');
 k2=full(mo2.K{3});
 r1=[k0./k2]-1; if norm(r1)>1e-10; error('Mismatch');end

 % split material:
 model2=fevisco('matsplit 100',model);
 K=fe_case(model,'assemble -matdes 2 1 5 -reset -cell');
 K2=fe_case(model2,'assemble -matdes 2 1 5 -reset -cell');
 if norm(K{3}-K2{3},'inf')/norm(K{3},'inf')>1e-12; error('Splitting failed'); end

 if 1==2   
  [C1,mo2.DOF]=fe_mknl('initnocon',mo2);
  i1=[1000 1000*100+[1:6]]; RO.dd=cell(1,length(i1));
  for j1=1:length(i1)
   [co,in,ID,r1]=p_solid('buildconstit 3 1',i1(j1)*[1 1],mo2.pl,mo2.il,mo2,C1);
   RO.dd{j1}=r1.dd;
  end
  dd=zeros(6); for j1=2:7;dd=dd+RO.dd{j1};end
  [RO.dd{1}-dd]
 end
 
 % Generate map using rst (not for JP yet)
data=struct('dir',{{'v1x','v1y','v1z'}}, ...
    'lab',{{'v1x','v1y','v1z','v2x','v2y','v2z'}}, ...
    'FieldFcn','data=field_REST2V123(node,elt,data,EC)');

elt=model.Elt;[EGroup,nGroup]=getegroup(elt);jGroup=1;
cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
ElemF= getegroup(model.Elt(EGroup(jGroup),:),jGroup);
C1=fe_mknl('init',model);EC=C1.GroupInfo{jGroup,end};
elt=feutil('addelt',evalin('caller','ElemF'), ...
       evalin('caller','elt(cEGI,:)'));
r1=elem0('VectFromDirAtNode',model,data,EC)

% Split into sub-energies
%  generate multiple elements with 01, 02, 03, ... 06
% sdtweb fe_caseg('enersub') : low level strategy for optimized energy
% For easier handling, of parametric models property split is needed
%  -fevisco SplitMat

% - - - -
% Some simple tests to check EnerK computation : - - - - - - - - - - - - - - - 
% isotropic:
mdl=femesh('testhexa8b');
def=fe_eig(mdl,2);
EnerK=fe_stress('enerk',mdl,def);
%figure;plot(def.data,(def.data*2*pi).^2,def.data,2*sum(EnerK.data,1))
if norm((def.data*2*pi).^2-2*sum(EnerK.data,1)')>1e-5; 
  error('pb in enerk computation'); 
end
% #orthotropic:
mdl=femesh('testhexa8b');
%mdl.pl=aubefan08('matblade41'); mdl.pl(1)=100;
mdl.pl=[100,0.41164696590954,80250000000,53330000000,8390000000,0.410622,0.059400049968847,0.027688,4730000000,8600000000,7500000000,1550,0,0,0,0,0];
def=fe_eig(mdl,2);
EnerK=fe_stress('enerk',mdl,def);
%figure;plot(def.data,(def.data*2*pi).^2,def.data,2*sum(EnerK.data,1))
if norm((def.data*2*pi).^2-2*sum(EnerK.data,1)')>1e-5; 
  error('pb in enerk computation'); 
end
% corresponding full anisotropic:
[co,in,ID,r1]=p_solid('buildconstit 3 1',[100 111],mdl.pl,mdl.il);
i1=[1  7 8  13 14 15   19 20 21 22   25 26 27 28 29  31 32 33 34 35 36]; %+36*(j1-1);
mdl.pl=[100 fe_mat('m_elastic','SI',3) r1.dd(i1) fe_mat('GetMat 100 rho',mdl)];
def=fe_eig(mdl,2); def0=def;
EnerK=fe_stress('enerk groupall',mdl,def);
%figure;plot(def.data,(def.data*2*pi).^2,def.data,2*sum(EnerK.data,1))
if norm((def.data*2*pi).^2-2*sum(EnerK.data,1)')>1e-5; 
  error('pb in enerk computation'); 
end
% test of splitting:
mdl=fevisco('matsplit 100',mdl);
def=fe_eig(mdl,2); 
i0=find(def0.data>1); % no rb for compare
if norm(abs(def.def(:,i0))-abs(def0.def(:,i0)))/norm(abs(def0.def(:,i0)))>1e-6...
   || norm(def.data(i0)-def0.data(i0))/norm(def0.data(i0))>1e-6;
 error('Splitting material for enerk computations')
end
EnerK=fe_stress('enerk groupall',mdl,def);
%figure;plot(def.data,(def.data*2*pi).^2,def.data,2*sum(EnerK.data,1))
if norm((def.data*2*pi).^2-2*sum(EnerK.data,1)')>1e-5; 
  error('pb in enerk computation'); 
end

%% Nastran orthotropic checking
cd(t_nast('wd','bulks_divers'))
b=nasread('mat9_mode.dat');b=fe_mknl(b);
a=nasread('mat9_mode.op2');
a.K={};[a.K{1},a.K{2}]=upcom(a,'assemble mk');a.Klab={'m','k'};
i1=fe_c(a.DOF,b.DOF,'ind');a.K{1}=a.K{1}(i1,i1);a.K{2}=a.K{2}(i1,i1);
[svd(full(a.K{2})) svd(full(b.K{2}))]


%% --------------------------------------------------------------------
else; error('%s unknown',CAM);
end

