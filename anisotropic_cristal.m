function [out,out1,out2,out3]=anisotropic_cristal(varargin)

% Illustration of static computations for two anisotropic cristals
% anisotropic_cristal('auto') % Nominal case
%
% Example contributed by Nicolas Ranc and Etienne Balmes
%
% See <a href="matlab: sdtweb _taglist anisotropic_cristal">TagList</a>


%       Etienne Balmes
%       Copyright (c) 2001-2025 by INRIA and SDTools, All rights reserved.
%       Use under OpenFEM trademark.html license and LGPL.txt library license

if nargin==0; [CAM,Cam]=comstr('auto',1);carg=1;
else;CAM=varargin{1}; carg=2; [CAM,Cam]=comstr(CAM,1);
end

%% #Auto SCRIPT part of the example - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'auto')
RO.Quad=1;

if RO.Quad;model=anisotropic_cristal('mesh -div12 -quad');  % Generate mesh
else;model=anisotropic_cristal('mesh  -div12');  % Generate mesh
end
model=anisotropic_cristal('matnominal',model); % Set anisotropic prop.

% Display, compute and post-process
cf=feplot(model);%fecom('colordatamat')
def=fe_simul('static',cf); %fecom colordataevalz
cf.def=def;
%fecom('colordataevaly') % color displacement
fecom('colordatastress AtNode -comp 3 -edgealpha.01')  % color sigma zz

 RC=struct('Origin',[0 0 0;0 0 .02],'steps',linspace(0,1,100));
 cut=fe_caseg('StressCut',RC);cut.type='observe';
 cut=fe_caseg('stresscut -radius .1 -SelOut',cut,model); % get output info
 cut.StressObs.CritFcn='';
 R1=fe_caseg('stressobserve',cut,def);
 R1.name='Stress';R1.DimPos=[2 1 3];R1.X{2}=R1.X{2}(:,7)*1000;
 R1.Xlab{2}='Pos [mm]';
 ci=iiplot;ci.Stack{'curve',R1.name}=R1;iicom(ci,'iixonly',R1.name)
 
 
% display the center pile of elements
st='innode {r<.6e-3 & z>.005 & z<.015}';
cf.sel={st,'ColorData stress -comp3 atcenter'}; % sigma zz


[r1,EC]=anisotropic_cristal('stressatgauss',model,def,st);
figure(1);plot(r1.X{2}(:,1),r1.Y(3,:))

r2=fe_stress('SchmidCubicCenteredFace',r1);
[r3,pos]=anisotropic_cristal('GroupAsLines',[],r2,EC);

figure(1); 
RO.comp=1;
plot(squeeze(pos.Y(3,:,:)),squeeze(r3.Y(RO.comp,:,:)),'+')
xlabel('z');ylabel(r3.X{1}{RO.comp});

% Central line is point 5 for 27 integration point
if RO.Quad;RO.line=5;else;RO.line=3;end
figure(11);plot(squeeze(pos.Y(3,:,RO.line)),squeeze(r3.Y(:,:,RO.line)),'+')
xlabel('z');ylabel('Schmid shear'); legend(r3.X{1});


%% #StressAtGauss curve=anisotropic_cristal('stressatgauss',model,def,sel)
elseif comstr(Cam,'stressatgauss');[CAM,Cam]=comstr(CAM,4);

mo1=varargin{carg};carg=carg+1;
def=varargin{carg};carg=carg+1;
if carg<=nargin; st=varargin{carg};carg=carg+1;
 mo1.Elt=feutil(['selelt innode ' st],mo1);
end  
% post-process all stress components
[eltid,mo1.Elt]=feutil('eltidfix',mo1);
mo1.DOF=feutil('getdof',mo1);
C1=fe_stress('stress -gstate',mo1,def);C1.jGroup=1;
r2=[];
for jGroup=1:size(C1.GroupInfo,1)
 C1.jGroup=jGroup;
 [opt,match]=elem0('GaussObserve',C1.GroupInfo{C1.jGroup,[8 3 4]},mo1,C1);
 r1=C1.GroupInfo{C1.jGroup,5};
 r1.Y=reshape(r1.Y,size(r1.Y,1),[]);
 r1.X={r1.X{1} opt.EltId};r1.Xlab(3)=[];r1.Xlab{2}='Gauss.Elt';
 if isempty(r2);r2=r1;
 else; r2.Y=[r2.Y;r1.Y];r2.X{2}=vertcat(r2.X{2},r1.X{2});
 end
end
out=r2;
out1=C1.GroupInfo{C1.jGroup,8};

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% reformat computed stresses along lines assuming a series of
% elements. From 'comp','Gauss','EltId' to 'Comp','pos','Line'
elseif comstr(Cam,'groupaslines')

match=varargin{carg};carg=carg+1;
r1=varargin{carg};carg=carg+1;
EC=varargin{carg};carg=carg+1; %#ok<*NASGU,*ASGLU>

[CAM,Cam,RO.dir]=comstr('-dir',[-25 1 1],CAM,Cam);
if isempty(RO.dir);RO.dir=3;end  % Direction of extrusion

[r2,i1,i2]=unique(EC.w(:,setdiff(1:3,RO.dir)),'rows');
[i2,i3]=sort(i2);
if ~isfield(match,'Node') % case with Gauss.EltId
 i4=size(r1.Y); i4=[i4(1) size(EC.w,1) i4(2)/size(EC.w,1) i4(3:end)];
 r1.Y=reshape(r1.Y,i4);
 match.Node=r1.X{2}(:,2:end);
end
     
if length(size(r1.Y))==3;r1.Y=r1.Y(:,i3,:);
else;r1.Y=r1.Y(:,i3,:,:);
end
i1=size(r1.Y);
out=r1; out.Y=reshape(r1.Y,[size(r1.Y,1)   ...
     size(EC.w,1)/length(r2) length(r2) i1(3:end)]);
out.Y=reshape(permute(out.Y,[1 2 4 3]),size(out.Y,1),[],size(out.Y,3));
out.Xlab(2:3)={'pos','line'};

% reorder positions
r3=reshape(match.Node,size(EC.w,1),[],3);r3=r3(i3,:,:);
r3=reshape(r3,[size(EC.w,1)/length(r2) length(r2) size(r3,2) 3]);
r3=permute(r3,[4 1 3 2]);

out1=struct('X',{{{'x';'y';'z'}, ...
      (1:size(out.Y,2))' (1:size(out.Y,3))' }}, ...
    'Xlab',{{'pos','index','line'}}, ...
    'Y',reshape(r3,3,[],length(r2)));

%% #Mat handle material properties for the example ---------------------------
elseif comstr(Cam,'mat');[CAM,Cam]=comstr(CAM,4);
model=varargin{carg};carg=carg+1;

%% Set the appropriate material properties
if comstr(Cam,'nominal')
   
 % The behavior of a material grain is assumed orthotropic
 C11=168.4e9; C12=121.4e9; C44=75.4e9; % GPa
 C=[C11 C12 C12 0 0 0;C12 C11 C12 0 0 0;C12 C12 C11 0 0 0;
   0 0 0 C44 0 0;    0 0 0 0 C44 0;    0 0 0 0 0 C44]; 

 model.pl=[m_elastic('formulaPlAniso 1',C,basis('bunge',[5.175 1.3071 4.2012]));
      m_elastic('formulaPlAniso 2',C,basis('bunge',[2.9208 1.7377 1.3921]))];
 model=feutil('setmat 1 rho 1',model);
 model=feutil('setmat 2 rho 1',model);
 model=p_solid('default;',model);
 
% Internal verification
%sp_util('diag',11);fe_mknl('init',model);
else; error('Mat%s unknown',CAM);
end
out=model;

%% #Mesh generate sample mesh with appropriate boundary conditions and loads
elseif comstr(Cam,'mesh')

[CAM,Cam,RO.quad]=comstr('-quad',[-25 3],CAM,Cam);
[CAM,Cam,RO.Div]=comstr('-div',[-25 2 1],CAM,Cam);
 
n1=[[0 0 0;5 0 0;5 5 0;0 5 0]; % bottom edge
      [0 0 0;5 0 5;5 5 5;0 5 0]+ones(4,1)*[0 0 7.5] % middle cut at 45deg
      [0 0 0;5 0 0;5 5 0;0 5 0]+ones(4,1)*[0 0 20]; % top
      ]; 
n1(:,1)=n1(:,1)-max(n1(:,1))/2;n1(:,2)=n1(:,2)-max(n1(:,2))/2;
n1=n1/1000; % Use SI (m)

model=struct('Node',[(1:size(n1,1))'*[1 0 0 0] n1],'Elt',[]);
model=feutil('addelt',model,'hexa8',[1:8 1 1;9 12 11 10 5 8 7 6 2 2]);
model.Elt=feutil('orient',model);

if isequal(RO.Div,1) % coarse refine
 dx=linspace(0,1,4);
 dz=[linspace(0,.89,3) .9:.025:1]; % vertical refinement
else % Fine division
 dx=linspace(0,1,RO.Div);
 dz=[linspace(0,.9,9) .91:.01:1]; % vertical refinement
end
model=feutil('divide',model,dx,dx,dz);

% Define boundary conditions and loads
data=struct('sel','z==.02', ...
'eltsel','selface & innode {z==.02}','def',-1,'DOF',.19);
model=fe_case(model,'fixdof','Clamp',1, ...%'z==0&r<=.001 -DOF 1 2 3', ...
     'fixdof','Bottom','z==0 -DOF 3', ...
     'Fsurf','Traction',data);
     %     'fixdof','ClampRot','x==0 & y==5e-3 & z==0 -DOF 1', ...
     %'fixdof','EdgeX','x==0 -DOF 1','fixdof','EdgeY','y==0 -DOF 2', ...
 
% Define material properties, this needs to be cleaned up
model.pl=m_elastic('dbval1 Steel','dbval2 Steel');

if RO.quad;model=feutil('lin2quad',model);end
out=model;

%% #DevEb development script -------------------------------------------------
% see also sdtweb t_stress('aniso')
elseif comstr(Cam,'deveb')

%cd D:\sdtdata\bug\ranc
model=anisotropic_cristal('mesh -div1'); 
model=anisotropic_cristal('matnominal',model); % Set anisotropic prop.

% Display, compute and post-process
cf=feplot(model);%fecom('colordatamat')
def=fe_simul('static',cf); %fecom colordataevalz
cf.def=def;

[r1,EC]=anisotropic_cristal('stressatgauss',model,def);
r2=anisotropic_cristal('Schmid',r1);
[r3,pos]=anisotropic_cristal('GroupAsLines',[],r2,EC);

RO.comp=3;%zz
amp=r1.Y(RO.comp,:)';
map=struct('vertex',reshape(pos.Y,3,[])', ...
    'normal',amp*full(sparse(1,RO.comp,1,1,3)), ...
    'opt',2, ... % MAP at point
    'color',amp);
fecom('showmap',map)

figure(1); 
RO.comp=3;
plot(squeeze(pos.Y(3,:,:)),squeeze(r3.Y(1,:,:)),'+')
xlabel('z');ylabel(r3.X{1}{RO.comp});


%% -------------------------------------------------------------------
else; error('%s unknown',CAM);
end