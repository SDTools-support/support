function [out,out1,out2,out3]=p_pml(varargin)

%P_PML support functions for PML implementation in SDT
%
%       [ProId type form pow a0 x0 Lx0 y0 Ly0 z0 Lz0]
%          remi(form,[],1) : 0 rectan, 1 cyl, 2 spherical, 
%                            4 euclidian distance to planes 
%          remi(form,[],2) : 0 stress at node, 1 stress at gauss
%          remi(form,[],3) : selection of implementation
%         Type fe_mat('p_pml','SI',1)
%         Dir place holder for variations of direction formulation
%         Attenuation of the form a0*((x-x0)/Lx0)^pow
% 
%       Split command generates 'mpc','PmlInt' 
% 
%       See sdtweb      fem (handling materials section), pl, fe_mat, p_shell
%       See also help   fe_mat


%       Etienne Balmes
%       Copyright (c) 2001-2021 by SDTools, All Rights Reserved.
%       For revision information use p_pml('cvs')


%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>

if nargin<1; help p_solid;return; end
if ischar(varargin{1}); [CAM,Cam]=comstr(varargin{1},1); il=[]; carg=2;
else;il=varargin{1};[CAM,Cam]=comstr(varargin{2},1); carg=3;
end

%% #Default ------------------------------------------------------------------
if comstr(Cam,'default') % 
 
 model=[];if carg<=nargin;model=varargin{carg};carg=carg+1;end
 if isempty(model); out=p_solid('database');out=out(1);
 else;              out=fe_mat(['defaultil' CAM(8:end)],model); % sdtweb fe_mat('default')
 end
 
%% #Info ---------------------------------------------------------------------
elseif comstr(Cam,'info')

 r1=p_solid('database');fprintf('%s\n',r1.name);

%% #Dbval --------------------------------------------------------------------
elseif comstr(Cam,'dbval')

 while 1==1
  i1=strfind(comstr(Cam,-27),'-unit'); out1={};
  if ~isempty(i1)
   [Unit,i2,i3,i4]=sscanf(CAM(i1+5:end),'%s',1);
   i4=i1+(0:4+i4);CAM(i4)=''; [CAM,Cam]=comstr(CAM,1);
  else;Unit='';
  end
  i2=strfind(comstr(Cam,-27),'-punit');
  if ~isempty(i2)
   [PUnit,i3,i4,i5]=sscanf(CAM(i2+6:end),'%s',1);
   i5=i2+[0:5+i5];CAM(i5)=''; [CAM,Cam]=comstr(CAM,1);
  else;PUnit='';
  end
  
  if ischar(CAM); [i1,CAM,Cam]=comstr(CAM,'dbval','%i');else; i1=[];end
  if isempty(CAM)&&carg<=nargin; st=varargin{carg};carg=carg+1;
  else;st=CAM;end
  if isempty(st);
  elseif ischar(st); mat=p_solid('database',st);
  elseif isnumeric(st)
   [typ,st1,i4]=fe_mat('typep',st(2));
   mat=struct('il',st,'name',sprintf('%i',st(1)),'type',typ,'unit',st1);
  end
  if ~isempty(PUnit)
   r1=fe_mat(sprintf('convert %s %s',mat.unit,PUnit),mat.il(1:2));
   mat.il(2)=r1(2); mat.unit=PUnit;
  end
  if ~isempty(Unit)
   mat.il=fe_mat(sprintf('convert %s %s',mat.unit,Unit),mat.il);mat.unit=Unit;
  end
  r1=mat.il; if length(i1)==1; r1(1)=i1;end
  if ~isempty(il); i2=find(il(:,1)==r1(1)); else;i2=[];end
  if isempty(i2); i2=size(il,1)+1;end
  il(i2,1:length(r1))=r1; %#ok<AGROW>
  if carg>nargin; break;end
  CAM=varargin{carg};carg=carg+1;if ischar(CAM);[CAM,Cam]=comstr(CAM,1);end
 end
 out=il;


%% #DataBase -----------------------------------------------------------------
elseif comstr(Cam,'database') 

  st=comstr(CAM,9);
  if isempty(st)&&carg<=nargin; st=varargin{carg}; carg=carg+1;end
  [MatId,i2,i3,i4]=sscanf(st,'%i',1); if i2~=1; MatId=1;end
  st=comstr(st,i4);
  
  %sdtw('_ewt','P_PML Not implemented yet')

  % MatId Typ                     CordM Integ Stres Isop Field
  out.il=[MatId fe_mat('p_solid','SI',1) 0 -3 ]; 
  out.name='Default for topology';
  out.type='p_solid';
  out.unit='SI';
  out1='0Thick';


%% #PropertyUnitType  --------------------------------------------------------
elseif comstr(Cam,'propertyunittype')

 if nargin==1;out=1; return; end % return subtypes ID
 i1=varargin{carg};
 out1={};
 switch i1 % PropertySubType
 case {1,2} % [ProId type Coordm In Stress Isop Fctn  ]
   st={ ...
   'ProId'   0  'sdtweb(''p_pml'')';
   'Type'    0  '';
   'Form'     0  'Integration/formulation';
   'pow'     0  'Attenuation power'
   'a0'      0  'Value at edge'
   'x0'      4  'edge position'
   'Lx0'     4  'PML thickness'
   'y0'      4  'edge position'
   'Ly0'     4  'PML thickness'
   'z0'      4  'edge position'
   'Lz0'     4  'PML thickness'
   };
 otherwise; st={'ProId' 0 'sdtweb(''p_pml'')'; 'Type', 0, ''};
 end
 if ~isempty(strfind(Cam,'cell')); out=st; else; out=[st{:,2}]; end

%% #SubTypeString ------------------------------------------------------------
elseif comstr(Cam,'subtypestring')

 i1=varargin{carg}; carg=carg+1;
 switch i1
 case 1;  out='p_solid';
 otherwise; out='p_solid';
 end

%% #const : defines EltConst based on model content --------------------------
% [EltConst,NDNDim]=p_solid('Const',ElemF,integ,constit,model,Case,cEGI,RunOpt);
elseif comstr(Cam,'const')

if nargin==1
  %% #Fomulas_for_Topology for 3 blocks of rho and 3 blocs of dd for 5 strains  -3
  ConstitLab=ConstitLab1;
  EC=struct('Nw',8,'NDNLabels',{{'',',x',',y',',z'}}, ...
      'N',zeros(1,8),'Nnode',8);
  EC=p_pml('const',EC,[],[-1 0 0],[],[],[]);
  
  r1=diag([14*ones(9,1);zeros(15,1)]);
  i2=[1 2 3 0 0 0;2 4 5 0 0 0;3 5 6 0 0 0;0 0 0 7 0 0;0 0 0 0 8 0;0 0 0 0 0 9];
  i2(i2~=0)=i2(i2~=0)+14;
  i3=[1:3 5 6];i4=10:14;r1(i4,i4)=i2(i3,i3);  
  i3=[1:3 4 6];i4=15:19;r1(i4,i4)=i2(i3,i3);  
  i3=[1:3 4 5];i4=20:24;r1(i4,i4)=i2(i3,i3);  
  st=cell(size(r1));st(:)={' '};st(r1~=0)=ConstitLab(r1(r1~=0));
  disp('Constit mass');disp(st);
  disp(comstr(int32(r1),-30,struct('NoClip',1)))
  %% Topology for damping 
  % off diagonal blocks of -1, dax day daz on diagonal, rax,ray,raz
  % start with (3 * 6)*3 strains and 3 * 6 stresses
  JJ=[4+(0:6:12) 18+5+(0:6:12)  36+6+(0:6:12) 54+[4 6+5 12+6]]; % strain columns that will be removed
  r1=zeros(18*4);II=eye(6)*54;
  r1(1:6,3*18   +(1:6))=II;% E._xx( = (Bx * qx)')' (bx=Bx'*Nx) sx % Compute of 
  r1(1:6,3*18+ 6+(1:6))=II;% E._xx '* sy
  r1(1:6,3*18+12+(1:6))=II;% E._xx '* sz
  r1(18+6+(1:6),3*18   +(1:6))=II;% e_yy '* sx
  r1(18+6+(1:6),3*18 +6+(1:6))=II;% e_yy '* sy
  r1(18+6+(1:6),3*18+12+(1:6))=II;% e_yy '* sz
  r1(36+12+(1:6),3*18   +(1:6))=II;% e_zz '* sx
  r1(36+12+(1:6),3*18+ 6+(1:6))=II;% e_zz '* sy
  r1(36+12+(1:6),3*18+12+(1:6))=II;% e_zz '* sz
  
  % Block with b_x'
  II=eye(6);
  r1(54+(1:6),  1:6)=II;% sx' * e_xx
  r1(54+(1:6), 7:12)=II;% sx' * e_xy
  r1(54+(1:6),13:18)=II;% sx' * e_xz
  r1(60+(1:6),19:24)=II;% sy' * e_yx
  r1(60+(1:6),25:30)=II;% sy' * e_yy
  r1(60+(1:6),31:36)=II;% sy' * e_yz
  r1(66+(1:6),37:42)=II;% sz' * e_zx
  r1(66+(1:6),43:48)=II;% sz' * e_zy
  r1(66+(1:6),49:54)=II;% sz' * e_zz
  r1(JJ,:)=[];r1(:,JJ)=[];
  for j1=1:3  % sx relaxation
   i2=[1 2 3 0 0 0;2 4 5 0 0 0;3 5 6 0 0 0;0 0 0 7 0 0;0 0 0 0 8 0;0 0 0 0 0 9];
   i2(i2~=0)=i2(i2~=0)+23+(j1-1)*9;% dax,day,daz
   if j1==1; i3=[1:3 5 6];elseif j1==2;i3=[1:3 4 6];else;i3=[1:3 4 5];end 
   i4=(46:50)+(j1-1)*5; % sx ' Dax sx, sy' Day sy ... 
   r1(i4,i4)=i2(i3,i3);
  end
  for j1=1:3 % qx ' rhoax qx
   i2=eye(3)*(50+j1); i3=60+(j1-1)*3+(1:3);r1(i3,i3)=i2;
  end

  st=cell(size(r1));st(:)={' '};st(r1~=0)=ConstitLab(r1(r1~=0));
  disp('Constit Damping');
  disp(st);
  disp(['EC.ConstitTopology{3}=' comstr(int32(r1),-30,struct('NoClip',1)) ';'])
  %integrules('matrixrule',EC);
  
  figure(1);
  s=EC.StrainDefinition{2};
  subplot(121);ii_plp('spy',sparse(s(:,1),s(:,3),s(:,2)));title('Mass');
  s=EC.StrainDefinition{3};
  subplot(122);ii_plp('spy',sparse(s(:,1),s(:,3),s(:,2)));title('Damping');
 

  return
end
 EC=varargin{carg};carg=carg+1;
 integ=varargin{carg};carg=carg+1;
 if carg<=nargin; constit=varargin{carg};carg=carg+1;else;constit=[];end
 if carg<=nargin; model=varargin{carg};carg=carg+1;else;model=[];end
 if carg<=nargin; Case=varargin{carg};carg=carg+1;else;Case=[];end
 if carg<=nargin; cEGI=varargin{carg};carg=carg+1;else;cEGI=[];end
 if carg<=nargin; RunOpt=varargin{carg};carg=carg+1;else;RunOpt=[];end
 
 if integ(3,1)/integ(4,1)==3 &&size(constit,1)==24
 %% #PML.2 : frequency domain formulation with cut -2
  if ~ischar(EC) % Allow for integrule switching here
  else;EC=integrules(EC,-3);% Classical init (at nodes)
  end
  SE3=model;SE3.Elt=feutil(sprintf('selelt group%i',Case.jGroup),model);
  SE3.Node=feutil('getnode groupall',SE3);
  [pro,i1]=stack_get(model,'pro'); 
  for j1=1:size(pro,1) % search for amplitude table
    try; 
        if isfield(pro{j1,3},'il')&&pro{j1,3}.il(1)==integ(2);
            pro=pro{j1,3};SE3.Stack(i1(j1),:)=[];
            break;
        end
    end
  end
  if isempty(pro); 
   warning('Missing pro.a field giving fe/fi(f)');
  end
  
  SE3.il=[];
  SE3=feutil('setpro',SE3,[integ(2) fe_mat('p_solid','SI',1) 0 -3 0 100]);% il(6) for gradient
  cut=fe_caseg('StressCut-SelOut',struct('type','Gauss'),SE3);
  cut=fe_caseg('stresscutToStrain',cut); obs=cut.StressObs;
  
   %st=feval(p_pml('@ConstitLab2'));st{13}

  obs.pow = constit(4);
  obs.xiL = constit(6:7);
  obs.yiL = constit(8:9);
  obs.ziL = constit(10:11);
  obs.C=zeros(6);obs.C([1 2 3 7 8 9 13 14 15 22 29 36])= ...
       constit([16 17 18   18 19 20  18 20 21 22 23 24]);
  obs.rho=constit(15);
  obs.Vs = sqrt(obs.C(6,6)/obs.rho);%/2/pi;
   
  obs=feutilb('placeindof',model.DOF,obs);
  obs.trans=feutilb('placeindof',model.DOF,obs.trans);
  obs.a=pro.a;
  obs.toK=@stretchD; obs.MatType=3; % Ready for Z assembly (MatType=3 since depends on s)
  EC.StrainDefinition=cell(1,5);
  EC.toK=vhandle.matrix('ZMatrix',obs);
  EC.ConstitTopology=cell(1,5);
  EC.MatrixIntegrationRule=cell(1,5);
  EC.VectMap=int32(reshape(1:3*EC.Nnode,3,EC.Nnode)'); 
  out2=[];out1=3; out=EC;
  return;
 end
 if constit(3);rule=fix(constit(3)/1000);else; rule=-2; end
 if ~ischar(EC) % Allow for integrule switching here
 else;EC=integrules(EC,rule);% Classical init (at nodes)
 end
 rule=[1 EC.Nw];
 EC.DofLabels=DofLabT;
 out2=[]; 

 if remi(constit(3),[],3)==0
 %% #PML.1 0xx time implementation, see Hadrien Pinault -2
 % define the deformation vector: row, NDN, DDL, NwStart, NwRule
 r1=repmat([1 1 4 rule],24,1);r1(:,3)=4:27;r1(:,1)=(1:24)';
 EC.StrainDefinition{2}=r1; % MatType=2 mass
 EC.StrainLabels{2}= ... % see DofLabT
     {'qx_x','qy_x','qz_x','qx_y','qy_y','qz_y','qx_z','qy_z','qz_z', ... % strain 1:9 DOF 13:21
      'Sxx_x','Syy_x','Szz_x','Szx_x','Sxy_x', ... % strain 10:14 DOF 22:26
      'Sxx_y','Syy_y','Szz_y','Syz_y','Sxy_y', ... % strain 15:19 DOF 27:31
      'Sxx_z','Syy_z','Szz_z','Syz_z','Szx_z'};    % strain 20:24 DOF 32:36
 EC.ConstitTopology{2}=int32([14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 15 16 17 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 16 18 19 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 17 19 20 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 22 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 23 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 15 16 17 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 16 18 19 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 17 19 20 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 21 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 23 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 15 16 17 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16 18 19 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 17 19 20 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 21 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 22]);
% MatType3=viscous damping
 EC.StrainLabels{3}={'Exx_xx','Eyy_xx','Ezz_xx','Ezx_xx','Exy_xx', ...%1:5 Bx qx e_xx
     'Exx_xy','Eyy_xy','Ezz_xy','Ezx_xy','Exy_xy', ...%6:10 Bx qy  = exy
     'Exx_xz','Eyy_xz','Ezz_xz','Ezx_xz','Exy_xz', ...%11:15 Bx qz = exz
     'Exx_yx','Eyy_yx','Ezz_yx','Eyz_yx','Exy_yx', ... % 16:20 By qx = eyx
     'Exx_yy','Eyy_yy','Ezz_yy','Eyz_yy','Exy_yy', ... % 21:25 By qy = eyy
     'Exx_yz','Eyy_yz','Ezz_yz','Eyz_yz','Exy_yz', ... % 26:30 By qz = eyz
     'Exx_zx','Eyy_zx','Ezz_zx','Eyz_zx','Ezx_zx', ... % 31:35 Bz qx = ezx
     'Exx_zy','Eyy_zy','Ezz_zy','Eyz_zy','Ezx_zy', ... % 36:40 Bz qy = ezy
     'Exx_zz','Eyy_zz','Ezz_zz','Eyz_zz','Ezx_zz', ... % 41:45 Bz qz = ezz
      'Sxx_x','Syy_x','Szz_x','Szx_x','Sxy_x', ... %46:50
      'Sxx_y','Syy_y','Szz_y','Syz_y','Sxy_y', ... %51:55
      'Sxx_z','Syy_z','Szz_z','Syz_z','Szx_z', ...    % 56:60
      'qx_x','qy_x','qz_x','qx_y','qy_y','qz_y','qx_z','qy_z','qz_z'}';%61:69
EC.ConstitTopology{3}=int32([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 54 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 54 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 54 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 54 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 54 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 54 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 54 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 54 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 54 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0 54 0 0 0 0 0 0 0 0 0;
1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 24 25 26 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 25 27 28 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 26 28 29 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 31 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 33 34 35 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 34 36 37 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 35 37 38 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 39 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 41 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 42 43 44 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 43 45 46 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 44 46 47 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 48 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 49 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 51 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 51 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 51 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 52 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 52 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 52 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 53 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 53 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 53]);

EC.StrainDefinition{3}= [
     1 2 4 rule;  4 2  6 rule; 5 2  5 rule% bsx qx (dof 4:6)
     6 2 7 rule;  9 2  9 rule;10 2  8 rule% bsx qy  (dof 7:9)
    11 2 10 rule;14 2 12 rule;15 2 11 rule% bsx qz  (dof 10:12)
    17 3 5 rule;19 3 6 rule;20 3 4 rule% bsy qx
    22 3 8 rule;24 3 9 rule;25 3 7 rule% bsy qy
    27 3 11 rule;29 3 12 rule;30 3 10 rule% bsy qz
    33 4  6 rule;34 4  5 rule;35 4  4 rule% bsz qx
    38 4  9 rule;39 4  8 rule;40 4  7 rule% bsz qy
    43 4 12 rule;44 4 11 rule;45 4 10 rule% bsz qz
    46 1 13 rule; 47 1 14 rule;48 1 15 rule;49 1 16 rule;50 1 17 rule;% strain 16 /DOF 13 Ns sx 
    51 1 18 rule; 52 1 19 rule;53 1 20 rule;54 1 21 rule;55 1 22 rule;% Ns sy 
    56 1 23 rule; 57 1 24 rule;58 1 25 rule;59 1 26 rule;60 1 27 rule;% Ns sz 
    61 1 4 rule; 62 1 5 rule; 63 1 6 rule; % strain 31 qx = Nq qx (dof4)
    64 1 7 rule; 65 1 8 rule; 66 1 9 rule; % strain 34 qy = Nq qy (dof7)
    67 1 10 rule;68 1 11 rule;69 1 12 rule; % strain 37 qz = Nq qz (dof10)
    ];

 EC=integrules('matrixrule',EC); % \int f v (see elem0 for format)
 EC.RhsDefinition=[];
 
 EC.VectMap=int32(reshape(1:(3*4+5*3)*EC.Nnode,(3*4+5*3),EC.Nnode)'); 
 out1=3; % Tell that BuildNDN rule is 3D
 % Build the InfoAtNode interpolation to interpolate alpha_x
 if isempty(model)
 else
  NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
  [EGroup,nGroup]=getegroup(model.Elt);
  cEGI = EGroup(Case.jGroup)+1:EGroup(Case.jGroup+1)-1;
  NodePos=model.Elt(cEGI,1:EC.Nnode)';
  [i1,i2,i3]=unique(NodePos);i3=reshape(i3,size(NodePos));
  IA=struct('data',zeros(30,length(i1)),'NodePos',int32(i3),'lab',{{ ...
  'd11x','d12x','d13x','d22x','d23x','d33x','d44x','d55x','d66x', ...
  'd11y','d12y','d13y','d22y','d23y','d33y','d44y','d55y','d66y', ...
  'd11z','d12z','d13z','d22z','d23z','d33z','d44z','d55z','d66z', ...
  'rax','ray','raz'}});
  r2=model.Node(NNode(i1),5:7);
  if constit(7)==0; ax=zeros(size(r2,1),1)';% z=ConstitLab1;z(6)
  else; ax=constit(5)*(abs(r2(:,1)-constit(6))'/constit(7)).^constit(4); 
  end
  if constit(9)==0; ay=zeros(size(r2,1),1)';% z=ConstitLab1;z(6)
  else; ay=constit(5)*(abs(r2(:,2)-constit(8))'/constit(9)).^constit(4); 
  end
  if constit(11)==0; az=zeros(size(r2,1),1)';% z=ConstitLab1;z(8)
  else; az=constit(5)*(abs(r2(:,3)-constit(10))'/constit(11)).^constit(4); 
  end
  if remi(constit(3),[],1)==4 
   % #ConstDir=xx4 distance but assuming rect box 3
   ax=sqrt(sum([ax.^2;ay.^2;az.^2],1));
   %ax(ax==0)=norm(ax,'inf')/100;dbstack
   ay=ax; az=ax; 
   %cf=feplot;cf.def=struct('def',ax,'DOF',i1+.19); fecom colordata19
  end
  %d=double(EC.ConstitTopology{3});d(d~=0)=constit(d(d~=0));
  IA.data(1:9,:)=constit(24:32)*ax;%D ax
  IA.data(10:18,:)=constit(33:41)*ay; % Day
  IA.data(19:27,:)=constit(42:50)*az; % Daz
  IA.data(28:30,:)=[constit(51)*ax;constit(52)*ay;constit(53)*az]; % rax,ray,raz
  EC=feval(p_solid('@ctableGen'),EC,IA,ConstitLab1);
  out2=IA;
 end
 
 elseif remi(constit(3),[],2)==1 % x1x stress at gauss tests
 % #ConstDirX -2
 % define the deformation vector: row, NDN, DDL, NwStart, NwRule
 EC.StrainDefinition{1}=[1 2 1 rule;5 2 3 rule;6 2 2 rule]; %exx+ shear
 EC.StrainLabels{1}=p_solid('ConvStrain3D');
 EC.ConstitTopology{1}=int32(12*ones(6)); % always zero since no stiffness
 out1=3; % Tell that BuildNDN rule is 3D

 % Kinetic energy for positions
 EC.StrainDefinition{2}= [1 1 1 rule;2 1 2 rule; 3 1 3 rule];
 EC.StrainLabels{2}={'u','v','w'};
 EC.ConstitTopology{2}=int32(eye(3)*12);% No field for now
 % Viscous for positions xxx Need Position dependence xxx
 EC.StrainDefinition{3}= [1 1 1 rule;2 1 2 rule; 3 1 3 rule];
 EC.StrainLabels{3}={'u','v','w'};
 EC.ConstitTopology{3}=EC.ConstitTopology{2};
  
 EC=integrules('matrixrule',EC);
 % \int f v (see elem0 for format)
 EC.RhsDefinition=[];
 
 EC.VectMap=int32(reshape(1:3*EC.Nnode,3,EC.Nnode)'); 
 else; error('Not a valid case');
 end
 if 1==2
     
 end
 
 EC=feutil('rmfield',EC,'material'); % Generic does not use material
 
 % integ(1)==-1 is done at the very beginning of the command
 
 if nargout==0
  integrules('texstrain',EC); % Does not work here due to repetition
  try; EC=integrules('stressrule',EC);integrules('texstress',EC);end
 else; out=EC;
 end


% -------------------------------------------------------------------------
% #BuildConstit Implementation of elastic constitutive law building 3D
%[constit,integ,Inits,Data]=p_solid('buildconstit Nfield Nnode',ID,pl,il, ...
%    model,Case); This is called by the element 'integinfo' command
elseif comstr(Cam,'buildconstit');

  RunOpt=struct('warn',{{}},'Dim',0,'ProStack',[],'DoElmap',1);
  try;RunOpt.Dim=comstr(Cam(13:end),-1);end
  if isempty(RunOpt.Dim);RunOpt.Dim=[3 8];end
   
  ID=varargin{carg};carg=carg+1;out1=int32(ID);out3=struct;
  pl=varargin{carg};carg=carg+1;   
  if isempty(pl);mat=[];else;mat=pl(pl(:,1)==ID(1),:);end
  il=varargin{carg};carg=carg+1; 
  if isempty(il);pro=[];else;pro=il(il(:,1)==ID(2),:);end
  if isempty(mat); 
   RunOpt.warn{end+1}=sprintf('MatId %i is not matched',ID(1));
   st='m_null';unit=1;typ=1;
  else;  [st,unit,typ]=fe_mat('typem',mat(2));
  end
  RunOpt.MatT=fe_mat('typemstring',mat(2));
  RunOpt.ProT=fe_mat('typepstring',pro(1,2));
  
  %% #p_pml.1  - - - - - - - - - - - - - - - - - - - - - - - - - -
  if strcmp(RunOpt.ProT,'p_pml.1')&&strcmp(RunOpt.MatT,'m_elastic.1')
   ID(5)=pro(3);  % Integrule
   if length(mat)<5;mat(5)=0;end
   if length(pro)<12;pro(12)=0;end
   if remi(pro(1,3),[],2)==1
    %% old non element definition
    ID(3:4)=RunOpt.Dim(2)*[3 1]; % NDOF,NNode
    out3.ConstitLab={'-1','il(2)','form','m','a0','x0','Lmpl','0','0', ...%1:9
        'pl(1)','pl(2)','rho', ...% 10:12
        'D11','D21','D31','D41','D51','D61', ...
        'D12','D22','D32','D42','D52','D62', ...
        'D13','D23','D33','D43','D53','D63', ...
        'D14','D24','D34','D44','D54','D64', ...
        'D15','D25','D35','D45','D55','D65', ...
        'D16','D26','D36','D46','D56','D66'};% 13:48 is DD
    if length(mat)<6||mat(6)==0;mat(6)=mat(3)/2/(1+mat(4));end
    dd=m_elastic('formulaENG2DD',mat([3 4 6]));
    
    constit=[-1;pro(2:12)';mat([1 2 5])';dd(:)]; % p_solid will call p_pml
    out3.dd=dd;
   else; 
    %% New symmetric formulation with DOFs
    ID(3:4)=RunOpt.Dim(2)*[27 1]; % NDOF,NNode
    if length(mat)<6||mat(6)==0;mat(6)=mat(3)/2/(1+mat(4));end
    dd=m_elastic('formulaENG2DD',mat([3 4 6]));dd=pinv(dd);
    out3.ConstitLab=ConstitLab1;
    i2=[1 7 13 8 14 15 22 29 36]';%a=zeros(6);a(i2)=i2
    constit=[-1;pro(2:11)';mat([1 2 5])';dd(i2);dd(i2);dd(i2);dd(i2);mat([5;5;5])';1]; 
    out3.dd=dd;
   end
  %% #p_pml.2 : ULB frequency formulation  - - - - - - - - - - - - - - - - - -
  elseif strcmp(RunOpt.ProT,'p_pml.2')&&strcmp(RunOpt.MatT,'m_elastic.1')
    ID(3:4)=RunOpt.Dim(2)*[3 1]; % NDOF,NNode
    if mat(6)==0;mat(6)=mat(3)/2/(1+mat(4));end % G
    dd=m_elastic('formulaENG2DD',mat([3 4 6]));%dd=pinv(dd);
    out3.ConstitLab=ConstitLab2;
    i2=[1 7 13 8 14 15 22 29 36]';%a=zeros(6);a(i2)=i2
    pro(1,end+1:12)=0;
    constit=[-1;pro(2:12)';mat([1 2 5])';dd(i2);]; 
    out3.dd=dd;
  else
      error('%s Not implemented',RunOpt.ProT)
  end
  out=constit(:);
  out1=int32(ID(:));
  out2=elem0('elmapmat_og',[double(out1(4)) double(out1(3))/double(out1(4))]);
  if ~isempty(RunOpt.warn);sdtw('_nb','%s\n',RunOpt.warn{:});end

% -------------------------------------------------------------------------
elseif comstr(Cam,'builddof')
%% #BuildDof comming from p_solid

  RunOpt=varargin{carg};carg=carg+1;
  pl=varargin{carg};carg=carg+1;
  il=varargin{carg};carg=carg+1;
  model=varargin{carg};carg=carg+1;
  st=fe_mat('typepstring',il(1,2));
  
  if strcmpi(st,'p_pml.2')
   RunOpt.FieldDofs=[1 2 3]; % 3 DOF case all the time
  elseif remi(il(1,3),[],2)==0 % New formulation with DOFs (stress at node)
   RunOpt.FieldDofs=[1 2 3 13:21 22:21+15]; % 3 DOF sum, qx,qy,qz, sx,sy,sz
  else; % Old implement with stress at gauss added at the end
   RunOpt.FieldDofs=[1 2 3]; % 3 DOF case all the time
  end

  out=RunOpt;

%% #assemble : adaptation of assembly operations
elseif comstr(Cam,'assemble');[CAM,Cam]=comstr(CAM,9);
%% #assemble : adaptation of assembly operations

model=varargin{carg};carg=carg+1;
if comstr(Cam,'exit')
%% #assembleexit : fe_time description of stress and NL 
    
Case=varargin{carg};carg=carg+1;
Load=varargin{carg};carg=carg+1;
RunOpt=varargin{carg};carg=carg+1;

%model.K=cellfun(@(x)x.GetData,model.K,'uni',0);
%figure(1);clf;ii_plp('spy',SE)
iq=fe_c(Case.DOF,(22:36)'/100,'ind',2);
is=fe_c(Case.DOF,(22:36)'/100,'ind');

% see sdtweb p_pml('nl_inout')
Ms=model.K{1}(is,is);Ms=Ms+spalloc(size(Ms,1),size(Ms,2),0);
NL=struct('type','nl_inout', 'c',Ms\model.K{2}(is,iq),'b',model.K{2}(iq,is));
%NL.Ca=full(diag(Ms\model.K{2}(is,is)));
if nnz(Ms-diag(diag(Ms)))>0
 NL.Ca=Ms\model.K{2}(is,is); % figure(1);spy(model.K{2}(is,is))
else; NL.Ca=diag(Ms\model.K{2}(is,is));
end
NL.Fu=p_pml('@nl_inout');
NL.unl=zeros(size(NL.c,1),1,3);
NL.c=v_handle('mklst',sparse(NL.c));
NL.b=v_handle('mklst',sparse(NL.b'));
NL.StoreFNL=0;

model.K=cellfun(@(x)x(iq,iq),model.K,'uni',0);
Case.DOF=Case.DOF(iq);
Case.T=Case.T(:,iq);
Load.def=Load.def(iq,:);Load.DOF=Load.DOF(iq);
%figure(1);ii_plp('spy',SE)

 model.NL={'NL','pml',NL};
 
 out={model,Case,Load,[]}; clear model
 if RunOpt.NeedSens;
     out{4}=evalin('caller','Sens');
 end
 % reset upper level matrix instances
 if ~isfield(RunOpt,'isMVR')||RunOpt.isMVR==0 % (out{1} handled by vout)
  assignin('caller','K',out{1}.K); evalin('caller','model.K=K;'); 
 else
  % do mkl when assembly procedure is shunted by MVR presence
  if 1==2%RunOpt.DoMkl;  
   k1=out{1}.K; out{1}.K=[]; 
   i1=cellfun(@(x) ~isa(x,'v_handle')&&isnumeric(x)&&~issparse(x),k1);
   if any(i1); k1(i1)=cellfun(@sparse,k1(i1),'uni',0); end
   out{1}.K=v_handle('mkls','k1');
  end
 end
 
%%    
else
 %% #AssembleOld adaptation of assembly 
[SE,CE]=fe_case(model,'assemble -matdes 2 3 1 -SE NoT');

RO.ind=find(cellfun(@(x)isfield(x,'type')&&strcmpi(x.type,'p_pml'),CE.GroupInfo(:,end)))';
RO.DOFs=cell(length(RO.ind),3);

for jGroup=RO.ind
 indGroup=find(RO.ind==jGroup);
 pro=feutil(sprintf('Getil %i -struct1',CE.GroupInfo{jGroup,3}(2)),model);
 if pro.a0==0;error('Not a valid case');end
 EC=CE.GroupInfo{jGroup,8};
 mo1=model;
 mo1.Elt=feutil('selelt proid',mo1,pro.ProId);
 %mo1.il=p_solid(sprintf('dbval %i d3 -3',pro.ProId));mo1.il(6)=100;% il(6) for gradient
 mo1.il=p_solid(sprintf('dbval %i d3 -2',pro.ProId));mo1.il(6)=100;% il(6) for gradient
 
 cut=fe_caseg('StressCut-SelOut',struct('type','Gauss'),mo1);
 cut=fe_caseg('stresscutToStrain',cut); 
 obs=cut.StressObs; obs.iq=fe_c(SE.DOF,obs.DOF,'ind'); obs.nd=0;% PML DOF only
 RO.DOFs{indGroup,3}=obs.iq;
 RO.Nsig=6;
 if ~isempty(pro.Lx0)&&pro.Lx0; obs.ax=@(x)pro.a0*(abs(x-pro.x0)/pro.Lx0).^pro.pow;
     obs.nd=obs.nd+1;
 else; obs.ax=[];
 end
 if ~isempty(pro.Ly0)&&pro.Ly0; obs.ay=@(x)pro.a0*(abs(x-pro.y0)/pro.Ly0).^pro.pow;
     obs.nd=obs.nd+1;
 else; obs.ay=[];
 end
 if ~isempty(pro.Ly0)&&pro.Lz0; obs.az=@(x)pro.a0*(abs(x-pro.z0)/pro.Lz0).^pro.pow;
     obs.nd=obs.nd+1;
 else; obs.az=[];
 end
 EC.pro=pro;EC.obs=obs;
  dd=reshape(CE.GroupInfo{jGroup,4}(16:51),6,6);%dd(:,2:4)=0;dd(2:4,:)=0;
  [II,JJ,KK]=find(dd); % scale inversion to avoid warning  
  i1=(1:size(obs.Elt,1)-1)*size(dd,1);
  i1=repmat(i1,size(dd,1),1)+repmat((-size(dd,1)+1:0)',1,size(i1,2));
  i2=(1:size(obs.Elt,1)-1)*size(dd,2);
  i2=repmat(i2,size(dd,2),1)+repmat((-size(dd,2)+1:0)',1,size(i2,2));
  
  II=i1(II,:);JJ=i2(JJ,:);KK=repmat(KK,1,size(II,2));
  EC.ddg=sparse(II(:),JJ(:),KK(:),size(dd,1)*size(obs.Node,1),size(dd,1)*size(obs.Node,1));
  
  st1={'u11','u13','u12','u22','u23','u21','u33','u32','u31'};
  [ia,ib]=ismember(st1,obs.X{1});ia=[1 5 6 2 4 6 3 4 5]; 
  % sdtweb sdt.pdf#subsection.6.1.1
  i3=[1 6 9];i1=ia(i3);i2=ib(i3);%st1(i3)
  Bx=sparse(repmat(i1,size(obs.Node,1),1)+(0:size(obs.Node,1)-1)'*ones(size(i1))*size(dd,2), ...
      repmat(i2,size(obs.Node,1),1)+(0:size(obs.Node,1)-1)'*ones(size(i2))*9, ...
      1,size(obs.Node,1)*size(dd,1),size(obs.Node,1)*9);
  i3=[3 4 8];i1=ia(i3);i2=ib(i3);%st1(i3)
  By=sparse(repmat(i1,size(obs.Node,1),1)+(0:size(obs.Node,1)-1)'*ones(size(i1))*size(dd,2), ...
      repmat(i2,size(obs.Node,1),1)+(0:size(obs.Node,1)-1)'*ones(size(i2))*9, ...
      1,size(obs.Node,1)*size(dd,1),size(obs.Node,1)*9);  
  i3=[2 5 7];i1=ia(i3);i2=ib(i3);%st1(i3)
  Bz=sparse(repmat(i1,size(obs.Node,1),1)+(0:size(obs.Node,1)-1)'*ones(size(i1))*size(dd,2), ...
      repmat(i2,size(obs.Node,1),1)+(0:size(obs.Node,1)-1)'*ones(size(i2))*9, ...
      1,size(obs.Node,1)*size(dd,1),size(obs.Node,1)*9);
  i1=ia;i2=ib;%i3=[2 5 7];i1=ia(i3);i2=ib(i3);%st1(i3)
  Ba=sparse(repmat(i1,size(obs.Node,1),1)+(0:size(obs.Node,1)-1)'*ones(size(i1))*size(dd,2), ...
      repmat(i2,size(obs.Node,1),1)+(0:size(obs.Node,1)-1)'*ones(size(i2))*9, ...
      1,size(obs.Node,1)*size(dd,1),size(obs.Node,1)*9);
  %[obs.X{1}';num2cell(full(Bx(1:6,1:9)+By(1:6,1:9)+Bz(1:6,1:9)))]
  EC.wj=diag(sparse(reshape(repmat(obs.wjdet',size(dd,1),1),[],1)));
  %k1=full(SE.K{3}(1:48,1:48));k2=full([obs.cta'*Ba'*wj*ddg*Ba*obs.cta]);norm(k1-k2,'inf')
  EC.Bx=Bx*obs.cta;EC.By=By*obs.cta;EC.Bz=Bz*obs.cta;  % PML strain observation
 
 %% Now extend DOF for each PML group q,qx,qy,qz,sx
   % qx 13:15, qy=16:18, qz=19:21,sx_y_z=22:24
   RO.DOFs{indGroup,1}=[obs.DOF+.12;obs.DOF+.15;obs.DOF+.18];% PML trans DOFs
   % PML stress DOF 3DOF per direction/per gauss point
   RO.DOFs{indGroup,2}=obs.nd*size(obs.Node,1)*RO.Nsig; 
 %% Save element constants
 CE.GroupInfo{jGroup,8}=EC;
end

%% Now deal with DOF augmentation  : add local translations, then stresses
%  starting at .21 using arbitrary nodes
DOFs=unique(fix(SE.DOF(vertcat(RO.DOFs{:,3}))));
i2=sum(vertcat(RO.DOFs{:,2}));
DOFs=feutil('getdof',(21+(1:ceil(i2/length(DOFs))))'/100,DOFs);
DOFs(i2+1:end)=[];i3=unique(round(vertcat(RO.DOFs{:,1})*100))/100;
i1=length(i3)+size(SE.DOF,1);
for j1=1:size(RO.DOFs,1); i1=i1(end)+(1:RO.DOFs{j1,2})';RO.DOFs{j1,2}=i1;end
DOFs=[i3;DOFs];
SE.DOF=[SE.DOF;DOFs];CE.DOF=[CE.DOF;DOFs];
CE.T(size(SE.DOF,1),size(CE.DOF,1))=0; 
for j1=1:length(SE.K);SE.K{j1}(size(CE.T,1),size(CE.T,1))=0;end

for jGroup=RO.ind
    % From Ni,j reconstruct Bx block diagonal
    %reshape(CE.GroupInfo{1,4}(3:end),6,6)-dd
  EC=CE.GroupInfo{jGroup,8};obs=EC.obs;
  csx=EC.ddg*EC.Bx; csy=EC.ddg*EC.By; csz=EC.ddg*EC.Bz;
  bsx=EC.Bx'*(EC.wj); bsy=EC.By'*EC.wj;bsz=EC.Bz'*EC.wj; % Low freq limit : -Kqs*q+sigma=0
  %k1=full(SE.K{3}(1:48,1:48));k2=full([(bsx+bsy+bsz)*(csx+csy+csz)]);norm(k1-k2,'inf')

  %ix=any(EC.Bx,2);iy=any(EC.By,2);iz=any(EC.Bz,2);
  ix=1:size(csx,1);iy=1:size(csy,1);iz=1:size(csz,1);%for RO.Nsig=6;
  csx=csx(ix,:);bsx=bsx(:,ix);csy=csy(iy,:);bsy=bsy(:,iy);csz=csz(iz,:);bsz=bsz(:,iz);
 
  if ~isequal(obs.trans.DOF,SE.DOF(obs.iq));error('problem');end
  rhowj=diag(sparse(reshape(repmat(obs.wjdet'*CE.GroupInfo{jGroup,4}(15),3,1),[],1)));
  M=obs.trans.cta'*(rhowj)*obs.trans.cta; csa=[];
  if ~isempty(obs.ax)
   Ax=diag(sparse(reshape(repmat(obs.ax(obs.Node(:,5)'),3,1),[],1)));
   Cax=obs.trans.cta'*(rhowj*Ax)*obs.trans.cta;%norm(SE.K{2}-C,'inf')
   Ax=diag(sparse(reshape(repmat(obs.ax(obs.Node(:,5)'),RO.Nsig,1),[],1)));
   csa=[csa;csx];
  else; Ax=[];Cax=[];
  end
  if ~isempty(obs.ay)
   Ay=diag(sparse(reshape(repmat(obs.ay(obs.Node(:,6)'),3,1),[],1)));
   Cay=obs.trans.cta'*(rhowj*Ay)*obs.trans.cta;%norm(SE.K{2}-C,'inf')
   Ay=diag(sparse(reshape(repmat(obs.ay(obs.Node(:,5)'),RO.Nsig,1),[],1)));
   csa=[csa;csy];
  else; Ay=[];Cay=[];
  end
  if ~isempty(obs.az)
   Az=diag(sparse(reshape(repmat(obs.az(obs.Node(:,7)'),3,1),[],1)));
   Caz=obs.trans.cta'*(rhowj*Az)*obs.trans.cta;%norm(SE.K{2}-C,'inf')
   Az=diag(sparse(reshape(repmat(obs.az(obs.Node(:,5)'),RO.Nsig,1),[],1)));
   csa=[csa;csz];
  else; Az=[];Caz=[];
  end
  is=RO.DOFs{find(RO.ind==jGroup),2};
  Css=speye(size(csx,1)+size(csy,1)+size(csz,1)); 
  Zsq=spalloc(size(Css,1),size(csa,2),0);
  Zss=spalloc(size(Css,1),size(Css,1),0);
  Zqq=spalloc(size(csa,2),size(csa,2),0);
  % Now extend the DOF : in pml : q,qx,qy,qz,sx
   % qx 13:15, qy=16:18, qz=19:21,sx_y_z=22:24

  % q are not active DOF in PML but sum of qx,qy,qs
  i3=fe_c(SE.DOF,obs.DOF,'ind');ind=1:length(i3);
  RO.DOFs{find(RO.ind==jGroup),3}=fe_c(CE.DOF,obs.DOF,'ind'); % q DOF inside a PML (eliminated at the end)
  RO.ifx=fe_c(SE.DOF,obs.DOF+.12,'ind');i2=fe_c(CE.DOF,obs.DOF+.12,'ind');
  CE.T(:,i2)=sparse(RO.ifx,ind,1,size(CE.T,1),length(ind))+ ...
      sparse(i3,ind,1,size(CE.T,1),length(ind));
  RO.ify=fe_c(SE.DOF,obs.DOF+.15,'ind');i2=fe_c(CE.DOF,obs.DOF+.15,'ind');
  CE.T(:,i2)=sparse(RO.ify,ind,1,size(CE.T,1),length(ind))+ ...
      sparse(i3,ind,1,size(CE.T,1),length(ind));
  RO.ifz=fe_c(SE.DOF,obs.DOF+.18,'ind');i2=fe_c(CE.DOF,obs.DOF+.18,'ind');
  CE.T(:,i2)=sparse(RO.ifz,ind,1,size(CE.T,1),length(ind))+ ...
      sparse(i3,ind,1,size(CE.T,1),length(ind));
 %     sparse(i2+3*ind(end)+size(Css,1),1:size(Css,1),1,size(Tp,1),size(Css,1))  ]; % I q_s
  %figure(1);spy(CE.T)

  if remi(pro.Dir,[],3)==1 % Dir=1xx
   %% check that static goes to nominal modes with merged fields
   dbstack; keyboard;
  z_ab=spalloc(size(SE.K{1},1),size(Zqq,2),0);
  SE.K{1}= [SE.K{1} z_ab
            z_ab' M];
  SE.K{3}= [SE.K{3} z_ab
            z_ab' (bsx+bsy+bsz)*(csx+csy+csz)];
  SE.K{2}(size(SE.K{1},1),size(SE.K{1},1))=0;
  CE.T=CE.T(1:end-nnz(ix)-2*size(M,1),1:end-nnz(ix)-2*size(M,1));
  SE.DOF=SE.DOF(1:end-nnz(ix)-2*size(M,1));
  CE.DOF=CE.DOF(1:end-nnz(ix)-2*size(M,1));
  elseif remi(pro.Dir,[],3)==2 % Dir=2xx
  %% No relaxation by field split
   dbstack; keyboard;
  z_ab=spalloc(size(SE.K{1},1),size(Zqq,2)*3,0);
  SE.K{1}= [SE.K{1} z_ab
            z_ab' [M    Zqq  Zqq ; 
                   Zqq  M    Zqq ; ...
                   Zqq Zqq   M ]
            ];
  SE.K{2}= [SE.K{2} z_ab
            z_ab' [Cax*0  Zqq Zqq; 
                   Zqq  Zqq Zqq ; ...
                   Zqq Zqq Zqq]
            ];
  SE.K{3}= [SE.K{3} z_ab
            z_ab' [bsx*[csy+csz+csx csy+csz+csx csy+csz+csx]; 
                   bsy*[csy+csz+csx csy+csz+csx csy+csz+csx];
                   bsz*[csy+csz+csx csy+csz+csx csy+csz+csx]]
            ];
  CE.T=CE.T(1:end-nnz(ix)-0*size(M,1),1:end-nnz(ix)-0*size(M,1));
  SE.DOF=SE.DOF(1:end-nnz(ix)-0*size(M,1));
  CE.DOF=CE.DOF(1:end-nnz(ix)-0*size(M,1));
  elseif remi(pro.Dir,[],3)==0 % Dir=0xx
  %% Now field split and relaxation
  SE.K{1}=addMat(SE.K{1},M,RO.ifx,RO.ifx);
  SE.K{1}=addMat(SE.K{1},M,RO.ify,RO.ify);
  SE.K{1}=addMat(SE.K{1},M,RO.ifz,RO.ifz);
  
  % Relaxation fields of qx,qy,qz
  if ~isempty(Cax); SE.K{2}=addMat(SE.K{2},Cax,RO.ifx,RO.ifx);end
  if ~isempty(Cay); SE.K{2}=addMat(SE.K{2},Cay,RO.ify,RO.ify);end
  if ~isempty(Caz); SE.K{2}=addMat(SE.K{2},Caz,RO.ifz,RO.ifz);end
  % Observation of stresses associated with qx+qy+qz
  SE.K{2}=addMat(SE.K{2},-csa,is,RO.ifx);
  SE.K{2}=addMat(SE.K{2},-csa,is,RO.ify);
  SE.K{2}=addMat(SE.K{2},-csa,is,RO.ifz);
  
  % I \dot sigma : identity for stress derivative
  SE.K{2}=addMat(SE.K{2},speye(length(is)),is,is);
  % Now stiffness contributions associated for directions not PML stresses
  if isempty(Cax)
   SE.K{3}=addMat(SE.K{3},[bsx*[csx csx csx];bsy*[csx csx csx];bsz*[csx csx csx]], ...
       [RO.ifx;RO.ify;RO.ifz],[RO.ifx;RO.ify;RO.ifz]);
  end
  if isempty(Cay)
   SE.K{3}=addMat(SE.K{3},[bsx*[csy csy csy];bsy*[csy csy csy];bsz*[csy csy csy]], ...
       [RO.ifx;RO.ify;RO.ifz],[RO.ifx;RO.ify;RO.ifz]);
  end
  if isempty(Caz)
   SE.K{3}=addMat(SE.K{3},[bsx*[csz csz csz];bsy*[csz csz csz];bsz*[csz csz csz]], ...
       [RO.ifx;RO.ify;RO.ifz],[RO.ifx;RO.ify;RO.ifz]);
  end
  
  % Stress contribution of PML + relaxation
  if ~isempty(Cax)
      SE.K{3}=addMat(SE.K{3},[bsx;bsy;bsz;Ax], ...
         [RO.ifx;RO.ify;RO.ifz;is],is);      
  end
  if ~isempty(Cay)
      SE.K{3}=addMat(SE.K{3},[bsx;bsy;bsz;Ay], ...
         [RO.ifx;RO.ify;RO.ifz;is],is);      
  end
  if ~isempty(Caz)
      SE.K{3}=addMat(SE.K{3},[bsx;bsy;bsz;Az], ...
         [RO.ifx;RO.ify;RO.ifz;is],is);      
  end
  else; error('Not a valid case');
  end
end % Loop on PML group
% add unit observe for PML stresses
i1=vertcat(RO.DOFs{:,2});
CE.T(i1,end+(-length(i1)+1:0))=speye(length(i1));
% remove q DOF used in PML from CE.DOF since qx,qy,qz
i1=unique(vertcat(RO.DOFs{:,3})); CE.T(:,i1)=[];CE.DOF(i1)=[];

%% Further elimination
i1=[];
if remi(pro.Dir,[],2)==2
   %% Fix qx23, qy13, qz12
   i1=fe_c(CE.DOF,[14 15 16 18 19 20]'/100,'ind');
end
CE.DOF(i1)=[];CE.T(:,i1)=[];

if 1==2 % Display matrix topology
 i1={feutil('getnode group1 & notin {group2}',model) % C
   feutil('getnode group1 & group2',model)           % I
    feutil('getnode group2 & notin{group1}',model)}; % qP
 i1=cellfun(@(x)fe_c(SE.DOF,feutil('getdof',x(:,1),(1:3)'/100),'ind'),i1,'uni',0);
 i1{4}=setdiff((1:length(SE.DOF))',vertcat(i1{:})); % qS
 figure(1);
 for j1=1:3
  subplot(1,3,j1);ii_plp('spyunsymm',struct('K',SE.K{j1},'ind',{i1}));
  title(SE.Klab{j1});
 end
 %figure(1);ii_plp('spy',struct('K',SE.K{2},'ind',{i1}));
end
out=SE;out1=CE;
end

elseif comstr(Cam,'solve');[CAM,Cam]=comstr(CAM,6);
%% #solve : PML associated solves 
if comstr(Cam,'dfrf')
 %% #SolveDfrf : compute direct frequency response
 model=varargin{carg};carg=carg+1;
 if carg>nargin;RO=struct;else;RO=varargin{carg};carg=carg+1;end
 il=feutil('getil',model);
 if isfield(RO,'Range')&&(isfield(RO.Range,'ncx')|| ...
     (isfield(RO.Range,'lab')&&any(strcmpi(RO.Range.lab,'ncx'))))
  %% Periodic model actually solved using fe_homo dftDfrf
  % called in Dynavoie learning phase of RedPML
  % xxx should be a direct call to dftDfrf
  %[SE,CE]=p_pml('assemble -NoT',model);Load=fe_load(SE,CE,'NoT');
  %SE.K=feutilb('tkt',CE.T,SE.K);SE.DOF=CE.DOF;

  if isfield(model,'K')&&~isempty(model.K)&&isfield(model,'Case')&&isfield(model,'Load')
    SE=model;CE=model.Case;Load=model.Load;
  else;[SE,CE,Load]=fe_case('assemble -SE -matdes 2 3 1 4 -NoT -load',model);
  end
  if ~isfield(RO,'AssembleCall'); % Bypass assembly in fe_homo dftdfrf
      RO.AssembleCall='[SE,Case,Load]=deal(mdl{:});mdl=SE;';
  end
  data=fe_case(SE,'getdata','Symmetry');% Do not apply periodicity on stress
  if ~isempty(data)
   data.IntDof=fe_c(SE.DOF, ...
      feutil('getdof',data.IntNodes(:,1),[1:3 13:21 22 25 26 28 30 31 34 35 36]'/100),'dof');
   i1=sparse(data.IntNodes(:,1),1,data.IntNodes(:,2));
   data.IntDof(:,2)=rem(data.IntDof,1)+i1(fix(data.IntDof));
   SE=fe_case(SE,'cyclic','Symmetry',data);
   CE=stack_set(CE,'cyclic','Symmetry',data);
   SE.Case=stack_set(SE.Case,'cyclic','Symmetry',data);
   d1=fe_homo('dftdfrf',{SE,CE,Load},RO);
   d1.LabFcn=['sprintf(''%s (%g Hz, ncx%.0f)'',def.lab_in{def.data(ch,2)},' ...
      'def.data(ch,1),def.Range.val(def.data(ch,3),1))'];
  end
  if isfield(RO,'defdouble')&&RO.defdouble
   d1=fe_cyclic('defdouble',d1);
  end
  
 elseif any(strcmpi(fe_mat('typepstring',il(:,2)),'p_pml.1'))
 %% Single model frequency domain format
 
 elseif any(strcmpi(fe_mat('typepstring',il(:,2)),'p_pml.3'))
 %% Single model old format
  error('Obsolete');   
  [SE,CE]=p_pml('assemble -NoT',model);
  Load=fe_load(SE,CE,'NoT');
  if size(Load.def,1)==size(CE.T,1);Load.def=CE.T'*Load.def;end
  K=feutilb('tkt',CE.T,SE.K);
  freq=stack_get(model,'info','Freq','get');
  d1=struct('def',zeros(length(SE.DOF),length(freq)),'DOF',SE.DOF,'data',freq);
  for j1=1:length(freq)
   w=freq(j1)*2*pi;
   Z=feutilb('sumkcoef',K,[-w^2 1i*w 1]);
   %[L,U]=lu(Z);
   d1.def(:,j1)=CE.T*(Z\Load.def)*(1i*w).^2;
  end
 else; error('Not a valid case');
 end
out=d1;
out1=SE;out1.Case=CE;
elseif comstr(Cam,'mode')
 %% #SolveModes : compute complex modes from reduced state space
 SE=varargin{carg};carg=carg+1;
 d1=varargin{carg};carg=carg+1;
 
 RP=stack_get(SE,'info','PML','get');i1=RP.nodeSet;
 i1=cellfun(@(x)fe_c(SE.DOF,feutil('getdof',x(:,1),(1:3)'/100),'ind'),i1,'uni',0);
 i1{4}=setdiff((1:length(SE.DOF))',vertcat(i1{:})); % qS
 RO.iset=i1; RO.iset(:,2)={'C';'I';'qP';'qS'};

%getdef=@(x)real(x.def(:,1));
getdef=@(x)[real(x.def) imag(x.def)];
resc=@(x)x*diag(sparse(1./sum(abs(x))));
Tc=getdef(d1);Tc(vertcat(RO.iset{2:4}),:)=0;Ti=getdef(d1);Ti(vertcat(RO.iset{[1 3:4]}),:)=0;
Tci=[Tc Ti];Tci=resc(Tci);

if 1==2
 [Tc,fc,mdr]=fe_norm(Tc,SE.K{1},SE.K{3});getdef=@(x)[real(x.def) imag(x.def)]*mdr;
 Tc=getdef(d1);Tc(vertcat(RO.iset{2:4}),:)=0;Ti=getdef(d1);Ti(vertcat(RO.iset{[1 3:4]}),:)=0;
 Tci=[Tc Ti];Tci=resc(Tci);
end

Tqp=getdef(d1);Tqp(vertcat(RO.iset{[1 2 4]}),:)=0;Tqp=resc(Tqp);
Tqs=getdef(d1);Tqs(vertcat(RO.iset{[1 2 3]}),:)=0;Tqs=resc(Tqs);
Z=@(x,y)zeros(size(x,2),sum(cellfun(@(z)size(z,2),y)));
I=@(x,y)eye(size(x,2));
Ta=[I(Tc) Z(Tc,{Ti,Tc,Ti,Tqp,Tqp,Tqs}) % Keep q_q
    Z(Ti,{Tc}) I(Ti) Z(Ti,{Tc,Ti}) I(Tqp) Z(Ti,{Tqp,Tqs}) % Keep q_i = q_p
    Z(Ti,{Tc,Ti}) I(Tc) Z(Ti,{Ti,Tqp,Tqp,Tqs}) % Keep q_q dot  
    Z(Ti,{Tc,Ti,Tc}) I(Ti) Z(Ti,{Tqp}) I(Tqp) Z(Ti,{Tqs}) % Keep dot q_i = dot q_p
    Z(Tqs,{Tc,Ti,Tc,Ti,Tqp,Tqp}) I(Tqs) % Keep q_s
    ]';

Cci=Tci'*SE.K{2}*Tci; Mci=Tci'*SE.K{1}*Tci; Kci=Tci'*SE.K{3}*Tci; 

Cqp=Tqp'*SE.K{2}*Tqp;Mqp=Tqp'*SE.K{1}*Tqp;Kqp=Tqp'*SE.K{3}*Tqp;
Cqs=Tqs'*SE.K{2}*Tqs;Mqs=Tqs'*SE.K{1}*Tqs;Kqs=Tqs'*SE.K{3}*Tqs;
cs=Tqs'*SE.K{2}*Tqp; bs=Tqp'*SE.K{3}*Tqs; 
% eig(Kqs,-Cqs)/2/pi
B=-[Cci Mci Z(Cci,{Cqp,Cqp,Cqs})
   Mci Z(Mci,{Mci,Cqp,Cqp,Cqs}) 
   Z(Cqp,{Mci,Mci}) Cqp Mqp Z(Cqp,{Cqs})
   Z(Cqp,{Mci,Mci})  Mqp Z(Mqp,{Cqp,Cqs})
   Z(Cqs,{Mci,Mci}) cs Z(Cqs,{Mqp}) Cqs
   ];
A=[Kci Z(Cci,{Kci,Cqp,Cqp,Cqs})
   Z(Cci,{Kci}) -Mci Z(Cci,{Cqp,Cqp,Cqs})
   Z(Kqp,{Kci,Kci}) Kqp Z(Kqp,{Cqp}) bs
   Z(Kqp,{Kci,Kci,Kqp}) -Mqp Z(Mqp,{Cqs}) 
   Z(Kqs,{Kci,Kci,Kqp,Kqp}) Kqs
   ];
if 1==1
 coef=1e-8;
 Tap=Ta;Tap(end+(-size(Kqs,1)+1:0),:)=Tap(end+(-size(Kqs,1)+1:0),:)/coef;
 TaT=Ta';TaT(:,end+(-size(Kqs,1)+1:0))=TaT(:,end+(-size(Kqs,1)+1:0))*coef;
 %Ta(end+(-size(Kqs,1):0),:)=Ta(end+(-size(Kqs,1):0),:)*1e-6;
 Ar=TaT*A*Tap;Br=TaT*B*Tap;
 %[svd(Ar) svd(Br)]
end
[psi,lambda]=eig(Ar,Br);psi=Tap*psi;lambda=diag(lambda);
if max(real(lambda))>0;error('Problem');end
 
% sqrt(eig(Kci,Mci))/2/pi
i1=[1:size(Tci,2) size(Tci,2)*2+(1:size(Tqp,2)) size(Tci,2)*2+size(Tqp,2)*2+(1:size(Tqs,2))];
def=struct('def',psi(i1,:), ...
    'data',[abs(lambda) -real(lambda)./abs(lambda)],'DOF',(1:length(i1))'+.99, ...
    'TR',struct('def',[Tci Tqp Tqs],'DOF',d1.DOF));
[def.data,i2]=sortrows(def.data);def.def=def.def(:,i2);

out=def;

%% #SolveEnd
else;error('Solve%s',CAM);
end

    
elseif comstr(Cam,'pcond')
%% #Pcond coef : rescales pressure DOFs for better conditionning -2
%  default coef is 1e8
% model=fe_case(model,'pcond','PML','p_pml(''Pcond'')');
[CAM,Cam,r1]=comstr('cond',[-25 2],CAM,Cam); if isempty(r1);r1=1e8;end 
 Case=evalin('caller','Case');model=evalin('caller','model');
 T=evalin('caller','T'); DOF=evalin('caller','DOF');
 pc=ones(size(Case.T,1),1);
if 1==1 
 %% Scale columns
 i1=fe_c(DOF,(22:36)'/100,'ind');
 pc=ones(length(DOF),1); pc(i1)=r1; pc=diag(sparse(pc)); 
 T=T*pc;
 if isfield(Case,'TIn'); warning('Case.TIn support not implemented for PML');end
 assignin('caller','T',T);return
elseif 1==1
 i1=fe_c(model.DOF,(22:36)'/100,'ind');
 pc(i1)=r1;
else
 i1=fe_c(model.DOF,([22:26 30])'/100,'ind'); pc(i1)=r1;
 i1=fe_c(model.DOF,([27:29 31:36])'/100,'ind'); pc(i1)=r1*10000;
end
 pc=diag(sparse(pc));
 
 T=pc*T;
 if isfield(Case,'TIn')
     Case.TIn=pc*Case.TIn; assignin('caller','Case',Case);
 end
 assignin('caller','T',T);
 if 1==2 % Field order for easier viewing
  DOF=evalin('caller','DOF');
  i1=[fe_c(DOF,(1:3)'/100,'ind');fe_c(DOF,(13:21)'/100,'ind');   
      fe_c(DOF,(22:36)'/100,'ind')];
  DOF=DOF(i1);T=T(:,i1);
  assignin('caller','DOF',DOF);assignin('caller','T',T);
 end
 
elseif comstr(Cam,'addpml')
%% #AddPML: meshing strategy -1
% model=femesh('testhexa8');mo1=p_pml('addPML',model,struct('Lp',[1 1 1  0 1 0],'pow',1,'a0',1000,'Lc',.2));

model=varargin{carg};carg=carg+1;
RO=varargin{carg};carg=carg+1; 
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
 'subtype(1#%g#"1 time, 2 frequency") ' ...
 'Form(0#%g#"PML form") ' ...
 'pow(2#%g#"defaut attenuation") ' ...
 'Np(5#%g#"number of layers if auto length") ' ...
 'MeshType(keep#%s#"keep standard or one") ' ...
 ],{RO,CAM}); Cam=lower(CAM);

if ~isfield(RO,'ProId');
    il=feutil('getil',model);RO.ProId=max(il(:,1)+1); 
end
if isfield(RO,'Dp') % Allow scalar Lp and vector Dp [1 1 1  0 1 0]
  RO.Lp=RO.Dp*RO.Lp;
end
if ~isfield(RO,'Lc')
    error('Missing mesh size check');
end
if ~isfield(model,'unit');error('Missing model.unit');end
% .Lp [Lx-,Ly-,Lz-,Lx+,Ly+,Lz+] 
if isfield(RO,'SetPro')
%% just set properties
n1=feutil(['getnode ' RO.SetPro],model);
RO.box=[min(n1(:,5:7)) max(n1(:,5:7))];

else
RO.box=[min(model.Node(:,5:7)) max(model.Node(:,5:7))];
RO.nElt=size(model.Elt,1); RO.OnePml=[];RO.PmlNode=[];
for j1=reshape(find(RO.Lp),1,[])
 RO.Elt=model.Elt;
 r1=feutil('mpid',RO.Elt(1:2,:)); % material used for extrusion
 r1=feutil(sprintf('getmat%i struct',r1(2,1)),model);
 cp=sqrt(r1.E.*(1-r1.Nu)./r1.Rho./(1+r1.Nu)./(1-2*r1.Nu));%pressure wave speed

 if j1==1;     model.Elt=feutil('selelt selface&innode {x==}',model,RO.box(j1)); r1=[-1 0 0];
 elseif j1==2; model.Elt=feutil('selelt selface&innode {y==}',model,RO.box(j1));r1=[0 -1 0];
 elseif j1==3; model.Elt=feutil('selelt selface&innode {z==}',model,RO.box(j1));r1=[0 0 -1];
 elseif j1==4; model.Elt=feutil('selelt selface&innode {x==}',model,RO.box(j1));r1=[1 0 0];
 elseif j1==5; model.Elt=feutil('selelt selface&innode {y==}',model,RO.box(j1));r1=[0 1 0];
 elseif j1==6; model.Elt=feutil('selelt selface&innode {z==}',model,RO.box(j1));r1=[0 0 1];
 end
 if strcmpi(RO.MeshType,'one')
  %% single coarse quad on edge
  if ~isempty(RO.PmlNode)
   [model.Elt,elt]=feutil('removeelt innode',model,RO.PmlNode);
  else; elt=[];
  end
  n1=feutil('getnode groupall',model);n1(:,5:7)=round(n1(:,5:7)*1e6);
  i1=find(any(diff(n1(:,5:7))))+4;
  r2=[min(n1(:,i1));max(n1(:,i1))];n2=r2([1 3;1 4;2 4;2 3]);
  [un1,i2]=ismember(n2,n1(:,i1),'rows');
  model.Elt=feutil('addelt','quad4',[n1(i2)' model.Elt(2,5:6)]);%
  r1b=(n1(:,i1)-repmat(n1(i2(1),i1),size(n1,1),1));
  r2=r1b(i2(2),:)';r2=r2/norm(r2)^2;r3=r1b(i2(4),:)';r3=r3/norm(r3)^2;
  r1b=r1b*[r2 r3];
  % now create the MPC for EdgeDisp
  n2=n1;n2(i2,:)=[];r1b(i2,:)=[]; r1b(:,4)=0; EC=integrules('quad4',r1b*2-1);
  z=EC.N*0;
  r2=struct('DOF',[feutil('getdof',(1:3)'/100,n2(:,1));
          feutil('getdof',(13:21)'/100,n1(i2))], ...
      'c',[eye(3*size(n2,1)) -[EC.N z z   EC.N  z z  EC.N  z z;
        z EC.N z z   EC.N  z z  EC.N  z;z z EC.N z z   EC.N  z z  EC.N ]], ...
      'slave',1:size(n2,1)*3);
   %[fe_c(r2.DOF) num2cell(r2.c')]
  if isempty(RO.OnePml);RO.OnePml=r2;
  else;
   i3=fe_c(r2.DOF(r2.slave),RO.OnePml.DOF,'ind');
   r2.c(i3,:)=[]; r2.slave(i3)=[]; RO.OnePml=fe_mpc('mpcMerge',RO.OnePml,r2);
  end
  if ~isempty(elt);
      model=feutil('addelt',model,elt);model=feutil('joinall',model);
  end
 end
 if isfield(RO,'Freq') % Frequency based meshing
  RO.Lp(j1)=sign(RO.Lp(j1))*cp/min(RO.Freq(:))/2;% RO.Lp is    
  r2=[0 logspace(log10(min(cp/max(RO.Freq(:)),RO.Lc)),...
        log10(cp/min(RO.Freq(:))/2),RO.Np)]*sign(RO.Lp(j1));
  while diff(r2(2:3))<r2(2)
   RO.Np=RO.Np-1;
   r2=[0 logspace(log10(min(cp/max(RO.Freq(:)),RO.Lc)),...
        log10(cp/min(RO.Freq(:))/2),RO.Np)]*sign(RO.Lp(j1));
  end
  
 else
  r2=feutil(sprintf('refineline %.15g',RO.Lc),[0 RO.Lp(j1)]);
 end
 model=feutil(sprintf('extrude 0 %g %g %g',r1),model,r2);
 RO.PmlNode=[RO.PmlNode;feutil('findnode groupall',model)];
 model=feutil('addelt',model,RO.Elt);
 
 if ~isfield(RO,'a0')||RO.a0<=0
 %% Default a0 see (4.30) in Hadrien's
   RO.a0=-log(.001)*3/2*max(cp)/norm(RO.Lp,'inf');
   fprintf('p_pml set a0=%.2f, cp =%g Lp= %g\n',RO.a0,cp,norm(RO.Lp,'Inf'));
   % R1.Lp=(0:.15:5*.15); 
 end
 
end
if size(model.Elt,1)==RO.nElt;
    error('No PML added : wrong box in model.Node ?');
end
end
if ~isfield(RO,'a0')
   RO.a0=1000; warning('Missing a0 init');
end
if isfield(RO,'OnePml')&&~isempty(RO.OnePml); 
  model=fe_case(model,'Mpc','OnePml',RO.OnePml);
end

mpid=feutil('mpid',model);
% 6 faces,12 edges,8 corners 
RO.il=RO.ProId+(1:26)';RO.il(:,2)=fe_mat('p_pml',model.unit,RO.subtype);
RO.il(:,4)=RO.pow; RO.il(:,5)=RO.a0; RO.il(:,3)=RO.Form;

st={'x<','y<','z<','x>','y>','z>'}; % 6 faces
for j1=1:length(st)
 mpid(feutil(sprintf('findelt withnode {%s}',st{j1}),model,RO.box(j1)),2)=RO.il(j1);
 ind=abs(st{j1}(1))-abs('x')+1;in2=ind+3*double(st{j1}(2)=='>');
 RO.il(j1,5+ind*2+[-1 0])=[RO.box(ind+3*double(st{j1}(2)=='>')) RO.Lp(in2)];
end
% 12 edges
st={'x<&z<',[1 3];'y<&z<',[2 3];'x>&z<',[4 3];'y>&z<',[5 3]; % bottom edges
    'x<&y>',[1 5];'x<&y<',[1 2];'x>&y<',[4 2];'x>&y>',[4 5]; % side edges
    'x<&z>',[1 6];'y<&z>',[2 6];'x>&z>',[4 6];'y>&z>',[5 6]; % top edges
};
for j1=1:size(st,1)
 mpid(feutil(sprintf('findelt withnode {%s}',st{j1,1}),model, ...
     RO.box(st{j1,2}(1)),RO.box(st{j1,2}(2))),2)=RO.il(j1+6);
 ind=abs(st{j1}(1))-abs('x')+1;in2=ind+3*double(st{j1}(2)=='>');
 RO.il(6+j1,5+ind*2+[-1 0])=[RO.box(ind+3*double(st{j1}(2)=='>')) RO.Lp(in2)];
 ind=abs(st{j1}(4))-abs('x')+1;in2=ind+3*double(st{j1}(5)=='>');
 RO.il(6+j1,5+ind*2+[-1 0])=[RO.box(ind+3*double(st{j1}(5)=='>')) RO.Lp(in2)];
end
% 8 corners
st={'x<&y<&z<';'x<&y>&z<';'x>&y<&z<';'x>&y>&z<';
    'x<&y<&z>';'x<&y>&z>';'x>&y<&z>';'x>&y>&z>';
};
for j1=1:size(st,1)
 ind=1:3;if ~isempty(strfind(st{j1},'x>'));ind(1)=4;end
 if ~isempty(strfind(st{j1},'y>'));ind(2)=5;end
 if ~isempty(strfind(st{j1},'z>'));ind(3)=6;end
 mpid(feutil(sprintf('findelt withnode {%s}',st{j1}),model, ...
     RO.box(ind(1)),RO.box(ind(2)),RO.box(ind(3))),2)=RO.il(j1+18);
 ind=1;RO.il(18+j1,5+ind*2+[-1 0])=[RO.box(ind+3*double(st{j1}(2)=='>')) RO.Lp(ind)];
 ind=2;RO.il(18+j1,5+ind*2+[-1 0])=[RO.box(ind+3*double(st{j1}(5)=='>')) RO.Lp(ind)];
 ind=3;RO.il(18+j1,5+ind*2+[-1 0])=[RO.box(ind+3*double(st{j1}(8)=='>')) RO.Lp(ind)];
end
model.Elt=feutil('mpid',model,mpid);
RO.il(~ismember(RO.il(:,1),mpid(:,2)),:)=[];
if ~isfield(model,'il');model.il=[];
elseif ~isempty(model.il)
 model.il(ismember(model.il(:,1),RO.il(:,1)),:)=[];
end
if RO.subtype==2% subtype 2 freq
 for j1=1:size(RO.il,1)
  r2=struct('name','pml','il',RO.il(j1),'unit','SI','type','p_pml','a',RO.a, ...
      'Form',RO.Form);
  model=stack_set(model,'pro',sprintf('PML%i',RO.il(j1)),r2);
 end
end
model.il(end+(1:size(RO.il,1)),1:size(RO.il,2))=RO.il;
model.Elt=feutilb('SeparatebyProp -max 2000000',model.Elt);

if isfield(RO,'quad')&&RO.quad==1;model=feutil('lin2quad',model);end

if RO.subtype==1
 %% Constraint on q=qx+qy+qz for time PML
 model=p_pml('DispMpc',model,RO); 
end
if nargout==0;cf=feplot(model);fecom(cf,'ShowFiPro');
else; out=model;
end

%% #DispMPC : add MPC for q=qx+qy+qz and possibly for Sigma ------------------
elseif comstr(Cam,'dispmpc');
    
model=varargin{carg};carg=carg+1;
if carg<=nargin;RO=varargin{carg};carg=carg+1;else; RO=struct;end
il=feutil('getil',model);
il=il(cellfun(@(x)strcmpi(fe_mat('typep',x),'p_pml'),num2cell(il(:,2))),:);
if ~isfield(RO,'il'); RO.il=il;end
n1=feutil('findnode proid',model,il(:,1)); % PML nodes


if 1==1
%% Now create a stress continuity constraint on interface
%  and displacement components on edge
st=sprintf('%i ',RO.il(:,1));
st1=sprintf('selelt proid~= %s & withnode {proid%s}',st,st);elt=feutil(st1,model);
n2=feutil(sprintf('getnode proid%s & proid~=%s',st,st),model);
mo1=model;mo1.Elt=elt;
mo3=fe_caseg('stresscut',struct('type','conform','sel',st1(7:end)),model);
mo2=model;mo2.Elt=mo3.Elt;mo3.Node=n2;mo3.Elt=feutil('addelt','node1',n2(:,1));
i1=setdiff(feutil('proid',mo2.Elt),0);
mo2=feutil('setpro',mo2,[i1 fe_mat('p_solid','SI',1) 0 -3 0 100]);% il(6) for gradient
i2=setdiff(feutil('matid',mo2.Elt),0);mo2.pl=m_elastic(sprintf('dbval %i Strain',i2(1)));
dd=feutil('getdd',[i2(1) i1 3 1],model);dd=dd.dd;
%cut=fe_caseg('stresscut -selout',mo3,mo2); % Observation elastic side
%Be=cut.StressObs;Be.cta(abs(Be.cta)<1e-10)=0; % Observation elastic side

% Observation PML side 
mo2=model; mo2.Elt=feutil(['selelt proid' st],mo2);
mo2.Elt=feutil('set groupall matid 1 proid 1',mo2);
mo2=feutil('setpro',mo2,[1 fe_mat('p_solid','SI',1) 0 -3 0 100]);% il(6) for gradient
mo2.pl=m_elastic('dbval 1 Strain'); ofact('silent');
cut=fe_caseg('stresscut -selout',mo3,mo2); Bp=cut.StressObs;Bp.cta(abs(Bp.cta)<1e-10)=0;

r1=unique(round([Bp.DOF;feutil('getdof',n1,[1:3 13:21 22:36]'/100)]*100))/100;
r1=struct('c',[],'DOF',r1);Bp=feutilb('placeindof',r1.DOF,Bp);
%Be=feutilb('placeindof',r1.DOF,Be);

Bx=sparse([1 5 6],[1 4 7],1,6,9); By=sparse([6 2 4],[2 5 8],1,6,9);
Bz=sparse([5 4 3],[3 6 9],1,6,9); 
%DBe=reshape(dd*(Bx+By+Bz)*reshape(Be.cta,9,[]),6*size(Be.Node,1),[]);
ix=[1:3 5 6];iy=[1:4 6]; iz=1:5;
DBx=reshape(dd(ix,:)*Bx*reshape(Bp.cta,9,[]),5*size(Bp.Node,1),[]);
DBy=reshape(dd(iy,:)*By*reshape(Bp.cta,9,[]),5*size(Bp.Node,1),[]);
DBz=reshape(dd(iz,:)*Bz*reshape(Bp.cta,9,[]),5*size(Bp.Node,1),[]);
%ind=reshape(1:size(DBe,1),6,[])';in2=1:size(n2,1);
% [U,s,v]=svd(DBy(10+(1:5),:),0);U(:,diag(s)/s(1)>1e-5)
if isfield(RO,'CycBuild'); % Support single direction periodicity
 model=fe_cyclic(RO.CycBuild,model);
 cyc=fe_case(model,'getdata','Symmetry');
 %in2=~ismember(n2(:,1),cyc.IntNodes(:));ind=ind(in2,:);
else; cyc=[];
end

r1.c={fe_c(r1.DOF,feutil('getdof',n1,(1:3)'/100))- ...
 fe_c(r1.DOF,feutil('getdof',n1,(13:15)'/100))- ...
 fe_c(r1.DOF,feutil('getdof',n1,(16:18)'/100))- ...
 fe_c(r1.DOF,feutil('getdof',n1,(19:21)'/100));% q-qx-qy-qz
-DBx+fe_c(r1.DOF,feutil('getdof',n2(:,1),(22:26)'/100)) % DBx q - sig_x
-DBy+fe_c(r1.DOF,feutil('getdof',n2(:,1),(27:31)'/100)) % DBy q - sig_y
-DBz+fe_c(r1.DOF,feutil('getdof',n2(:,1),(32:36)'/100))}; % DBz q - sig_y

%  -DBe(ind(:,1:3),:)+fe_c(r1.DOF,feutil('getdof',n2(in2,1),(22:24)'/100))+...
%    fe_c(r1.DOF,feutil('getdof',n2(in2,1),(27:29)'/100))+...
%    fe_c(r1.DOF,feutil('getdof',n2(in2,1),(32:34)'/100)) % sigma_e - sigma x- sigma y-sigma z
%  -DBe(ind(:,4),:)+fe_c(r1.DOF,feutil('getdof',n2(in2,1),.30))+...
%    fe_c(r1.DOF,feutil('getdof',n2(in2,1),.35)) % sigma_e - sigma x- sigma y-sigma z
%  -DBe(ind(:,5),:)+fe_c(r1.DOF,feutil('getdof',n2(in2,1),.25))+...
%    fe_c(r1.DOF,feutil('getdof',n2(in2,1),.36)) % sigma_e - sigma x- sigma y-sigma z
%  -DBe(ind(:,6),:)+fe_c(r1.DOF,feutil('getdof',n2(in2,1),.26))+...
%    fe_c(r1.DOF,feutil('getdof',n2(in2,1),.31)) % sigma_e - sigma x- sigma y-sigma z
%     };
 %%  continuity and elimination of stress on interface
 
 rd=r1;[rd.c,rd.slave]=feutil('fixMpcMaster',vertcat(r1.c{1}));%fe_c(r1.DOF(rd.slave))
 % 2. eliminate interface stresses s_xyz using constitutive law since no delay there
 if 1==2
  r3.c=vertcat(r1.c{2:4}); r3.c=r3.c-sparse(r3.c(:,rd.slave))*rd.c; r3.c(:,rd.slave)=0;
  [r3.c,r3.slave]=feutil('fixMpcMaster',r3.c);%fe_c(r1.DOF(r3.slave))
  r2=r1;r2.c=[rd.c;r3.c]; r2.slave=[rd.slave;r3.slave];%sparse( r2.c(:,r2.slave)-speye(length(r2.slave)))
  n1=setdiff(n1,n2(:,1)); % Then eliminate null strains only not on interface
 else;r2=rd;
 end
 model=fe_case(model,'mpc','PmlDispMpc',r2); 
  
 % e_yy^x = inv(dd)(yy,:)*sx =0 
 % e_zz^x = inv(dd)(zz,:)*sx =0 
 r1=inv(dd(1:3,1:3));  % Interior Nodes
 i2=2:3; r3=r1(i2,i2)\r1(i2,:);r3(:,i2)=eye(length(i2));r3=r3(1);
 %i2=[1 3]; r2=r1(i2,i2)\r1(i2,:);r2(:,i2)=eye(length(i2));r2
 %i2=[1 2]; r2=r1(i2,i2)\r1(i2,:);r2(:,i2)=eye(length(i2));r2
 
 r1=struct('c',[], ...
    'DOF',feutil('getdof',n1,[22:24   28 27 29   34 32 33  ...
        25 36 26 31 30  35]'/100));% DofLabT(1+[25 36 26 31 30  35]'/100)
 r1.c=[fe_c(r1.DOF,n1+.23)+r3*fe_c(r1.DOF,n1+.22)
       fe_c(r1.DOF,n1+.24)+r3*fe_c(r1.DOF,n1+.22)
       fe_c(r1.DOF,n1+.27)+r3*fe_c(r1.DOF,n1+.28)
       fe_c(r1.DOF,n1+.29)+r3*fe_c(r1.DOF,n1+.28)
       fe_c(r1.DOF,n1+.32)+r3*fe_c(r1.DOF,n1+.34)
       fe_c(r1.DOF,n1+.33)+r3*fe_c(r1.DOF,n1+.34)
       ];
 if 1==1 % No shear
    r1.c=[r1.c;
       fe_c(r1.DOF,n1+.25);fe_c(r1.DOF,n1+.36) % ax e_xy^x + ay e_xy^y
       fe_c(r1.DOF,n1+.26);fe_c(r1.DOF,n1+.31)
       fe_c(r1.DOF,n1+.30);fe_c(r1.DOF,n1+.35)
       ];   
 elseif 1==2 % Equal shear
    r1.c=[r1.c;
       fe_c(r1.DOF,n1+.25)-fe_c(r1.DOF,n1+.36) % ax e_xy^x + ay e_xy^y
       fe_c(r1.DOF,n1+.26)-fe_c(r1.DOF,n1+.31)
       fe_c(r1.DOF,n1+.30)-fe_c(r1.DOF,n1+.35)
       ];   
 elseif 1==2 % Oposite shear
    r1.c=[r1.c;
       fe_c(r1.DOF,n1+.25)-fe_c(r1.DOF,n1+.36) % ax e_xy^x + ay e_xy^y
       fe_c(r1.DOF,n1+.26)-fe_c(r1.DOF,n1+.31)
       fe_c(r1.DOF,n1+.30)-fe_c(r1.DOF,n1+.35)
       ];   
 end
 [r1.c,r1.slave]=feutil('fixMpcMaster',r1.c);
 model=fe_case(model,'mpc','PmlStressMpc',r1); 
 
 if  ~isempty(cyc) % Force choice of interface nodes
  cyc.IntDof=fe_c(fe_c(r2.DOF,r2.DOF(r2.slave),'dof',2),cyc.IntNodes(:,1),'dof');
  nind=sparse(cyc.IntNodes(:,1),1,cyc.IntNodes(:,2));
  cyc.IntDof(:,2)=full(nind(fix(cyc.IntDof(:,1))))+rem(cyc.IntDof(:,1),1);
  model=fe_case(model,'stack_set',{'cyclic','Symmetry',cyc}); 
 end
else
 % Disp MPC within the PML
 i1=[1 13 16 19 2 14 17 20 3 15 18 21]'/100;
 r1=struct('c',sparse(repmat(1:length(n1)*3,4,1), ... % row
    1:length(n1)*12,repmat([1;-1;-1;-1],1,length(n1)*3)), ...
    'DOF',feutil('getdof',n1,i1),'slave',1:4:length(n1)*12);
 model=fe_case(model,'Mpc','PmlDispMPC',r1);
end

out=model;

%% #View -------------------------------------------------------------------------
elseif comstr(Cam,'view');[CAM,Cam]=comstr(CAM,5);
if comstr(Cam,'topo');
%% #ViewTopo -------------------------------------------------------------------------

SE=varargin{carg};carg=carg+1;
CE=varargin{carg};carg=carg+1;

K=feutilb('tkt',CE.T,SE.K);
i1={1:3,13:15,16:18,19:21,22:26,27:31,32:36};
st={'q','qx','qy','qz','s..x','s..y','s..z'};
for j1=1:length(i1)
 i1{j1}=fe_c(CE.DOF,i1{j1}(:)/100,'ind');
end
i2=cellfun(@isempty,i1);disp(st(~i2));i1(i2)=[];
figure(1);clf;
for j1=1:3
  subplot(1,3,j1);ii_plp('spy',struct('K',K{j1},'ind',{i1}));
  title(SE.Klab{j1});
end
elseif comstr(Cam,'lublocks')
%% #ViewLuBlocks : diagnostics using LU blocks -------------------------
RO=varargin{carg};carg=carg+1;
[i2,i1]=sort(rem(RO.DOF,1));
if isfield(RO,'K');
 [L,U,p]=lu(RO.K);figure(11);semilogy(abs(diag(U(i1,i1))));
 d=abs(diag(U)); iz=RO.DOF(d<1e-2);fecom('shownodemark',iz)
elseif isfield(RO,'Ti')
  % Check of fields in fe_homo DftRedP2set
  % T2.DOF=SE.DOF(i2);p_pml('viewLuBlocks',T2)
  r1=abs(RO.Ti(i1,1:3)); r1(r1<.001*norm(r1(:),'inf'))=0;
  figure(11);semilogy(r1);  
  iz=RO.DOF(any(r1,2));fecom('shownodemark',iz)
end
 axis('tight');
 i3=find(diff(round(i2*100)));ii_plp(i3*[1 0]);%set(gca,'xtick',i3) 
 x=(i3+[0;i3(1:end-1)])/2;
 set(gca,'xtick',x,'xticklabel', ...
      cellfun(@(x)sprintf('%i',x),num2cell(unique(round(i2*100))),'uni',0)) 
  1;
elseif comstr(Cam,'redsvd')
%% #ViewRedSvd : analyze zero subspace

if carg<=nargin; SE=varargin{carg};carg=carg+1; 
    ki=SE.K;
else;
 PA=dyn_ui('paramvh');SE=PA.mt.Stack{1,3};
 ki=feutilb('sumkcoef',SE.K,[1e6 1 1e3]);
end
[u,s]=svd(full(ki)); d1=struct('def',fliplr(u),'DOF',SE.DOF,'data',flipud(diag(s)));
d1.data(1:10)
if isfield(SE,'TR'); d1.TR=SE.TR; cf=feplot(SE,d1);
 [i2,i1]=sort(rem(SE.TR.DOF,1));d2=fe_def('subdef',d1,1:4);
 d2.def=d1.TR.def*d2.def;d2.DOF=d1.TR.DOF;
 d2=feutilb('placeindof',SE.TR.DOF(i1),d2);d2=feutil('rmfield',d2,'TR');
 r2=abs(d2.def);r2(abs(r2)<1e-10)=0;r2=r2*diag(1./max(r2));
 figure(11);clf;semilogy(abs(r2))
 cf.def=d2; fecom('colordataA');fecom('coloralpha');% Not EvalA to see stress 
 i3=find(diff(round(i2*100)));ii_plp(i3*[1 0]);%set(gca,'xtick',i3) 
 x=(i3+[0;i3(1:end-1)])/2;
 set(gca,'xtick',x,'xticklabel', ...
      cellfun(@(x)sprintf('%i',x),num2cell(unique(round(i2*100))),'uni',0))

 iz=d2.DOF(any(r2>.5,2));fecom('textnode',unique(fix(iz)));fecom('shownodemark',iz);
 feval(p_pml('@DofLabT'),sort(iz)) 

else; cf=feplot;cf.def=d1;
end


 1;
elseif comstr(Cam,'q')
%% #ViewQ : view the various translation fields

cf=feplot;d2=cf.Stack{'def'};
if isempty(d2);stack_set(cf,'curve','def',cf.def);d2=cf.def; end
d2=feutil('rmfield',d2,'scale','opt');

% Qx
if Cam(2)=='x'
 RO.d=.13; RO.e=([1:3 16:21 27:36])'/100; 
 RO.f=(13:15)'/100; RO.fi=.12; RO.g=(22:26)'/100; RO.gi=0;
elseif Cam(2)=='y'
 RO.d=.16; RO.e=([1:3 13:15 19:21 22:26 32:36])'/100; 
 RO.f=(16:18)'/100; RO.fi=.15; RO.g=(22:26)'/100;   RO.gi=.05; 
else
 % Qz
 RO.d=.19; RO.e=([1:3 13:18 22:31])'/100; 
 RO.f=(19:21)'/100; RO.fi=.18; RO.g=(32:36)'/100;RO.gi=.1;
end

i2=fix(fe_c(d2.DOF,RO.d,'dof'));
d3=fe_def('subdofind',d2,fe_c(d2.DOF,feutil('getdof',i2,RO.e),'ind',2));
i3=fe_c(d3.DOF,feutil('getdof',i2,RO.f),'ind');d3.DOF(i3)=d3.DOF(i3)-RO.fi;
i3=fe_c(d3.DOF,feutil('getdof',i2,RO.g),'ind');d3.DOF(i3)=d3.DOF(i3)-RO.gi;

cf.def=d3;


else; error('View%s',CAM);

end

% -------------------------------------------------------------------------
%% #Test basic test : p_pml('testElt')
elseif comstr(Cam,'test');[CAM,Cam]=comstr(CAM,5);

if comstr(Cam,'elt')    
%% #TestElt Do not integrate example into femesh -2
 
end
%% #End ----------------------------------------------------------------------
elseif comstr(Cam,'tablecall');out='';
elseif comstr(Cam,'@');out=eval(CAM);
elseif comstr(Cam,'cvs');
 out=sdtcheck('revision');
 %out='$Revision: 527 $  $Date: 2020-10-21 19:10:15 +0200 (Wed, 21 Oct 2020) $'; return;
else;sdtw('''%s'' not known',CAM);
end
end

%% #SubFunc

%% #DfrfBIN Full/Reduced DRF with enforced shape - ---------------------------
function def1=DfrfBIn(MVR,mdl,oProp)

if 1==2
 % This is obsolete but could be tested with 
 mo1=d_pml('MeshUlb1D SW',struct('v',1,'quad',1,'Lc',2));feutilb('_write',mo1)
 RF=struct('DfrfBIn',p_pml('@DfrfBIn'));
 d1=fe_simul('dfrf',stack_set(mo1,'info','Freq',[10;100]),RF);d_pml('View1DPS',mo1,d1);
 % Pressure wave wave
mo2=d_pml('MeshUlb1D',struct('quad',0,'Lc',2));
RF=struct('DfrfBIn',p_pml('@DfrfBIn'));
d2=fe_simul('dfrf',stack_set(mo2,'info','Freq',[10;50]),RF);
d_pml('View1DPS',mo2,d2);

end

  warning('Obsolete use implicit matrix now'); v=3;
  Case=evalin('caller','Case');EC=Case.GroupInfo{1,8};

  na=size(MVR.BIN{1},2);if ~isfield(MVR,'br');MVR.br=[];end
  na=na+size(MVR.br,2);
  BIN=MVR.BIN;K=MVR.K;w=MVR.w;
  if isfield(MVR,'TIn')&&size(MVR.TIn,2)<na; MVR.TIn(1,na)=0;end

  def1=struct('def',zeros(size(MVR.DOF,1),length(w)*na), ...
   'DOF',MVR.DOF,'data',w/2/pi,'Xlab',{{'DOF','Freq'}});
  if na>1; 
    def1.data=[reshape(repmat(def1.data',na,1),[],1) ...
     repmat((1:na)',length(w),1)];def1.Xlab{2}={'Freq';'In'};
   end
  [zCoef,K]=feval(fe_simul('@safe_dfrf_zcoef'),w,'K',mdl);
  if isfield(MVR,'lab_in');def1.lab_in=MVR.lab_in;end
  if v==3;  elseif v==2 % xxx Will need to implement tkt
    obs=EC.obs;obs.cta=obs.cta*Case.T;obs.DOF=Case.DOF;
    obs.trans.cta=obs.trans.cta*Case.T;obs.trans.DOF=Case.DOF;
    K{1,end+1}=vhandle.matrix('ZMatrix',obs);zCoef(:,end+1)=w; 
  else
   obs=EC.obs;obs.cta=obs.cta*Case.T;obs.DOF=Case.DOF;
   obs.trans.cta=obs.trans.cta*Case.T;obs.trans.DOF=Case.DOF;
  zCoef(:,end+1)=w; K{1,end+1}=@(x)obs.toK(obs,x);
  end

  for j1=1:length(w)
   %  zCoef=fix_loss(mdl,zCoef);
   %[K{length(zCoef)+1},K{length(zCoef)+2}]=obs.StretchD(obs,w(j1)/2/pi);
   %Z=feutilb('sumkcoef',K,[zCoef(j1,:) -w(j1)^2 1]);% append pml m,k
   Z=feutilb('sumkcoef',K,zCoef(j1,:));% append pml m,k
   if ~isempty(MVR.br);b=[feutilb('sumkcoef',BIN,zCoef(j1,:)) MVR.br];
   else;b=feutilb('sumkcoef',BIN,zCoef(j1,:));
   end
   if MVR.full
    r1=ofact(Z,b);
   else;r1=(Z\b);
   end
   if ~isempty(MVR.T);r1=MVR.T*r1+MVR.TIn; % Restitution on all DOFs
   else; r1(end+1)=1; r1=MVR.Sens.cta*r1;
   end
   def1.def(:,(j1-1)*na+(1:na)) = r1;
  end

end

%% Add a matrix in proper location
function K=addMat(K,M,ix,iy);

[II,JJ,KK]=find(M);
K=K+sparse(ix(II),iy(JJ),KK,size(K,1),size(K,2));
end

%% #ConstitLab1 case of time PML (developped with Hadrien Pinault)
function out=ConstitLab1
out={'-1','il(2)','form','pow','a0','x0','Lx0','y0','Ly0','z0','Lz0', ...
        'pl(1)','pl(2)','rho', ...% 12:14
        'd11','d12','d13','d22','d23','d33','d44','d55','d66', ... % Non attenuated
        'd11x','d12x','d13x','d22x','d23x','d33x','d44x','d55x','d66x', ... % atten x : 24:32
        'd11y','d12y','d13y','d22y','d23y','d33y','d44y','d55y','d66y', ... % atten y (33:41)
        'd11z','d12z','d13z','d22z','d23z','d33z','d44z','d55z','d66z', ... % atten z (42:50)
        'rax','ray','raz','1' %51:54
        };
end    
    
%% #ConstitLab2 case of freq PML (developped with Arnaud Deraemaeker)
function out=ConstitLab2
out={'-1','il(2)','form','pow','a0','x0','Lx0','y0','Ly0','z0','Lz0', ...
        'b0','pl(1)','pl(2)','rho', ...% 12:15
        'd11','d12','d13','d22','d23','d33','d44','d55','d66'}; % 16:24 Non attenuated
end
    
function out=ViscousDampingBoundary(model,RM);
%% #ViscousDampingBoundary #VDB Method INFINITE DAMPING (Viscous Damping Boundary Method)

%% #both top and side
%%%%%%%%%%%%%%%%%%%
% cyl <=r i nx 

data=struct('sel',RM.SelNode, ... 
             'eltsel',['Withnode{ ' RM.SelNode '}'],'def',1,'DOF',.19);
model=fe_case(model,'Fsurf','Surface load',data);
Load = fe_load(model);

% Arnaud - to check if nodes are identified correctly
if 1==2
  nd=unique(floor(Load.DOF(Load.def~=0)));feplot(model);fecom('textnode',nd)
 % @cd this is a strange load distribution no ?
 feplot(model,Load); %
 %fecom('showdefarrow') % to see the load as vectors
 % pause
end

% Node to apply cbush
[~,NR] = feutil(['findnode ' RM.SelNode ],model); % Ok the nodes are fine

Elfact = [10 0];
elem_data=[];

model.bas=[]; 
if ~isfield(RM,'ki');RM.ki=1e2*ones(1,6);end % Default spring stiffness

for i_n =1:length(NR);
    % Position/direction actual spring
%     NR(i_n,1);
    ind = fe_c(Load.DOF,NR(i_n,1),'ind');
%     Load.DOF(ind);
    F1 = Load.def(ind);
    
    if abs(sum(F1))~=0
        if size(F1,1)==2
        dir = [(F1')./norm(F1) 0]; % Added a zero because only x and y load for some elements ...
        else
        dir = (F1')./norm(F1);
        end
            
        fact = -sign(NR(i_n,5:7)*dir')*norm(F1);
        %CorID Type 0     Ax Ay Az       Ux Uy Uz Vx Vy Vz Wx Wy Wz s
        p =  basis(dir,[0 0 0],1);
        model.bas(end+1,:) = [i_n 1 0 NR(i_n,5:7) [p(:,1)' p(:,2)' p(:,3)'] 1] ;
        Elfact(end+1,:) = [Elfact(end,1)+1 fact];
        ProId = Elfact(end,1);
        % spring element definition
        elem_data(end+1,:)=[[NR(i_n,1) 0] [ProId ProId 0] 0 0 0 i_n 1];
        % Spring properties
        model.il(end+1,1:14)=[ProId fe_mat('p_spring','SI',2) RM.ki fact*RM.ci];
        %     end
    end
end
clear load NR 
model.Elt=feutil('AddElt',model.Elt,'cbush',elem_data);
model=fe_case(model,'Remove','Surface load');
out=model; 

end
%% #stretch : evaluate stretch function -2
function lx=stretch(x,xiL,pow,fei,Vslw)
    if xiL(2)==0; % No stretch direction
        lx=1;
    else
     rpos=abs((x-xiL(1))/xiL(2))^pow;% Relative position within PML
     lx=1+complex(real(fei),Vslw*imag(fei))*rpos;
    end
end
%% #stretchD : vectorized streching
% xxx missing ref to PML.tex equations
% should actually be a single Z matrix
function [m,k]=stretchD(obs,freqRad) 

freq=freqRad/2/pi;
fei=interp1(obs.a.X{1},obs.a.Y,freq,'linear','extrap')*[1;-1i];
% build the block-diagonal DD
DK=zeros(81,size(obs.Node,1));IK=DK;JK=DK;
ik=(1:9)'*ones(1,9);jk=reshape(ik',[],1);ik=ik(:);
DM=zeros(3,size(obs.Node,1));
im=(1:3)';jm=im;
for jw=1:size(obs.Node,1)
  if size(obs.xiL,2)==1;i3=1;else;i3=jw;end
  lx=stretch(obs.Node(jw,5),obs.xiL(:,i3),obs.pow,fei,obs.Vs/obs.xiL(2)/freqRad);
  ly=stretch(obs.Node(jw,6),obs.yiL(:,i3),obs.pow,fei,obs.Vs/obs.yiL(2)/freqRad);
  lz=stretch(obs.Node(jw,7),obs.ziL(:,i3),obs.pow,fei,obs.Vs/obs.ziL(2)/freqRad);
  Lambda=[1/lx 0 0   0 0 0      0 0 0;
           0 0 0      0 1/ly 0   0 0 0;
           0 0 0      0 0 0      0 0 1/lz;
           0 0 0      0 0 1/lz  0 1/ly 0;
           0 0 1/lz   0 0 0     1/lx 0 0;
           0 1/ly 0   1/lx 0 0  0 0 0;
           ];
  if jw==1&&any(~isfinite(Lambda(:))); error('Problem with data');end
  Ce=(Lambda.'*obs.C*Lambda)*lx*ly*lz*obs.wjdet(jw);
  DK(:,jw)=Ce(:);IK(:,jw)=ik(:);JK(:,jw)=jk(:);ik=ik+9;jk=jk+9;
  DM(:,jw)=obs.rho*lx*ly*lz*obs.wjdet(jw)*[1;1;1];
  IM(:,jw)=im(:);JM(:,jw)=jm(:);im=im+3;jm=jm+3;
end
DK=sparse(IK,JK,DK);
DM=sparse(IM,JM,DM);
m=obs.trans.cta'*DM*obs.trans.cta;
k=obs.cta'*DK*obs.cta;
if nargout==1; % Return Z for implicit/zmatrix implementation
 m=k-(freqRad)^2*m; 
end
end

%% #DofLabT
function out=DofLabT(i1);

ind=[1:6 13:36];
%nind=sparse(nind,1,1:length(nind))
st=[{'u','v','w','rx','ry','rz'}  ...
     {'qx_x','qy_x','qz_x','qx_y','qy_y','qz_y','qx_z','qy_z','qz_z' ...
      'Sxx_x','Syy_x','Szz_x','Szx_x','Sxy_x', ... % strain 10:14 DOF 22:26
      'Sxx_y','Syy_y','Szz_y','Syz_y','Sxy_y', ... % strain 15:19 DOF 27:31
      'Sxx_z','Syy_z','Szz_z','Syz_z','Szx_z'}];
if nargin==0
 if nargout==1
   out=st;
 else %feval(p_pml('@DofLabT'));
  disp([st(:) num2cell(ind(:))])   
 end
elseif all(i1<1) % Used labels feval(p_pml('@DofLabT'),i2/100)
 i1=round(i1*100);
 nind=sparse(ind,1,1:length(ind));
 out=[st(nind(i1))' num2cell(i1(:))];
else;
 i2=round(remi(i1,1)*100);
 nind=sparse(ind,1,1:length(ind));
 out=[num2cell(fix(i1)) st(nind(i2))' num2cell(i2(:))];
 if nargout==0; out=out';fprintf('%i %s(%i)\n',out{:});clear out;end
end
end


%% #nl_inout : euler integration of stresses
function out=nl_inout(NL,fc,model,u,v,a,opt,Case)

% M_ss (sig_(n+1) - sig_n)/dt + C_sq q_n + C_ss sig_n = 0 
% sig_(n+1)= sig_n - (M_ss \ ( Csq q_n + C_ss sig_n ))*dt
% sig_(n+1)= sig_n - (c q_n)*dt - (Ca sig_n )*dt 
% unl(:,1,1) contains (M_ss/dt) \ ( Csq q_n )
% unl(:,1,3) contains sig_n 

% Ms=model.K{1}(is,is);Ms=Ms+spalloc(size(Ms,1),size(Ms,2),0);
% NL=struct('type','nl_inout', 'c',Ms\model.K{2}(is,iq),'b',model.K{2}(iq,is));
% %NL.Ca=full(diag(Ms\model.K{2}(is,is)));
% if nnz(Ms-diag(diag(Ms)))>0
%  NL.Ca=Ms\model.K{2}(is,is); % figure(1);spy(model.K{2}(is,is))
% else; NL.Ca=diag(Ms\model.K{2}(is,is));
% end

unlj1=NL.unl(:,:,3);
cunl=NL.unl(:,:,1); % M_ss \ ( Csq q_n)

if size(NL.Ca,2)==1
 sig = unlj1 + ( cunl + NL.Ca.*unlj1)*(-opt.Opt(4));
else
 sig = unlj1 + (cunl + NL.Ca*unlj1)* (-opt.Opt(4));
end

i0=of_time(-1,NL.unl,sig*[1 0 1]);% return sigma and shift unlj1 
out=struct('snl',sig);

end
