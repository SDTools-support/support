function out=fe_honeycomb(varargin)

% FE_HONEYCOMB support function of honeycomb modeling
%
% fe_honeycomb('matcore') % orthotropic core model (Hashbi formulas) 


%       Etienne Balmes, Jean-Michel Leclere
%       Copyright (c) 1990-2025 by SDTools, All Rights Reserved.
%       For revision information use fe_honeycomb('cvs')

if nargin==0; CAM=''; carg=1; else;CAM=varargin{1}; carg=2;  end
[CAM,Cam]=comstr(CAM,1);

% -------------------------------------------------------------------------
%    mat300=fe_honeycomb('matcore',model)
%    mat300=fe_honeycomb('matcore',RO)
if comstr(Cam,'matcore');[CAM,Cam]=comstr(CAM,8);%

if isempty(Cam)&&nargin==1 % Test
  CAM=['MatCore a .0027 b .0027 t 7.62e-5 t2 1.524e-4' ...
      'E 3e9 nu .3 rho 1380'];
  %fe_honeycomb(CAM)
  [CAM,Cam]=comstr(CAM,8);
end

    opt={'a(#%g#"Cell L wall length") ', ...
         'b(#%g#"Cell W wall length") ', ...
         'theta(30#%g#"W wall angle (deg)") ', ...
         't (#%g#"simple wall thickness")', ...
         't2 (#%g#"double wall thickness")', ...
         'E (#%g#"double wall thickness")', ...
         'nu (#%g#"double wall thickness")', ...
         'rho (#%g#"double wall thickness")', ...
         'G (#%g#"double wall thickness")', ...
         'MatId (200#%i#"Material identifier")', ...
         'Split (0#%i#"Split components")', ...
      };

    
    model=[];if carg<=nargin;model=varargin{carg};carg=carg+1;end
    if isfield(model,'Elt');
        if carg>nargin;RO=stack_get(model,'info','CoreProp','getdata');
        else;RO=varargin{carg};carg=carg+1;
        end
    else; RO=model;model=[];
    end
    RO=cingui('paramedit-given',sprintf('%s ',opt{:}),{CAM,RO});

    RO=fe_def('CleanEntry',RO); st=fieldnames(RO); 
    for j1=1:length(st);eval(sprintf('%s=RO.%s;',st{j1},st{j1}));end

    %a=alveole(1);b=alveole(2);theta=alveole(3);
    theta=theta*pi/180;
    
    E1ref=E*(t/b)^3*((a/b+sin(theta))/(cos(theta))^3);
    E2ref=E*(t/b)^3*cos(theta)/((a/b+sin(theta))*...
        (sin(theta))^2);
    E3ref=E*(2*t+a/b*t2)/(a+b*sin(theta))/cos(theta);
    nu12ref=0.8; nu31ref=0.3; nu32ref=0.3;
    G12ref=E*(a/b+sin(theta))/(2*(a)^2*cos(theta))*1/(b/...
        (2*(t2)^3)+a/t^3);
    G1z_inf=E/2/1.3*(a*t2+b*t*sin(theta))^2/(b*cos(theta)...
        *(a+b*sin(theta))*(2*a*t2+b*t));
    G1z_sup=E/2/1.3*(a*t2+2*b*t*(sin(theta))^2) ...
        /(2*b*cos(theta)*(a+b*sin(theta)));
    G1zref=(G1z_inf+G1z_sup)/2;
    G2zref=E/2/1.3*cos(theta)/(a/b+sin(theta))*t/b;
    rhoref=rho*(t2*a+2*t*b)/(2*b*cos(theta)*(a+b*sin(theta)));
    coef=1;
    %  size_eltref=10e-3
    if ~isfield(RO,'E1'); RO.E1=[E1ref]; end
    if ~isfield(RO,'E2'); RO.E2=[E2ref]; end
    if ~isfield(RO,'E3'); RO.E3=[E3ref]; end
    if ~isfield(RO,'nu12'); RO.nu12=[nu12ref]; end
    if ~isfield(RO,'nu31'); RO.nu31=[nu31ref]; end
    if ~isfield(RO,'nu32'); RO.nu32=[nu32ref]; end
    if ~isfield(RO,'G12'); RO.G12=[G12ref]; end
    if ~isfield(RO,'G1z'); RO.G1z=[G1zref]; end
    if ~isfield(RO,'G2z'); RO.G2z=coef*[G2zref]; end
    if ~isfield(RO,'rho'); RO.rho=[rhoref]; end
      
    RO.nu21=RO.E2/RO.E1*RO.nu12;
    RO.nu13=RO.E1/RO.E3*RO.nu31;
    RO.nu23=RO.E2/RO.E3*RO.nu32;

    % Construction du tenseur de souplesse du nida homogene
    C=[1/RO.E1 -RO.nu21/RO.E2 -RO.nu31/RO.E3 0 0 0;...
        -RO.nu12/RO.E1 1/RO.E2 -RO.nu32/RO.E3 0 0 0;...
        -RO.nu13/RO.E1 -RO.nu23/RO.E2 1/RO.E3 0 0 0;...
        0 0 0 1/RO.G2z 0 0;...
        0 0 0 0 1/RO.G1z 0;...
        0 0 0 0 0 1/RO.G12];
    S=inv(C);fprintf('G1z %.2f (ref %.2f) G2z %.2f (ref %.2f) MPa\n', ...
        [S(5,5),G1zref,S(4,4),G2zref]/1e6);
    if any(eig(S)<0)
        disp(S); error('Not a positive constitutive law');
        %mat2=[];
    end
    % Orthotropic core material property 
    mat2=zeros(1,24);
    mat2([1,2,3,4,5,6,7,8,12,17,23,24])=[200,fe_mat('m_elastic','SI',3),...
        S(1,1),S(1,2),S(2,2),S(1,3),S(2,3),S(3,3),S(4,4),S(5,5),S(6,6),RO.rho];
    if RO.Split
     mat2=[[301;302;303],mat2([1;1;1],2),mat2([1 1 1],3:end)/3]; % Gij/3
     mat2(1,[12 17])=[S(4,4)/4 S(5,5)/4];% Mat 201 (17) G1z=G1z0 (12) G2z=G2z0
     mat2(2,[12 17])=[2*S(4,4)/4 S(5,5)/4];% Mat 202 G1z0 2*G2z0
     mat2(3,[12 17])=[S(4,4)/4 2*S(5,5)/4];% Mat 203 2*G1z0 G2z0
    end
 
    if ~isempty(model);   
       for j1=1:size(mat2,1)
           i1=[];
           try; i1=find(model.pl(:,1)==mat2(j1,1));end
           if isempty(i1);i1=size(model.pl,1)+1;end
           model.pl(i1,1:size(mat2,2))=mat2(j1,:);
       end
       out=model;
    else; out=mat2;
    end
 
    
% ---------------------------------------------------------------------
elseif comstr(Cam,'test');[CAM,Cam]=comstr(CAM,5);
    
% place a piezo patch using fe_fmesh
if comstr(Cam,'patch')
    
 mdl=femesh('test quad4 -divide 30 30');
 fe_fmesh patch
 mdl=fe_fmesh('patch Ori .15 .15 0 width .2 height .1 Thick .01 Divx 2 DivY3',mdl)

% - - - - - - - - - - - - - - - - - - - 
else;error('Test%s unknown',CAM);
end
    
% -----------------------------------------------------------------
elseif comstr(Cam,'cvs')
 out=sdtcheck('revision');
 %out='$Revision: 492 $  $Date: 2020-02-26 12:10:10 +0100 (Wed, 26 Feb 2020) $';
else;error('%s unknown',CAM);
end
