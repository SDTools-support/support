function [out,out1,out2]=t_subspace(varargin)

% Tests used for the development of subspace methods
%
% t_subspace('stab',ID) expects a structure with fields
% .po [freq damp order]
% .res resisues
%
% See <a href="matlab: sdtweb _taglist t_subspace">TagList</a>
% See example with t_subspace('poly')
% 
% Contributed by E. Balmes

%       Etienne Balmes
%       Copyright (c) 1990-2025 by SDTools, All Rights Reserved.


%#ok<*NOSEM,*ASGLU>
if nargin==0; Cam='';
else; [CAM,Cam]=comstr(varargin{1},1);carg=2; 
end

if isempty(Cam)  % These are tests that run automatically

[XF,Solution]=t_subspace('testone');

[XF,Solution]=t_subspace('testsdt_id');
% Frequency domain polyreference
ID=t_subspace('polyref',XF,20);t_subspace('stab',ID);
figure(1);ii_plp(Solution.po)

% Frequency domain orthogonal polynomial
RO=struct('ind',3:15);
ID=t_subspace('poly',XF,RO);
t_subspace('stab',ID);

if 1==2 
 %% #Cross Identification from cross spectra
 % edit(sdeb('wdToDo\corus\t_corus.m'))
 [XF,Solution]=t_subspace('testone',struct('GxxRef',1));Solution.F=[6 16];
 %[XF,Solution]=t_subspace('testsdt_idAll',struct('GxxRef',1));Solution.F=[30 230];
 XF.Edit=struct('GenFcn','process_r(''@PosSpec'')', ...
     'Fmin',struct('type','double','value',Solution.F(1)), ...
     'Fmax',struct('type','double','value',Solution.F(2)), ...
     'Generate',struct('type','push','value','Test', ...
     'callback',{{'ii_mmif'}}));
 ci=iiplot;iicom(ci,'curveInit','Cross',XF);iicom('curtabStack','Cross');
 iicom('curveinit','Test',ii_mmif(XF));iicom submagpha
 ci=iiplot;ci.Stack{'IdMain'}=Solution;idcom est
 %idcom(sprintf('e .01 %g',Solution.po(1)));idcom('ea');idcom('est');
 %idcom('e .01 10');idcom('ea');idcom('est');
 ci.Stack{'IdMain'}
 
 
end

elseif comstr(Cam,'corus')
   
% - Auto-mac added to menu
% - Fmin cross use the java value of idopt
% - Fmax cross use the java value of idopt
% - Recompute
% - GfSvd
% - buttons

iicom('curveloadGartid');
clear('sdtroot');UI=sdtroot('PARAMUI');UI.DefBut.Id
 UI.DefBut.Id.mtable={'Id.ObjectInfo','','',''
  'Id.w', ...
  {'Id.feplot';'Id.vf20';'Id.svd';'Id.stab';'Id.eup';'Id.eopt';'Id.est'} ...
  {            'Id.vf20';'Id.vsvd';'Id.vstab';'Id.veup';'Id.veopt';'Id.vest'} ''
  'Id.w0' 'Id.wmo','',''
  'Id.iipo1',{'Id.er';'Id.vf20';'Id.ea'},'Id.iipo',''
  'Id.e','','Id.ve',''
  };
UI.DefBut.Id.w.sdim='[2 1 165 -1 5  120 -1 5 1 1]';
UI.DefBut.Id.stab.callback={'diag_stab'};

sdtroot('setUI',UI);
delete(ancestor(clean_get_uf('iiplotpro',2),'figure'))
iicom('pro');idcom;iicom('curtab Ident');

%% 
elseif comstr(Cam,'test');[CAM,Cam]=comstr(CAM,5); 

if carg<=nargin;RO=varargin{carg};carg=carg+1;else;RO=struct;end

if comstr(Cam,'one')  
%% #TestOneDof

f=linspace(0,20,1024)';
XF=struct('w',f,'xf',qbode(1,[1 2*.01*10*2*pi (10*2*pi)^2],f*2*pi), ...
    'dof',[1.01 1.01]);
XF.idopt=idopt;XF.idopt(1:6)=[3 0 11 1 length(f) 1];
out1=struct('po',[10 .01],'res',1);

elseif comstr(Cam,'sdt_id') 
%% #TestSdt_Id Second example with 4 poles

iicom curveload sdt_id
r1=load('sdt_id.mat');
ci=iiplot;  %idcom('polyroot 10 10')
XF=ci.Stack{'Test'};
if isempty(strfind(Cam,'all'));ind=2;% example with a single channel
 XF.xf=XF.xf(:,ind); 
 XF.dof=XF.dof(ind,:); % example with a single channel
end
out1=struct('po',r1.def.data,'res',[]);out1.po(:,2)=.01;

else; error('Test%s unknown');
end

if isfield(RO,'GxxRef') % Generate cross correlation to ref
 XF.xf=XF.xf.*conj(repmat(XF.xf(:,RO.GxxRef),1,size(XF.xf,2)));
 % xxx idopt
end
out=XF;

%% -----------------------------------------------------------------------
elseif comstr(Cam,'polyref')
%-- Poly reference Least Square Complex Frequency --%
%-- Common denominator model --%
%-- discrete time domain / complex coefficients --%
%-- from P. Verboven and B. Cauberghe -- PhD theses --%
%--  http://mech.vub.ac.be/avrg/phd.htm
%-- po=diag_stab(XF,max_order);
%--

XF=varargin{carg};carg=carg+1;
if carg<=nargin;RunOpt.n=varargin{carg};carg=carg+1;
else;RunOpt.n=10;
end

% select a range (see iiplot)
IIw=XF(1).w;IIxf=XF(1).xf;
ind_w=XF(1).idopt(4):XF(1).idopt(5);
ind_xf=1:size(IIxf,2); xf=IIxf(ind_w,ind_xf);

%-- reduce data with svd --%

[u,s,v]=svd(xf,0);
s=diag(s);
sr=0;
ds=diff(log(s));
for i1=2:length(ds);
    sr(i1)=sr(i1-1)+ds(i1).^2;
end;
sr=sr/sum(ds.^2);
ind_s=max([find(sr<=0.9) RunOpt.n/2]);ind_s=1:ind_s;
if max(ind_s)>size(IIxf,2)
    ind_s=1:size(IIxf,2);
end;
xf=u(:,ind_s)*diag(s(ind_s))*v(ind_s,ind_s)';

%--

[nw,nxf]=size(xf); %xf=xf.*(w(:,ones(nxf,1)).^(-2));
xf=[xf ; conj(xf(end:-1:1,:))];


w=1i*IIw(ind_w);
Z_trans_coeff=1;
z=exp(1i*2*pi*[1:2*length(w)]'/(2*length(w))*Z_trans_coeff);

%-- coefficient for discrete time to continus time conversion --%

AA=[w(1) 1 ; w(end) 1];
tt=AA\[log(z(1)) ; log(z(end/2))];
gamma=tt(1);  delta=tt(2);


Poly=zeros(length(z),RunOpt.n+1);
for j1=0:RunOpt.n; Poly(:,j1+1)=z.^j1; end


R=Poly'*Poly;
%[u,s,v]=svd(R); possibly pseudo inverse
Rm1=inv(R);


Ups=[];
S=[];
T=[];
M=zeros(RunOpt.n+1);
%h = waitbar(0,'Computation in progress...');
for i1=1:nxf
   % waitbar(i1/nxf,h);
%-- Store data if residues are to be computed --%  
%   Ups=xf(:,i1*ones(nd+1,1)).*Poly;
%   S=[S ; Poly(:,1:nn+1)'*Ups];
%   T=[T ; Ups'*Ups];
%   M=M + T((nd+1)*(i1-1)+(1:nd+1),:)-S((nn+1)*(i1-1)+(1:nn+1),:)'*Rm1*S((nn+1)*(i1-1)+(1:nn+1),:);

  Ups=xf(:,i1*ones(RunOpt.n+1,1)).*Poly;
  S=Poly(:,1:RunOpt.n+1)'*Ups;
  T=Ups'*Ups;
  M=M + T-S'*Rm1*S;
end;
%close(h);

% Now loop on pole estimation
out=struct('po',[],'res',[]);

for j1=2:RunOpt.n
  % build polynomial and find roots
  a=[1 ; -M(2:j1+1,2:j1+1)\M(2:j1+1,1) ];
  uu=roots(flipud(a));  
  % convert poles from discrete to continus time domain 
  po=(log(uu)-delta.*sign(imag(uu)))/gamma;
  
  po(imag(po)<0)=[]; % keep positive frequencies
  po=[abs(po) -real(po)./abs(po) ones(size(po,1),1)*j1];
  out.po=[out.po;po];   
end

%% -------------------------------------------------------------------------
% stabilisation diagram for a frequency domain polynomial
elseif comstr(Cam,'poly')

XF=varargin{carg};carg=carg+1;
if carg<=nargin;RunOpt=varargin{carg};carg=carg+1;
    if ~isstruct(RunOpt);RunOpt=struct('ind',RunOpt);end
else;RunOpt.ind=3:15;
end

out=struct('po',[],'res',[],'dof',XF.dof,'idopt',XF.idopt);

for j2=RunOpt.ind(:)' % compute the values for a variable order
  [num,den]=id_poly(XF.xf,XF.w,j2*2,j2*2,XF.idopt);
  po=ii_pof(roots(den)/2/pi,3); 
  [res,po,xe]=id_rcopt(XF.xf,po,XF.w,XF.idopt);
  po(:,3)=j2;out.po=[out.po;po];out.res=[out.res;res];
end

% -------------------------------------------------------------------------
elseif comstr(Cam,'stab')

ID=varargin{carg};carg=carg+1;

% build the stabilization diagram
RO.level=unique(ID.po(:,3));
RO.ftol=.01;
RO.dtol=.1;
RO.rtol=.1;

for j1=2:length(RO.level)
  if j1==2; iold=find(ID.po(:,3)==RO.level(1));
  else; iold=ind;
  end
  ind=find(ID.po(:,3)==RO.level(j1));
  pold=ID.po(iold,:); 
  if size(ID.res,1)>iold(end);resold=ID.res(iold,:);else;resold=[];end
  
  for j2=1:length(ind) % Poles of current level
   po=ID.po(ind(j2),:);
   i1=find(abs(po(1)-pold(:,1))<po(1)*RO.ftol);
   po(4)=0; if ~isempty(i1); po(4)=1;end % frequency stabilized
   if po(4)
    [i2,i3]=min(abs(po(2)-pold(i1,2)));i1=i1(i3);
    if abs(po(2)-pold(i1,2))<po(2)*RO.dtol; po(4)=2; end% damping stabilized
    if isempty(resold)
    elseif norm(ID.res(ind(j2),:)-resold(i1,:))< ...
            norm(ID.res(ind(j2),:))*RO.rtol;
        po(4)=3; % residue stabilized
    end
    if po(2)<0;po(4)=-1; end % unstable pole
   end
   ID.po(ind(j2),1:4)=po;
  end
end

figure(1); clf;
RO.seq='xo+*';
RO.cseq=[0 0 1;0 .5 0;1 .5 .5;1 0 0];
RO.wseq=[.5 1 2 2];
for j1=0:3
  ind=ismember(ID.po(:,4),j1);
  if ~any(ind); continue; end
  h(j1+1)=line(ID.po(ind,1),ID.po(ind,3),'linestyle','none', ...
      'marker',RO.seq(j1+1),'color',RO.cseq(j1+1,:), ...
      'linewidth',RO.wseq(j1+1));
end
xlabel('Frequency [Hz]');ylabel('Order');
st={'none','freq','damp','res'};
legend(h,st{1:length(h)},'location','EastOutside')

%plot(IIw,db(ii_mmif(IIxf,IDopt,'sum')));

% ------------------------------------------------------------
elseif comstr(varargin{1},'cvs')
 out=sdtcheck('revision');
  %out='$Revision: 490 $  $Date: 2020-02-26 09:49:43 +0100 (Wed, 26 Feb 2020) $'; return;
% ------------------------------------------------------------
else; error('%s unknown',CAM);
end

