function [out,out1,out2]=fe_homo(varargin)

% FE_HOMO: finite element utilities for homogeneization and periodic 
%   computations
%
% EXPERIMENTAL utilities that are not included in SDT. 
% SDTools thus does not guarantee that support will be provided.
%
% See <a href="matlab: sdtweb _taglist fe_homo">TagList</a>


%       E. Balmes
%       Copyright (c) 1990-2025 by SDTools, All Rights Reserved.
%       For revision information use fe_homo('cvs')

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*NOSEM>
if nargin<1; CAM=''; Cam='';
elseif ~ischar(varargin{1}); 
  obj=varargin{1};evt=varargin{2};[CAM,Cam]=comstr(varargin{3},1); carg=24;
else; [CAM,Cam]=comstr(varargin{1},1); carg=2;
end

if comstr(Cam,'@');out=eval(CAM);
%% #DFT : periodic computations assuming propagating waves
% theory currently in rotor.tex and tex/dft.tex
% Load fixed in space
%  DofLoad sec(i), curve(t)
%   def.k_spectrum xxx Wave discrete
%   def.om_spectrum 
%   Generate a range grid with
%    'lcx','lcy','freq'
%    for cyclic symmetry 'dia'
%    sdtweb t_cyclic('basicHarmonic')
elseif comstr(Cam,'dft'); [CAM,Cam]=comstr(CAM,4);
%% #DftBuild fe_homo('DftBuild',mdl,dir) 
if comstr(Cam,'build');[CAM,Cam]=comstr(CAM,6);
 mdl0=varargin{carg};carg=carg+1;
 dir=varargin{carg};carg=carg+1;
 r0=[];
 if ~isempty(dir)&&size(dir,1)~=3
  if dir(1)==-1;dir(1)=max(mdl0.Node(:,5))-min(mdl0.Node(:,5));end
  if length(dir)>1&&dir(2)==-1;dir(2)=max(mdl0.Node(:,6))-min(mdl0.Node(:,6));end
  if length(dir)>2&&dir(3)==-1;dir(3)=max(mdl0.Node(:,7))-min(mdl0.Node(:,7));end
  dir=diag(dir);dir(end+1:3,1)=0;
 end
 if ~isempty(dir)&&any(dir(:,1))
  mdl0=fe_cyclic(sprintf('build -1 %.15g %.15g %.15g %s',dir(:,1),CAM),mdl0);
  r0=fe_case('GetDataSymmetry',mdl0);
 end
 if size(dir,2)>1&&any(dir(:,2))
  % builds periodicity /y :
  if isempty(strfind(Cam,'epsl'));CAM=[ ' epsl .001' CAM];end
  mdl0=fe_cyclic(sprintf('build -1 %.15g %.15g %.15g %s', ...
      dir(:,2),CAM),mdl0);
  r1=fe_case('GetDataSymmetry',mdl0);
  if isempty(r0);r0=r1;else; r0.IntYNodes=r1.IntNodes;end
 end
 if size(dir,2)>2&&any(dir(:,3))
 % builds periodicity /z :
  mdl0=fe_cyclic(sprintf('build -1 %.15g %.15g %.15g %s',dir(:,3),CAM),mdl0);
  r1=fe_case(mdl0,'getdata','Symmetry');
  r0.IntZNodes=r1.IntNodes;
 end
 r0.CellDir=dir;r0=feutil('rmfield',r0,'trans');
 mdl0=fe_case(mdl0,'stack_set','cyclic','Symmetry',r0);
 if isfield(mdl0,'Case');mdl0.Case=stack_set(mdl0.Case,'cyclic','Symmetry',r0);end
 out=mdl0;
 
elseif comstr(Cam,'ref') 
%% #DFTRef Verify that a wave propagating at skewed angle gives a single harmonic
% should be moved to t_ function
xk=0:49; yk=0:99;
p=repmat(xk(:),1,length(yk));
p(:,:,2)=repmat(yk(:)',length(xk),1);

dir=2*pi.*[1/10 1/25];
u=sin(dir(1)*p(:,:,1)+dir(2)*p(:,:,2));

U=fft2(u);U(abs(U)<norm(U,'inf')*1e-10)=0;sparse(U)

%figure(1);mesh(abs(U))
figure(1);clf;
subplot(211);pcolor(abs(U));%ii_plp('spy',abs(U))
subplot(212);pcolor(u)

% Write boundary conditions on cell edges
% Equivalent procedures 
% compute Z q = F  for variable frequency (E,G(omega)) and wavelength
% Energy for given volume, find 

elseif comstr(Cam,'eig')||comstr(Cam,'ceig') % 
%% #DftEig : compute frequency response with harmonic
% fe_homo('dftEig',RO,model) RO.nc : cell periods, RO.FixReal
%  'info','EigOpt' should be set externally

 RO=varargin{carg};carg=carg+1;
 if ~isfield(RO,'UseLong');RO.UseLong=0;end
 model=varargin{carg};carg=carg+1;
 out1=[];
 if ~isfield(model,'K') % Force assembly
  model.K={};model=fe_case('assemble NoT -matdes 2 1 -SE',model);
  model=fe_def('zcoef-default',model);
 end
 if comstr(Cam,'ceig');
 elseif isfield(RO,'kCoef')
  model.K={model.K{1} feutilb('sumkcoef',model.K,RO.kCoef)};
 end
 if isfield(model,'Case');Case=model.Case;
 else;Case=fe_case(model,'getcase');
 end
 RO=dftu('clGen',RO,model);
 if isfield(RO,'EigOpt');eigopt=RO.EigOpt;
 elseif isfield(RO,'list'); eigopt=[]; % No need for eigopt in fe2xf Build
 else; eigopt=fe_def('defeigopt',model);
 end
 if ischar(eigopt);eigopt=eval(eigopt);end
 if isfield(RO,'Load')  
  %% #DftEigReduce -2, called from fe_homo('dftdisp') - - - - - - - - - - - - 
  RO.b=RO.Load.def;RO.ak=1;RO.u=1;% 1i];
  [K2,T,RO]=dftu('buildK',RO,model,Case);%sdtweb fe_homo buildk
  K2(1:2)={(K2{1}+K2{1}')/2,(K2{2}+K2{2}')/2};
  if ~isreal(K2{2}); dbstack; keyboard;
      K2=[K2{1} real(K2{2}) imag(K2{2})];
   if length(K2)>2;SE.Opt(2,3)=1;SE.Klab{3}='Ki';end
  end
  SE=struct('Node',[1 0 0 0  0 0 0],'Elt',[], ...
      'K',{K2},'Klab',{model.Klab},'Opt',model.Opt,'DOF', ...
      (1:size(K2{1},1))'+.01,'Stack',{ ...
       {'info','EigOpt',eigopt;'info','Range',1;'info','Freq',1}});
  %xxx if isfield(model,'Opt');SE.Opt=model.Opt;SE.Klab=model.Klab;end
  if size(T,1)==2*size(RO.b,1)
    dbstack; keyboard
  else
    RO.Load=struct('DOF',SE.DOF,'def',T'*RO.b);
  end
  CE=struct('T',speye(length(SE.DOF)),'DOF',SE.DOF,'Stack',{{ ...
      'DofLoad','IN',RO.Load}});
  SE=stack_set(SE,'case','Case 1',CE);
  if size(SE.DOF,1)<20||(~isempty(eigopt)&&eigopt(1)==2)% Keep all
   d1=fe_eig({SE.K{SE.Opt(2,:)==2},SE.K{SE.Opt(2,:)==1},SE.DOF},2);   
   %struct('def',eye(size(SE.br,1)),'DOF',SE.DOF,'data',(1:length(SE.DOF))');
   SE.K=feutilb('tkt',d1.def,SE.K);
   SE.br=d1.def'*RO.Load.def;SE.lab_in=fe_c((1:size(SE.br,2))'+.99);
   zCoef=stack_get(model,'info','zCoef');
   if ~isempty(zCoef);SE=stack_set(SE,zCoef);end
  elseif ~isfield(RO,'list'); 
    SE=fe_reduc(['free-SE -matdes' sprintf('%i ',SE.Opt(2,:))],SE);
    zCoef=stack_get(model,'info','zCoef');
    if ~isempty(zCoef);SE=stack_set(SE,zCoef);end
    SE.br=SE.TR.def'*RO.Load.def;SE.lab_in=fe_c((1:size(SE.br,2))'+.99);
    d1=SE.TR;
    %SE=fe_reduc('free-SE',SE); d1=SE.TR;
  else; SE=fe2xf('Build;',SE,RO); d1=SE.TR;
    zCoef=stack_get(model,'info','zCoef');
    if ~isempty(zCoef);SE=stack_set(SE,zCoef);end
    SE.br=SE.TR.def'*RO.Load.def;SE.lab_in=fe_c((1:size(SE.br,2))'+.99);
  end
  out1=SE;
 elseif comstr(Cam,'ceig')
  %% Complex dispersion with no load - - - - - - - - - - - - - - - - - - - - - - 
  [K2,T,RO]=dftu('buildK',RO,model,Case);
  K3=cell(1,4); RO.UseLong=1;
  i1=model.Opt(2,:)==1;K3{1}=(K2{i1}+K2{i1}')/2;
  i1=model.Opt(2,:)==2;K3{3}=(K2{i1}+K2{i1}')/2;
  i1=model.Opt(2,:)==4;if any(i1);K3{3}=K3{3}+1i*K2{i1};end
  i1=model.Opt(2,:)==3;if any(i1);K3{2}=K2{i1};end
  K2=K3; clear K3; K2{4}=(1:size(K2{1},1))'+.01;
  d1=fe_ceig(K2,eigopt);
 else
  %% Real Dispersion with no load - - - - - - - - - - - - - - - - - - - - - - 
  if length(model.K)>2;error('Should have had RO.kCoef');end
  [K2,T,RO]=dftu('buildK',RO,model,Case);
  K2={(K2{1}+K2{1}')/2,(K2{2}+K2{2}')/2};
  d1=fe_eig({K2{1},K2{2},(1:size(K2{1},1))'+.01},eigopt);
 end
 N=size(model.K{1},1);md1=full(T*d1.def); d1.DOF=model.DOF; 
 if ~RO.isReal
  if RO.UseLong; d1.def=md1;d1.DOF=[d1.DOF;d1.DOF+.5];
  else; md1=complex(md1(1:N,:),md1(N+1:2*N,:));d1.def=md1;
  end
 elseif RO.UseLong;% Force long for real
     md1(N*2,1)=0; d1.def=md1; d1.DOF=[d1.DOF;d1.DOF+.5];
 else; d1.def=md1;
 end
 % d2=fe_def('subdofind',d1,1:N);d2.TR=model.TR;cf=feplot;cf.def=d2;
 data=stack_get(Case,'cyclic','Symmetry','get');
 d1.data(:,2:length(RO.dataOpt))= ...
     repmat(RO.dataOpt(2:end),size(d1.data,1),1);
 d1.Xlab={'DOF',{'Freq';'nx';'ny';'nz'}};
 d1.Xlab{2}(length(RO.dataOpt)+1:end)=[];
 st=sprintf('%.1f,',RO.dataOpt(2:end));
 d1.label=sprintf('nc(%s)',st(1:end-1));
 if isfield(RO,'jpar') % Allow parametric study
     d1.data(:,2)=RO.jpar; d1.data(:,3:end)=[];d1.Xlab{2}={'Freq','jPar'};
 end
 out=d1; if ~isempty(out1);out1.TR=out;end
 
elseif comstr(Cam,'disp') 
%% #DFTDisp dispersion curve computation - - - - - - - - - - - - 

if 1==2
 %% test, see also t_cyclic('DftIfft')
 model=comp12('rve1fiber',r1);cf=feplot(model);
 model.pl(:,5)=1.2e-9; % nominal density very high freq
 model=stack_set(model,'info','EigOpt',[5 10 -1e9]);
 % RO=struct('nc',[0 1e3 0],'FixReal',0,'EigOpt','vf',r1.vf);
 Range=struct('val',logspace(log10(2),3,20)','lab',{{'kcx'}});
 def=fe_homo('dftDisp',model,Range);
 fh=reshape(def.data(:,1),[],size(def.Range.val,1));
 ph=reshape(def.data(:,2),[],size(def.Range.val,1)); 
 figure(1);loglog(1./max(def.Range.val(ph(1,:),:),[],2),fh')
end

model=varargin{carg};carg=carg+1;RO=struct;
if carg<=nargin; Range=varargin{carg};carg=carg+1;
  if isfield(Range,'Range');RO=Range;Range=RO.Range;end
else  % Default wave in x direction
   Range=struct('val',logspace(log10(2.01),2,20)'*[1 0 0], ...
       'lab',{{'ncx','ncy','ncz'}});
   data=fe_case(model,'getdata','Symmetry');
   if ~isfield(data,'IntZNodes');Range.val(:,3)=[];Range.lab(3)=[];end
   if ~isfield(data,'IntYNodes');Range.val(:,2)=[];Range.lab(2)=[];end
end
    
if carg<=nargin&&isstruct(varargin{carg}); RO=varargin{carg};carg=carg+1;end
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   'UseLong(#3#"Use long storage") ' ...
   'FRF(0#3#"Compute FRF") ' ...
   ],{RO,CAM});Cam=lower(CAM); %#ok<NBRAK>

def=[];
if isfield(model,'TR')&&isfield(model.TR,'adof')
    RO.MatDes=model.Opt(2,:);
    SE=model;CE=fe_case(SE,'getcase');
    CE.T=speye(length(SE.DOF));
    if RO.FRF
     SE=fe_sens('br&lab',SE,'');
     Load=struct('def',SE.br,'DOF',SE.DOF,'lab',{SE.lab_in});
     SE=fe_sens('br&lab',SE,'');  
     RO.Freq=stack_get(model,'info','Freq','get');RO.Freq=RO.Freq(:);
     if ~isempty(fe_case(SE,'getdata','SensDof'))
      SE=fe_sens('cr&lab',SE,''); 
     end
    end
elseif RO.FRF
 if ~isfield(RO,'MatDes'); RO.MatDes=[2 1 3 4];end % Vary uses 5 (done top level)
 st=sprintf('%i ',RO.MatDes);
 [SE,CE,Load,Sens]=fe_case(sprintf('assemble -matdes %s -NoT -SE -load -sens', ...
     st),model);
 RO.Freq=stack_get(model,'info','Freq','get');RO.Freq=RO.Freq(:);
 SE=fe_def('zcoef-default',SE);sdtkey('cvsnum>1.333;','fe_def')
else;
 if ~isfield(RO,'MatDes');RO.MatDes=[2 1];end 
 [SE,CE]=fe_case(sprintf('assemble -matdes %s -NoT -SE', ...
     sprintf(' %i',RO.MatDes)),model);
end
SE.Case=CE;
RO.cyc=stack_get(CE,'cyclic','Symmetry','get');
RO=dftu('RangeNc',RO,Range); out2=[];
if isfield(model,'TR')&&RO.FRF==0; 
    RO.RangeNc.kCoef=double(RO.MatDes==1|RO.MatDes==5);
end
RO.nM=0;
RO.follow=zeros(1,2);
if 1==2
cingui('timerStart', ...
   struct('caller','fe_homo', ...
   'DispFcn','if ishandle(uo.h);waitbar(uo.data,uo.h);end', ...
   'h',waitbar(0,'PerDfrf'),'StartDelay',5, ...
   'StopFcn','if ishandle(uo.h);delete(uo.h);end'),RO.follow);
end
    
if RO.FRF 
 %% #DftDispFrf Case where shapes are computed too -3
 ofact silent
 if isempty(stack_get(SE,'','zCoef'));SE=fe_def('zcoef-default',SE);end
 for jpar=1:size(RO.RangeNc.val,1)
  R1=RO.RangeNc; R1.val=R1.val(jpar,:);
  sp_util('setinput',RO.follow,jpar/size(RO.RangeNc.val,1),zeros(1));
  RO.EigOpt=stack_get(SE,'info','EigOpt','get');
  if nnz(Load.def)==0;error('Zero load');end
  R1.Load=Load; 
  if isfield(RO,'BuildList'); R1.list=RO.BuildList;end
  if ~isfield(RO,'MifDes');RO.MifDes={'frf',[]};end
  R1.jpar=jpar;R1.UseLong=1;
  [d1,MVR]=fe_homo('dftEig',R1,SE);d1.label=['mode ',d1.label];% d1 expected in long format
  
  MVR=stack_set(MVR,{'info','Freq',RO.Freq;'info','MifDes',{RO.MifDes}});
  MVR.cr=eye(size(MVR.br,1));MVR.lab_out=fe_c((1:size(MVR.TR.def,2))'+.01);
  R2=fe2xf('frfzr',MVR);
  if isfield(RO,'cpxmode');% Store complex mode (not very robust)
      d2=fe_ceig(MVR);i4=size(d2.def,2)-size(d1.data,1)+1:size(d2.def,2);
      i4(i4<1)=[];d1.data(1:length(i4),[1 3])=d2.data(i4,:);
      d1.def=MVR.TR.def*d2.def(:,i4);
  end
  d2=reshape(permute(R2.Y,[2 1 3]),size(R2.Y,2),[]);
  if isfield(RO,'Sens')
   Sens=fe_case(model,'Sens');Sens=feutilb('placeindof',MVR.TR.DOF,Sens);
   d2=struct('def',(Sens.cta*MVR.TR.def)*d2, ...
      'DOF',Sens.tdof,'data',R2.X{1}, ...
      'Xlab',{{'DOF',{'Freq';'jPar'}}},'name','RESP');d2.data(:,2)=jpar;
   
  else;
   d2=struct('def',MVR.TR.def*d2, ...
      'DOF',MVR.TR.DOF,'data',R2.X{1}, ...
      'Xlab',{{'DOF',{'Freq';'jPar'}}},'name','RESP');d2.data(:,2)=jpar;
  end 
  if size(R2.Y,3)~=1;d2.data=repmat(d2.data,size(R2.Y,3),1);end% Multi input
  if isfield(RO,'SubDof'); % Allow extraction but keep order
   if ~isfield(RO,'iSubDof');RO.iSubDof=sort(fe_c(d1.DOF,RO.SubDof,'ind'));end
   d1=fe_def('subdofind',d1,RO.iSubDof);
   d2=fe_def('subdofind',d2,RO.iSubDof);
  end
  if isempty(out2); out2=d2; 
      if isfield(RO,'Range')&&isfield(RO.Range,'val')
       out2.Range=RO.Range;out2.Range.CellDir=RO.cyc.CellDir;
      else;% Just a loop on periodicity
       out2.Range=RO.RangeNc;out2.Range.CellDir=RO.cyc.CellDir;
      end
      RO.i2de=[0 numel(out2.def)];RO.i2da=[0 numel(out2.data)];
      out2.def(1,size(RO.RangeNc.val,1)*size(out2.def,2))=0;
      out2.data(size(RO.RangeNc.val,1)*size(out2.data,1),1)=0;out2.data=out2.data';
  else; % out2=fe_def('appenddef',out2,d2);
      sp_util('setinput',out2.def,d2.def,RO.i2de);
      sp_util('setinput',out2.data,d2.data',RO.i2da);
  end
  if isempty(def) % Deal with reduced basis output
  elseif ~isreal(d1.def)&&R1.val(strcmpi(R1.lab,'ncx'))~=1&&size(d1.def,2)>=2*RO.nM
      d1=fe_def('subdef',d1,1:2:2*RO.nM); % One vect
  elseif size(d1.def,2)<RO.nM;d1.def(1,RO.nM)=0;d1.data(end+1:RO.nM,1)=NaN;
  else; d1=fe_def('subdef',d1,1:RO.nM); 
  end 
  if ~RO.nM; RO.nM=size(d1.def,2);end
  if isempty(def);def=d1;
  else; def=fe_def('appendDef',def,d1);
  end
 end
%% Disp but not FRF
else
 if isfield(RO.RangeNc,'EigFcn');RO.EigFcn=RO.RangeNc.EigFcn;
 elseif ~isfield(RO,'EigFcn');RO.EigFcn={'fe_homo','dftEig'};
 end
 for jpar=1:size(RO.RangeNc.val,1)
  sp_util('setinput',RO.follow,jpar/size(RO.RangeNc.val,1),zeros(1));
  R1=RO.RangeNc; R1.val=R1.val(jpar,:);
  if isfield(RO,'UseLong');R1.UseLong=RO.UseLong;end
  d1=feval(RO.EigFcn{:},R1,SE);%fe_homo('dftEig')
  if isfield(R1,'ModalFilter')&&R1.ModalFilter==R1.val(1)
    m=SE.K{1};z=spalloc(size(m,1),size(m,2),1);k=SE.K{2};
    m=[m z;z m];k=[k z;z k];%feutilb('dtkt',[real(d1.def);imag(d1.def)],{m,k})    
    d1.Filter=d1;
    if rem(d1.DOF(end),1)<.5;r1=[real(d1.def);imag(d1.def)]; else;r1=d1.def;end
    d1.Filter.pdef=pinv(r1'*m*r1)*(r1'*m); d1.Filter.def=r1; 
  end
  d1.data(:,2)=jpar; d1.data(:,3:end)=[];d1.Xlab{2}={'Freq','jPar'};
  if R1.val(1)~=1 % Not ncx==1
   d1=fe_def('subdef',d1,1:2:size(d1.def,2)); % On vector of the pair
  end
  if size(d1.def,2)<RO.nM; 
    d1.def(:,end+1:RO.nM)=NaN;d1.data(end+1:RO.nM,1)=NaN;
  elseif ~RO.nM&&~isempty(def); 
      RO.nM=size(d1.def,2);def=fe_def('subdef',def,1:RO.nM);
  end
  def=fe_def('appendDef',def,d1);
 end
end

try
 i1=reshape(1:size(def.def,2),[],size(RO.RangeNc.val,1))';
 def=fe_def('subdef',def,reshape(i1,[],1));
 ph=reshape(def.data(:,2),size(RO.RangeNc.val,1),[]); 
 RO.NeedHist=1;
catch
 sdtw('_nb','Failed hist reshape');
 RO.NeedHist=0;
end
if ~isempty(out2);out2.data=out2.data';end

% Clean up output
def.Range=RO.RangeNc;  
if isfield(RO.cyc,'trans');def.Range.CellDir=RO.cyc.trans(:);end
def.LabFcn='fe_range(''LabDef'',def,ch);';
def.Range.param=struct;
def.Range.param.Freq.LabFcn='sprintf(''%.2f Hz'',def.data(ch,1));';

%if any(strcmp(def.Range.lab,'inx'))
% def.Range.param.icx=struct('LabFcn', ...
%     'sprintf('' %i@kx=%.2f m^{-1}'',nnz(def.data(1:ch,2)==def.data(ch,2)),1/val)');
%end
i1=find(strcmp(def.Range.lab,'iny'));
if ~isempty(i1)&&~any(def.Range.val(:,i1))
 def.Range.param.iny=struct('LabFcn','''''');
end
def=dftu('addrange',def,'kx');
if any(strcmpi(def.Range.lab,'ncx'))
def.Range.param.ncx=struct('LabFcn', ...
     'sprintf('' %i@lx=%.2f m'',nnz(def.data(1:ch,2)==def.data(ch,2)),val)');
end

out=def;
if RO.NeedHist
 data=fe_case(model,'getdata','Symmetry');
 i1=strcmpi(def.Range.lab,'kx');
 hist=struct('X',{{ def.Range.val(:,i1) ... %*diag(data.CellDir(1:length(RO.i_nc))))...
     []}}, ... % Kappa=1/(nc)
     'Xlab',{{'kx','Mode'}}, ...
     'Y',reshape(def.data(:,1),size(def.Range.val,1),[]));
 if isfield(SE,'dof_in') % OutIn
   cr=fix(SE.dof_in)+(1:3)'/100;
   if isfield(SE,'TR');
       st=fe_c(SE.TR.DOF,cr,'dofs');cr=fe_c(SE.TR.DOF,cr)*SE.TR.def;
   else; st=fe_c(SE.DOF,cr,'dofs');cr=fe_c(SE.DOF,cr);
   end
   hist.X{3}=[{'Freq'};st];hist.X{2}=(1:size(hist.Y,2))';
   hist.Xlab{3}={'Freq/DOF'}; 
   cr=reshape(cr*def.def(1:size(cr,2),:),length(st),size(def.Range.val,1),[]);
   hist.Y(:,:,1+(1:length(st)))=permute(cr,[2 3 1]);
 end
 i2=all(hist.X{1}==0|~isfinite(hist.X{1})|hist.X{1}==1); % 0=1=Infinite period
 if any(i2); hist.X{1}(:,i2)=[];end
 i2=find(hist.X{1}==0); % Real
 if ~isempty(i2)&&size(hist.Y,3)==1; % duplicate frequencies at k=0
  r1=reshape(hist.Y([i2;i2],:),1,[]);
  hist.Y(i2,:)=r1(1:size(hist.Y,2));
 end
 if isfield(RO,'Range')&&isfield(RO.Range,'val');hist.Range=RO.Range;
 else;hist.Range=RO.RangeNc;end
 hist.Range.CellDir=RO.cyc.CellDir;
 hist.X{2}=(1:size(hist.Y,2))'; % Mode indices
 out1=hist;
else;out1=[];
end

elseif comstr(Cam,'dfrf') 
%% #DFTDfrf : direct frequency response at multiple k/eta

[mdl,RunOpt,SE,Case,Load,w,ft,carg]=safeInitDfrf(varargin,carg);

 if isfield(Case,'TIn') % support DOFSet entries - - - - - - - - - - - - 
     error('DofSet not supported for DFTDfrf')
 else % standard load input - - - - - - - - - - - - - - - - - - - - - - -
  na=size(Load.def,2);
  if isfield(RunOpt,'loadcombine');na=1;end
  if isempty(Load.def)||nnz(Load.def)==0; error('Load building failed');
  elseif isequal(Load.DOF,SE.DOF);b=Load.def;
  elseif isequal(Load.DOF,Case.DOF);b=Load.def;
  else;error('Load building failed');
  end
  if size(SE.K{1},1)~=size(b,1) 
   error('Dimension mismatch between Load and matrix');
  end
  if ~isempty(Load.def) && isfield(Case,'TIn')
   sdtw('_nb','when load and DofSet, load is ignored'); end
  def1=struct('def',zeros(size(SE.DOF,1),length(w)*na), ...
   'DOF',SE.DOF,'data',w/2/pi,'Xlab',{{'DOF',{'Freq';'IndDef';'jPar'}}});
  if isfield(Load,'lab');  def1.lab_in=Load.lab;
  else; def1.lab_in=cellfun(@(x) sprintf('in_%i',x),...
          num2cell((1:na)'),'UniformOutput',0);
  end
  w=w(:);
  if isfield(Load,'curve')&&~isempty(Load.curve)
  end
  [zCoef,K]=feval(fe_simul('@safe_dfrf_zcoef'),w,SE,SE);
  if size(Case.T,2)<1e3; ofact('silent');end
  RunOpt=dftu('clGen',RunOpt,mdl); 
  RunOpt.icol=0;if ~isfield(RunOpt,'Post');RunOpt.Post='';end
  for j1=1:length(w)
  model=SE;model.K={feutilb('sumkcoef',K,zCoef(j1,:))};
  %if j1==length(w);clear K;end
  for jk=1:size(RunOpt.RangeNc.val,1) % may not be smart to loop on k internally
   RunOpt.nc=RunOpt.RangeNc.val(jk,1:end-1);
   if isfield(RunOpt,'Range')&&isfield(RunOpt.Range,'lab')&&isequal(RunOpt.Range.lab,{'ncx','Freq'})
    if ~ismember(round([RunOpt.nc w(j1)/2/pi]*1000),round(RunOpt.Range.val*1000),'rows')
        continue;
    end
    fprintf('Nc=%i w=%.1fHz\n',[RunOpt.nc w(j1)/2/pi]);
   end  
   if isempty(ft);u=eye(na);
   elseif isfield(RunOpt,'loadcombine');u=ft(j1,:).'; % Use combined loads
   else;u=diag(ft(j1,:));
   end
   RunOpt.u=u;
   RunOpt.ak=RunOpt.RangeNc.val(jk,strcmpi(RunOpt.RangeNc.lab,'ak'));
   if isempty(RunOpt.ak);error('Missing ak in range definition');end
   RunOpt.b=b; %b load on cell 1, decomposed on [r(bu) /im(bu)]
   %RunOpt.mpcond={(22:36)'/100,1e8};dbstack
   if size(RunOpt.RangeNc.val,1)==1
    [K2,T,RO]=dftu('buildK',RunOpt,'model',Case);
   else
    [K2,T,RO]=dftu('buildK',RunOpt,model,Case);% Can't be 'model' for multiple
   end
   if isfield(RO,'blab');def1.lab_in=RO.blab;end 
   if isfield(RO,'iter');error('Report EB should be oProp');
    kd=K2{1};u=ofact(kd,T'*RO.b,struct('iter',RO.iter,'iterOpt',{RO.iterOpt}));
    u=T*u;
   else
    %figure(1);r1=feutilb('sumkcoef',K,zCoef(j1,:));r1=mean(abs(r1));[r2,i2]=sortrows([round(rem(model.DOF(:,1),1)*100) fix(model.DOF)]);semilogy(r1(i2))   
    %z=model.DOF(r1<1e-7);z(fe_c(z,(1:3)'/100,'ind'))=[];fecom('shownodemark',z)
    if 1==2%exist('K','var')
     Z=Case.T'*feutilb('sumkcoef',K,zCoef(j1,:))*Case.T;p_pml('ViewLuBlocks',struct('K',Z,'DOF',Case.DOF));
    end
    if 1==2
     [r3,i1]=sort(round(rem(Case.DOF,1)*100));r1=r3;
     for j2=1:3;r1(:,j2)=abs(diag(Case.T'*K{j2}*Case.T))*mean(abs(zCoef(:,j2)));end
     figure(1);semilogy(r1(i1,:));ii_plp(find(diff(r3))*[1 0]);
     Z=full(K2{1}); [L,U,p]=lu(Z);figure(11);semilogy(abs(diag(U)))
     Z=feutilb('sumkcoef',K,zCoef(j1,:));p_pml('ViewLuBlocks',struct('K',Z,'DOF',RunOpt.DOF));
     u=T*(kd\(T'*RO.b)); ofact('clear',kd);% [real;imag]* cpx(omega)
    end
    %i3=find(any(T(:,~any(K2{1})),2));
    %if ~isempty(i3);r2=[model.DOF;model.DOF+.5];disp(fe_c(r2(i3)));end 
    
    kd=K2{1}; kd=ofact(kd,RunOpt.oProp{:});u=T*(kd\(T'*RO.b)); ofact('clear',kd);% [real;imag]* cpx(omega)
    if 1==2
     % diagnostic for PML
     z=struct('K',Case.T'*feutilb('sumkcoef',K,zCoef(j1,:))*Case.T,'DOF',Case.DOF);p_pml('viewlublocks',z)
     d1=struct('def',T(1:size(T,1)/2,:)*fliplr(a),'DOF',SE.DOF,'data',flipud(diag(b)));
     ic=1;fe_c(d1.DOF(any(abs(d1.def(:,ic))>.5*norm(d1.def(:,ic),'inf'),2)));
     cf.def=d1; fecom('colordata22')
    % [a,b]=svd(full(K2{1}));d1=struct('def',T(1:size(T,1)/2,:)*a,'DOF',SE.DOF,'data',diag(b),'TR',SE.TR)
      cf=feplot;def=struct('def',reshape(u,[],2),'DOF',SE.DOF);cf.def=def;
      cyc=fe_case(model,'getdata','Symmetry');    
     RO.ds=fe_c(def.DOF,(22:36)'/100,'dof');
     RO.de=fe_c(RO.ds,cyc.IntDof(:),'dof');
     RO.di=fe_c(RO.ds,RO.de,'dof',2);
     fe_c(def.DOF,RO.de)*def.def;sum(abs(ans))
     fe_c(def.DOF,RO.di)*def.def;sum(abs(ans))
     % check oscilation that happens within interior
     r1=full(fe_c(def.DOF,feutil('getdof',[36 116 121 189]',[22 25 26 28 30 31 34 35 36]'/100))*def.def);r1(abs(r1)<norm(r1,'inf')*.001)=0;r1=reshape(r1,9,[]);[real(r1);imag(r1)]
     [U,s]=svds(K2{1},5,'smallest');def=struct('def',reshape(T*U(:,1),[],2),'DOF',SE.DOF);
     
    end
   end
   %  cf=feplot;def=struct('def',reshape(u,[],2),'DOF',SE.DOF);cf.def=def;
   % fecom colordataevala
   RunOpt.icol=RunOpt.icol(end)+(1:size(u,2));na=size(u,2);
   def1.def(1:size(u,1),RunOpt.icol) =u;
   def1.data(RunOpt.icol,1)=w(j1)/2/pi;
   def1.data(RunOpt.icol,3)=jk;
   def1.data(RunOpt.icol,2)=(1:na)';
   eval(RunOpt.Post);
   % jar
  end
  end
  def1.LabFcn='sprintf(''%s (%g Hz, jk%i)'',def.lab_in{def.data(ch,2)},def.data(ch,[1 3]))';
  def1.Range=RunOpt.RangeNc;def1.Range.CellDir=RunOpt.CellDir;
 end % DofSet or not
 if max(def1.data(:,1))>1e4
  def1.Range.param.Freq.LabFcn='sprintf(''%.0f kHz'',def.data(ch,1)/1000);';
 else
  def1.Range.param.Freq.LabFcn='sprintf(''%.0f Hz'',def.data(ch,1));';
 end

 if size(def1.DOF,1)==size(def1.def,1)/2
   def1.DOF=[def1.DOF;def1.DOF+.5]; % Use long format with DOF shifted by 50
 end
 def1=dftu('addrange',def1,{'kx'});

 out=def1;
elseif comstr(Cam,'redp2') 
%% #DFTRedP2set : SE generation with left/right (phase2) 
 Case=[];SE=[]; def=[];
 eval(iigui({'SE','def','RO','Case','Load'},'MoveFromCaller'))

 % Inspired from sdtweb dyn_solve RedV2Phase2 
 cyc=stack_get(Case,'cyclic','Symmetry','get');
RO.EdgeDof=fe_c(def.DOF,cyc.IntNodes(:,1),'ind');
RO.EdgeDof(:,2)=fe_c(def.DOF,cyc.IntNodes(:,2),'ind');
if ~isfield(RO,'EdgeTol');RO.EdgeTol=1e-5;end
RO.DofIn=def.DOF;RO.DofIn(RO.EdgeDof(:))=[];

 % Now gradually build the basis
zCoef=stack_get(SE,'','zCoef','get');
if isempty(zCoef)
 k=SE.K(ismember(SE.Opt(2,:),[1 5]));if length(k)>1;error('Mismatch');end
 m=SE.K(ismember(SE.Opt(2,:),2));if length(m)>1;error('Mismatch');end
 k=k{1};m=m{1};
else; 
 [zCoef,SE]=feval(fe2xf('@get_zCoef'),SE); 
 m=feutilb('sumkcoef',SE.K,vertcat(zCoef{2:end,2}));
 k=feutilb('sumkcoef',SE.K,vertcat(zCoef{2:end,3}));
end
%RO.pml=fe_c(SE.DOF,(13:36)','ind');
% keep independent mecha shapes
RO.mech=find(diag(k)~=0);%fe_c(SE.DOF,unique(fix(SE.DOF(RO.pml))),'ind',2);
RO.pml=setdiff((1:length(SE.DOF))',RO.mech);
%i1=~ismember(RO.EdgeDof(:,1),RO.mech);
%RO.EdgePml=RO.EdgeDof(i1,:);RO.EdgeDof(i1,:)=[];

%cg=feplot(10);cg=feplot(SE,struct('def',[r1.Tl+r1.Tr],'DOF',SE.DOF))
%% Possibly LR shapes by DOF subsets
%if ~isfield(RO,'P2Sets');RO.P2Sets={SE.DOF,'fe_norm'};end
%if ~isfield(RO,'P2Sets');RO.P2Sets={SE.DOF,'lrisvd'};end
% fecom('shownodemark',RO.P2Sets(:,1),'marker','o')
if ~isfield(RO,'P2Sets'); % DftRedP2LU 
  RO.P2Sets={SE.DOF(RO.EdgeDof(:)),'svd lrilu static','edge';
             SE.DOF,struct('type','fe_norm','noEdge',1),'interior'
             };
elseif ischar(RO.P2Sets{1}); %{name,DOF,type} -> {DOF,type,name};
  RO.P2Sets=RO.P2Sets(:,[2 3 1]);
elseif size(RO.P2Sets,2)<3; 
 RO.P2Sets(:,3)=cellfun(@(x)sprintf('Set%i',x),num2cell(1:size(RO.P2Sets,1)),'uni',0);
end
for j1=1:size(RO.P2Sets,1) % Robust format RC
 if ~isstruct(RO.P2Sets{j1,2});
        RO.P2Sets{j1,2}=struct('type',RO.P2Sets{j1,2});
 end
 RC=struct('EdgeDof',[],'EdgeTol',RO.EdgeTol,'DOF',[],'k',[],'Active',0);
 RC=sdth.sfield('AddMissing',RO.P2Sets{j1,2},RC);
 RO.P2Sets{j1,2}=RC;
end
if ~isfield(RO,'fe_coor');RO.fe_coor='lrilu';end
if ~isfield(RO,'SvdTol');RO.SvdTol=1e-8;end
RO.usedIndDof=find(~any(def.def,2));kd=[];
if ~isreal(def.def);def.def=[real(def.def) imag(def.def)];end

for j1=1:size(RO.P2Sets,1)
 RC=RO.P2Sets{j1,2};i2=RO.P2Sets{j1,1};if ischar(i2);eval(i2);end
 i2=fe_c(SE.DOF,i2,'ind');i2=setdiff(i2,RO.usedIndDof);
 RO.usedIndDof=[RO.usedIndDof;i2];
 ci2=setdiff(1:size(def.def,1),i2);
 T3=sparse(i2,1:length(i2),1,size(def.def,1),length(i2));
 RO.curActive=0;
 if RC.Active==1
  i2(fe_c(SE.DOF(i2),Case.DOF,'ind',2))=[]; RO.curActive=1;% Only keep active DOFs
  T3=Case.T(:,fe_c(Case.DOF,SE.DOF(i2),'ind'));T2=def.def(i2,:);
  RC.k=T3'*k*T3;RC.DOF=SE.DOF(i2);
 else; 
   T2=def.def(i2,:);k2=k(i2,i2);RC.k=k(i2,i2);RC.DOF=SE.DOF(i2);
 end
 RO.curCoor=RO.fe_coor;RO.curSetType=RC.type;if isempty(T2); continue;end
 [T2,RC,RO]=genT2(T2,RC,RO,i2);

 % fecom('shownodemark',{RC.DOF(RC.EdgeDof(:,1)),RC.DOF(RC.EdgeDof(:,2)))})
 % d_dft('viewdebugT2')
 %% Build a basis that uses all given DOF and span the learning subspace
 if isfield(RC,'NoAdof')&&ischar(RC.NoAdof); 
    RC.NoAdof=fe_c(RC.DOF,feutil(['findnode' RC.NoAdof],SE),'ind');
 end
 if ~isnumeric(T2)
 elseif isempty(RC.EdgeDof)&&~isfield(RC,'noEdge');
   sdtw('_nb','P2Set %s no EdgeDof : if not an error use .noEdge=1',RO.P2Sets{j1,3});
   break;
  % else; T2=struct('Tl',[],'Tr',[],'Ti',T2);
  % end
 else;
  if 1==2
   curCoor=RO.curCoor;T2ref=def.def(i2,:);eval(iigui({'T2','RC','curCoor','T2ref'},'SetInBaseC'));
   T3=fe_coor(curCoor,T2,RC);subspace([T3.Tl,T3.Tr,T3.Ti],T2ref)
  end
  T0=T2;
  T2=fe_coor(RO.curCoor,T0,RC);% sdtweb fe_coor('lrisvd')
 end
 % d_dft('viewdebugT2') p_pml('viewqz'); 
 %  % d1=fe_eig({m,k,T3*[T2.Tl T2.Tr T2.Ti],SE.DOF},2);cf=feplot(10,';');cf.def=d1
 
 if isempty(T2.Tl);T2.Tl=zeros(size(T3,1),0);T2.Tr=zeros(size(T3,1),0);
 elseif ~isempty(strfind(RO.curSetType,'static'));
   %% Static restitution of existing
   T2.Tl=T3*T2.Tl; T2.Tr=T3*T2.Tr;
   if isempty(kd);
     RO.Stat.ic=fe_c(Case.DOF,SE.DOF(RO.EdgeDof(:)),'ind',2);
     RO.Stat.ia=fe_c(Case.DOF,SE.DOF(i2),'ind');
     RO.Stat.if=fe_c(SE.DOF,Case.DOF(RO.Stat.ia),'ind');
     RO.Stat.Tc=Case.T(:,RO.Stat.ic); kd=ofact(feutilb('tkt',RO.Stat.Tc,k));
   end 
   T2.Tl=T2.Tl+ RO.Stat.Tc* (kd\(-(RO.Stat.Tc'*k*T2.Tl)));
   T2.Tr=T2.Tr+ RO.Stat.Tc* (kd\(-(RO.Stat.Tc'*k*T2.Tr)));
   %Removed already used subspace
   def.def=def.def-T2.Tl*(T2.Tl(i2,:)\def.def(i2,:));
   def.def=def.def-T2.Tr*(T2.Tr(i2,:)\def.def(i2,:));
 elseif strcmpi(RO.curSetType,'CaseT')
 else
   T2.Tl=T3*T2.Tl; T2.Tr=T3*T2.Tr;
 end
 if strcmpi(RO.curSetType,'CaseT');
 elseif ~isempty(T2.Ti); T2.Ti=T3*T2.Ti;
 else; T2.Ti=zeros(size(T3,1),0);
 end
 % fecom('shownodemark',vertcat(T2.adof{:}),'marker','o','color','r')
 RO.P2info(j1,1:3)=[size(T2.Tl,2) size(T2.Tr,2) size(T2.Ti,2)];
 if j1==1;r1=T2;
 else
  r1.Tl=[r1.Tl T2.Tl];r1.Tr=[r1.Tr T2.Tr];
  r1.Ti=[r1.Ti T2.Ti];
  r1.adof=cellfun(@(x,y)[x;y],r1.adof,T2.adof,'uni',0);
 end
end % loop on sets
 % d1=fe_eig({m,k,[r1.Tl r1.Tr r1.Ti],SE.DOF},2);cf.def=d1; 
 % cf=feplot;cf.def=struct('def',[r1.Tl r1.Tr r1.Ti],'DOF',SE.DOF,'adof',vertcat(r1.adof{:}))
if ~isempty(kd);ofact('clear',kd);end
   
% fecom('shownodemark',SE.TR.adof,'marker','o','color','r')
% cf=feplot(10);feplot(cf,SE,TR);fecom('showfimdef');
TR=struct('def',[r1.Tl r1.Ti r1.Tr],'DOF',SE.DOF, ...
    'adof',[r1.adof{1};r1.adof{3}; r1.adof{2}], ...
    'info',[size(r1.Tl,2) size(r1.Ti,2) size(r1.Tl,2)],'P2info',RO.P2info);
%% fix NoMass % fecom('shownodemark',TR.DOF(diag(m)==0))
mr=feutilb('tkt',TR.def,m);i1=find(diag(mr)==0);
if ~isempty(i1)
 i2={fe_c(r1.adof{1},TR.adof(i1),'ind') fe_c(r1.adof{2},TR.adof(i1),'ind') ...
     fe_c(r1.adof{3},TR.adof(i1),'ind')};
 i2{1}=unique(vertcat(i2{1:2}));i2{2}=i2{1};
 if ~isempty(i2{3});
     fprintf('Clipping interior DOF %s\n',sdtm.toString(fe_c(r1.adof{3}(i2{3}))'))
 end
 r1.Tl(:,i2{1})=[];r1.Tr(:,i2{2})=[];r1.Ti(:,i2{3})=[];
 r1.adof{1}(i2{1})=[];r1.adof{2}(i2{2})=[];
 if ~isempty(i2{3});r1.adof{3}(i2{3})=[];end

 TR=struct('def',[r1.Tl r1.Ti r1.Tr],'DOF',SE.DOF, ...
    'adof',[r1.adof{1};r1.adof{3}; r1.adof{2}], ...
    'info',[size(r1.Tl,2) size(r1.Ti,2) size(r1.Tl,2)]);
end
if 1==2
 z=SE;z.TR=TR;
 z.K=cellfun(@(x)z.TR.def'*x*z.TR.def,z.K,'uni',0);
 s=2i*pi; Z=feutilb('sumkcoef',z.K,[s^2 1 s 1i]);diag(Z)
 [u,s]=svd(Z);s=diag(s);[un1,i1]=max(abs(u(:,s<1e-3)));fe_c(TR.adof(i1))
end
if 1==2 % Check matrix diagonal and topology
  figure(11);clf;semilogy(feutilb('dtkt',TR.def,SE.K));
  ii_plp(cumsum([RO.P2info(:,1);RO.P2info(:,3);RO.P2info(:,2)])*[1 0]); 
  legend(SE.Klab)
  % check dynamic stiffness
  dt=1e-3;Z=feutilb('sumkcoef',feutilb('tkt',TR.def,SE.K),[dt^-2 1 1/dt]);
  [U,s]=svd(full(Z));s=flipud(diag(s));U=fliplr(U);
  d1=struct('def',U,'DOF',TR.adof,'data',s,'TR',TR);feplot(SE,d1);
  ic=1:4;fe_c(d1.DOF(any(abs(d1.def(:,ic))>.1*norm(d1.def(:,ic),'inf'),2)))
  ic=1;d2=fe_def('subdef',TR,fe_c(d1.DOF,d1.DOF(any(abs(d1.def(:,ic))>.1*norm(d1.def(:,ic),'inf'),2)),'ind'));
  d2=fe_def('subdof',d2,(21:23)'/100);d2.DOF=d2.DOF-.2;cf=feplot;cf.def=d2;
  % check topology
  figure(10);ii_plp('spy',struct('Node',[],'Klab',{SE.Klab},'K',{feutilb('tkt',TR.def,SE.K)}));
  figure(10);ii_plp('spy',SE);
end
%disp('zzz');TR=struct('def',Case.T,'adof',Case.DOF,'DOF',TR.DOF)

SE.TR=TR; SE.K=feutilb('tkt',SE.TR.def,SE.K);SE.DOF=SE.TR.adof;
% d1=fe_eig(SE);d1.TR=SE.TR;cf=feplot; cf.def=d1;
[Case,name]=fe_case(SE,'getcase');
for j1=find(strcmpi(Case.Stack(:,1),'dofload'))'% Make sure load propagated
 d1=feutilb('placeindof',TR.DOF,Case.Stack{j1,3});
 d1.def=TR.def'*d1.def; d1.DOF=TR.adof;
 Case.Stack{j1,3}=d1;
end
for j1=find(strcmpi(Case.Stack(:,1),'sensdof'))'% Make sure sensors propagated
 Case.Stack{j1,3}.UseSE=2;
end

Case.T=speye(length(SE.DOF));Case.DOF=SE.DOF; Case.mDOF=SE.DOF; % Save as SE
Case=feutil('rmfield',Case,'cGL','bas','GroupInfo','jGroup','DofPerElt');
SE=stack_set(SE,'case',name,Case);

assignin('caller','out',SE);

elseif comstr(Cam,'redlist') 
%% #DFTRedList : learning phase for reduction 

 model=varargin{carg};carg=carg+1;
 RO=varargin{carg};carg=carg+1;
 if ~isfield(RO,'BuildList');RO.BuildList=RO.peig;end % DV compatibility
 if ~isfield(RO,'Phase2');RO.Phase2={'fe_homo','DftRedP2'};end
 if ~isfield(RO,'matdes');RO.matdes=[2 1];end
 if isfield(RO,'Assemble')&&~isempty(RO.Assemble)
  [SE,Case,Load]=fe_case(model,RO.Assemble);
 else
 [SE,Case,Load]=fe_case(model, ...
     sprintf('assemble -matdes %s -SE NoTload',sprintf('%i ',RO.matdes)));
 end
 def=[];
 % Keep_SleeperBandGap mode Anti resonance 'Eig2 -pcond=1',[5 3 0]
 for j2=1:size(RO.BuildList,1) % Peig as list see dv15 cfg

     [R2,st,CAM]=cingui('paramedit -DoClean',[ ...
      'fmax(#%g#"handling of callback for each step")' ...
       ],{RO,RO.BuildList{j2,1}}); Cam=lower(CAM);
     SE3=SE;SE3.Case=Case;SE3.Load=Load;
     val=RO.BuildList{j2,2};
     if ~isfield(val,'zCoef'); % Deal with matrix building
      fprintf('\n%s ---\n',RO.BuildList{j2,1});
     elseif iscell(val.zCoef) % Possibly deal with zCoef matrix definitions
       r1=[{''},SE.Klab(:)'];
       for j3=1:size(val.zCoef,1);
        r1(j3+1,:)=[val.zCoef(j3) num2cell(val.zCoef{j3,3})];
       end
       fprintf('\n%s --- \n%s \n',RO.BuildList{j2,1},comstr(r1,-30));
        SE3.Klab=val.zCoef(:,1)';SE3.Opt=1;SE3.Opt(2,1:size(val.zCoef,1))=vertcat(val.zCoef{:,2})';
        SE3.K=cellfun(@(x)feutilb('sumkcoef',SE.K,x),val.zCoef(:,3),'uni',0)';
        val=val.val;
     else; error('Not expected');
     end
     if comstr(Cam,'tstat')% {'tstat',linspace(0,.5,RO.Static)'}
      def2=fe_cyclic(CAM,SE3,val); def2.data(:,2)=0;
      def2.LabFcn='sprintf(''%i@ %.1f rad/m'',ch,2*pi*def.data(ch,1)/.6)';
     elseif comstr(Cam,'feval')%  [d1,SE]=p_pml('SolveDfrf',mo1,RD);
      i1=find(cellfun(@(x)ischar(x)&&strncmpi(x,'@',1),val));
      r1=val;% Allow callback to use r1 variable name rather than val
      for j1=i1(:)';val{j1}=eval(val{j1}(2:end));end
      def2=feval(val{:});  % external solver. Example p_pml SolveDfrf
     else % {eig2 -pcond=1',[5 3 0]}
      def2=fe_cyclic(CAM,stack_set(SE3,'info','EigOpt',val));
     end
     if ~isempty(R2.fmax);
         i1=def2.data(:,1)<R2.fmax; if nnz(i1)==0;i1=1:4;end
         def2=fe_def('subdef',def2,i1);
     end
     def=fe_def('appenddef',def,def2); clear def2
     def.LabFcn='sprintf(''%i@ %.1f %i'',ch,def.data(ch,1:2))';
 end

 %% Now do phase 2 : reduction fe_homo('DftRedP2');
 % cf=feplot(SE,def);fecom showfimdef;fecom('coloralpha');
 out1=def; % Allow export of reference
 if size(SE.K{1},1)*2==size(def.def,1) % long vector go to double format
  def.def=reshape(def.def,size(def.def,1)/2,[]); def.DOF=def.DOF(1:length(def.DOF)/2);
 end
 feval(RO.Phase2{:});
 
elseif comstr(Cam,'redrest') 
%% #DFTRedRest : restitution on unreduced DOF
 SE=varargin{carg};carg=carg+1;
 def=varargin{carg};carg=carg+1;
 i1=1:size(SE.TR.def,2);
 def.def=[SE.TR.def*def.def(i1,:);SE.TR.def*def.def(i1+i1(end),:)];
 def.DOF=[SE.TR.DOF;SE.TR.DOF+.5];
 out=def;
    
elseif comstr(Cam,'ak') 
%% #DFTAk [Ck,Cn]=fe_homo('dftAk -nc(i) -ik (k)',cx) wave domain amplitude

if carg<=nargin&&isfield(varargin{carg},'Y'); cx=varargin{carg};carg=carg+1;
else; % dirac on first point
 cx=struct('X',{{0}},'Xlab',{{'inx'}},'Y',ones(1,1)); 
end
if carg<=nargin&&isstruct(varargin{carg}); RO=varargin{carg};carg=carg+1;
else; RO=struct;
end
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   'nx(#%g#"Number of cells for x spatial transform") '...
   'skx(#%g#" subsamples of kx spectrum")'...
   'ny(#%g#"Number of cells for y spatial transform") '...
   'nz(#%g#"Number of cells for y spatial transform") '...
   'iky(#%g#" subsamples of ky spectrum")'...
   'inx(#%g#" subsamples of slices")'...
   ],{RO,CAM});Cam=lower(CAM);
if isempty(RO.nx)&&~isfield(RO,'kcx'); error('You must at least give nx');end
if ~isempty(RO.skx) % Place ik points in the half spectrum
 RO.nx=ceil(RO.nx/(RO.skx-1)/2)*(RO.skx-1)*2+1;
 RO.ikx=2:RO.skx:RO.nx/2;
else; % Keep all points
 if rem(RO.nx,2)==0;RO.nx=RO.nx+1;end % odd
end

if isfield(RO,'kcx') % Given possibly irrregular Kcx
 
 if ~isfield(RO,'inx'); RO.inx=(0:RO.nx-1)';end
 out=dftu('Enk',RO.inx,RO.kcx);

elseif 1==2 % old/test
 kcx=RO.kcx(:);
 if any(kcx>pi);error('Cannot have kcx>pi');end
 
 
 % One seeks to have a finite response with spectrum independent of N
 % adding zeros should not change the response

  dk=diff([0;kcx])+diff([kcx;pi]);
  RO.Enk=reshape([cos(inx*kcx');-sin(inx*kcx')]*diag(dk),size(inx,1),[])/2/pi;  

 
 an=zeros(RO.nx,size(cx.Y,2));
 an(cx.X{1}+1,:)=cx.Y; % Place non-zero at values given in cx (0 numbering)
 Ak=fft(an); %N=length(Ak);
 
 X=RO.kcx;if X(1)~=0;X=[0;X];end % Add kcx=0 if needed
 kc=X; X(kc~=0)=2*pi./kc(kc~=0);X(kc==0)=1;
 X=[X kc];
 
 % Create symetric kappa in [0 Pi] for kappa > Pi
 ix=find(X(:,2)>pi);
 Y=X(ix(1):length(X),2);
 Y=(Y-pi);
 X(:,2)=unique([X(1:ix(1)-1,2);Y]);
 X(:,1)=[1;sort(2*pi./X(2:end,2),'descend')];
  Ak=Ak(1:RO.nx);
 out=cx; out.X{1}=X;  out.Xlab{1}={'ncx','kcx'};out.Y=Ak;
 % Create matrix E %E=fe_cyclic_modif('inttransform',X(:,2)); %[Re(0) Re(1) Im(1) ...]
 kappa=transpose(X(:,2));N=length(kappa)*2;
 if ~isfield(RO,'inx'); RO.inx=(0:RO.nx-1)';end
 r1=zeros(length(RO.inx),N+1); % harmonic to physical transform
 for j1=1:length(RO.inx)
  %k=0:N/2; 
  r1(j1,[1 2:2:(end-2)])=cos(kappa*RO.inx(j1)).*[1 repmat(2,1,length(kappa)-2) 1];
  r1(j1,3:2:end)=-2*sin(kappa*RO.inx(j1));
 end
E=r1;
A=zeros(2*length(X)+1,1);
for i=1:length(X)-1
    A(2*i+1)=X(i+1,2)-X(i,2);
    A(2*(i+1))=X(i+1,2)-X(i,2);
end
A(1)=X(1,2);
A(2)=X(1,2);
A(end-1)=pi-X(end,2);
A(end)=pi-X(end,2);
B=diag(A);
E=E*B;
out.Enk=[E(:,1)*[1 0] E(:,2:end)]; %[Re(0) 0, Re(1) Im(1), ...]
out.inx=RO.inx; 
else % old case for regular spacing (should be removed)
 an=zeros(RO.nx,size(cx.Y,2));
 an(cx.X{1}+1,:)=cx.Y; % Place non-zero at values given in cx (0 numbering)
 Ak=fft(an); N=length(Ak);
 X=[[1;N./(1:(N-1)/2)']  ...  % Ncx from (long to short wavelength)
    (0:(N-1)/2)'/N]; % (/m) not rad/m
 Ak=Ak(1:size(X,1));
 out=cx; out.X{1}=X*diag([1 2*pi]); out.Y=Ak; out.Xlab{1}={'ncx','kcx'};
 E=fe_cyclic('reftransform',length(an))/length(an); %[Re(0) Re(1) Im(1) ...]
 out.Enk=[E(:,1)*[1 0] E(:,2:end)]; %[Re(0) 0, Re(1) Im(1), ...]
 if ~isempty(RO.inx); out.Enk=out.Enk(RO.inx,:); out.inx=RO.inx; end
end


if isfield(RO,'ikx')&&~isempty(RO.ikx) % Reinterpolate from subsampled k values
 out.X{1}=X(RO.ikx,:); out.Y=out.Y(RO.ikx,:);
 r1=of_time('lininterp',[X(RO.ikx,2) eye(length(RO.ikx))], ...
    X(:,2),zeros(1,3));
 out.interp=r1;
 % Just a few harmonics
 r2=[1;1]*sum(r1);i2=[RO.ikx*2-1;RO.ikx*2];
 out.ikx=RO.ikx;
 out.Ensk=out.Enk(:,i2(:))*diag(r2(:)); % subsampled
 % Reinterpolate integration rule
 r2=zeros(size(r1)*2);
 r2(1:2:end,1:2:end)=r1; % real interp
 r2(2:2:end,2:2:end)=r1; % imag interp
 out.Enik=out.Enk*r2; % Back transform from [Re I] (interpolated)
end
if nargout>1 % Provide the spatial input
 out1=cx;out1.X{1}=(0:size(an,1)-1)'; out1.Y=an;
end

elseif comstr(Cam,'homo') 
%% #DFTHomo : basic strategy for homogeneisation
 
model=varargin{carg};carg=carg+1;
RO=varargin{carg};carg=carg+1;
RO.nc=[0 1e3 0];d1=fe_homo('dft eig',RO,model);
RO.nc=[1e3 0 0];d2=fe_homo('dft eig',RO,model);
RO.nc=[0 0 1e3];d3=fe_homo('dft eig',RO,model);

d1=fe_def('AppendDef',fe_def('subdef',d1,1:2:6), ...
    fe_def('subdef',d2,1:2:6),fe_def('subdef',d3,1:2:6));

%i1=abs(mean(fe_c(d1.DOF,.02)*d1.def))./mean(abs(d1.def));
%d1=fe_def('subdef',d1,i1>1); % shear yz, shear xz, axial yy

%model.K={};model=fe_case('assemble NoT -matdes 2 1 -SE',model);
% For uniform properties the strain state is uniaxial
% The same is not true for non-uniaxial properties
RO.pl1=RO.pl; % First iteration
rho=RO.pl(1,12); [r1,unit,r3]=fe_mat('type',RO.pl(1,2));
iter=3;
while iter;iter=iter-1;
 mou=model;mou.pl=RO.pl1;
 mou.Elt=feutil('setgroupall matid 1 proid1',mou);
 RO.nc=[0 1e3 0];du1=fe_homo('dft eig',RO,mou);
 RO.nc=[1e3 0 0];du2=fe_homo('dft eig',RO,mou);%comp12('rveDefFlex',cf,d2);
 RO.nc=[0 0 1e3];du3=fe_homo('dft eig',RO,mou);%comp12('rveDefFlex',cf,d3);
 du1=fe_def('AppendDef',fe_def('subdef',du1,1:6), ...
    fe_def('subdef',du2,1:6),fe_def('subdef',du3,1:6));

 % Adjust phase of corner node
 data=fe_case(model,'getdata','Symmetry');
 i1=intersect(intersect(data.IntNodes(:,1),data.IntYNodes(:,1)), ...
    data.IntZNodes(:,1));
 %r1=fe_c(d1.DOF,i1)*d1.def;r1=phaseb(r1)
 %d1.def=d1.def*
 a=ii_mac(d1,du1,'mac paira -combine .001');
 i1=vertcat(a.GroupB{:});i1=i1(:,1); % Matched modes
 mou=fevisco('matsplit -type"ortho" matid1',mou);
 mou.K={};mou=fe_case(mou,'assemble -matdes 1 -1 -SE -NoT');
 mou.K(1:2)=[];mou.Klab(1:2)=[];
 r3=stack_get(mou,'info','MatSplit','get');
 r1=feutilb('dtkt',real(du1.def(:,i1)),mou.K)+ ...
    feutilb('dtkt',imag(du1.def(:,i1)),mou.K);
 r2=[sum(r1,2) ((du1.data(i1,1)*2*pi).^2) ((d1.data(:,1)*2*pi).^2)];

 r4=pinv(r1,1e-3)*r2(:,3);r4(abs(r4)>10)=1;
 dd1=reshape(r3.dd1*r4,6,6);
 comp12('matinfo',dd1)
 idd=[1  7 8  13 14 15   19 20 21 22   25 26 27 28 29  31 32 33 34 35 36]; %+36*(j1-1);
 RO.pl1=[1 fe_mat('m_elastic',unit,3) dd1(idd) rho];
end
fprintf('Initial\n');
[a,b,c,dd]=p_solid('buildconstit',[1;1],RO.pl,mou.il);dd=dd.dd;
comp12('matinfo',dd)

elseif comstr(Cam,'vg');[CAM,Cam]=comstr(CAM,5);
%% #DFTVg : estimate group velocity from a dispersion curve

%  'xxx see nida13 xxx'
% Group velocity is d omega/ d k
% physical wave length is 1/(n_c dx) = k/dx

hist=varargin{carg};carg=carg+1;
RO=struct('ind',1:3,'unit','SI');
 % 1./(hist.X{1}/2/pi)/RO.Lc(1) % check Nc 
 %deriv=feval(ii_mmif('@getOp'),struct('type','velc'),hist);

if carg<=nargin;RO=sdth.sfield('addmissing',varargin{carg},RO);carg=carg+1;end
if isfield(hist,'dx');dx=hist.dx;
else; dx=norm(hist.Range.CellDir(:,1));% xxx check
    if strcmpi(hist.Xlab{1},'kx'); dx=1;end
end

[Fu,dFu]=feval(nlutil('@cubeInterp'),'pppchip',hist);
 h2=hist;
 %h2.Y=2*pi*diag(sparse(dx./diff(hist.X{1})))*diff(h2.Y);
 %h2.X{1}(end)=[];h2.Y(:,:,2)=hist.Y(1:end-1,:); h2.DimPos=[1 3 2];
 h2.Y=dFu(hist.X{1})*2*pi;h2.Y(:,:,2)=hist.Y;h2.DimPos=[1 3 2];
 h2.X{3}={'Vg','m/s',[];'Freq','Hz',[]}; h2.Xlab{3}={'info'};
 %iicom('curveinit','V_g',h2); iicom('polar-1'); 
 coef=fe_mat(sprintf('convert%sSI',RO.unit),5);
 if isfield(RO,'rel')
  if ~isfinite(RO.rel(1)); RO.rel=squeeze(h2.Y(end,:,1));end
  h2.Y(:,:,1)=h2.Y(:,:,1)*diag(sparse(1./RO.rel));
  h2.X{3}(1,1:3)={sprintf('Vg/(%.2g m/s)',RO.rel*coef),'',[]};coef=1; 
 end
 figure(1);
 if isfield(RO,'hold');hold(RO.hold);end
 semilogx(h2.Y(:,RO.ind,2)/1000, ...
     h2.Y(:,RO.ind,1)*coef)
 xlabel('Frequency (kHz)');ylabel(sprintf('%s ',h2.X{3}{1,1:2}));axis tight
 if isfield(RO,'hold');hold('off');end
 if nargout>0; out=RO;end

elseif comstr(Cam,'deflong');[CAM,Cam]=comstr(CAM,5);
%% #DFTdeflong : transform to long format with Imag in DOF >.50
 def=varargin{carg};carg=carg+1;
 i1=(rem(def.DOF,1)>.5);
 if any(i1); % already long
 else
   def.DOF=[def.DOF;def.DOF+.5];
   def.def=[real(def.def);imag(def.def)];
 end
 out=def;

elseif comstr(Cam,'rk');[CAM,Cam]=comstr(CAM,3);
%% #DFTRK : typical wave number ranges
if carg<=nargin;RO=varargin{carg};else;RO=struct;end
if comstr(Cam,'lx')
  % line along x
  [CAM,Cam,i1]=comstr('lx',[-25 1],CAM,Cam);if isempty(i1);i1=30;end
  out=struct('val',linspace(0,pi,i1+2)','lab',{{'kcx'}});
  out.val([1 end])=[];
elseif comstr(Cam,'circ') 
%% %DftRkCircle in the wave domain
% used to be DispK RO=nida13('DispK',RO) : build range
% .nc_lim : 20 points at angle 0
% .angle = linspace(0,90,10);
% Range=struct('val',k_chos,'lab',{{'ncx','ncy'}});

if ~isfield(RO,'nc_lim');error('Expecting nc_lim range in [2 Inf]');end
if ~isfield(RO,'angle');RO.angle=0;warning('Using angle=0');;end

if ischar(RO.nc_lim)&&strncmpi(RO.nc_lim,'wangle',6);
 % Basic angular study for propagation at high frequencies
  RO.angle=[1 5:5:85 89]; 
  [CAM,Cam,nc]=comstr('nc',[-25 2],RO.nc_lim,lower(RO.nc_lim));
  if isempty(nc);nc=[100 101];
  elseif length(nc)<2;nc(2)=nc(1)*1.01;
  end
  if ~isfield(RO,'MatPl');RO.MatPl='';end
  RO.HistName=sprintf('%s %s lambda=%.0f mm',RO.MatPl,'Wangle',nc(1)*RO.lc(1));
elseif length(RO.nc_lim)==2&&diff(RO.nc_lim); 
  nc=logspace(log10(RO.nc_lim(1)),log10(RO.nc_lim(2)),20);
else; nc=RO.nc_lim;
end

nc=struct('val',nc(:),'lab',{{'nc'}},'type','vect');
ang=struct('val',RO.angle(:),'lab',{{'angle'}},'type','vect');
if min(nc.val(:))<2; warning('expecting nc>2 to avoid aliasing');end

Range=fe_range('grid',{nc;ang});

% Use ncx,ncy as in DFT commands
r1=Range.val;

if 1==1
 % Go to physical wavelength, assuming nc given along lc(1)
 % Then define cell periodicities
 kappa=2*pi./(r1(:,1)*RO.lc(1)); % Physical |k|
 r2=kappa(:,[1 1]).*[cosd(r1(:,2)) sind(r1(:,2))]; % physical k->
 i2=(r2~=0);r2(i2)=2*pi./r2(i2); % physical Ld
 r2=r2*diag(1./RO.lc(1:2)); % Cell periods
 Range.val=r2; Range.lab={'ncx','ncy'};

elseif 1==2
 Range.val=Range.val(:,[1 1]).*[cosd(Range.val(:,2)) sind(Range.val(:,2))];
 Range.lab={'ncx','ncy'};
else % Generate a circle in the wavenumber domain
 % Periodic wave number
 r2=(r1(:,[1 1]).^(-1)).*[RO.lc(1)*cosd(Range.val(:,2)) RO.lc(2)*sind(Range.val(:,2))];
 % Verify physical wavenumber : kx=2pi/l_x=2pi k/dx should be on a circle
 %kxyz=r2*2*pi*diag(1./RO.lc(1:2));
 % figure(1);plot(kxyz(:,1),kxyz(:,2));axis equal
 r2(r2~=0)=1./r2(r2~=0); r2(r2==0)=1; % k=0 <-> nc =1; 
 Range.val=r2; Range.lab={'ncx','ncy'};

end

% Angular wavenumber rad/[len] = 2*pi*kappa/dx
% Wavelength [len] in cells
% xxx lambda=Range.val; %lambda(lambda==1)=0;lambda=lambda*diag(RO.lc(1:2));
% ld=sqrt(sum(lambda.^2,2))
% RO.X=[Range.val r1(:,2) 2*pi./lambda lambda];
lambda=Range.val; i1=lambda==1|lambda==0;
lambda=lambda*diag(RO.lc(1:2));lambda(i1)=0; % physical wavelength
%kappa=2*pi./sqrt(sum(lambda.^2,2)); % norm of physical wave vector
kappa_xyz=lambda;kappa_xyz(~i1)=2*pi./lambda(~i1); kappa=sqrt(sum(kappa_xyz.^2,2));
RO.X=[Range.val r1(:,2) kappa 2*pi./kappa];
r2='l';
if isfield(RO,'unit');
    r2=fe_mat(['convertSI' RO.unit],'length');r2=r2{4};
end
RO.Xlab={'ncx','ncy','angle',sprintf('kappa [rad/%s]',r2), ...
    sprintf('lambda [%s]',r2)};
RO.Range=Range;
out=RO;
    
elseif comstr(Cam,'dva')
 % DynavoieA : uneven spacing of wavelength
 [CAM,Cam,RO.div]=comstr('dva',[-25 1],CAM,Cam);if isempty(RO.div);RO.div=30;end
 [CAM,Cam,RO.ncxlim]=comstr('nclim',[-25 2],CAM,Cam);
 [CAM,Cam,RO.ncy]=comstr('ncy',[-25 2],CAM,Cam);
 
 r1=[300 10 9.95 2.01];RO.ncxlim(end+1:4)=r1(length(RO.ncxlim)+1:end);
 Range=struct('val',[1 logspace(log10(RO.ncxlim(1)),log10(RO.ncxlim(2)),RO.div),  ...
    logspace(log10(RO.ncxlim(3)),log10(RO.ncxlim(4)),RO.div-1)]','lab',{{'ncx'}});
 if ~isempty(RO.ncy)
  Range.lab{2}='ncy';Range.val(:,2)=RO.ncy; % simple periodicity
 end
 out=Range;
  
else;error('DftRK%s',CAM);
end
elseif comstr(Cam,'converge');[CAM,Cam]=comstr(CAM,5);
%% DftConverge : packaging of convergence test
mo1=varargin{carg};carg=carg+1;
RO=varargin{carg};carg=carg+1;
mo1=stack_set(mo1,'info','EigOpt',2);

if ~isfield(RO,'name');RO.name='Wave';end
if ~isfield(RO,'ncx');RO.ncx=1./linspace(1e-1,.49,30)';end
r1=fe_case(mo1,'getdata','Symmetry'); 
if isfield(r1,'IntYNodes'); Range.ncy=1; end

for jpar=1:size(RO.scale,1)
 mo2=mo1; r3=RO.scale(jpar,:);r3(end+1:3)=1;
 mo2.Node(:,5:7)=mo2.Node(:,5:7)*diag(r3(1:3));
 Range.ncx=RO.ncx*max(RO.scale(:,1))/RO.scale(jpar);
 r2=r1;r2.CellDir=r2.CellDir*diag(RO.scale(jpar,1:size(r2.CellDir,2)));
 mo2=fe_case(mo2,'cyclic','Symmetry',r2);
 [def,hist]=fe_homo('dftDisp -UseLong',mo2,struct('Range',Range));
 if jpar==1
  RD=struct('mno',(0:10)','cf',102,'ci',112,'ViewHist',hist,'ktype','kcx');
  fe_homo('dfpInitSelDef',mo2,def,RD);fecom('ShowFiCEvalz')
 end
 if jpar==1; hp=hist;
     hp.Xlab{3}={'lc','mm',[],'#st3{j2}=sprintf(''%.1fmm'',r2(j2)*1e3);'}; % xxx unit
 end
 %hist.X{1}(1)
 hp.X{3}=RO.scale(1:jpar,:);hp.Y(:,:,jpar)=hist.Y;
end

if nargout==0
 hp.Xlab{2}={'Mode','',[],'Mode %i'};hp.Ylab={'Frequency','Hz',[]};
 hp.DimPos=[1 3 2]; hp.X{1}=2*pi./hp.X{1}*1000; hp.Xlab{1}={'\lambda_x','mm',[]};
 ci=iiplot(111); iicom('curveInit',RO.name,hp); 
 iicom(ci,'ch',{'Mode','All'});iicom(ci,';xlog;ylog');
end



elseif comstr(Cam,'test');[CAM,Cam]=comstr(CAM,5);
%% #DFTTest : sample models
if carg<=nargin;RO=varargin{carg};else;RO=struct;end
if isempty(Cam)
 %% Solid cube
 if ~isfield(RO,'lc');RO.lc=[.1 .1 1];end
 if ~isfield(RO,'divide');RO.divide=[1 1 4];end
 model=femesh(['testhexa20 divide', sprintf(' %i',RO.divide)]);
 model.Node(:,5:7)=model.Node(:,5:7)*diag(RO.lc); % Cell dims
 model=fe_case(model,'fixdof','plane',.02);
 model.pl=m_elastic('dbval 100 steel -unit MM');
 model=fe_homo('DftBuild',model,RO.lc(1:2));
 model=stack_set(model,'info','EigOpt',[5 50 -1e-12]);
 model=fe_case(model,'DofLoad','In',8.03,'SensDof','Out',8.03);

 if isfield(RO,'RemY') 
  % Just in x direction
  data=fe_case(model,'getdata','Symmetry');
  data=feutil('rmfield',data,'IntYNodes');data.CellDir=data.CellDir(1);
  model=fe_case(model,'stack_set','cyclic','Symmetry',data);
 end
 if isfield(RO,'In')
  % Peridiodic frequency response
  if ~isfield(RO,'Freq'); RO.Freq=linspace(3,150,10)'*1e5; end;
  RO.InF=1.03;
  model=fe_case(model,'remove','In', ...
      'DofLoad','In',struct('def',1,'DOF',1.03),'fixdof','top','z==1');
  model=stack_set(model,'info','Freq',RO.Freq);
 end
 
 out=model; 
elseif comstr(Cam,'bar')
 if ~isfield(RO,'lc');RO.lc=.1;end
 model=femesh('testbar1 divide 1');
 model.Node(:,5:7)=model.Node(:,5:7)*diag(RO.lc); % Cell dims
 model=fe_case(model,'fixdof','plane',[.02;.03]);
 model.pl=m_elastic('dbval 100 steel -unit MM');
 model.il=p_beam('ConvertTo1',p_beam('dbval 112 rectangle .1 .1'));
 model=fe_homo('DftBuild',model,RO.lc(1));
 model=stack_set(model,'info','EigOpt',[5 50 1e3]);
 model=fe_case(model,'DofLoad','In',1.01,'SensDof','Out',1.02);
 model=feutil('setmat 100 eta .1',model);

 out=model; 
 
else; error('Test%s unknown',CAM);
end
else; error('DFT%s Unknown',CAM);
end

elseif comstr(Cam,'dfp');[CAM,Cam]=comstr(CAM,4);
%% #Dfp : generic feplot utilities associated with DFT computations
if comstr(Cam,'initseldef');[CAM,Cam]=comstr(CAM,11);
%% #DfpInitSelDef : initialize for viewing of deformations

model=varargin{carg};carg=carg+1;
def=varargin{carg};carg=carg+1;
if carg<=nargin;RO=varargin{carg};carg=carg+1;
else; RO=struct;
end
if nargout>0; 
elseif isfield(RO,'cf');
 if ~isa(RO.cf,'sdth');RO.cf=comgui('guifeplot -reset',RO.cf);end
else; RO.cf=feplot;
end

if isfield(RO,'cf')&&~isempty(RO.cf)&&isa(model,'struct');
   if isa(RO.cf,'sdth');elseif ishandle(RO.cf);RO.cf=get(RO.cf,'userdata');
   else; RO.cf=feplot(RO.cf);
   end
   feplot(RO.cf,'initmodel-back',model);
end
data=fe_case(model,'getdata','Symmetry');

sel=sdth.GetData(model);
if ~isempty(CAM);sel.Elt=feutil(['selelt ' CAM],sel);
elseif isfield(RO,'sel');sel.Elt=feutil(['selelt ' RO.sel],sel);
end
if ~isfield(RO,'LinFace'); RO.LinFace='-linface';end
sel=feutil(['getpatch new' RO.LinFace],sel);sel.selelt=CAM;
% sel.cna
if ~isfield(def,'Range');warning('Missing def.Range');return;
elseif isfield(def.Range,'Ck');RO.Ck=def.Range.Ck;
end

if ~isfield(RO,'mno'); % Fill in volume in all dirs
  if isfield(RO,'N'); % Give N 
  elseif isfield(RO,'Ck')
   if isfield(RO.Ck,'inx');RO.N=length(RO.Ck.inx);
   else; RO.N=size(RO.Ck.Enk);RO.N(end)=[];
   end
  elseif ~isfield(RO,'N'); 
    RO.mno=(-10:10)';RO.N=max(RO.mno)-min(RO.mno)+1;
  end
  st=cell(length(RO.N),1);
  for j1=1:length(st); 
   st{j1}=sprintf('lab "%s" min %i max %i NPoints %i', ...
       char(abs('x')-1+j1),-fix(RO.N(j1)/2)+[0 RO.N(j1)],RO.N(j1)+1);
  end
  r1=fe_range('Grid',st);RO.mno=r1.val;
else;% Locations given in mno
  if isstruct(RO.mno);r2=fe_range('buildgrid',RO.mno);RO.mno=r2.val;end
  RO.N=max(RO.mno)-min(RO.mno)+1;
end 
if isfield(model,'TR')% place observation and restitution
 d1=def;d1.TR=model.TR; RO.TR='mdl';
 sel=fe_def('iselcna',[],d1,sel,'new'); 
 r2=spalloc(size(sel.cna{1},1),size(sel.cna{1},2),0);
 sel.cna{1}=[sel.cna{1} r2;r2 sel.cna{1}]; % re/im
elseif isfield(def,'TR') 
 % Recover shapes here
 if size(def.def,1)==size(def.TR.def,2);def.def=def.TR.def*def.def;
     def.def=[real(def.def);imag(def.def)];
 else;def.def=[def.TR.def*def.def(1:size(def.TR.def,2),:); 
         def.TR.def*def.def(size(def.TR.def,2)+1:end,:)];
 end
 def.DOF=def.TR.DOF;def=feutil('rmfield',def,'TR');
 def.DOF=[def.DOF;def.DOF+.5];
end
if isfield(RO,'TR')&&isequal(RO.TR,'mdl')
elseif rem(def.DOF(end),1)<.5 % Used in sdtweb fe_homo('mono-harmonic')
 error('Not reimplemented case with re/im in same DOF'); 
else % place observation from full
 sel=fe_def('iselcna',[],def,sel,'new'); 
 i1=1:size(def.DOF,1)/2;
 sel.cna{1}=[sel.cna{1};sel.cna{1}(:,[i1+i1(end) i1])]; % re/im
end

if ~isempty(RO.N)
data.CellDir=getCellDir('vect',data);

RO.Nnode =size(sel.vert0,1);RO.Nmax=max(sel.Node);
sel.vert0=repmat({sel.vert0},size(RO.mno,1),1);
sel.Node=repmat({sel.Node},size(RO.mno,1),1);
sel.fs=repmat({sel.fs},length(sel.vert0),1);
sel.f2=repmat({sel.f2},length(sel.vert0),1);
%RO.shift=ones(size(RO.mno,1),length(RO.N)); % done on expand
for j1=1:size(RO.mno,1); % Pave with Nx,Ny,Nz
    r2=sel.vert0{j1}; r3=sel.Node{j1};
    for j2=1:min(size(data.CellDir,2),size(RO.mno,2)); % Offset for node position and id
      r2=r2+ones(size(r2,1),1)*((RO.mno(j1,j2))*data.CellDir(:,j2)');
      r3=r3+j1*RO.Nmax;%RO.mno(j1,j2)*RO.N(j2)*RO.Nmax;'xxx negative cells'
      %RO.shift(j1,j2)=exp(2i*pi*RO.mno(j1,j2)/RO.N(j2));
    end 
    sel.vert0{j1}=r2;sel.Node{j1}=r3;
    sel.fs{j1}=int32(double(sel.fs{j1})+(j1-1)*RO.Nnode);
    if ~isempty(sel.f2{1})
     sel.f2{j1}=int32(double(sel.f2{j1})+(j1-1)*RO.Nnode);
    end       
    %i2=II;
    %i2(:,1)=i2(:,1)+RO.mno(j1,1)*RO.Nnode; % xshift
    %i2(:,2)=i2(:,2)+(j1-1+RO.N-1)*RO.Nnode; % yshift
    %i2(:,3)=i2(:,3)+(j1-1+2*RO.N-2)*RO.Nnode; % zshift
end
sel.vert0=vertcat(sel.vert0{:});sel.Node=vertcat(sel.Node{:});
sel.fs=vertcat(sel.fs{:}); %sel.opt(1)=2;
sel.f2=vertcat(sel.f2{:}); %sel.opt(1)=2;
if ~isfield(RO,'Ck');RO.Ck=[];end
sel.StressObs=struct('Ck',RO.Ck, ... %'Shift',RO.shift, ...
    'DOF',feutil('getdof',[.01;.02:.03],sel.Node));
%sel.fvcs=fe_homo('@DftShow');

%% #DfpInitDef : initialize expansion data
% step 1 : dvert(@nvert0*(x,y,z),ik) : p,d,k
% step 2 : d_vert * E = dvert(@vert0*(x,y,z),n) : p,d,n
% step 3 : permute to obtain dvert([@vx,n @vy,n @vz,n] 

%cf.SelF{1}.opt(1)=3; % Color interp at node
 d1=def;
 d1.DOF=feutil('getdof',[.01;.02;.03],sel.Node);
 fun=fe_homo('@DftRest'); 
 def.Rest=struct('mno',RO.mno,'Ck',RO.Ck, ...%'Shift',RO.shift, ...
    'DOF',d1.DOF,'cna',sel.cna{1},'CellDir',[]);
 sel.cna{1}=speye(length(d1.DOF));
 d1.def=curvemodel('Source',def,'yRef',fun,'getXFcn',{fun,fun});
 d1.LabFcn='feval(fe_homo(''@DftRest''),''LabFcn'',def,ch)';
else;% No mno
 d1=def;
 sel=feval(sdtroot('@iselcna'),[],d1,sel,'new');
end
if isfield(RO,'type')&&comstr(RO.type,'dfrf')
%% RO.type='dfrf' % recombined FRF with (freq x k) order
 r2=d1.def.Source; r2.Rest.type='dfrfb';
 r2.Rest.kcx=dftu('getx',d1,struct('lab',{{'kcx'}}));
 r1=feval(fe_def('@omethods'),'xvec',d1,[],'jPar');
 if diff(r1(1:2)); % sort kappa x freq
  r1=reshape(d1.data(:,1),size(r2.Range.val,1),[])';r1=r1(:,1);
  r2.Rest.type='dfrf';
 else
  r1=reshape(d1.data(:,1),[],size(r2.Range.val,1));r1=r1(:,1);r2.type='dfrfb';
 end
 d1.data=r1; d1.Xlab{2}='Freq';% One freq
 fun=fe_homo('@DftRest'); r2.data=d1.data;
 d1.def=curvemodel('Source',r2,'yRef',fun,'getXFcn',{fun,fun});
elseif (~isfield(RO,'type')||isempty(RO.type))&& ...
        length(def.Xlab{2})>1&&strcmpi(def.Xlab{2}{2},'IndDef')
    RO.type='dfrf';
elseif isfield(RO,'type')
elseif isfield(RO,'sens')&&size(def.data,1)/length(unique(def.data(:,1)))>1.9 
  RO.type='dfrf';% non recombined FRF
else; RO.type='hist';
end

if nargout==0;
    if isfield(RO,'cf');cf=RO.cf;else;cf=feplot;end
    cf.DefF={};cf.SelF{1}=sel;
    cf.def=d1;clf(cf.opt(1));feplot
else;out=sel;
end
if isfield(RO,'ci')&&isfield(RO,'sens')&&strcmpi(RO.type,'dfrf')
 %% Display 2DFRF
 ci=comgui('guiiiplot',RO.ci);
 fe_homo('dfpSensObserve -dimpos 2 3 1',model,RO.sens,def, ...
     struct('ktype',RO.ktype,'PlotInfo','surface','ci',ci))
elseif isfield(RO,'ViewHist') % Also display dispersion diag 
 fe_homo('dfpViewHist',RO,RO.ViewHist);
end



elseif comstr(Cam,'viewhist');
%% #DfpViewHist : view history and possibly
% fe_homo('dfpViewHist',struct('cf',feplot),hist,def)
%  dyn_post(sprintf('plotDisp -jx=%i',RO.jx),hist,struct('cf',RO.ci)); 

RO =varargin{carg};carg=carg+1;
hist=varargin{carg};carg=carg+1;
if carg<=nargin;def=varargin{carg};carg=carg+1;else;def=[];end

if ~isfield(RO,'ci');RO.ci=1;end; figure(RO.ci);
RO.hist=hist;RO=ktypeCheck(RO,[]);hist=RO.hist;
r1=hist.X{1};st=hist.Xlab{1};
if iscell(st);st(~cellfun(@ischar,st))=[];
    if length(st)>1;st=sprintf('%s [%s]',st{1:2});else; st=st1{1};end
end
if size(hist.Y,3)==2 % If markers defined loop on markers
 r2=unique(hist.Y(:,:,2)); marker={'or','+g','db','sk','^m','<c','vy'};
 li={}; 
 for j1=1:length(r2);
     r3=hist.Y(:,:,1);r3(hist.Y(:,:,2)~=r2(j1))=NaN;
     li(j1*3+[-2:0])={r3,r1,marker{remi(j1,length(marker))}};
 end
 plot(li{:});xlabel('Omega [Hz]');ylabel(st);
elseif norm(hist.Y(isfinite(hist.Y)),'inf')>1e4&&size(hist.Y,3)==1
 plot(hist.Y/1000,r1);xlabel('Omega [kHz]');ylabel(st);
elseif size(hist.Y,3)>1;warning('Forced response not expected');
 plot(squeeze(hist.Y(:,:,strncmpi(hist.X{3},'Freq',4))),r1);
 xlabel('Omega [Hz]');ylabel(st);
else
 plot(hist.Y,r1);xlabel('Omega [Hz]');ylabel(st);
end
if isfield(hist,'axProp'); cingui('objset',gca,hist.axProp); end% Should move to fe_homo

if isfield(RO,'cf')&&sdtdef('isinteractive'); 
   fe_homo('dfpDatatip',RO.ci,RO); % datatip with update
end


elseif comstr(Cam,'viewcurve');
%% #DfpViewCurve : curve model for Q(n,omega) viewing in slice
%  For Q(k,omega) viewing use fe_homo('dfpSensObserve')
% Build a curvemodel for freq at DOF observation
% see dyn_post('PerFEII')
def=varargin{carg};carg=carg+1;%d2=cf.def.def;def=d2.Source; 
if carg<=nargin; RO=varargin{carg};carg=carg+1;
else; RO=struct;end
if isfield(RO,'CurInd');
elseif isfield(RO,'DOF');
  RO.CurInd=fe_c(def.DOF,RO.DOF,'ind');
  RO.CurInd(rem(def.DOF(RO.CurInd),1)>.5)=[];% Not the long DOF
elseif isfield(RO,'sens');def.CurInd=1:size(RO.sens.cta,1);
else; RO.CurInd=1:3;
end
if ~isfield(RO,'mno');RO.mno=0;end

if isfield(RO,'type')&&strncmpi(RO.type,'cfrf',4)
 %% Cfrfb : recompose all kappa
 if isfield(RO,'sens') % Observe as pseudo DOF to allow recombine
   def.def=[RO.sens.cta*def.def(1:size(RO.sens.cta,2),:);
       RO.sens.cta*def.def(size(RO.sens.cta,2)+1:end,:)];
   def.DOF=(1:size(RO.sens.cta,1))'+.01; def.DOF=[def.DOF;def.DOF+.5];
   def.CurInd=1:size(RO.sens.cta,1);
 else; def.CurInd=RO.CurInd;
 end
 fun=fe_homo('@DftRest'); 
 def.Rest=struct('mno',RO.mno,'kcx',dftu('getx',def,struct('lab',{{'kcx'}})), ...
     'type',RO.type);
 r1=feval(fe_def('@omethods'),'xvec',def,[],'jPar');
 if diff(r1(1:2)); % sort kappa x freq
  r1=reshape(def.data(:,1),size(def.Range.val,1),[])';r1=r1(:,1);
  def.Rest.type='cfrf';
 else
  r1=reshape(def.data(:,1),[],size(def.Range.val,1));r1=r1(:,1);def.type='cfrfb';
 end
 def.X={r1,[]};def.Xlab={'Frequency [Hz]','DOF'}; 
 if isfield(RO,'sens');def.X{2}=RO.sens.lab;end
 % C1.X={Freq,Kappa,DOF}
 if size(RO.mno,1)>1; def.Xlab{3}='Slice';def.X{3}=RO.mno;end
 def.Ylab='Disp(omega)';
 C1=curvemodel('Source',def,'yRef',fun, ...
     'getXFcn',repmat({fun},1,length(def.X)),'DimPos',1:length(def.X));
else % Basic case with extraction of periodic mode at specific node
 def.CurInd=RO.CurInd;
 fun=fe_homo('@DftRest'); 
 % C1.X={Freq,Kappa,DOF}
 r1=reshape(def.data(:,1),[],size(def.Range.val,1));
 [RO,def]=ktypeCheck(RO,def);
 def.X={def.data(1:size(def.data,1)/size(def.Range.val,1),1), ...
    def.Range.val,[]}; 
 def.Xlab={'Frequency [Hz]',def.Range.lab{1},'DOF'};
 def.Rest=struct('mno',RO.mno,'Ck',[],'type','q(omega,k,dof)'); def.Ylab='Disp';
 C1=curvemodel('Source',def,'yRef',fun,'getXFcn',{fun,fun,fun},'DimPos',1:3);
 if isfield(RO,'hist') % Add PoleLines from hist
  C1.Source.ID={histToId(RO,def)};
 end
 %C1.Source.CurInd=fe_c(C1.Source.DOF,42+(1:3)'/100,'ind');
end
if isfield(RO,'doGetData');C1=C1.GetData;end
if isfield(RO,'PlotInfo');
  r1=RO.PlotInfo;r1{cellfun(@(x)isequal(x,'@C1'),r1)}=C1; C1=feval(r1{:});
end
if nargout==0; ci=iiplot;iicom(ci,'curveinit','Disp',C1);
else; out=C1;
end

elseif comstr(Cam,'datatip');
%% #Dfpdatatip : initalize datatip cursor : fe_homo('dfpDatatip',ci)
 cf=varargin{carg};carg=carg+1;
 
 if ~isa(cf,'sdth') % Datatip for dispersion curve
  ga=get(cf,'CurrentAxes');
  r1=get(ga,'children');r1(strcmpi(get(r1,'type'),'text'))=[];
  R1=struct('X',1);R1.PostFcn={@dftTip,'RangeCh'}; 
  iimouse('datatipnew',r1(1),R1)
  if carg<=nargin&&isfield(varargin{carg},'cf'); % Associate feplot from RO.cf
    RO=varargin{carg};carg=carg+1;
    ua=v_handle('uo',ga);ua.cf=RO.cf;   
  end
  
 elseif cf.opt(1,3)==1 % iiplot
  cf.ua.DataTipFmt='out=sprintf(''xxx=%.2f\n%.2f Hz'',pos);out1=pos;';
  cf.ua.PostFcn={@dftTip,'RangeCh'}; 
  cf.ua.cf=feplot;
  %RO=struct('X',1);
  %iimouse('datatipnew',cf.ga,RO)
 end
 if 1==2
   evd=struct('Character','','Modifier',{[]},'Key','rightarrow');
   gf=11;
   hThis=getuimode(11,'KeyPressFcn');
   hThis.modeWindowKeyPressFcn(gf,evd,hThis,[]);
 end
elseif comstr(Cam,'cursor');[CAM,Cam]=comstr(CAM,7);
%% #DfpCursor 
 if comstr(Cam,'iihist')
  % fe_homo('dfpCursorIIHist',12) : fecom ch cursor display
  gf=varargin{carg};carg=carg+1;
  iimouse('cursoron',gf);
  go=findall(gf,'tag','cursor');uo=get(go,'userdata');
  uo(1).HistoryFcn={@fe_homo,'DfpCursorHistory'};
  set(go,'userdata',uo);
 elseif comstr(Cam,'history')
 %% #DfpCursorHistory : display appropriate channels
  xi=evt(1).XHistory;
  x=evt(1).xdata(:);x=[x (1:length(x))'];
  xi(:,2)=fix(of_time('lininterp',x,xi,zeros(1,3)));
  for j1=1:size(xi,1)
   sdtweb('_link',sprintf('fecom(''ch%i'')',xi(j1,2)),sprintf('%.1f Hz',xi(j1,1)))   
  end
 else; error('DFPCursor%s Unknown',CAM);
 end
elseif comstr(Cam,'sensobserve');
 %% #DfpSensObserve C1=fe_homo('dfpSensObserve -dimpos 2 3 1',SE,sens,dFrf,struct('ktype','kc'))
 SE=varargin{carg};carg=carg+1;
 sens=varargin{carg};carg=carg+1;
 def=varargin{carg};carg=carg+1;
 if carg>nargin;RO=struct;else; RO=varargin{carg};carg=carg+1;end
 
 [RO,def]=ktypeCheck(RO,def);
 if 1==2% Does not work for mode tracking
    RO.aFreq=def.data(:,1); def.data(:,1)=repmat((1:size(def.data,1)/size(def.Range.val,1))',size(def.Range.val,1),1);
 end
 if rem(def.DOF(end),1)>.5;def=fe_cyclic('defdouble',def);end% Long
 if ~isfield(SE,'TR')
 elseif size(SE.TR.def,2)==size(def.def,1);def.TR=SE.TR;
 else; error('Inconsistent TR');
 end
 C1=fe_case(CAM,SE,sens,def);% Problem with mode tracking 
 i1=cellfun(@(x)isequal(x,'jPar'),C1.Xlab);
 C1.X{i1}=def.Range.val(:,1);C1.Xlab{i1}={RO.defKtype{2:3},[]};
 if isequal(C1.X{1},'Freq');C1.X{1}={'Frequency','Hz',[]};end
 i1=1:length(C1.X); 
 i2=[find(strncmpi(C1.Xlab,'Fr',2)); find(strncmpi(C1.Xlab,'k',1))];
 i3=find(strncmpi(C1.Xlab,'ReIm',4));
 i1=[i2(:)' setdiff(i1,[i2;i3]) i3];
 C1.Y=permute(C1.Y,i1);C1.Xlab=C1.Xlab(i1);C1.X=C1.X(i1);
 if isfield(RO,'hist') % Add PoleLines from hist
  C1.ID={histToId(RO)};
 end
 if isfield(RO,'PlotInfo')
   C1.PlotInfo=ii_plp('PlotInfo 2D -type "surface" ',C1);
 end
 if isfield(RO,'PostObs')
  %RO.PostObs='C2=ii_mmif(''acc -struct'',C1);C2.ID=C1.ID;C2.name=''EdgeZ'';C1=C2;';
  eval(RO.PostObs); 
 end
 if nargout==0
  if ~isfield(RO,'ci');RO.ci=[];end;ci=comgui('guiiiplot',RO.ci);
  iicom(ci,'curveinit',C1.name,C1);
 else
  out=C1;
 end
 
else; error('DFP%s Unknown',CAM);
end

elseif comstr(Cam,'vf');[CAM,Cam]=comstr(CAM,3);
%% #VF utilities associated with virtual field method

if comstr(Cam,'gen')
%% #VF gen : numeric generation of virtual fields
% The interace where non-zero edge loads are applied should be clamped
% tested in xxx
% publis : SDT_trackers\publi_MCV
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
 'MatId(#%i#"Material to apply split")' ...
 'splitFcn("Ortho"#%s#"split function to call")' ...
 '-visc(#3#"to activate viscoelasticity")' ... forwarded
 ],{struct,CAM}); Cam=lower(CAM);
model=varargin{carg}; carg=carg+1;
data=varargin{carg}; carg=carg+1;
if carg<=nargin; RO.sel=varargin{carg}; carg=carg+1; else; RO.sel=''; end
%[vkh,vb,dof]=fe_homo('VFgen',mo1,data,'withnode{z<=15}');
if isempty(RO.MatId); pl=feutil('getpl',model); RO.MatId=pl(1); end

mo1=model; % working model, maybe split
if ~isempty(RO.sel)
 mo1=sdth.GetData(mo1); mo1.Elt=feutil(sprintf('selelt %s',RO.sel),mo1);
 mo1.Node=feutil('getnodegroupall',mo1);
end
if RO.visc; matdes=[2 1 4];else; matdes=[2 1];end
mo1=fe_case(mo1,sprintf('assemble -matdes %s -SE -NoT -reset',sprintf('%i ',matdes)));

% define model with split properties
mo2=fe_homo(sprintf('vfSplit%s Matid%i %s',RO.splitFcn,RO.MatId,CAM),mo1);
% Generate the virtual fields
if ischar(data); data=fe_homo(sprintf('vfdb %s',data),model);end
for j1=1:numel(data)
 vf1=elem0('VectFromDirAtDof',model,data(j1),mo2.DOF);
 if j1==1;vf=vf1; else; vf=fe_def('appenddef',vf,vf1); end
end
if isfield(data,'lab'); vf.lab={data.lab};
else; vf.lab=num2cell((1:numel(data))');
end
vf=feutilb('placeindof',mo2.DOF,vf);

% RHS
if RO.visc; vb=vf.def'*(mo1.K{2}+1i*mo1.K{3}); 
else; vb=vf.def'*mo1.K{2}; 
end
% LHS: virtual field applied to homogeneous stiffness parts
vkh=struct('DOF',mo2.DOF,'vf',vf,'def',{cell(1,length(mo2.K)-1)},'k0',mo1.K{2});
if RO.visc; vkh.k0=vkh.k0+1i*mo1.K{3}; end
for j1=1:min(size(vf.def,2),length(mo2.K)-1)
 vkh.def{1,j1}=vf.def(:,j1)'*mo2.K{j1+1};
end

% output
out=vkh; out1=vb;
if nargout>2; out2=mo2; end


elseif comstr(Cam,'solve'); [CAM,Cam]=comstr(CAM,6);
 %% #VFSolve - - -
 vkh=varargin{carg}; carg=carg+1;
 vb=varargin{carg}; carg=carg+1;
 def=varargin{carg}; carg=carg+1;
 
 def=feutilb('placeindof',vkh.DOF,def);
 
 if isreal(def.def)
  out=zeros(size(def.def,2),length(vkh.def));
%   for j1=find(reshape(~cellfun(@isempty,vkh.def),1,[]))
%    out(:,j1)=diag((vkh.def{j1}(:,j1:3:end)*def.def(j1:3:end,:)))\...
%     (vb(j1,j1:3:end)*def.def(j1:3:end,:))' ;
%   end
  for j1=find(reshape(~cellfun(@isempty,vkh.def),1,[]))
   out(:,j1)=diag((vkh.def{j1}*def.def))\...
    (vb(j1,:)*def.def)' ;
  end
  
 else
  out=zeros(size(def.def,2),length(vkh.def),2);
  for j1=find(reshape(~cellfun(@isempty,vkh.def),1,[]))
   for j2=1:size(out,1)
    r1=vb(j1,:)*def.def(:,j2);
    out(j2,j1,:)=...
     [ vkh.def{j1}*real(def.def(:,j2))  -vkh.def{j1}*imag(def.def(:,j2));
       vkh.def{j1}*imag(def.def(:,j2))  +vkh.def{j1}*real(def.def(:,j2))]\...
     [real(r1);imag(r1)];
    %      [((real(vb(j1,:))*qr(:,j2)-imag(vb(j1,:))*qi(:,j2))');
    %       ((real(vb(j1,:))*qi(:,j2)+imag(vb(j1,:))*qr(:,j2))')];
   end
  end
 end

elseif comstr(Cam,'db'); [CAM,Cam]=comstr(CAM,3);
%% #VFdb : database of virtual fields functions
model=varargin{carg}; carg=carg+1;
out=struct('sel','groupall','DOF',(1:3)'/100,...
 'dir',{{'x.*(max(x)-x)','y*0','z*0'};
 {'x*0','y.*(max(y)-y)','z*0'};
 {'x*0','y*0','z.*(max(z)-z)'}});

elseif comstr(Cam,'orth')
%% #VFOrth : numeric orthonormalization 
model=varargin{carg};carg=carg+1;
T=varargin{carg};carg=carg+1;
if carg<=nargin;RO=varargin{carg};carg=carg+1;else;RO=struct('v',1);end

if fix(RO.v)==2
%% Find mode with specific matrix set to zero

mo2=model;C2=mo2.Case; 

% Compute null space of target matrix
zCoef=ones(1,length(mo2.K));zCoef(1)=0; zCoef(1,RO.jmat+1)=0;
d1=fe_eig({mo2.K{1},feutilb('sumkcoef',mo2.K,zCoef),C2.T,mo2.DOF});
fe_homo('VFCheckT',d1.def,mo2)

% Keep most interesting
zCoef(1,RO.jmat+1)=1;
[T,fr]=fe_norm(d1.def,mo2.K{1},feutilb('sumkcoef',mo2.K,zCoef));
fe_homo('VFCheckT',T,mo2)

elseif fix(RO.v)==1
%% First attempt at making orthogonal
  k=feutilb('sumkcoef',model.K(2:end),ones(1,length(model.K)-1));
  K=feutilb('tkt',T,model.K);
  RunOpt.T=cell(length(model.K),2);
  for jmat=2:length(model.K);
  kj=K{jmat};
  dkj=feutilb('sumkcoef', ...
      K(setdiff(1:length(K),[1 jmat])),ones(1,length(K)-2));
  [U,s,v]=svd(full(dkj));s=diag(s); % null space of dkj
  % now max response in null space of dkj
  i1=find(s<s(1)*1e-7); %if isempty(i1);keyboard;end 
  [U2,s2,v2]=svd(U(:,i1)'*kj*U(:,i1));U2=U(:,i1)*U2;s2=diag(s2);
  T2=T*U2; % This is the subspace
  z=feutilb('dtkt',T2,model.K(2:end));z(z<0)=0;z=diag(1./max(z,[],2))*z;
  i3=find(z(:,jmat-1)>.9);T2=T2(:,i3);z=z(i3,:);
  if round(rem(RO.v,1)*10)==1 
    % 1.1%~isempty(T2);% order within acceptable
    T2=fe_norm(T2,model.K{1},k);
    z=feutilb('dtkt',T2,model.K(2:end));z(z<0)=0;z=diag(1./max(z,[],2))*z;
  end
  r3=cellfun(@(x)sprintf('%s z=%.1e',model.Klab{jmat},x),num2cell(s2(i3)),'uni',0);
  RunOpt.T(jmat,1:3)={T2 r3 z};
  end
  
  cf=feplot;cf.def=struct('def',horzcat(RunOpt.T{:,1}), ...
      'DOF',model.DOF,'lab',{vertcat(RunOpt.T{:,2})});
  r3=[model.Klab(2:end);num2cell(vertcat(RunOpt.T{:,3}))];
  r3=[[{''};num2cell((1:size(r3,1)-1)')] r3];
  comstr(r3,-17,'tab',...
      struct('fmt',{[{'%i'} repmat({'%.2f'},1,9)]},'HasHead',1, ...
      'name',sprintf('VfOrh %i',RO.v)))
  % missing table selection / feplot link
elseif RO.v==0 
%% Just display T and associated component repartition
  z=feutilb('dtkt',T,model.K(2:end));z(z<0)=0;z=diag(1./max(z,[],2))*z;
  d1=struct('def',T,'DOF',model.DOF,'lab',{cell(size(T,2),1)});
  for j1=1:size(T,2);
    [r3,i3]=sort(z(j1,:),'descend');
    d1.lab{j1}=sprintf('%s %.2f, %s %.2f ',model.Klab{i3(1)+1},r3(1),...
        model.Klab{i3(2)+1},r3(2));
  end
  cf=feplot;cf.def=d1;
  r3=[model.Klab(2:end);num2cell(z)];
  r3=[[{''};num2cell((1:size(r3,1)-1)')] r3];
  comstr(r3,-17,'tab',...
      struct('fmt',{[{'%i'} repmat({'%.2f'},1,9)]},'HasHead',1, ...
      'name',sprintf('VfOrh %i',RO.v)))
   
end


elseif comstr(Cam,'split'); [CAM,Cam]=comstr(CAM,6);
%% VFSplit - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if comstr(Cam,'@') % #VFSPlit@ (callback)-3
 stFun=regexp(Cam,'^.*@([a-zA-Z0-9_]*).*$','tokens'); stFun=stFun{1}{1};
 [CAM,Cam]=comstr(CAM,length(stFun)+2);
 [mo2,C2]=feval(stFun,CAM,varargin{:});
 
elseif comstr(Cam,'ortho') % #VFSplitOrtho -3
 mo2=varargin{carg};carg=carg+1;
 [CAM,Cam,RO.MatId]=comstr('matid',[-25 1],CAM,Cam);
 
 mo2.Elt=feutil(sprintf('set group all matid %i proid %i', ...
  RO.MatId,RO.MatId),mo2); % Homog
 mo3=fevisco(sprintf('matsplit MatId %i -partype1 -type"orthou" %s', ...
  RO.MatId,CAM),mo2);
 [mo2,C2]=fe_case(mo3,'assemble -matdes 2 -1 -SE -NoT');
 
else; error('fe_homo VFSplit %s unknown',CAM);
end
out=mo2; if nargout==1; out.Case=C2;else; out1=C2;end
 
elseif comstr(Cam,'checkt')
%% VFCheckT

 T=varargin{carg};carg=carg+1;
 model=varargin{carg};carg=carg+1;
 i2=2:length(model.K);
 if isfield(T,'def'); TR=T; T=TR.def;end
 r1=feutilb('dtkt',T,model.K(2:end));r1(r1<0)=0;r1=diag(1./max(r1,[],2))*r1;
 disp(svd(r1)')
 out=[model.Klab(i2);num2cell(r1)];
 if nargout==0 % Display fractions
  out=[[{''};num2cell((1:size(out,1)-1)')] out];
  comstr(out,-17,'tab',...
      struct('fmt',{[{'%i'} repmat({'%.2f'},1,9)]},'HasHead',1, ...
      'name',sprintf('VfOrh %i',0)))
  clear out;
 end
 
elseif comstr(Cam,'est'); [CAM,Cam]=comstr(CAM,6);
%% VFEst estimate modules based on given VF and frequency responses - - - - -

SE=varargin{carg};carg=carg+1;
dvirt=varargin{carg};carg=carg+1;
d1=varargin{carg};carg=carg+1;

if ~isequal(d1.DOF,SE.DOF);d1=feutilb('placeindof',SE.DOF,d1);end
if ~isequal(dvirt.DOF,SE.DOF);error('Mismatch');end

G=cellfun(@(x)dvirt.def'*x*d1.def,SE.K,'uni',0);
G{1}=dvirt.def'*SE.K{1}*d1.def*diag((d1.data(:,1)*2*pi).^2);
i2=cellfun(@(x)normest(x),G);i2=find(i2>1e-10*max(i2)); %disp(SE.Klab(i2));

r2=zeros(length(i2)-1,size(d1.def,2));
for jmode=1:size(d1.def,2)
  r3=cellfun(@(x)x(:,jmode),G(i2(2:end)),'uni',0);wK=horzcat(r3{:});
  %[{0} SE.Klab(i2(2:end));[num2cell(1:size(dvirt.def,2))' num2cell(r3)]]
  wM=G{1}(:,jmode);
  r2(:,jmode)=pinv(wK)*wM;
  %i3=abs(wM)>.01*norm(wM,'inf');wK(i3,:)\wM(i3)
end
  %[dvirt.lab(:) num2cell(G{1})]
ua=struct('table',{[SE.Klab(i2(2:end))' num2cell(r2)]}, ...
    'ColumnName',{num2cell(0:size(r2,2))});ua.ColumnName{1}='Mode';
ua.ColumnName(3,2:end)={'0.0'};

if nargout==0; comstr(ua.table,-17,'tab',ua);
else;out=ua;
end

else; error('VF%s Unknown',CAM);
end


elseif comstr(Cam,'rve');[CAM,Cam]=comstr(CAM,4);
%% #Rve : simple periodic homogeneization computations

if comstr(Cam,'pbc')||comstr(Cam,'simpleload')||comstr(Cam,'kubc')||comstr(Cam,'mubc')
%% #RveKubc #RvePcb #RveSimpleLoad: implementation of KUBC homogeneization method
 model=varargin{carg};carg=carg+1;
 if carg<=nargin;RO=varargin{carg};carg=carg+1;else; RO=struct;end
 if ~isfield(RO,'btype');RO.btype=CAM;end
 if ~isfield(RO,'alpha');RO.alpha=-100;end % Mass shift
 %DirDofInd={1:3,1:3,[1 2 0]};
 if ~isfield(RO,'epsl');RO.epsl=sp_util('epsl');end
 if ~isfield(RO,'silent');RO.silent=';';end
 if ~isfield(RO,'oProp');RO.oProp={};end
 if ~isfield(model,'nmap');model.nmap=vhandle.nmap;end
 model.nmap('NoLoadWarn')='on';model.nmap('NoMaterialWithLoss')='on';

 [z,RO,mo1,C1,Load,w,ft,carg]=safeInitDfrf([{CAM,model,RO},varargin(carg:end)],2);

 mo2=mo1;fe_case(mo2,'stack_rm','DofSet'); [Case,NNode,DOF]=fe_case(mo2,'gett');
 mo1.Case=C1; [R2,d2,mo1]=BuildUb(mo1,RO);c=R2.c;% Possibly deal with visco

 % cg=feplot(10);cg.model=model;cg.def=d2;
 % Now enforced motion but use MPCs
 % cg=feplot;cg.sel={['selface & innode' sprintf('%i ',data.IntNodes)],'colorfacew'}
 
 r2=fe_case(model,'stack_get','Dofset');
 if ~isempty(r2);r2(~ismember(lower(r2(:,2)),lower(RO.Load)),:)=[];end
 if ~isempty(r2); % Externally enforced motion
  for j1=1:size(r2,1)
   c1=fe_c(d2.DOF,r2{j1,3}.DOF);
   c=[c;c1];d2.lab(end+1)=r2(j1,2); %#ok<AGROW>
  end
  d2.def(:,end+(1:size(mo1.Case.TIn,2)))=mo1.Case.TIn;
 end
 % Now compute the forced response
 if strcmpi(RO.btype,'mubc') 
 elseif any(strcmpi(RO.btype,{'kubc'})) % Enforce motion on full edge (not for piezo SimpleLoad)
  c=find(sum(abs(c),1));c=sparse(1:length(c),c,1,length(c),length(d2.DOF));
  %[C1.T,C1.TIn] = fe_coor(c,[4 1 2],(c*d2.def));
  T = fe_coor('lu',struct('c',c,'TIn',d2.def));C1.T=T.T;C1.TIn=T.TIn;
  1;
  % [T,TIn] = fe_coor(c(:,fe_c(d2.DOF,Case.DOF,'ind')),[4 1 2],(c*d2.def));
  
 else % K=PBC : left disp=right disp
  if 1==2%~isfield(RO,'lut')||~RO.lut
   [T,TIn] = fe_coor(c(:,fe_c(d2.DOF,Case.DOF,'ind')),[4 1 2],(c*d2.def));
  else
   i2=fe_c(d2.DOF,Case.DOF,'ind');
   T = fe_coor('lu',struct('c',c(:,i2),'TIn',d2.def(i2,:)));
   if normest(c(:,i2)*T.T)>1e-10 % Problem if DOF only in TIn
    %sdtw('_ewt','Report lu problem');
    [T,TIn] = fe_coor(c(:,i2),[4 1 2],(c*d2.def));
   else;   TIn=T.TIn;T=T.T;   
   end
  end
  mo1.DOF=R2.DOF; C1.Stack=Case.Stack; 
  C1.T=Case.T*T;C1.TIn=Case.T*TIn;
  if size(C1.TIn,2)<size(d2.def,2);% Problem with TIn only DOF (piezo MFC) 
      C1.TIn(:,end+1:size(d2.def,2))=d2.def(:,size(C1.TIn,2)+1:size(d2.def,2));
  end
  C1.DOF=(1:size(T,2))'+.99;mo1.Stack(strcmp(mo1.Stack(:,2),'model.Case'),:)=[];
  mo1.Stack{strcmpi(mo1.Stack(:,1),'case'),3}=C1;
  d1=d2; d1.name='homog';d2.name='hete';% Save homogeneous cases
 end
 if isempty(fe_c(mo1.DOF,.21,'ind'))
  %% standard case with a loop on properties
  if ~isfield(RO,'toFun');RO.toFun=@toOrtho;C2=[];end
  % xxxREMOVE_FC  if isfield(model,'Visco') && ~isfield(RO,'Range'); model=stack_set(model,'info','MatSplit',RO.MatSplit);RO.Range=comp18('StepRveVisco', model); end

  if ~isfield(RO,'Range') % Use default range
   if strcmpi(RO.btype,'mubc')
    for j1=1:length(R2.li)
     d2.def(:,j1)=R2.li{j1}.TIn-ofact({mo1.K{strcmpi(mo1.Klab,'k')},R2.li{j1}.T},mo1.K{strcmpi(mo1.Klab,'k')}*R2.li{j1}.TIn,RO.oProp{:});
    end
   elseif isempty(C1.T);d2.def=C1.TIn;
   else
     d2.def=C1.TIn-ofact({mo1.K{strcmpi(mo1.Klab,'k')},C1.T},mo1.K{strcmpi(mo1.Klab,'k')}*C1.TIn,RO.oProp{:});
   end
   r1=d2.def'*mo1.K{ismember(lower(mo1.Klab),{'k','1','5'})}*d2.def/R2.V;
   if isfield(RO,'pl')
    if isfield(mo1,'unit');RO.unit=mo1.unit;end
    C2=RO.toFun(r1,[],sdth.sfield('addselected',struct,RO,{'type','unit','pl'}));% toOrtho
    if isfield(RO,'volume')
     C2.pl(strcmpi(m_elastic('propertyunittypecell',6),'rho'))= ...
      full(feutilb('tkt',sum(fe_c(mo1.DOF,.01))',mo1.K{1})/RO.volume);
    end
   else
    C2=RO.toFun(r1);% toOrtho
   end
  else
   if isfield(RO,'LearnVal')
   %% #RveKubcLearnVal start by doing a reduced model -3
    Range=RO.Range;Range.val=RO.LearnVal;
    d2.def=zeros(size(C1.TIn,1),size(C1.TIn,2)*size(Range.val,1));
    d2.LearnDD=[];
    for jpar=1:length(Range.val)
          Range.jPar=jpar;r2=Range.param.iVisco.data{Range.val(jpar)};
          r2(:,2)=cellfun(@real,r2(:,2),'uni',0);% Learning is elastic
          r2=feutilb('sumkcoef',mo1.K,r2); 
          def=C1.TIn-ofact({r2.Kh,C1.T},r2.Kh*C1.TIn,RO.oProp{:});
          r1=(def'*r2.Kh*def)/R2.V;
          if jpar==1;r20=r2;r20.r1=r1; disp(r1);
          end
          % qh_cc = - Kh_cc \ Kh_ci q_d
          % Note that C1.TIn not included since fixed
          d2.def(:,size(C1.TIn,2)*(jpar-1)+(1:size(C1.TIn,2)))=def-C1.TIn;
          d2.LearnDD=RO.toFun(d2.LearnDD,Range,r1);
    end
    [T,fr]=fe_norm(d2.def,mo1.K{1},r2.Kh);
    mo1.TR=struct('def',[C1.TIn T],'DOF',mo1.DOF,'data',fr/2/pi,'adof',(1:size(T,2))'+.99);
    C1.TIn=speye(size(mo1.TR.def,2),size(C1.TIn,2));
    C1.T=speye(size(C1.TIn,1));C1.T(:,1:size(C1.TIn,2))=[];
    mo1.K=feutilb('tkt',mo1.TR.def,mo1.K);
    % see Florian Conejos Thesis section 3.4.2 
    Range=RO.Range;     
    for jpar=1:size(Range.val,1)
      Range.jPar=jpar;      
      r2=feutilb('sumkcoef',mo1.K,Range.param.iVisco.data{Range.val(jpar)});
      kr=feutilb('tkt',C1.T,r2.Kh); 
      def=C1.TIn-C1.T* (kr\ ( C1.T'*(r2.Kh*C1.TIn))); % qh_cc = - Kh_cc \ Kh_ci q_d     
      r1=(def'*r2.Kh*def)/R2.V;
      C2=RO.toFun(C2,Range,r1);
    end
    1;
   else;
    %% #without_LearnVal_model reduction phase -4
    Range=RO.Range; C2=[] ; 
    for jpar=1:size(Range.val,1)
      Range.jPar=jpar;
      r2=feutilb('sumkcoef',mo1.K,Range.param.iVisco.data{jpar});
      if isfield(mo1,'TR')
       r3=feutilb('sumkcoef',mo1.KIn,Range.param.iVisco.data{jpar});
       def=C1.TIn-mo1.TR.def* (r2.Kh\r3.Kh); % qh_cc = - Kh_cc \ Kh_ci q_d
      elseif any(strcmp(RO.btype,{'subc','mubc'}))
        for j1=1:length(R2.li)
               def(:,j1)=R2.li{j1}.TIn-ofact({r2.Kh,R2.li{j1}.T},r2.Kh*R2.li{j1}.TIn,RO.oProp{:});
        end
      else
        def=C1.TIn-ofact({r2.Kh,C1.T},r2.Kh*C1.TIn,RO.oProp{:}); % qh_cc = - Kh_cc \ Kh_ci q_d
      end
      %out1=struct('dd',d2.def'*K*d2.def/R2.V,'lab',{d2.lab});
      %welas=@(x,k)real(x)'*real(k)*real(x)+imag(x)'*real(k)*imag(x);
      %welas(d2.def,r2.Kh)./real(d2.def'*r2.Kh*d2.def)
      %wdiss=@(x,k)real(x)'*imag(k)*real(x)+imag(x)'*imag(k)*imag(x);
      %wdiss(d2.def,r2.Kh)./imag(d2.def'*r2.Kh*d2.def)
      r1=def'*r2.Kh*def/R2.V;          
      C2=RO.toFun(C2,Range,r1);
    end
   end
  end
  out=C2;
  if nargout>1; 
     rb=feutilb('geomrb',mo1,[0 0 0],d2.DOF);
     d2.def=d2.def-rb.def* (rb.def\d2.def);out1=d2; 
  end
  out2=mo1; 
  return
 else % For piezo use preconditionning
  if ~isfield(RO,'alpha');RO.alpha=0;end
  K=feutilb('sumkcoef',mo1.K,[RO.alpha 1 0 0]);
  pcond=diag(K);
  pcond(pcond<0)=sqrt(mean(pcond(pcond>0))/-mean(pcond(pcond<0)))/100;
  pcond(pcond>0)=1;
  pcond=diag(sparse(pcond));
  %K1=pcond*K*pcond;
  d2.def=C1.TIn-ofact({K,pcond*C1.T},K*C1.TIn);
  out=d2;
 end
 
elseif comstr(Cam,'subc')
%% #RveSUBC : implementation of SUBC method
  
model=varargin{carg};carg=carg+1;
if carg<=nargin; RO=varargin{carg};carg=carg+1;else;RO=struct;end

if ~isfield(RO,'alpha');RO.alpha=-1000;end % Mass shift
[R2,d2,model]=BuildUb(model,RO);c=R2.c; d2.name='q_dK'; % d2 homogenous 
data=fe_case(model,'getdata','Symmetry');
RO.DOF=R2.DOF;
 r1=[1 0 0;0 0 0;0 0 0;0 0 0;0 0 1;0 1 0]';
model=fe_case(model,'reset','FSurf','F1', ...
    struct('sel',data.IntNodes(:,1), ...
    'def',repmat(r1,size(data.IntNodes,1),1), ...
    'DOF',feutil('getdof',data.IntNodes(:,1),(1:3)'/100)));
model=fe_case(model,'FSurf','F2', ...
    struct('sel',data.IntNodes(:,2), ...
    'def',repmat(-r1,size(data.IntNodes,1),1), ...
    'DOF',feutil('getdof',data.IntNodes(:,2),(1:3)'/100)));
if isfield(data,'IntYNodes')
 r1=[0 0 0;0 1 0;0 0 0;0 0 1;0 0 0;1 0 0]';
 model=fe_case(model,'FSurf','F3', ...
    struct('sel',data.IntYNodes(:,1), ...
    'def',repmat(r1,size(data.IntYNodes,1),1), ...
    'DOF',feutil('getdof',data.IntYNodes(:,1),(1:3)'/100)));
 model=fe_case(model,'FSurf','F4', ...
    struct('sel',data.IntYNodes(:,2), ...
    'def',repmat(-r1,size(data.IntYNodes,1),1), ...
    'DOF',feutil('getdof',data.IntYNodes(:,2),(1:3)'/100)));
end
if isfield(data,'IntZNodes')
 r1=[0 0 0;0 0 0;0 0 1;0 1 0;1 0 0;0 0 0]';
 model=fe_case(model,'FSurf','F5', ...
    struct('sel',data.IntZNodes(:,1), ...
    'def',repmat(r1,size(data.IntZNodes,1),1), ...
    'DOF',feutil('getdof',data.IntZNodes(:,1),(1:3)'/100)));
 model=fe_case(model,'FSurf','F6', ...
    struct('sel',data.IntZNodes(:,2), ...
    'def',repmat(-r1,size(data.IntZNodes,1),1), ...
    'DOF',feutil('getdof',data.IntZNodes(:,2),(1:3)'/100)));
end

RO.NeedHomog=1;
[z,RO,mo1,C1,Load,w,ft,carg]=safeInitDfrf([{CAM,model,RO},varargin(carg:end)],2);

RO.dd0=(d2.def'*mo1.K{RO.ihomog}*d2.def)/R2.V;% q_dK^T K_h q_dK = dd0*V 


if 1==2 % Check virtual work of homogeneous
 dd=reshape(C2.GroupInfo{4}(C2.GroupInfo{end}.ConstitTopology{1},1),6,6)
 norm(dd-(d2.def'*mo2.K{1}*d2.def)/R2.V)
end

K=feutilb('sumkcoef',mo1.K,[RO.alpha 1]);
% Compute qh
if isequal(Load.DOF,mo1.DOF) 
 def=struct('def',ofact({K,C1.T},Load.def),'DOF',mo1.DOF);
else
 def=struct('def',ofact(feutilb('tkt',C1.T,K),Load.def),'DOF',C1.DOF);
end
% qdS
out=def;out1=struct('dd',((d2.def'*mo1.K{strcmpi(mo1.Klab,'homog')}*def.def)/R2.V)\RO.dd0,'lab',{d2.lab});
%cf=feplot;cf.def=def; 

elseif comstr(Cam,'perdef')
%% #RvePerDef : periodic displacement

if 1==2
 % Signal periodic over Nx integer number of cells of the form
 % cos((n*dx+x0)/Nx*dx*2*pi)
 R2=struct('nx',51,'ny',500,'nz',500);
 cx=struct('X',{{(0:R2.nx-1)',(0:R2.ny-1)',(0:R2.nz-1)'}}, ...
    'Xlab',{{'ncx','ncy','ncz'}}, ... % ncx cell number x direction
    'Y',[]); 
 R2=struct('nx',51,'ikx',1:10);
 cx=struct('X',{{(0:R2.nx-1)'}},'Xlab',{{'ncx'}},'Y',[]);
 cx.Y=cos((cx.X{1}+1)/R2.nx*2*pi)/R2.nx*2;
 RO.Ck=fe_homo('dftAk',cx,R2);RO.Ck.Y(abs(RO.Ck.Y)<1e-12)=0;RO.Ck.Y
end


if carg>nargin &&isfield(varargin{carg},'ncx')
 %% Call to initialize RO 
 RO=varargin{carg};carg=carg+1;
 RO=struct('Range',RO,'bFcn',fe_homo('@RvePerDef'));
 RO.Post='def1.data(RunOpt.icol,4+(1:6))=RO.b''*u;def1.X{2}(5:10)={''vw''};';
 RO.Post='fe_homo(''RveHarmPos'')';
 RD=struct('mno',struct('m',0:9,'n',0),'cf',2);
 if isfield(RO.Range,'ncz');RD.mno.o=0;end
 out=RO;out1=RD;
 return;
end
[model,RO,SE,Case,Load,w,ft,carg]=safeInitDfrf(varargin,carg);
RO=dftu('RangeNc',RO);  
if ~isfield(RO,'oProp');RO.oProp={};end
if isempty(w);w=0;elseif length(w)>1;error('Not expected');end

RO.CellDir=getCellDir('vect',model);RO.RangeNc.CellDir=RO.CellDir;
RO.k=dftu('getx',struct('RangeNc',RO.RangeNc,'CellDir',RO.CellDir), ...
    struct('lab',{{'kx','ky','kz'}}));
out1=[];
RO.iVisco=find(strcmpi(RO.RangeNc.lab,'iVisco'));
SE1=SE;
for jpar=1:size(RO.RangeNc.val,1)
 %RO.cx=sprintf('exp(1i*x*%.15g)',2*pi/(RO.nx*data.CellDir(1)));
 % 2*pi/10/RO.CellDir(end)
 % Build Ud(kappa,x)
 
 if isempty(RO.iVisco) % Detect need to reassemble or not
   error('Need renew r2.Kh');
 elseif jpar==1 || diff(RO.RangeNc.val(jpar+[-1 0],RO.iVisco))
     sdtw('_ewt','detect example')
   r3=RO.RangeNc.param.iVisco.data{RO.RangeNc.val(jpar,RO.iVisco)};
   r2=feutilb('sumkcoef',SE.K,r3);
   SE1.K=struct2cell(r2)';SE1.Klab=fieldnames(r2)';
 else; r2=SE1.Klab(:)';r2(2,:)=SE1.K;r2=struct(r2{:});
 end
 if jpar==1
  out=struct('X',{{(1:size(RO.RangeNc.val,1))',{'u';'v';'w'}, ...
    {'Kh';'Km';'Ed'} }},'Xlab',{{'jpar','dir','comp'}}, ...
    'Y',zeros(size(RO.RangeNc.val,1),3,3));
 end
 st=sprintf('exp(1i*(x*%.15g+y*%.15g+z*%.15g))',RO.k(jpar,:));
 d1={'u',struct('dir',{{st,'0','0'}},'DOF',[.01;.02;.03])
     'v',struct('dir',{{'0',st,'0'}},'DOF',[.01;.02;.03])
     'w',struct('dir',{{'0','0',st}},'DOF',[.01;.02;.03])};
 ud=[];
 RB=dftu('clGen',struct,SE); 
 zCoef=stack_get(fe_def('zcoef-default;',SE),'info','zCoef','get');

 for j1=1:size(d1,1)
  vf1=elem0('VectFromDirAtDof',SE,d1{j1,2},SE.DOF);
  ud=fe_def('appenddef',ud,vf1); 
 end
 ud=feutil('rmfield',ud,'dir');ud.lab=d1(:,1); % missing data
 ud.data=[1:3;ones(1,3)]';ud.data(:,3)=jpar;ud.Xlab={'DOF',{'Index','d','jpar'}};
 if ~isfield(RO,'PerLoad');RO.PerLoad='mhqd';end
 switch lower(RO.PerLoad)
 case 'kmqd'% Generate load using homogeneous stiffness
   fp=ud; fp.def=r2.Km*ud.def; 
 case 'khqd'% Generate load using heterogeneous stiffness
  fp=ud; fp.def=r2.Kh*ud.def; %
 case 'mhqd'% Generate load using heterogeneous mass
  fp=ud; fp.def=r2.Mh*ud.def; 
 otherwise; error('Not a valid case');
 end
 RO.ed=diag(real(ud.def'*fp.def));%diag(real(Fd)'*real(fd.def)+imag(Fd)'*imag(fd.def));
 fp=fe_cyclic('defLong',fp);fp.name='kd*u';
 fp.Range=RO.RangeNc;fp.Range.val=RO.RangeNc.val(jpar,:);% save cur experiment param
 RB.nc=RO.RangeNc.val(jpar,1:3);
 RB.Load=fp; [K2,T,RB]=dftu('buildK',RB,SE1,Case);
 r2=SE1.Klab(:)';r2(2,:)=K2;r2=struct(r2{:});
 q=ofact(r2.Kh,T'*fp.def,RO.oProp{:});% This now the hetero response to hetero Fh
 RO.eh=feutilb('dtkt',q,r2.Kh); % Energy on hetero resp
 % real(|qd|_khetero) = RO.Ed ==> real(RO.eh)./RO.ed = 1 ==> OK :
 % field build with case RO.PerLoad='khqd' 
 q=ofact(r2.Km,T'*fp.def);% This now the homog response to hetero Fh
 RO.em=feutilb('dtkt',q,r2.Km); % the energy of homog response
 out.Y(jpar,:,:)=[RO.eh RO.em RO.ed];
 if nargout>1
   d1=struct('def',T*q,'DOF',fp.DOF,'data',[(1:3)' jpar*ones(size(q,2),1)], ...
     'Xlab',{{'DOF',{'dir';'jpar'}}});
   out1=fe_def('appenddef',out1,d1);
 end
 %dtkt=@(q,k)diag(real(q)'*k*real(q)+imag(q)'*k*imag(q));dtkt(q,K2{2})
end
out1.Range=RO.RangeNc;

elseif comstr(Cam,'harm')
%% #RveHarmPost : post-processing to have harmonic response

RO=[];RunOpt=[];SE=[];Case=[];def1=[];model=[];
eval(iigui({'model','RO','RunOpt','SE','Case','def1'},'GetInCaller')) 
L=sqrt(sum(RunOpt.CellDir.^2));

r1=RunOpt.nc;r1(r1~=0)=1./r1(r1~=0);
K=RunOpt.CellDir(:,1:length(r1))*r1(:);
u=def1.def(:,RunOpt.icol);
VW=u.'*RO.b*sum(RunOpt.nc.^2)/pi^2;
u=reshape(u,[],2); E=u'*model.K{end}*u;E=sum(diag(E));
loss=@(x)imag(x)./real(x);%loss(E)
def1.data(RunOpt.icol,4)=loss(E);def1.Xlab{2}{4,1}='Loss';
assignin('caller','def1',def1);
%dd=RO.dd{1}+RO.dd{2};r1=VW/dd(1);r1(abs(r1)<1e-10)=0;
%disp(r1)
return

elseif comstr(Cam,'khomo')
%% #RveKhomo : assemble the homogeneous material

model=varargin{carg};carg=carg+1;
Case=varargin{carg};carg=carg+1;
[Case,model.DOF]=fe_mknl('init',fe_case(model,'reset'));
RO=varargin{carg};carg=carg+1;
r1=Case.GroupInfo(:,4);RO.K_homo=cell(1,2);
for j1=1:2
 for jgroup=1:size(Case.GroupInfo,1)
  % For now only deal with uniform material
  constit=r1{jgroup}(:,1);EC=Case.GroupInfo{jgroup,8};
  dd=constit(EC.ConstitTopology{1}); 
  if j1==1;dd(4:6,:)=0;dd(:,4:6)=0;
  else; dd(1:3,:)=0;dd(:,1:3)=0;
  end
  constit(3:38)=dd(:);RO.dd{j1}=dd;
  Case.GroupInfo{jgroup,4}=repmat(constit,1,size(Case.GroupInfo{jgroup,4},2));
 end
 model.name='';
 RO.K_homo{j1}=fe_mknl('assemble NoT',model,Case,1);
end
out=RO;

elseif comstr(Cam,'displc')
%% #RveDispLC : dispersion curve to mesh convergence
model=varargin{carg};carg=carg+1;
RO=varargin{carg};carg=carg+1;
RO=sdth.sfield('addmissing',RO,struct('nc',[5:10 20:10:100]','gf',1));
gf=RO.gf; figure(gf);clf; 
%  compression speed sqrt(model.pl(3)/model.pl(5))
% E=model.pl(3);rho=model.pl(5);nu=model.pl(4); sqrt(E*(1-nu)/(1+nu)/(1-2*nu)/rho)

RO.rel=NaN;
for j1=1:length(RO.Lc)
 mo1=feutil(sprintf('objecthexa %i %i',model.pl(1),model.il(1)),[0 0 0;eye(3)], ...
     [0 RO.Lc(j1)],[0 RO.Lc(j1)],[0 RO.Lc(j1)]);
 if isfield(RO,'quad')&&RO.quad; mo1=feutil('lin2quad',mo1);end
 mo1.pl=model.pl;mo1.il=model.il;
 mo1=fe_case(mo1,'fixdof','z','z==0 -DOF3','fixdof','y','y==0 -DOF2');
 if isfield(RO,'epsl');st=sprintf('dftbuild epsl%g;',RO.epsl);else; st='dftbuild;';end
 mo1=fe_homo(st,mo1,RO.Lc(j1)*[1 1 1]);
 mo1=stack_set(mo1,'info','EigOpt',[5 10 0]);
 mo1.name='';
 RT=sdth.sfield('addselected',struct('Range',struct('ncx',RO.nc,'ncy',1,'ncz',1)),RO,'MatDes');
 [def,hist]=fe_homo('dftDisp',mo1,RT); hist.Y=hist.Y(:,1);
 r1=fe_homo('dftVg',hist,struct('unit',model.unit,'ind',1,'hold','on','cf',gf,'rel',RO.rel));
 RO.rel=r1.rel;
 %hold on; plot(hist.Y(:,1)/1e3,pi*hist.X{1},'-+');hold off;
 %xlabel('Omega (kHz)');ylabel(hist.Xlab{1});
end
cingui('objset',gf,{'@OsDic',{'ImGrid','ImSw40'},'@axes',{'xscale','log'},'@line',{'marker','.'}})
st=cellfun(@(x)comstr(x,-30),num2cell(RO.Lc),'uni',0);
legend(st,'location','best');
if isfield(RO,'Report')&&RO.Report
 cingui('plotwd',gf,'@OsDic(SDT Root)',{'WrW49c'});comgui('imwrite',gf);
end

else;error('Rve%s unknown',CAM);
end

elseif comstr(Cam,'ass');[CAM,Cam]=comstr(CAM,4);
%% #Ass generic utilities associated with assembly

if comstr(Cam,'gaussobs')
%% #AssGaussObs : strain at gauss observation fe_homo('AssGaussObs',model)
mo1=varargin{carg};carg=carg+1;
if carg<=nargin; st=varargin{carg};carg=carg+1;
 mo1.Elt=feutil(['selelt innode ' st],mo1);
end  

% post-process all stress components
% Assume one material per group
mpid=unique(feutil('mpid',mo1),'rows'); 
if size(mpid,1)-1~=nnz(~isfinite(mo1.Elt))
   warning(['Using feutilb(''SeparatebyProp'') to divide in groups ' ...
    'with unique MatId/ProId']);
end
[eltid,mo1.Elt]=feutil('eltidfix;',mo1);
sel=fe_caseg('stresscut -selout',struct('type','Gauss'),mo1); 
obs=sel.StressObs;
obs.X(2:3)={obs.Lambda{1,4}.w,eltid(eltid~=0)};
obs.Xlab={'StrainComp','Gauss','EltId'};
obs.cta=getConst('ild',obs)*obs.cta;
out=obs;

else; error('Ass%s Unknown',CAM);
end


%% #end ------------------------------------------------------------------- -2
elseif comstr(Cam,'@');out=eval(CAM);
elseif comstr(Cam,'cvs')
 out=sdtcheck('revision');
 %out='$Revision: 530 $  $Date: 2020-12-11 11:10:34 +0100 (Fri, 11 Dec 2020) $'; 
else; error('%s unknown',CAM);
end

%% #subFunc
%% #getConst Build assembled constitutive law
function out=getConst(CAM,obs)

[CAM,Cam]=comstr(CAM,1);

if comstr(Cam,'ild')  % strain = ILD * stress
 i1=getConst('MatPos',obs); 
 % Get back to strains
 ld=obs.Lambda(:,1);  for j1=1:size(ld,1);ld{j1}=inv(ld{j1});end
 ld=reshape(ld(i1),1,[]);ld=horzcat(ld{:}); % All constitutive
 II=repmat((1:size(ld,1))',1,size(ld,2))+ ... % Row in Lambda
  reshape(repmat((0:length(i1)-1)*size(ld,1),size(ld,1)^2,1),size(ld,1),[]);
 JJ=repmat(repmat(1:size(ld,1),size(ld,1),1),1,length(i1)) + ... %col in Ld
  reshape(repmat((0:length(i1)-1)*size(ld,1),size(ld,1)^2,1),size(ld,1),[]);
 out=sparse(II,JJ,ld);
elseif comstr(Cam,'ldj')  % wjdet * LD symmetric
 
 if any(obs.wjdet<=0);error('Non positive jdet');end
 out=getConst(['ld' CAM(4:end)],obs);
 [CAM,Cam,i2]=comstr('mat',[-25 1],CAM,Cam); if isempty(i2);i2=1;end
 if i2==1;i2=size(obs.X{1},1);elseif i2==2;i2=size(obs.trans.X{1},1);
 else; error('Not an expected case');
 end
 jdet=diag(sparse(reshape(repmat(sqrt(obs.wjdet'),i2,1),[],1)));
 out=jdet*out*jdet; % Wait

elseif comstr(Cam,'ld')  % stress = LD * strain, -mat i to specify matrix
    
 [CAM,Cam,i2]=comstr('mat',[-25 1],CAM,Cam); if isempty(i2);i2=1;end
 ld=obs.Lambda(:,i2);  i1=getConst('MatPos',obs); 
 ld=reshape(ld(i1),1,[]);ld=horzcat(ld{:}); % All constitutive
 II=repmat((1:size(ld,1))',1,size(ld,2))+ ... % Row in Lambda
  reshape(repmat((0:length(i1)-1)*size(ld,1),size(ld,1)^2,1),size(ld,1),[]);
 JJ=repmat(repmat(1:size(ld,1),size(ld,1),1),1,length(i1)) + ... %col in Ld
  reshape(repmat((0:length(i1)-1)*size(ld,1),size(ld,1)^2,1),size(ld,1),[]);
 out=sparse(II,JJ,ld);
 
elseif comstr(Cam,'matpos')  % strain = ILD * stress
 i2=cellfun(@(x)x.ID(1),obs.Lambda(:,end)); % Get the MatID
 i1=feutil('proid',feutil('selelt eltname node',obs));
 i1(i1==0)=[]; % ProId associated with gausspoint
 [i1,out]=ismember(i1,i2);
else;error('%s unknown',CAM);
end  



%% #DftRest : curvemodel implementation of spatial restitution
function out=DftRest(varargin) 
 if ischar(varargin{1});[CAM,Cam]=comstr(varargin{1},1);carg=2;end
 persistent omethods
 if isempty(omethods); omethods=fe_def('@omethods');end
 
 if comstr(Cam,'getx')
 %% #getX -3
   RO=varargin{2}.Source;ind=varargin{3};typ=varargin{4};
   if isfield(RO,'X')
    if ind==1;out=RO.X{1};
    elseif ind==2&&~isequal(RO.Xlab{2},'DOF');out=RO.X{2}; % k
    elseif isempty(RO.X{ind}); out=RO.DOF(RO.CurInd);
    else; out=RO.X{ind};
    end
    if ~comstr(Cam,'getxstring') 
    elseif iscell(out);out=out{typ};% Called in iiplot
    else; out=sprintf('%.2f',out(typ));
    end
   else % Case for feplot (DOFs)
    if ind==1;out=RO.Rest.DOF; % DOF
    elseif ind==2;out=RO.data(:,1);  % channels
    end
   end
 elseif comstr(Cam,'()'); 
  %% #GetY Generate Y and return values -3
  
  obj=varargin{2};def=obj.Source;S=varargin{3};
  r2=def.Rest;
  if length(S)>1;error('Not expected');end
  ch=S.subs{2};
  if isfield(r2,'type')&&strncmpi(r2.type,'dfrf',4)
  %% #Multi-harmonic case Source.Rest.type='dfrf'; -3
  if length(ch)>1; error('Expecting single channel');
  elseif ischar(ch)&&isequal(ch,':') % multiple ch restitution
   r1=zeros(size(def.Rest,1),size(def.data,1));
   for ch=1:size(def.data,1)
     S1=S;S1.subs{2}=ch; r2=DftRest('()',obj,S1);
     S1.subs{2}=':';r2=subsref(r2,S1);
     r1(:,ch)=r2;
   end
   out=r1; return
  end
  
  if strcmpi(r2.type,'dfrf')
  %% #Dfrf case with [k x omega storage] -4
   if ~isfield(r2,'Ck')||isempty(r2.Ck)
    ind=(ch-1)*length(r2.kcx)+(1:length(r2.kcx));
   else
    ind=1:size(r2.Ck.Y,1); ind=(ch-1)*ind(end)+ind;
   end
   r1=def.def(:,ind); 
  else % #Dfrfb if strcmpi(r2.type,'dfrfrb') % case [omega x k] -4
   ind=ch:size(def.data,1):size(def.def,2);
   r1=def.def(:,ind); 
  end
  r3=dftu('enk',r2.mno,r2.kcx);r3=r3.Enk;% May have problems with mno cols
  if size(r1,1)==size(r2.cna,2)*2
   error('Not reimplemented')
   %r1=r2.cna*reshape(r1,size(r2.cna,2),[]); % (point,dir) x (re im)
  elseif size(r1,1)==size(r2.cna,2)
   r1=reshape(r2.cna*r1,size(r2.cna,1)/2,[]);
   %r1=r2.cna*r1;r1=reshape([real(r1);imag(r1)],size(r1,1),[]);
   r1=r1*r3';
  else;  r1=r1*r3';

  end
  % Back-transform of monoharmonic
  % Re(U(k)e(j (kx m+ky n + kz o))
  % Cnk (for cells in ind), 20 allows negative mno for 20 blocks
  %r3=r2.Ck.Enk(remi(r2.mno+1+size(r2.Ck.Enk,1)*20,size(r2.Ck.Enk,1)),:); 
  
  % xxx deal with complex feplot scaling
  % Step3 permute 
  r1=reshape(r1,[],3,size(r1,2)); %dof,dir,cell
  r1=permute(r1,[1 3 2]);r1=reshape(r1,[],1); % Give column
  out=r1;
  
  elseif ~isfield(def,'X')
  %% #Mono-harmonic case, single channel for feplot -3
  % step 1 : dvert(@nvert0*(x,y,z),ik) : p,d,k
  % step 2 : d_vert * E = dvert(@vert0*(x,y,z),n) : p,d,n
  % step 3 : permute to obtain dvert([@vx,n @vy,n @vz,n] 
  r1=def.def(:,ch); 
  if size(r1,2)>1; error('Expecting single channel')
  elseif size(r1,1)==size(r2.cna,2)*2
   error('Not reimplemented')
   %r1=r2.cna*reshape(r1,[],2); % (point,dir) x (re im)
  elseif size(r1,1)==size(r2.cna,2)
   r1=reshape(r2.cna*r1,[],2);
  end
  r3=dftu('getx',def,struct('lab',{{'kcx','kcy','kcz'}}),ch);
  if isscalar(r3); Enk=exp(1j*r2.mno(:,1)*r3(1)); % single dir for now
  else;
     Enk=exp(1j*r2.mno*reshape(r3,size(r2.mno,2),1));
  end
  % Back-transform of monoharmonic
  % Re(U(k)e(j (kx m+ky n + kz o))
  Enk=[real(Enk) -imag(Enk)];
  if size(r1,2)==1; error('Expecting different form');end
  r1=r1*Enk';

  % xxx deal with complex feplot scaling
  % Step3 permute 
  r1=reshape(r1,[],3,size(r1,2)); %dof,dir,cell
  r1=permute(r1,[1 3 2]);r1=reshape(r1,[],1); % Give column
  out=r1;
  elseif strncmp(def.Rest.type,'cfrfb',4)
  %% #curve_on_first cell with recombined harmonics
   i1=[def.CurInd(:);size(def.def,1)/2+def.CurInd(:)];% long format
   if def.Rest.type(end)=='b' % [(dof x freq) x (r/i * k)
     out=reshape(def.def(i1,:),length(def.CurInd),2,[],size(def.Rest.kcx,1));
     out=reshape(permute(out,[3 1 2 4]),[],2*size(def.Rest.kcx,1)); 
   else % Sort with (k*freq)
    out=reshape(def.def(i1,:),length(def.CurInd),2,size(def.Rest.kcx,1),[]);
    out=reshape(permute(out,[4 1 2 3]),[],2*size(def.Rest.kcx,1));     
   end
   r3=dftu('enk',def.Rest.mno(:,1),def.Rest.kcx);
   out=reshape(out*r3.Enk',[],length(def.CurInd),size(r3.Enk,1)); 
   out=subsref(out,S); % freq,DOF
   
  else
  %% #curve_on_first cell with specific harmonic
  
   out=reshape(def.def(def.CurInd,:),[],size(def.X{1},1),size(def.X{2},1));
   out=permute(out,[2 3 1]);
   out=subsref(out,S); % freq,kappa,DOF
   if length(S.subs)==3&&isscalar(S.subs{2})&&isfield(def,'Disp') 
    % Single kappa
    try
     i1=S.subs{2}; 
     if ischar(i1)
     else; obj.Source.ID=struct('po',def.Disp.Y(i1,:)');
     end
    catch; obj.Source.ID=[];
    end
   end
   
  end
 elseif comstr(Cam,'labfcn'); % generate LabFcn
   def=varargin{carg};carg=carg+1;
   ch=varargin{carg};carg=carg+1;
   i2=find(strcmpi(def.Xlab{2},'jpar'));
   st=fe_range('labdef',def,ch);
   r2=def.Xlab{2};
   if iscell(r2)&&strcmpi(r2{1},'kx')&&i2==2 % clean kx
    out=st;
   elseif size(def.data,2)==1 % Just frequency
    out=sprintf('%.1f Hz',def.data(ch,1));
   elseif isfield(def,'Range')&&isfield(def.Range,'param') && ...
       isfield(def.Range.param,'Freq')
       out=st;
   else
    if ischar(st);st={st};end
    out=sprintf('%.1f Hz load %i, %s',def.data(ch,1:2),st{1});
   end
 elseif comstr(Cam,'xlab');
     obj=varargin{2};def=obj.Source;
     if isfield(def,'Xlab')
      out=def.Xlab;
     else
        dbstack;
      out={'Time [s]','channel'};
     end
 else; error('% unknown',Cam);
 end
%% #ktypeCheck 
function  [RO,def]=ktypeCheck(RO,def);

if ~isfield(RO,'ktype');
 try; RO.ktype=def.Range.lab{1};
 catch;RO.ktype='kcx';
 end
end
if isfield(RO,'hist');hist=RO.hist; else;hist=[];end 
if isfield(def,'Range')&&isfield(def.Range,'CellDir'); 
    RO.CellDir=def.Range.CellDir;
else; RO.CellDir=RO.hist.Range.CellDir;
end
if length(RO.CellDir)>1;RO.CellDir=norm(RO.CellDir);end
%% edit the def range

if ~isfield(RO,'unit')||strncmpi(RO.unit,'us',2);
  r1=fe_mat(sprintf('convertSIUS')); r1=sdtm.toStruct(r1(:,[4 4]));
else
 r1=fe_mat(sprintf('convert%sUS',RO.unit(1:2))); r1=sdtm.toStruct(r1(:,[4 2]));
end
r1={'ncxkcx',@(x)2*pi./x,'kcx','rad/cell';
    'ncxkc',@(x)1./x,'kc','1/cell';
   'ncxkx',@(x)2*pi/RO.CellDir./x,'kx',['rad/' r1.length]
   'kxncx',@(x)2*pi./x./RO.CellDir,'ncx','cell'
   'kxlx',@(x)2*pi./x,'lx',r1.length
   'lxkx',@(x)2*pi./x,'kx','rad/len xxx'
   'lxkcx',@(x)2*pi*RO.CellDir./x,'kcx','rad/cell xxx'
   'kxkcx',@(x)x*RO.CellDir,'kcx','rad/cell';
   'kxkc',@(x)x*RO.CellDir/(2*pi),'kc','1/cell';
   };
for j1=1:2
 if j1==1; if isempty(def);continue;end;st=def.Range.lab{1};
 elseif isempty(hist); break;
 else; st=hist.Xlab{1};
 end
 if iscell(st);st=st{1};end
 if strcmpi(RO.ktype,st) % no change
 else
  if strncmpi(st,'lambda',6);st='lx';end
  r2=r1(strcmpi(r1(:,1),[st,RO.ktype]),:);
  if isempty(r2);
      error('%s not in %s',[st,RO.ktype],sdtm.toString(r1(:,1)'));
  elseif j1==1; RO.defKtype=r2(2:4);
     if ~isempty(def);
         def.Range.lab{1}=RO.ktype;
         def.Range.val(:,1)=r2{2}(def.Range.val(:,1));
     end
  elseif j1==2; 
    hist.Xlab{1}={r2{3:4} []};
    hist.X{1}=r2{2}(hist.X{1});
  end
 end
end
if ~isempty(hist);RO.hist=hist;end

 
%% #DFTu utilitiles
function [out,out1,out2]=dftu(CAM,varargin)

[CAM,Cam]=comstr(CAM,1);
if comstr(Cam,'rangenc')
 %% #RangeNC extract the Nc range -3
 % dftu('RangeNc',RO,Range)
   RO=varargin{1};
   if nargin>2;Range=varargin{2};
   elseif isfield(RO,'Range');Range=RO.Range;
   else % Use amplitude curves
    if ~isfield(RO,'a');r1=1;
    else; r1=RO.a;
    end
    st=r1.Xlab{1}; if ~iscell(st);st={st};end
    i1=strcmpi(st,'kx');
    if any(i1);
      r2=r1.X{1}(:,i1);r2(r2~=0)=1./r2(r2~=0);
      Range=struct('val',[r2 r1.Y],'lab',{{'inx','ak'}});
    else; error('not implemented')
    end
   end
   if ~isfield(Range,'lab')
    Range=feutil('rmfield',fe_range('buildgrid',Range),'level','type','VectFcn','Ngrid');
   end
   [i1,i2]=ismember({'ncx','ncy','ncz','kx','ky','kz', ...
       'kcx','kcy','kcz','inx','iny','inz','ak','iVisco','zCoef'},Range.lab);
   RO.i_nc=i2(i1); % Colums used for wave number spec
   r2=Range;r2.val=r2.val(:,i2(i1));r2.lab=r2.lab(i2(i1));
   [i1,i2,i3]=unique(r2.val,'rows','first');r2.val=r2.val(sort(i2),:);
   if ~any(strcmpi(r2.lab,'ak'));r2.val(:,end+1)=1;r2.lab{end+1}='ak';end
   if isfield(RO,'cyc'); % Propagate periodicity information
    if isfield(RO.cyc,'trans');  RO.cyc.CellDir=RO.cyc.trans(:);
    else;r2.CellDir=RO.cyc.CellDir; 
    end
   end
   RO.RangeNc=r2;
   out=RO;
   
elseif comstr(Cam,'rangefreq')
 %% #RangeFreq build -3
   RO=varargin{1};
   if nargin==3;Range=varargin{2};else;Range=RO.Range;end
   
   [i1,i2]=ismember({'Freq'},Range.lab);
   RO.i_freq=i2(i1);
   r2=Range;r2.val=r2.val(:,i2(i1));r2.lab=r2.lab(i2(i1));
   [i1,i2,i3]=unique(r2.val,'rows','first');r2.val=r2.val(sort(i2),:);
   RO.RangeFreq=r2;
   out=RO;
   
elseif comstr(Cam,'clgen')
 %% #Clgen build RO=feval(fe_homo('@dftu'),'clGen',RO,model) adds cl cr to RO-3
 
 RO=varargin{1};model=varargin{2};
    
 if ~isfield(RO,'cLx')
  if ~isfield(model,'DOF');model.DOF=feutil('getdof',model);end
  if ~isfield(RO,'DOF');RO.DOF=model.DOF;end
  r2=fe_cyclic('solveclcr',model,RO);
  if isfield(r2,'cL');
     r2.cLx=r2.cL;r2.cRx=r2.cR; % Robust 1D
     r2=feutil('rmfield',r2,'cL','cR');
  end
  st=fieldnames(r2);
  for j1=1:length(st);RO.(st{j1})=r2.(st{j1});end
 end
 if isfield(RO,'DirDofInd') 
  %% DirDofInd={1:4,1:4,[1 2 0 0]}; eliminate certain constraints
  for j1=1:length(RO.DirDofInd)
    st=char(abs('x')-1+j1);
    r2=RO.(['cL',st]);
    i2=RO.DirDofInd{j1}(:);i2=i2(:);
    i2=repmat(i2,size(r2,1)/length(i2),1);i2=(i2~=0);
    r2=r2(i2,:);RO.(['cL',st])=r2;
    RO.(['cR',st])=RO.(['cR',st])(i2,:);
  end
 end
 RO.CellDir=getCellDir('vect',model);
 out=RO;
 
elseif comstr(Cam,'buildk')
 %% #BuildK that verifies periodicity conditions -3
 RO=varargin{1};model=varargin{2};Case=varargin{3};
 
 if ischar(model); st=model;model=evalin('caller',st);evalin('caller',['clear ' st]);end
 data=fe_case(model,'getdata','Symmetry');
 if data.opt(1)~=-1; error('Not a periodic problem');end
 if ~isfield(RO,'nc') % Build Nc from range
   [RO.nc,r2]=dftu('getx',struct('Range',RO),struct('lab',{{'ncx','ncy','ncz'}}));
 end
 if isempty(RO.nc);error('Empty .nc not a valid case');end
 data.opt(1+(1:length(RO.nc)))=RO.nc; RO.dataOpt=data.opt;
 c=cell(6,1);  RO.isReal=1;
 for j1=1:length(data.opt)-1 % Loop on xyz
  Ld=data.opt(j1+1);
  if Ld==0; Ld=1; cma=1; sma=0; 
  else;    %2pi/(cell_periodicity) : exp(k x) with k=2pi/(Nx Lc) 
    ma=2*pi/Ld; cma=cos(ma); sma=sin(ma);    
  end
  cL=RO.(['cL' char(119+j1)]);cR=RO.(['cR' char(119+j1)]);
  z=spalloc(size(cL,1),size(cL,2),0);
  if Ld==2  % 
   c(2*j1+[-1 0])={[cL+cR z];[z cL+cR]};
  elseif Ld==1 % Simple periodicity Re(qL-qR)=0  Im(qL-qR)=0 
   % fecom('shownodemark',model.DOF(any(RO.cRz)))
   c(2*j1+[-1 0])={[cL-cR z];[z cL-cR]};
   if isfield(RO,'FixReal')&&RO.FixReal % No trans for real direction
     i1=find(any(cL));i1=fix(model.DOF(i1(1)));
     c{2*j1-1}=[c{2*j1-1};fe_c(model.DOF,i1+j1/100) z(1,:)];
     c{2*j1}=[c{2*j1};z(1,:) fe_c(model.DOF,i1+j1/100)];
   end
  else % Non real case
   RO.isReal=0;
   c(2*j1+[-1 0])={ % verified 2013
       [cL-cR*cma -sma*cR];
       [sma*cR cL-cma*cR]};
  end
 end
 % Fix one node since rigid body modes are periodic
 % xxx does not work
 if isfield(RO,'FixNode')
   r2=fe_c(model.DOF,RO.FixNode);
   c{7}=[r2 z(1:size(r2,1),:)];
   c{8}=[z(1:size(r2,1),:) r2];
 end
 K2=model.K;
 %{feval(fe_reduc('@get_mat'),model.K,2,model,struct) ...
 %      feval(fe_reduc('@get_mat'),model.K,1,model,struct)};
 if RO.isReal
  c=vertcat(c{1:2:end,1});c(:,size(c,2)/2+1:end)=[];
  T=Case.T*fe_coor(c*Case.T); 
  K2=feutilb('tkt',T,K2);
  if isfield(RO,'bFcn') % [Real;imag(B*ak)]*(complex(u(omega))
    dbstack; keyboard;
  elseif isfield(RO,'b');
    if ~isreal(RO.ak);error('ak should be real but u can be complex');end
    RO.b=real(RO.b*RO.ak)*RO.u;
  end
 else
  z=spalloc(size(model.K{1},1),size(model.K{1},2),1);
  for j1=1:length(K2);K2{j1}=[K2{j1} z;z K2{j1}];end
  c=vertcat(c{:}); 
  z=spalloc(size(Case.T,1),size(Case.T,2),0);T=[Case.T z;z Case.T];
  if ~isfield(RO,'coor');RO.coor='lut';end
  if strcmpi(RO.coor,'1')
   T2=T*fe_coor(c*T,[1 1]);
  elseif strcmpi(RO.coor,'3')
   T=T*fe_coor(c*T,[3 3]); % separate columns, LU
  else % if strcmpi(RO.coor,'lut');%exist('luq','file');
   %T=fe_coor(RO.coor,c*T,'T'); % replace t
   T1=fe_coor(RO.coor,struct('c',c,'T',T,'cDOF',[Case.DOF;Case.DOF+.5]));
   T=T1.T;
  end
  if 1==2%isfield(RO,'mpcond')% Scale model DOF
   i1=fe_c(model.DOF,RO.mpcond{1},'ind');T(i1,:)=T(i1,:)*RO.mpcond{2};
  end
  if isfield(RO,'bFcn') % [Real;imag(B*ak)]*(complex(u(omega))
    RO=feval(RO.bFcn,model,RO,data);
  elseif isfield(RO,'b') % [Real;imag(B*ak)]*(complex(u(omega))
    b=RO.b*RO.ak;
    RO.b=[real(b);imag(b)]*RO.u; 
  end
  clear model 
  K2=feutilb('tkt',T,'K2');
 end
 clear c
 out=K2; out1=T;out2=RO;
 
elseif comstr(Cam,'getx');
 %% #dftu_GetX (find nc,ld,Kx)  -3
 % feval(fe_homo('@dftu'),'getx lx',def,RO,ch)
 %  RO=struct('lab',{{'ak'}})
 def=varargin{1}; if nargin>2; RO=varargin{2};else; RO=struct;end
 if isfield(def,'Range');  Range=def.Range;
   RO.CellDir=getCellDir('vect',Range);
 elseif isfield(def,'RangeNc'); Range=def.RangeNc;
   RO.CellDir=getCellDir('vect',def);
 else; error('GetX implemented for range'); 
 end
 
 r1=Range.val(:,strncmpi(Range.lab,'ncx',3));
 if isfield(RO,'lab'); RO.Needed=RO.lab; else;RO.Needed=Range.lab;end
 if isfield(Range,'lc')&&isscalar(Range.lc)  % Monodim
    error('Update EB');
    if isempty(r1);
     r1=Range.val(:,strncmpi(Range.lab,'lambda',5))/Range.lc;
    end
    out=struct('X',[r1 Range.lc*r1 2*pi/Range.lc./r1], ...
        'Xlab',{{'nc';'Lambda';'Kx'}});
    if isfield(Range,'unit')
      r2=fe_mat(['convertSI',Range.unit]);r2=r2{4,4};
      out.Xlab={'nc',sprintf('Lambda [%s]',r2),sprintf('Kx [rd/%s]',r2)};
    end
 else; % By default define nc.
   ind=find(~cellfun('isempty',regexp(Range.lab,'nc[xyz]')));
   out=struct('X',Range.val(:,ind),'Xlab',{Range.lab(ind)});
   RO.nc=out.X;
   if isempty(RO.nc)&&size(RO.CellDir,2)==1 
    % single dir with just kx
    i1=strncmpi(Range.lab,'kx',2);
    if any(i1); 
      r1=Range.val(:,i1); % kx=2*pi/kcx=2*pi*(dx*kx)
      RO.nc=r1;RO.nc(r1~=0)=2*pi/norm(RO.CellDir)./r1(r1~=0);
      RO.nc(r1==0)=1; % one periodicity for
    end
   end
 end
 st=setdiff(RO.Needed,out.Xlab);st=setdiff(st,{'ncy','ncz'});
 if ~isempty(st) % Fill any desired vector
  if size(RO.CellDir,1)==1; % r1 cols contains periodicity directions
    r1=diag(RO.CellDir);if size(r1,1)<3;r1(3,1)=0; end
  else;r1=RO.CellDir;
  end
  if ~isempty(RO.nc);
   r2=RO.nc; r2(r2==0)=1; % nc=0 same as nc=1 or infinity
   r3=2*pi./(r2*diag((sum(r1(:,1:size(r2,2)).^2,1)))); % 2*pi/(ncx |dx|^2)
   r3(r2==1)=0; % wave number of 2*pi should be zero;
   if size(RO.CellDir,2)==1;  
     RO.k_xyz=r3*norm(r1); % Single dir
   else; RO.k_xyz=r3*r1(:,1:size(r3,2))'; % kappa xyz
   end
  else; RO.k_xyz=[];
   i2=find(strcmpi(Range.lab,'kcx')); % Given in kx
   if ~isempty(i2); RO.k_xyz=Range.val(:,i2)/norm(r1(:,1));end 
  end
  for j1=1:length(st)
   st1=sscanf(st{j1},'%s',1);
   switch lower(st1) % sdtweb nida13('dispk')
   case 'lambda'; 
       r1=sqrt(sum(RO.k_xyz.^2,2)); % amplitude of wave number rad/l
       r1(r1~=0)=2*pi./r1(r1~=0); % lambda [l]
   case 'kx';  % physical wave number in x direction
      r1=RO.k_xyz(:,1);
   case 'ky';  r1=RO.k_xyz(:,2); % physical wave number in y direction
   case 'kz';  r1=RO.k_xyz(:,3); % physical wave number in z direction
   case 'kappa';r1=sqrt(sum(RO.k_xyz.^2,2)); % amplitude of wave vector rad/l
   % Ku,Kv,Kw    
   case 'kcx';  % cell wave number in x direction rad
     if isempty(RO.nc)
       r1=Range.val(:,strncmpi(Range.lab,'kappa',5))*norm(RO.CellDir(:,1));
     else; 
       r2=RO.nc(:,1); r2(r2==0)=1; r1=2*pi./r2;r1(r2==1)=0; 
     end
   case 'kcy';  % cell wave number in y direction rad
    if size(RO.CellDir,2)<2;% Possibly ignore (no periodicity in y)
      RO.lab(strcmpi(RO.lab,st{j1}))=[];st{j1}='';
      continue;
    elseif size(RO.nc,2)>1
       r2=RO.nc(:,2); r2(r2==0)=1; r1=2*pi./r2;r1(r2==1)=0; 
    else
     error('Not implemented');
    end 
   case 'kcz';  % cell wave number in y direction rad
    if size(RO.CellDir,2)<3||size(RO.nc,2)<3;% Possibly ignore
      RO.lab(strcmpi(RO.lab,st{j1}))=[];st{j1}='';
      continue;
    elseif size(RO.nc,2)>2
       r2=RO.nc(:,3); r2(r2==0)=1; r1=2*pi./r2;r1(r2==1)=0; 
    else;    error('Not implemented');
    end 
   case 'ak';  % reuse value in Range
     r1=Range.val(:,strcmpi(Range.lab,'ak'));
   case 'ncx';  % find ncx from kcx
     r1=Range.val(:,strcmpi(Range.lab,'kcx'));
     if isempty(r1);% xxx needs checking
      r1=Range.val(:,strcmpi(Range.lab,'kx'));% physical wave kx=2pi/(ncx*dx)
      i1=(r1~=0); r1(i1)=2*pi./r1(i1)/norm(RO.CellDir(:,1));%ncx=2*pi/kx/dx
      
     else;
       i1=(r1~=0); r1(i1)=2*pi./r1(i1);
       r1(~i1)=1; %ncx=1 for kcx=0 (simple periodicity) 
     end
   otherwise;
     dbstack; keyboard
   end
   if isempty(r1); % Possibly reuse value given in range
     i1=find(strncmpi(Range.lab,st1,length(st1)));
     if ~isempty(i1);r1=Range.val(:,i1);end
   end
   out.X(:,end+1)=r1;out.Xlab{end+1}=st{j1};
  end
  st(cellfun('isempty',st))=[]; % y,z possibly ignored
 end
 if isfield(RO,'lab') % Return specific units
     ind=zeros(1,length(RO.lab));
     for j1=1:length(ind);
         i2=find(strncmpi(out.Xlab,RO.lab{j1},length(RO.lab{j1})));
         if ~isempty(i2);ind(j1)=i2(1);end
     end
     i2=(ind==0&strncmpi(RO.lab,'nc',2));
     ind(i2)=[]; % ignore nc. not present (used to detect ncy,ncz)
     out1=out.Xlab(ind);out=out.X(:,ind);
 end
 if nargin>3; % Allow ch selection based on jpar value
   ch=varargin{3};
   i1=strcmpi(def.Xlab{2},'jpar');
   if any(i1);ch=def.data(ch,i1);end
   out=out(ch,:);
 end 
elseif comstr(Cam,'addrange')
 %% #AddRange feval(fe_homo('@dftu'),'addrange',def,{'kx'}) -3
 def=varargin{1};
 st=varargin{2};if ~iscell(st);st={st};end
 for j1=1:length(st)
  st1=sscanf(st{j1},'%s',1);
  i1=strncmpi(def.Range.lab,st1,length(st1));
  if ~any(i1);def.Range.lab{end+1}=st{j1};end
 end
 r1=dftu('getx',def);
 
 for j1=1:length(r1.Xlab)
  i1=strcmpi(def.Range.lab,r1.Xlab{j1});
  def.Range.val(:,i1)=r1.X(:,j1);
 end
 out=def;
 
 
elseif comstr(Cam,'enk')
 %% #Enk feval(fe_homo('@dftu'),'enk',inx,kcx) -3
 inx=varargin{1}(:); kcx=varargin{2}(:);
 dk=diff([0;kcx])+diff([kcx;kcx(end)]);
 if pi-kcx(end)>=dk(end)/2; dk(end)=dk(end)*2;end
 
 % if abs(kcx(end)-pi)>1e-5;dk(end)=dk(end-1);end % No mid spectrum adjust for DFT 
 out=struct('inx',inx,'kcx',kcx,'Enk', ...
   reshape([cos(inx*kcx');-sin(inx*kcx')]*diag(dk/2/pi), ...
      length(inx),[])); % ' xxx delta x'
  
 if 1==2 % Manual verification
  N=24;kcx=(0:floor(N/2))'/N*2*pi;
  inx=(0:N-1)';
  dk=diff([0;kcx])+diff([kcx;kcx(end)]);
  if pi-kcx(end)>=dk(end)/2; dk(end)=dk(end)*2;end
  %dk=diff([0;kcx])+diff([kcx;pi]);
  %if abs(kcx(end)-pi)>1e-5;dk(end)=dk(end-1);end % No mid spectrum
  RO.Enk=reshape([cos(inx*kcx');-sin(inx*kcx')]*diag(dk),size(inx,1),[])/2/pi;
  RO.Ekn=reshape([cos(inx*kcx');-sin(inx*kcx')],size(inx,1),[])';
  r2=fe_cyclic('reftransform',N)/N;norm(r2-RO.Enk(:,[1 3:size(r2,2)+1]))
  
  % Now check with part spectrum
  ik=1:N/2;ik(5:2:end)=[]; 
  kcx=kcx(ik); q0=RO.Enk;%reshape(RO.Enk,size(RO.Enk,1)*2,[]); % reference shapes
  dk=diff([0;kcx])+diff([kcx;pi]);
  RO.Enk=reshape([cos(inx*kcx');-sin(inx*kcx')]*diag(dk),size(inx,1),[])/2/pi;
  RO.Enk*RO.Ekn;
 end

elseif comstr(Cam,'tospace')
 %% #ToSpace feval(fe_homo('@dftu'),'toSpace',E,def) -3
 E=varargin{1};
 r1=varargin{2};
 if isstruct(r1)
  r2=reshape(r1.def,size(r1.def,1)/2,[]);
  out=r1;
  if ~isempty(E)
   out.def=r2*E';
   out.data=(0:size(out.def,2)-1)';
  else; % Just return real and imaginary columns
    out.def=r2; 
    try;out.data=reshape([out.data';out.data'],size(out.data,2),[])';end
  end
 else
  r1=reshape([real(r1);imag(r1)],size(r1,1),[]);
  out=r1*E';
 end
 
 
else; error('%s',CAM);
end


function [T2,RC,RO]=genT2(T2,RC,RO,i2);
 %% #genT2 : 
 if strcmpi(RO.curSetType,'keep'); 
  T2=speye(length(i2));
 elseif strcmpi(RO.curSetType,'CaseT'); % Keep full bases on some DOF
  i3=fe_c(def.DOF(i2),Case.DOF,'ind',2); RO.usedIndDof=[RO.usedIndDof;i2(i3)];
  i2(i3)=[];
  i3={intersect(i2,RO.EdgeDof(:,1)),intersect(i2,RO.EdgeDof(:,2)), ...
      setdiff(i2,RO.EdgeDof(:))};
  i3=cellfun(@(x)fe_c(Case.DOF,def.DOF(x),'ind'),i3,'uni',0);
  T2=struct('Tl',Case.T(:,i3{1}),'Tr',Case.T(:,i3{2}),'Ti',Case.T(:,i3{3}),...
      'adof',{cellfun(@(x)Case.DOF(x),i3,'uni',0)});
  fprintf('%20s: CaseT L%i : R%i : I%i\n',RO.P2Sets{j1,3},cellfun(@length,i3));
 elseif strcmpi(RO.curSetType,'di');
   di=RO.di;%fprintf('Using finite learning');dbstack; 
    di.def=reshape(di.def,size(di.def,1),[]);
    di=feutilb('placeindof',SE.DOF,di);
    [u,s]=svd(di.def(i2,:),0);s=diag(s); 
    [T2,fr]=fe_norm(u,m(i2,i2),RC.k);
 elseif strncmpi(RO.curSetType,'fe_norm',7);
   if RO.curActive
     m=evalin('caller','T3''*m*T3');
     k=evalin('caller','T3''*k*T3');
     [T2,fr]=fe_norm(T2(:,any(T2)),m,k);
   else
     m=evalin('caller','m(i2,i2)');
     [T2,fr]=fe_norm(T2(:,any(T2)),m,RC.k);
   end
   if contains(RO.curSetType,'lrik');RO.curCoor='lrik';end
 else; 
   if isfield(RC,'coef'); 
       T2=T2*RC.coef;% Pcond scale coef for SvdTol
   % cf=feplot(10);cf.def=struct('def',T2,'DOF',SE.DOF(i2));fecom(';colordataA;coloralpha');
   end 
   %[T2,s]=svd(T2,0,'vector');T2=T2(:,s>RO.SvdTol*s(1));
   if ~sdtm.Contains(RO.curSetType,'lrilu');RO.curCoor='lrilu';end
   RC.SvdTol=RO.SvdTol; % How delayed to fe_coor
   %fprintf('%20s: svd_1=%.1e nkept=%i, %s\n',RO.P2Sets{j1,3},s(1),size(T2,2),RO.curCoor);
 end
 nind=sparse(i2,1,1:length(i2));
 RC.EdgeDof=full(nind(RO.EdgeDof(ismember(RO.EdgeDof(:,1),i2),:)));


%% #RvePerDef_gen rve periodic disp
function [out,out1,out2]=RvePerDef(varargin)

model=varargin{1};RO=varargin{2};data=varargin{3};
Case=evalin('caller','Case');

data.opt(end+1:4)=0;

ang=zeros(size(model.Node,1),1);
for j1=1:3 % Loop on xyz
  Ld=data.opt(j1+1);
  if Ld==0||RO.nc(j1)==0
  else; ang=ang+model.Node(:,4+j1)/(norm(RO.CellDir(:,j1))*RO.nc(j1));
  end
end
Z=zeros(size(ang));ex=exp(2i*pi.*ang);
d1=struct('def',[[ex;Z;Z] [Z;ex;Z] [Z;Z;ex] ], ...%[Z;ex;ex] [ex;Z;ex] [ex;ex;Z]
    'DOF',feutil('getdof',(1:3)'/100,model.Node(:,1)));
d1=feutilb('placeindof',model.DOF,d1);
if ~isfield(RO,'K_homo')
  RO=fe_homo('RveKHomo',model,Case,RO);
end
RO.b=[real(d1.def);imag(d1.def)]; % B with spatial real/spatial imag
RO.b=[RO.K_homo{1}*real(d1.def) RO.K_homo{2}*real(d1.def)
      RO.K_homo{1}*imag(d1.def) RO.K_homo{2}*imag(d1.def)]; % B with spatial real/spatial imag
if isfield(RO,'PerDefi');RO.b=RO.b(:,RO.PerDefi);end
RO.blab=cellfun(@(x)sprintf('u%i',x),num2cell(1:size(RO.b,2)),'uni',0);
out=RO;

%% #dftTip : datatip handling
function out=dftTip(varargin) %#ok<STOUT>

obj=varargin{1}; evt=varargin{2};[CAM,Cam]=comstr(varargin{3},1);

if comstr(Cam,'rangech')
 % clean seek of feplot ch   
 go=evt.Target; ga=ancestor(go,'axes');ua=get(ga,'userdata');

 if ~isfield(ua,'cf');
  return; %warning('No ua.cf : cannot update channel');return;
  %cf=feplot;
 else; cf=ua.cf;
 end
 if ~isa(cf,'sdth');cf=get(cf,'userdata');end
 r1=cf.vfields.DefF;
 if isempty(r1)||isempty(r1{1});return;else;r1=r1{1};end
 def=r1;r1=r1.data; 
 if size(r1,2)<2; warning('Inconsistent def in feplot(%i)',cf.opt(1));return;end
 if isa(evt,'matlab.graphics.internal.DataTipEvent');pos=evt.Position;
 elseif isfield(evt,'pos');pos=evt.pos;
 else; pos=evt.Position;
 end
 ch=find(r1(:,1)==pos(1)); % seek frequency
 if isempty(ch)&&~isempty(strfind(ga.XLabel.String,'kHz'));
  ch=find(abs(r1(:,1)-pos(1)*1000)<1e-6);
 end
 if length(ch)>1; ch=ch(ua.ch);end % Instance of freq at correct k
 if ~isempty(ch) % dispersion curve
 elseif strcmpi(go.Type,'surface')
  r2=go.XData; [r3,i1]=min(abs(r2-pos(1)));pos(:,1)=r2(i1);
  r2=go.YData; [r3,i2]=min(abs(go.YData-pos(1,2)));pos(:,2)=r2(i2);
  ch=find(r1(:,1)==pos(1)&r1(:,2)==i2);
 elseif isfield(ua,'PlotType')&&any(strcmpi(ua.PlotType,{'surface','contour'}))
  go=handle(ua.ob(1));
  r2=go.XData; [r3,i1]=min(abs(r2-pos(1)));pos(:,1)=r2(i1);
  r2=go.YData; [r3,i2]=min(abs(go.YData-pos(1,2)));pos(:,2)=r2(i2);
  ch=find(r1(:,1)==pos(1)&r1(:,strcmp(def.Xlab{2},'jPar'))==i2);
 else;
  ch=find(r1(:,1)==pos(2)); % seek frequency will generate jpar
 end
 
 gf1=handle(cf.opt(1));
 go=findall(gf1,'tag','TipCh'); 
 if isempty(go);go=uicontrol('parent',gf1,'tag','TipCh','visible','off');end
 if length(dbstack)<15
   fecom(cf,sprintf('ch%i',ch));
 else
  set(go,'userdata',ch);  % Update tip channel but do nothing else
 end
end

%% #BuildUb : build analytic strains
function [R2,d2,model]=BuildUb(model,RO);
 data=fe_case(model,'getdata','Symmetry');
 if ~isfield(RO,'CellDir')||isempty(RO.CellDir)
  % KUBC with non periodic mesh
  data.IntNodes= ...
   feutil('findnode x==|x==',model,min(model.Node(:,5)),max(model.Node(:,5)));
  data.IntYNodes= ...
   feutil('findnode y==|y==',model,min(model.Node(:,6)),max(model.Node(:,6)));
  data.IntZNodes= ...
   feutil('findnode z==|z==',model,min(model.Node(:,7)),max(model.Node(:,7)));
  R2=struct('DOF',feutil('getdof',model));
  R2.c=fe_c(R2.DOF,[data.IntNodes;data.IntYNodes;data.IntZNodes]);
 else
  if isempty(data)
   model=fe_homo(sprintf('dftBuild epsl %.15g%s',RO.epsl,RO.silent),model,RO.CellDir);
   data=fe_case(model,'getdata','Symmetry');
   model=fe_case(model,'cyclic','Symmetry',data);
   if isfield(model,'Case');
      model.Case=stack_set(model.Case,'cyclic','Symmetry',data);
   end
  elseif isfield(model,'Case')&&isempty(stack_get(model.Case,'cyclic','Symmetry'))
      dbstack; keyboard;
  end
  % Remove edges to avoid circular repetition (wrong not needed)
  %if length(RO.CellDir)>1
  % data.IntYNodes(ismember(data.IntYNodes(:,1),data.IntNodes(:,1)),:)=[];
  %end
  R2=dftu('clGen',RO,model);RO.DOF=R2.DOF;
  c=R2.cLx-R2.cRx;
  if isfield(R2,'cLy');  c=[c;R2.cLy-R2.cRy];end
  if isfield(R2,'cLz');  c=[c;R2.cLz-R2.cRz];end
  R2.c=c;
 end
if ~isfield(RO,'btype');RO.btype='kubc';end
switch lower(RO.btype)
case {'pbc','kubc','simpleload'}
 % Select reference homogeneous field
 % [a,b]=feval(fe_homo('@BuildUb'),model,struct('btype','kubc'))
 d1={'e11',struct('dir',{{'x','0','0'}},'DOF',[.01;.02;.03])
     'e22',struct('dir',{{'0','y','0'}},'DOF',[.01;.02;.03])
     'e33',struct('dir',{{'0','0','z'}},'DOF',[.01;.02;.03])
     'e12',struct('dir',{{'y','0','0'}},'DOF',[.01;.02;.03])
     'e23',struct('dir',{{'0','0','y'}},'DOF',[.01;.02;.03])
     'e13',struct('dir',{{'0','0','x'}},'DOF',[.01;.02;.03])
     'exx',struct('dir',{{'x','0','0'}},'DOF',[.01;.02;.03])
     'eyy',struct('dir',{{'0','y','0'}},'DOF',[.01;.02;.03])
     'ezz',struct('dir',{{'0','0','z'}},'DOF',[.01;.02;.03])
     'eyz',struct('dir',{{'0','z/2','y/2'}},'DOF',[.01;.02;.03])
     'ezx',struct('dir',{{'z/2','0','x/2'}},'DOF',[.01;.02;.03])
     'exy',struct('dir',{{'y/2','x/2','0'}},'DOF',[.01;.02;.03])
     };
 if isfield(RO,'Load'); d1=d1(ismember(d1(:,1),lower(RO.Load)),:);
 else; d1=d1(7:end,:); RO.Load=d1; 
 end
 %% Build analytic solution on full volume
 d2=[];%d2=struct('def',zeros(length(model.DOF),size(d1,1)),'DOF',model.DOF);
 %d2.lab=d1(:,1);  
 for j1=1:size(d1,1)
  d3=elem0('VectFromDirAtDof',model,d1{j1,2},R2.DOF);d3.lab=d1(j1,1);
  d2=fe_def('appenddef',d2,d3);if j1>1;d2.dir(j1,:)=d3.dir;end
 end
 
case 'mubc'
  %%CL INNER FACE NORMAL X
  NodesFace=sort(data.IntNodes(:));
  d1={'exx',struct('dir',{{'x'}},'DOF',[.01]), ...
   struct('dir',{{'0'}},'DOF',[.02]), ...
   struct('dir',{{'0'}},'DOF',[.03])
   'eyy',struct('dir',{{'0'}},'DOF',[.01]), ...
   struct('dir',{{'y'}},'DOF',[.02]), ...
   struct('dir',{{'0'}},'DOF',[.03])
   'ezz',struct('dir',{{'0'}},'DOF',[.01]), ...
   struct('dir',{{'0'}},'DOF',[.02]), ...
   struct('dir',{{'z'}},'DOF',[.03])
   'eyz',struct('dir',{{'0'}},'DOF',[.01]), ...
   struct('dir',{{'0','y/2'}},'DOF',[.01;.03]), ...
   struct('dir',{{'0','z/2'}},'DOF',[.01;.02])
   'ezx',struct('dir',{{'0','x/2'}},'DOF',[.02;.03]), ...
   struct('dir',{{'0'}},'DOF',[.02]), ...
   struct('dir',{{'z/2','0'}},'DOF',[.01;.02])
   'exy',struct('dir',{{'x/2','0'}},'DOF',[.02;.03]), ...
   struct('dir',{{'y/2','0'}},'DOF',[.01;.03]), ...
   struct('dir',{{'0'}},'DOF',[.03])
   };
  d2=struct('def',zeros(length(model.DOF),size(d1,1)),'DOF',model.DOF);
  d2.lab=d1(:,1); 
  for j1=1:size(d1,1)
   d3=elem0('VectFromDirAtDof',model,d1{j1,2},feutil('getdof',data.IntNodes(:),d1{j1,2}.DOF));
   d4=elem0('VectFromDirAtDof',model,d1{j1,3},feutil('getdof',data.IntYNodes(:),d1{j1,3}.DOF));
   d3=feutilb('placeindof',unique(round([d3.DOF;d4.DOF]*100))/100,d3);
   d4=feutilb('placeindof',d3.DOF,d4);d3.def=d3.def+d4.def;
   d4=elem0('VectFromDirAtDof',model,d1{j1,4},feutil('getdof',data.IntZNodes(:),d1{j1,4}.DOF));
   d3=feutilb('placeindof',unique(round([d3.DOF;d4.DOF]*100))/100,d3);
   d4=feutilb('placeindof',d3.DOF,d4);d3.def=d3.def+d4.def;
   d1{j1,5}=d3;
   TIn=fe_c(model.DOF,d3.DOF,d3.def')';T=fe_c(model.DOF,fe_c(model.DOF,d3.DOF,'dof',2))';
   d1{j1,1}=struct('T',T,'TIn',TIn,'DOF',model.DOF,'name',d1{j1,1});
  end
  R2.li=d1(:,1);
  
otherwise 
        error('Not implemented');
end
R2.V=prod(max(model.Node(:,5:7))-min(model.Node(:,5:7)));% Not very good

 %% #histToId
function  r2=histToId(RO,def)
  r2=struct('marker','xy','po',RO.hist);r2.po.Y(:,:,2:end)=[];
  if nargin==2
   r2.po.Y(r2.po.Y>def.X{1}(end))=NaN;r2.po.Y=r2.po.Y;
  end
  r2.po.X{1}(r2.po.X{1}==1|r2.po.X{1}==0)=NaN;
  r2.zfun=@(xa,ya,ga)ii_plp('plzcrop',ga,{'linewidth',1,'linestyle','-'});
  r2.po.MainDim='y';

 function   r1=toDD(C1,Range,r1);
%% #toDD : return DD

 r1=struct('dd',C1); 

 function   r1=toOrtho(C1,Range,r1);
%% #toOrtho : extract orthotropic moduli

if nargin==1||(nargin==3&&isfield(r1,'pl'))
 %dd stiffness matrix, c softness 
 c=inv(C1);
M(1,1)=1./c(1,1);
M(1,9)=conj(-c(1,2)*M(1,1));% nu12
M(1,8)=-c(1,3)*M(1,1); % nu13
M(1,2)=1./c(2,2);
M(1,7)=-c(2,3)*M(2); % nu12
M(1,3)=1./c(3,3);
M(1,6)=1./c(6,6);
M(1,5)=1./c(5,5);
M(1,4)=1./c(4,4);
if nargout==0||(nargin==3&&isfield(r1,'pl'))
  %% set model property keep rho
  if nargin<3||~isfield(r1,'pl')||isempty(r1.pl);r1=struct('pl',1);end
  if nargin<2;RO=struct;else; RO=Range;end
  if ~isfield(r1,'unit');r1.unit='US';end
  M(10)=-c(3,1)*M(3); % SDT uses nu31 rather than 13
  if ~isscalar(r1.pl)
     mat=feutil(sprintf('getpl %i -struct1',r1.pl(1)),r1);
  elseif isfield(r1,'Rho');
    mat=sdth.sfield('addselected',struct,r1,{'Rho','Eta'});
    r1=sdtm.rmfield(r1,{'Rho','Eta'});
  else; mat.Rho=1e-10;
  end
  r1.pl=[r1.pl(1) fe_mat('m_elastic;',r1.unit,6) M([1 2 3 7 10 9 4 5 6]) mat.Rho];
  r1.type='m_elastic';
  if nargout==0
      feutilb('_writepl',r1);clear r1;
  end
else
    r1=struct('X',{{{'E1', 'E2', 'E3','G23','G13','G12','nu23','nu13','nu12'}'}}, ...
     'Xlab',{{'Comp'}},'Y',M(:));
end
else
 %% current value to extract certain components
 ind=[1 8 15 22 29 36 7 13 14];% reshape(1:36,6,6)
 if isempty(C1);
   C1=struct('X',{{{'c11','c22','c33','c44','c55','c66','c12','c13','c23'}', ...
       (1:size(Range.val,1))'}},'Xlab',{{'Comp','jPar'}}, ...
        'Y',zeros(9,size(Range.val,1)));
   
 end
 C1.Y(:,Range.jPar)=r1(ind);
 r1=C1;
end

 
%% #safeInitDfrf : performs assembly
function [mdl,RunOpt,SE,Case,Load,w,ft,carg]=safeInitDfrf(varg,carg);
    
 mdl=varg{carg};carg=carg+1;if iscell(mdl);mo1=mdl{1};else;mo1=mdl; end
 if carg<=length(varg); RunOpt=varg{carg};carg=carg+1;else; RunOpt=struct;end
 r1=stack_get(mo1,'info','oProp','get');if ~isempty(r1);RunOpt.oProp=r1;end
 [w,carg]=fe_def('DefFreq;',mo1,carg,varg);RunOpt.f=w;w=w*2*pi;

 if isfield(RunOpt,'MatSplit')
  %% Actually use a MatSplit call
  mdl=fevisco('MatSplit',mdl,RunOpt.MatSplit); 
  if isfield(RunOpt,'AssembleCall'); eval(RunOpt.AssembleCall);
  else;
   if ~isfield(RunOpt,'matdes')||isempty(RunOpt.matdes); RunOpt.matdes=[2 -1];end
   st=sprintf('assemble -matdes %s NoT -se gett load',sprintf('%i ',RunOpt.matdes));
   [SE,Case,Load]=fe_case(mdl,st);
  end
  
 elseif ~isfield(RunOpt,'NeedHomog')
  if isfield(RunOpt,'Range');RunOpt=dftu('RangeNc',RunOpt);end
  if ~isfield(RunOpt,'oProp');RunOpt.oProp={};end

  % assemble (if necessary) taking boundary conditions into account
  %[w,carg]=fe_def('DefFreq',mdl,carg,varargin);RunOpt.f=w;w=w*2*pi;
  if isfield(RunOpt,'AssembleCall'); eval(RunOpt.AssembleCall);
  else;
   if ~isfield(RunOpt,'matdes')||isempty(RunOpt.matdes); RunOpt.matdes=[2 1 3 4];end
   st=sprintf('assemble -matdes %s NoT -se gett load',sprintf('%i ',RunOpt.matdes));
   [SE,Case,Load]=fe_case(mdl,st);
  end
 elseif RunOpt.NeedHomog==1
 %% NeedHomog==1 : decompose Kh in real and imaginary (no split)
 %  Florian (3.90)
 %  assemble {'m','Khr','Khi', 'Km(Em)'} 
 %  assemble {'m','Kh(E'')','Kv(E'''')', 'Km(Em)'} 
  if isfield(RunOpt,'AssembleCall'); eval(RunOpt.AssembleCall);
  else;
   if ~isfield(RunOpt,'matdes')||isempty(RunOpt.matdes); RunOpt.matdes=[2 1 3 4];end
   st=sprintf('assemble -matdes %s NoT -se gett load',sprintf('%i ',RunOpt.matdes));
   [SE,Case,Load]=fe_case(mdl,st);
  end
  SE.Klab{1}='Mh';SE.Klab{strcmp(SE.Klab,'k')}='Khr';
  SE.Klab(strcmp(SE.Klab,'4'))={'Khi'};
  r1=stack_get(mdl,'info','oProp','get');if ~isempty(r1);RunOpt.oProp=r1;end

  if size(Load.def,2)==36&&all(ismember({'F1';'F2';'F3';'F4';'F5';'F6'},Case.Stack(:,2)))
   Load.def=-(Load.def(:,1:6)+Load.def(:,7:12)+Load.def(:,13:18)+ ...
     Load.def(:,19:24)+Load.def(:,25:30)+Load.def(:,31:36));
   Load=feutil('rmfield',Load,'bset');
  end
  i1=cellfun(@isempty,SE.K);SE.K(i1)=[];SE.Klab(i1)=[];SE.Opt(:,i1)=[];
  % mo1 (heterogeneous), mo2 homogeneous
  %mo1=SE;mo1.K=[];mo1.DOF=[];mo1=fe_case(mo1,'reset','DofLoad','Ref',Load);
  % mo2 Kd = homogenous
  mo2=SE;%pl=m_elastic('dbval 1 Strain');
  mo2.il(:,3)=0; %homogenized in the global coordinates system xxx may not be robust
  if isfield(RunOpt,'pl');pl=RunOpt.pl(1,:);else;pl=mo2.pl(1,:);end
  if ~isreal(pl); error('Not valid use NeedHomog2');end
  mo2.pl=[mo2.pl(:,1) repmat(pl(2:end),size(mo2.pl,1),1)];
  mo2.K=[];mo2.DOF=[];mo2.name='homog';[mo2,C2]=fe_case(mo2,'assemble -matdes 1 NoT -SE');
  RunOpt.ihomog=length(SE.Klab)+1;
  SE.Klab{end+1}='homog';SE.K{length(SE.Klab)}=mo2.K{1};SE.Opt(2,length(SE.K))=1;
  w=[];ft=[];
 elseif RunOpt.NeedHomog==2
  %% Uses a prior MatSplit to generate matrices
  if ~isfield(RunOpt,'matdes')||isempty(RunOpt.matdes); RunOpt.matdes=[2 -1];end
  st=sprintf('assemble -matdes %s NoT -se gett load',sprintf('%i ',RunOpt.matdes));
  [SE,Case,Load]=fe_case(mdl,st);
  w=[];ft=[];
    
 end
 [Load,ft]=fe_load('builduNoBuildtime',mdl,Case,Load,w,struct('Opt',0));

 
 
function out=getCellDir(varargin);
%% #getCellDir
[CAM,Cam]=comstr(varargin{1},1);
data=varargin{2};

if comstr(Cam,'vect') % Return CellDir column vectors
  if isfield(data,'Elt');data=fe_case(data,'getdata','Symmetry');end
  if isfield(data,'CellDir')
  elseif isfield(data,'trans');data.CellDir=data.trans(:);
  end
  if size(data.CellDir,1)==1;
   r1=diag(data.CellDir); if size(r1,1)<3;r1(3,1)=0;end      
  else; r1=data.CellDir;
  end
  out=r1;
elseif comstr(Cam,'dx') % length of x
  data.CellDir=getCellDir('vect',data);
  out=norm(data.CellDir(:,1));
else; error('%s unknown',CAM);
end
