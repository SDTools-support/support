function [out,out1,out2]=d_dft(varargin); %#ok<*STOUT>

% D_DFT Support for demonstrations related to periodic computations in SDT
%
% See <a href="matlab: sdtweb _taglist d_dft">TagList</a>

%       Etienne Balmes, SDTools
%       Copyright (c) 1990-2025 by SDTools, All Rights Reserved.
%       For revision information use d_dft('cvs')

if nargin==0
 if sdtdef('isinteractive');d_dft('tuto');return;
 else;% Non interactive test
  d_dft('ScriptNida');
  d_dft('ScriptRve');
 end
end

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*DEFNU,*NOSEM>
[CAM,Cam]=comstr(varargin{1},1);carg=2;
if carg>nargin||~isstruct(varargin{2});RO=struct('info',1);
else;RO=varargin{carg};carg=carg+1;
end

%% #Script -------------------------------------------------------------------
if comstr(Cam,'script');[CAM,Cam]=comstr(CAM,7);

if comstr(Cam,'nida')
%% #TutoNida : reduced dispertion diagram for honeycomb model

%% Step1 : compute a reduced model and display interactive dispersion diagram
 mo1=d_mesh('rvenida',struct('lc',3));% feplot(mo1)
 mo1=feutil('setmat 1 eta .02',mo1); mo1=feutil('setmat 2 eta .02',mo1);%2% loss
 mo1=fe_case(mo1,'DofLoad','In',struct('def',1,'DOF',1.03));
 sens=struct('cta',1,'DOF',1.03,'tdof',1.03,'lab',{{'1z'}});
 mo1=fe_case(mo1,'SensDof','Out',sens);
 %mo1=d_mesh('rvenida',struct('lc',1));% feplot(mo1)
 
 mo1=fe_cyclic('BuildEqualY',mo1);
 
 R1=struct('BuildList',{{'Eig1',[5 20 0];'Eig50',[5  20 0]}},'matdes',[2 1 3 4]);
 R1.Phase2={'fe_homo','DftRedp2LU'};R1.EdgeTol=1e-5;
 SE=fe_homo('dftRedList',mo1,R1);
  
 Range=struct('ncx',1./linspace(1e-3,.49,30)');
 SE=stack_set(SE,'info','EigOpt',[5 30 1e3]);
 [def,hist]=fe_homo('dftDisp -UseLong',SE,struct('Range',Range,'fmax',1e7));

 RD=struct('mno',(0:4)','cf',102,'ci',112,'ViewHist',hist,'ktype','kcx');
 fe_homo('dfpInitSelDef',SE,def,RD);fecom('ShowFiCEvalz')

%% Step2 : Build a 2D force response and use interactive display

 % Now compute response with reduced model
freq=1000*linspace(0,40,100)';SE=stack_set(SE,'info','Freq',freq);
RR=struct('Range',Range,'fmax',1e7);

[dmode,hist,dFrf]=fe_homo('dftDisp-FRF -UseLong',SE,RR);%Reduced basis Dfrf

% Interactivity between 2D FRF and shapes (uses DfpSensObserve)
RE=struct('mno',(0:4)','cf',102,'ci',112,'ktype','kc','sens',sens,'PlotInfo','surface');
fe_homo('dfpInitSelDef',SE,dFrf,RE);fecom('ShowFiCEvalz')
iicom(';xlin;ylin');ci=iiplot;fe_homo('dfpDataTip',ci);

%% Step3 Check full response of initial model
jf=[4 10 30];jk=[1 10 20];RF=RR;RF.Range.ncx=RF.Range.ncx(jk,:);
df=fe_homo('dftdfrf',stack_set(mo1,'info','Freq',freq(jf)),RF);

dfr=fe_homo('dftredRest',SE,fe_homo('dftdfrf',stack_set(SE,'info','Freq',freq(jf)),RF));
for j1=1:size(df.def,2); disp(dfr.def(:,j1)\df.def(:,j1));end

 C1=fe_homo('dfpSensObserve',mo1,sens,df,struct('ktype','kc'));
 C1b=fe_homo('dfpSensObserve',mo1,sens,dfr,struct('ktype','kc'));
 r3=[reshape(C1b.Y,2,[]);reshape(C1.Y,2,[])];
 r3(1:2,:)./r3(3:4,:) % This is not perfect be kind of expected since not exactly on freq


%% Step4 :  Recombine harmonics xxx needs further checks
 C1=struct('type','cfrf');C1.sens=sens;C1=fe_homo('DfpViewCurve',dFrf,C1); C2=C1.GetData;

 % Initialize multi-cell display
 RB=struct('mno',(0:5)','cf',102,'ktype','kc');
 fe_homo('dfpInitSelDef',SE,df,RB);fecom('ShowFiCEvalz')

 
%% EndTuto

elseif comstr(Cam,'rve')
%% #TutoRve : homogeneization testing

mo1=d_mesh('RveLayered',struct('x',.1,'y',.1,'z',2,'Per',3,'nz',3,'nx',2,'ny',2));
mo1=fe_case(mo1,'remove','1D');
feutilb('_write',mo1)

[dd1,d1]=fe_homo('RveKubc',mo1);
[d2,dd2]=fe_homo('RveSubc',mo1);

RP=struct('Range',struct('ncx',1,'ncy',1,'ncz',[3 10 100 1000]));
%dd1.Y([1 6 5;6 2 4;5 4 3])
RP.pl=m_elastic('FormulaDDtoOrtho',dd2.dd);RP.matdes=[2 1 4 -1];
[d1,C1]=fe_homo('rvePerDef',mo1,RP);

squeeze(C1.Y(:,:,2))./squeeze(C1.Y(:,:,strcmpi(C1.X{3},'Ed')))
RE=struct('mno',(0:4)'*[0 0 1],'cf',102);
fe_homo('dfpInitSelDef',mo1,d1,RE);fecom('ShowFiCEvalZ')%sdtweb fe_homo dfpSensObs

% EndTuto

elseif comstr(Cam,'converge')
%% #TutoConverge : mesh convergence analysis

sdtweb src18 converge

% EndTuto


else; error('Script%s unknown',CAM);
    
end

%% #Solve -------------------------------------------------------------------
elseif comstr(Cam,'solve');[CAM,Cam]=comstr(CAM,6);

if comstr(Cam,'reducebase+f1')
%% #Base+f1 Enforce base acceleration with correction for first frequency -2

elseif comstr(Cam,'fmaxstudy')
%% #SolveFmaxStudy -2

model=RO; RO=varargin{carg};carg=carg+1;

for jpar=1:length(RO.fmax)
 model.il(2,4:5)=[2 RO.fmax(jpar)*2*pi]; 
 [d1,SE]=p_pml('SolveDfrf',model);
 r2=fe_def('subDof',d1,RO.In+.01);
 if jpar==1;C1=struct('X',{{r2.data}},'Xlab',{{'Freq','fmax'}},'Y',r2.def(:));end
 C1.X{2}=RO.fmax(1:jpar)*2*pi;C1.Y(:,jpar)=r2.def(:);C1.Ylab='Resp';
 if isfield(RO,'jdisp')&&RO.jdisp==jpar
  cf=comgui('guifeplot -reset',2);cf.model=SE;cf.def=d1;
  fecom('colordataEvalA');
 end
end
ci=comgui('guiiiplot -reset',3);
iicom(ci,'curveinit','Test',C1);iicom(ci,';submagpha;chall;cax1;xlog;cax2;xlog');


else;error('Solve%s unknown',CAM);
end
elseif comstr(Cam,'view')

if comstr(Cam,'viewdebugt2')
 %% debug views for fe_homo P2Sets, fe_coor LriLU  
 eval(iigui({'T3','SE','T2','RC'},'GetInCaller')) 
 c10=feplot(10,';');
 if comstr(Cam,'T2')
   if isfield(T2,'adof')&&iscell(T2.adof{:})
    c10.def=struct('def',T3*[T2.Tl T2.Tr T2.Ti],'DOF',SE.DOF,'adof',vertcat(T2.adof{:}));
   end
   % t_cyclic('debugCoorView')
   % cf=feplot(10,';');cf.def=struct('def',T3*Ti,'DOF',SE.DOF);
 end
else; error('%s',CAM)
end

%% #Tuto: recover model from a specific tuto step -3
elseif comstr(Cam,'tuto'); 
 eval(sdtweb('_tuto',struct('file','d_dft','CAM',CAM)));
 if nargout==0; clear out; end

elseif comstr(Cam,'cvs')
 out=sdtcheck('revision','$Revision: a244ccc $  $Date: 2021-04-15 20:25:44 +0200 $');
else;error('%s',CAM);

%% #End function
end




