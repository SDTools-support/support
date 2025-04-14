function [out,out1,out2]=pak2sdt(varargin)

% 

%       Etienne Balmes, Jean-Philippe Bianchi
%       Copyright (c) 1990-2025 by SDTools, All Rights Reserved.
%       $Revision: 1.29 $  $Date: 2015/10/06 16:36:20 $

obj=[]; evt=[];
if nargin==0;[CAM,Cam]=comstr('init',1);carg=1;
elseif ischar(varargin{1})
  [CAM,Cam]=comstr(varargin{1},1);carg=2;
else; obj=varargin{1};evt=varargin{2}; [CAM,Cam]=comstr(varargin{3},1);carg=4;
end

persistent checked UI
if isempty(checked)  % Verify updates compared to SDT 6.5
 UI=sdtroot('PARAMUI');
end

%#ok<*ASGLU,*NASGU,*NOSEM,*NBRAK>

% sdtweb('_tracker','pje',704)

%% #Init : GUI tab generation ------------------------------------
if comstr(Cam,'init'); [CAM,Cam]=comstr(CAM,5);
 
PARAM=sdtroot('PARAMVh;');
if comstr(Cam,'ftree')
 r1=PARAM.PakFile;
 table={'DataSet'};level=1;
 r2=struct2cell(r1.Measurements.DataSets)';
 r2(:,4)=cellfun(@(x)sprintf('%i x %i',length(x(1).Data),length(x)), ...
      r2(:,4),'uni',0);
 table(end+(1:size(r2,1)),1:3)=r2(:,[4 2 1]);level(end+1:size(table,1),1)=2;
 ua=struct('table',{table},'level',level, ...
     'name','PAK','ToolTip', ...
     'Content of pak file','NeedClose',1,'setSort',4, ...
     'ColWidth',[100 40 -1]); 
  uo=feval('cinguj','tabbedpaneAdd',UI.gf,ua); 
  set(UI.gf,'visible','on');
  
else; % xxx sdtroot
end

%% #Init : GUI tab generation ------------------------------------
elseif comstr(Cam,'set'); [CAM,Cam]=comstr(CAM,4);
    
    
elseif comstr(Cam,'test'); [CAM,Cam]=comstr(CAM,5);
 % cd 'D:\balmes\Dropbox\PAK2sdtools\PAK2Mat Example\ExampleMIMO'
 PARAM=sdtroot('PARAMVh');
 PARAM.PakFile=load('Test_01.mat');
 pak2sdt('initftree')
 
 fe_range('uiRO2Table',r1.Measurements,[],struct('MaxCol',10,'MaxLevel',3))
 
elseif comstr(Cam,'cvs')
 out=sdtcheck('revision');
else;error('%s unknown',CAM);
end
