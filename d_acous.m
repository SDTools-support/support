function [out,out1]=d_acous(varargin);

% D_ACOUS Support for demonstrations related to acoustics
%
% See <a href="matlab: sdtweb _taglist d_acous">TagList</a>
%
% Etienne Balmes, SDTools


%       Copyright (c) 1990-2025 by SDTools, All Rights Reserved.
%       For revision information use d_acous('cvs')

if nargin==0

 d_acous('MeshBoxWithFelt')
 return

end

%#ok<*NASGU,*ASGLU,*CTCH,*TRYNC,*DEFNU>
[CAM,Cam]=comstr(varargin{1},1);carg=2;
if carg>nargin||~isstruct(varargin{2});RO=struct('info',1);
else;RO=varargin{carg};carg=carg+1;
end

%% #Script -------------------------------------------------------------------
if comstr(Cam,'script');[CAM,Cam]=comstr(CAM,7);

if comstr(Cam,'manual')
%% #ScriptManual : compute things by hand

error('not implemneted')

else; error('Script%s unknown');
    
end
%% #Mesh -------------------------------------------------------------------
elseif comstr(Cam,'mesh');[CAM,Cam]=comstr(CAM,5);

if comstr(Cam,'boxwithfelt')
%% #MeshBoxWithFelt simple example of a box with an absorbing wall
 
% Parameter handling
if carg<=nargin&&isstruct(varargin{carg});RO=varargin{carg};carg=carg+1;
else; RO=struct;end
[RO,st,CAM]=cingui('paramedit -DoClean',[ ...
   'lx(0.505#%g#"box length")'...
   'ly(0.505#%g#"box width")'...
   'h(0.43#%g#"patch height")'...
   'lc(0.01#%g#"desired length in fluid")'...
   'ys(.40#%g#"speaker y")'...
   'zs(.30#%g#"spearker z")'...
   'ws(.05#%g#"speaker width")'...
   'hs(.05#%g#"spearker height")'...
   'Quad(#3#"linear or quadratic elements")'...
   ],{RO,CAM});

%% Mesh

 % Length based split in x
 mo1=feutil('objectbeam 1 1',[0 0 0;RO.lx 0 0],1);
 mo1=feutil(sprintf('refinebeam%.15g',RO.lc),mo1);
 RunOpt.lx=sort(mo1.Node(:,5));

 % Proper split in y direction
 mo2=struct('Node', ...
    feutil('addnode',[],[0 RO.ly-[RO.ys RO.ys-RO.ws 0]]'*[0 1 0]));
 mo2.Elt=feutil('objectbeamline',1:size(mo2.Node,1));
 mo2=feutil(sprintf('refinebeam%.15g',RO.lc),mo2);
 RunOpt.ly=sort(mo2.Node(:,6));
 
 % Proper split in z direction
 mo2=struct('Node', ...
    feutil('addnode',[],[0 RO.h-[RO.zs RO.zs-RO.hs 0]]'*[0 1 0]));
 mo2.Elt=feutil('objectbeamline',1:size(mo2.Node,1));
 mo2=feutil(sprintf('refinebeam%.15g',RO.lc),mo2);
 RunOpt.lz=sort(mo2.Node(:,6));

 mo1=feutil('objecthexa 1 1',[0 0 0;1 0 0;0 1 0;0 0 1], ...
     RunOpt.lx,RunOpt.ly,RunOpt.lz);
 
 % Identify speaker as surface elements
 i1=feutil('findnode y>= & y<= & z>= & z<= & x==0',mo1, ...
     RO.ys,RO.ys+RO.ws,RO.zs,RO.zs+RO.hs);
 elt=feutil('selelt selface & innode',mo1,i1);
 mo1=feutil('addelt',mo1,elt);
 mo1.Elt=feutil('set group2 matid 2 proid 2',mo1);
 
 % Add absorbing wall
 i1=feutil('findnode x==',mo1,RO.lx);
 elt=feutil('selelt selface & innode',mo1,i1);
 mo1=feutil('addelt',mo1,elt);
 mo1.Elt=feutil('set group3 matid 3 proid 3',mo1);
 
 feplot(mo1); fecom colordatamat-alpha.4-edgealpha.05

%%%% Define material properties
mo1.pl=m_elastic('dbval -unit SI Air');
mo1.unit='SI';


% Integration rules for volumes
mo1=p_solid('default;',mo1);

if RO.Quad;mo1=feutil('lin2quad',mo1);end

out=mo1;  % Send output to out variable
%% MeshEnd
else;error('Mesh%s unknown');
end
    
%% clean end
elseif comstr(Cam,'@');out=eval(CAM);
elseif comstr(Cam,'cvs')
 out=sdtcheck('revision');
 %out='$Revision: 507 $  $Date: 2020-05-13 08:54:00 +0200 (Wed, 13 May 2020) $';
else; error('%s unknown',CAM);
end 
%% #End function
end
