function [out,out1,out2,out3] = esi2sdt(varargin)
% esi2sdt : Interface between esi FEM files and SDT
%
% model=esi2sdt('read',FileName)
%
%       Supported files for reading are
%         .pc   pam crash files
%
% See <a href="matlab: sdtweb _taglist esi2sdt">TagList</a>

%
%       E. Balmes, P. Lorong, G. Martin
%       Copyright (c) 1990-2025 by SDTools, All Rights Reserved.
%       Use esi2sdt('cvs') for revision information

if nargin==0
 sdtweb _taglist esi2sdt % see structure of esi2sdt file
 return
end

% Input parsing
obj=[]; evt=[];
if ~ischar(varargin{1}) % From GUI
 obj=varargin{1}; evt=varargin{2}; [CAM,Cam]=comstr(varargin{3},1); carg=4;
else;[CAM,Cam]=comstr(varargin{1},1); carg=2; % From script
end

if comstr(Cam,'read'); [CAM,Cam]=comstr(CAM,5);
 %% #Read-----------------------------------------------------------------1
 % Read an file describing a FE model in ESI pam crash format (*.pc)
 % an output the corresponding open fem model.
 % If second argment given and equal to 1 a modal analysis is performed.

 % Next developments :
 % - passer ?? text_scan pour l'analyse des lignes
 % - changer de nom : esi2sdt.m
 % - ajouter les fonction esiread, esiwrite (voir nasread ...)
 % - ajouter des tests : t_esi
 % - Utilisation de dictionnaires : container.map

 % Input parsing
 FileName=varargin{carg}; carg=carg+1; % Second arg = FileName
 if carg>nargin||~isstruct(varargin{carg});RO=struct; % No option in third arg
 else;RO=varargin{carg};carg=carg+1; % Option structure in third arg
 end
 % Then parse CAM for command options,
 % and assign default values to unspecified options.
 % values declared prior in the RO structure are not overriden
 [RO,st,CAM]=cingui('paramedit -DoClean',[ ...
  'def(#3#"Compute modal analysis and display in feplot")' ...
  ],{RO,CAM}); Cam=lower(CAM);

 disp(['> Reading: ',FileName])

 % Read keys :
 %  - NODE   : nb_nod, node_id(nb_node), node_xyz(nb_node,3)
 %  - TETR4  : nb_T4, tet4(nb_T4,4), tet4_id(nb_T4), tet4_mat(nb_T4)
 %  - TETR10 : nb_T10, tet10(nb_T10), tet10_id(nb_T10), tet10_mat(nb_T10)
 %  - MATER  :
 %      Type 1  : 3D, ?lastique, isotrope
 %                nb_mat_1, mat_1_E(nb_mat_1), mat_1_nu(nb_mat_1), mat_1_rho(nb_mat_1)
 %  - SPRING : nb_spring,
 %  - MASS   : nb_mass
 %  - PART   :
 %  - FRAME  :
 %  - BOUNC  :
 %  - MTOCO  : à finir

 nb_Part_Solid = 0 ;
 nb_Part_Spring = 0 ;
 nb_Part_Beam = 0 ;

 nb_nod = 0 ;

 nb_T4 = 0 ;
 nb_T10 = 0 ;
 nb_beam = 0 ;
 nb_spring = 0 ;
 nb_mass = 0 ;

 nb_mat_1 = 0 ;
 nb_mat_225 = 0 ;
 nb_mat_213 = 0 ;

 nb_frame = 0 ;
 nb_BOUNC = 0 ;
 nb_MTOCO = 0 ;

 mdl = struct('Node',[],'Elt',[],'pl',[],'il',[],'name',' ','unit','') ;

 %% First reading
% Init empty lists
  node_id = [];
  node_xyz = [];

  tet4 = [];
  tet4_id = [];
  tet4_part = [];
  tet4_mat = [];

  tet10 = [];
  tet10_id = [];
  tet10_part = [];
  tet10_mat = [];

  beam = [];
  beam_id = [];
  beam_part = [];
  beam_mat = [];

  spring_Id = [];
  spring_IDPRT = [];
  spring_IDNOD1 = [];
  spring_IFRA = [];

  mass_IDNOD = [];
  mass_Mixy = [];
  mass_J = [];

  frame_id = [];
  frame_inod = [];
  %frame_U = zeros(nb_frame,3) ;
  %frame_V = zeros(nb_frame,3) ;
  frame_basis = [];

  bounc_XYZUVW = [];
  bounc_IFRA = [];
  bounc_names = {};
  bounc_lst_nodes = {} ;

  mtoco_id = [];
  mtoco_IDNODi = [];
  mtoco_XYZUVW = [];
  mtoco_IFRA1 = [];
  mtoco_names ={};
  mtoco_lst_nodes = {};

  Part_Beam_IDPRT=[];
  Part_Beam_IDMAT=[];
  Part_Solid_IDPRT=[];
  Part_Solid_IDMAT=[];
  Part_Spring_IDPRT=[];
  Part_Spring_IDMAT=[];

 %% Reading
 tic ;
 fid = fopen(FileName,'r');
 txt=fileread(FileName);
 fclose(fid);

 % Get header indiced in char array
 headind=strfind(txt,'$#');headind(end+1)=length(txt)+1;
 % First header
 mdl=pc_header(mdl,txt,headind);
 % Loop on headers
 jhead=0;
 while jhead<length(headind)-1
  [st,jhead,head]=nextBlock(txt,jhead,headind);
  if isempty(st); continue; end
  i1=regexp(st(1:7),'....../'); % Detect data type : NODE, TETR4, ...
  if i1; typ=st(1:6);
  else; typ='';
  end
  switch typ
   case 'NODE  '
    %% #Node-------------------------------------------------------------2
    % Get data from first line
    i2=0:65:length(st); % Line break indices (16*4+1=65)
    i3=i2+(1:8)'; % indice containing NODE  / 
    st(i3)=' '; % => replace NODE  / by empty char
    st(i2(2:end))=''; % Remove line breaks
    % Reshape by block of 16 charaters
    st=reshape(st,16,[]);
    st(end+1,:)=' '; % Add whitespace between each block of character
    % Now use sscanf to convert into numeric value
    n1=reshape(sscanf(st,'%16d%16f%16f%16f'),4,[])';
    %n1=sscanf(st(i1+1:end),'NODE  / %8d%16f%16f%16f\n');

    % Store data (indice in col1 and  x y z in column 5 6 7)
    mdl.Node(nb_nod+1:nb_nod+size(n1,1),[1 5:7])=n1;
    % Increment node indice
    nb_nod = nb_nod+size(n1,1);

   case {'TETR4 ','TETR10','BEAM  ','SPRING','MASS  '}
    %% #Elt--------------------------------------------------------------2
    % Add Elt with list of node ids and ESI partid as matid and proid
    % Conversion from ESI partid to matid and proid is done in SDT model
    % filling section
    % For each element type : fill elt lines and eltname
    % elt line format : [nodes_id partid partid eltid]
    if strncmpi(typ,'tetr4',5)
     %% #Tetr4------------------------------------------------------------3
     error('Need revise, contact Guillaume Martin');
     % Get data from first line
     n1=sscanf(st,'TETR4 / %8d%8d%8d%8d%8d%8d'); % eltid partid nodes_id
     % Attempt to read all following TETR4 lines
     n2=fscanf(fid,'TETR4 / %8d%8d%8d%8d%8d%8d\n',[6 Inf]);
     % Store data [node_list partid partid eltid]
     elt=[n1 n2]'; elt=elt(:,[3:6 2 2 1]);
     eltname='tetra4'; % SDT eltname
     nb_T4=nb_T4+size(elt,1);
    elseif strncmpi(typ,'tetr10',6)
     %% #Tetr10-----------------------------------------------------------3
     % Get first element skipping headers and get remaining elements
     [st2,jhead,head]=nextBlock(txt,jhead,headind);
     [st3,jhead,head]=nextBlock(txt,jhead,headind);
     st=sprintf('%s\n%s\n%s',st,st2,st3);
     % 17 blocks of 8 digits + 3 line return = 139 characters
     i2=0:139:length(st);
     i3=i2+([1:8 25:41 106:122])'; % indice containing TETR10  / and empty 8 spaces
     st(i3)='';
     i2=8*12+1:8*12+1:length(st); % line return indices
     st(i2)='';
     % Reshape by block of 8 charaters
     st=reshape(st,8,[]);
     st(end+1,:)=' '; % Add whitespace between each block of character
     % Now use sscanf to convert into numeric value
     elt=reshape(sscanf(st,repmat('%8d',1,12)),12,[])';
     % Store data [node_list partid partid eltid]
     elt=elt(:,[3:12 2 2 1]);
     eltname='tetra10'; % SDT eltname
     nb_T10=nb_T10+size(elt,1);
    elseif strncmpi(typ,'beam',4)
     %% #Beam-------------------------------------------------------------3
     % 10 blocks of 8 digits + 1 line return = 81 characters
     i2=0:81:length(st);
     i3=i2+([1:8 49:80])'; % indice containing TETR10  / and unused columns
     st(i3)=''; 
     i2=8*5+1:8*5+1:length(st); % line return indices
     st(i2)='';
     % Reshape by block of 8 charaters
     st=reshape(st,8,[]);
     st(end+1,:)=' '; % Add whitespace between each block of character
     % Now use sscanf to convert into numeric value
     elt=reshape(sscanf(st,repmat('%8d',1,5)),5,[])';
     % Store data [node_list partid partid node_bending_plane 0 0 eltid]
     elt=[elt(:,[3 4 2 2 5]) zeros(size(elt,1),2) elt(:,1)];
     eltname='beam1'; % SDT eltname
     nb_beam=nb_beam+size(elt,1);
    elseif strncmpi(typ,'spring',6)
     %% #Spring-----------------------------------------------------------3
     error('Need revise, contact Guillaume Martin');
     if nb_spring==0; warning('Should speed up using fscanf, contact Guillaume Martin'); end
     nb_spring = nb_spring + 1 ;
     IDEL = str2double(st(9:16));
     IDPRT = str2double(st(17:24));
     IDNOD1 = str2double(st(25:32));
     IDNOD2 = str2double(st(33:40));
     IFRA = str2double(st(41:48));
     if IFRA ~= 0
      warning('SPRING : IFRA ~= 0, not yet used')
     end
     if IDNOD2 ~= 0
      error('SPRING : IDNOD2 ~= 0 not yet handled')
     end
     spring_Id(nb_spring) = IDEL ;
     spring_IDPRT(nb_spring) = IDPRT ;
     spring_IDNOD1(nb_spring) = IDNOD1 ;
     spring_IFRA(nb_spring) = IFRA ;

    elseif strncmpi(typ,'mass',4)
     %% #Mass-------------------------------------------------------------3
     error('Need revise, contact Guillaume Martin');
     if nb_mass==0; warning('Should speed up using fscanf, contact Guillaume Martin'); end
     nb_mass = nb_mass+ 1 ;
     IDNOD = str2double(st(9:16));
     IFRA = str2double(st(17:24));
     if IFRA ~= 0
      error('MASS : IFRA ~= 0 not yet handled')
     end
     st=nextLine(fid) ;
     %NAME = line(6:length(line)) ;
     st=nextLine(fid) ;
     Mx = str2double(st(9:24));
     My = str2double(st(25:40));
     Mz = str2double(st(41:56));
     st=nextLine(fid) ;
     Ix = str2double(st(9:24));
     Iy = str2double(st(25:40));
     Iz = str2double(st(41:56));
     st=nextLine(fid) ;
     Ixy = str2double(st(9:24));
     Iyz = str2double(st(25:40));
     Ixz = str2double(st(41:56));
     mass_IDNOD(nb_mass) = IDNOD ;
     mass_Mixy(nb_mass) = [Mx, My, Mz] ;
     mass_J(nb_mass) = [Ix, Iy, Iz, Ixy, Iyz, Ixz] ;
    end
    mdl.Elt=feutil('AddElt',mdl.Elt,eltname,elt);
   case 'TITLE '
    mdl.name = menage_str(st(8:length(st))) ;

   case 'MATER '
    %% #MATER-------------------------------------------------------------2
    IDMAT = str2double(st(9:16));
    MATYP = str2double(st(17:24));
    RHO = str2double(st(25:40));
    [st,jhead,head]=nextBlock(txt,jhead,headind) ;
    % Skip name ? For now yes
    [st,jhead,head]=nextBlock(txt,jhead,headind) ;
    if MATYP == 1
     nb_mat_1 = nb_mat_1 + 1 ;
     G = str2double(st(1:10));
     [st,jhead,head]=nextBlock(txt,jhead,headind) ;
     K = str2double(st(1:10));
     mat_1_id(nb_mat_1) = IDMAT ;
     mat_1_E(nb_mat_1) = 9*K*G/(3*K+G);
     mat_1_nu(nb_mat_1)= (3*K-2*G)/(2*(3*K+G));

     %                     mat_1_E(nb_mat_1) = K ;
     %                     mat_1_nu(nb_mat_1)= K/(2*G)-1;
     mat_1_rho(nb_mat_1)= RHO;

    elseif MATYP == 225
     error('need revise. Contact GM')
     nb_mat_225 = nb_mat_225 + 1;
     mat_225_id(nb_mat_225) = IDMAT ;

     MASS = str2double(st(11:20));
     INERTIA = str2double(st(21:30));
     [st,jhead,head]=nextBlock(txt,jhead,headind) ;
     STIFTR = str2double(st(1:10));
     DAMVTR = str2double(st(11:20));
     KSITR  = str2double(st(21:30));
     [st,jhead,head]=nextBlock(txt,jhead,headind) ;
     STIFTS = str2double(st(1:10));
     DAMVTS = str2double(st(11:20));
     KSITS  = str2double(st(21:30));
     [st,jhead,head]=nextBlock(txt,jhead,headind) ;
     STIFTT = str2double(st(1:10));
     DAMVTT = str2double(st(11:20));
     KSITT  = str2double(st(21:30));
     [st,jhead,head]=nextBlock(txt,jhead,headind) ;
     STIFRR = str2double(st(1:10));
     DAMVRR = str2double(st(11:20));
     KSIRR  = str2double(st(21:30));
     [st,jhead,head]=nextBlock(txt,jhead,headind) ;
     STIFRS = str2double(st(1:10));
     DAMVRS = str2double(st(11:20));
     KSIRS  = str2double(st(21:30));
     [st,jhead,head]=nextBlock(txt,jhead,headind) ;
     STIFRT = str2double(st(1:10));
     DAMVRT = str2double(st(11:20));
     KSIRT  = str2double(st(21:30));
     [st,jhead,head]=nextBlock(txt,jhead,headind) ;
     TTTR3 = str2double(st(1:10));
     RRRR3 = str2double(st(11:20));
     % Pour eviter les raideur nulles
     mat_225_stiff(nb_mat_225,:) = [ STIFTR, STIFTS, STIFTT, STIFRR, STIFRS, STIFRT] ;

    elseif MATYP == 213
     % Beam elasto-plastic behaviour
     nb_mat_213 = nb_mat_213 + 1;
     mat_213_id(nb_mat_213) = IDMAT ;
     mat_213_E(nb_mat_213) = str2double(st(1:10));
     mat_213_nu(nb_mat_213)= str2double(st(11:20));
     mat_213_rho(nb_mat_213)= RHO;

    else
     disp('MATER, MATYP:')
     disp(MATYP)
     error('MATER: MATYP not yet handled')
    end

   case 'PART  '
    %% #PART--------------------------------------------------------------2
    st2=strsplit(st,'\n');
    st3=st2{1};
    IDPRT = str2double(st3(9:16));
    ATYPE = menage_str(st3(17:24));
    IDMAT = str2double(st3(25:32));
    st3=st2{2};
    NAME = st3(6:length(st3)) ;
    %disp(['IDPRT = ',num2str(IDPRT),'  ATYPE = ',ATYPE, ...
    %      '  IDMAT = ',num2str(IDMAT), '  NAME = ',NAME])
    if strcmp(ATYPE(1:4),'BEAM')
     nb_Part_Beam = nb_Part_Beam +1;
     Part_Beam_IDPRT(nb_Part_Beam) = IDPRT ;
     Part_Beam_IDMAT(nb_Part_Beam) = IDMAT ;
     % Skip 7 headers
     [st,jhead,head]=nextBlock(txt,jhead,headind,7) ;
     Part_Beam_IDSEC(nb_Part_Beam) = str2double(st(1:5));
     Part_Beam_CA(nb_Part_Beam) = str2double(st(11:20));
     if Part_Beam_IDSEC(nb_Part_Beam) ~= 2
      Part_Beam_B(nb_Part_Beam) = str2double(st(21:30));
      if Part_Beam_IDSEC(nb_Part_Beam) == 3
       Part_Beam_C(nb_Part_Beam) = str2double(st(31:40));
      end
     end
    elseif strcmp(ATYPE(1:5),'SOLID') || strcmp(ATYPE(1:5),'TETRA')
     nb_Part_Solid = nb_Part_Solid + 1 ;
     Part_Solid_IDPRT(nb_Part_Solid) = IDPRT;
     Part_Solid_IDMAT(nb_Part_Solid) = IDMAT;

    elseif strcmp(ATYPE(1:6),'SPRING')
     nb_Part_Spring = nb_Part_Spring + 1;
     Part_Spring_IDPRT(nb_Part_Spring) = IDPRT ;
     Part_Spring_IDMAT(nb_Part_Spring) = IDMAT ;
    else
     disp('PART, ATYPE:')
     disp(ATYPE)
     error('PART: ATYPE not yet handled')
    end

   case 'FRAME '
    %% #FRAME-------------------------------------------------------------2
    nb_frame = nb_frame + 1;
    frame_id(nb_frame) = str2double(st(9:16));
    IFRATY = str2double(st(17:24));
    IAXIS = str2double(st(25:32));
    if IFRATY ~= 1
     error('Frames : IFRATY ~= 1 not yet handled')
    end
    if IAXIS ~= 0
     error('Frames : IAXIS ~= 0 not yet handled')
    end
    st2=strsplit(st,'\n');
    st3=st2{3};
    NAME = menage_str(st3(5:length(st3)));
    [st,jhead,head]=nextBlock(txt,jhead,headind);
    Ux = str2double(st(9:24));
    Uy = str2double(st(25:40));
    Uz = str2double(st(41:56));
    [st,jhead,head]=nextBlock(txt,jhead,headind);
    Vx = str2double(st(9:24));
    Vy = str2double(st(25:40));
    Vz = str2double(st(41:56));
    frame_inod(nb_frame) = str2double(st(57:64));
    frame_basis(nb_frame,:,:) = basis_UV([ Ux, Uy, Uz ],[ Vx, Vy, Vz ]) ;

   case 'BOUNC '
    %% #BOUNC-------------------------------------------------------------2
    nb_BOUNC = nb_BOUNC + 1 ;
    IDNOD = str2double(st(9:16));
    bounc_XYZUVW(nb_BOUNC) = str2double(st(19:24));
    bounc_IFRA(nb_BOUNC) = str2double(st(25:32));
    ISENS = str2double(st(33:40));
    if ISENS ~= 0
     warning('BOUNC : ISENS ~= 0 not yet used')
    end
    st2=strsplit(st,'\n');
    st3=st2{2};
    bounc_names{nb_BOUNC} = st3(6:length(st3)) ;
    %NAME = strcat(num2str(nb_BOUNC), ':',NAME) ;
    if IDNOD == 0
     lines=st2(3:end);
     bounc_lst_nodes{nb_BOUNC} = read_nodesID(lines) ;
    else
     bounc_lst_nodes{nb_BOUNC} = IDNOD ;
    end

   case 'MTOCO '
    %% #MTOCO-------------------------------------------------------------2
    nb_MTOCO = nb_MTOCO + 1 ;
    mtoco_id(nb_MTOCO) = str2double(st(9:16));
    mtoco_IDNODi(nb_MTOCO) = str2double(st(17:24));
    mtoco_XYZUVW(nb_MTOCO) = str2double(st(27:32));
    mtoco_IFRA1(nb_MTOCO) = str2double(st(33:40));
    st2=strsplit(st,'\n');
    st3=st2{3};
    mtoco_names{nb_MTOCO} = st3(6:length(st3)) ;
    st3=st2(5:end);
    %lst_nodes = read_nodesID(line,fid) ;
    mtoco_lst_nodes{nb_MTOCO} = read_nodesID(st3);

    %sdtBC = esi2sdt_BC(XYZUVW) ;
    %mdl=fe_case(mdl,'rigid',NAME,[IDNODi sdtBC lst_nodes]);

   otherwise
    if isempty(typ) % Header without NODE  / ; TETR10 / ;..... => skip
    else; disp(['    Untreated key: ',typ]);
    end
  end
 end

 disp('> Summary of read data:')
 if nb_Part_Solid  > 0 ; disp(['    nb Part Solid  : ',num2str(nb_Part_Solid)]) ; end
 if nb_Part_Beam   > 0 ; disp(['    nb Part Beam   : ',num2str(nb_Part_Beam)]) ; end
 if nb_Part_Spring > 0 ; disp(['    nb Part Spring : ',num2str(nb_Part_Spring)]) ; end
 if nb_nod         > 0 ; disp(['    nb nodes       : ',num2str(nb_nod)])  ; end
 if nb_T4          > 0 ; disp(['    nb T4          : ',num2str(nb_T4)])  ; end
 if nb_T10         > 0 ; disp(['    nb T10         : ',num2str(nb_T10)]) ; end
 if nb_beam        > 0 ; disp(['    nb beam        : ',num2str(nb_beam)]) ; end
 if nb_spring      > 0 ; disp(['    nb spring      : ',num2str(nb_spring)]); end
 if nb_mass        > 0 ; disp(['    nb mass        : ',num2str(nb_mass)]) ; end
 if nb_mat_1       > 0 ; disp(['    nb mat_1       : ',num2str(nb_mat_1)]) ; end
 if nb_mat_225     > 0 ; disp(['    nb mat_225     : ',num2str(nb_mat_225)]) ; end
 if nb_mat_213     > 0 ; disp(['    nb mat_213     : ',num2str(nb_mat_213)]) ; end
 if nb_BOUNC       > 0 ; disp(['    nb BOUNC       : ',num2str(nb_BOUNC)]) ; end
 if nb_MTOCO       > 0 ; disp(['    nb MTOCO       : ',num2str(nb_MTOCO)]) ; end
 if nb_frame       > 0 ; disp(['    nb FRAME       : ',num2str(nb_frame)]) ; end


 %% #SDT_model_filling----------------------------------------------------2
 disp('> SDT model filling')

 %% Node declaration
 usr2idx = containers.Map(mdl.Node(:,1),1:size(mdl.Node,1)) ;

 %% Association element / materiau (esi uses PART)
 % Assign matid from partid
 mpid=feutil('mpid',mdl.Elt);
 partid=[Part_Solid_IDPRT Part_Beam_IDPRT Part_Solid_IDPRT];
 matid=[Part_Solid_IDMAT Part_Beam_IDMAT Part_Solid_IDMAT];
 [i1,i2]=ismember(mpid(:,1),partid);
 mpid(i1~=0,1)=matid(i2(i1~=0)); % Replace each partid by matid
 mdl.Elt=feutil('mpid',mdl.Elt,mpid); % Apply back

 %% Element declaration
 %  [Element_type_declaration ;
 %   NodeNumbers MatId ProId EltId OtherInfo ]
 %  MatId material ID number in model.pl
 %  ProId element property (e.g. section area) ID number in model.il

 if 1==2 % To propagate at reading section when revising spring and mass
  if nb_spring > 0
   mdl.Elt(end+1,1:6)=[Inf abs('celas')];
   for iElt=1:nb_spring
    IDPRT = spring_IDPRT(iElt) ;
    IDMAT_spring = -1 ;
    for iPart=1:nb_Part_Spring
     if Part_Spring_IDPRT(iPart) == IDPRT
      IDMAT_spring = Part_Spring_IDMAT(iPart) ;
      break
     end
    end
    if IDMAT_spring == -1
     disp('spring_IDPRT:')
     disp(spring_IDPRT)
     disp('Part_Spring_IDPRT:')
     disp(Part_Spring_IDPRT)
     error('SPRING: IDMAT not found')
    end
    iMat225_OK = -1 ;
    for iMat225 = 1:nb_mat_225
     if mat_225_id(iMat225) == IDMAT_spring
      %disp(mat_225_stiff(iMat225,:))
      if mat_225_stiff(iMat225,1) ~= 0
       mdl.Elt(end+1,1:7)=[spring_IDNOD1(iElt) 0 -1 1 0 0 mat_225_stiff(iMat225,1)];
      end
      if mat_225_stiff(iMat225,2) ~= 0
       mdl.Elt(end+1,1:7)=[spring_IDNOD1(iElt) 0 -2 2 0 0 mat_225_stiff(iMat225,1)];
      end
      if mat_225_stiff(iMat225,3) ~= 0
       mdl.Elt(end+1,1:7)=[spring_IDNOD1(iElt) 0 -3 3 0 0 mat_225_stiff(iMat225,3)];
      end
      if mat_225_stiff(iMat225,4) ~= 0
       mdl.Elt(end+1,1:7)=[spring_IDNOD1(iElt) 0 -4 4 0 0 mat_225_stiff(iMat225,4)];
      end
      if mat_225_stiff(iMat225,5) ~= 0
       mdl.Elt(end+1,1:7)=[spring_IDNOD1(iElt) 0 -5 5 0 0 mat_225_stiff(iMat225,5)];
      end
      if mat_225_stiff(iMat225,6) ~= 0
       mdl.Elt(end+1,1:7)=[spring_IDNOD1(iElt) 0 -6 6 0 0 mat_225_stiff(iMat225,6)];
      end
      iMat225_OK = iMat225 ;
      break
     end
    end
    if iMat225_OK == -1
     disp('spring_IDPRT:')
     disp(spring_IDPRT)
     disp('Part_Spring_IDPRT:')
     disp(Part_Spring_IDPRT)
     disp('Part_Spring_IDMAT:')
     disp(Part_Spring_IDMAT)
     disp('mat_225_id:')
     disp(mat_225_id)
     error('SPRING: Mat 225 not found')
    end
   end
  end
 
  if nb_mass > 0
   mdl.Elt(end+1,1:6)=[Inf abs('mass1')];
   for imass = 1:nb_mass
    mdl.Elt(end+1,1:7) = [mass_IDNOD(imass), ...
     mass_Mixy(imass,1), mass_Mixy(imass,2), mass_Mixy(imass,3), ...
     mass_J(imass,1), mass_J(imass,2), mass_J(imass,3), ...
     imass] ;
    if (mass_J(imass,4) ~= 0) || (mass_J(imass,5) ~= 0) || (mass_J(imass,6) ~= 0)
     error('MASS: rectangular inertia product not yet handled')
    end
   end
  end
 end


 %% Materials declaration

 for iMat_1 = 1:nb_mat_1
  %disp(['iMat_1 ',num2str(iMat_1)])
  mdl.pl = [ mdl.pl ; ...
   mat_1_id(iMat_1), ...
   fe_mat('m_elastic','SI',1), ...
   mat_1_E(iMat_1), ...
   mat_1_nu(iMat_1), ...
   mat_1_rho(iMat_1)] ;
 end
 for iMat_213 = 1:nb_mat_213
  %disp(['iMat_213 ',num2str(iMat_213)])
  mdl.pl = [ mdl.pl ; ...
   mat_213_id(iMat_213), ...
   fe_mat('m_elastic','SI',1), ...
   mat_213_E(iMat_213), ...
   mat_213_nu(iMat_213), ...
   mat_213_rho(iMat_213)] ;
 end

 %% Properties declaration
 nbColil = 0 ;
 if nb_Part_Solid > 0
  nbColil = max(nbColil,6) ;
 end
 if nb_Part_Beam > 0
  nbColil = max(nbColil,8) ;
 end
 nb_Part = nb_Part_Solid + nb_Part_Beam;
 mdl.il = zeros(nb_Part,nbColil) ;

 iPart = 0 ;
 for iPartS = 1:nb_Part_Solid
  iPart = iPart + 1 ;
  mdl.il(iPart,1:6) = [ Part_Solid_IDPRT(iPartS), ...
   fe_mat('p_solid','SI',1), ...
   0, -3, 0, 0] ; % -3 = default element integ rule
 end

 for iPartB = 1:nb_Part_Beam
  iPart = iPart + 1 ;
  if Part_Beam_IDSEC(iPartB) == 1
   mdl.il(iPart,1:6) = [ Part_Beam_IDPRT(iPartB), ...
    fe_mat('p_beam','SI',3), ...
    0, ...
    comstr('TUBE',-32), ...
    Part_Beam_CA(iPartB)+Part_Beam_B(iPartB)/2,...
    Part_Beam_CA(iPartB)-Part_Beam_B(iPartB)/2] ;
  elseif Part_Beam_IDSEC(iPartB) == 2
   mdl.il(iPart,1:5) = [ Part_Beam_IDPRT(iPartB), ...
    fe_mat('p_beam','SI',3), ...
    0, ...
    comstr('ROD',-32), ...
    Part_Beam_CA(iPartB)] ;
  elseif Part_Beam_IDSEC(iPartB) == 3
   mdl.il(iPart,1:8) = [ Part_Beam_IDPRT(iPartB), ...
    fe_mat('p_beam','SI',3), ...
    0, ...
    comstr('BOX',-32), ...
    Part_Beam_B(iPartB)+Part_Beam_C(iPartB)/2,...
    Part_Beam_CA(iPartB)+Part_Beam_C(iPartB)/2,...
    Part_Beam_C(iPartB),...
    Part_Beam_C(iPartB)] ;
  elseif Part_Beam_IDSEC(iPartB) == 4
   mdl.il(iPart,1:6) = [ Part_Beam_IDPRT(iPartB), ...
    fe_mat('p_beam','SI',3), ...
    0, ...
    comstr('BAR',-32), ...
    Part_Beam_B(iPartB),...
    Part_Beam_CA(iPartB)] ;

  end
 end


 %% Boundary condition (disp) and MPC

 for i_bounc=1:nb_BOUNC
  lst_nodes = cell2mat(bounc_lst_nodes(i_bounc)) ;
  nb_nodes = length(lst_nodes) ;

  % Si parmi les noeuds de lst_nodes certains sont des noeuds maitres
  % de MTOCO il faut les supprimer.
  % C'est provisoire ....
  nb_master_nodes = 0 ;

  lst_nodes_F = [] ;
  for inod=1:nb_nodes
   node_bc = lst_nodes(inod) ;
   for i_mtoco=1:nb_MTOCO
    if node_bc == mtoco_IDNODi(i_mtoco)
     nb_master_nodes = 1 ;
    end
   end
   if nb_master_nodes == 0
    lst_nodes_F = [lst_nodes_F, node_bc] ;
   end
  end

  if length(lst_nodes_F) > 0
   if bounc_IFRA(i_bounc) == 0
    basis = eye(3) ;
   else
    % Need to use frame_id !
    ind=find(frame_id==bounc_IFRA(i_bounc)); % correction by GM 04/07/2023
    basis = reshape(frame_basis(ind,:,:),3,3) ;
   end
   mdl=add_BOUNC(mdl, bounc_names{i_bounc}, bounc_XYZUVW(i_bounc), ...
    lst_nodes_F, basis);
  end
 end

 for i_mtoco=1:nb_MTOCO
  if mtoco_IFRA1(i_mtoco) == 0
   basis = eye(3) ;
  else
   % Need to use frame_id !
   ind=find(frame_id==mtoco_IFRA1(i_mtoco)); % correction by GM 04/07/2023
   basis = reshape(frame_basis(ind,:,:),3,3) ;
  end
  mdl=add_MTOCO(mdl, mtoco_names{i_mtoco}, mtoco_IDNODi(i_mtoco), mtoco_XYZUVW(i_mtoco), ...
   cell2mat(mtoco_lst_nodes(i_mtoco)), basis, ...
   usr2idx);
 end

 disp('> Model filling completed')
 t_second = toc;
 disp(['> Reading time: ',num2str(t_second),' s'])

 %% Plotting
 %disp('pc2sdt: SDt plot')
 if RO.def

  nb_modes = 12;

  %analyse_modale(mdl, nb_modes)

  % Voir avec Guillaume Martin la signification de ce qui suit ! (01/10/2021)
  %ofact('method mklserv_utils -silent');
  %def=fe_eig(mdl,[5 nb_modes 1000]);

  %cf0 = feplot(mdl) ;
  %cf0.def = def;

  def=fe_eig(mdl,[5 nb_modes 0]);
  cf=feplot(mdl,def);

 end

 out=mdl;
elseif comstr(Cam,'write'); [CAM,Cam]=comstr(Cam,6);
 %% #Write----------------------------------------------------------------1
 error('Writing not supported yet');

 %% classical commands
elseif comstr(Cam,'@'); out=eval(CAM);
elseif comstr(Cam,'cvs')
 out=sdtcheck('Revision');
else; error('Unknown command %s',CAM);
end %commands
end % baseFunction
%% #SubFunc --------------------------------------------------------------1
function mdl=pc_header(mdl,txt,headind)
%% #pc_header : store header in model, interprete unit--------------------2
 st=txt(1:headind(1)-1);
 st2=regexp(st,'UNIT\s*(\w*)\s*(\w*)\s*(\w*)\s*(\w*)\s*\n','tokens');
 if isempty(st2); sdtw('_nb','Model unit interpretation failed');
 else;
  st2=st2{1};
  if strcmpi(st2{1},'MM')&&strcmpi(st2{2},'KG')
   mdl.unit='MM';
  else; warning('Deal with unit system (%s,%s)',st2{1},st2{2});
  end
 end
 mdl=stack_set(mdl,'info','header',st);
end

function [st,jhead,head]=nextBlock(txt,jhead,headind,nb_line)
%% #nextBlock-------------------------------------------------------------2
if nargin == 3
 nb_line = 1 ;
end
jhead=jhead+nb_line;
if jhead>=length(headind); st=''; head=''; return; end % Header in last line, return empty 
st=txt(headind(jhead):headind(jhead+1)-1);
i1=regexp(st,'\n','once'); 
head=st(1:i1-1);
st=st(i1+1:end);
i2=regexp(st,'\n *$'); % Remove last line return and empty spaces
if ~isempty(i2); st(i2:end)=''; end

end
function lst_nodesID = read_nodesID(st)
%% #read_nodesID----------------------------------------------------------2
lst_nodesID = [];
j1=0;
while j1<length(st)
 j1=j1+1; st2=st{j1};
 if strcmp(st2(1:11),'        NOD')
  n1=str2num(st2(12:end));
  lst_nodesID=[lst_nodesID n1]; 
 end
end
end

function ligne_prop = menage_str(ligne)
%% #menage_str------------------------------------------------------------2
deb = 1 ;
while strcmp(ligne(deb),' ') && (deb ~= length(ligne))
 deb = deb + 1 ;
end
fin = length(ligne) ;
while strcmp(ligne(fin),' ') && (fin ~= 1)
 fin = fin - 1 ;
end
if (fin==1) || (deb==length(ligne))
 ligne_prop = ' ';
else
 ligne_prop = ligne(deb:fin) ;
end
end

function sdtBC = esi2sdt_BC(esiBC)
%% #esi2sdt_BC------------------------------------------------------------2
% Convert Boundary Condition (disp) frome ESI (VPS) format to SDT format :
%   Example of ESI format    : 111111, 011001
%   Corresponding SDT format : 123456,    236
% input:
%   esiBC: integer
% output:
%   sdtBC: integer

if esiBC == 0
 error('esiBC == 0 : such a bc is useless')
end

BC = zeros(6,1) ;

if esiBC >= 100000
 BC(1) = 1;
 esiBC = esiBC - 100000;
end
if esiBC >= 10000
 BC(2) = 1;
 esiBC = esiBC - 10000;
end
if esiBC >= 1000
 BC(3) = 1;
 esiBC = esiBC - 1000;
end
if esiBC >= 100
 BC(4) = 1;
 esiBC = esiBC - 100;
end
if esiBC >= 10
 BC(5) = 1;
 esiBC = esiBC - 10 ;
end
if esiBC >= 1
 BC(6) = 1;
end

sdtBC = 0 ;
expo = -1 ;
for ii=6:-1:1
 if BC(ii) == 1
  expo = expo + 1;
  sdtBC = sdtBC + ii*10^expo ;
 end
end

end

function mdl=add_BOUNC(mdl, NAME, XYZUVW, lst_nodes, basis)
%% #add_BOUNC-------------------------------------------------------------2
% XYZUVW : integer, 6 digit to specify bc (111111)
% lst_nodes : 1 row array with a list of user id of nodes
% basis : 3x3 array, ROWS are vectors components

debug = 0 ;

XYZUVW_list = XYZUVW_split(XYZUVW) ;

nb_nodes = length(lst_nodes) ;
lst_dof = reshape([lst_nodes+0.01 ; lst_nodes + 0.02 ; lst_nodes + 0.03 ;
 lst_nodes+0.04 ; lst_nodes + 0.05 ; lst_nodes + 0.06],[],1) ;
c=cell(1,6); % Constraint matrix for each direction
if debug
 disp(strcat("BOUNC : ",NAME))
end
for ii=1:6
 if XYZUVW_list(ii)
  if debug
   disp(['  dof : ',num2str(ii)])
  end
  const_mat = zeros(nb_nodes,6*nb_nodes) ;
  if ii <= 3
   iiv = ii;
  else
   iiv = ii - 3;
  end
  for inod=1:nb_nodes
   if ii <= 3 % translation dofs
    const_mat(inod,6*(inod-1)+1:6*(inod-1)+3) = basis(iiv,:) ;
   else % rotation dofs
    const_mat(inod,6*(inod-1)+4:6*(inod-1)+6) = basis(iiv,:) ;
   end
   if debug
    disp(lst_dof)
    disp(const_mat)
   end
  end
  c{ii}=const_mat;
  % mdl=fe_case(mdl,'mpc',[char(NAME),'_bc_dof',num2str(ii)],data);
 end
end
% Merge data and store in Case (one case entry for all directions)
data=struct('DOF',lst_dof,'c',vertcat(c{:}));
r1=fe_case('stack_get',mdl,[NAME,'_bc'],'g');
if isempty(r1) % Add new case
 mdl=fe_case(mdl,'mpc',[NAME,'_bc'],data);
else; error('Merging several BOUNC not handled yet, contact Guillaume Martin');
end

%data=struct('DOF',[1.01;1.03],'c',[1 -1]);
%mo2=fe_case(mo1,'mpc','z=x',data);
end

function mdl=add_MTOCO(mdl, NAME, IDNODi, XYZUVW, lst_nodes, basis, usr2idx)
%% #add_MTOCO-------------------------------------------------------------2
debug = 0 ;
if debug
 disp(strcat("MTOCO : ",NAME))
end

nb_nodes = length(lst_nodes) ;
XYZUVW_list = XYZUVW_split(XYZUVW) ;
nodes = mdl.Node ;

lst_nd2 = [IDNODi, lst_nodes ] ;

lst_dof = reshape([lst_nd2+0.01 ; lst_nd2 + 0.02 ; lst_nd2 + 0.03 ; ...
 lst_nd2+0.04 ; lst_nd2 + 0.05 ; lst_nd2 + 0.06],[],1) ;
if debug
 disp(lst_dof')
end

master_idx = usr2idx(IDNODi);
xyz_master = nodes(master_idx,5:7);
c=cell(1,6); % Constraint matrix for each direction
for ii=1:6
 %disp(['ii : ',num2str(ii)])
 if XYZUVW_list(ii)
  if debug
   disp(['  dof : ',num2str(ii)])
  end
  const_mat = zeros(nb_nodes,length(lst_dof)) ;
  if ii <= 3
   % d : direction vector
   % M, U_M, theta : master node, location and displacement
   % and rotation
   % S, U_S : slave node, location and displacement
   % -d.U_M + d.U_S + cross(d,MS).theta_M = 0
   iiv = ii;
   for inod=1:nb_nodes
    node_idx = usr2idx(lst_nodes(inod));
    xyz_n = nodes(node_idx,5:7);
    const_mat(inod,1:3) = -basis(iiv,:) ;
    const_mat(inod,4:6) = ...
     cross(basis(iiv,:),xyz_n-xyz_master) ;
    const_mat(inod,6*inod+1:6*inod+3) = basis(iiv,:) ;
   end
  else
   iiv = ii - 3;
   for inod=1:nb_nodes
    const_mat(inod,4:6) = -basis(iiv,:) ;
    const_mat(inod,6*inod+4:6*inod+6) = basis(iiv,:);
   end
  end

  if debug
   disp(const_mat)
  end
  c{ii}=const_mat;
  %mdl=fe_case(mdl,'mpc',[char(NAME),'_mpc_dof',num2str(ii)],data);
 end
end
% Merge data and store in Case (one case entry for all directions)
data=struct('DOF',lst_dof,'c',vertcat(c{:}));%,'slave',(1:6)');
mdl=fe_case(mdl,'mpc',[NAME,'_mpc'],data);
% % Add mass1 to master node (allows to keep mpc with feutilb submodel)
% i1=feutil('findnode inelt{eltname mass1};',mdl);
% if ~ismember(IDNODi,i1)
%  mdl=feutil('AddElt',mdl,'mass1',IDNODi);
% end
%mdl=feutil('AddElt',mdl,'mass1',IDNODi);
%data=struct('DOF',[1.01;1.03],'c',[1 -1]);
%mo2=fe_case(mo1,'mpc','z=x',data);

end

function XYZUVW_list = XYZUVW_split(XYZUVW)
%% #XYZUVW_split----------------------------------------------------------2
if XYZUVW == 0
 error('XYZUVW == 0 : such a bc is useless')
end

XYZUVW_list = zeros(6,1) ;

if XYZUVW >= 100000
 XYZUVW_list(1) = 1;
 XYZUVW = XYZUVW - 100000;
end
if XYZUVW >= 10000
 XYZUVW_list(2) = 1;
 XYZUVW = XYZUVW - 10000;
end
if XYZUVW >= 1000
 XYZUVW_list(3) = 1;
 XYZUVW = XYZUVW - 1000;
end
if XYZUVW >= 100
 XYZUVW_list(4) = 1;
 XYZUVW = XYZUVW - 100;
end
if XYZUVW >= 10
 XYZUVW_list(5) = 1;
 XYZUVW = XYZUVW - 10 ;
end
if XYZUVW >= 1
 XYZUVW_list(6) = 1;
end

end

function basis = basis_UV(U,V)
%% #basis_UV--------------------------------------------------------------2


X = U / sqrt(sum(U.*U)) ;
Z = cross(X,V) ;
Z = Z / sqrt(sum(Z.*Z)) ;
Y = cross(Z,X) ;

basis = zeros(3) ;
basis(1,:) = X ;
basis(2,:) = Y ;
basis(3,:) = Z ;

end