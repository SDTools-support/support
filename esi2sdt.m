function [out,out1,out2,out3] = esi2sdt(varargin)
% esi2sdt : Interface between esi FEM files and SDT
%
% model=esi2sdt('read',FileName)
%
%       Supported files for reading are
%         .pc   pam crash files
%
%       E. Balmes, P. Lorong, G. Martin
%       Copyright (c) 1990-2023 by SDTools, All Rights Reserved.
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

 mdl = struct('Node',[],'Elt',[],'pl',[],'il',[],'name',' ','unit','SI') ;

 %% First reading
 disp('> First reading')

 tic ;
 fid = fopen(FileName,'r');
 while 1
  line = fgetl(fid);
  if ~ischar(line), break, end
  if length(line) > 7
   if strcmp(line(1:7),'NODE  /')
    nb_nod = nb_nod+1 ;
   elseif strcmp(line(1:7),'TETR4 /')
    nb_T4 = nb_T4 + 1 ;
   elseif strcmp(line(1:7),'TETR10/')
    nb_T10 = nb_T10 + 1 ;
   elseif strcmp(line(1:7),'BEAM  /')
    nb_beam = nb_beam + 1 ;
   elseif strcmp(line(1:7),'SPRING/')
    nb_spring = nb_spring + 1;
   elseif strcmp(line(1:7),'MASS  /')
    nb_mass = nb_mass + 1;
   elseif strcmp(line(1:7),'FRAME /')
    nb_frame = nb_frame + 1;
   elseif strcmp(line(1:7),'BOUNC /')
    nb_BOUNC = nb_BOUNC + 1;
   elseif strcmp(line(1:7),'MTOCO /')
    nb_MTOCO = nb_MTOCO + 1;
   end
  end
 end
 fclose(fid);
 t_first = toc;
 disp(['> Reading time: ',num2str(t_first),' s'])

 node_id = zeros(nb_nod,1) ;
 node_xyz = zeros(nb_nod, 3);

 tet4 = zeros(nb_T4, 4);
 tet4_id = zeros(nb_T4,1);
 tet4_part = zeros(nb_T4,1);
 tet4_mat = zeros(nb_T4,1);

 tet10 = zeros(nb_T10,10);
 tet10_id = zeros(nb_T10,1);
 tet10_part = zeros(nb_T10,1);
 tet10_mat = zeros(nb_T10,1);

 beam = zeros(nb_beam,3) ;
 beam_id = zeros(nb_beam,1);
 beam_part = zeros(nb_beam,1);
 beam_mat = zeros(nb_beam,1);

 spring_Id = zeros(nb_spring,1);
 spring_IDPRT = zeros(nb_spring,1);
 spring_IDNOD1 = zeros(nb_spring,1);
 spring_IFRA = zeros(nb_spring,1);

 mass_IDNOD = zeros(nb_mass,1) ;
 mass_Mixy = zeros(nb_mass,3) ;
 mass_J = zeros(nb_mass,6) ;

 frame_id = zeros(nb_frame,1) ;
 frame_inod = zeros(nb_frame,1) ;
 %frame_U = zeros(nb_frame,3) ;
 %frame_V = zeros(nb_frame,3) ;
 frame_basis = zeros(nb_frame,3,3) ;

 bounc_XYZUVW = zeros(nb_BOUNC, 1) ;
 bounc_IFRA = zeros(nb_BOUNC, 1) ;
 bounc_names = strings(nb_BOUNC, 1) ;
 bounc_lst_nodes = cell(nb_BOUNC, 1) ;

 mtoco_id = zeros(nb_MTOCO, 1) ;
 mtoco_IDNODi = zeros(nb_MTOCO, 1) ;
 mtoco_XYZUVW = zeros(nb_MTOCO, 1) ;
 mtoco_IFRA1 = zeros(nb_MTOCO, 1) ;
 mtoco_names = strings(nb_MTOCO, 1) ;
 mtoco_lst_nodes = cell(nb_BOUNC, 1) ;


 %% Second reading
 disp('> Second reading')
 tic ;
 fid = fopen(FileName,'r');

 i_nod = 0 ;
 i_T4 = 0 ;
 i_T10 = 0 ;
 i_beam = 0 ;
 i_spring = 0 ;
 i_mass = 0 ;
 i_frame = 0 ;
 i_bounc = 0 ;
 i_mtoco = 0 ;

 % Reading lines
 while 1
  line = fgetl(fid);
  if ~ischar(line), break, end
  if (length(line) > 7) && strcmp(line(7),'/')
   %disp(line)
   switch line(1:6)

    case 'NODE  '
     i_nod = i_nod+1 ;
     node_id(i_nod) = sscanf(line(9:16),'%8d') ;
     node_xyz(i_nod,1) = sscanf(line(17:32),'%16f') ;
     node_xyz(i_nod,2) = sscanf(line(33:48),'%16f') ;
     node_xyz(i_nod,3) = sscanf(line(49:64),'%16f') ;
     %cellfun(@str2double,textscan(tline,'%8c%16c%16c%16c'))

    case 'TETR4 '
     i_T4 = i_T4 + 1 ;
     tet4_id(i_T4) = sscanf(line(9:16),'%8d');
     tet4_part(i_T4) = sscanf(line(17:24),'%8d');
     tet4(i_T4,1) = sscanf(line(25:32),'%8d');
     tet4(i_T4,2) = sscanf(line(33:40),'%8d');
     tet4(i_T4,3) = sscanf(line(41:48),'%8d');
     tet4(i_T4,4) = sscanf(line(49:56),'%8d');

    case 'TETR10'
     i_T10 = i_T10 + 1 ;
     tet10_id(i_T10) = sscanf(line(9:16),'%8d');
     tet10_part(i_T10) = sscanf(line(17:24),'%8d');

     line=nextLine(fid) ;
     tet10(i_T10,1) = sscanf(line(17:24),'%8d');
     tet10(i_T10,2) = sscanf(line(25:32),'%8d');
     tet10(i_T10,3) = sscanf(line(33:40),'%8d');
     tet10(i_T10,4) = sscanf(line(41:48),'%8d');
     tet10(i_T10,5) = sscanf(line(49:56),'%8d');
     tet10(i_T10,6) = sscanf(line(57:64),'%8d');
     tet10(i_T10,7) = sscanf(line(65:72),'%8d');
     tet10(i_T10,8) = sscanf(line(73:80),'%8d');

     line=nextLine(fid) ;
     tet10(i_T10,9)  = sscanf(line(17:24),'%8d');
     tet10(i_T10,10) = sscanf(line(25:32),'%8d');

    case 'BEAM  '
     i_beam = i_beam + 1 ;
     beam_id(i_beam) = sscanf(line(9:16),'%8d');
     beam_part(i_beam) = sscanf(line(17:24),'%8d');
     beam(i_beam,1) = sscanf(line(25:32),'%8d');
     beam(i_beam,2) = sscanf(line(33:40),'%8d');
     beam(i_beam,3) = sscanf(line(41:48),'%8d');

    case 'SPRING'
     i_spring = i_spring + 1 ;
     IDEL = str2double(line(9:16));
     IDPRT = str2double(line(17:24));
     IDNOD1 = str2double(line(25:32));
     IDNOD2 = str2double(line(33:40));
     IFRA = str2double(line(41:48));
     if IFRA ~= 0
      warning('SPRING : IFRA ~= 0, not yet used')
     end
     if IDNOD2 ~= 0
      error('SPRING : IDNOD2 ~= 0 not yet handled')
     end
     spring_Id(i_spring) = IDEL ;
     spring_IDPRT(i_spring) = IDPRT ;
     spring_IDNOD1(i_spring) = IDNOD1 ;
     spring_IFRA(i_spring) = IFRA ;

    case 'MASS  '
     i_mass = i_mass+ 1 ;
     IDNOD = str2double(line(9:16));
     IFRA = str2double(line(17:24));
     if IFRA ~= 0
      error('MASS : IFRA ~= 0 not yet handled')
     end
     line=nextLine(fid) ;
     %NAME = line(6:length(line)) ;
     line=nextLine(fid) ;
     Mx = str2double(line(9:24));
     My = str2double(line(25:40));
     Mz = str2double(line(41:56));
     line=nextLine(fid) ;
     Ix = str2double(line(9:24));
     Iy = str2double(line(25:40));
     Iz = str2double(line(41:56));
     line=nextLine(fid) ;
     Ixy = str2double(line(9:24));
     Iyz = str2double(line(25:40));
     Ixz = str2double(line(41:56));
     mass_IDNOD(i_mass) = IDNOD ;
     mass_Mixy(i_mass) = [Mx, My, Mz] ;
     mass_J(i_mass) = [Ix, Iy, Iz, Ixy, Iyz, Ixz] ;

    case 'TITLE '
     mdl.name = menage_str(line(8:length(line))) ;

    case 'MATER '
     IDMAT = str2double(line(9:16));
     MATYP = str2double(line(17:24));
     RHO = str2double(line(25:40));
     line=nextLine(fid) ;
     line=nextLine(fid) ;
     %NAME = line(6:length(line)) ;
     line=nextLine(fid) ;
     if MATYP == 1
      nb_mat_1 = nb_mat_1 + 1 ;
      G = str2double(line(1:10));
      line=nextLine(fid) ;
      K = str2double(line(1:10));
      mat_1_id(nb_mat_1) = IDMAT ;
      mat_1_E(nb_mat_1) = 9*K*G/(3*K+G);
      mat_1_nu(nb_mat_1)= (3*K-2*G)/(2*(3*K+G));

      %                     mat_1_E(nb_mat_1) = K ;
      %                     mat_1_nu(nb_mat_1)= K/(2*G)-1;
      mat_1_rho(nb_mat_1)= RHO;

     elseif MATYP == 225
      nb_mat_225 = nb_mat_225 + 1;
      mat_225_id(nb_mat_225) = IDMAT ;

      MASS = str2double(line(11:20));
      INERTIA = str2double(line(21:30));
      line=nextLine(fid) ;
      STIFTR = str2double(line(1:10));
      DAMVTR = str2double(line(11:20));
      KSITR  = str2double(line(21:30));
      line=nextLine(fid) ;
      STIFTS = str2double(line(1:10));
      DAMVTS = str2double(line(11:20));
      KSITS  = str2double(line(21:30));
      line=nextLine(fid) ;
      STIFTT = str2double(line(1:10));
      DAMVTT = str2double(line(11:20));
      KSITT  = str2double(line(21:30));
      line=nextLine(fid) ;
      STIFRR = str2double(line(1:10));
      DAMVRR = str2double(line(11:20));
      KSIRR  = str2double(line(21:30));
      line=nextLine(fid) ;
      STIFRS = str2double(line(1:10));
      DAMVRS = str2double(line(11:20));
      KSIRS  = str2double(line(21:30));
      line=nextLine(fid) ;
      STIFRT = str2double(line(1:10));
      DAMVRT = str2double(line(11:20));
      KSIRT  = str2double(line(21:30));
      line=nextLine(fid) ;
      TTTR3 = str2double(line(1:10));
      RRRR3 = str2double(line(11:20));
      % Pour eviter les raideur nulles
      mat_225_stiff(nb_mat_225,:) = [ STIFTR, STIFTS, STIFTT, STIFRR, STIFRS, STIFRT] ;

     elseif MATYP == 213
      % Beam elasto-plastic behaviour
      nb_mat_213 = nb_mat_213 + 1;
      mat_213_id(nb_mat_213) = IDMAT ;
      mat_213_E(nb_mat_213) = str2double(line(1:10));
      mat_213_nu(nb_mat_213)= str2double(line(11:20));
      mat_213_rho(nb_mat_213)= RHO;

     else
      disp('MATER, MATYP:')
      disp(MATYP)
      error('MATER: MATYP not yet handled')
     end

    case 'PART  '
     IDPRT = str2double(line(9:16));
     ATYPE = menage_str(line(17:24));
     IDMAT = str2double(line(25:32));
     line=nextLine(fid) ;
     NAME = line(6:length(line)) ;
     %disp(['IDPRT = ',num2str(IDPRT),'  ATYPE = ',ATYPE, ...
     %      '  IDMAT = ',num2str(IDMAT), '  NAME = ',NAME])
     if strcmp(ATYPE(1:4),'BEAM')
      nb_Part_Beam = nb_Part_Beam +1;
      Part_Beam_IDPRT(nb_Part_Beam) = IDPRT ;
      Part_Beam_IDMAT(nb_Part_Beam) = IDMAT ;
      for ii=1:7
       line=nextLine(fid) ;
      end
      Part_Beam_IDSEC(nb_Part_Beam) = str2double(line(1:5));
      Part_Beam_CA(nb_Part_Beam) = str2double(line(11:20));
      if Part_Beam_IDSEC(nb_Part_Beam) ~= 2
       Part_Beam_B(nb_Part_Beam) = str2double(line(21:30));
       if Part_Beam_IDSEC(nb_Part_Beam) == 3
        Part_Beam_C(nb_Part_Beam) = str2double(line(31:40));
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
     i_frame = i_frame + 1;
     frame_id(i_frame) = str2double(line(9:16));
     IFRATY = str2double(line(17:24));
     IAXIS = str2double(line(25:32));
     if IFRATY ~= 1
      error('Frames : IFRATY ~= 1 not yet handled')
     end
     if IAXIS ~= 0
      error('Frames : IAXIS ~= 0 not yet handled')
     end
     line=nextLine(fid) ;
     NAME = menage_str(line(5:length(line)));
     line=nextLine(fid) ;
     Ux = str2double(line(9:24));
     Uy = str2double(line(25:40));
     Uz = str2double(line(41:56));
     line=nextLine(fid) ;
     Vx = str2double(line(9:24));
     Vy = str2double(line(25:40));
     Vz = str2double(line(41:56));
     frame_inod(i_frame) = str2double(line(57:64));
     frame_basis(i_frame,:,:) = basis_UV([ Ux, Uy, Uz ],[ Vx, Vy, Vz ]) ;

    case 'BOUNC '
     i_bounc = i_bounc + 1 ;
     IDNOD = str2double(line(9:16));
     bounc_XYZUVW(i_bounc) = str2double(line(19:24));
     bounc_IFRA(i_bounc) = str2double(line(25:32));
     ISENS = str2double(line(33:40));
     if ISENS ~= 0
      warning('BOUNC : ISENS ~= 0 not yet used')
     end
     line=nextLine(fid) ;
     bounc_names(i_bounc) = line(6:length(line)) ;
     %NAME = strcat(num2str(i_bounc), ':',NAME) ;
     if IDNOD == 0
      line=nextLine(fid) ;
      bounc_lst_nodes(i_bounc) = {read_nodesID(line,fid)} ;
     else
      bounc_lst_nodes(i_bounc) = {[IDNOD]} ;
     end

    case 'MTOCO '
     i_mtoco = i_mtoco + 1 ;
     mtoco_id(i_mtoco) = str2double(line(9:16));
     mtoco_IDNODi(i_mtoco) = str2double(line(17:24));
     mtoco_XYZUVW(i_mtoco) = str2double(line(27:32));
     mtoco_IFRA1(i_mtoco) = str2double(line(33:40));
     line=nextLine(fid) ;
     mtoco_names(i_mtoco) = line(6:length(line)) ;
     line=nextLine(fid) ;
     %lst_nodes = read_nodesID(line,fid) ;
     mtoco_lst_nodes(i_mtoco) = {read_nodesID(line,fid)} ;

     %sdtBC = esi2sdt_BC(XYZUVW) ;
     %mdl=fe_case(mdl,'rigid',NAME,[IDNODi sdtBC lst_nodes]);

    otherwise
     disp(['    Untreated key: ',line(1:7)])
   end
  end
 end

 fclose(fid);
 t_second = toc;

 disp(['> Reading time: ',num2str(t_second),' s'])
 disp('> Summary of read data:')
 if nb_Part_Solid  > 0 ; disp(['    nb Part Solid  : ',num2str(nb_Part_Solid)]) ; end
 if nb_Part_Beam   > 0 ; disp(['    nb Part Beam   : ',num2str(nb_Part_Beam)]) ; end
 if nb_Part_Spring > 0 ; disp(['    nb Part Spring : ',num2str(nb_Part_Spring)]) ; end
 if nb_nod         > 0 ; disp(['    nb nodes       : ',num2str(nb_nod)])  ; end
 if nb_T4          > 0 ; disp(['    nb T4          : ',num2str(nb_T4)])  ; end
 if nb_T10         > 0 ; disp(['    nb T10         : ',num2str(nb_T10)]) ; end
 if nb_spring      > 0 ; disp(['    nb spring      : ',num2str(nb_spring)]); end
 if nb_mass        > 0 ; disp(['    nb mass        : ',num2str(nb_mass)]) ; end
 if nb_mat_1       > 0 ; disp(['    nb mat_1       : ',num2str(nb_mat_1)]) ; end
 if nb_mat_225     > 0 ; disp(['    nb mat_225     : ',num2str(nb_mat_225)]) ; end
 if nb_mat_213     > 0 ; disp(['    nb mat_213     : ',num2str(nb_mat_213)]) ; end
 if nb_BOUNC       > 0 ; disp(['    nb BOUNC       : ',num2str(nb_BOUNC)]) ; end
 if nb_MTOCO       > 0 ; disp(['    nb MTOCO       : ',num2str(nb_MTOCO)]) ; end
 if nb_frame       > 0 ; disp(['    nb FRAME       : ',num2str(nb_frame)]) ; end


 %% SDT model filling
 disp('> SDT model filling')

 %% Node declaration
 mdl.Node = zeros(nb_nod,7) ;
 for iNod=1:nb_nod
  mdl.Node(iNod,1) = node_id(iNod) ;
  mdl.Node(iNod,5:7) = node_xyz(iNod,1:3);
 end
 idx = 1:nb_nod;
 usr2idx = containers.Map(mdl.Node(:,1),idx) ;


 %% Association element / materiau (esi uses PART)
 % Computation of nbColElt
 nbColElt = 0 ;
 if nb_T4 > 0
  nbColElt = max(nbColElt,4+3) ;
  for iElt=1:nb_T4
   for iPart=1:nb_Part_Solid
    if tet4_part(iElt) == Part_Solid_IDPRT(iPart)
     tet4_mat(iElt) = Part_Solid_IDMAT(iPart) ;
    end
   end
  end
 end
 if nb_T10 > 0
  nbColElt = max(nbColElt,10+3) ;
  for iElt=1:nb_T10
   for iPart=1:nb_Part_Solid
    if tet10_part(iElt) == Part_Solid_IDPRT(iPart)
     tet10_mat(iElt) = Part_Solid_IDMAT(iPart) ;
    end
   end
  end
 end
 if nb_beam > 0
  nbColElt = max(nbColElt,9) ;
  for iElt=1:nb_beam
   for iPart=1:nb_Part_Beam
    if beam_part(iElt) == Part_Beam_IDPRT(iPart)
     beam_mat(iElt) = Part_Beam_IDMAT(iPart) ;
    end
   end
  end
 end

 %% Element declaration
 %  [Element_type_declaration ;
 %   NodeNumbers MatId ProId EltId OtherInfo ]
 %  MatId material ID number in model.pl
 %  ProId element property (e.g. section area) ID number in model.il

 row_elt = 0 ;

 if nb_T4 > 0
  mdl.Elt=[mdl.Elt; zeros(nb_T4+1,nbColElt)] ;
  row_elt = row_elt + 1 ;
  mdl.Elt(row_elt,1:7)=[Inf  abs('tetra4') ] ;
  for iElt=1:nb_T4
   mdl.Elt(row_elt+iElt,1:(4+3)) = ...
    [ tet4(iElt,:), ...
    tet4_mat(iElt), ...    % MatId <= IDMAT
    tet4_part(iElt), ...   % ProId <= IDPRT
    tet4_id(iElt) ];
  end
  row_elt = row_elt+nb_T4 ;
 end

 if nb_T10 > 1
  mdl.Elt=[mdl.Elt; zeros(nb_T10+1,nbColElt)] ;
  row_elt = row_elt + 1 ;
  mdl.Elt(row_elt,1:8)=[Inf  abs('tetra10') ] ;
  for iElt=1:nb_T10
   mdl.Elt(row_elt+iElt,1:(10+3)) = ...
    [ tet10(iElt,:), ...
    tet10_mat(iElt), ...    % MatId <= IDMAT
    tet10_part(iElt), ...   % ProId <= IDPRT
    tet10_id(iElt)];
  end
  row_elt = row_elt+nb_T10 ;
 end

 if nb_beam > 0
  mdl.Elt=[mdl.Elt; zeros(nb_beam+1,nbColElt)] ;
  row_elt = row_elt + 1 ;
  mdl.Elt(row_elt,1:6)=[Inf  abs('beam1') ] ;
  for iElt=1:nb_beam
   mdl.Elt(row_elt+iElt,1:(4+3)) = ...
    [ beam(iElt,1:2), ...
    beam_mat(iElt), ...    % MatId <= IDMAT
    beam_part(iElt), ...   % ProId <= IDPRT
    beam(iElt,3), ...      % node defining bending plane 1
    0, 0 ];
  end
  %row_elt = row_elt+nb_beam ;
 end
 %disp(mdl.Elt)


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
    basis = reshape(frame_basis(bounc_IFRA(i_bounc),:,:),3,3) ;
   end
   mdl=add_BOUNC(mdl, bounc_names(i_bounc), bounc_XYZUVW(i_bounc), ...
    lst_nodes_F, basis);
  end
 end

 for i_mtoco=1:nb_MTOCO
  if mtoco_IFRA1(i_mtoco) == 0
   basis = eye(3) ;
  else
   basis = reshape(frame_basis(mtoco_IFRA1(i_mtoco),:,:),3,3) ;
  end
  mdl=add_MTOCO(mdl, mtoco_names(i_mtoco), mtoco_IDNODi(i_mtoco), mtoco_XYZUVW(i_mtoco), ...
   cell2mat(mtoco_lst_nodes(i_mtoco)), basis, ...
   usr2idx);
 end

 disp('> Model filling completed')

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
function line=nextLine(fid, nb_line)
%% #nextLine--------------------------------------------------------------2
if nargin == 1
 nb_line = 1 ;
end
for il=1:nb_line
 while 1
  line = fgetl(fid);
  if line(1) ~= '$', break, end
 end
end
end

function lst_nodesID = read_nodesID(line, fid)
%% #read_nodesID----------------------------------------------------------2
lst_nodesID = [];
nbNodes = 0 ;
while strcmp(line(1:11),'        NOD')
 first_split = strsplit(line(12:length(line)),' ');
 for i1 = 1:length(first_split)
  split_elem = char(first_split(i1)) ;
  if length(split_elem) > 0
   sgnd_split=strsplit(split_elem,{':'});
   if length(sgnd_split) == 1
    nbNodes = nbNodes + 1 ;
    lst_nodesID(nbNodes) = str2double(char(sgnd_split(1)));
   else
    for inod=str2double(char(sgnd_split(1))):str2double(char(sgnd_split(2)))
     nbNodes = nbNodes + 1 ;
     lst_nodesID(nbNodes) = inod ;
    end
   end

  end
 end
 line=nextLine(fid) ;
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
mdl=fe_case(mdl,'mpc',[char(NAME),'_bc'],data);

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
data=struct('DOF',lst_dof,'c',vertcat(c{:}));
mdl=fe_case(mdl,'mpc',[char(NAME),'_mpc'],data);
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