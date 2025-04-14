function sdt_anim(block)

%
%
%
% Based on Level-2 MATLAB file S-Function for continuous time variable step demo.
%   Copyright 1990-2025 The MathWorks, Inc.

  setup(block);
  
%endfunction

function setup(block)
  
  %% Register number of input and output ports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  ref=evalin('base','ref');
  block.InputPort(1).Dimensions        = size(ref.c,1);
  block.InputPort(1).DirectFeedthrough = true;
  
  block.OutputPort(1).Dimensions       = size(ref.b,2)-1;
  
  %% Set block sample time to fixed sample time
  block.SampleTimes = [.01 0];
  
  %% Set the block simStateCompliance to default (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  %% Register methods
  block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
  block.RegBlockMethod('InitializeConditions',    @InitConditions); 
  block.RegBlockMethod('Outputs',                 @Fnl);  
  %block.RegBlockMethod('Update',                  @Update);  

%endfunction

function DoPostPropSetup(block)

  %% Setup Dwork
  ref=evalin('base','ref');
  block.NumDworks = 1;
  block.Dwork(1).Name = 'tu'; 
  block.Dwork(1).Dimensions      = size(ref.c,1);
  block.Dwork(1).DatatypeID      = 0;
  block.Dwork(1).Complexity      = 'Real';
  block.Dwork(1).UsedAsDiscState = true;

  
%endfunction

function InitConditions(block)

  %% Initialize state to zero
  block.Dwork(1).Data = zeros(size(block.Dwork(1).Data));
  
%endfunction

function Fnl(block)

  ref=evalin('base','ref');
  if isfield(ref,'fnl')
    ref.fnl.c=reshape(ref.b(:,2),[],2)';
    y=ref.fnl.c*block.InputPort(1).Data; % Obersvation collocated with input1
    fnl=ref.fnl.fun(y); % Simple case of a cubic load
    block.OutputPort(1).Data=fnl;
  end

  %block.OutputPort(1).Data=
  % block.Dwork(1).Data;
  
  %% Set the next hit for this block 
  %block.NextTimeHit = block.CurrentTime + block.InputPort(1).Data(2);
  
%endfunction

function Update(block)

% dbstack; keyboard
% block.Dwork(1).Data = block.InputPort(1).Data(1);
  
%endfunction

