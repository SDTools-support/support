(tuto:ExpGUI)=
# Expansion GUI
This tutorial aims at illustrating how to implement and analyze expansion
 technique using GUI in SDT.
## Step 1 : Load required data for expansion

:::{dropdown} <a href="matlab:d_cor('tutoExpGUI-s1;')"><img src="../_images/run16.png" >Run step </a>. Expand to display source code.
```matlab

[model,def,id]=demosdt('DemoGartDataCoShape');

```
:::
Execute the line below to load all required data for expansion :
 ```MATLAB
 [model,def,id]=demosdt('DemoGartDataCoShape');
 ```
 - `model` : GARTEUR model (very light FEM model classically used in SDT tutorials)
that also contains the sensor definition and wireframe geometry in the 
`SensDof` Stack entry
 - `def` : FEM modes
 - `id` :  identified modes
 :::{note}
 :class: dropdown
 The `SensDof` is used to build the observation matrix of the sensors (tdof) 
 on the FEM (DOF).
 You can access its structure in `wire` and build the observation matrix in
 `Sens`  with the command below
 ```MATLAB
 wire=fe_case('getsensdof',model);
 Sens=fe_case('Sens',model);
 ```
 :::
## Step 2 : Define a model parameter from script

:::{dropdown} <a href="matlab:d_cor('tutoExpGUI-s2;')"><img src="../_images/run16.png" >Run step </a>. Expand to display source code.
```matlab
model=fe_case(model,'ParAdd k 1.0 0.1 10','ConstrainedLayer','pro4');

```
:::
:::{warning}
 :class: dropdown
 Note that this step often needs the module `PARAM` and custom strategies
 when dealing with industrial models. This is a very simple
 parameterization example to show how it can be used later on.
 :::
 Execute the line below to parameterize the ContrainedLayer stiffness in
 the model.
 ```MATLAB
 model=fe_case(model,'ParAdd k 1.0 0.1 10','ConstrainedLayer','pro4');
 ```
 This parameter is introduced here to illustrate how to deal with it during
 expansion.
## Step 3 : Open the dock `CoShape`

:::{dropdown} <a href="matlab:d_cor('tutoExpGUI-s3;')"><img src="../_images/run16.png" >Run step </a>. Expand to display source code.
```matlab
RO=struct('model',model,'va',id,'vb',def);
VC=iicom('DockCoShape',RO);

```
:::
To open the dock `CoShape` with the required data, store them in a `RO` 
 (running option) structure with the fields `model`, `va` and `vb`. 

 `va` reads **V**ector set **A** (here the identified modes) 

 `vb` reads **V**ector set **B** (here the FEM modes)
 ```MATLAB
 RO=struct('model',model,'va',id,'vb',def);
 ```
 Then execute the command below to open the dock `CoShape`
 ```MATLAB
 VC=iicom('DockCoShape',RO);
 ```
## Step 4 : Build the reduced parameterized model for expansion

:::{dropdown} <a href="matlab:d_cor('tutoExpGUI-s4;')"><img src="../_images/run16.png" >Run step </a>. Expand to display source code.
```matlab
% Compute reduced model for expansion
ii_mac(VC,'SetCoExp',struct('sparse',1,'Reduce','Mode+Sens','do','Reduce'));

```
:::
Open the tab `CoExp` to build the reduced parameterized model
 To performe `Minimum Dynamic Residual Expansion`, we will use the
 reduction method `Mode+Sens` which uses mode shapes and an enrichment with
 static response at each sensor.
 Click on ![](../_icons/run16.png)
## Step 5 : Specify design points for expansion computation

:::{dropdown} <a href="matlab:d_cor('tutoExpGUI-s5;')"><img src="../_images/run16.png" >Run step </a>. Expand to display source code.
```matlab
ii_mac(VC,'SetCoExp',struct('jw','1:5','lgamma','1:5',...
'param1','@log(-1,1,11)','MDRE','do'));

```
:::
## Step 6 : Switch between post-expansion analyzes

:::{dropdown} <a href="matlab:d_cor('tutoExpGUI-s6;')"><img src="../_images/run16.png" >Run step </a>. Expand to display source code.
```matlab
ii_mac(VC,'SetMDRE',struct('ShowExp','do'));
ii_mac(VC,'SetMDRE',struct('ShowExpTest','do'));
ii_mac(VC,'SetMDRE',struct('ShowErrMod','do'));
ii_mac(VC,'SetMDRE',struct('ShowErrModModesOnly','do'));
ii_mac(VC,'SetMDRE',struct('ShowErrModSensOnly','do'));
%ii_mac(VC,'SetMDRE',struct('ShowErrTest','do')); % Not useful as it is
%shown in feplot3

```
:::
## Step 7 : Show model error vs. parameter for all $\gamma$

:::{dropdown} <a href="matlab:d_cor('tutoExpGUI-s7;')"><img src="../_images/run16.png" >Run step </a>. Expand to display source code.
```matlab
ii_mac(VC,'SetMDRE',struct('ShowExpTest','do'));
ii_mac(VC,'SetMDRE',struct('Xaxis3',1,'Dim3',1,'Dim2','all','Dim4','ErrMod','Dim5',1));

```
:::