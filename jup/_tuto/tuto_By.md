(tuto:By)=
# Experimental Modal Analysis `By` block
The previous tutorial [](tutoEMA) has focused on the identification of 
 SIMO or MIMO transfers. 
 
 For some test cases, it is necessary to perform several identification by
 blocks because pole shifts are expected between these blocks.
 
 Here are examples of typical parameters defining the blocks :
 - Temperature
 - Input location
 - Input level
 - Sensor batches
 - Pre-loading
 - ...
 
 SDT provides tools to easily navigate between blocks, perfom the 
 identification for each block and finally post-treat the result 
 (pole tracking, clustering, main shape extraction,...)
 
 The following tutorial highlights the procedure through the parametric 
 identification of the `DemoGartBy` example in which each block of 
 transfers depends on the temperature.
## Step 1 : Load wireframe and parametric transfers

:::{dropdown} <a href="matlab:gartid('tutoBy-s1;')"><img src="../../_images/run16.png" >Run step </a>. Expand to display source code.
```matlab
[XF,wire]=demosdt('DemoGartBy');
[ci,cf]=iicom('dockid',struct('model',wire,'XF',XF));

```
:::
Load the curve containing all the measurement blocks and eventually the 
 wireframe. 

 Store them in the dock Id with the command
 ```matlab
 % Example of measurement blocks with GARTEUR
 [XF,wire]=demosdt('DemoGartBy'); 
 % Open dock with measurement blocks in curve 'Test'
 [ci,cf]=iicom('dockid',struct('model',wire,'XF',XF)); 
 ```
:::{figure} ../_images/TutoBy_1.png
:width: 100%
:name: TutoBy_1
:::
## Step 2 : Switch to By mode

:::{dropdown} <a href="matlab:gartid('tutoBy-s2;')"><img src="../../_images/run16.png" >Run step </a>. Expand to display source code.
```matlab
% Simulate click on button "Init By Mode"
iicom(ci,'SetIdent',struct('ByMode','do'))



```
:::
To initialize the identification by block, on tab `Ident`, expand 
 `IDopt` and click on `Init By Mode`.
:::{figure} ../_images/TutoBy_2a.png
:width: 50%
:name: TutoBy_2a
:::
In the `Ident` tab, a block of buttons `By` has appeared.
:::{figure} ../_images/TutoBy_2b.png
:width: 50%
:name: TutoBy_2b
:::
These button help in the navigation between block through all the 
 procedure of identification by blocks. All the blocks are listed in the
 pop button on the top right (here 11 temperatures). You can increase or 
 decrease the selected block by clicking on the + or - buttons or directly 
 selecting the wanted one in the list. When one block is selected, 
 measurement corresponding to this block only is stored in the curve Test 
 and eventually previous identification results.