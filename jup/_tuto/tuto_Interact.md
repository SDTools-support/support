(tuto:Interact)=
# `Feplot` interactivity
This tutorial aims at illustrating the main ways to interact with `feplot`
 figures.
 <iframe width="560" height="315" 
 src="https://www.youtube.com/embed/xET78h-wNPY?si=BXWMqd0CVAqH149m" 
 title="YouTube video player" frameborder="0" allow="accelerometer; 
 autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; 
 web-share" allowfullscreen></iframe>

## Step 1 : Initalize feplot with GARTEUR test case

:::{dropdown} <a href="matlab:d_feplot('tutoInteract-s1;')"><img src="../_images/run16.png" >Run step </a>. Expand to display source code.
```matlab
cf=demosdt('DemoGartDataFEM -plot');


```
:::
Run this step to load in a feplot figure the GARTEUR test case.
 This test case contains a very light FEM model and the computed mode
 shapes. 
 We will use this test case to illustrate
 - camera controls
 - shape animation and navigation between shape channels
 - images and movies generation
 - navigation into the feplot property figure
 
 **Note that you can open the list of all interaction with the `feplot` 
 figure with the key `?` when the `feplot` figure is the current figure**
## Step 2 : Camera controls

:::{dropdown} <a href="matlab:d_feplot('tutoInteract-s2;')"><img src="../_images/run16.png" >Run step </a>. Expand to display source code.
```matlab

```
:::
```{list-table}
 :header-rows: 1
 :widths: 15 20 40
 * - Purpose
   - Direction
   - Action / Shortcut
 * - Rotation
   - Around axis perpendicular to the sight axis and the mouse displacement
   vector
   - `Ctrl + Drag`
 * - 
   - Around u (horizontal axis)
   - `u + Drag` `u + Scroll` `u` `U`
 * - 
   - Around v (vertical axis)
   - `v + Drag` `v + Scroll` `v` `V`
 * - 
   - Around w (sight axis)
   - `w + Drag` `w + Scroll` `w` `W`
 * - Zoom
   - Along the sight
   - `Ctrl + Scroll` `a + Drag` `a + Scroll` `a` `A`
 * - 
   - To the area in the box
   - `Drag`    
 * - Translation
   - Following the mouse displacement
   - `Shift + Drag`
 * - 
   - Through u (horizontal axis)
   - `x + Drag` `x + Scroll` `x` `X`
 * - 
   - Through v (vertical axis)
   - `y + Drag` `y + Scroll` `y` `Y`
 * - 
   - Through w (sight axis)
   - `z + Drag` `z + Scroll` `z` `Z`    
 * - Default views
   - (XZ) plane
   - `1`
 * - 
   - (XY) plane
   - `2`
 * - 
   - "diagonal" view
   - `3`
 * - 
   - (YZ) plane
   - `4`
 * - ResetView
   -
   - `DoubleClick` `i`
 * - Flip view
   - From behinds
   - `f`
 ``` 
## Step 3 : Shape animation + channel navigation

:::{dropdown} <a href="matlab:d_feplot('tutoInteract-s3;')"><img src="../_images/run16.png" >Run step </a>. Expand to display source code.
```matlab

```
:::
```{list-table}
 :header-rows: 1
 :widths: 15 60
 * - Purpose
   - Action / Shortcut
 * - Animation
   - ![an.png](../_icons/an.png) `Shift + Scroll` `UpArrow` `DownArrow` 
 * - Channel
   - ![bplus_plus.png](../_icons/bplus_plus.png) 
     ![bplus_minus.png](../_icons/bplus_minus.png)
     `Scroll` `+` `-` `RightArrow` `LeftArrow` `Home` `End`
 * - Deformation scale
   - `l + Scroll` `l + Drag` `l` `L` 
 ``` 
## Step 4 : Feplot capture

:::{dropdown} <a href="matlab:d_feplot('tutoInteract-s4;')"><img src="../_images/run16.png" >Run step </a>. Expand to display source code.
```matlab

```
:::
```{list-table}
 :header-rows: 1
 :widths: 15 60
 * - Purpose
   - Action / Shortcut
 * - Capture image
   - ![an.png](../_icons/snapshot.png)
 * - Generate movie
   - ![bplus_plus.png](../_icons/movie16.png) 
 ``` 
 The captured image and movie are shown below :
:::{figure} ../_images/TutoFeplotInteract_1.png
:width: 50%
:name: TutoFeplotInteract_1
:::
:::{figure} ../_images/TutoFeplotInteract_2.gif
:width: 50%
:name: TutoFeplotInteract_2.gif
:::
## Step 5 : Interaction with the model / analysis

:::{dropdown} <a href="matlab:d_feplot('tutoInteract-s5;')"><img src="../_images/run16.png" >Run step </a>. Expand to display source code.
```matlab
%%

```
:::
The list below provides the entry points to analyze the model and
 customize the rendering.
 ```{list-table}
 :header-rows: 1
 :widths: 15 60
 * - Purpose
   - Action / Shortcut
 * - Figure menu
   - `File` : mainly to save/load data

     `Feplot` : rendering options (color, animation type,...)
 * - Context menu
   - `RightClick` on the axes background
 * - Node cursor
   - `n`
 * - Model properties
   - ![mprop.png](../_icons/mprop.png)