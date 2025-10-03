(TMW)=
# Tracking of support / enhancement requests to The MathWorks


(TMW_JavaJS)=
## Tracking of the transition of the SDT GUI from java to JS based rendering (MATLAB>=2025a)


(TMW_uitable)=
## `vhandle.tab` rendering as `uitable`


Regressions/roadblocks
 - 😠 How to implement cell based tooltips ? 
 - How to analyze rendering performance on the JS side ?
- UseCase 
  - [TagList](https://www.sdtools.com/helpcur/base/sdtweb.html#_taglist), ChannelTab 

TODO

 - TableCellRenderer
   - get TMW to explicit strategies for efficient partial update table 
 - `vhandle.tab.jsSafeCell`
   - `pop` @TMW says it is possible, but I only find how to do pop on a column by column basis (not a cell by cell)
   - `push` improve decoration
   - `level` column of `TreeTable` is shown as decorations in the first column. Current implementation is in `vhandle.tab.level2iconkey`. `reshape(char(9660+[-99:100]),10,[])` char(9660) <a style="font-family:Courier;font-size:20px;">:▼:☹:☺:</a> 
   - Microsoft uses private area U+F04A
   - Monospaced font selection for part of the row. 
   - Need to clean all reference to javaframe 
   -  consider switching icons to svg : `'@matlabroot\ui\icons\16x16\yAxisView.svg'`

DONE

 - In Java [`TableCellRenderer`](https://docs.oracle.com/javase/8/docs/api/javax/swing/table/TableCellRenderer.html). Writing a Matlab based `cellRenderer` is a major performance killer and will not occur. 
   - `vhandle.tab.asUitable` will deal with filling a `uitable.Data` peer containing the view model
   - since the `uitable` will contain a view of the data rendering of buttons will be done by html styling in `[text,style]=vhandle.tab.jsSafeCell` 
 - `vhandle.tab.jsSafeCell` need to render the basic button types
   - `text/double` : do nothing. @TMW Is there a performance issue if `ta.Data` a contains mixed type cell array ? 
   - `push` do a decoration around text. Currently MATLAB does not support inline reference to icons (they are placed in the style file for an unknown HTML rendering limitation @TMW it would be useful to have a clear answer of whether this will stay or go) -> use unicode/emoji 
   - `pop` @TMW said it was possible to have cell based popup menus, waiting for details 
  

 uit.DisplayDataChangedFcn = @(src,event) updatePlot(src,ax);


(TMW_tabbedpane)=
## Rendering of `TabbedPane` as `tabgroup`


Regressions/roadblocks
 - @TMW Missing closing button 
 
TODO
 - port use of `uitab` with `javacomponent` for the viewport rather than a `TabbedPane`

DONE 
 - URN for tabbed pane group inclusion `figure(1).Tab`
 - `uf=sdth.urn('figure(ui.Tree).uitab{Tree}');`
 - `t_js('uitree')` is basic test 


 (TMW_dock)=
## Recovering a docking framework `SplitPane` as `uigridlayout`

Regressions/roadblocks
 - no documentation of multiple docks (dock name associated with groups )
 - documentation of programmatic splitting ? 

TODO
 - continue revising `sdtu.ui.dock`


DONE
 - added compatibility layer for >=2025a



 (TMW_splitpane)=
## Rendering of `SplitPane` as `uigridlayout`

Regressions/roadblocks
 - Missing interactive resize of grid width

TODO
 - AlwaysOnTop (bypass is probably possible with modal ? )
 - Dock management (more than the single Figure dock is necessary for us (we need them to guide users in processes that occur). I will work on our transition prototype and come back to you for guidance. 



(TMW_errors)=
## Miscellaneous annoyance ☹ / errors 😠

- ☹ `quit` does not work for us in either 2025a or 2025b, we need to use `quit('force')`. Suggestion on how to diagnose the problem ? 
- Annoying caught error

```
figure(1);plot([0 1]); dbstop if caught error;print(1,'-noui','-dpng','xxx.png')

Caught-error breakpoint was hit in printing\private\setGraphicsProperty at line 0. The error was:
Error using matlab.graphics.axis.Axes/set
Unrecognized property YColorMode_I for class Axes.

```

- ☹ in the web browser there is no limit on the text length shown for a tab (makes using multiple tabs a pain). 
- ☹ setFigureIcon (a cosmetic issue for our customers selling products developed with our toolbox)
- how to use git to track a live script ? 
- ☹ named HTML color support (not a good idea that we are the only ones to support it)
- ☹ why to we have to do bypasses to open file at specific line `matlab.desktop.editor.Document.goToLine` ?