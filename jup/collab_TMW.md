(TMW)=
# Tracking of support / enhancement requests to The MathWorks


(TMW_JavaJS)=
## Tracking of the transition of the SDT GUI from java to JS based rendering (MATLAB>=2025a)

Smiley convention : annoyance ☹ / errors-roadblocks 😠 / fine 😊

This was discussed in [ ref:!00Di00Ha1u.!500UU0V14ir:ref ] and on the Matlab exchange https://fr.mathworks.com/matlabcentral/answers/2179688-treetable-in-matlab-2025a. 
 
(TMW_uitable)=
## `vhandle.tab` rendering as `uitable`


Regressions/roadblocks
 - 😠 How to implement cell based tooltips ? 
   - @TMW.mt may occur in the future
   - Current bypass thought use cell extention on focus capability (but will fail for pop/categorical/edit)
 - ☹ How to analyze rendering performance on the JS side ?
- Use cases 
  - [TagList](https://www.sdtools.com/helpcur/base/sdtweb.html#_taglist), ChannelTab, [FEMLink](https://www.sdtools.com/helpcur/base/sdttab.html#tabfemlink)

TODO

 - TableCellRenderer
   - get @TMW.rm or @TMW.mt to explicit strategies for efficient partial update table 
 - `vhandle.tab.jsSafeCell` deals with rendering of [CinCell](https://www.sdtools.com/helpcur/base/gui_data.html#CinCell)
   - `push` improve decoration 😊. ❓😠 is it possible to use css styles to decorate part of text in a cell ?
   -  consider switching icons to svg : `'@matlabroot\ui\icons\16x16\yAxisView.svg'`
   -  ❓@TMW.mt why support it in toolbars and not uitable? 
-  [uistyle Interpreter=html ](https://fr.mathworks.com/help/releases/R2025b/matlab/ref/uistyle.html#mw_e57c82dd-e4a5-4f48-a8db-0561c0ee0341), tested in {m}`t_js('TutoTMW251001 -s{componentRenderer_proto} -show4')`
   -  😊 overal
   -  need progress on unicode. 
   -  ❓😠 how to do cell based tooltip ?
   -  ❓😠 how to mix fonts ? it seems possible in the latex side. Monospaced font selection for part of the row. 
   -  ❓😠 : do we need to trash the idea that inline SVG will be supported in table cells ? 
   - ❓😠 how can we get out of the loop that TMW takes a long time to answer : see code of uifigure https://matlabthoughts.com/2022/11/10/web-figures-uifigures-and-uihtml/
   


DONE

 - In Java [`TableCellRenderer`](https://docs.oracle.com/javase/8/docs/api/javax/swing/table/TableCellRenderer.html), transition to a SDT view model
   - @TMW.rm writing a Matlab based `cellRenderer` is a major performance killer and will not occur. 
   - {m}`vhandle.tab.asUitable` will deal with filling a `uitable.Data` peer containing the view model. `vhandle.tab` stored as 'sdt' appdata in uitable. Row model is stored in `vh.GHandle.irow` 
   - since the `uitable` will contain a view of the data rendering of buttons done by html styling in `[text,style]=vhandle.tab.jsSafeCell` 
   - Need to work on non-synchronous filling, sdtu.logger.do
 - `TreeTable` implementation
   - ⚠️ `level` column of `TreeTable` is shown as decorations in the first column. Current implementation is in `vhandle.tab.level2iconkey`.
   - ⛔ : using icons trashed since TMW does not give a time frame
 - {m}`sdtweb t_js vhandle.tab.jsSafeCell` tests rendering of the basic button/cell types
   - `text/double` : do nothing. @TMW Is there a performance issue if `ta.Data` a contains mixed type cell array ? 
   - `push` do a decoration around text. Currently MATLAB does not support inline reference to icons (they are placed in the style file for an unknown HTML rendering limitation @TMW it would be useful to have a clear answer of whether this will stay or go) -> use unicode/emoji 
   - `pop` example of categorical  
  
  
```matlab
 % Possibly useful code 
 uit.DisplayDataChangedFcn = @(src,event) updatePlot(src,ax);

%% @TMW.mt Cell font (☹ which is not multi fonts in Cell)
t = uitable(uifigure, 'Data', {'First','Second'});

% Define different styles for each cell
style1 = uistyle('FontName', 'Arial', 'FontColor', 'red');
style2 = uistyle('FontName', 'Courier New', 'FontColor', 'blue');

% Apply styles to individual cells
addStyle(t, style1, 'cell', [1 1]); % Cell (1,1)
addStyle(t, style2, 'cell', [1 2]); % Cell (1,2)

%% Pop>Categorical 
% sdtweb t_js categorical

Colors = categorical({'Red';'Blue';'Green';}); 
Numbers = categorical({'1';'2';'3';}); 
Names = categorical({'John';'Jane';'Ace'}); 
Letters =  categorical({'A';'B';'C'}); 

gf=uifigure;
t = uitable(gf, 'Data', table({Colors(1); Numbers(1)},{ false;Letters(1)},{ Names(1);1.3}), 'ColumnEditable', true);
% Note strangely : multiple columns is not the same as cell array
ta = uitable(gf, 'Data', table([{Colors(1); Numbers(1)},{ false;Letters(1)},{ Names(1);1.3}]), 'ColumnEditable', true);


```

(TMW_tabbedpane)=
## Rendering of `TabbedPane` as `tabgroup`


Regressions/roadblocks
 - ☹ @TMW Missing closing button 
 - Use case {m}`t_js('uitree')` 

TODO
 - port use of `uitab` with `javacomponent` for the viewport rather than a `TabbedPane`
 - Report issues to TMW [ ref:!00Di00Ha1u.!500UU0WS8zX:ref ]
 - check ability to capture content of a single pane. 

DONE 
 - URN for tabbed pane group inclusion `figure(1).Tab`
 - {m}`uf=sdth.urn('figure(ui.Tree).uitab{Tree}');`
 - 
 

 (TMW_dock)=
## Recovering a docking framework 

Regressions/roadblocks
 - ☹ no documentation of multiple docks (dock name associated with groups )
 - documentation of programmatic splitting ? 

TODO
 - continue revising {m}`sdtu.ui.dock`
 - ☹❓Port [robot capture methods](https://www.sdtools.com/helpcur/base/comgui.html#ImWrite)


DONE
 - added compatibility layer for >=2025a
 - AlwaysOnTop (bypass is probably possible )



 (TMW_splitpane)=
## Rendering of `SplitPane` as `uigridlayout`

Regressions/roadblocks
 - 😠Missing interactive resize of grid width

TODO
 - Dock management (more than the single Figure dock is necessary for us (we need them to guide users in processes that occur). I will work on our transition prototype and come back to you for guidance. 

DONE


(TMW_errors)=
## Miscellaneous annoyance ☹ / errors 😠

- ☹ `quit` does not work for us in either 2025a, 2025b, 2026a, we need to use `quit('force')`. Suggestion on how to diagnose the problem ? 

- ☹ in the web browser there is no limit on the text length shown for a tab (makes using multiple tabs a pain). 
- ☹ setFigureIcon (a cosmetic issue for our customers selling products developed with our toolbox)
- how to use git to track a live script ? 
- ☹ named HTML color support (not a good idea that we are the only ones to support it)
- ☹ inability to have a stable link to [uistyle/interpreter/html](https://fr.mathworks.com/help/releases/R2025b/matlab/ref/uistyle.html#mw_e57c82dd-e4a5-4f48-a8db-0561c0ee0341). ToDo SDTools implement `@MatlabHelp`


- Error using fopen The 'all' option will be removed in a future release. For a list of all open file identifiers, use the openedFiles function
instead. You are saying will be removed when it aready has. 


Solved : 
- ☹ `h=matlab.desktop.editor.Document.findEditor(which('evt.WhichName');` is not always the same as `h=matlab.desktop.editor.Document.findEditor(evt.WhichName);`. Now switched to using `opentoline`. ☹ why do we have to do bypasses to open file at specific line {m}`matlab.desktop.editor.Document.goToLine` ?
- ☹ Annoying caught error (solved 2026a)

```
figure(1);plot([0 1]); dbstop if caught error;print(1,'-noui','-dpng','xxx.png')

Caught-error breakpoint was hit in printing\private\setGraphicsProperty at line 0. The error was:
Error using matlab.graphics.axis.Axes/set
Unrecognized property YColorMode_I for class Axes.

```
