```{include} ../header.md
```

```{tableofcontents}
```

(collab-install)=
# Installation / upgrade

(collab-install-MATLAB+SDT)=
## MATLAB+SDT
With MATLAB you need 
- SDT : www.sdtools.com/sdtcur
- {m}`sdtcheck 'patchupgrade ../download/daily/sdtdaily_dis.p'` is used to update your installation based on the daily post of SDT. Use this if you are running a project with SDTools only. 
- {m}`sdtcheck 'patchupgrade'` update from the latest beta release located at https://www.sdtools.com/distrib/beta/sdtcur.zip  

(collab-install-Python+MATLAB)=
## Python/MATLAB integration for use with SDT
- Install a Python compatible with MATLAB  https://www.mathworks.com/support/requirements/python-compatibility.html,  {m}`doc 'configure your system to use python'`  
Latest windows binary compatible with MATLAB (Python 3.9) is found here : [Python Release Python 3.9.13 | Python.org](https://www.python.org/downloads/release/python-3913/)  
Assign which python executable MATALB should use with command  
```MATLAB
pyenv('Version','_full_path_to_python.exe_')
```

(collab-install-JupyterBook)=
## Build Jupyter book documentation

If using Python / Jupyter book to generate this documentation

- First setup Python for MATLAB following instructions in 
- Install Jupyter book for Python by executing the commands listed by the call below in a Windows batch  
{m}`sdtpy.install('jup')`
- Update the matlab lexer used by `Pygments` to color code to add SDT flavor.  
Use command {m}`sdtm.jup('sdtlexer')`
- Create symbolic links to organize your directories
  - {m}`mklink /d sdt.cur\_jup   target` where {m}`target` is typically tempdir\_jup 
  - {m}`mklink /d target\base sdt.cur\support\jup` to link files of the SDT support as base. 
  - {m}`mklink /d target\proj sdt.git\proj\jup` to add any project to the building tree. 
  - within the  `sdt.git\proj\jup directory` use {m}`mklink /d _images ..\tex\plots`  to link the expected jup/_images to actual location
 `mklink /d _icons ..\tex\icons`  if appropriate


Other Python installs of interest

- {m}`sdtpy.install('scipy')` lists command to install SciPy. Verify `import scipy; scipy.signal.butter(1,.1,"lowpass")`
- {m}`sdtpy.install('gmsh')` lists command to install GMSH. Verify using `import gmsh; gmsh.initialize(); gmsh.fltk.run(); gmsh.finalize()`
	
(collab-install-VSCode)=
## VSCode as Markdown editor

- VSCode (Visual Studio Code) is a modern open source code editor
- `Ctrl+k  v` opens the preview editor side by side 
- The explorer allows nice navigation with multiple files if you define a workspace. For example sdt.code-workspace contains

```
{"folders": [
		{"path": "support/jup"},{"path": "test2/jup"}
	]}
```

(collab-install-GitHub)=
## GitHub desktop as git interface

- Github desktop is an intuitive interface for Git access. It is the choice used by SDTools to explain git to users that are not familiar with this versioning system.   

