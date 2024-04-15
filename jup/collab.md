```{include} ../header.md
```
(markdown)=
# Collaboration with SDTools

```{tableofcontents}
```

## Installation 

With MATLAB you need 
- SDT : www.sdtools.com/sdtcur
- {m}`sdtcheck 'patchupgrade ../download/daily/sdtdaily_dis.p'` is used to update your installation based on the daily post of SDT. Use this if you are running a project with SDTools only. 
- {m}`sdtcheck 'patchupgrade'` update from the latest beta release located at https://www.sdtools.com/distrib/beta/sdtcur.zip  

If using Python / Jupyter book to generate this documentation

- Install a Python compatible with MATLAB  https://www.mathworks.com/support/requirements/python-compatibility.html,  {m}`doc 'configure your system to use python'`
- {m}`sdtpy.install('jup')` lists command to install Jupyter book. 
- Create symbolic links to organize your directories
  - {m}`mklink /d sdt.cur\_jup   target` where {m}`target` is typically tempdir\_jup 
  - {m}`mklink /d target\base sdt.cur\support\jup` to link files of the SDT support as base. 
  - {m}`mklink /d target\proj sdt.git\proj\jup` to add any project to the building tree. 
  - within the  `sdt.git\proj\jup directory` use {m}`mklink /d _images ..\tex\plots`  to link the expected jup/_images to actual location
 `mklink /d _icons ..\tex\icons`  if appropriate


Other Python installs of interest

- {m}`sdtpy.install('scipy')` lists command to install SciPy. Verify `import scipy; scipy.signal.butter(1,.1,"lowpass")`
- {m}`sdtpy.install('gmsh')` lists command to install GMSH. Verify using `import gmsh; gmsh.initialize(); gmsh.fltk.run(); gmsh.finalize()`
	
