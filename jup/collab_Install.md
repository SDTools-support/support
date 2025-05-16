(collab_install)=
# Installation / upgrade

(collab_install_MATLAB_SDT)=
## MATLAB+SDT
With MATLAB you need 
- SDT : www.sdtools.com/sdtcur
- {m}`sdtcheck 'patchupgrade ../download/daily/sdtdaily_dis.p'` is used to update your installation based on the daily post of SDT. Use this if you are running a project with SDTools only (since daily versions may have issues). 
- {m}`sdtcheck 'patchupgrade ../download/daily/nlsimdaily_dis.p'` may be used if you also have an SDT-nlsim license.  
- {m}`sdtcheck 'patchupgrade'` update from the latest beta release located at https://www.sdtools.com/distrib/beta/sdtcur.zip  

(collab-install-Python+MATLAB)=
## Python/MATLAB integration for use with SDT
- Install a Python compatible with MATLAB  https://www.mathworks.com/support/requirements/python-compatibility.html,  {m}`doc 'configure your system to use python'`  
Latest windows binary compatible strongly depends on MATLAB revision.Assign which python executable MATLAB should use with command  
```MATLAB
pyenv('Version','_full_path_to_python.exe_')
```

Other Python installs of interest
- {m}`sdtpy.install('scipy')` , shows the system commands needed for this installation. You must execute those manually. 
- {m}`sdtpy.install('scipy')` lists command to install SciPy. Verify `import scipy; scipy.signal.butter(1,.1,"lowpass")`
- {m}`sdtpy.install('gmsh')` lists command to install GMSH. Verify using `import gmsh; gmsh.initialize(); gmsh.fltk.run(); gmsh.finalize()`
	
(collab_install.JupyterBook)=
## Build Jupyter book documentation

If using Python / Jupyter book to generate this documentation
- First setup Python for MATLAB following instructions in 
- Install Jupyter book for Python by executing the commands listed by the call below in a Windows CMD/unix terminal   
{m}`sdtpy.install('jup')`
- Update the matlab lexer used by `Pygments` to color code to add SDT flavor.  
Use command {m}`sdtm.jup('sdtlexer')`
- {m}`sdtpy.install(''jup'')`) provides download instructions 
- {m}`sdtm.jup('build{support}')` generates the documentation of a specific book

(collab_install_VSCode)=
## VSCode as Markdown editor

- VSCode (Visual Studio Code) is a modern open source code editor
- It integrates with {m}`sdtu.f.open` if you set it as the editor for Markdown
  {m}`sdtdef('SDTools.MdEdt-setpref',{'VScode','c:\Users\$USER$\AppData\Local\Programs\Microsoft VS Code\Code.exe'})`
- You can possibly define this as text editor instead of the SDTools standard NotePad++ sdtdef('SDTools.TxtEdt-setpref',{'notepad++','D:\APP\win64\npp\notepad++.exe'})
- `Ctrl+k  v` opens the preview editor side by side 
- The explorer allows nice navigation with multiple files if you define a workspace. For example `sdt.code-workspace` contains

```
{"folders": [
		{"path": "support/jup"},{"path": "test2/jup"}
	]}
```

- Suggested extensions 
  - LaTeX workshop (from James Yu) 
  - LTEX for spelling / grammar compatible with both LaTeX and Markdown.
    - Note that add to dictionary is in Quick fix:Add to dictionary 
    - In settings search for ltex
    - Change item=dictionary to value=userExternal file
    - Change ltex:Dictionary (opens the `settings.json`) where you should have 

```
    "ltex.enabled": true,
    "ltex.dictionary": {
          "en-US": [":path_to_be_edited/sdt.cur/tex/en-US.usr"]
            },
```

Note that an execution bug may require to manually edit the @user/.vscode/extensions/valentjn.vscode-ltex-13.1.0/dist/extensions.js file to contain `shell:true`. See https://github.com/valentjn/vscode-ltex/issues/886.
```
        const executableOptions = {
                encoding: 'utf-8',
                timeout: 15000,
                shell: true
              };
```


VSCode is now the preferred editor for MarkDown (used to open files when calling {m}`sdtu.f.open(''File.md'')`). 

 Preview side-by-side (Ctrl+K V) 
 Install jupyter notebook support


 python3 

!"C:\Program Files (x86)\Microsoft Visual Studio\Shared\Python39_64\python.exe" -m pip install jupyter-book

Note that VSCode uses https://katex.org/docs/supported  for default Math preview which is not yet compatible with SDT math macros. 


(collab-install-GitHub)=
## GitHub desktop as git interface

- [GitHub desktop](https://github.com/apps/desktop) is an intuitive interface for Git access. It is the choice used by SDTools to explain git to users that are not familiar with this versioning system.   

(xxxtesth3)=
### xxx test


(xxxtesth2)=
## xxxtest2