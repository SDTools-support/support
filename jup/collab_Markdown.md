```{include} ../header.md
```
(markdown)=
# Markdown examples (as used in documentation)

```{tableofcontents}
```

(markdown-SDT)=
## SDT 

- MATLAB code generated with ``{m}`function('hello','python')`:``

{m}`function('hello','python')`

- MATLAB command **S**tring generated with ``{s}`MATLAB` ``+ `_static\sdt.css` 

{s}`MATLAB`

- Code block
``` matlab
%% Step1 : Initialize model
model=struct('Node',[1 0 0 0  0   0   0;     2 0 0 0 0   0   0.15;
                     3 0 0 0  0.4 1.0 0.176; 4 0 0 0 0.4 0.9 0.176],...
             'Elt',[],'unit','SI','name','GARTEUR');
%% Step2 Fuselage
model.Elt=feutil('ObjectBeamLine 1 2',model);
model=feutil('Extrude 0  1.0 0.0 0.0',model,...
             [linspace(0,.55,5) linspace(.65,1.4,6) 1.5]);
```

- Custom directives and roles (hello world level...)

```{helloworlddirective}
```

{helloworldrole}`test`

```{role} bluecolor(raw)
:format: html latex
```
{bluecolor}`<font color="blue">BlueColor custom role</font>`

(collab-cross)=
## Cross references

https://myst-parser.readthedocs.io/en/latest/syntax/cross-referencing.html

````{list-table} SDT math macros
:widths: auto
:header-rows: 1

* - Insert command
  - Insert result
  - Reference command
  - Reference result
* - `{{fig % ("id.png","FigID","caption")}}`
  - {{fig % ("id.png","FigID","caption")}}
  - `{refnum}`FigID``  
  `{ref}`FigID``
  - {refnum}`FigID`  
  {ref}`FigID`

````

## Tables 

Using pipes and - 

 | a | b | c |
 | ---: | :--- | :---: |
 | $+-$ right | left | center |
 | $ -$ right | left | center |

Using jupyter :::

xxx need documentation


(markdown-VS)=
## VS code shortcuts 

 Preview side-by-side (Ctrl+K V) 
 Install jupyter notbook support


 python3 

!"C:\Program Files (x86)\Microsoft Visual Studio\Shared\Python39_64\python.exe" -m pip install jupyter-book
   

(markdown-Math)=
## Math

````{list-table} SDT math macros
:widths: auto
:header-rows: 1

* - Command
  - Mathjax rendering	
* - `\ma{A}`
  - $\ma{A}$
* - `\mam{a & b \\ c & d}`
  - $\mam{a & b \\ c & d}$
* - `\vem{a & b \\\\ c & d}`
  - $\vem{a & b\\c & d}$
* - `\pam{a & b \\\\ c & d}`
  - $\pam{a & b\\c & d}$	
* - `\ve{A}`
  - $\ve{A}$	
* - `\pa{A}`
  - $\pa{A}$		
* - `\pa{-\omega^2\ma{M}+j\omega\ma{C}+\ma{K}}\ve{x}=0`
  - $\pa{-\omega^2\ma{M}+j\omega\ma{C}+\ma{K}}\ve{x}=0$
* - `\diag{A}`
  - $\diag{A}$
* - `\norm{A}`
  - $\norm{A}$
* - `\su{A}`
  - $\su{A}$
* - `\du{A}`
  - $\du{A}$	
* - `\oft`
  - $\oft$
* - `\ofw`
  - $\ofw$
````

$$
  w_{t+1} = (1 + r_{t+1}) s(w_t) + y_{t+1}
$$
`$$(my_other_label)` for jupyter


Ref to equation {eq}`my_other_label`. Has failed;

Next 

