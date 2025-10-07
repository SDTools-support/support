```{include} ../header.md
```
# Basic navigation in code/documentation 

- Learn about SDT {m}`sdtweb start`
- Learn about usual variables and coding style {m}`sdtweb syntax`


## Navigate within SDT code and documentation 


- {m}`sdtweb('_taglist')` (see {m}`sdtweb('sdtweb#_taglist')`)
  - xxx document tags   -3  at end of line without trailing spaces gives level of tag. 
  - Add `t` favorites in main MATLAB window 
     ![](./_images/Fav_taglist.png)

- {m}`sdtweb('d_doe','nmap')` (identical to {m}`sdtweb d_doe nmap`) opens the MATLAB editor at the proper tag within the `d_doe.m` file. 
- {m}`sdtweb('_mtag','nmap')` searches tags within all `.m` files.
   - This requires that you have either `grep` function in your system path (Cygwin for example) or that you have a `grep.m` function
   - Search path is `sdtroot`, `sdtroot/*`, `sdtroot/../sdt.git/*/m/*`
- {m}`sdtweb('_textag','nmap')` searches tags within all `.tex` files.
   - Search path is `sdtroot/tex`, `sdtroot/../sdt.git/*/tex/*`


## Build your tag list

- https://www.sdtools.com/helpcur/syntax.html#syntCode
