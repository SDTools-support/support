```{include} ../header.md
```

```{contents}
```

(collab-git)=
# Directory tree and GIT access

(collab-tree)=
## Standard directory structure

Each user should place (either physically or through symbolic links)
 - {m}`sdt.daily` xxx . {m}`sdtcheck patchupgrade` 
 - {m}`sdtcheck patchmkl` xxx / sdtdef problem if exist
 -  clone by url https://github.com/SDTools-support/support

xxx clone_Support.png

 - For internal SDTools processes this is {m}`sdt.cur`
 - {m}`sdt.git` contains all the project directories. For example 
   - {m}`sdt.git\support`  should be a clone from https://github.com/SDTools-support/support
 - {m}`sdt.git/repo` should then contain the standard directories 
   -  {m}`sdt.git\repo\m` should contain `.m` source files. Standard file names are  
      -  {m}`d_*.m` correspond to demos and are never pcoded. 
      -  {m}`t_*.m` correspond to test files.
      -  {m}`*(year).m` {m}`dev24.m` correspond to a project and the associated year. 
-  {m}`sdt.git\repo\tex` should contain `.tex` source files.  
-  {m}`sdt.git\repo\tex\plots` can contain `.tex` source files but **do not add binary figures (png/pdf/jpg) in the directory**. 
-  {m}`sdt.git\repo\jup` should contain `.md` or `.rst` source files for documentation generation. 

(collab-GitHub)=
## GitHub desktop as git interface

- [GitHub desktop](https://github.com/apps/desktop) is an intuitive interface for Git access. It is the choice used by SDTools to explain git to users that are not familiar with this versioning system.   

