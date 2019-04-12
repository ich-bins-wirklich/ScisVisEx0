# Exercise 0 - Framework Introduction

## Folder structure

Your folder should look like this:

    exercise
        |- build        --> everything will be build in this folder
        |- 0-intro      --> source code of current exercise
        |- cgv          --> clone of cgv framework repository


### workflow for windows

- clone the [cgv framework](https://github.com/sgumhold/cgv.git)
- drag *build* directory onto *cgv/define_system_variables.bat*
- drag *0-intro/intro.pj* onto *cgv/bin/generate_makefiles.bat*
    - use our VS version (for VS2017 use option *g*)
- start solution from *build/vs141/scivis-intro/scivis-intro.sln*
- in Visual Studio build solution either as Debug DLL or Release DLL
    - (for VS2017 either change *Windows SDK-version* in project properties or install Windows SDK-version 10.0.17134.0)
