# Exercise 2 - Particle

## folder structure

All exercises are independent of each other.

    sci-vis-exercises
        |- build        --> solutions, libs, binaries will be build in this folder
        |- 0-intro
        |- 1-stereo
        |- 2-particle   --> source code of recent exercise
        |- cgv          --> clone of cgv framework repository
        |- data         --> example data like obj-files


## work-flow for windows

### optional: set up cgv framework

You can skip this step, if you already used the framework for the previous task.

- clone the [cgv framework](https://github.com/sgumhold/cgv.git)
- drag *build* directory onto *cgv/define_system_variables.bat*

When you initially use the *cgv/bin/generate_makefiles.bat* in the following steps, choose your VS version (for VS2017 use option *g*).


#### create VS solution and compile

- pull recent master branch from cgv framework!
- drag *particle.pj*-file onto *cgv/bin/generate_makefiles.bat*
- start solution from *build/* directory (e.g. *build/vs141/scivis-particle/scivis-particle.sln*)
- in Visual Studio: build solution either as *Debug DLL* or *Release DLL*
    - for loading meshes the Release DLL is the preferred option since it will take very long to load it in Debug mode


### known issues

- for VS2017 either change *Windows SDK-version* in project properties or install Windows SDK-version 10.0.17134.0