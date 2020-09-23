

# Adaption of Salome Geometry module  into an independent library

to compile with [OpenCASCADE](https://www.opencascade.com/), an industrial quality open source CAD kernel, without the whole [Salome](https://www.salome-platform.org/) , an open source full-featured CAE platform with GUI and python scripting), just like [`SMesh`](https://github.com/LaughlinResearch/SMESH), an extraction of Salome mesh module.

adapted by Qingfeng Xia, 2019, [UKAEA]()

License: LGPL v2.1 as Salome, OpenCASCADE

## Salome Geometry module

### Introduction of Geometry module

Salome 9.x has the `Shaper` module, more dev may go into that new module

### adapted features from Salome v9.3

+ GEOMAlgo:  BoundingSphere, BndBoxTree, Gluer, classifier, FinderShapeOn

+ ShHealOper: fillHoles

  

  what is the difference between OCCT's shape healing algo in `TKShHealing` toolket

Todo:

+ OCC2VTK:  may not quite necessary, as OCCT has such features
+ ShapeRecognition: depends on OpenCV, but it does not output TopoDS_Shape

### Tested occt version

 It should work with OCCT 7.x, but probably not with the community edition (OCE)



## Process of adaption

### modification on `CMakeFiles.txt`
After commenting out "GetInPlaceAPI.cxx" in CMakeFiles.txt, GEOMAlgo compiled, without link to `GeomUtils` which is a Kernel module of Salome.

### OCC version macro 
replacing OCC version macro header and code in `GEOMAlgo_Gluer.hxx` with OCCT official header and macro `<Standard_Version.hxx>`

```cpp
//#include <Basics_OCCTVersion.hxx>
#include <Standard_Version.hxx>

#if OCC_VERSION_HEX >= 0x070200

// use OCC_VERSION_HEX instead of OCC_VERSION
```

### Add some extra cmake modules

several cmake modules such as "FindOpenCascade.cmake" added into the `cMake` folder in this repo

`SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cMake)`


## tool to automatically track salome development

It is a feature in the todo list: 

+ git generate patch after comparison

+ CI test


