

# Adaption of Salome Geometry module  into an independent library

to compile with [OpenCASCADE](https://www.opencascade.com/), an industrial quality open source CAD kernel, without the whole [Salome](https://www.salome-platform.org/) , an open source full-featured CAE platform with GUI and python scripting), just like [`SMesh`](https://github.com/LaughlinResearch/SMESH), an extraction of Salome mesh module.

adapted by Qingfeng Xia, 2019, [UKAEA]()

License: LGPL as Salome, OpenCASCADE

## Salome Geometry module

### Introduction of Geometry module

Salome 9.x has the `Shaper` module, more dev may go into that new module

### adapted features of Salome v9.3

+ GEOMAlgo

Todo:

+ OCC2VTK:  may not quite necessary, as OCCT has such features
+ ShHealOper: what is the difference between OCCT's shape healing algo
+ ShapeRecognition:

### Tested occt version

 It should work with OCCT 7.x, but probably not with the community edition (OCE)



## Process of adaption

### modification on `CMakeFiles.txt`
After commenting out "GetInPlaceAPI.cxx" in CMakeFiles.txt, GEOMAlgo compiled, without link to GeomUtils, Kernel module of Salome.

### OCC version macro 
replacing OCC version macro header and code in `GEOMAlgo_Gluer.hxx`

```cpp
//#include <Basics_OCCTVersion.hxx>
#include <Standard_Version.hxx>

#if OCC_VERSION_HEX >= 0x070200

OCC_VERSION
```

cmake support




## tool to automatically track salome development

```bash
#git generate patch after comparison

#CI test
```



### examples for CI