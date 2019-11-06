

# Adaption of Salome Goem module  into an independent library

to compile with OpenCASCADE (Open source CAD kernel) without the whole salome, just like `SMesh`, an extraction of Salome mesh module 

by Qingfeng Xia, 2019, UKAEA

License: LGPL as Salome

## Salome Geom module 9.3

### Introduction of Geom module

Salome 9.x has the Shaper module, more dev may go into that new module

### adapted features

+ GEOMAlgo

### tested occt version

 it should works with OCCT 7.x, proabably not community edition (OCE)



## process of adaption

### modification on `CMakeFiles.txt`
after comment out "GetInPlaceAPI.cxx" in CMakeFiles.txt

this module compiled, without link to GeomUtils. kernel module

### OCC version macro 
replacing OCC version macro header and code in `GEOMAlgo_Gluer.hxx`

```cpp
//#include <Basics_OCCTVersion.hxx>
#include <Standard_Version.hxx>

#if OCC_VERSION_HEX >= 0x070200

OCC_VERSION
```






## tool to automatically track salome development

```bash
#git generate patch after comparison

#CI test
```



### examples for CI