// This script allows to choose which physical model to setup interactively:
//
// 1) Create a geometry with Gmsh, or load an existing geometry (".geo" file)
//    in Gmsh with `File->Open'
// 2) Merge this file ("Interactive.pro") with `File->Merge'
// 3) You will be prompted to setup your materials, sources and boundary
//    conditions for each physical group, interactively
// 4) Press "Run" to solve the model
//
// Everytime "Run" is pressed a ".pro" file is created (with the same prefix as
// the geometry file) with all the choices made interactively, for later non-
// interactive use.

Include "pile_common.pro";