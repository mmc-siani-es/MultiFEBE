# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com), and this
project adheres to [Semantic Versioning](https://semver.org).


## [2.0.2] - 2025-09-21
### Added
- Smoothing of shell FE stress resultants by simple averaging at nodes.
- Choose degshell FE (standard or MITC) formulation via input data file, although MITC remains as the default.
- Change in the documentation folder structure. Now instead of "tutorials", all cases are "examples", some serve as validation cases, other as tutorial cases and others as simple examples of application. New examples are included.

### Changed
- MITC3 formulation is updated to MITC3i new scheme.
- Bug in shape functions of tri4 (triangular element with centroid node with cubic bubble function).
- Manual has been updated.
- Include BE body loads including point loads in time harmonic analysis of fluid & elastic solid. 
- Update calculation of internal points, including how results are exported. 
- Update stress tensor calculation along boundaries.

## [2.0.1] - 2023-02-24
### Added
- Handling of BE bodyload elements when coincident to symmetry planes

### Changed
- Path handling is now completely general, user can use relative paths or absolute paths. It is compiler/platform independent.

### Removed
- Internally force straight strbeam elements. Now is responsability of the user.
