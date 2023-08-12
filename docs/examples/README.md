MultiFEBE - Examples
====================

This document lists all the included examples for users, and it also
describes the guidelines for developing and including new examples.

Examples are organized by the type of problem `PP`, type of analysis
`AA`, modeling `MM`, and number `NNN`, i.e. a generic example is encoded
as `PP-AA-MM-NNN`. The examples whose title starts with "\[TUTORIAL\]"
include a document with a description of the modeling approach, input
files, output files and a comparison of the obtained results against
other reference results (analytical or otherwise), which serves as
validation.

List of examples for users
==========================

-   `ME` Linear elastic mechanics

    -   `ST` Static

        -   `SR` Structures

            -   `001` [\[TUTORIAL\] Cantilever beam (FEM
                model)](ME-ST-SR-001/).

        -   Elastostatics (continuum)

            -   `001` [\[TUTORIAL\] Plain strain square (BEM
                model)](ME-ST-EL-001/).

            -   `002` [\[TUTORIAL\] Cube (BEM model)](ME-ST-EL-002/).

        -   `CO` Coupled

            -   `001` [\[TUTORIAL\] Floating pile response under head
                force or moment (BEM-FEM model)](ME-ST-CO-001/).

    -   `TH` Time Harmonic

        -   `SR` Structures

            -   `001` [\[TUTORIAL\] Cantilever beam (FEM
                model)](ME-TH-SR-001/).

        -   `EL` Elastodynamics (continuum)

            -   `001` [\[TUTORIAL\] Cube (BEM model)](ME-TH-EL-001/).

            -   `002` [\[TUTORIAL\] Cantilever wall (BEM
                model)](ME-TH-EL-002/).

            -   `003` [\[TUTORIAL\] Soria arch dam: fixed-base (BEM
                model).](ME-TH-EL-003/)

        -   `PO` Poroelastodynamics (continuum)

            -   `none`

        -   `AC` Acoustics (continuum)

            -   `none`

        -   `CO` Coupled

            -   `001` [\[TUTORIAL\] Impedances of inclined pile
                foundations (BEM-FEM model)](ME-TH-CO-001/).

            -   `002` [\[TUTORIAL\] Seismic response of a single
                inclined pile (BEM-FEM model)](ME-TH-CO-002/).

            -   `003` [\[TUTORIAL\] Impedances of a suction caisson /
                bucket foundation (BEM-FEM model)](ME-TH-CO-003/).

            -   `004` [\[TUTORIAL\] Seismic response of an Offshore Wind
                Turbine (BEM-FEM model)](ME-TH-CO-004/).

            -   `005` [\[TUTORIAL\] Soria arch dam: compliant base (BEM
                model)](ME-TH-CO-005/).

            -   `006` [\[TUTORIAL\] Soria arch dam: compliant base with
                112 m water depth (BEM model)](ME-TH-CO-006/).

Notes for developers
====================

Examples should be categorized within the following tree:

-   `ME` Linear elastic mechanics

    -   `ST` Static

        -   `SR` Structures. Examples using only structural and discrete
            finite elements.

        -   `EL` Elastostatics (continuum). Examples using either
            boundary elements or finite element, where the domain is
            discretized without structural simplifications, i.e. a
            continuum or solid model.

        -   `CO` Coupled. Examples where structures are coupled with the
            surrounding media.

    -   `TH` Time Harmonic

        -   `SR` Structures. Examples using only structural and discrete
            finite elements.

        -   `EL` Elastodynamics (continuum). Examples using either
            boundary elements or finite element, where the domain is
            discretized without structural simplifications, i.e. a
            continuum or solid model. Only elastic solids.

        -   `PO` Poroelastodynamics (continuum). The same as above, but
            when using poroelastic media.

        -   `AC` Acoustics (continuum). The same as above, but when
            using inviscid fluids.

        -   `CO` Coupled. Examples where structures are coupled with the
            surrounding media, and/or where different material models
            are interacting.

-   `LA` Laplace (to be developed)

Within the `examples` folder, the example must have a dedicated folder
with the name `PP-AA-MM-NNN` according to its categorization and
numbering. The structure within the example folder must be:

-   `PP-AA-MM-NNN.pdf` - It is a mandatory document describing example.
    It can be as simple as a readme type of document, where a brief
    description of the example is given, or it can be a fully documented
    example in a tutorial format.

-   `doc_src` - It is a mandatory folder containing the sources for
    creating `PP-AA-MM-NNN.pdf`. It can be a plain text file
    `PP-AA-MM-NNN.txt`, a Microsoft Word file `PP-AA-MM-NNN.doc` or
    `PP-AA-MM-NNN.docx`, a Open Document Text file `PP-AA-MM-NNN.odt`,
    or a LaTeX file `PP-AA-MM-NNN.tex`. The preferred document source is
    LaTeX.

-   `case_files` - It is a mandatory folder containing all input files,
    including MultiFEBE, mesh and other files required for running the
    example.
