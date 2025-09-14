# Macroelement Doubletwist PFC7

This is the repository for the PFC7 implementation of the double-twisted hexagonal mesh wire model, implemented using remote contact bonds between discrete particles, within the Discrete Element framework.

The example cases folder include the code to setup and run the model in PFC3Dv7.0. They employ a pre-compiled version of the contact model, compiled as 'Release' dll.

If you want to compile it yoursef, make modifications and run the debugger, the source code is included in the respective folder. You will still need the PFC7 libraries, Qt and some version of Visual Studio, as described in the manual.
https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/optional_features/plugins/plugins.html#plugins-pfc

I had previously made a video showing how to use the VS debugger with PFC, which some people found useful: https://youtu.be/TZUDU_Ggcxg

# Notes & additional context
The code was compiled successfully using Visual Studio Community 2019 and 2022. The last modifications to the code occurred in May 2025. At the time, PFC9 was still under development, with alpha and beta releases, and not all features and files were available. For this reason, it was never ported to PFC9. In order to use the Von Mises constitutive model in a coupled FDM-DEM model, the FLAC9 model was ported to FLAC7 https://github.com/marcoprevitali/vonMises_FLAC7.
