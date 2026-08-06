# EphemerisPc

The EphemerisPc tool uses ephemeris tables to calculate the statistically expected number of collisions (Nc), along with bounding collision probabilities
(PcMax and PcMin). It nominally evaluates the entire ephemeris overlap span and generates summary plots for the full assessment span as well as for each
encounter segment.

For more information on the methodology, read "Ephemeris-Based Satellite Collision Rates and Probabilities," by Doyle Hall in '../../References'.
The paper also covers example scenarios, including low relative velocity LEO and GEO conjunctions, and a high relative velocity LEO conjunction.
For a quicker overview of the tool's outputs and usage, see 'EphemerisPc_Example.pptx' in the 'Example' folder.

The primary function is EphemerisPc.m. The header contains all of the information needed to run the tool, including a detailed description of inputs and outputs.
The function ingests either the CCSDS OEM ephemerides or BFMC-generated ephemeris/minima/VCM files, and calculates the Nc, PcMax, and PcMin.
While Pc represents the probability of at least one collision, Nc can sum expected collisions from multiple encounters. 
For a single independent encounter Nc≈Pc, but they can diverge when encounter segments blend together, which can occur with low relative velocity conjunctions.
The tool does not calculate cumulative Pc, but provides cumulative Pc bounds. 
Where the true cumulative Pc lies within those bounds depends on how independent the encounter segments are.
PcMin assumes all encounters are fully dependent, while PcMax assumes they are fully independent. 
Reference the paper for more details.

The function currently supports the following input modes:
1. OEM Mode: Provide the path to a folder containing standard CCSDS OEM ephemerides. Public users will run EphemerisPc in this mode.
2. BFMC Mode: Provide the path to an output folder from a BFMC run. Public users will not have the inputs to run EphemerisPc in this mode. It is intended for CARA analysis team internal use.

OEM ephemerides can be generated from raw ephemeris data using 'src/CCSDSWriter.m'.

The following folders are contained within the EphemerisPc directory:

1. Analysis\_Team\_Driver\_Scripts: Contains driver scripts for running EphemerisPc in OEM or BFMC mode and for a single ephemeris pair or a series of ephemeris
pairs. This folder is intended for the internal use by the CARA analysis team.
2. Documentation: Contains an EphemerisPc user guide
2. Example: Contains a saved example EphemerisPc run, a script for running the example scenario, and a slide package outlining it.
3. src: Provides supporting functions for EphemerisPc.m.
4. TestEphOffset: Test code intended for internal CARA analysis team use.
5. UnitTest: Contains a unit test for EphemerisPc

