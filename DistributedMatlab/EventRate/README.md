##
EventRate is a semi-empirical risk analysis tool that estimates the risk for a prospective mission by 
analyzing historical conjunction data from existing satellites in similar orbits. The primary risk metrics
estimated by EventRate are:

- Rate of conjunctions in which the calculated Pc exceeds some user-defined threshold (High-Risk Events or HREs),
- Cumulative Pc over the total duration of a mission.

EventRate can also produce other outputs that may be of use to mission designers, like altitude-latitude
plots to spatially visualize where conjunctions and risk are concentrated. A more detailed overview of 
EventRate outputs can be found in doc/EventRate_Description.


##
DIRECTORIES:

- data/ - A directory to store OCMDB/PcTable files used to run EventRate. A PcTable file containing conjunction data from
			a curated list of surrogate missions will be provided in this directory. If the data in this distribution does
			not meet their mission needs, users may contact CARA via an issue request on the GitHub site.
	
- doc/ - A directory containing documentation files:
	- EventRate_Description.pdf and EventRate_Description.pptx - A detailed overview of EventRate, its functionality, and its outputs
		> For more information on the theory behind the EventRate tool, see the document "Hall_2019_EventRateTheory_AAS_19-631.pdf" in the ../../References directory
	- EventRate Curated Satellite List.xlsx - A list of satellites included in the public PcTable file and their orbital parameters. This list can be used to choose suitable surrogates for mission analysis.
	- EventRate Orbit Regime Definitions.xlsx - A spreadsheet defining the various orbital regimes referenced in the example mission parameter files.
	- EventRate Truncated OCMDB Field Definitions.xlsx - A detailed overview of the OCMDB table format
	
	
- output/ - A directory where EventRate will store plots and log files. EventRate runs are each assigned a unique output folder named based on the configuration and timestamp at the start of the run	
- params/ - A directory to store EventRate parameter files
- src/ - The location of all EventRate source code and functions
	
##
USAGE:

To run EventRate, the user must supply at least an OCMDB extract file or a PcTable file derived from an OCMDB
extract and at least one satellite ID included as a primary in that extract file to serve as the mission surrogate. 
Mission information should be put into an EventRate parameter structure, examples of which may be found in the params/ directory. 
The parameter structure is also used to control the functionality of the EventRate tool. Unset parameters will be set to default values (See src/EventRate_default_params.m for more details.)

Once a parameter structure script is written, EventRate can be run from the MATLAB command line as follows:

- output = EventRate('test_mission')

where the parameter file being used in this example is params/test_mission.m. 
Refer to doc/EventRate_Description.pdf for further information regarding the output structure.