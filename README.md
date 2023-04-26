# trial-invitation-algorithm

An algorithm to decide which patients to invite to join a trial.

The code is demonstrated using a model where patients are organised into age/sex groups by GP practice

- **prog**: python code for the invitation algorithm
- **inputgen**: R code and associated demonstration data. The code generates input files for the algorithm using the demonstration data

## Overview 

To run the algorithm the following steps need to be done

1. Provide a list of GPs in the same format as the demonstration (see **inputgen/updatefiles**)
1. If a second round invitation list (wave 2+), add update data in the same format as the demonstration (see **inputgen/updatefiles**)
1. The R script in `inputgen/ ` will create input files for the algorithm 
	* `Rscript trial-invitation-geninput.r /path/to/inputfile ## round 1`
	* `Rscript trial-invitation-geninput.r /path/to/inputfile /path/to/bookingfile ##round2/3`
1. There is an example site at Kings College London provided. eg. To generate input files
	* `Rscript trial-invitation-geninput.r updatefiles/site_distance_KingsCollegeLondon.csv # round 1`
	* `Rscript trial-invitation-geninput.r updatefiles/site_distance_KingsCollegeLondon.csv updatefiles/bookinvites_KingsCollegeLondon.csv ##round 2`
1. The input files will be saved in genfiles with date and unit name in the filename, and also copied to `../prog/input/` for the algorithm
1. To run the algorithm, call `python3 prog/trial-invitation-frontend.py` and edit slots, proportion of remaining capacity to aim to book, uptake assumption
1. Output that is inflated for national opt outs is in `prog/output/nhsd-invite-inflated.csv`
 
# License

GNU GPL v3

# Last updated
26th April 2023, Adam Brentnall

