# ebi_project
Files to merge indelible and pyvolve functions to provided simulated data with different omegas and indels based on custom matrices.

#some files on here for personal reference:
- project.py, first file made, indelible and pyvolve simulations
- site.py, second file made, indelible, pyvolve, sites extracted from indelible dna_RATES
- sites_added.py, third file made, indelible, pyvolve, and sites then added to pyvolve simulations

- merge_test.py = merges the pyvolve simualtions and applies the indel patterns from indelible dna_TRUE.

- indel_pyvolve_simulator.py = the final rendition so far that will be editied and cleaned, does all the same as the one above but ACTUALLY merges based on indel patterns, and applies the indels in a vectorised way. INPUT = tree file path, ,output_file path. output = simualtions from indelible, pyvolve, merged sequences, sequence with indels and sequence with indels removed. 

- need to sort out files made and removed/ make sure they parse to eachother
- need to change hard coded paths so it can run without me
- sync git with both terminals ebi/me
