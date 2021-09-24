# glider-sampling-simulation
Simulated underwater glider sampling for Antarctic krill

To run the simulation, download all files in the 'R' and 'Data' directories and the rmarkdown files ('.rmd') from the 'Main' directory into the same working directory. Render the rmd scripts using 'knitr'. The results for one random simulation of 9 replicate samples are shown in the pdfs in the 'Main' directory for the corresponding '.rmd' files. The resolution of the pdfs is better when they are downloaded than when they are viewed directly from GitHub. 

This simulation produces about 8 GB of output files and takes several minutes. More replicates will produce better statistical results but require more time and produce more output. 100 replicates can require several days on a multicore computer and produce about 40 GB of output files.

The data to be sampled in the simulation are from AMLR research ship surveys in two sampling strata around the Antarctic Peninsula, the Southern Area ('SA', also called Bransfield Strait) and the Western Area ('WA', also called Cape Shirreff). 'Glider_simulation_SA.rmd' and 'Glider_simulation_WA.rmd' should be rendered first to produce the glider sampling files from each strata. Then 'Figures1-7_SA.rmd' and 'Figures1-7_WA.rmd' can be rendered to produce the Figures for each of the two strata. 'Figure8.rmd' combines glider samples from both strata.
