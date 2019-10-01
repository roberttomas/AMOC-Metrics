# AMOC-Metrics
This is a first repository for AMOC Metrics work I began during January/February 2019.  To date I've forcused on the RAPID
(26.5N) and MOVE (16N) regions.  I started using code available from Roberts (RapidMoc, 201?) on github written in Python.  
To use the code with LR and HR POP, I wrote code in NCL to modify POP history files so that they are compatible with the 
RapidMoc code.  After initial work with the RapidMoc code, I started looking at the MOVE region and wrote my own code using
NCL.  I started using model velocities and later, calculated geostrophic velocities using monthly averaged T, S and SSH, and
bottom up and top down methodologies.
