# yeast_virulance_project
Yeast virulance selection genomic data analysis with the Murphy Lab

2026.11.04
First script used a python enviroment to implement picard. 
Now attempting to use picard built into gatk version 4 instead.

To deal with storage in the future use I made a symbolic link between the yeast_virulance project directory and the scratch directory with:
"mv ./ output /sciclone/scr10/jbusbin/yeast_vir_proj_output/"
and 
"ln -s /sciclone/scr10/jbusbin/yeast_vir_proj_output ~/murphy_lab/yeast_virulance_project/output" 

The script could be made for more reproducable with using a few more bash variables for the directories strucutre. And more personalized documentation. 
 
