
----------------------------------------------------------------------------------
Hawkins, BL. 2019. Exploring the role of temperature and possible alternative 
	stable states in brook trout (Salvelinus fontinalis) population structure. 
	M.S. Thesis. Univeristy of California San Diego.
----------------------------------------------------------------------------------

To recreate my analysis, here’s what to do:

1) Download R and RStudio.
2) Create a new R project called FishLengths in RStudio.
3) Create three directories in the FishLengths directory.
	Code
	Data
	Results
4) Move hawkins_thesis_analysis.R and hawkins_thesis_models.R to the Code directory.

The hawkins_thesis_analysis.R script uses the lake and fish length datasets from 
Rob Grasso and Roland Knapp, extracts population size structure features, and then 
explores their relationship with environmental factors. The hawkins_thesis_models.R 
script includes the normal mixture models used in step 3 of the analysis script.

5) Get the fish and environmental datasets (courtesy of Rob Grasso and Roland Knapp, 
   contact me or Shurin Lab for these files) and move them to the Data directory.
6) Check that your directories and files are structured the following way.

	FishLengths
	   Code
	      hawkins_thesis_analysis.R
	      hawkins_thesis_models.R
	   Data
	      RobGrassoYosemite
	         meta.data.csv
	         length.data.csv
	      RolandKnappYosemite directory
	         Lake ID and elevation.csv
	         LakeData_ToShurin_10Apr2019.csv
	   Results

7) Install libraries in the setup section of hawkins_thesis_analysis.R.
8) Source hawkins_thesis_analysis.R.

Once you source the script, you will end up with a new directory in Results, named with
the date and time you sourced the script. It will include seven files:

Tables from thesis
	1) Table1.csv
	2) Table2.csv
	3) Table4.csv
Figures from thesis
	4) Figure2.png
	5) Figure3.png
Novel datasets describing population size structure
	6) population.csv
	7) sizeclass.csv

Note that this script uses deterministic and stochastic processes in four steps to extract 
features of a population's size structure. The first two steps, calculating average fish size
and finding peaks, are deterministic. Therefore, the results regarding mean fish size and 
the number of peaks for each population will be consistent every time you run the script. 
However, fitting Bayesian normal mixture models in step 3 uses MCMC simulations, which are 
inherently stochastic. Therefore, the number of size classes and evenness of size classes 
will be slightly different each time you run the script, as will the novel datasets for the 
population and size class features.

