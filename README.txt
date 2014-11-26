NCR_Code
========

This package contains the Matlab code implementing the NCR algorithm described in ECCV 2014 paper "Consistent Re-identification in a Camera Network".

================ CONTENTS OF THE FILE =================
The package contains the following 3 '.m' files.
1. NCR.m - The code for the case when all the persons are present in all the cameras. The input to the code is the camera pairwise similarity socres between all the persons. The format of the pairwise similarity score matrix is described as comments inside the 3rd cell titled - "Load pairwise similarity scores and generate optimization cost function and constraints". Other configurable parameters are described inside the 2nd cell titled - "Parameter initialization and addition of paths". The output are the CMC plots with and without using the NCR algorithm.
2. NCR_Generalized.m - The code for the generalized NCr when every person may not be present in every camera. The input for this case also, is the camera pairwise similarity socres. The format of the pairwise similarity score matrix is described as comments inside the 3rd cell titled - "Load pairwise similarity scores and generate optimization cost function and constraints". Other configurable parameters are described inside the 2nd cell titled - "Parameter initialization and addition of paths". The output are the individual and average (across different trials) accuracy value (as defined in the paper) for different 'k' values for each camera pair.
3. set_label_style.m - This is a function to format the fig files.

The package also contains two folders as follows.
1. Pairwise_Similarity - This is the placeholder for the mat files containing the pairwise similarity scores. As examples pairwise similarity scores for the test peoples for both WARD and RAiD dataset are provided. Note that these pairwise similarity scores are generated using the FT method as described in the paper. Using NCR.m on the WARD data will generate Fig 3 (Only NCR on FT and FT) while Using NCR_Generalized.m on 'Pairwise_sim_Less_People_RAiD_6a' and 'Pairwise_sim_Less_People_RAiD_6b' will generate data for Figures 6(a) and 6(b) respectively.
2. Results - This is a folder to store intermediate and final results so that the code can reuse data for repeated running. We have provided the code for saving these data for 'NCR.m' only. Please incorporate the same if you need to save them for 'NCR_Generalized.m' too. The placeholder for this purpose (GeneralizedNCR) is already provided.



================ PREREQS (CPLEX) =================
To solve the binary integer program we use "IBM ILOG CPLEX Optimization Studio". IBM ILOG CPLEX Optimization Studio is an optimization software package by IBM. This proprietary software is frrely available to the academicians through 'IBM Academic Initiative'. Please visit the following website to get more details about it - http://www-304.ibm.com/ibm/university/academic/pub/page/membership .
After installing CPLEX please connect CPLEX with Matlab by following the procedure.
a) Open MATLAB and type 'pathtool' (without quotes).
b) Add the CPLEX matlab folder which generally is 'C:\Program Files\IBM\ILOG\CPLEX_StudioXXXX\cplex\matlab\x64_win64' in Windows (XXXX denotes the version number).
c) (optional) Repeat the step (b) for the following folder to use example file of the help window. 'C:\Program Files\IBM\ILOG\CPLEX_StudioXXXX\cplex\examples\src\matlab' (XXXX denotes the version number).


================ PLEASE CITE =================
@inproceedings{Das2014,
	title     = {Consistent Re-identification in a Camera Network},
	author    = {Abir Das and Anirban Chakraborty and Amit K. Roy-Chowdhury},
	booktitle = {European Conference on Computer Vision},
	year      = {2014},
	publisher = {Springer},
	location  = {Zurich},
	series    = {Lecture Notes in Computer Science},
	volume    = {8690},
	pages     = {330--345}
}

For any questions, feel free to email at abir.das@email.ucr.edu
