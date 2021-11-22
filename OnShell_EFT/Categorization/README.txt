This file has the newest categorization script. Uses root_numpy package to quickly remove branches and operate on and make new branches.

To Run this code call the funciton ./gammaH_Categorization.py -i <input tree> -o <Output_Folder> -b <branchlist>

where branchlist is the list of branches that you want to have saved. 
There is also a script for submitting to condor to run categorization on all trees in the alltrees.txt file

To add your own categorization you can perform the categorization during the initial loop over all entries in the input tree. 
Save any calculated values to an array (This represents a new branch you want to add)
Convert that array into a numpy array 
Then use the array2tree function to add the branch to the newest tree. 
