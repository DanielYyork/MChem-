The data collation script is designed to be used with orca output files.
To run properly it requires a file structure as shown below:

-Master folder
    - Biocrude chitin model
        - Molecule #1 (this is the folder generated during the orca calculation)
        - Molecule #2
        - Molecule #3
        - etc.
    - Biocrude seaweed model
        - Molecule #1
        - Molecule #2
        - Molecule #3
        - etc.
    - etc.
        
This code will create a .csv file containing a range of data acquired from the orca output file.

USAGE:
- This code will create 1 csv at a file for each folder within the master folder,
    the csv file will contain data for each molecules (i.e. molecule 1, 2, 3, etc.).
- It can currently only handle 1 folder at once,
    i.e. 2 runs of code will be required to generate the csv files for seaweed and chitin biocrude data.
- Before using, check the path is correct (i.e. to the Biocrude chitin model folder),
    and that the .csv file produced has the name you desire.