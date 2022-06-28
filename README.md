# CRTBP Atlas
MATLAB companion codes for globally searching for connecting orbits in the Circular Restricted Three Body Problem (CRTBP). This is companion code which includes Atlas and Chart classes which inherit from the base classes in the [IMP library](https://github.com/skepley/IMP). The atlases and mining algorithms are both adapted to dynamically perform Levi-Civita regularization during the analytic continuation of the manifolds. 

## Searching for connections
For simplicity the code can be run with only a few steps. 

* Clone the IMP library and the CRTBP Atlas repository. 
2. Save a copy of the */CRTBP Atlas/production runs/template file with a new name. Save it **without** any file extension (e.g. "shane\_template" is ok but not "shane\_template.txt").
3. In your new template, edit the following code block with your own information.

```
17 %   Author: NAME
18 %   email: EMAIL
```
* In your new template edit the following code block with the path to your local copy of the IMP library, the CRTBP Atlas repo, and the path where you would like the code to save your data. The atlases are often 10GB or more so I suggest you do not save to Dropbox or iCloud, etc. 

```
27 addpath(genpath('PATH/TO/IMP/LIBRARY')) % this should point to the base folder named "IMP"
28 addpath(genpath('PATH/TO/CRTBP ATLAS'))  % this should point to the base folder named "CRTBP Atlas"
29 % path to use for saving/loading manifold and connection data. MUST END IN BACKSLASH! DO NOT USE DROPBOX FOLDER!
30 savePath = 'PATH/TO/SAVE/DATA/';  
```

Save and close the new template

* Run the MATLAB function */CRTBP Atlas/production runs/new\_run\_script.m. The first argument is the name of the template created above. The second is the desired name for the script which will execute the connection search. The remaining choices are related to parameters for the CRTBP. Example code

```
templateFilename = 'shane_template'
saveAsFilename = 'connection_search_1'
mu = 0.5
energy = 3
time = [-5, 5]
source = 'newL4_equalmass_local_manifolds.mat'
target = 'P1'
new_run_script(templateFilename, saveAsFilename, mu, energy, time, source, target)

```
The output of the script is another MATLAB script named 'connection\_search_1.m' saved in */CRTBP Atlas/production runs. 

* Run the script just created. 
