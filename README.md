# Landsat-derived Arctic-optimised Sea Surface Temperature
**If you use this code and method please cite the following reference:**

Ing, R., Nienow, P., Slater, D., Medina-Lopez, E. (in review). Investigating the use of Landsat-derived Sea Surface Temperatures as a proxy for changes in ocean forcing of Greenland’s marine-terminating outlet glaciers. Journal of Glaciology.

The method used to derive and optimise the Landsat SST algorithm for Greenland coastal waters is explained in the above paper. The code used to derive the algorithm is available at: https://doi.org/10.5281/zenodo.20610025. The algorithm was trained on radiative transfer simulations for atmospheric profiles over Greenlandic coastal waters, however it has been shown through validation to work well for Arctic coastal waters in general. 

This repository provides example code for the python version utilising the Google Earth Engine Python API, with the main SST algorithm pipeline in "landsat_sst". **"Arctic_Landsat_SST_Demo.ipynb"** provides and example juypter notebook of how to view Landsat scenes and calculate SST.

A more accessible, user-friendly, Google Earth Engine app is detailed below.

## Google Earth Engine App:
A web app has been designed to allow users to view and create timeseries of Landsat-derived SST optimised for the Arctic. This web app is available here: (https://ryan-ing.projects.earthengine.app/view/arctic-landsat-sst-viewer)

### How to use the app:

**Step 1: Define geometry**
Draw a rectangle over your chosen study area. A smaller geometry is preferred, as it avoids scenes that only partially cover your area of interest from being included. The clear geometry button allows the user to remove all drawn geometries.

**Step 2: Define cloud cover**
Set max cloud cover (%) of scenes. A smaller number (5 %) is preferred.

**Step 3: Choose Landsat scene**
The app will now show all available scenes (1982 to present) with the applied geometry and cloud filters. Select an image to display it on the scene. 

**Step 4: Export to drive**
Unfortunatelty Google Earth Engine does not allow users to export to their google drive in web apps. To use this feature, see the section below. 

**Step 5: Visualisation**
Set the visualisation SST range for the scene displayed on the screen.

**Step 6: Inspector**
Click anywhere on the screen to view the SST value of the choosen pixel

![image](https://github.com/ry-ing/arctic-landsat-sst-GEE/blob/main/img/gee_screenshot1.png)

**Step 7 & 8. Time-series Analysis**
The step allows the user to generate a timeseries of SST for a chosen year range. A figure of monthly aggregates is also generated. This data can be directly downloaded as a CSV file. 

![image](https://github.com/ry-ing/arctic-landsat-sst-GEE/blob/main/img/gee_screenshot2.png)

### How to use the app (with google drive export feature):
To export scenes of Landsat-derived SST, the app needs to be run in the Google Earth Engine console. To do this users will need to sign up for a Google Earth Engine account and a Google cloud project. This can be done here: https://earthengine.google.com/signup/. 

After signing up for a GEE account, you should be able to access the console-version of the web app here: https://code.earthengine.google.com/f96545829b6e11b929a07eb3cf13901b. After clicking on the link, press "**Run**", and the app will load. The "Export Scene to Drive" button will now work and export scenes to your Google Drive. 

![image](https://github.com/ry-ing/arctic-landsat-sst-GEE/blob/main/img/gee_screenshot3.png)