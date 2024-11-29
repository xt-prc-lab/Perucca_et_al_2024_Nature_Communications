With these codes, the intensity profiles of Fibroblasts and ECM images are computed.
The orientations of Fibroblast and ECM toward the Fibrobast edge are computed.
This document describes the overall process. Each script include more details about 
how it works and how to use it.

* I/ System requirements *
The scripts have been used in Windows 10 and Linux 22 PCs. The Matlab version should be R2020b or higher.



* II/ Installation guide * 
MATLAB R2020b or higher, with the following toolboxes : 
	- Image Processing
	- Signal Processing
	- Statistics and Machine Learning
	- Parallel Computing

CT-FIRE (https://eliceirilab.org/software/ctfire/)


* III/ Demo and instructions to run on data *
1. Create a main folder (here "MATLAB4ECM_example_data") including one folder per position to analyze.
Each fold to analyze should include : 
	- CAF.tif image 8bits, 515*512 pixels
	- If applicable ECM images (FN, COL_IV) : ECM.tif images 8bits, 512*512 pixels
Before running ctFIRE the folder "MATLAB4ECM_example_data" should have the following structure : 
|MATLAB4ECM_example_data/
├── folder1
│   ├── CAF.tif
│   ├── FN.tif
│   └── ctFIRE
│      ├── CAF.tif
│      └── FN.tif
│  
└── folder2
    ├── CAF.tif
    ├── COL_IV.tif
    └── ctFIRE
        ├── CAF.tif
        └── COL_IV.tif
        
2. Open FIJI.
 Open each image and register the background pixels values in a cell-free region (doing a mean of several points if necessary).

3. Open CT-FIRE and load your single / batch images.

4. Click "update" and edit the background pixel value measure in FIJI in the "tresh_im2" window.

5. Run CT-FIRE. After running ctFIRE the folder "MATLAB4ECM_example_data" should have the following structure :
|MATLAB4ECM_example_data/
├── folder1
│   ├── CAF.tif
│   ├── ctFIRE
│   │   ├── CAF.tif
│   │   ├── ctFIREout
│   │   │   ├── ctFIREout_CAF.mat
│   │   │   ├── ctFIREout_FN.mat
│   │   │   ├── ctfParam_CAF.tif.csv
│   │   │   ├── ctfParam_FN.tif.csv
│   │   │   ├── CTRimg_CAF.tif
│   │   │   ├── CTRimg_FN.tif
│   │   │   ├── currentP_CTF.mat
│   │   │   ├── HistANG_ctFIRE_CAF.csv
│   │   │   ├── HistANG_ctFIRE_FN.csv
│   │   │   ├── HistLEN_ctFIRE_CAF.csv
│   │   │   ├── HistLEN_ctFIRE_FN.csv
│   │   │   ├── HistSTR_ctFIRE_CAF.csv
│   │   │   ├── HistSTR_ctFIRE_FN.csv
│   │   │   ├── HistWID_ctFIRE_CAF.csv
│   │   │   ├── HistWID_ctFIRE_FN.csv
│   │   │   ├── OL_ctFIRE_CAF.tif
│   │   │   └── OL_ctFIRE_FN.tif
│   │   └── FN.tif
│   └── FN.tif
│  
└── folder2
    ├── CAF.tif
    ├── COL_IV.tif
    └── ctFIRE
        ├── CAF.tif
        ├── COL_IV.tif
        └── ctFIREout
            ├── ctFIREout_CAF.mat
            ├── ctFIREout_COL_IV.mat
            ├── ctfParam_CAF.tif.csv
            ├── ctfParam_COL_IV.tif.csv
            ├── CTRimg_CAF.tif
            ├── CTRimg_COL_IV.tif
            ├── currentP_CTF.mat
            ├── HistANG_ctFIRE_CAF.csv
            ├── HistANG_ctFIRE_COL_IV.csv
            ├── HistLEN_ctFIRE_CAF.csv
            ├── HistLEN_ctFIRE_COL_IV.csv
            ├── HistSTR_ctFIRE_CAF.csv
            ├── HistSTR_ctFIRE_COL_IV.csv
            ├── HistWID_ctFIRE_CAF.csv
            ├── HistWID_ctFIRE_COL_IV.csv
            ├── OL_ctFIRE_CAF.tif
            └── OL_ctFIRE_COL_IV.tif

6. Run the MATLAB script "MAIN.mat" in order to  :
	- Draw the Fibroblast edgen export a matrix where the value of each pixel is its distance to the closest piece of edge
	- Compute intensity profiles of the given images, as function of the distance to the edge
	- Compute the width of the ECM barrier
	- Compute the width of the CAF barrier
	- Use the ctFIRE fibers to find their angle toward the edge
	- Compute the angle profile of the given images, as a function of the distance to the edge
	- Merge angle and intensity profiles
	- Convert the 'MEASURE.mat' to '.csv'
	- Merge the profiles across different positions

9. With 'Patch_export_profiles_07_02_22.mat', the intensity profiles of several positions are gathered


On the example data provided, the whole process takes around five minutes, because the angles maps are already 
smoothed (./Matlab_data/Vecfield.mat). The computation and smoothing of these maps, from the previously drawn Fibroblast edge, 
take from 2 to 5 hour per matrix, depending on their size (512*512, 1024*1024, 2048*2048). A GPU accelerated version of this filter 
implemented in Python is available on demand. 
