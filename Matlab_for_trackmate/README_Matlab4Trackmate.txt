This repository contains the code for the paper Micro Immune Response On-chip (MIRO) models the tumour-stroma interface for immunotherapy testing, by Alice Perucca, Andrea Gómez Llonin, Oriol Mañé Benach, Clement Hallopeau, Elisa I. Rivas, Jenniffer Linares, Marta Garrido, Anna Sallent-Aragay, Tom Golde, Julien Colombelli, Eleni Dalaka, Judith Linacero, Marina Cazorla, Teresa Galan, Jordi Pastor Viel, Xavier Badenas, Alba Recort-Bascuas, Laura Comerma, Patricia Fernandez-Nogueira, Ana Rovira, Pere Roca-Cusachs, Joan Albanell, Xavier Trepat, Alexandre Calon and Anna Labernadie

Codes details

1. System requirements
The scripts have been used in Windows 10 and Linux 22 PCs. The Matlab version used was: R2020b (tested on higher versions)

2. Installation guide
The following toolboxes are necessary to run this code:
	- Image Processing
	- Signal Processing
	- Statistics and Machine Learning
	- Parallel Computing
It will be necessary to add to MATLAB path the FIJI folder 'scripts'.

3. Demo and instructions to run on data
The published Matlab codes are aimed to compute immune cell track characteristics previously obtained with FIJI plugin Trackmate.
The code to run is "MAIN.m" and here we provide a description of the structure of the input folder and the type, characteristics and specific naming of the files.
|Matlab4Trackmate_example_data/
├── Frames : Frames folder / For each position, 1 binary mask of the immune cells channel for each timepoint exported directly from FIJI (FIJI -> Save as -> Image Sequence); tif image 8bits 1024x1024 pixels
│   ├── n2_Pos7_C20000.tif
│   ├── n2_Pos7_C20001.tif
|   ├── ...
│   ├── n2_Pos7_C20089.tif
│   ├── n3_Pos5_C20000.tif
|   ├── ...
│   └── n3_Pos5_C20089.tif
├── n2_Pos7_C2.tif : Image of Fibroblasts; 1 tif file per position corresponding to the first frame of the fibroblasts channel; RGB image, 1024x1024 pixels
├── n2_Pos7_C2_Tracks.xml : Track file; xml file exported from the FIJI Trackmate plugin (FIJI -> Trackmate -> Export tracks to XML file); named according to the corresponding fibroblast image
├── n3_Pos5_C2.tif
└── n3_Pos5_C2_Tracks.xml

The output of "MAIN.m" includes:
- Radius folder: the distance matrix for each fied of view. For each pixel, the value is the distance to the closest piece of edge, in micrometers.
- Morphological parameters: cell tracks, with the morphological parameters, after selection, for each field of view. 
- Matlab_data: visual validation of the tracks
- CSV file: dataset containing information about the immune cell tracking, such as X and Y coordinates, track ID (ID), immune cell area (Area), minor and major axis (Majax, Minax), circularity, distance to the edge (dist), max area per track (MAX_Area), position as respect to the edge ('IN', 'EDGE', 'OUT') (bin), displacement angle towards the edge (Theta), Instant and Mean displacement (Di, Dm), Instant and Mean Velocity (Vi, Vm), Roundness and  Mean roundness (Ri, Rm), Aspect Ratio and Mean aspect ratio (AFi, AFM)

Additional useful information:
- in the section 'Draw the edge' the user will be asked to draw a poligon along the CAF edge and then to validate the drawing.
- in the section 'Measure and save morphological parameters' the user will be asked to manually validate each track (between Trackmate and Matlab???)
- in the section 'create dataset' the user will be asked to input the time between 2 frames (seconds) for each position.

Time required: few minutes.

4. Instructions for use
The codes can be downloaded from the online repository with the demo data. Explanations on the tasks performed by each line are contained in the codes.