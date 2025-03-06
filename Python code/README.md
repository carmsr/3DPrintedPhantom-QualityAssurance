# 3D-printed phantom Quality Assurance: Python Original Code
@Author: __Carmen Salvador Ribés__

## First steps
<p align="justify">The code is developed in Python 3.7. To be able to use the code, you have to download the folder <em>3DPrintedPhantom-QualityAssurance/Python code</em> and install the packages specified in the file <em>requirements.txt</em> manually or automatically running the following command on pip:</p>

`pip install -r requirements.txt`

<p align="justify"><em>3DPrintedPhantom-QualityAssurance/Python code</em> folder contains two subfolders. Code scripts are in <em>Code</em> subfolder, and additional required files are in <em>Required files</em> subfolder.</p>


## Overview
<p align="justify">The 3D-printed phantom Quality Assurance code enables the analysis of different quality parameters of PET, CT and MR images (quantification, resolution, co-registration, distortion and radiomics) with 3D-printed phantom acquisitions.</p>
<p align="justify">The code is made up of 7 scripts that are shown in Figure 1. To start the analysis, the script <em>3DPrintedPhantom_QA.py</em> should be run, which allows to choose the type of quality analysis to be performed (Figure 2). <em>Auxiliar_functions.py</em> script contains all necessary functions used in the different analyses. The rest of the scripts perform the analysis corresponding to its name. You can also launch each analysis individually by running the corresponding script.</p>

<div align="center"><img src="https://github.com/user-attachments/assets/eac58c31-2c3e-42e0-9891-dab36b568134" width="25%"></div>
<div align="center"><strong>Figure 1. QA 3D-printed phantom code scripts.</strong></div>
<p>&nbsp;</p>

<div align="center"><img src="https://github.com/user-attachments/assets/c4838fbd-1f25-4534-9136-10640157f652" width="40%"></div>
<div align="center"><strong>Figure 2. Initial window that appears when running the script 3DPrintedPhantom_QA.py.</strong></div>
<p>&nbsp;</p>

<p align="justify"><strong>General prerequisites</strong>: It is necessary to have the PET, CT and MR images to be analysed in DICOM format.     
A specific folder organisation is required for the code to work (Figure 3). The working folder is created (in this case the <em>Analysis folder</em>), where we create a subfolder called Images in which we will place the folders with the different dates of the studies. Within these folders, we will place the DICOM series corresponding to an image in a folder with the name of the modality (CT, PET and MR). As an exception, in the MR folder there shall be subfolders corresponding to the possible existing sequences (as an example, two sequences T1 and T2), which shall contain the corresponding DICOM series. It is not necessary that the folders have these two names, they can have any name, the code parses the sequence type from the DICOM header. The <em>Output</em> folder will be automatically generated in the analysis process within the work folder. In it, different output files will be generated within the corresponding date folder.</p>

<div align="center"><img src="https://github.com/user-attachments/assets/58c8b942-7a7f-4e77-ace2-a5cab78f6ffb" width="70%"></div>
<div align="center"><strong>Figure 3. Folder organisation for QA 3D-printed phantom analysis.</strong></div>
<p>&nbsp;</p>

## QUANTIFICATION analysis description
<p align="justify"><strong>Purpose</strong>: For signal quantification analysis home-made fillable inserts with radioactive solution are employed for PET, the commercial CIRS density inserts for CT and home-made inserts for MR. In PET analysis, the comparison with respect to the expected activity concentrations is done by computation of the recovery coefficient (RC) for each radioactive insert. In CT analysis, the calibration curve (HU vs. Density) is displayed from different inserts.</p>

<div align="center"><img src="https://github.com/user-attachments/assets/c160273b-4873-46b5-a0e4-7b9db03d5654" width="45%"></div>
<div align="center"><strong>Figure 4. Axial slices of PET (left) and CT (right) images showing radioactive and density inserts.</strong></div>
<p>&nbsp;</p>

<p align="justify"><strong>Prerequisites</strong>: For PET studies, it is required to know the <strong>initial activity (Bq)</strong>, <strong>time when activity was measured (HHMMSS)</strong>, <strong>the volume (ml)</strong> into which the activity has been injected, and the <strong>radiopharmaceutical employed (18F-FDG or Other)</strong>. For CT/MR studies, it’s required to specify the <strong>inserts</strong> used and its <strong>disposition</strong> (if it is different from default).</p>

<p align="justify"><strong>How it works</strong>: When the code is executed, a window appears (Figure 5). It asks us to select the work folder (<em>Analysis folder</em>) with the “Browse directory”. Once chosen, we click on "Run".</p>

<div align="center"><img src="https://github.com/user-attachments/assets/067302d4-9fee-416f-b922-3adc292692f0" width="35%"></div>
<div align="center"><strong>Figure 5. Initial window of Quantification analysis.</strong></div>
<p>&nbsp;</p>

<p align="justify">Within the folder of each measurement, the code goes through each image (PET, CT or MR) and a window where we can visualise the different slices of the image appears (Figure 6). We choose the initial slice where the 7 active inserts (not density inserts, which are the four central inserts) start to be seen by pressing "enter" on the keyboard / "Select" button and then closing the window. Thus, the initial slice of the image analysis region will have been chosen. </p>

- <p align="justify">If a <strong>CT or MR image</strong> folder is found, it asks us which solid inserts were used (for Calibration Curve or for Materials Characterization).</p>

<div align="center"><img src="https://github.com/user-attachments/assets/526b5435-ec61-4457-9db7-84e2c1f0390c" width="45%"></div>
<div align="center"><strong>Figure 6. Window to choose the type of solid inserts used for CT/MR quantification analysis.</strong></div>
<p>&nbsp;</p>

<p align="justify">Once chosen, you will have to indicate whether the disposition of the inserts is the one shown by default or if it is another one. In the latter case, you will be able to modify the disposition, as shown in Figure 7. If you enter a new position-insert pair but you don’t click in “Save” button, the modification will not take effect. When the desired layout appears, click on “Ok” button.</p>

<div align="center"><img src="https://github.com/user-attachments/assets/76cd8c60-7b77-4f7d-ad3b-3dcfb103db29" width="80%"></div>
<div align="center"><strong>Figure 7. Window showing the default disposition of the inserts (left). If it is different from the default, another window will appear allowing to modify the default disposition (right).</strong></div>
<p>&nbsp;</p>

-  <p align="justify">If a <strong>PET image</strong> folder is found, it asks us to specify the parameters of the PET study, and if we choose "Other" radiopharmaceutical that is not 18F-FDG, we will be asked to specify its semi decay period (seconds), as shown in Figure 8.</p>

<div align="center"><img src="https://github.com/user-attachments/assets/93276323-e4fd-4208-8b88-69127eec2a7f" width="55%"></div>
<div align="center"><strong>Figure 8. Windows to specify PET study parameters.</strong></div>
<p>&nbsp;</p>

<p align="justify">The analysis <strong>OUTPUT</strong> is an Excel file with the results of the analysis (average intensities of the inserts from CT, MR and PET images, and recovery coefficients for PET image), as well as a plot of the PET RCs of each active insert and a plot of CT calibration curve.</p>


## RESOLUTION analysis description
<p align="justify"><strong>Purpose</strong>: PET, CT and MR images resolution analysis using 3D-printed phantom with an insert with 6-cylinder sections of different diameters. Resolution analysis of each section by calculating contrast, line profiles and cylinders’ number.</p>

<div align="center"><img src="https://github.com/user-attachments/assets/51a41d04-3700-45ac-ad8b-36d84035433a" width="85%"></div>
<div align="center"><strong>Figure 9. Axial slices of PET (left), CT (centre-left), MR T2 (centre-right) and MR T1 (right) images showing resolution insert.</strong></div>
<p>&nbsp;</p>

<p align="justify"><strong>How it works</strong>: When running the code, a window appears asking us to select the work folder (<em>Analysis folder</em>) and to specify the name of the output excel file. Once chosen, we click on "Run".</p>

<div align="center"><img src="https://github.com/user-attachments/assets/1751f13c-10fe-4b4a-8755-8f1a9043ef0d" width="35%"></div>
<div align="center"><strong>Figure 10. Initial window of Resolution analysis.</strong></div>
<p>&nbsp;</p>

<p align="justify">For both modalities it is necessary to specify the central slice of the region of interest, i.e. the central slice of the resolution insert. When the code passes through each image, a window appears where we can visualise the different slices of the image. We choose the central resolution insert slice by pressing "enter" on the keyboard / "Select" button and then closing the window. Thus, the central slice of the analysis region of the image will have been chosen.</p>
<p align="justify">Once the corresponding image slice has been chosen, 6 windows open to choose the vertices of the 6 triangular sections to be analysed (Figure 11). In each window you have to choose the 3 vertices by clicking on the point where you want to place them (they will be shown in red in the slice). Once the 3 vertices have been chosen, close the window, and the window corresponding to the next section will appear.</p>

<div align="center"><img src="https://github.com/user-attachments/assets/edbcd29c-d6b3-4467-801f-bdbd5476799f" width="85%"></div>
<div align="center"><strong>Figure 11. Selection of the triangular vertices of each section.</strong></div>
<p>&nbsp;</p>

<p align="justify">The analysis <strong>OUTPUT</strong> is an Excel file with the results of the analysis of all input images. In each Excel there is an analysis of an image. Figure 12 shows the results of one of these excels, corresponding to PET image.

<div align="center"><img src="https://github.com/user-attachments/assets/42cc9642-03e6-4f94-9825-8c482ef36d0e" width="70%"></div>
<div align="center"><strong>Figure 12. Output resolution excel results (up) and output intensity profile of one section (down) of a PET image.</strong></div>
<p>&nbsp;</p>


## CO-REGISTRATION analysis description
<p align="justify"><strong>Purpose</strong>: PET/CT and PET/MR images co-registration analysis using 3D-printed phantom with a radiactive cylindrical insert. Dice-Simmilarity-Coefficient (DSC) obtaining between modalities for each axial slice.</p>

<div align="center"><img src="https://github.com/user-attachments/assets/b46636b8-2e52-4144-9338-9fd97d3b7040" width="80%"></div>
<div align="center"><strong>Figure 13. PET/CT images (left) and PET/MR images (right) showing co-registration cylindrical insert in red.</strong></div>
<p>&nbsp;</p>

<p align="justify"><strong>How it works</strong>: When the code is executed, a window appears asking to select the work folder (<em>Analysis folder</em>) and the system type (PET/CT or PET/MR).</p>

<div align="center"><img src="https://github.com/user-attachments/assets/5cf480d0-a610-4d92-b89f-6333158bfca9" width="35%"></div>
<div align="center"><strong>Figure 14. Initial window of Co-registration analysis.</strong></div>
<p>&nbsp;</p>

<p align="justify">Once the folder is selected, we run the code, which goes through the folders of the different dates with the pairs of PET/CT or PET/MR images.</p>
<p align="justify">The anaylisis <strong>OUTPUT</strong> is a png file showing the DSC value for each slice and the mean, as shown in Figure 15.</p>

<div align="center"><img src="https://github.com/user-attachments/assets/f00467d8-7967-4e7f-bf7e-f04a80ffb1d5" width="65%"></div>
<div align="center"><strong>Figure 15. Output DSC plot for all insert slices.</strong></div>
<p>&nbsp;</p>


## DISTORTION analysis description
<p align="justify"><strong>Purpose</strong>: Distortion analysis in CT and MR images using a phantom with a cylindrical distortion insert and a panal distortion insert. Analysis of the constancy of the sides of the cylindrical insert meshes, the constancy of the axial distance between the meshes and the constancy of the panal insert walls.</p>

<div align="center"><img src="https://github.com/user-attachments/assets/2b3de31f-7358-4456-a488-32a1db08e7ce" width="85%"></div>
<div align="center"><strong>Figure 16. Cylindrical insert (left) and panal insert (centre-left) in CT image, cylindrical insert for MR T2 image (centre-right) and MR T1 image (right).</strong></div>
<p>&nbsp;</p>

<p align="justify"><strong>How it works</strong>: When the code is executed, a window appears asking to select the work folder (<em>Analysis folder</em>).</p>

<div align="center"><img src="https://github.com/user-attachments/assets/33bfc9b3-c341-4396-a3d3-514761da5f78" width="35%"></div>
<div align="center"><strong>Figure 17. Initial window of Distortion analysis.</strong></div>
<p>&nbsp;</p>

<p align="justify">Once the data has been entered, we click on “Run” and the code starts to work. First, a window appears asking if the cylindrical distortion insert has been filled.</p>

<div align="center"><img src="https://github.com/user-attachments/assets/167a64dc-31b8-4b3a-b8a7-7c40f721867c" width="40%"></div>
<div align="center"><strong>Figure 18. Cyindrical insert filled window.</strong></div>
<p>&nbsp;</p>

<p align="justify">It is necessary to specify the initial and final slice of the cylindrical insert mesh; and the panal insert slice for each image to analyse. For <strong>CT images</strong>, in the first window, the code asks us to choose the two slices of the cylindrical insert. To do this, we can go through the different slices and once we find the desired slice, click on the "Select" button (twice: one for the initial and one for the final) and then close the window. Then, a second window appears in which we choose in the same way the slice that includes the panal insert.</p>

<p align="justify">For <strong>MR images</strong>, the first and second windows ask us for the same thing as for CT: the initial and final slices of the distortion cylinder, and the central slice of the panal insert. But a third window appears that asks you to specify by click, the walls of the distortion panal in 6 different directions (right, left, up-right diagonal, down-left diagonal, up-left diagonal and down-right diagonal), as shown in Figure 19. Once you have selected the walls in the corresponding direction, close the window, and the window corresponding to the next direction will appear.</p>

<div align="center"><img src="https://github.com/user-attachments/assets/4eef788d-a1c5-4d2c-9800-1285081c3b93" width="85%"></div>
<div align="center"><strong>Figure 19. Window for choosing panal walls in different directions in MR T1 image.</strong></div>
<p>&nbsp;</p>

<p align="justify">The analysis <strong>OUTPUT</strong> are 3 excel files per image with the data of: distortion of the cylindrical insert in the axial axis (DistortionModality_AxialDistances.xlsx), in the xy-plane (DistortionModality_MeshSides.xlsx) and transaxial distortion of the panal insert in different directions (DistortionModality_PanalWalls.xlsx).</p>

<p align="justify">The following figure shows a representation of the results of the distortion analysis of the cylindrical insert for CT image:</p>

<div align="center"><img src="https://github.com/user-attachments/assets/b2a725db-6500-4dd5-b569-320225a29906" width="75%"></div>
<div align="center"><strong>Figure 20. Results of the distortion analysis of the cylindrical insert (CT image).</strong></div>
<p>&nbsp;</p>

<p align="justify">The following figure shows a graphical representation of the results of the distortion analysis of the panal insert for CT image:</p>

<div align="center"><img src="https://github.com/user-attachments/assets/85731ae9-c6d5-49f0-a897-250c92fa8416" width="75%"></div>
<div align="center"><strong>Figure 21. Wall distances of the panal insert in different directions (CT image).</strong></div>
<p>&nbsp;</p>


## RADIOMICS analysis description
<p align="justify"><strong>Purpose</strong>: Radiomic features analysis in CT, PET and MR images using the resolution insert of the 3D-printed phantom.</p>

<div align="center"><img src="https://github.com/user-attachments/assets/708c2a1b-2694-486a-8d5e-12d84828255d" width="75%"></div>
<div align="center"><strong>Figure 22. CT (left), PET (centre) and MR (right) reference images showing radiomic reference segmentations in red.</strong></div>
<p>&nbsp;</p>

<p align="justify"><strong>Prerequisites</strong>: The same folder organisation as in Figure 3, with the difference that in this case we also create <em>Segmentations</em> folder, in which we place different subfolders with images and segmentations of reference. In this folder, we will also introduce a <em>Params.yaml</em> file with some configuration parameters for PET radiomic analysis, as shown in Figure 23.</p>

<div align="center"><img src="https://github.com/user-attachments/assets/695f5cba-1114-4590-9712-116585440478" width="85%"></div>
<div align="center"><strong>Figure 23. Folder organisation for Radiomics analysis.</strong></div>
<p>&nbsp;</p>

<p align="justify"><strong>How it works</strong>: When running the code, a window appears asking us to select the work folder (<em>Analysis folder</em>). Once chosen, we click on "Run".</p>

<div align="center"><img src="https://github.com/user-attachments/assets/8a53e339-4c77-4893-b0bb-7284e1abe491" width="35%"></div>
<div align="center"><strong>Figure 24. Initial window of Radiomics analysis.</strong></div>
<p>&nbsp;</p>

<p align="justify">For all image types, when segmentation registration to analysis image finishes, a window appears showing registered segmentations and then, at command line, it asks us if it is all okay with registered segmentations or if we want to provide our own segmentations (provide the folder path).</p>

<div align="center"><img src="https://github.com/user-attachments/assets/2aa92e25-7055-4285-86e3-c8c12d12c597" width="75%"></div>
<div align="center"><strong>Figure 25. Registered segmentations plot (up) and question to define the segmentations to be used (down).</strong></div>
<p>&nbsp;</p>

<p align="justify">The analysis <strong>OUTPUT</strong> is an excel document for each image with radiomic features in each segmentation.</p>
<p>&nbsp;</p>

**Contact:** <carmen_salvador@iislafe.es>
