# 3D-printed phantom Quality Assurance: Python Original Code
@Author: __Carmen Salvador Ribés__

## Overview
<p align="justify">The 3D-printed phantom Quality Assurance code enables the analysis of different quality parameters of PET, CT and MR images (quantification, resolution, co-registration, distortion and radiomics) with 3D-printed phantom acquisitions.
The code is made up of 7 scripts that are shown in Figure 1. To start the analysis, the script "3DPrintedPhantom_QA.py" must be run, which allows to choose the type of quality analysis to be performed (Figure 2). The script "Auxiliar_functions.py" contains different functions that are used in the different analyses. The rest of the scripts perform the analysis corresponding to its name.</p>

![image](https://github.com/user-attachments/assets/eac58c31-2c3e-42e0-9891-dab36b568134)

 **Figure 1. QA 3D-printed phantom code scripts.**
 
 ![image](https://github.com/user-attachments/assets/c4838fbd-1f25-4534-9136-10640157f652)
 
**Figure 2. Initial window that appears when running the script 3DPrintedPhantom_QA.py.**

<p align="justify"><strong>General prerequisites</strong>: It is necessary to have the PET, CT and MR images to be analysed in DICOM format.     
A specific folder organisation is required for the code to work (Figure 3). The working folder is created (in this case the Analysis folder), where we create a subfolder called Images in which we will place the folders with the different dates of the studies. Within these folders, we will place the DICOM series corresponding to an image in a folder with the name of the modality (CT, PET and MR). As an exception, in the MR folder there shall be subfolders corresponding to the possible existing sequences (as an example, two sequences T1 and T2), which shall contain the corresponding DICOM series. It is not necessary that the folders have these two names, they can have any name, the code parses the sequence type from the DICOM header. The Output folder will be automatically generated in the analysis process within the work folder. In it, different output files will be generated within the corresponding date folder.</p>

![image](https://github.com/user-attachments/assets/58c8b942-7a7f-4e77-ace2-a5cab78f6ffb)

**Figure 3. Folder organisation for QA 3D-printed phantom analysis.**

## QUANTIFICATION analysis description
<p align="justify"><strong>Purpose</strong>: For signal quantification analysis home-made fillable inserts with radioactive solution are employed for PET, the commercial CIRS density inserts for CT and home-made inserts for MR.</p> 
<p align="justify">In PET analysis, the comparison with respect to the expected concentrations is done by computation of the recovery coefficient (RC) for each radioactive insert.</p>
In CT analysis, the calibration curve (HU vs. Density) is displayed.</p>

![image](https://github.com/user-attachments/assets/c160273b-4873-46b5-a0e4-7b9db03d5654)

**Figure 4. Axial slices of PET (left) and CT (right) images showing radioactive and density inserts.**

<p align="justify"><strong>Prerequisites</strong>: For PET studies, it is required to know the <strong>initial activity (Bq)</strong>, <strong>time when activity was measured (HHMMSS)</strong>, <strong>the volume (ml)</strong> into which the activity has been injected, and the <strong>radiopharmaceutical employed (18F-FDG or Other)</strong>. For CT/MR studies, it’s required to specify the <strong>inserts</strong> used and its <strong>disposition</strong> (if it is different from default).</p>

<p align="justify"><strong>How it works</strong>: When the code is executed, a window appears (Figure 5). It asks us to select the work folder (Analysis folder) with the “Browse directory”. Once chosen, we click on "Run".</p>

![image](https://github.com/user-attachments/assets/067302d4-9fee-416f-b922-3adc292692f0)

**Figure 5. Initial window of Quantification analysis.**

<p align="justify">Within the folder of each measurement, the code goes through each image (PET, CT or MR) and a window where we can visualise the different slices of the image appears (Figure 6). We choose the initial slice where the 7 active inserts (not density inserts, which are the four central inserts) start to be seen by pressing "enter" on the keyboard / "Select" button and then closing the window. Thus, the initial slice of the image analysis region will have been chosen. </p>

- <p align="justify">If a <strong>CT or MR image</strong> folder is found, it asks us which solid inserts were used (for Calibration Curve or for Materials Characterization).</p>

![image](https://github.com/user-attachments/assets/526b5435-ec61-4457-9db7-84e2c1f0390c)

**Figure 6. Window to choose the type of solid inserts used for CT/MR quantification analysis.**

<p align="justify">Once chosen, you will have to indicate whether the disposition of the inserts is the one shown by default or if it is another one. In the latter case, you will be able to modify the disposition, as shown in Figure 8. If you enter a new position-insert pair but you don’t click in “Save” button (right window Fig. 8), the modification will not take effect. When the desired layout appears, click on “Ok” button (left window Fig. 7).</p>

![image](https://github.com/user-attachments/assets/76cd8c60-7b77-4f7d-ad3b-3dcfb103db29)

**Figure 7. Window showing the default disposition of the inserts (left). If it is different from the default, another window will appear allowing to modify the default disposition (right).**

-  <p align="justify">If a <strong>PET image</strong> folder is found, it asks us to specify the parameters of the PET study, and if we choose "Other" radiopharmaceutical that is not 18F-FDG, we will be asked to specify its semi decay period (seconds), as shown in Figure 8.</p>

![image](https://github.com/user-attachments/assets/93276323-e4fd-4208-8b88-69127eec2a7f)

**Figure 8. Windows to specify PET study parameters.**

<p align="justify">The analysis <strong>OUTPUT</strong> is an Excel file with the results of the analysis (average intensities of the inserts from CT, MR and PET images, and recovery coefficients for PET image), as well as a plot of the PET RCs of each active insert and a plot of CT calibration curve.</p>


## RESOLUTION analysis description
<p align="justify"><strong>Purpose</strong>:PET, CT and MR images resolution analysis using 3D-printed phantom with an insert with 6-cylinder sections of different diameters. Resolution analysis of each section by calculating contrast, line profiles and cylinders’ number.</p>
 
![image](https://github.com/user-attachments/assets/51a41d04-3700-45ac-ad8b-36d84035433a)

**Figure 9. Axial slices of PET (left), CT (centre-left), MR T2 (centre-right) and MR T1 (right) images showing resolution insert.**

<p align="justify"><strong>How it works</strong>: When running the code, a window appears asking us to select the work folder (Analysis folder) and to specify the name of the output excel file. Once chosen, we click on "Run".</p>

![image](https://github.com/user-attachments/assets/1751f13c-10fe-4b4a-8755-8f1a9043ef0d)

**Figure 10. Initial window of Resolution analysis.**

<p align="justify">For both modalities it is necessary to specify the central slice of the region of interest, i.e. the central slice of the resolution insert. When the code passes through each image, a window appears where we can visualise the different slices of the image. We choose the central resolution insert slice by pressing "enter" on the keyboard / "Select" button and then closing the window. Thus, the central slice of the analysis region of the image will have been chosen.</p>
<p align="justify">Once the corresponding image slice has been chosen, 6 windows open to choose the vertices of the 6 triangular sections to be analysed (Figure 11). In each window you have to choose the 3 vertices by clicking on the point where you want to place them (they will be shown in red in the slice). Once the 3 vertices have been chosen, close the window, and the window corresponding to the next section will appear.</p>

![image](https://github.com/user-attachments/assets/edbcd29c-d6b3-4467-801f-bdbd5476799f)

**Figure 11. Selection of the triangular vertices of each section.**

<p align="justify">The analysis OUTPUT is an Excel file with the results of the analysis of all input images. In each Excel there is an analysis of an image. Figure 12 shows the results of one of these excels, corresponding to PET image.

![image](https://github.com/user-attachments/assets/42cc9642-03e6-4f94-9825-8c482ef36d0e)

**Figure 12. Output resolution excel results (up) and output intensity profile of one section (down) of a PET image.**


## CO-REGISTRATION analysis description
<p align="justify"><strong>Purpose</strong>: PET/CT and PET/MR images co-registration analysis using 3D-printed phantom with a radiactive cylindrical insert.</p>

![image](https://github.com/user-attachments/assets/b46636b8-2e52-4144-9338-9fd97d3b7040)

**Figure 13. PET/CT images (left) and PET/MR images (right) showing co-registration cylindrical insert in red.**

<p align="justify"><strong>How it works</strong>: When the code is executed, a window appears asking to select the work folder (Analysis folder) and the system type (PET/CT or PET/MR).</p>

![image](https://github.com/user-attachments/assets/5cf480d0-a610-4d92-b89f-6333158bfca9)

**Figure 14. Initial window of Co-registration analysis.**

<p align="justify">Once the folder is selected, we run the code, which goes through the folders of the different dates with the pairs of PET/CT or PET/MR images.</p>
<p align="justify">The final result is a png file showing the Dice-Simmilarity-Coefficient (DSC) for each slice, as shown in Figure 15.</p>

**Figure 15. Output DSC plot for all insert slices.**


## DISTORTION analysis description
<p align="justify"><strong>Purpose</strong>: Distortion analysis in CT and MR images using a phantom with a cylindrical distortion insert and a panal distortion insert. Analysis of the constancy of the sides of the cylindrical insert meshes, the constancy of the axial distance between the meshes and the constancy of the panal insert walls.</p>

![image](https://github.com/user-attachments/assets/2b3de31f-7358-4456-a488-32a1db08e7ce)

**Figure 16. Cylindrical insert (left) and panal insert (centre-left) in CT image, cylindrical insert for MR T2 image (centre-right) and MR T1 image (right).**

<p align="justify"><strong>How it works</strong>: When the code is executed, a window appears asking to select the work folder (Analysis folder).</p>

![image](https://github.com/user-attachments/assets/33bfc9b3-c341-4396-a3d3-514761da5f78)

**Figure 17. Initial window of Distortion analysis.**



**Contact:** <carmen_salvador@iislafe.es>
