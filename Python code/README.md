# 3D-printed phantom Quality Assurance: Python Original Code
@Author: __Carmen Salvador Ribés__

## Overview
<p style="text-align: justify;">The 3D-printed phantom Quality Assurance code enables the analysis of different quality parameters of PET, CT and MR images (quantification, resolution, co-registration, distortion and radiomics) with 3D-printed phantom acquisitions.
The code is made up of 7 scripts that are shown in Figure 1. To start the analysis, the script "3DPrintedPhantom_QA.py" must be run, which allows to choose the type of quality analysis to be performed (Figure 2). The script "Auxiliar_functions.py" contains different functions that are used in the different analyses. The rest of the scripts perform the analysis corresponding to its name.</p>

![image](https://github.com/user-attachments/assets/eac58c31-2c3e-42e0-9891-dab36b568134)
 **Figure 1. QA 3D-printed phantom code scripts.**
 
 ![image](https://github.com/user-attachments/assets/c4838fbd-1f25-4534-9136-10640157f652)
**Figure 2. Initial window that appears when running the script 3DPrintedPhantom_QA.py.**

<u>**General prerequisites**</u>: It is necessary to have the PET, CT and MR images to be analysed in DICOM format.     
A specific folder organisation is required for the code to work (Figure 3). The working folder is created (in this case the Analysis folder), where we create a subfolder called Images in which we will place the folders with the different dates of the studies. Within these folders, we will place the DICOM series corresponding to an image in a folder with the name of the modality (CT, PET and MR). As an exception, in the MR folder there shall be subfolders corresponding to the possible existing sequences (as an example, two sequences T1 and T2), which shall contain the corresponding DICOM series. It is not necessary that the folders have these two names, they can have any name, the code parses the sequence type from the DICOM header. The Output folder will be automatically generated in the analysis process within the work folder. In it, different output files will be generated within the corresponding date folder


Contact: <carmen_salvador@iislafe.es>
