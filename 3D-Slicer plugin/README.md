# 3D-printed phantom Quality Assurance: QAHybrid Plugin
@Authors: __Tobias Fechter__, __Carmen Salvador Ribés__ and __Montserrat Carles Fariña__

<div align="center"><img src="https://github.com/user-attachments/assets/36e86fcb-5353-4b82-81da-0bf4d38d1def" width="45%"></div>

## First steps
<p align="justify">The <strong>QAHybrid</strong> plugin was validated in 3D-Slicer version 5.6.2. To install it:</p>
<p align="justify">1. Download the plugin folder <a href="https://github.com/carmsr/3DPrintedPhantom-QualityAssurance/tree/main/3D-Slicer%20plugin/QAHybrid" target="_blank"> <em>3DPrintedPhantom-QualityAssurance/3D-Slicer plugin/QAHybrid</em></a>.</p>
<p align="justify">2. Open the 3D-Slicer program, and go to the menu ‘Edit - Application Settings’.</p>
<div align="center"><img src="https://github.com/user-attachments/assets/9c71d87b-483c-40db-954c-bee3507c6d3e" width="35%"></div>

<p align="justify">3. In the ‘Modules’ section, add the path to the plugin folder in the 'Additional module paths' section, and restart 3D-Slicer.</p>
<div align="center"><img src="https://github.com/user-attachments/assets/bee5582b-4ed0-4aa2-8a84-0859a2365191" width="85%"></div>

<p align="justify">4. The ‘QAHybrid’ module will appear in the ‘Quantification’ category.</p>
<div align="center"><img src="https://github.com/user-attachments/assets/1dc8c46b-4c1e-4ec3-9ba3-8c0e12932279" width="55%"></div>

## Analysis
### General
1.	In 3D-Slicer, import the DICOM series to be analysed (PET, CT or MR T1/T2).
2.	Once imported into the DICOM database, upload them into the workspace.
3.	Go to the category ‘Modules - Quantification - QAHybrid’ to access to the analysis module.
4.	Once in the **QAHybrid** module, select the images to be analysed as input.
<div align="center"><img src="https://github.com/user-attachments/assets/bb92db42-904b-4354-9f0c-3c7045d99538" width="45%"></div>
5.	Select the output folder where the results will be saved.
<div align="center"><img src="https://github.com/user-attachments/assets/6dceceee-139e-41d3-af8e-c671cab032db" width="45%"></div>

### Quantification analysis
1. In the ‘Inputs’ section, select the analysis type ‘Quantification’ under ‘Analysis-Type’. The quantification section will open automatically, and the parameters that still need to be specified in order to start the analysis will appear in red.
<div align="center"><img src="https://github.com/user-attachments/assets/b799e0f0-67fd-41e3-b03f-5e8c9e527fbd" width="45%"></div>

2. The first step is the determination of the initial axial slice for the analysis.

   - For **CT and PET imaging**, it can be determined automatically by means of the ‘detect’ button and also manually. It is recommended to check and adjust (if necessary) the results of the automatic detection to ensure correct analysis.
   - For **MR imaging**, it is manual only.
<div align="center"><img src="https://github.com/user-attachments/assets/cc6008f8-d3b6-4298-835b-503b6dc83f1f" width="45%"></div>


3. For <strong>PET quantification</strong>, the parameters are entered manually:
    - Volume of radioactive solution (in mL)
    - Radiotracer used: 18F-FDG or Other
    - Half-life of the radiotracer (in s)
    - Initial activity of the radioactive solution (in Bq)
    - Time of measurement of the initial activity (in HHMMSS format)
<div align="center"><img src="https://github.com/user-attachments/assets/6b1b263d-8eb8-486a-9be9-3db30b468434" width="45%"></div>

4. For <strong>CT and MR quantification</strong>, ‘Select Insert Type’ selects the insert types used: for calibration curve (for CT imaging only) and for material characterisation (for MR or CT imaging). Within each option, the insert arrangement can be specified/modified.
<div align="center"><img src="https://github.com/user-attachments/assets/16e4a41a-78ec-41bf-ab42-c904134a70f8" width="45%"></div>

5. Once all the necessary parameters have been determined/entered, click on the ‘Apply’ button. The analysis will be carried out automatically and the outputs will be obtained. For each image: table with quantification information. For PET image, figure with activity concentration RC values. For CT image, figure with calibration curve. The figures will be automatically saved in the directory initially specified, while the tables will be displayed in the 3D-Slicer interface and can be saved in different formats.

### Resolution analysis
1. In the ‘Inputs’ section, select the analysis type ‘Resolution’ under ‘Analysis-Type’. The resolution section will open automatically, and the parameters that still need to be specified in order to start the analysis will appear in red.
<div align="center"><img src="https://github.com/user-attachments/assets/903bf71e-a587-40fe-8fb1-da743e9170a2" width="45%"></div>

2. Automatic detection of the central axial slice of the resolution insert and automatic generation of the vertices of each triangular region by means of the ‘detect’ button. For all imaging modalities.
<div align="center"><img src="https://github.com/user-attachments/assets/86c276b4-8817-4bf1-9ec3-188b4cbd5163" width="45%"></div>

3. Once the vertices have been generated, it is recommended to check and adjust manually, avoiding including the central cylinder of the insert in the triangular regions. For manual modification, it is necessary that it is first detected automatically. The vertices generated for different images (from left to right: CT, PET and MR T2) are shown below.
<div align="center"><img src="https://github.com/user-attachments/assets/1e7d3274-46d6-4382-9c70-1c5df501e7e8" width="70%"></div>

4. Click on the ‘Apply’ button. The analysis will be carried out automatically and the outputs will be obtained. For each image: figures of the intensity profiles in the different triangular regions and tables with the contrast values. The figures will be automatically saved in the directory initially specified, while the tables will be displayed in the 3D-Slicer interface and can be saved in different formats.

### Co-registration analysis
1. In the ‘Inputs’ section, select the analysis type ‘Co-registration’ under ‘Analysis-Type’.
<div align="center"><img src="https://github.com/user-attachments/assets/77c792ff-4fcb-43dd-86f5-504eff54ddfd" width="45%"></div>

2. Click on the ‘Apply’ button, as the analysis is fully automatic. A figure with DSC values evaluating system co-registration will be generated and automatically saved in specified directory.

### Distortion analysis
1. In the ‘Inputs’ section, select the analysis type ‘Distortion’ under ‘Analysis-Type’. The distortion section will open automatically, and the parameters that still need to be specified in order to start the analysis will appear in red.
<div align="center"><img src="https://github.com/user-attachments/assets/cd202c68-04e5-495f-91c4-e44e914a8788" width="45%"></div>

2. Automatic detection of the initial and final axial slices of distortion cylinder and central axial slice of panal (all modalities), and automatic generation of panal wall points (MR imaging only) by means of the ‘detect’ button.
<div align="center"><img src="https://github.com/user-attachments/assets/bbdb6e19-1c08-4336-a2dc-9eb26f372487" width="45%"></div>

3. For MR images, once the panal wall points have been generated, it is recommended to check and adjust manually. The panal walls points generated for two directons (Diagonal-Up-Left and Diagonal-Down-Right) for MR T1 image are shown below.
<div align="center"><img src="https://github.com/user-attachments/assets/97ba052d-51d9-4d81-8388-651438790ed7" width="25%"></div>

4. Click on the ‘Apply’ button. The analysis will be carried out automatically and the outputs will be obtained. Figures and tables with distortion evaluation will be generated. The figures will be automatically saved in the directory initially specified, while the tables will be displayed in the 3D-Slicer interface and can be saved in different formats.

### Radiomics analysis
1. In the ‘Inputs’ section, select the analysis type ‘Radiomics’ under ‘Analysis-Type’. The radiomics section will open automatically, and the parameters that still need to be specified in order to start the analysis will appear in red.
<div align="center"><img src="https://github.com/user-attachments/assets/70229686-c657-4c4c-b51e-77715bdfe585" width="45%"></div>

2. The segmentations in which the analysis is to be carried out can be provided. Otherwise, the reference segmentations contained in the software shall be used, adapting them to the analysis image space.
<div align="center"><img src="https://github.com/user-attachments/assets/8182eb84-ab87-4058-b2a7-e3b8dd8a9e08" width="45%"></div>

3. Click on the ‘Apply’ button. The analysis will be carried out automatically and the outputs will be obtained. For each image, a table with the radiomics features in each segmentation. The table could be saved in different formats.
<p>&nbsp;</p>

**Contact**: <carmen_salvador@iislafe.es>
