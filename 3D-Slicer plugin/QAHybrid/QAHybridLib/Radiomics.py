# Radiomics features extraction for PET, CT and MR DICOM images: 3D Printed Phantom.

import tkinter as tk
from tkinter import filedialog


window = tk.Tk()


class Path:
    def __init__(self):
        self.directory = "None"
        self.directory_path = tk.StringVar()

    def browse_directory(self):
        self.directory = filedialog.askdirectory()
        self.directory_path.set(self.directory)


main = Path()


# Main code
def main_code():
    import os
    import numpy as np
    import SimpleITK as sitk
    import pandas as pd
    import pydicom
    import Auxiliar_functions

    global window, main

    # Folders' paths
    main_path = main.directory  # main folder path
    images_main_path = os.path.join(main_path, 'Images')  # images' folder path
    seg_ref_main_path = os.path.join(main_path, 'Segmentations')  # segmentations' folder path
    params_path = os.path.join(seg_ref_main_path, 'Pruebas/Params.yaml')    # radiomic parameter file path

    # Create an Output folder
    output_path = os.path.join(main_path, 'Output')
    try:
        os.makedirs(output_path)
    except FileExistsError:
        pass

    # Obtains the names of the date folders contained in Images' folder
    sub_dirs = [d for d in os.listdir(images_main_path) if os.path.isdir(os.path.join(images_main_path, d))]

    # Goes through date folders
    for dir in sub_dirs:
        path_date = os.path.join(images_main_path, dir)
        sub_folders = [d for d in os.listdir(path_date) if os.path.isdir(os.path.join(path_date, d))]  # images folders names

        # Create an Output folder for each date
        output_date_path = os.path.join(output_path, dir)
        try:
            os.makedirs(output_date_path)
        except FileExistsError:
            pass

        # Create an Output folder for the type of analysis
        output_date_analysis_path = os.path.join(output_date_path, 'Radiomics')
        try:
            os.makedirs(output_date_analysis_path)
        except FileExistsError:
            pass

        # Goes through the image folders contained in the date folder
        for sub_folder in sub_folders:
            if sub_folder == 'PET' or sub_folder == 'CT' or sub_folder == 'MR':
                if sub_folder == 'PET':
                    print(f'Radiomics analysis of {sub_folder} image dated {dir}...')
                    path_image = os.path.join(path_date, sub_folder)
                    seg_ref_path = os.path.join(seg_ref_main_path, 'PET')
                    seg_ref_paths = [os.path.join(seg_ref_path, d) for d in os.listdir(seg_ref_path) if d != 'PET.nrrd']
                    image_ref_path = os.path.join(seg_ref_path, 'PET.nrrd')

                    Auxiliar_functions.radiomics_analysis(path_image, image_ref_path, seg_ref_paths, sub_folder, dir, output_date_analysis_path, params_path)

                elif sub_folder == 'CT':
                    print(f'Radiomics analysis of {sub_folder} image dated {dir}...')
                    path_image = os.path.join(path_date, sub_folder)
                    seg_ref_path = os.path.join(seg_ref_main_path, 'CT')
                    seg_ref_paths = [os.path.join(seg_ref_path, d) for d in os.listdir(seg_ref_path) if d != 'CT.nrrd']
                    image_ref_path = os.path.join(seg_ref_path, 'CT.nrrd')

                    Auxiliar_functions.radiomics_analysis(path_image, image_ref_path, seg_ref_paths, sub_folder, dir, output_date_analysis_path, params_path)

                elif sub_folder == 'MR':
                    path_images = os.path.join(path_date, sub_folder)
                    images_path = [os.path.join(path_images, d) for d in os.listdir(path_images)]
                    seg_ref_path = os.path.join(seg_ref_main_path, 'MR')
                    seg_ref_paths = [os.path.join(seg_ref_path, d) for d in os.listdir(seg_ref_path) if d != 'MR.nrrd']
                    image_ref_path = os.path.join(seg_ref_path, 'MR.nrrd')

                    for image_path in images_path:
                        # Read image type (T1, T2...) from DICOM header
                        rm_file_names = os.listdir(image_path)
                        first_file_path = os.path.join(image_path, rm_file_names[0])
                        ds = pydicom.filereader.dcmread(first_file_path)
                        image_type = ds["SeriesDescription"].value
                        print(f'Radiomics analysis of {sub_folder + image_type} image dated {dir}...')
                        if 'T1' in image_type:
                            img_name = 'MRT1'
                        elif 'T2' in image_type:
                            img_name = 'MRT2'
                        Auxiliar_functions.radiomics_analysis(image_path, image_ref_path, seg_ref_paths, img_name, dir, output_date_analysis_path, params_path)

    print('Radiomics analysis completed.')

    window.destroy()


window.title('Radiomics analysis')

e1 = tk.Button(window, text="Browse directory", command=main.browse_directory)
lbl1 = tk.Label(window, textvariable=main.directory_path)

tk.Label(window, text='Select analysis directory:').grid(row=0, column=0)
e1.grid(row=0, column=1)
lbl1.grid(row=0, column=2)

tk.Button(window, text='Run', command=main_code).grid(row=1, column=1, sticky=tk.W, pady=4)

tk.mainloop()
