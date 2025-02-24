# Resolution analysis for PET, CT and MR DICOM images: 3D-Printed Phantom.

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
    import pandas as pd
    import SimpleITK as sitk
    import pydicom
    import Auxiliar_functions


    global window, main

    # Folders' paths
    main_path = main.directory  # main folder path
    images_main_path = os.path.join(main_path, 'Images')  # images' folder path

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
        sub_folders = [d for d in os.listdir(path_date) if os.path.isdir(os.path.join(path_date, d))]   # names of the images folders

        # Create an Output folder for each date
        output_date_path = os.path.join(output_path, dir)
        try:
            os.makedirs(output_date_path)
        except FileExistsError:
            pass

        # Create an Output folder for the type of analysis
        output_date_analysis_path = os.path.join(output_date_path, 'Resolution')
        try:
            os.makedirs(output_date_analysis_path)
        except FileExistsError:
            pass

        # Choose central slice of PET, CT and MR analysis area, and the vertices of triangular sections
        for sub_folder in sub_folders:
            if sub_folder == 'PET':
                coord_pet = Auxiliar_functions.select_slices(os.path.join(path_date, sub_folder), sub_folder, 'Choose central slice for resolution analysis', dir)[0]
                vertices_pet = []
                for i in range(6):
                    coordinates_pet = Auxiliar_functions.choose_vertices(os.path.join(path_date, sub_folder), coord_pet, sub_folder, 'Choose triangle vertices for all sections', str(i + 1))
                    vertices_pet.append(coordinates_pet)

            elif sub_folder == 'CT':
                coord_ct = Auxiliar_functions.select_slices(os.path.join(path_date, sub_folder), sub_folder, 'Choose central slice for resolution analysis', dir)[0]
                vertices_ct = []
                for i in range(6):
                  coordinates_ct = Auxiliar_functions.choose_vertices(os.path.join(path_date, sub_folder), coord_ct, sub_folder, 'Choose triangle vertices for all sections', str(i + 1))
                  vertices_ct.append(coordinates_ct)

            elif sub_folder == 'MR':
                mr_images_path = [os.path.join(os.path.join(path_date, sub_folder), d) for d in os.listdir(os.path.join(path_date, sub_folder))]
                for mr_image_path in mr_images_path:
                    rm_file_names = os.listdir(mr_image_path)
                    first_file_path = os.path.join(mr_image_path, rm_file_names[0])  # first dcm path
                    ds = pydicom.filereader.dcmread(first_file_path)
                    type_mr = ds["SeriesDescription"].value
                    if "T1" in type_mr:
                        coord_mrt1 = Auxiliar_functions.select_slices(mr_image_path, sub_folder + "_" + type_mr, 'Choose central slice for resolution analysis', dir)[0]
                        vertices_mrt1 = []
                        for i in range(6):
                            coordinates_mrt1 = Auxiliar_functions.choose_vertices(mr_image_path, coord_mrt1, sub_folder + "_" + type_mr, 'Choose triangle vertices for all sections', str(i + 1))
                            vertices_mrt1.append(coordinates_mrt1)
                    elif "T2" in type_mr:
                        coord_mrt2 = Auxiliar_functions.select_slices(mr_image_path, sub_folder + "_" + type_mr, 'Choose central slice for resolution analysis', dir)[0]
                        vertices_mrt2 = []
                        for i in range(6):
                            coordinates_mrt2 = Auxiliar_functions.choose_vertices(mr_image_path, coord_mrt2, sub_folder + "_" + type_mr, 'Choose triangle vertices for all sections', str(i + 1))
                            vertices_mrt2.append(coordinates_mrt2)

        # Goes through the image folders contained in the date folder
        for sub_folder in sub_folders:
            path_image_type = os.path.join(path_date, sub_folder)

            if sub_folder == 'PET':
                print(f'Resolution analysis of {sub_folder} image dated {dir}...')

                pet_image_path = path_image_type  # DICOM serie folder path

                # Extract image metadata and convert to NIFTI format
                reader_pet = sitk.ImageSeriesReader()
                dicom_names_pet = reader_pet.GetGDCMSeriesFileNames(pet_image_path)
                reader_pet.SetFileNames(dicom_names_pet)
                serie_pet = reader_pet.Execute()

                nifti_path = os.path.join(output_date_analysis_path, sub_folder + '_image.nii')
                sitk.WriteImage(serie_pet, nifti_path)

                im = sitk.ReadImage(nifti_path)

                # Resolution analysis #
                pet_contrast, pet_holes, pet_rc_holes, pet_sections = Auxiliar_functions.resol_pet(im, coord_pet, vertices_pet, output_date_analysis_path, serie_pet.GetSpacing())

                # Save Excel results file
                for n in range(0, len(pet_contrast)):
                    if len(pet_sections) != 0 and len(pet_contrast) != 0:
                        dic_pet = {'Section': pet_sections}
                        dic_pet['Diameter (mm)'] = [5, 7.5, 9, 11, 12, 15]
                        dic_pet['Contrast'] = pet_contrast
                        dic_pet['Holes number'] = pet_holes
                        dic_pet['RC holes'] = pet_rc_holes
                        df_pet = pd.DataFrame(data=dic_pet)
                        writer = pd.ExcelWriter(os.path.join(output_date_analysis_path, f'ResolutionPET.xlsx'))
                        df_pet.to_excel(writer, sheet_name='Resolution', index=False)
                        dic_pet.clear()
                        writer.save()

            elif sub_folder == 'CT':
                print(f'Resolution analysis of {sub_folder} image dated {dir}...')

                ct_image_path = path_image_type

                # Extract image metadata and convert to NIFTI format
                reader = sitk.ImageSeriesReader()
                dicom_names = reader.GetGDCMSeriesFileNames(ct_image_path)
                reader.SetFileNames(dicom_names)
                serie = reader.Execute()
                nifti_path = os.path.join(output_date_analysis_path, sub_folder + '_image.nii')
                sitk.WriteImage(serie, nifti_path)

                im = sitk.ReadImage(nifti_path)

                # Resolution analysis #
                ct_contrast, ct_holes, ct_rc_holes, ct_sections = Auxiliar_functions.resol_ct(im, coord_ct, vertices_ct, output_date_analysis_path, serie.GetSpacing())

                # Save Excel results file
                for n in range(0, len(ct_contrast)):
                    if len(ct_sections) != 0 and len(ct_contrast) != 0:
                        dic_ct = {'Section': ct_sections}
                        dic_ct['Diameter (mm)'] = [5, 7.5, 9, 11, 12, 15]
                        dic_ct['Contrast'] = ct_contrast
                        dic_ct['Holes number'] = ct_holes
                        dic_ct['RC holes'] = ct_rc_holes
                        df_ct = pd.DataFrame(data=dic_ct)
                        writer = pd.ExcelWriter(os.path.join(output_date_analysis_path, f'ResolutionCT.xlsx'))
                        df_ct.to_excel(writer, sheet_name='Resolution', index=False)
                        dic_ct.clear()
                        writer.save()

            elif sub_folder == 'MR':
                print(f'Resolution analysis of {sub_folder} images dated {dir}...')

                mr_path = path_image_type  # MR images' folder path
                mr_images_path = [os.path.join(mr_path, d) for d in os.listdir(mr_path)]  # DICOM series folder path
                for path_mr_images in mr_images_path:
                    # Extract image metadata and convert to NIFTI format
                    reader_mr = sitk.ImageSeriesReader()
                    dicom_names_mr = reader_mr.GetGDCMSeriesFileNames(path_mr_images)
                    reader_mr.SetFileNames(dicom_names_mr)
                    serie_mr = reader_mr.Execute()

                    rm_file_names = os.listdir(path_mr_images)
                    first_file_path = os.path.join(path_mr_images, rm_file_names[0])  # first dcm path
                    ds = pydicom.filereader.dcmread(first_file_path)
                    type_mr = ds["SeriesDescription"].value
                    nifti_path = os.path.join(output_date_analysis_path, sub_folder + str(type_mr) + '_image.nii')
                    sitk.WriteImage(serie_mr, nifti_path)

                    im = sitk.ReadImage(nifti_path)

                    # Resolution analysis #
                    if "T1" in type_mr:
                        # N4BiasFieldCorrection
                        #im = sitk.ReadImage(nifti_path, sitk.sitkFloat32)
                        #corrector = sitk.N4BiasFieldCorrectionImageFilter()
                        #corrector.SetMaximumNumberOfIterations([1000])
                        #corrector.SetSplineOrder(200)
                        #im_filtered = corrector.Execute(im)
                        #im_filtered.SetDirection(im.GetDirection())
                        #im_filtered.SetSpacing(im.GetSpacing())
                        #im_filtered.SetOrigin(im.GetOrigin())

                        #mr_contrast, mr_holes, mr_sections = Auxiliar_functions.resol_mr_t1(im_filtered, coord_mrt1, vertices_mrt1, output_date_analysis_path, serie_mr.GetSpacing()[1])
                        img_name = 'MRT1'
                        mr_contrast, mr_holes, mr_rc_holes, mr_sections = Auxiliar_functions.resol_mr_t1(im, coord_mrt1, vertices_mrt1, output_date_analysis_path, serie_mr.GetSpacing())

                    elif "T2" in type_mr:
                        img_name = 'MRT2'
                        mr_contrast, mr_holes, mr_rc_holes, mr_sections = Auxiliar_functions.resol_mr_t2(im, coord_mrt2, vertices_mrt2, output_date_analysis_path)

                    for n in range(0, len(mr_contrast)):
                        if len(mr_sections) != 0 and len(mr_contrast) != 0:
                            dic_mr = {'Section': mr_sections}
                            dic_mr['Diameter (mm)'] = [5, 7.5, 9, 11, 12, 15]
                            dic_mr['Contrast'] = mr_contrast
                            dic_mr['Holes number'] = mr_holes
                            dic_mr['RC holes'] = mr_rc_holes
                            df_mr = pd.DataFrame(data=dic_mr)
                            writer = pd.ExcelWriter(os.path.join(output_date_analysis_path, f'Resolution{img_name}.xlsx'))
                            df_mr.to_excel(writer, sheet_name='Resolution', index=False)
                            dic_mr.clear()
                            writer.save()

    print('Resolution analysis completed.')

    window.destroy()


window.title('Resolution analysis')

e1 = tk.Button(window, text="Browse directory", command=main.browse_directory)
lbl1 = tk.Label(window, textvariable=main.directory_path)

tk.Label(window, text='Select analysis directory:').grid(row=0, column=0)
e1.grid(row=0, column=1)
lbl1.grid(row=0, column=2)

tk.Button(window, text='Run', command=main_code).grid(row=1, column=1, sticky=tk.W, pady=4)

tk.mainloop()
