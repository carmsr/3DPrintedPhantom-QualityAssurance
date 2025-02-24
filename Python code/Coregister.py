# Co-registration Analysis for PET/CT and PET/MR systems with the 3D-Printed Phantom.

import os
import tkinter as tk
from tkinter import filedialog, ttk


window = tk.Tk()


class Path:
    def __init__(self):
        self.directory = "None"
        self.directory_path = tk.StringVar()

    def browse_directory(self):
        self.directory = filedialog.askdirectory()
        self.directory_path.set(self.directory)


class OptionDialog(tk.Toplevel):
    def __init__(self, parent, title, prompt, options):
        super().__init__(parent)
        self.title(title)
        self.prompt = prompt
        self.options = options
        self.result = None

        self.label = ttk.Label(self, text=prompt)
        self.label.pack(padx=5, pady=5)

        label_width = self.label.winfo_reqwidth()
        window_width = label_width + 50

        self.geometry(f"{window_width}x100")

        self.option_var = tk.StringVar(self)
        self.option_var.set(options[0])

        self.option_menu = ttk.OptionMenu(self, self.option_var, options[0], *options)
        self.option_menu.pack(padx=5, pady=5)

        self.button_ok = ttk.Button(self, text="Ok", command=self.ok)
        self.button_ok.pack(padx=5, pady=5)

    def ok(self):
        self.result = self.option_var.get()
        self.destroy()


main = Path()


# Main code
def main_code():

    import statistics
    import pydicom
    import SimpleITK as sitk
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as cls
    import Auxiliar_functions

    global window, main

    system_type = str(w1.get())     # PET/CT or PET/MR system

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
        sub_folders = [d for d in os.listdir(path_date) if os.path.isdir(os.path.join(path_date, d))]  # names of the images folders

        # Create an Output folder for each date
        output_date_path = os.path.join(output_path, dir)
        try:
            os.makedirs(output_date_path)
        except FileExistsError:
            pass

        # Create an Output folder for the type of analysis
        output_date_analysis_path = os.path.join(output_date_path, 'Co-registration')
        try:
            os.makedirs(output_date_analysis_path)
        except FileExistsError:
            pass

        # Goes through the image folders contained in the date folder
        for sub_folder in sub_folders:
            path_image_type = os.path.join(path_date, sub_folder)
            if sub_folder == 'PET':
                pet_dicom_path = path_image_type
            elif sub_folder == 'CT':
                ct_dicom_path = path_image_type
            elif sub_folder == 'MR':
                mr_path = path_image_type  # MR images' folder path
                mr_images_path = [os.path.join(mr_path, d) for d in os.listdir(mr_path)]  # DICOM series folder path
                for mr_image_path in mr_images_path:
                    rm_file_names = os.listdir(mr_image_path)
                    first_file_path = os.path.join(mr_image_path, rm_file_names[0])  # first dcm path
                    ds = pydicom.filereader.dcmread(first_file_path)
                    type_mr = ds["SeriesDescription"].value
                    if "T1" in type_mr:
                        mr_sequence = "T1"
                        mr_dicom_path = mr_image_path
                    elif "T2" in type_mr:
                        mr_sequence = "T2"
                        mr_dicom_path = mr_image_path

        # Convert DICOM series to NIFTI format
        # PET image
        reader_pet = sitk.ImageSeriesReader()
        dicom_names_pet = reader_pet.GetGDCMSeriesFileNames(pet_dicom_path)
        reader_pet.SetFileNames(dicom_names_pet)
        serie_pet = reader_pet.Execute()
        voxelvol_pet = serie_pet.GetSpacing()[0] * serie_pet.GetSpacing()[1] * serie_pet.GetSpacing()[2]

        pet_image_path = os.path.join(output_date_analysis_path, 'PET_image.nii')
        sitk.WriteImage(serie_pet, pet_image_path)
        pet_im = sitk.ReadImage(pet_image_path)

        # CT or MR image
        if system_type == "PET/CT":
            reader_ct = sitk.ImageSeriesReader()
            dicom_names_ct = reader_ct.GetGDCMSeriesFileNames(ct_dicom_path)
            reader_ct.SetFileNames(dicom_names_ct)
            serie_ct = reader_ct.Execute()
            ct_image_path = os.path.join(output_date_analysis_path, 'CT_image.nii')
            sitk.WriteImage(serie_ct, ct_image_path)
            other_im = sitk.ReadImage(ct_image_path)
        elif system_type == "PET/MR":
            reader_mr = sitk.ImageSeriesReader()
            dicom_names_mr = reader_mr.GetGDCMSeriesFileNames(mr_dicom_path)
            reader_mr.SetFileNames(dicom_names_mr)
            serie_mr = reader_mr.Execute()
            mr_image_path = os.path.join(output_date_analysis_path, 'MR_' + str(mr_sequence) + '_image.nii')
            sitk.WriteImage(serie_mr, mr_image_path)
            other_im = sitk.ReadImage(mr_image_path)
            system_type += str(mr_sequence)

        # PET threshold segmentation: VOI_PET, Ref_VOI_PET
        voi_pt, ref_voi_pt = Auxiliar_functions.pet_th(pet_im, voxelvol_pet)

        # CT/MR threshold segmentation: Threshold_CT/Threshold_MR
        threshold = Auxiliar_functions.ct_mr_th(other_im, system_type)

        # Ref_VOI_PET resampling to CT/MR dimensions: Ref_VOI_CT/Ref_VOI_MR
        ref_voi = Auxiliar_functions.resampling(other_im, ref_voi_pt)

        # Intersection Threshold_CT/Threshold_MR and Ref_VOI_CT/Ref_VOI_MR: VOI_CT/VOI_MR
        th = sitk.GetArrayFromImage(threshold)
        th = th.astype(bool)
        rf_voi = sitk.GetArrayFromImage(ref_voi)
        rf_voi = rf_voi.astype(bool)

        th[~rf_voi] = 0  # intersection
        th = th.astype(int)

        voi = sitk.GetImageFromArray(th)
        voi.SetDirection(threshold.GetDirection())
        voi.SetSpacing(threshold.GetSpacing())
        voi.SetOrigin(threshold.GetOrigin())

        # Intersection VOI_PET and VOI_CT_resampled/VOI_MR_resampled
        # Resample VOI_CT/VOI_MR to PET image dimensions
        voi_resamp = Auxiliar_functions.resampling(pet_im, voi)

        # Intersection
        voi_pet = sitk.GetArrayFromImage(voi_pt)
        voi_resampled = sitk.GetArrayFromImage(voi_resamp)
        intersection = voi_pet * voi_resampled  # intersection

        # DSC calculation #
        # Slices number calculation for DSC analysis
        li_pet = 0
        li_other = 0
        lf_pet = len(voi_pet[:, 0, 0]) - 1
        lf_other = len(voi_resampled[:, 0, 0]) - 1

        # Choose initial and final slice for PET image
        for s in range(0, len(voi_pet[:, 0, 0])):
            if np.mean(voi_pet[s, :, :]) == 0:
                li_pet = li_pet + 1
            else:
                break
        for s in range(li_pet, len(voi_pet[:, 0, 0])):
            if np.mean(voi_pet[s, :, :]) != 0:
                lf_pet = s
            else:
                break

        # Choose initial and final slice for CT image
        for s in range(0, len(voi_resampled[:, 0, 0])):
            if np.mean(voi_resampled[s, :, :]) == 0:
                li_other = li_other + 1
            else:
                break
        for s in range(li_other, len(voi_resampled[:, 0, 0])):
            if np.mean(voi_resampled[s, :, :]) != 0:
                lf_other = s
            else:
                break

        # Choose inferior and superior slice for DSC calculation
        if li_pet < li_other:
            li = li_pet
        else:
            li = li_other
        if lf_pet > li_other:
            lf = lf_pet
        else:
            lf = lf_other

        # DSC calculation
        dsc = []
        z = []
        for j in range(li, lf + 1):
            slice = voi_resampled[j, :, :]
            pet_slice = voi_pet[j, :, :]
            intersect_slice = intersection[j, :, :]
            area_ct = np.sum(slice == 1)
            area_pet = np.sum(pet_slice == 1)
            area_intersect = np.sum(intersect_slice == 1)

            dsccoeff = (2 * area_intersect) / (area_ct + area_pet)
            dsc.append(dsccoeff)
            z.append(j)

        # Get VOI_PET and VOI_CT_resampled/VOI_MR_resampled to calculate distance between CT/MR and PET
        # Load position of origin voxel
        PixelPosition_other = voi_resamp.GetOrigin()
        PixelPosition_PET = voi_pt.GetOrigin()
        PixelSpacing_other = voi_resamp.GetSpacing()
        PixelSpacing_PET = voi_pt.GetSpacing()

        # Difference between first VOI_CT/VOI_MR and VOI_PET slices in mm
        pos_other = PixelPosition_other[2] + PixelSpacing_other[2] * li_other  # Position in mm of the first VOI_CT_resampled/VOI_MR_resampled slice
        pos_PET = PixelPosition_PET[2] + PixelSpacing_PET[2] * li_pet  # Position in mm of the first VOI_PET slice
        DIF_CORTES = abs(pos_PET - pos_other)

        # DSC's plot according to the slice
        mean_dsc = statistics.mean(dsc)
        std_dsc = statistics.stdev(dsc)
        min_dsc = np.min(dsc)
        max_dsc = np.max(dsc)

        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111)
        ax1.plot(z, dsc, color=cls.to_rgba('C7', 0.7))
        ax1.axhline(y=mean_dsc, color=cls.to_rgba('C9', 0.8))
        ax1.axhline(y=0.91, color=cls.to_rgba('C10', 0.8), linestyle='-')
        ax1.fill_between(z, dsc, color=cls.to_rgba('C7', 0.3))
        ax1.axhline(y=1, color='k', linestyle='dashed')
        ax1.set_xlim([min(z), max(z)])
        ax1.set_ylim([0, 1.05])
        ax1.set_xlabel('Slice', fontsize=15)
        ax1.set_ylabel('DSC', fontsize=15)

        if round(mean_dsc, 2) >= 0.91:
            if DIF_CORTES > 2:
                color = ('black', 'green', 'black')
                ax1.legend([f'DSC coefficient: The {system_type} system is not properly co-registered', f'Mean DSC= {round(mean_dsc, 2)} \u00B1 {round(std_dsc, 2)} and $\Delta$z={round(DIF_CORTES, 2)} mm', f'DSC threshold value: 0.91'], labelcolor=color, title=f'Co-registration properties')
            else:
                color = ('black', 'green', 'black')
                ax1.legend([f'DSC coefficient: The {system_type} system is well co-registered', f'Mean DSC= {round(mean_dsc, 2)} \u00B1 {round(std_dsc, 2)} and $\Delta$z={round(DIF_CORTES, 2)} mm', f'DSC threshold value: 0.91'], labelcolor=color, title=f'Co-registration properties')
        else:
            color = ('black', 'red', 'black')
            ax1.legend([f'DSC coefficient: The {system_type} system is not properly co-registered', f'Mean DSC= {round(mean_dsc, 2)} \u00B1 {round(std_dsc, 2)} and $\Delta$z={round(DIF_CORTES, 2)} mm', f'DSC threshold value: 0.91'], labelcolor=color, title=f'Co-registration properties')

        loc_min = []
        loc_max = []
        for i in range(0, len(dsc)):
            if dsc[i] == min_dsc:
                loc_min.append(i)
            if dsc[i] == max_dsc:
                loc_max.append(i)

        if len(loc_min) == 1:
            ax1.plot(loc_min[0] + min(z), min_dsc, marker=".", color="red")
            ax1.text(loc_min[0] + min(z), min_dsc - 0.03 * min_dsc, round(min_dsc, 4))

        if len(loc_min) != 1 and len(loc_min) != 0:
            ax1.plot(loc_min[0] + min(z), min_dsc, marker=".", color="red")
            ax1.text(loc_min[0] + min(z), min_dsc - 0.03 * min_dsc, round(min_dsc, 4))
            for i in range(1, len(loc_min)):
                ax1.plot(loc_min[i] + min(z), min_dsc, marker=".", color="red")

        if len(loc_max) == 1:
            ax1.plot(loc_max[0] + min(z), max_dsc, marker=".", color="green")
            ax1.text(loc_max[0] + min(z), max_dsc - 0.03 * max_dsc, round(max_dsc, 4))

        if len(loc_max) != 1 and len(loc_max) != 0:
            ax1.plot(loc_max[0] + min(z), max_dsc, marker=".", color="green")
            ax1.text(loc_max[0] + min(z), max_dsc - 0.03 * max_dsc, round(max_dsc, 4))
            for i in range(1, len(loc_max)):
                ax1.plot(loc_max[i] + min(z), max_dsc, marker=".", color="green")

        if system_type == "PET/CT":
            plt.savefig(os.path.join(output_date_analysis_path, f'Coregistration_DSC_{dir.replace("-", "")}_PETCT.png'))
        elif "PET/MR" in system_type:
            plt.savefig(os.path.join(output_date_analysis_path, f'Coregistration_DSC_{dir.replace("-", "")}_PETMR{mr_sequence}.png'))

    print('Co-registration analysis completed.')

    window.destroy()


window.title('Co-register analysis')

e1 = tk.Button(window, text="Browse directory", command=main.browse_directory)
lbl1 = tk.Label(window, textvariable=main.directory_path)

tk.Label(window, text='Select analysis directory:').grid(row=0, column=0)
e1.grid(row=0, column=1)
lbl1.grid(row=0, column=2)

choices1 = ['PET/CT', 'PET/MR']
w1 = ttk.Combobox(window, values=choices1)
lbl2 = tk.Label(window, text='Choose the system:')
w1.grid(row=1, column=1)
lbl2.grid(row=1, column=0)

tk.Button(window, text='Run', command=main_code).grid(row=2, column=1, sticky=tk.W, pady=4)

tk.mainloop()
