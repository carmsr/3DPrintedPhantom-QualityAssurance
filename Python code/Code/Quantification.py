# Signal Quantification Analysis for PET/CT and PET/MR systems with the 3D-Printed Phantom.

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


class OptionPet(tk.Toplevel):
    def __init__(self, parent, title):
        super().__init__(parent)
        self.title(title)

        self.result_vol = None
        self.result_act = None
        self.result_time = None
        self.result_radiopharmaceutical = None
        self.result_semidecay_period = None

        self.label1 = ttk.Label(self, text="Initial activity (Bq)")
        self.label1.grid(row=0, column=0, padx=10, pady=10)
        self.entry1 = ttk.Entry(self)
        self.entry1.insert(tk.END, "")
        self.entry1.grid(row=0, column=1, padx=10, pady=10)

        self.label2 = ttk.Label(self, text="Time activity (HHMMSS)")
        self.label2.grid(row=1, column=0, padx=10, pady=10)
        self.entry2 = ttk.Entry(self)
        self.entry2.insert(tk.END, "")
        self.entry2.grid(row=1, column=1, padx=10, pady=10)

        self.label3 = ttk.Label(self, text="Dissolution volume (ml)")
        self.label3.grid(row=2, column=0, padx=10, pady=10)
        self.entry3 = ttk.Entry(self)
        self.entry3.insert(tk.END, "1200")
        self.entry3.grid(row=2, column=1, padx=10, pady=10)

        self.label4 = ttk.Label(self, text="Radiopharmaceutical employed")
        self.label4.grid(row=3, column=0, padx=10, pady=10)
        self.combo = ttk.Combobox(self, values=["18F-FDG", "Other"])
        self.combo.grid(row=3, column=1, padx=10, pady=10)
        self.combo.bind("<<ComboboxSelected>>", self.on_combobox_selected)

        self.button_ok = ttk.Button(self, text="Ok", command=self.ok)
        self.button_ok.grid(row=4, column=1, padx=10, pady=10)

    def on_combobox_selected(self, event):
        selected_option = self.combo.get()
        if selected_option == "Other":
            self.open_semidecay_period_entry()

    def open_semidecay_period_entry(self):
        self.semidecay_period_window = tk.Toplevel(self)
        self.semidecay_period_window.title("Semidecay Period")

        self.label5 = ttk.Label(self.semidecay_period_window, text="Semidecay Period (s)")
        self.label5.grid(row=0, column=0, padx=10, pady=10)
        self.entry4 = ttk.Entry(self.semidecay_period_window)
        self.entry4.grid(row=0, column=1, padx=10, pady=10)

        self.button_ok = ttk.Button(self.semidecay_period_window, text="Ok", command=self.save_semidecay_period)
        self.button_ok.grid(row=1, column=1, padx=10, pady=10)

    def save_semidecay_period(self):
        self.result_semidecay_period = self.entry4.get()
        self.semidecay_period_window.destroy()

    def ok(self):
        self.result_act = self.entry1.get()
        self.result_time = self.entry2.get()
        self.result_vol = self.entry3.get()
        self.result_radiopharmaceutical = self.combo.get()

        if self.result_radiopharmaceutical == "Other":
            self.open_semidecay_period_entry()

        self.destroy()


class OptionCT(tk.Toplevel):
    def __init__(self, parent, title):
        super().__init__(parent)
        self.title(title)

        self.geometry("450x100")  # Ancho x Alto

        self.result_inserts_type = None
        self.result_inserts = None
        self.result_modification = "No"

        self.insert_types = {
            "Calibration Curve": [(0, 'D10'), (1, 'D28a'), (2, 'D28b'), (3, 'D28c'), (4, 'D28d'), (5, 'D20'),
                                  (6, 'D12'), (7, 'Trabecular bone'), (8, 'Dense bone'), (9, 'Muscle'), (10, 'Lung (exhale)')],
            "Materials Characterization": [(0, 'D10'), (1, 'D28a'), (2, 'D28b'), (3, 'D28c'), (4, 'D28d'), (5, 'D20'),
                                           (6, 'D12'), (7, 'DAP0102'), (8, 'SolidWater'), (9, 'HW04'), (10, 'DAP0203')]}

        self.label_source = ttk.Label(self, text="Select Inserts Type")
        self.label_source.grid(row=0, column=0, padx=10, pady=10)
        self.combo_source = ttk.Combobox(self, values=list(self.insert_types.keys()))
        self.combo_source.grid(row=0, column=1, padx=10, pady=10)

        self.button_select = ttk.Button(self, text="Select", command=self.select_inserts_type)
        self.button_select.grid(row=1, column=0, columnspan=2, padx=10, pady=10)

        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)

    def select_inserts_type(self):
        self.result_inserts_type = self.combo_source.get()
        if self.result_inserts_type in self.insert_types:
            self.default_inserts = self.insert_types[self.result_inserts_type]
            self.show_inserts_window()

    def show_inserts_window(self):
        self.inserts_window = tk.Toplevel(self)
        self.inserts_window.title("Position-Insert")
        self.inserts_window.geometry("400x550")

        self.inserts_window.grid_columnconfigure(0, weight=1)
        self.inserts_window.grid_columnconfigure(1, weight=1)

        self.label_inserts_info = ttk.Label(self.inserts_window, text="Position-Insert:")
        self.label_inserts_info.grid(row=0, column=0, columnspan=2, padx=10, pady=10)

        self.insert_labels = []
        for idx, (pos, name) in enumerate(self.default_inserts):
            label = ttk.Label(self.inserts_window, text=f"{pos}: {name}")
            label.grid(row=idx + 1, column=0, columnspan=2, padx=10, pady=5)
            self.insert_labels.append(label)

        self.button_modify = ttk.Button(self.inserts_window, text="Modify", command=self.modify_inserts)
        self.button_modify.grid(row=len(self.default_inserts) + 1, column=0, padx=10, pady=10)

        self.button_ok = ttk.Button(self.inserts_window, text="Ok", command=self.ok)
        self.button_ok.grid(row=len(self.default_inserts) + 1, column=1, padx=10, pady=10)

    def modify_inserts(self):
        self.modify_window = tk.Toplevel(self.inserts_window)
        self.modify_window.title("Modify Position-Insert")
        self.modify_window.geometry("500x120")

        self.modify_window.grid_columnconfigure(0, weight=1)
        self.modify_window.grid_columnconfigure(1, weight=1)
        self.modify_window.grid_columnconfigure(2, weight=1)
        self.modify_window.grid_columnconfigure(3, weight=1)

        ttk.Label(self.modify_window, text="Enter new Position-Insert:").grid(row=0, column=0, columnspan=4, padx=10, pady=10)

        self.new_inserts_entries = []
        self.add_entry()  # Add the first entry row

        self.save_button = ttk.Button(self.modify_window, text="Save", command=self.save_new_inserts)
        self.save_button.grid(row=2, column=0, columnspan=4, padx=10, pady=10)

    def add_entry(self):
        row = len(self.new_inserts_entries) + 1

        new_pos_entry = tk.Entry(self.modify_window)
        new_name_entry = tk.Entry(self.modify_window)
        self.new_inserts_entries.append((new_pos_entry, new_name_entry))

        ttk.Label(self.modify_window, text="Position:").grid(row=row, column=0, padx=5, pady=5, sticky=tk.W)
        new_pos_entry.grid(row=row, column=1, padx=5, pady=5)

        ttk.Label(self.modify_window, text="Name:").grid(row=row, column=2, padx=5, pady=5, sticky=tk.W)
        new_name_entry.grid(row=row, column=3, padx=5, pady=5)

    def save_new_inserts(self):
        new_inserts = [(int(pos_entry.get()), name_entry.get()) for pos_entry, name_entry in self.new_inserts_entries if pos_entry.get().strip() and name_entry.get().strip()]
        for new_pos, new_name in new_inserts:
            for i, (pos, name) in enumerate(self.default_inserts):
                if pos == new_pos:
                    self.default_inserts[i] = (new_pos, new_name)
                    break
            else:
                self.default_inserts.append((new_pos, new_name))

        self.result_modification = "Yes"
        self.result_inserts = self.default_inserts

        # Update the inserts window with new inserts
        for label in self.insert_labels:
            label.destroy()
        self.insert_labels.clear()

        for idx, (pos, name) in enumerate(self.default_inserts):
            label = ttk.Label(self.inserts_window, text=f"{pos}: {name}")
            label.grid(row=idx + 1, column=0, columnspan=2, padx=10, pady=5)
            self.insert_labels.append(label)

        self.modify_window.destroy()

    def ok(self):
        self.result_inserts = self.default_inserts
        self.inserts_window.destroy()
        self.destroy()


main = Path()


# Main code
def main_code():

    import os
    import sys
    import pandas as pd
    import numpy as np
    import math
    import SimpleITK as sitk
    import pydicom
    import Auxiliar_functions
    import matplotlib.pyplot as plt
    from scipy.stats import linregress

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
        output_date_analysis_path = os.path.join(output_date_path, 'Quantification')
        try:
            os.makedirs(output_date_analysis_path)
        except FileExistsError:
            pass

        intensity_ct = []  # list with the mean CT intensities of one image
        err_intensity_ct = []  # list with the mean intensity CT errors of one image
        intensity_pet = []  # list with the mean PET intensities of one image
        err_intensity_pet = []  # list with the mean intensity PET errors of one image
        rc_ac = []  # list with the concentration RC of one image
        intensity_mr_total = []  # list with the mean intensities of various MR images
        err_intensity_mr_total = []  # list with the mean intensity errors of various MR images
        diameter_pet = []  # list with inserts diameter for PET image of one image
        diameter_ct = []  # list with inserts diameter for CT image of one image
        diameter_mr_total = []  # list with inserts diameter of various MR images
        type_mr_total = []  # list with images names of various MR images
        names_mr_total = []  # list with inserts names of various MR images

        # Specify PET parameters and inserts information
        ctmr_executed = False

        # Choose initial slice of PET, CT analysis area, and specify PET parameters and inserts information
        for sub_folder in sub_folders:
            if sub_folder == 'PET':
                # Choose initial slice
                coord_pet = Auxiliar_functions.select_slices(os.path.join(path_date, sub_folder), sub_folder, 'Choose initial slice for quantification analysis', dir)[0]

                # Choose PET acquisition information
                dialog_pet = OptionPet(window, f'Parameters for determining the activity concentration of the solution used for PET {dir}')
                window.wait_window(dialog_pet)
                volume = dialog_pet.result_vol
                activity = dialog_pet.result_act
                time_activity = dialog_pet.result_time
                radioph = dialog_pet.result_radiopharmaceutical
                if radioph == 'Other':
                    other_semidecay_period = dialog_pet.result_semidecay_period
                    if other_semidecay_period:
                        other_semidecay_period = float(other_semidecay_period)
                    elif not other_semidecay_period:
                        print('The semi-decay period has not been entered, it is not possible to continue.')
                        sys.exit()
                if activity:
                    activity = float(activity)
                elif not activity:
                    print('The activity has not been entered, it is not possible to continue.')
                    sys.exit()
                if not time_activity:
                    print('The activity measurement time has not been entered, it is not possible to continue.')
                    sys.exit()
                if volume:
                    volume = float(volume)
                elif not volume:
                    print('The volume has not been entered, it is not possible to continue.')
                    sys.exit()

            elif sub_folder == 'CT':
                # Choose initial slice
                coord_ct = Auxiliar_functions.select_slices(os.path.join(path_date, sub_folder), sub_folder, 'Choose initial slice for quantification analysis', dir)[0]

                # Choose inserts disposition
                if not ctmr_executed:
                    dialog_ctmr = OptionCT(window, f'Disposition of quantification inserts dated {dir}')
                    window.wait_window(dialog_ctmr)
                    inserts_type = dialog_ctmr.result_inserts_type
                    modification = dialog_ctmr.result_modification
                    if modification == 'Yes':
                        new_inserts_disposition = dialog_ctmr.result_inserts
                    ctmr_executed = True

            elif sub_folder == 'MR':
                # Choose initial slice
                mr_images_path = [os.path.join(os.path.join(path_date, sub_folder), d) for d in os.listdir(os.path.join(path_date, sub_folder))]
                for mr_image_path in mr_images_path:
                    rm_file_names = os.listdir(mr_image_path)
                    first_file_path = os.path.join(mr_image_path, rm_file_names[0])  # first dcm path
                    ds = pydicom.filereader.dcmread(first_file_path)
                    type_mr = ds["SeriesDescription"].value
                    if "T1" in type_mr:
                        coord_mrt1 = Auxiliar_functions.select_slices(mr_image_path, sub_folder + "_" + type_mr, 'Choose initial slice for quantification analysis', dir)[0]
                    elif "T2" in type_mr:
                        coord_mrt2 = Auxiliar_functions.select_slices(mr_image_path, sub_folder + "_" + type_mr, 'Choose initial slice for quantification analysis', dir)[0]

                # Choose inserts disposition
                if not ctmr_executed:
                    dialog_ctmr = OptionCT(window, f'Disposition of quantification inserts dated {dir}')
                    window.wait_window(dialog_ctmr)
                    inserts_type = dialog_ctmr.result_inserts_type
                    modification = dialog_ctmr.result_modification
                    if modification == 'Yes':
                        new_inserts_disposition = dialog_ctmr.result_inserts
                    ctmr_executed = True

        # Goes through the image folders contained in the date folder
        for sub_folder in sub_folders:
            if sub_folder == 'CT' or sub_folder == 'MR' or sub_folder == 'PET':
                path_image_type = os.path.join(path_date, sub_folder)

                if sub_folder == 'PET':
                    print(f'Quantifying {sub_folder} image dated {dir}...')

                    pet_image_path = path_image_type  # DICOM serie folder path
                    pet_file_names = os.listdir(pet_image_path)

                    # Extract image metadata: obtain Acquisition Time
                    min_acquisition_time = None
                    min_acquisition_time_file = None
                    for file_name in pet_file_names:
                        file_path = os.path.join(pet_image_path, file_name)
                        ds = pydicom.filereader.dcmread(file_path)
                        time_pet = ds["AcquisitionTime"].value
                        acquisition_time_pet = time_pet.split(".")[0]

                        if min_acquisition_time is None or acquisition_time_pet < min_acquisition_time:
                            min_acquisition_time = acquisition_time_pet
                            min_acquisition_time_file = file_path

                    ds = pydicom.filereader.dcmread(min_acquisition_time_file)
                    time_pet = ds["AcquisitionTime"].value
                    acquisition_time_pet = time_pet.split(".")[0]

                    # Convert to NIFTI format
                    reader_pet = sitk.ImageSeriesReader()
                    dicom_names_pet = reader_pet.GetGDCMSeriesFileNames(pet_image_path)
                    reader_pet.SetFileNames(dicom_names_pet)
                    serie_pet = reader_pet.Execute()
                    nifti_path = os.path.join(output_date_analysis_path, sub_folder + '_image.nii')
                    sitk.WriteImage(serie_pet, nifti_path)

                    im = sitk.ReadImage(nifti_path)
                    image = sitk.GetArrayFromImage(im)  # PET image as a numpy array

                    # Create the interest analysis area
                    area_pet, final_coord_pet = Auxiliar_functions.segmentation_region(im, coord_pet, float(serie_pet.GetSpacing()[2]), sub_folder, 'RC')

                    # CONCENTRATION analysis #
                    # Radiactive inserts segmentation
                    names_pet, centers_pet, inserts_pet = Auxiliar_functions.inserts_pet(im, area_pet, coord_pet, final_coord_pet, serie_pet.GetSpacing())
                    theor_diameter_pet = ['D10', 'D28a', 'D28b', 'D28c', 'D28d', 'D20', 'D12']

                    # Expected activity concentration calculation at image acquisition time
                    ac_time = Auxiliar_functions.time_converter(acquisition_time_pet)  # image acquisition time
                    in_time = Auxiliar_functions.time_converter(time_activity)  # activity measure time
                    diff_time = abs(ac_time - in_time)
                    if radioph == '18F-FDG':
                        activity_decay = (activity / volume) * math.exp(-diff_time * ((math.log(2)) / 6586.200165))
                    elif radioph == 'Other':
                        activity_decay = (activity / volume) * math.exp(-diff_time * ((math.log(2)) / other_semidecay_period))

                    # Go through each insert
                    for i in range(len(inserts_pet)):
                        # Create a reduced analysis region inside the insert for activity quantification
                        segm_activity, diameter_insert = Auxiliar_functions.reduced_segmentation_cyl(inserts_pet[i], float(serie_pet.GetSpacing()[0]))

                        if not np.isnan(diameter_insert):
                            diameter_pet.append(diameter_insert)

                            # ACTIVITY quantification
                            in_pet, err_pet = Auxiliar_functions.quantification(image, sitk.GetArrayFromImage(segm_activity))
                            intensity_pet.append(in_pet)
                            err_intensity_pet.append(err_pet)

                            rc = in_pet / activity_decay  # recovery coefficient
                            rc_ac.append(rc)

                    # Plot the recovery coefficient (RC_activity) of each active insert
                    data_points = list(zip(theor_diameter_pet, rc_ac))
                    data_points_sorted = sorted(data_points, key=lambda point: point[1])
                    theor_diameter_pet_sorted, rc_ac_sorted = zip(*data_points_sorted)
                    fig_rc = plt.figure(num='RC')
                    plt.plot(theor_diameter_pet_sorted, rc_ac_sorted, 'o', label='Activity concentration RC')
                    plt.xlabel('Insert diameter (mm)')
                    plt.ylabel('Recovery Coefficient')
                    plt.ylim([0.3, 1.3])
                    plt.title('Quantification PET')
                    plt.legend()
                    plt.savefig(os.path.join(output_date_analysis_path, f'QuantificationPET_RC.png'))
                    plt.close(fig_rc)

                    # Save Excel results file
                    if len(names_pet) != 0 and len(intensity_pet) != 0:
                        dic_pet = {'Position': names_pet}
                        dic_pet['Insert'] = theor_diameter_pet
                        dic_pet['Diameter (mm)'] = diameter_pet
                        dic_pet['Mean (Bq/ml)'] = intensity_pet
                        dic_pet['Standard Deviation (Bq/ml)'] = err_intensity_pet
                        dic_pet['RC Concentration'] = rc_ac
                        df_pet = pd.DataFrame(data=dic_pet)
                        df_pet = df_pet.sort_values(by='Diameter (mm)')
                        writer = pd.ExcelWriter(os.path.join(output_date_analysis_path, 'QuantificationPET.xlsx'))
                        df_pet.to_excel(writer, sheet_name='Quantification', index=False)
                        dic_pet.clear()
                        writer.save()

                elif sub_folder == 'CT':
                    print(f'Quantifying {sub_folder} image dated {dir}...')

                    ct_image_path = path_image_type  # DICOM serie folder path

                    # Extract image metadata and convert to NIFTI format
                    reader = sitk.ImageSeriesReader()
                    dicom_names = reader.GetGDCMSeriesFileNames(ct_image_path)
                    reader.SetFileNames(dicom_names)
                    serie = reader.Execute()
                    nifti_path = os.path.join(output_date_analysis_path, sub_folder + '_image.nii')
                    sitk.WriteImage(serie, nifti_path)

                    im = sitk.ReadImage(nifti_path)
                    image = sitk.GetArrayFromImage(im)  # CT image as a numpy array

                    # Create the interest analysis area
                    area_ct, final_coord = Auxiliar_functions.segmentation_region(im, coord_ct, float(serie.GetSpacing()[2]), sub_folder, inserts_type)

                    # HU analysis #
                    # Select inserts
                    names_ct, centers_ct, inserts_ct = Auxiliar_functions.inserts_ct(im, coord_ct, final_coord, serie.GetSpacing(), inserts_type)

                    # Go through each insert
                    for i in range(len(inserts_ct)):
                        # Create a reduced analysis region inside the insert for HU quantification
                        segm_den, diameter_in_ct = Auxiliar_functions.reduced_segmentation_cyl(inserts_ct[i], float(serie.GetSpacing()[0]))
                        diameter_ct.append(diameter_in_ct)

                        intensity, err = Auxiliar_functions.quantification(image, sitk.GetArrayFromImage(segm_den))
                        intensity_ct.append(intensity)
                        err_intensity_ct.append(err)

                    if inserts_type == 'Calibration Curve':
                        if modification == 'No':
                            theor_diameter_ct = ['D10', 'D28a', 'D28b', 'D28c', 'D28d', 'D20', 'D12', 'Trabecular bone', 'Dense bone', 'Muscle', 'Lung (exhale)']
                        elif modification == 'Yes':
                            theor_diameter_ct = list()
                            for pos, ins in new_inserts_disposition:
                                theor_diameter_ct.append(ins)
                    elif inserts_type == 'Materials Characterization':
                        if modification == 'No':
                            theor_diameter_ct = ['D10', 'D28a', 'D28b', 'D28c', 'D28d', 'D20', 'D12', 'DAP0102', 'SolidWater', 'HW04', 'DAP0203']
                        elif modification == 'Yes':
                            theor_diameter_ct = list()
                            for pos, ins in new_inserts_disposition:
                                theor_diameter_ct.append(ins)

                    # Calibration curve calculation: HU vs density
                    if inserts_type == 'Calibration Curve':
                        insert_density = {'Muscle': 1.06, 'Lung (inhale)': 0.20, 'Liver': 1.07, 'Breast': 0.99, 'Water': 1.004, 'Adipose': 0.97, 'Trabecular bone': 1.16, 'Lung (exhale)': 0.50, 'Dense bone': 1.61}
                        insert_position = []
                        density_list = []
                        density_list_all = []
                        for insert_name in theor_diameter_ct:
                            if insert_name in insert_density:
                                insert_position.append(theor_diameter_ct.index(insert_name))
                                density_list.append(insert_density[insert_name])
                                density_list_all.append(insert_density[insert_name])
                            else:
                                density_list_all.append(np.NaN)

                        intensity_ct_list = np.array([intensity_ct[indice] for indice in insert_position])

                        regression = np.polyfit(density_list, intensity_ct_list, 1)
                        x = np.linspace(np.min(density_list), np.max(density_list), 2000)
                        y = regression[0] * x + regression[1]

                        slope, intercept, r_value, _, _ = linregress(density_list, intensity_ct_list)
                        R2 = r_value ** 2

                        plt.figure()
                        plt.plot(density_list, intensity_ct_list, 'bo')
                        plt.plot(x, y, 'g--', alpha=0.5, label='Linear regression')
                        if intercept > 0:
                            plt.title(f'Linear regression: {slope:.0f}x + {intercept:.0f}. R square value = {R2:.4f}')
                        else:
                            plt.title(f'Linear regression: {slope:.0f}x - {-intercept:.0f}. R square value = {R2:.4f}')
                        plt.xlabel('Density(g/cm3)')
                        plt.ylabel('Hounsfield Units')
                        plt.savefig(os.path.join(output_date_analysis_path, f'QuantificationCT_CalibrationCurve.png'), dpi=600, format='png', bbox_inches='tight')
                        plt.close()

                    # Save Excel results file
                    if len(names_ct) != 0 and len(intensity_ct) != 0:
                        dic_ct = {'Position': names_ct}
                        if inserts_type == 'Calibration Curve':
                            dic_ct['Density (g/cm3)'] = density_list_all
                        dic_ct['Insert'] = theor_diameter_ct
                        dic_ct['Mean (HU)'] = intensity_ct
                        dic_ct['Standard Deviation (HU)'] = err_intensity_ct

                        df_ct = pd.DataFrame(data=dic_ct)
                        if inserts_type == 'Calibration Curve':
                            df_ct = df_ct.sort_values(by='Density (g/cm3)')
                        else:
                            df_ct = df_ct.sort_values(by='Mean (HU)')
                        writer = pd.ExcelWriter(os.path.join(output_date_analysis_path, 'QuantificationCT.xlsx'))
                        df_ct.to_excel(writer, sheet_name='Quantification', index=False)
                        dic_ct.clear()
                        writer.save()

                elif sub_folder == 'MR':
                    mr_path = path_image_type  # MR images' folder path
                    mr_images_path = [os.path.join(mr_path, d) for d in os.listdir(mr_path)]  # DICOM series folder paths
                    for mr_image_path in mr_images_path:
                        intensity_mr = []
                        err_intensity_mr = []
                        diameter_mr = []

                        # Extract image metadata and convert to NIFTI format
                        rm_file_names = os.listdir(mr_image_path)
                        first_file_path = os.path.join(mr_image_path, rm_file_names[0])  # first dcm path
                        ds = pydicom.filereader.dcmread(first_file_path)
                        type_mr = ds["SeriesDescription"].value
                        name_mr = str(sub_folder + type_mr)
                        print(f'Quantifying {name_mr} image dated {dir}...')

                        reader_mr = sitk.ImageSeriesReader()
                        dicom_names_mr = reader_mr.GetGDCMSeriesFileNames(mr_image_path)
                        reader_mr.SetFileNames(dicom_names_mr)
                        serie_mr = reader_mr.Execute()
                        nifti_path = os.path.join(output_date_analysis_path, sub_folder + str(type_mr) + '_image.nii')
                        sitk.WriteImage(serie_mr, nifti_path)

                        im = sitk.ReadImage(nifti_path)
                        image = sitk.GetArrayFromImage(im)

                        if 'T1' in type_mr:
                            coord_mr = coord_mrt1
                            type_mr_total.append('MRT1')
                        elif 'T2' in type_mr:
                            coord_mr = coord_mrt2
                            type_mr_total.append('MRT2')

                        # Create the interest analysis area
                        area_mr, final_coord = Auxiliar_functions.segmentation_region(im, coord_mr, float(serie_mr.GetSpacing()[2]), sub_folder, inserts_type)

                        # MR analysis #
                        # Inserts segmentation
                        names_mr, centers_mr, inserts_mr = Auxiliar_functions.inserts_mr(im, coord_mr, final_coord, serie_mr.GetSpacing())

                        # Go through each insert
                        for i in range(len(inserts_mr)):
                            # Create a reduced analysis region inside the insert for intensity quantification
                            segm_mr, diameter_in_mr = Auxiliar_functions.reduced_segmentation_cyl(inserts_mr[i], float(serie_mr.GetSpacing()[0]))
                            diameter_mr.append(diameter_in_mr)

                            intensity, err = Auxiliar_functions.quantification(image, sitk.GetArrayFromImage(segm_mr))
                            intensity_mr.append(intensity)
                            err_intensity_mr.append(err)

                        intensity_mr_total.append(intensity_mr)
                        err_intensity_mr_total.append(err_intensity_mr)
                        diameter_mr_total.append(diameter_mr)
                        names_mr_total.append(names_mr)

                        if inserts_type == 'Calibration Curve':
                            if modification == 'No':
                                theor_diameter_mr = ['D10', 'D28a', 'D28b', 'D28c', 'D28d', 'D20', 'D12', 'Trabecular bone', 'Dense bone', 'Muscle', 'Lung (exhale)']
                            elif modification == 'Yes':
                                theor_diameter_mr = list()
                                for pos, ins in new_inserts_disposition:
                                    theor_diameter_mr.append(ins)
                        elif inserts_type == 'Materials Characterization':
                            if modification == 'No':
                                theor_diameter_mr = ['D10', 'D28a', 'D28b', 'D28c', 'D28d', 'D20', 'D12', 'DAP0102', 'SolidWater', 'HW04', 'DAP0203']
                            elif modification == 'Yes':
                                theor_diameter_mr = list()
                                for pos, ins in new_inserts_disposition:
                                    theor_diameter_mr.append(ins)

                    # Save Excel results file
                    for n in range(0, len(names_mr_total)):
                        if len(names_mr_total[n]) != 0 and len(intensity_mr_total[n]) != 0:
                            dic_mr = {'Position': names_mr_total[n]}
                            dic_mr['Insert'] = theor_diameter_mr
                            dic_mr['Mean (ms)'] = intensity_mr_total[n]
                            dic_mr['Standard Deviation (ms)'] = err_intensity_mr_total[n]

                            df_mr = pd.DataFrame(data=dic_mr)
                            df_mr = df_mr.sort_values(by='Mean (ms)')
                            writer = pd.ExcelWriter(os.path.join(output_date_analysis_path, f'Quantification{type_mr_total[n]}.xlsx'))
                            df_mr.to_excel(writer, sheet_name=f'Quantification', index=False)
                            dic_mr.clear()
                            writer.save()

    print('Quantification analysis completed.')

    window.destroy()


window.title('Quantification analysis')
window.geometry("290x70")

e1 = tk.Button(window, text="Browse directory", command=main.browse_directory)
lbl1 = tk.Label(window, textvariable=main.directory_path)

tk.Label(window, text='Select analysis directory:').grid(row=0, column=1)
e1.grid(row=0, column=2)
lbl1.grid(row=0, column=3)

tk.Button(window, text='Run', command=main_code).grid(row=1, column=2, sticky=tk.W, pady=4)

tk.mainloop()
