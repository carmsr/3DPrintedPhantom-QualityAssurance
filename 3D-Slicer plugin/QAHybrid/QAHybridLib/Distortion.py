# Distortion analysis for CT and MR DICOM images: 3D-Printed Phantom.

import math
import os
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog


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
    import SimpleITK as sitk
    import numpy as np
    import pydicom
    import matplotlib.pyplot as plt
    import Auxiliar_functions
    from scipy.ndimage.morphology import binary_erosion
    from scipy.spatial import distance_matrix
    import pandas as pd

    global window, main

    # Folders' paths
    main_path = main.directory  # main folder path
    images_main_path = os.path.join(main_path, "Images")  # images' folder path

    # Create an Output folder
    output_path = os.path.join(main_path, "Output")
    try:
        os.makedirs(output_path)
    except FileExistsError:
        pass

    # Obtains the names of the date folders contained in Images' folder
    sub_dirs = [d for d in os.listdir(images_main_path) if os.path.isdir(os.path.join(images_main_path, d))]

    # Goes through date folders
    for dir in sub_dirs:
        path_date = os.path.join(images_main_path, dir)
        sub_folders = [
            d for d in os.listdir(path_date) if os.path.isdir(os.path.join(path_date, d))
        ]  # names of the images folders

        # Create an Output folder for each date
        output_date_path = os.path.join(output_path, dir)
        try:
            os.makedirs(output_date_path)
        except FileExistsError:
            pass

        # Create an Output folder for the type of analysis
        output_date_analysis_path = os.path.join(output_date_path, "Distortion")
        try:
            os.makedirs(output_date_analysis_path)
        except FileExistsError:
            pass

        # Check whether the cylindrical distortion insert has been filled
        choices = ["yes", "no"]
        dialog = OptionDialog(
            window,
            "Insert filled",
            f"Has the cylindrical distortion insert been filled for the measurement dated {dir}?",
            choices,
        )
        window.wait_window(dialog)
        filled_cyl = dialog.result

        # Choose initial and final slices of cylindrical insert mesh and the panal insert slice
        for sub_folder in sub_folders:
            if sub_folder == "CT":
                in_sl, fin_sl = Auxiliar_functions.select_slices(
                    os.path.join(path_date, sub_folder),
                    sub_folder,
                    "Choose the initial and final slices of cylindrical insert",
                    dir,
                )
                panal_slice = Auxiliar_functions.select_slices(
                    os.path.join(path_date, sub_folder), sub_folder, "Choose the distortion panal slice", dir
                )[0]
            elif sub_folder == "MR":
                mr_images_path = [
                    os.path.join(os.path.join(path_date, sub_folder), d)
                    for d in os.listdir(os.path.join(path_date, sub_folder))
                ]
                for mr_image_path in mr_images_path:
                    rm_file_names = os.listdir(mr_image_path)
                    first_file_path = os.path.join(mr_image_path, rm_file_names[0])  # first dcm path
                    ds = pydicom.filereader.dcmread(first_file_path)
                    type_mr = ds["SeriesDescription"].value
                    if "T1" in type_mr:
                        in_sl_mr_t1, fin_sl_mr_t1 = Auxiliar_functions.select_slices(
                            mr_image_path,
                            sub_folder + "_" + type_mr,
                            "Choose the initial and final slices of cylindrical insert",
                            dir,
                        )
                        panal_slice_mr_t1 = Auxiliar_functions.select_slices(
                            mr_image_path, sub_folder + "_" + type_mr, "Choose the distortion panal slice", dir
                        )[0]
                        # Choose panal insert points
                        panal_points_t1 = []
                        directions = [
                            "Right",
                            "Left",
                            "Diag. up-right",
                            "Diag. down-left",
                            "Diag. up-left",
                            "Diag. down-right",
                        ]
                        for i in range(6):
                            coordinat = Auxiliar_functions.choose_points(
                                mr_image_path, panal_slice_mr_t1, sub_folder + "_" + type_mr, directions[i]
                            )
                            panal_points_t1.append(coordinat)
                    elif "T2" in type_mr:
                        in_sl_mr_t2, fin_sl_mr_t2 = Auxiliar_functions.select_slices(
                            mr_image_path,
                            sub_folder + "_" + type_mr,
                            "Choose the initial and final slices of cylindrical insert",
                            dir,
                        )
                        panal_slice_mr_t2 = Auxiliar_functions.select_slices(
                            mr_image_path, sub_folder + "_" + type_mr, "Choose the distortion panal slice", dir
                        )[0]
                        # Choose panal insert points
                        panal_points_t2 = []
                        directions = [
                            "Right",
                            "Left",
                            "Diag. up-right",
                            "Diag. down-left",
                            "Diag. up-left",
                            "Diag. down-right",
                        ]
                        for i in range(6):
                            coordinat = Auxiliar_functions.choose_points(
                                mr_image_path, panal_slice_mr_t2, sub_folder + "_" + type_mr, directions[i]
                            )
                            panal_points_t2.append(coordinat)

        # Goes through the image folders contained in the date folder
        for sub_folder in sub_folders:
            path_image_type = os.path.join(path_date, sub_folder)

            if sub_folder == "CT":
                print(f"Distortion analysis of {sub_folder} image dated {dir}...")

                ct_dicom_path = path_image_type  # DICOM serie folder path

                # Read DICOM CT serie
                reader = sitk.ImageSeriesReader()
                dicom_names = reader.GetGDCMSeriesFileNames(ct_dicom_path)
                reader.SetFileNames(dicom_names)
                serie = reader.Execute()

                # Read SliceThickness and PixelSpacing metadata from CT image
                files = [os.path.join(ct_dicom_path, f) for f in os.listdir(ct_dicom_path)]
                if not files:
                    raise ValueError(f"No DICOM files found in folder {ct_dicom_path}")
                ds = pydicom.dcmread(files[0])
                ct_spacingslices = abs(ds.SpacingBetweenSlices)
                ct_spacing = serie.GetSpacing()

                # Convert DICOM CT serie to NIFTI format
                ct_image_path = os.path.join(output_date_analysis_path, sub_folder + "_image.nii")
                sitk.WriteImage(serie, ct_image_path)

                im = sitk.ReadImage(ct_image_path)

                # Obtain phantom center in xy plane
                center_phantom = Auxiliar_functions.phantom_center(im, panal_slice, ct_spacingslices, sub_folder)

                ## CYLINDRICAL INSERT ##
                # Threshold containing cylindrical insert
                thresh_cyl_sitk = Auxiliar_functions.th_cylinsert(
                    im, in_sl, fin_sl, ct_spacing, center_phantom, filled_cyl, sub_folder
                )
                thresh_cyl = sitk.GetArrayFromImage(thresh_cyl_sitk)

                z_planes = []  # list with cyl insert slices (float)
                z_slices = []  # list with cyl insert slices (int)
                # Iterate along dimension z to identify contiguous insert's slices
                for z in range(thresh_cyl.shape[0] - 1):
                    current_slice = thresh_cyl[z, :, :]
                    next_slice = thresh_cyl[z + 1, :, :]
                    # Check for non-zero values in both slices that exceed the threshold
                    if np.count_nonzero(current_slice) > 300 and np.count_nonzero(next_slice) > 300:
                        # Combines both slices in binary form
                        combined_slice = current_slice + next_slice
                        new_slice = np.where(combined_slice >= 1, 1, combined_slice)
                        # Calculate the average of z
                        averaged_z = int((z + (z + 1)) / 2)
                        z_planes.append(((z + (z + 1)) / 2))
                        z_slices.append(averaged_z)
                        # Update the segmentation values
                        thresh_cyl[z, :, :] = 0
                        thresh_cyl[z + 1, :, :] = 0
                        thresh_cyl[averaged_z, :, :] = new_slice
                    elif np.count_nonzero(current_slice) < 300 and np.count_nonzero(next_slice) < 300:
                        thresh_cyl[z, :, :] = 0
                        thresh_cyl[z + 1, :, :] = 0
                    elif np.count_nonzero(current_slice) > 300 and np.count_nonzero(next_slice) < 300:
                        z_planes.append(z)
                        z_slices.append(z)

                c_sitk = sitk.GetImageFromArray(thresh_cyl)
                c_sitk.SetDirection(thresh_cyl_sitk.GetDirection())
                c_sitk.SetSpacing(thresh_cyl_sitk.GetSpacing())
                c_sitk.SetOrigin(thresh_cyl_sitk.GetOrigin())

                # OBTAINS MESH VERTICES (xy coordinates)
                mesh = np.zeros_like(thresh_cyl, dtype=int)
                fused_vertices = []
                for z in z_slices:
                    z_thresh = thresh_cyl[z, :, :]  # image slice
                    # Proceeds if there is any non-null pixel
                    if np.any(z_thresh != 0):
                        # Erosion of the segmentation to remove the wall thickness of the squares
                        eroded_segmentation = binary_erosion(z_thresh, iterations=3)
                        # plt.imshow(eroded_segmentation, cmap='gray')
                        # plt.show()

                        # Get coordinates of pixels that are part of the eroded segmentation
                        y_coords, x_coords = np.nonzero(eroded_segmentation)

                        if len(y_coords) == 0 or len(x_coords) == 0:
                            continue
                        else:
                            # Get the coordinates of the vertices of each square of the mesh
                            vertices = [(z, y, x) for y, x in zip(y_coords, x_coords)]

                            # Calculate distance matrix between vertices
                            dist_matrix = distance_matrix(vertices, vertices)

                            n_vertices = len(vertices)
                            visited_vertices = [False] * n_vertices

                            for i in range(n_vertices):
                                if visited_vertices[i]:
                                    continue

                                # Identify the vertices near this vertex
                                close_vertices = np.where(dist_matrix[i] < 3)[0]

                                # Take the average value of the coordinates of nearby vertices
                                fused_vertex = np.mean(np.array([vertices[j] for j in close_vertices]), axis=0)

                                # Check if the merged coordinate is already present in fused_vertices or if any of the contiguous coordinates are present
                                contiguous_vertices = [
                                    vertex for vertex in fused_vertices if np.linalg.norm(vertex - fused_vertex) < 2
                                ]
                                if not contiguous_vertices:
                                    fused_vertices.append(fused_vertex)  # Add merged coordinate to the list
                                else:
                                    averaged_vertex = np.mean(
                                        np.concatenate((contiguous_vertices, np.array([fused_vertex]))), axis=0
                                    )
                                    fused_vertices = [
                                        v
                                        for v in fused_vertices
                                        if not any(np.array_equal(v, vertex) for vertex in contiguous_vertices)
                                    ]  # Delete contiguous_vertices from fused_vertices
                                    fused_vertices.append(averaged_vertex)  # Add the average to the list fused_vertices

                                # Mark visited vertices as merged
                                close_vertices = np.array(close_vertices).astype(int)
                                visited_vertices = np.array(visited_vertices).astype(int)
                                visited_vertices[close_vertices] = True

                # Convert the merged list of coordinates into a NumPy array and generate the point grid
                vertices_array_fused = np.array(fused_vertices).astype(int)
                mesh[vertices_array_fused[:, 0], vertices_array_fused[:, 1], vertices_array_fused[:, 2]] = 1
                mesh_sitk = sitk.GetImageFromArray(mesh)
                mesh_sitk.SetDirection(thresh_cyl_sitk.GetDirection())
                mesh_sitk.SetSpacing(thresh_cyl_sitk.GetSpacing())
                mesh_sitk.SetOrigin(thresh_cyl_sitk.GetOrigin())
                # sitk.WriteImage(mesh_sitk, "C:\\Users\\20970341E\\Documents\\Mesh_CT.nrrd")

                # AXIAL DISTANCE: Analyse the distance between the planes (z_planes) of the mesh in the axial axis
                distance_between_planes = np.diff(z_planes)
                axial_distance = distance_between_planes * ct_spacingslices
                mean_axial_distance = np.mean(axial_distance)
                max_axial_distance = np.max(axial_distance)
                min_axial_distance = np.min(axial_distance)
                std_value = np.std(axial_distance)

                # Plots axial distances
                plt.figure()
                plt.scatter(range(len(axial_distance)), axial_distance, marker="o")
                plt.axhline(y=mean_axial_distance, color="green", label="Mean distance")
                plt.ylim([min_axial_distance - 2, max_axial_distance + 3])
                plt.legend(
                    [
                        f"Distances between axial planes",
                        f"Mean distance in axial axis: {round(mean_axial_distance, 2)} \u00B1 {round(std_value, 2)} mm",
                    ]
                )
                plt.ylabel("Measured distance (mm)")
                plt.title(
                    f"Distances along superior-inferior direction of cylindrical insert \n Slice spacing: {round(ct_spacingslices, 2)} mm"
                )
                plt.xticks([])
                plt.savefig(os.path.join(output_date_analysis_path, f"DistortionCT_AxialDistances.png"))
                plt.close()

                # MESH SIDES: Analyse the sides of the mesh squares in each slice.
                mesh_point = sitk.GetArrayFromImage(mesh_sitk)

                right_sides = []
                left_sides = []
                up_sides = []
                down_sides = []
                z_slices_filtered = []
                mesh_centers = []
                centers_differences = []
                for i in z_slices:
                    plane = mesh_point[i, :, :]
                    y_nonzero, x_nonzero = np.where(plane == 1)
                    if len(y_nonzero) >= 7 and len(y_nonzero) <= 9:
                        z_slices_filtered.append(i)
                        nonzero_coords = np.where(plane == 1)
                        combined_nonzero_coords = list(zip(nonzero_coords[0], nonzero_coords[1]))
                        y_nonzero_mean, x_nonzero_mean = np.mean(y_nonzero), np.mean(x_nonzero)

                        distances = [
                            np.sqrt((coord[0] - y_nonzero_mean) ** 2 + (coord[1] - x_nonzero_mean) ** 2)
                            for coord in combined_nonzero_coords
                        ]
                        y_central, x_central = combined_nonzero_coords[np.argmin(distances)]
                        mesh_centers.append([y_central, x_central])

                        # Calculates deviation of (y,x) from first
                        first_y, first_x = mesh_centers[0]
                        difference_yx = [y_central - first_y, x_central - first_x]
                        centers_differences.append(difference_yx)

                        # Choose the neighbours closest to the central point
                        neigh_coords = [
                            coord
                            for coord in combined_nonzero_coords
                            if np.sqrt(
                                ((coord[0] - y_central) * ct_spacing[1]) ** 2
                                + ((coord[1] - x_central) * ct_spacing[0]) ** 2
                            )
                            <= 16
                        ]

                        # Calculate horizontal sides (x-axis) and vertical sides (y-axis):
                        for n_coord in neigh_coords:
                            y_diff = n_coord[0] - y_central
                            x_diff = n_coord[1] - x_central

                            # Horizontal sides
                            if abs(y_diff) <= (3 / ct_spacing[1]) and abs(x_diff) >= (3 / ct_spacing[0]):
                                if x_diff > 0:
                                    right_sides.append(abs(x_diff * ct_spacing[0]))
                                elif x_diff < 0:
                                    left_sides.append(abs(x_diff * ct_spacing[0]))
                            # Vertical sides
                            if abs(x_diff) <= (3 / ct_spacing[0]) and abs(y_diff) >= (3 / ct_spacing[1]):
                                if y_diff < 0:
                                    up_sides.append(abs(y_diff * ct_spacing[1]))
                                elif y_diff > 0:
                                    down_sides.append(abs(y_diff * ct_spacing[1]))

                        all_sides = [right_sides, left_sides, up_sides, down_sides]
                        max_size = max(len(side_list) for side_list in all_sides)
                        if len(right_sides) < max_size:
                            right_sides.extend([None] * (max_size - len(right_sides)))
                        if len(left_sides) < max_size:
                            left_sides.extend([None] * (max_size - len(left_sides)))
                        if len(up_sides) < max_size:
                            up_sides.extend([None] * (max_size - len(up_sides)))
                        if len(down_sides) < max_size:
                            down_sides.extend([None] * (max_size - len(down_sides)))

                # Plots sides values
                sides_fig = plt.figure()
                sides_colors = ["red", "blue", "green", "orange"]
                sides_labels = ["Right", "Left", "Up", "Down"]
                markers = ["o", "s", "D", "v"]
                ax = sides_fig.add_subplot(1, 1, 1)
                # Iterate over each type of side and its respective values
                for i, side_list in enumerate([right_sides, left_sides, up_sides, down_sides]):
                    side_list_filtered = [value for value in side_list if value is not None]
                    mean_value = np.mean(side_list_filtered)
                    mean_std = np.std(side_list_filtered) / np.sqrt(len(side_list_filtered))
                    ax.scatter(
                        z_slices_filtered, side_list, color=sides_colors[i], label=sides_labels[i], marker=markers[i]
                    )
                    plt.axhline(
                        y=mean_value,
                        color=sides_colors[i],
                        label=f"Mean {sides_labels[i]} side: {round(mean_value, 2)} \u00B1 {round(mean_std, 2)} mm",
                    )

                ax.set_ylabel("Measured distance (mm)")
                ax.set_xlabel("Slice")
                ax.set_title(
                    f"Distances along transaxial plane of the cylindrical insert \n Pixel dimension: {round(ct_spacing[0], 2)}x{round(ct_spacing[1], 2)} mm\u00B2"
                )
                ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc="upper left")
                plt.ylim([6, 14])
                sides_fig.savefig(
                    os.path.join(output_date_analysis_path, f"DistortionCT_MeshSides.png"), bbox_inches="tight"
                )
                plt.close(sides_fig)

                # Plots (y,x) mesh central deviation
                central_fig = plt.figure()
                central_colors = ["blue", "green"]
                central_labels = ["y-axis", "x-axis"]
                markers = ["o", "s"]
                ax = central_fig.add_subplot(1, 1, 1)

                ax.scatter(
                    z_slices_filtered,
                    [coord[0] for coord in centers_differences],
                    color=central_colors[0],
                    label=central_labels[0],
                    marker=markers[0],
                )
                ax.scatter(
                    z_slices_filtered,
                    [coord[1] for coord in centers_differences],
                    color=central_colors[1],
                    label=central_labels[1],
                    marker=markers[1],
                )

                max_y = max(coord[0] for coord in centers_differences)
                min_y = min(coord[0] for coord in centers_differences)
                max_x = max(coord[1] for coord in centers_differences)
                min_x = min(coord[1] for coord in centers_differences)
                upper_limit = max(max_y, max_x)
                lower_limit = min(min_y, min_x)
                margin = 10
                upper_limit += margin
                lower_limit -= margin

                ax.set_ylabel("Mesh centres deviation")
                ax.set_xlabel("Slice")
                ax.set_title(
                    f"Mesh centres deviation \n Pixel dimension: {round(ct_spacing[0], 2)}x{round(ct_spacing[1], 2)} mm\u00B2"
                )
                ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc="upper left")
                ax.set_ylim([lower_limit, upper_limit])
                central_fig.savefig(
                    os.path.join(output_date_analysis_path, f"DistortionCT_MeshCentres.png"), bbox_inches="tight"
                )
                plt.close(central_fig)

                # Excel document creation
                writer_sides = pd.ExcelWriter(os.path.join(output_date_analysis_path, "DistortionCT_MeshSides.xlsx"))
                writer_ax = pd.ExcelWriter(os.path.join(output_date_analysis_path, "DistortionCT_AxialDistances.xlsx"))

                # Dictionary and DataFrame creation
                for n in range(0, len(z_slices)):
                    if len(right_sides) != 0 and len(left_sides) != 0 and len(up_sides) != 0 and len(down_sides) != 0:
                        sides = {"Slice": z_slices_filtered}
                        sides["Right side (mm)"] = right_sides
                        sides["Left side (mm)"] = left_sides
                        sides["Up side (mm)"] = up_sides
                        sides["Down side (mm)"] = down_sides
                        sides["Mesh y-axis center deviation"] = [coord[0] for coord in centers_differences]
                        sides["Mesh x-axis center deviation"] = [coord[1] for coord in centers_differences]

                        df = pd.DataFrame(data=sides)
                        df.to_excel(writer_sides, sheet_name="Mesh sides distortion", index=False)

                for p in range(0, len(axial_distance)):
                    if len(axial_distance) != 0:
                        dist = {"Axial distances (mm)": axial_distance}
                        df_dist = pd.DataFrame(data=dist)
                        df_dist.to_excel(writer_ax, sheet_name="Axial distortion", index=False)

                writer_sides.save()
                writer_ax.save()

                ## PANAL INSERT ##
                # Threshold of panal insert
                thresh_panal_sitk = Auxiliar_functions.ct_th_panalinsert(
                    im, panal_slice, center_phantom, serie.GetSpacing()
                )
                thresh_panal = sitk.GetArrayFromImage(thresh_panal_sitk)
                z_thresh_panal = thresh_panal[panal_slice, :, :]  # axial slice with the panal segmentation

                # Distortion analysis
                y_panal, x_panal = z_thresh_panal.shape
                y_ini = round(center_phantom[0])  # panal center y coordinate
                x_ini = round(center_phantom[1])  # panal center x coordinate

                # First, we check distances in the left-right direction
                pixel_right = []
                coord_right = []
                pixel_left = []
                coord_left = []
                # We move to the right
                for x in range(x_ini, x_panal, 1):
                    pixel_right.append(z_thresh_panal[y_ini, x])
                    coord_right.append([y_ini, x])
                # We move to the left
                for x in range(x_ini, -1, -1):
                    pixel_left.append(z_thresh_panal[y_ini, x])
                    coord_left.append([y_ini, x])

                # Finally, we check distances in the two diagonal direction
                # Secondary diagonal
                pixel_ds_up = []
                coord_ds_up = []
                pixel_ds_down = []
                coord_ds_down = []
                prev_x = x_ini
                prev_y = y_ini
                for i in range(0, 125):
                    new_x = round(prev_x + 2 * math.cos(math.radians(60)))
                    new_y = round(prev_y - 2 * math.sin(math.radians(60)))
                    pixel_ds_up.append(z_thresh_panal[new_y, new_x])
                    coord_ds_up.append([new_y, new_x])
                    prev_x = new_x
                    prev_y = new_y

                prev_x = x_ini
                prev_y = y_ini
                for i in range(0, 125):
                    new_x = round(prev_x - 2 * math.cos(math.radians(60)))
                    new_y = round(prev_y + 2 * math.sin(math.radians(60)))
                    pixel_ds_down.append(z_thresh_panal[new_y, new_x])
                    coord_ds_down.append([new_y, new_x])
                    prev_x = new_x
                    prev_y = new_y

                # Main diagonal
                pixel_dm_up = []
                coord_dm_up = []
                pixel_dm_down = []
                coord_dm_down = []
                prev_x = x_ini
                prev_y = y_ini
                for i in range(0, 125):
                    new_x = round(prev_x + 2 * math.cos(math.radians(60)))
                    new_y = round(prev_y + 2 * math.sin(math.radians(60)))
                    pixel_dm_down.append(z_thresh_panal[new_y, new_x])
                    coord_dm_down.append([new_y, new_x])
                    prev_x = new_x
                    prev_y = new_y

                prev_x = x_ini
                prev_y = y_ini
                for i in range(0, 125):
                    new_x = round(prev_x - 2 * math.cos(math.radians(60)))
                    new_y = round(prev_y - 2 * math.sin(math.radians(60)))
                    pixel_dm_up.append(z_thresh_panal[new_y, new_x])
                    coord_dm_up.append([new_y, new_x])
                    prev_x = new_x
                    prev_y = new_y

                # Once we have the lists with the pixel values in different directions, we obtain the distances of the
                # walls from the center of the panal in each direction (in mm)
                # Distances in list position between walls and panal origin
                wall_coord_r = Auxiliar_functions.panal_walls(pixel_right, coord_right)
                distances_coord_r_mm = Auxiliar_functions.calculate_mm_distances(
                    wall_coord_r, ct_spacing[1], ct_spacing[0]
                )

                wall_coord_l = Auxiliar_functions.panal_walls(pixel_left, coord_left)
                distances_coord_l_mm = Auxiliar_functions.calculate_mm_distances(
                    wall_coord_l, ct_spacing[1], ct_spacing[0]
                )

                wall_coord_dsu = Auxiliar_functions.panal_walls(pixel_ds_up, coord_ds_up)
                distances_coord_dsu_mm = Auxiliar_functions.calculate_mm_distances(
                    wall_coord_dsu, ct_spacing[1], ct_spacing[0]
                )

                wall_coord_dsd = Auxiliar_functions.panal_walls(pixel_ds_down, coord_ds_down)
                distances_coord_dsd_mm = Auxiliar_functions.calculate_mm_distances(
                    wall_coord_dsd, ct_spacing[1], ct_spacing[0]
                )

                wall_coord_dmu = Auxiliar_functions.panal_walls(pixel_dm_up, coord_dm_up)
                distances_coord_dmu_mm = Auxiliar_functions.calculate_mm_distances(
                    wall_coord_dmu, ct_spacing[1], ct_spacing[0]
                )

                wall_coord_dmd = Auxiliar_functions.panal_walls(pixel_dm_down, coord_dm_down)
                distances_coord_dmd_mm = Auxiliar_functions.calculate_mm_distances(
                    wall_coord_dmd, ct_spacing[1], ct_spacing[0]
                )

                # Plots walls distances
                walls_fig = plt.figure()
                walls_colors = ["#FF4444", "#55CC55", "#FFCC33", "#3366FF", "#CC6699", "#6666CC"]
                walls_labels = [
                    "Right",
                    "Left",
                    "Diag. up-right",
                    "Diag. down-left",
                    "Diag. up-left",
                    "Diag. down-right",
                ]
                markers = ["o", "s", "D", "v", "h", "p"]
                ax = walls_fig.add_subplot(1, 1, 1)

                # Iterate over each type of side and its respective values
                for i, walls_list in enumerate(
                    [
                        distances_coord_r_mm,
                        distances_coord_l_mm,
                        distances_coord_dsu_mm,
                        distances_coord_dsd_mm,
                        distances_coord_dmu_mm,
                        distances_coord_dmd_mm,
                    ]
                ):
                    walls_list_filtered = [value for value in walls_list if value is not None]
                    mean_value = np.mean(walls_list_filtered)
                    mean_std = np.std(walls_list_filtered) / np.sqrt(len(walls_list_filtered))

                    ax.scatter(
                        range(len(walls_list)),
                        walls_list,
                        color=walls_colors[i],
                        label=walls_labels[i],
                        marker=markers[i],
                    )
                    plt.axhline(
                        y=mean_value,
                        color=walls_colors[i],
                        label=f"Mean {walls_labels[i]} side: {round(mean_value, 2)} \u00B1 {round(mean_std, 2)} mm",
                    )

                ax.set_ylabel("Measured distance (mm)")
                ax.set_title(
                    f"Distances between panal insert walls \n Pixel dimension: {round(ct_spacing[0], 2)}x{round(ct_spacing[1], 2)} mm\u00B2"
                )
                plt.ylim([23, 30])
                ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc="upper left")
                ax.set_xticks([])
                walls_fig.savefig(
                    os.path.join(output_date_analysis_path, f"DistortionCT_PanalWalls.png"), bbox_inches="tight"
                )
                plt.close(walls_fig)

                # Save distances in each direction in an Excel file
                writer_panal = pd.ExcelWriter(os.path.join(output_date_analysis_path, "DistortionCT_PanalWalls.xlsx"))

                # Dictionary and DataFrame creation
                series_r = pd.Series(distances_coord_r_mm, name="Right (mm)")
                series_l = pd.Series(distances_coord_l_mm, name="Left (mm)")
                series_dsu = pd.Series(distances_coord_dsu_mm, name="Diag. up-right (mm)")
                series_dsd = pd.Series(distances_coord_dsd_mm, name="Diag. down-left (mm)")
                series_dmu = pd.Series(distances_coord_dmu_mm, name="Diag. up-left (mm)")
                series_dmd = pd.Series(distances_coord_dmd_mm, name="Diag. down-right (mm)")

                df = pd.concat([series_r, series_l, series_dsu, series_dsd, series_dmu, series_dmd], axis=1)
                df.to_excel(writer_panal, sheet_name="Panal sides distortion", index=False)
                writer_panal.save()

            elif sub_folder == "MR":
                print(f"Distortion analysis of {sub_folder} images dated {dir}...")

                mr_path = path_image_type  # MR images' folder path
                mr_images_path = [os.path.join(mr_path, d) for d in os.listdir(mr_path)]  # DICOM series folder path
                for mr_image_path in mr_images_path:
                    # Extract image metadata and convert to NIFTI format
                    reader_mr = sitk.ImageSeriesReader()
                    dicom_names_mr = reader_mr.GetGDCMSeriesFileNames(mr_image_path)
                    reader_mr.SetFileNames(dicom_names_mr)
                    serie_mr = reader_mr.Execute()

                    rm_file_names = os.listdir(mr_image_path)
                    first_file_path = os.path.join(mr_image_path, rm_file_names[0])  # first dcm path
                    ds = pydicom.filereader.dcmread(first_file_path)
                    type_mr = ds["SeriesDescription"].value

                    # Convert DICOM CT serie to NIFTI format
                    mr_nifti_path = os.path.join(output_date_analysis_path, sub_folder + str(type_mr) + "_image.nii")
                    sitk.WriteImage(serie_mr, mr_nifti_path)

                    im = sitk.ReadImage(mr_nifti_path)
                    spacing = im.GetSpacing()
                    mr_spacingslices = spacing[2]

                    if "T1" in type_mr:
                        panal_slice_mr = panal_slice_mr_t1
                        in_sl_mr = in_sl_mr_t1
                        fin_sl_mr = fin_sl_mr_t1
                        panal_points = panal_points_t1
                        iter = 2
                        mr_img = "MRT1"

                        # N4BiasFieldCorrection
                        # ima = im_filtered
                        ima = im

                    elif "T2" in type_mr:
                        panal_slice_mr = panal_slice_mr_t2
                        in_sl_mr = in_sl_mr_t2
                        fin_sl_mr = fin_sl_mr_t2
                        panal_points = panal_points_t2
                        iter = 3
                        ima = im
                        mr_img = "MRT2"

                    # Obtain phantom center in xy plane
                    center_phantom = Auxiliar_functions.phantom_center(
                        ima, panal_slice_mr, mr_spacingslices, sub_folder
                    )

                    # CYLINDRICAL INSERT #
                    # Threshold containing cylindrical insert
                    thresh_cyl_sitk = Auxiliar_functions.th_cylinsert(
                        ima,
                        in_sl_mr,
                        fin_sl_mr,
                        serie_mr.GetSpacing(),
                        center_phantom,
                        filled_cyl,
                        sub_folder + str(type_mr),
                    )
                    thresh_cyl = sitk.GetArrayFromImage(thresh_cyl_sitk)

                    # sitk.WriteImage(thresh_cyl_sitk, "C:\\Users\\20970341E\\Documents\\CylMRAxT1_121223.nrrd")

                    z_planes = []  # list with cyl insert slices (float)
                    z_slices = []  # list with cyl insert slices (int)
                    # Iterate along dimension z to identify contiguous slices
                    for z in range(thresh_cyl.shape[0] - 1):
                        current_slice = thresh_cyl[z, :, :]
                        next_slice = thresh_cyl[z + 1, :, :]
                        # Check for non-zero values in both slices that exceed the threshold
                        if np.count_nonzero(current_slice) > 50 and np.count_nonzero(next_slice) > 50:
                            # Combines both slices in binary form
                            combined_slice = current_slice + next_slice
                            new_slice = np.where(combined_slice > 1, 1, combined_slice)
                            # Calculate the average of z
                            averaged_z = int((z + (z + 1)) / 2)
                            z_planes.append(((z + (z + 1)) / 2))
                            z_slices.append(averaged_z)
                            # Update the segmentation values
                            thresh_cyl[z, :, :] = 0
                            thresh_cyl[z + 1, :, :] = 0
                            thresh_cyl[averaged_z, :, :] = new_slice
                        elif np.count_nonzero(current_slice) < 50 and np.count_nonzero(next_slice) < 50:
                            thresh_cyl[z, :, :] = 0
                            thresh_cyl[z + 1, :, :] = 0
                        elif np.count_nonzero(current_slice) > 50 and np.count_nonzero(next_slice) < 50:
                            z_planes.append(z)
                            z_slices.append(z)

                    # OBTAINS MESH VERTICES
                    mesh = np.zeros_like(thresh_cyl, dtype=int)
                    fused_vertices = []
                    for z in z_slices:
                        z_thresh = thresh_cyl[z, :, :]  # image slice
                        # Proceeds if there is any non-null pixel
                        if np.any(z_thresh != 0):
                            # Erosion of the segmentation to remove the wall thickness of the squares
                            eroded_segmentation = binary_erosion(z_thresh, iterations=iter)

                            # Get coordinates of pixels that are part of the eroded segmentation
                            y_coords, x_coords = np.nonzero(eroded_segmentation)

                            if len(y_coords) == 0 or len(x_coords) == 0:
                                continue
                            else:
                                # Get the coordinates of the vertices of each square of the mesh
                                vertices = [(z, y, x) for y, x in zip(y_coords, x_coords)]

                                # Calculate distance matrix between vertices
                                dist_matrix = distance_matrix(vertices, vertices)

                                n_vertices = len(vertices)
                                visited_vertices = [False] * n_vertices

                                for i in range(n_vertices):
                                    if visited_vertices[i]:
                                        continue

                                    # Identify the vertices near this vertex
                                    close_vertices = np.where(dist_matrix[i] < 3)[0]

                                    # Take the average value of the coordinates of nearby vertices
                                    fused_vertex = np.mean(np.array([vertices[j] for j in close_vertices]), axis=0)

                                    # Check if the merged coordinate is already present in fused_vertices or if any of the contiguous coordinates are present
                                    contiguous_vertices = [
                                        vertex for vertex in fused_vertices if np.linalg.norm(vertex - fused_vertex) < 2
                                    ]
                                    if not contiguous_vertices:
                                        fused_vertices.append(fused_vertex)  # Add merged coordinate to the list
                                    else:
                                        averaged_vertex = np.mean(
                                            np.concatenate((contiguous_vertices, np.array([fused_vertex]))), axis=0
                                        )
                                        fused_vertices = [
                                            v
                                            for v in fused_vertices
                                            if not any(np.array_equal(v, vertex) for vertex in contiguous_vertices)
                                        ]  # Delete contiguous_vertices from fused_vertices
                                        fused_vertices.append(
                                            averaged_vertex
                                        )  # Add the average to the list fused_vertices

                                    # Mark visited vertices as merged
                                    close_vertices = np.array(close_vertices).astype(int)
                                    visited_vertices = np.array(visited_vertices).astype(int)
                                    visited_vertices[close_vertices] = True

                    # Convert the merged list of coordinates into a NumPy array and generate the point grid
                    vertices_array_fused = np.array(fused_vertices).astype(int)
                    mesh[vertices_array_fused[:, 0], vertices_array_fused[:, 1], vertices_array_fused[:, 2]] = 1
                    mesh_sitk = sitk.GetImageFromArray(mesh)
                    mesh_sitk.SetDirection(thresh_cyl_sitk.GetDirection())
                    mesh_sitk.SetSpacing(thresh_cyl_sitk.GetSpacing())
                    mesh_sitk.SetOrigin(thresh_cyl_sitk.GetOrigin())
                    # sitk.WriteImage(mesh_sitk, "C:\\Users\\20970341E\\Documents\\MeshMRAxT1_121223.nrrd")

                    # AXIAL DISTANCE: Analyse the distance between the planes (z_planes) of the mesh in the axial axis.
                    distance_between_planes = np.diff(z_planes)
                    axial_distance = distance_between_planes * mr_spacingslices
                    mean_axial_distance = np.mean(axial_distance)
                    max_axial_distance = np.max(axial_distance)
                    min_axial_distance = np.min(axial_distance)
                    std_value = np.std(axial_distance)

                    # Plots axial distances
                    plt.figure()
                    plt.scatter(range(len(axial_distance)), axial_distance, marker="o")
                    plt.gca().axes.get_xaxis().set_visible(False)
                    plt.axhline(y=mean_axial_distance, color="green", label="Mean distance")
                    plt.ylim([min_axial_distance - 2, max_axial_distance + 3])
                    plt.legend(
                        [
                            f"Distances between axial planes",
                            f"Mean distance in axial axis: {round(mean_axial_distance, 2)} \u00B1 {round(std_value, 2)} mm",
                        ]
                    )
                    plt.ylabel("Measured distance (mm)")
                    plt.title(
                        f"Distances along superior-inferior direction of cylindrical insert \n Slice spacing: {round(mr_spacingslices, 2)} mm"
                    )
                    plt.xticks([])
                    plt.savefig(os.path.join(output_date_analysis_path, f"Distortion{mr_img}_AxialDistances.png"))
                    plt.close()

                    # MESH SIDES: Analyse the sides of the mesh squares in each slice.
                    mesh_point = sitk.GetArrayFromImage(mesh_sitk)

                    right_sides = []
                    left_sides = []
                    up_sides = []
                    down_sides = []
                    z_slices_filtered = []
                    mesh_centers = []
                    centers_differences = []
                    for i in z_slices:
                        plane = mesh_point[i, :, :]
                        y_nonzero, x_nonzero = np.where(plane == 1)
                        if len(y_nonzero) == 9:
                            z_slices_filtered.append(i)
                            nonzero_coords = np.where(plane == 1)
                            combined_nonzero_coords = list(zip(nonzero_coords[0], nonzero_coords[1]))
                            y_nonzero_mean, x_nonzero_mean = np.mean(y_nonzero), np.mean(x_nonzero)

                            distances = [
                                np.sqrt((coord[0] - y_nonzero_mean) ** 2 + (coord[1] - x_nonzero_mean) ** 2)
                                for coord in combined_nonzero_coords
                            ]
                            y_central, x_central = combined_nonzero_coords[np.argmin(distances)]
                            mesh_centers.append([y_central, x_central])

                            # Calculates deviation of (y,x) from first
                            first_y, first_x = mesh_centers[0]
                            difference_yx = [y_central - first_y, x_central - first_x]
                            centers_differences.append(difference_yx)

                            # Choose the neighbours closest to the central point
                            neigh_coords = [
                                coord
                                for coord in combined_nonzero_coords
                                if np.sqrt(
                                    ((coord[0] - y_central) * spacing[1]) ** 2
                                    + ((coord[1] - x_central) * spacing[0]) ** 2
                                )
                                <= 12
                            ]

                            # Calculate horizontal sides (right, left) and vertical sides (up, down)
                            for n_coord in neigh_coords:
                                y_diff = n_coord[0] - y_central
                                x_diff = n_coord[1] - x_central

                                # Right
                                if x_diff > 0 and abs(x_diff) >= (3 / spacing[0]) and abs(y_diff) <= (3 / spacing[1]):
                                    x_diff = abs(x_diff * spacing[0])
                                    right_sides.append(x_diff)
                                # Left
                                elif x_diff < 0 and abs(x_diff) >= (3 / spacing[0]) and abs(y_diff) <= (3 / spacing[1]):
                                    x_diff = abs(x_diff * spacing[0])
                                    left_sides.append(x_diff)
                                # Up
                                if y_diff < 0 and abs(y_diff) >= (3 / spacing[1]) and abs(x_diff) <= (3 / spacing[0]):
                                    y_diff = abs(y_diff * spacing[1])
                                    up_sides.append(y_diff)
                                # Down
                                elif y_diff > 0 and abs(y_diff) >= (3 / spacing[1]) and abs(x_diff) <= (3 / spacing[0]):
                                    y_diff = abs(y_diff * spacing[1])
                                    down_sides.append(y_diff)

                            all_sides = [right_sides, left_sides, up_sides, down_sides]
                            max_size = max(len(side_list) for side_list in all_sides)
                            if len(right_sides) < max_size:
                                right_sides.extend([None] * (max_size - len(right_sides)))
                            if len(left_sides) < max_size:
                                left_sides.extend([None] * (max_size - len(left_sides)))
                            if len(up_sides) < max_size:
                                up_sides.extend([None] * (max_size - len(up_sides)))
                            if len(down_sides) < max_size:
                                down_sides.extend([None] * (max_size - len(down_sides)))

                    # Plots sides values
                    sides_fig = plt.figure()
                    sides_colors = ["red", "blue", "green", "orange"]
                    sides_labels = ["Right", "Left", "Up", "Down"]
                    markers = ["o", "s", "D", "v"]
                    ax = sides_fig.add_subplot(1, 1, 1)
                    # Iterate over each type of side and its respective values
                    for i, side_list in enumerate([right_sides, left_sides, up_sides, down_sides]):
                        side_list_filtered = [value for value in side_list if value is not None]
                        mean_value = np.mean(side_list_filtered)
                        mean_std = np.std(side_list_filtered) / np.sqrt(len(side_list_filtered))
                        ax.scatter(
                            z_slices_filtered,
                            side_list,
                            color=sides_colors[i],
                            label=sides_labels[i],
                            marker=markers[i],
                        )
                        plt.axhline(
                            y=mean_value,
                            color=sides_colors[i],
                            label=f"Mean {sides_labels[i]} side: {round(mean_value, 2)} \u00B1 {round(mean_std, 2)} mm",
                        )

                    ax.set_ylabel("Measured distance (mm)")
                    ax.set_xlabel("Slice")
                    ax.set_title(
                        f"Distances along transaxial plane of the cylindrical insert \n Pixel dimension: {round(spacing[0], 2)}x{round(spacing[1], 2)} mm\u00B2"
                    )
                    ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc="upper left")
                    plt.ylim([6, 14])
                    sides_fig.savefig(
                        os.path.join(output_date_analysis_path, f"Distortion{mr_img}_MeshSides.png"),
                        bbox_inches="tight",
                    )
                    plt.close(sides_fig)

                    # Plots (y,x) mesh central coordinates
                    central_fig = plt.figure()
                    central_colors = ["blue", "green"]
                    central_labels = ["y-axis", "x-axis"]
                    markers = ["o", "s"]
                    ax = central_fig.add_subplot(1, 1, 1)

                    ax.scatter(
                        z_slices_filtered,
                        [coord[0] for coord in centers_differences],
                        color=central_colors[0],
                        label=central_labels[0],
                        marker=markers[0],
                    )
                    ax.scatter(
                        z_slices_filtered,
                        [coord[1] for coord in centers_differences],
                        color=central_colors[1],
                        label=central_labels[1],
                        marker=markers[1],
                    )

                    max_y = max(coord[0] for coord in centers_differences)
                    min_y = min(coord[0] for coord in centers_differences)
                    max_x = max(coord[1] for coord in centers_differences)
                    min_x = min(coord[1] for coord in centers_differences)
                    upper_limit = max(max_y, max_x)
                    lower_limit = min(min_y, min_x)
                    margin = 10
                    upper_limit += margin
                    lower_limit -= margin

                    ax.set_ylabel("Mesh centres deviation")
                    ax.set_xlabel("Slice")
                    ax.set_title(
                        f"Mesh centres deviation \n Pixel dimension: {round(spacing[0], 2)}x{round(spacing[1], 2)} mm\u00B2"
                    )
                    ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc="upper left")
                    ax.set_ylim([lower_limit, upper_limit])
                    central_fig.savefig(
                        os.path.join(output_date_analysis_path, f"Distortion{mr_img}_MeshCentres.png"),
                        bbox_inches="tight",
                    )
                    plt.close(central_fig)

                    # Excel document creation
                    writer_sides = pd.ExcelWriter(
                        os.path.join(output_date_analysis_path, f"Distortion{mr_img}_MeshSides.xlsx")
                    )
                    writer_ax = pd.ExcelWriter(
                        os.path.join(output_date_analysis_path, f"Distortion{mr_img}_AxialDistances.xlsx")
                    )

                    # Dictionary and DataFrame creation
                    for n in range(0, len(z_slices_filtered)):
                        if (
                            len(right_sides) != 0
                            and len(left_sides) != 0
                            and len(up_sides) != 0
                            and len(down_sides) != 0
                        ):
                            sides = {"Slice": z_slices_filtered}
                            sides["Right side (mm)"] = right_sides
                            sides["Left side (mm)"] = left_sides
                            sides["Up side (mm)"] = up_sides
                            sides["Down side (mm)"] = down_sides
                            sides["Mesh y-axis center deviation"] = [coord[0] for coord in centers_differences]
                            sides["Mesh x-axis center deviation"] = [coord[1] for coord in centers_differences]

                            df = pd.DataFrame(data=sides)
                            df.to_excel(writer_sides, sheet_name="Mesh sides distortion", index=False)

                    for p in range(0, len(axial_distance)):
                        if len(axial_distance) != 0:
                            dist = {"Axial distances (mm)": axial_distance}
                            df_dist = pd.DataFrame(data=dist)
                            df_dist.to_excel(writer_ax, sheet_name="Axial distortion", index=False)

                    writer_sides.save()
                    writer_ax.save()

                    ## PANAL INSERT ##
                    # Calculates the distance between the different points of the user-defined walls (in mm)
                    distances = []
                    for panal_points_dir in panal_points:
                        distance_dir = []
                        for i in range(len(panal_points_dir) - 1):
                            y1, x1 = panal_points_dir[i]
                            y2, x2 = panal_points_dir[i + 1]
                            distance = math.sqrt(
                                ((y2 - y1) * serie_mr.GetSpacing()[1]) ** 2
                                + ((x2 - x1) * serie_mr.GetSpacing()[0]) ** 2
                            )
                            distance_dir.append(distance)
                        distances.append(distance_dir)

                    # Plots walls distances
                    walls_fig = plt.figure()
                    walls_colors = ["#FF4444", "#55CC55", "#FFCC33", "#3366FF", "#CC6699", "#6666CC"]
                    walls_labels = [
                        "Right",
                        "Left",
                        "Diag. up-right",
                        "Diag. down-left",
                        "Diag. up-left",
                        "Diag. down-right",
                    ]
                    markers = ["o", "s", "D", "v", "h", "p"]
                    ax = walls_fig.add_subplot(1, 1, 1)

                    # Iterate over each type of side and its respective values
                    for i, walls_list in enumerate(distances):
                        walls_list_filtered = [value for value in walls_list if value is not None]
                        mean_value = np.mean(walls_list_filtered)
                        mean_std = np.std(walls_list_filtered) / np.sqrt(len(walls_list_filtered))

                        ax.scatter(
                            range(len(walls_list)),
                            walls_list,
                            color=walls_colors[i],
                            label=walls_labels[i],
                            marker=markers[i],
                        )
                        plt.axhline(
                            y=mean_value,
                            color=walls_colors[i],
                            label=f"Mean {walls_labels[i]} side: {round(mean_value, 2)} \u00B1 {round(mean_std, 2)} mm",
                        )

                    plt.xticks([])
                    ax.set_ylabel("Measured distance (mm)")
                    ax.set_title(
                        f"Distances between panal insert walls \n Pixel dimension: {round(spacing[0], 2)}x{round(spacing[1], 2)} mm\u00B2"
                    )
                    plt.ylim([15, 40])
                    ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc="upper left")
                    ax.set_xticks([])
                    walls_fig.savefig(
                        os.path.join(output_date_analysis_path, f"Distortion{mr_img}_PanalWalls.png"),
                        bbox_inches="tight",
                    )
                    plt.close(walls_fig)

                    # Save distances in each direction in an Excel file
                    writer_panal = pd.ExcelWriter(
                        os.path.join(output_date_analysis_path, f"Distortion{mr_img}_PanalWalls.xlsx")
                    )

                    # Dictionary and DataFrame creation
                    series_r = pd.Series(distances[0], name="Right (mm)")
                    series_l = pd.Series(distances[1], name="Left (mm)")
                    series_dsu = pd.Series(distances[2], name="Diag. up-right (mm)")
                    series_dsd = pd.Series(distances[3], name="Diag. down-left (mm)")
                    series_dmu = pd.Series(distances[4], name="Diag. up-left (mm)")
                    series_dmd = pd.Series(distances[5], name="Diag. down-right (mm)")

                    df = pd.concat([series_r, series_l, series_dsu, series_dsd, series_dmu, series_dmd], axis=1)
                    df.to_excel(writer_panal, sheet_name="Panal sides distortion", index=False)
                    writer_panal.save()

    print("Distortion analysis completed.")

    window.destroy()


window.title("Distortion analysis")

e1 = tk.Button(window, text="Browse directory", command=main.browse_directory)
lbl1 = tk.Label(window, textvariable=main.directory_path)

tk.Label(window, text="Select analysis directory:").grid(row=0, column=0)
e1.grid(row=0, column=1)
lbl1.grid(row=0, column=2)

tk.Button(window, text="Run", command=main_code).grid(row=1, column=1, sticky=tk.W, pady=4)

tk.mainloop()
