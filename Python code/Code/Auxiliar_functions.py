# Quality Assurance analysis for PET/CT and PET/MR systems with the 3D-Printed Phantom: Auxiliar Functions.

import os
import numpy as np
import SimpleITK as sitk
import math
import scipy.ndimage.measurements
import nrrd
import yaml
import pandas as pd
from scipy import ndimage
from skimage import measure
from mahotas.labeled import bwperim
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from radiomics import featureextractor
from skimage.measure import regionprops
from scipy.optimize import curve_fit
from tkinter import filedialog


# Choose axial slices in 3D DICOM image
def select_slices(dicom_path, image_name, message, date):
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(dicom_path)
    reader.SetFileNames(dicom_names)
    im = reader.Execute()
    image = sitk.GetArrayFromImage(im)

    axial_slice_index = 0
    selected_slice = None
    indices = []

    def on_select_button_clicked(event):
        nonlocal selected_slice
        nonlocal indices
        selected_slice = image[axial_slice_index]
        indices.append(axial_slice_index)

    def update_axial_slice():
        nonlocal axial_slice_index
        axial_slice = image[axial_slice_index, :, :]
        ax.imshow(axial_slice, cmap='gray', origin='lower')
        ax.set_title(f'Axial Slice {axial_slice_index} \n of {image_name} image dated {date}')
        ax.invert_yaxis()
        fig.canvas.draw_idle()

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)

    axial_slice = image[axial_slice_index, :, :]
    ax.imshow(axial_slice, cmap='gray', origin='lower')
    ax.set_title(f'Axial Slice {axial_slice_index} \n of {image_name} image dated {date}')
    ax.invert_yaxis()

    slider_axial = Slider(ax=plt.axes([0.15, 0.1, 0.65, 0.03]), label='Axial Slice', valmin=0, valmax=image.shape[0]-1, valinit=axial_slice_index, valstep=1)

    select_button_axial = plt.Button(ax=plt.axes([0.5, 0.02, 0.1, 0.04]), label='Select')
    select_button_axial.on_clicked(on_select_button_clicked)

    def on_slider_change(val):
        nonlocal axial_slice_index
        axial_slice_index = int(val)
        update_axial_slice()

    slider_axial.on_changed(on_slider_change)

    fig.canvas.set_window_title(message)
    plt.show(block=True)

    if selected_slice is not None:
        slices = indices
    else:
        slices = []

    return slices


# QUANTIFICATION analysis #
# Creates analysis region from a slice for QUANTIFICATION analysis (PET, CT and MR)
def segmentation_region(im, coord, spacing, modality, inserts_type):
    image = sitk.GetArrayFromImage(im)

    distance_pet = 60

    if inserts_type == 'Calibration Curve':
        distance_ct = 60
    elif inserts_type == 'Materials Characterization':
        distance_ct = 21.5
        distance_mr = 21.5

    segm = np.zeros_like(image)
    segm[coord, :, :] = 1

    seg = segm
    if modality == 'PET':
        desplaz = int(round((distance_pet / spacing) - 1))
        for i in range(1, desplaz + 1):
            seg_roll = np.roll(segm, i, axis=0)
            seg = seg + seg_roll
    elif modality == 'CT':
        desplaz = int(round((distance_ct / spacing) - 1))
        for i in range(1, desplaz + 1):
            seg_roll = np.roll(segm, i, axis=0)
            seg = seg + seg_roll
    elif 'MR' in modality:
        desplaz = int(round((distance_mr / spacing) - 1))
        for i in range(1, desplaz + 1):
            seg_roll = np.roll(segm, i, axis=0)
            seg = seg + seg_roll

    segm_sitk = sitk.GetImageFromArray(seg)
    segm_sitk.SetDirection(im.GetDirection())
    segm_sitk.SetSpacing(im.GetSpacing())
    segm_sitk.SetOrigin(im.GetOrigin())

    final_coord = coord + desplaz

    return segm_sitk, final_coord


# Generate cylindrical segmentation from its centre
def generate_cylinder_segmentation(image, y_center, x_center, radius, ini, fin, voxeldim):
    seg = np.zeros_like(image)

    z_indices, y_indices, x_indices = np.indices(image.shape, dtype=float)
    x_indices -= x_center
    y_indices -= y_center

    distance_to_center_mm = np.sqrt((x_indices*voxeldim[0]) ** 2 + (y_indices*voxeldim[1]) ** 2)
    mask = (distance_to_center_mm <= radius) & np.logical_and(ini <= z_indices, z_indices <= fin)
    seg[mask] = 1

    return seg


# Obtain the centers of the CT inserts and insert segmentations geometrically
def inserts_ct(im, initial_coord, final_coord, voxeldim, inserts_type):
    image = sitk.GetArrayFromImage(im)

    # Calculate the phantom center in xy plane
    if inserts_type == 'Calibration Curve':
        centre_slice = int(round(final_coord + (15 / float(voxeldim[2]))))
    elif inserts_type == 'Materials Characterization':
        centre_slice = int(round(final_coord + (10 / float(voxeldim[2]))))
    slice = image[centre_slice, :, :]

    thresh = np.where(slice >= -644, 1, 0)

    labeled_array, num_features = ndimage.label(thresh)
    volumes = []
    for i in range(1, num_features + 1):
        single_segmentation = np.zeros_like(thresh)
        single_segmentation[labeled_array == i] = 1
        numvoxels = len(single_segmentation[single_segmentation == 1])
        volumes.append(numvoxels)

    sorted_volumes = sorted(volumes, reverse=True)
    wanted_volume = sorted_volumes[0]
    for i, volume in enumerate(volumes):
        if volume == wanted_volume:
            single_segmentation = np.zeros_like(thresh)
            single_segmentation[labeled_array == i + 1] = 1
            panal = ndimage.binary_fill_holes(single_segmentation).astype(np.uint8)

    center_phantom = scipy.ndimage.measurements.center_of_mass(panal)   # phantom's centre in xy plane

    # Obtain CT inserts centres (z,y,x) from the centre of the phantom
    central_coord = int(round((initial_coord + final_coord) / 2))
    center_1 = np.array((central_coord, round(center_phantom[0] - ((18 * 1.17188)/float(voxeldim[1]))), round(center_phantom[1] + ((83 * 1.17188)/float(voxeldim[0])))))
    center_2 = np.array((central_coord, round(center_phantom[0] + ((54 * 1.17188)/float(voxeldim[1]))), round(center_phantom[1] + ((67 * 1.17188)/float(voxeldim[0])))))
    center_3 = np.array((central_coord, round(center_phantom[0] + ((86 * 1.17188)/float(voxeldim[1]))), round(center_phantom[1] + ((0 * 1.17188) / float(voxeldim[0])))))
    center_4 = np.array((central_coord, round(center_phantom[0] + ((54 * 1.17188)/float(voxeldim[1]))), round(center_phantom[1] - ((66 * 1.17188)/float(voxeldim[0])))))
    center_5 = np.array((central_coord, round(center_phantom[0] - ((17 * 1.17188)/float(voxeldim[1]))), round(center_phantom[1] - ((83 * 1.17188)/float(voxeldim[0])))))
    center_6 = np.array((central_coord, round(center_phantom[0] - ((75 * 1.17188)/float(voxeldim[1]))), round(center_phantom[1] - ((37 * 1.17188)/float(voxeldim[0])))))
    center_7 = np.array((central_coord, round(center_phantom[0] - ((75 * 1.17188)/float(voxeldim[1]))), round(center_phantom[1] + ((37 * 1.17188)/float(voxeldim[0])))))
    center_8 = np.array((central_coord, round(center_phantom[0] + ((1 * 1.17188)/float(voxeldim[1]))), round(center_phantom[1] + ((34 * 1.17188)/float(voxeldim[0])))))
    center_9 = np.array((central_coord, round(center_phantom[0] + ((1 * 1.17188)/float(voxeldim[1]))), round(center_phantom[1] - ((34 * 1.17188)/float(voxeldim[0])))))
    center_10 = np.array((central_coord, round(center_phantom[0] - ((33 * 1.17188)/float(voxeldim[1]))), round(center_phantom[1])))
    center_11 = np.array((central_coord, round(center_phantom[0] + ((36 * 1.17188)/float(voxeldim[1]))), round(center_phantom[1])))

    center_inserts = [center_1, center_2, center_3, center_4, center_5, center_6, center_7, center_8, center_9, center_10, center_11]
    radius_inserts = [5, 14, 14, 14, 14, 10, 6, 14, 14, 14, 14]

    # Numbers the inserts as they have been defined
    name_all_inserts = []
    for i in range(1, 12):
        name_all_inserts.append(i)

    # Generates cylindrical inserts' segmentations
    inserts = []
    num = 0
    for center in center_inserts:
        # Cylindrical segmentation
        cyl = generate_cylinder_segmentation(image, center[1], center[2], radius_inserts[num], initial_coord, final_coord, voxeldim)
        cyl_sitk = sitk.GetImageFromArray(cyl)
        cyl_sitk.SetDirection(im.GetDirection())
        cyl_sitk.SetSpacing(im.GetSpacing())
        cyl_sitk.SetOrigin(im.GetOrigin())

        inserts.append(cyl_sitk)
        num += 1

    return name_all_inserts, center_inserts, inserts


# Obtain the centres of the PET inserts, segmentations around them and insert segmentations
def inserts_pet(im, segm, initial_coord, final_coord, voxeldim):
    image = sitk.GetArrayFromImage(im)

    # Calculate the phantom center in xy plane
    centre_slice = int(round(final_coord + (20 / float(voxeldim[2]))))
    print(f'Centre slice: {centre_slice}')
    slice = image[centre_slice, :, :]

    thresh = np.where(slice >= 1000, 1, 0)

    center_phantom = scipy.ndimage.measurements.center_of_mass(thresh)  # phantom's centre in xy plane

    # Obtain PET inserts centres (z,y,x) from the centre of the phantom
    central_coord = int(round((initial_coord + final_coord) / 2))
    center_1 = np.array((central_coord, round(center_phantom[0] - ((18 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] + ((83 * 1.17188) / float(voxeldim[0])))))
    center_2 = np.array((central_coord, round(center_phantom[0] + ((54 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] + ((66 * 1.17188) / float(voxeldim[0])))))
    center_3 = np.array((central_coord, round(center_phantom[0] + ((87 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] + ((2 * 1.17188) / float(voxeldim[0])))))
    center_4 = np.array((central_coord, round(center_phantom[0] + ((55 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] - ((66 * 1.17188) / float(voxeldim[0])))))
    center_5 = np.array((central_coord, round(center_phantom[0] - ((16 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] - ((83 * 1.17188) / float(voxeldim[0])))))
    center_6 = np.array((central_coord, round(center_phantom[0] - ((74 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] - ((38 * 1.17188) / float(voxeldim[0])))))
    center_7 = np.array((central_coord, round(center_phantom[0] - ((75 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] + ((36 * 1.17188) / float(voxeldim[0])))))
    center_inserts = [center_1, center_2, center_3, center_4, center_5, center_6, center_7]

    # Numbers the inserts as they have been defined
    name_all_inserts = []
    for i in range(1, 8):
        name_all_inserts.append(i)

    # Generates cylindrical segmentations containing each insert and inserts' segmentations
    inserts = []
    for center in center_inserts:
        # Cylindrical segmentation
        cyl = generate_cylinder_segmentation(image, center[1], center[2], 23, initial_coord, final_coord, voxeldim)
        cyl_sitk = sitk.GetImageFromArray(cyl)
        cyl_sitk.SetDirection(im.GetDirection())
        cyl_sitk.SetSpacing(im.GetSpacing())
        cyl_sitk.SetOrigin(im.GetOrigin())

        # Insert's segmentation
        insert = region_growing(im, segm, cyl_sitk, center, 4, 'PET')
        inserts.append(insert)

    return name_all_inserts, center_inserts, inserts


# Obtain the centers of the MR inserts and insert segmentations geometrically
def inserts_mr(im, initial_coord, final_coord, voxeldim):
    image = sitk.GetArrayFromImage(im)

    # Obtain phantom's centre in xy plane
    centre_slice = int(round(final_coord + (16 / float(voxeldim[2]))))  # slice to calculate phantom's centre
    otsu_filter = sitk.OtsuThresholdImageFilter()
    otsu_filter.SetInsideValue(0)
    otsu_filter.SetOutsideValue(1)
    otsu_result = otsu_filter.Execute(im[:, :, centre_slice])
    otsu_threshold = otsu_filter.GetThreshold()
    thresh = np.where(image[centre_slice, :, :] >= otsu_threshold, 1, 0)
    labeled_array, num_features = ndimage.label(thresh)
    sizes = ndimage.sum(thresh, labeled_array, range(num_features + 1))
    largest_label = np.argmax(sizes[1:]) + 1
    water_region = np.zeros_like(thresh)
    water_region[labeled_array == largest_label] = 1
    panal = ndimage.binary_fill_holes(water_region).astype(np.uint8)
    center_phantom = ndimage.measurements.center_of_mass(panal)  # phantom's centre (y,x)

    # Obtain MR inserts centres (z,y,x) from the centre of the phantom
    central_coord = int(round((initial_coord + final_coord) / 2))
    center_1 = np.array((central_coord, round(center_phantom[0] - ((18 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] + ((83 * 1.17188) / float(voxeldim[0])))))
    center_2 = np.array((central_coord, round(center_phantom[0] + ((54 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] + ((67 * 1.17188) / float(voxeldim[0])))))
    center_3 = np.array((central_coord, round(center_phantom[0] + ((86 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] + ((0 * 1.17188) / float(voxeldim[0])))))
    center_4 = np.array((central_coord, round(center_phantom[0] + ((54 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] - ((66 * 1.17188) / float(voxeldim[0])))))
    center_5 = np.array((central_coord, round(center_phantom[0] - ((17 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] - ((83 * 1.17188) / float(voxeldim[0])))))
    center_6 = np.array((central_coord, round(center_phantom[0] - ((75 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] - ((37 * 1.17188) / float(voxeldim[0])))))
    center_7 = np.array((central_coord, round(center_phantom[0] - ((75 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] + ((37 * 1.17188) / float(voxeldim[0])))))
    center_8 = np.array((central_coord, round(center_phantom[0] + ((1 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] + ((34 * 1.17188) / float(voxeldim[0])))))
    center_9 = np.array((central_coord, round(center_phantom[0] + ((1 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1] - ((34 * 1.17188) / float(voxeldim[0])))))
    center_10 = np.array((central_coord, round(center_phantom[0] - ((33 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1])))
    center_11 = np.array((central_coord, round(center_phantom[0] + ((36 * 1.17188) / float(voxeldim[1]))), round(center_phantom[1])))

    center_inserts = [center_1, center_2, center_3, center_4, center_5, center_6, center_7, center_8, center_9, center_10, center_11]
    radius_inserts = [5, 14, 14, 14, 14, 10, 6, 14, 14, 14, 14]

    # Numbers the inserts as they have been defined
    name_all_inserts = []
    for i in range(1, 12):
        name_all_inserts.append(i)

    # Generates cylindrical inserts' segmentations
    inserts = []
    num = 0
    for center in center_inserts:
        # Cylindrical segmentation
        cyl = generate_cylinder_segmentation(image, center[1], center[2], radius_inserts[num], initial_coord, final_coord, voxeldim)
        cyl_sitk = sitk.GetImageFromArray(cyl)
        cyl_sitk.SetDirection(im.GetDirection())
        cyl_sitk.SetSpacing(im.GetSpacing())
        cyl_sitk.SetOrigin(im.GetOrigin())

        inserts.append(cyl_sitk)
        num += 1

    return name_all_inserts, center_inserts, inserts


# Creates a reduced cylindrical segmentation from an input segmentation
def reduced_segmentation_cyl(segm, pixel_x):
    segmentation = sitk.GetArrayFromImage(segm)
    centre_mass = scipy.ndimage.measurements.center_of_mass(segmentation)
    center = np.array(centre_mass)  # center of mass of segmentation

    if not np.any(segmentation == 1):
        cylinder_segmentation = np.zeros_like(segmentation)
        diameter = np.nan

    else:
        znum = np.sum(segmentation, axis=0)
        xnum = np.sum(segmentation, axis=2)
        diameter = np.max(xnum) * pixel_x

        radius_cylinder = np.max(xnum) / 4
        height_cylinder = np.max(znum) / 2

        shape = segmentation.shape

        z, y, x = np.indices(shape)
        z = z - center[0]
        y = y - center[1]
        x = x - center[2]

        distance_to_center = np.sqrt(x**2 + y**2)

        # Generates a boolean array indicating the positions inside the cylinder
        cylinder_mask = (distance_to_center <= radius_cylinder) & (z >= (-height_cylinder/2)) & (z <= (height_cylinder/2))

        # Creates a new cylindrical segmentation with the same shape as the original segmentation
        cylinder_segmentation = np.zeros_like(segmentation)

        # Generates cylindrical segmentation and saves it
        cylinder_segmentation[cylinder_mask] = 1

    reduced_segm_sitk = sitk.GetImageFromArray(cylinder_segmentation)
    reduced_segm_sitk.SetDirection(segm.GetDirection())
    reduced_segm_sitk.SetSpacing(segm.GetSpacing())
    reduced_segm_sitk.SetOrigin(segm.GetOrigin())

    return reduced_segm_sitk, diameter


# Quantification analysis in a segmentation (all modalities): mean value and its stdev
def quantification(img, segm):
    img = np.copy(img).astype('float')
    seg = np.copy(segm)

    if not np.any(seg == 1):
        mean_value = np.nan
        dev_value = np.nan
    else:
        seg = np.copy(seg).astype(bool)
        img = img[seg]  # defines image in the segmentation region only
        mean_value = np.mean(img)  # mean value in segmentation area
        dev_value = np.std(img)  # deviation value
        #num_voxels = np.size(img)  # number of voxels in the segmentation
        #mean_dev_value = dev_value / math.sqrt(num_voxels)  # mean deviation value

    return mean_value, dev_value


# Obtains the image max value in a segmentation
def maxvalue(image, seg, percent):
    img = np.copy(image).astype('float')
    seg = np.copy(seg).astype('bool')
    img[~seg] = np.nan

    max_value = np.nanmax(img)
    max_value_percent = np.nanmean(img[img >= (percent/100) * max_value])  # mean of the voxels with value >= percent*max

    return max_value, max_value_percent


def get_nbhd(pt, checked, dims):
    nbhd = []

    if (pt[0] > 0) and not checked[pt[0]-1, pt[1], pt[2]]:
        nbhd.append((pt[0]-1, pt[1], pt[2]))
    if (pt[1] > 0) and not checked[pt[0], pt[1]-1, pt[2]]:
        nbhd.append((pt[0], pt[1]-1, pt[2]))
    if (pt[2] > 0) and not checked[pt[0], pt[1], pt[2]-1]:
        nbhd.append((pt[0], pt[1], pt[2]-1))

    if (pt[0] < dims[0]-1) and not checked[pt[0]+1, pt[1], pt[2]]:
        nbhd.append((pt[0]+1, pt[1], pt[2]))
    if (pt[1] < dims[1]-1) and not checked[pt[0], pt[1]+1, pt[2]]:
        nbhd.append((pt[0], pt[1]+1, pt[2]))
    if (pt[2] < dims[2]-1) and not checked[pt[0], pt[1], pt[2]+1]:
        nbhd.append((pt[0], pt[1], pt[2]+1))

    return nbhd


# Region growing function
def region_growing(im, ar, ins, start_coord, t, modality):
    img = sitk.GetArrayFromImage(im).astype(float)
    area = sitk.GetArrayFromImage(ar).astype(bool)
    insert = sitk.GetArrayFromImage(ins).astype(bool)
    img[~area] = np.nan

    seg = np.zeros(img.shape, dtype=np.bool)
    checked = np.zeros_like(seg)

    start_coord = tuple(map(int, start_coord))

    seg[start_coord] = True
    checked[start_coord] = True
    needs_check = get_nbhd(start_coord, checked, img.shape)

    if modality == 'CT':
        if img[start_coord] >= -200:
            seg[start_coord] = True
        else:
            seg[start_coord] = False

        checked[start_coord] = True
        needs_check = get_nbhd(start_coord, checked, img.shape)

    elif modality == 'PET':
        seg[start_coord] = True
        checked[start_coord] = True
        needs_check = get_nbhd(start_coord, checked, img.shape)

        max_insertpet, maxperc_insertpet = maxvalue(img, insert, 70)
        th40 = 0.4 * maxperc_insertpet
        th = th40

    while len(needs_check) > 0:
        pt = needs_check.pop()
        # It is possible that the point was already checked and was put in the needs_check stack multiple time
        if checked[pt]:
            continue

        checked[pt] = True

        # Handle borders
        imin = max(pt[0] - t, 0)
        imax = min(pt[0] + t, img.shape[0] - 1)
        jmin = max(pt[1] - t, 0)
        jmax = min(pt[1] + t, img.shape[1] - 1)
        kmin = max(pt[2] - t, 0)
        kmax = min(pt[2] + t, img.shape[2] - 1)

        if modality == 'PET':
            if img[pt] >= th:
                seg[pt] = True
                needs_check += get_nbhd(pt, checked, img.shape)
            else:
                if img[imin:imax + 1, jmin:jmax + 1, kmin:kmax + 1].mean() > th:
                    # Include the voxel in the segmentation and add its neighbors to be checked
                    seg[pt] = True
                    needs_check += get_nbhd(pt, checked, img.shape)

        elif modality == 'CT':
            if img[pt] >= -200:
                seg[pt] = True
                needs_check += get_nbhd(pt, checked, img.shape)
            else:
                if img[imin:imax + 1, jmin:jmax + 1, kmin:kmax + 1].mean() > -200:
                    seg[pt] = True
                    needs_check += get_nbhd(pt, checked, img.shape)

    seg = np.copy(seg).astype(int)
    segm_sitk = sitk.GetImageFromArray(seg)
    segm_sitk.SetDirection(im.GetDirection())
    segm_sitk.SetSpacing(im.GetSpacing())
    segm_sitk.SetOrigin(im.GetOrigin())

    return segm_sitk


# Segmentation's volume calculator
def seg_vol(seg, voxeldim):
    segm = sitk.GetArrayFromImage(seg)
    numvoxels = np.sum(segm)
    vol = numvoxels * (voxeldim / 1000)     # segmentation's volume in ml

    return vol


# Convert HHMMSS string to seconds
def time_converter(time_string):
    hhmmss = time_string
    time_in_seconds = int(hhmmss[-2:]) + int(hhmmss[2:4]) * 60 + int(hhmmss[:2]) * 3600

    return time_in_seconds



# RESOLUTION analysis #
# Choose triangular vertices
def choose_vertices(dicom_path, coord, image_name, message, section):
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(dicom_path)
    reader.SetFileNames(dicom_names)
    im = reader.Execute()
    image = sitk.GetArrayFromImage(im)

    axial_slice = image[coord, :, :]

    fig, ax = plt.subplots()
    ax.imshow(axial_slice, cmap='gray', origin='lower')
    ax.set_title(f'Select Vertices Coordinates of \n {image_name} image Section {section}')
    ax.invert_yaxis()

    coordinates = []

    def on_image_click(event):
        if event.inaxes == ax:
            x = int(event.xdata)
            y = int(event.ydata)
            coordinates.append([y, x])
            ax.plot(x, y, 'ro', markersize=3)
            plt.draw()

    fig.canvas.mpl_connect('button_press_event', on_image_click)
    fig.canvas.set_window_title(message)
    plt.show(block=True)

    return coordinates


# Resolution analysis for PET images
def resol_pet(im, coord, vertices, output, voxeldim):
    image = sitk.GetArrayFromImage(im)
    shape = image.shape

    # Creates 6 triangular segmentations and analyses the resolution in each section #
    # Defines vertices coordinates of each triangular segmentation in the xy plane (z=coord)
    vertices6 = np.array(vertices[5])
    vertices5 = np.array(vertices[4])
    vertices4 = np.array(vertices[3])
    vertices3 = np.array(vertices[2])
    vertices2 = np.array(vertices[1])
    vertices1 = np.array(vertices[0])
    list_vertices = [vertices1, vertices2, vertices3, vertices4, vertices5, vertices6]
    list_diameters = [5, 7.5, 9, 11, 12, 15]

    # Generates each triangular section
    contrast_total = []
    section_names = []
    holes_total = []
    rc_holes_total = []
    theoretical_holes_total = [55, 31, 20, 15, 10, 8]
    # SECTION-BY-SECTION loop
    for n in range(len(list_vertices)):
        tri_segm = np.zeros_like(image)
        vertices = list_vertices[n]

        # Obtain the (y,x) coordinates of the triangle vertices
        A, B, C = vertices
        y1, x1 = A
        y2, x2 = B
        y3, x3 = C

        # Calculate original triangle area
        area = abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))

        # Iterate over all coordinates in the image
        for i in range(shape[0]):
            if i == coord:
                for j in range(shape[1]):
                    for k in range(shape[2]):
                        # Calculate the area of each of the triangles formed by the point and each of the sides of the original triangle
                        sub_area1 = abs((x1 - k) * (y2 - j) - (x2 - k) * (y1 - j))
                        sub_area2 = abs((x2 - k) * (y3 - j) - (x3 - k) * (y2 - j))
                        sub_area3 = abs((x3 - k) * (y1 - j) - (x1 - k) * (y3 - j))

                        # Check whether the sum of the areas of the triangles formed by the point and the sides of the original triangle is equal to the area of the original triangle
                        if abs(sub_area1 + sub_area2 + sub_area3 - area) < 1e-6:
                            tri_segm[i, j, k] = 1

        # Expand triangular segmentation to previous and next slice
        tri = tri_segm
        for l in range(-1, 0):
            seg_roll = np.roll(tri_segm, l, axis=0)
            tri = tri + seg_roll

        for l in range(1, 2):
            seg_roll = np.roll(tri_segm, l, axis=0)
            tri = tri + seg_roll

        # Analyse the resolution in triangular section #
        tri_segme = np.copy(tri).astype(bool)
        image_copy = np.copy(image).astype(float)
        image_copy[~tri_segme] = np.nan  # image only in triangular section
        max_value = np.nanmax(image_copy)  # max intensity value in triangular section
        # max_value70 = np.nanmean(image_copy[image_copy >= (70/100) * max_value])
        # th40 = 0.40 * max_value70
        th40 = 0.46 * max_value  # 46% threshold of the maxvalue

        # Create holes segmentation
        thresh_phantom = np.where((image >= 0.1 * max_value), 1, 0)
        thresh_phantom_slice_filled = ndimage.binary_fill_holes(thresh_phantom[coord, :, :])
        thresh_phantom_erosion = ndimage.binary_erosion(thresh_phantom_slice_filled, iterations=6)
        thresh_phantom_er = np.zeros_like(thresh_phantom, dtype=int)
        thresh_phantom_er[coord, :, :] = thresh_phantom_erosion
        thresh_phantom_er[coord + 1, :, :] = thresh_phantom_erosion
        thresh_phantom_er[coord - 1, :, :] = thresh_phantom_erosion

        thresh_holes = np.where((image_copy < th40), 1, 0).astype(int)
        th_holes = thresh_holes & thresh_phantom_er

        # Contrast calculation
        c_u40 = np.nanmean(image_copy[image_copy > th40])
        image_copy_c = np.copy(image_copy).astype(float)
        thresh_phantom_er = np.copy(thresh_phantom_er).astype(bool)
        image_copy_c[~thresh_phantom_er] = np.nan
        c_l40 = np.nanmean(image_copy_c[image_copy_c < th40])
        contrast = c_u40 / c_l40
        if contrast > 10:
            contrast = round(contrast)
        else:
            contrast = round(contrast, 1)
        contrast_total.append(contrast)
        section_names.append(str(n + 1))

        # Separate each hole of the section for INTENSITY PROFILE creation
        labeled_array, num_features = measure.label(th_holes, background=0, return_num=True, connectivity=None)
        holes_total.append(num_features)
        rc_holes_total.append(num_features / theoretical_holes_total[n])

        if num_features != 0:
            dist_total = []
            centroid_total = []
            for i in range(1, num_features + 1):
                single_hole = np.zeros_like(th_holes)
                single_hole[labeled_array == i] = 1

                z_coords, y_coords, x_coords = np.where(single_hole == 1)
                centroid_z, centroid_y, centroid_x = np.mean(z_coords), np.mean(y_coords), np.mean(x_coords)
                centroid = int(centroid_z), int(centroid_y), int(centroid_x)
                dist = np.sqrt(((centroid_x - x1)*voxeldim[0]) ** 2 + ((centroid_y - y1)*voxeldim[1]) ** 2)
                dist_total.append(dist)
                centroid_total.append(centroid)

            min_value_distance = min(dist_total)
            if min_value_distance <= 95:
                min_index = dist_total.index(min_value_distance)
                initial_centroid = centroid_total[min_index]  # we choose the hole closest to vertex A
                starting_coord = initial_centroid

                # Create intensity profile
                yx_vector = B - A
                yx_vector[0] = int(yx_vector[0])
                yx_vector[1] = int(yx_vector[1])
                xyz_vector = np.insert(yx_vector, 0, 0)

                xyz_final = xyz_vector + initial_centroid
                end_coord = tuple(xyz_final)

                # Calculate the original distance between the two points
                dist_original = math.sqrt((end_coord[1] - starting_coord[1]) ** 2 + (end_coord[2] - starting_coord[2]) ** 2)

                # Reduce the line by 30%
                dist_new = dist_original * 0.7

                # Calculate the coordinates of the new end point
                x_new = starting_coord[2] + (end_coord[2] - starting_coord[2]) * dist_new / dist_original
                y_new = starting_coord[1] + (end_coord[1] - starting_coord[1]) * dist_new / dist_original

                new_end_coord = tuple([end_coord[0], y_new, x_new])

                z_lin, y_lin, x_lin = np.linspace(starting_coord[0], new_end_coord[0], num=60), np.linspace(starting_coord[1], new_end_coord[1], num=60), np.linspace(starting_coord[2], new_end_coord[2], num=60)

                intensity_profile = [(image[int(round(z_lin[i])), int(round(y_lin[i])), int(round(x_lin[i]))]) for i in range(60)]
                max_intensity = max(intensity_profile)
                normalized_profile = [value / max_intensity for value in intensity_profile]

                fig, ax = plt.subplots(nrows=1, ncols=2)
                fig.suptitle(f'Resolution PET Section {n+1} (D={list_diameters[n]}mm)')
                ax[0].plot(normalized_profile)
                ax[0].set_title(f"Intensity Profile (From Image Center to Edges)")
                ax[0].set_ylabel("Normalized Intensity")

                ax[1].imshow(image[starting_coord[0], :, :], cmap='gray')
                ax[1].imshow(th_holes[starting_coord[0], :, :], cmap='Blues', alpha=0.3)
                ax[1].plot([starting_coord[2], new_end_coord[2]], [starting_coord[1], new_end_coord[1]], color='red', linewidth=1)
                ax[1].set_title(f"Profile Line")

                fig.savefig(os.path.join(output, 'ResolutionPET_IntensityProfile_Section' + str(n + 1) + '.png'))
                plt.close(fig)

    return contrast_total, holes_total, rc_holes_total, section_names


# Resolution analysis for CT images
def resol_ct(im, coord, vertices, output, voxeldim):
    image = sitk.GetArrayFromImage(im)

    # Obtain phantom's centre
    thresh_all = np.where(-250 <= image[coord, :, :], 1, 0).astype(int)
    labeled_array, num_features = ndimage.label(thresh_all)
    sizes = ndimage.sum(thresh_all, labeled_array, range(num_features + 1))
    largest_label = np.argmax(sizes[1:]) + 1
    all_region = np.zeros_like(thresh_all)
    all_region[labeled_array == largest_label] = 1
    center_phantom = ndimage.measurements.center_of_mass(all_region)  # phantom's centre (y,x)

    # Select all HOLES of the resolution insert
    # Obtain holes threshold
    shape = image.shape
    z, y, x = np.indices(shape)
    y = y - (center_phantom[0])
    x = x - (center_phantom[1])
    distance_to_center = np.sqrt((x * voxeldim[0]) ** 2 + (y * voxeldim[1]) ** 2)   # distance in mm
    max_radius = 100  # maximum radius in mm
    mask = distance_to_center <= max_radius
    mask = np.copy(mask).astype(bool)
    image_masked = np.zeros_like(image).astype(float)
    image_masked[mask] = image[mask]
    image_masked[~mask] = np.nan
    valid_pixels = image_masked[coord, :, :][~np.isnan(image_masked[coord, :, :])]

    # Double gaussian function definition to obtain rods threshold
    def double_gaussian(x, a1, b1, c1, a2, b2, c2):
        return (a1 * np.exp(-(x - b1) ** 2 / (2 * c1 ** 2)) + a2 * np.exp(-(x - b2) ** 2 / (2 * c2 ** 2)))

    hist, bin_edges = np.histogram(valid_pixels, bins=400, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    initial_guesses = [1, np.mean(valid_pixels), np.std(valid_pixels), 1, np.mean(valid_pixels) + 70, np.std(valid_pixels) / 2]
    params, covariance = curve_fit(double_gaussian, bin_centers, hist, p0=initial_guesses)
    a1, b1, c1, a2, b2, c2 = params
    centroid1 = b1
    centroid2 = b2
    sigma2 = abs(c2)  # Cylinders' gaussian standard deviation
    num_std = 1
    th_cyl = centroid2 - num_std * sigma2   # cylinders' calculated threshold

    # Show double gaussian function
    # x_fit = np.linspace(min(valid_pixels), max(valid_pixels), 1000)
    # y_fit = double_gaussian(x_fit, *params)
    #
    # plt.figure(figsize=(10, 6))
    # plt.hist(valid_pixels, bins=100, density=True, alpha=0.6, color='g', label='Datos')
    # plt.plot(x_fit, y_fit, label='Ajuste Doble Gaussiana', color='blue')
    # plt.axvline(centroid1, color='red', linestyle='dashed', linewidth=1, label=f'Centroide Agua: {centroid1:.2f}')
    # plt.axvline(centroid2, color='orange', linestyle='dashed', linewidth=1, label=f'Centroide Cilindros: {centroid2:.2f}')
    # plt.axvline(th_cyl, color='orange', linestyle='dotted', linewidth=1, label=f'Límite Inferior Cilindros: {th_cyl:.2f}')
    # plt.legend()
    # plt.xlabel('Intensidad de píxel')
    # plt.ylabel('Densidad')
    # plt.title('Ajuste de Doble Gaussiana')
    # plt.show()

    thresh_insert = np.where(-250 <= image, 1, 0).astype(int)
    thresh_holes = np.where(image >= th_cyl, 1, 0).astype(int)
    thresh_water = thresh_insert - thresh_holes
    thresh_water = np.copy(thresh_water).astype(float)

    distance_to_center = np.sqrt((x*voxeldim[0]) ** 2 + (y*voxeldim[1]) ** 2)
    max_radius = 130  # maximum radius in mm
    mask = distance_to_center <= max_radius
    thresh_water[~mask] = np.nan

    th_w = np.zeros_like(image)
    th_water_filled = np.zeros_like(image)
    th_w[coord, :, :] = thresh_water[coord, :, :]

    # th_w_dil = ndimage.binary_dilation(th_w, iterations=2)
    # th_w_dil_erosion = ndimage.binary_erosion(th_w_dil, iterations=1).astype(int)

    th_w_dil_fill = ndimage.binary_fill_holes(th_w[coord, :, :])
    th_w_dil_fill_er = ndimage.binary_erosion(th_w_dil_fill, iterations=3)

    th_water_filled[coord, :, :] = th_w_dil_fill_er
    th_water_filled[coord + 1, :, :] = th_w_dil_fill_er
    th_water_filled[coord - 1, :, :] = th_w_dil_fill_er  # water segmentation

    th_holes = thresh_holes & th_water_filled
    th_holes = np.copy(th_holes).astype(float)  # holes' segmentation

    # Creates 6 triangular segmentations and analyses the resolution in each section #
    # Defines vertices coordinates of each triangular segmentation in the xy plane (z=coord)
    vertices6 = np.array(vertices[5])
    vertices5 = np.array(vertices[4])
    vertices4 = np.array(vertices[3])
    vertices3 = np.array(vertices[2])
    vertices2 = np.array(vertices[1])
    vertices1 = np.array(vertices[0])
    list_vertices = [vertices1, vertices2, vertices3, vertices4, vertices5, vertices6]
    list_diameters = [5, 7.5, 9, 11, 12, 15]

    contrast_total = []
    section_names = []
    holes_total = []
    rc_holes_total = []
    theoretical_holes_total = [55, 31, 20, 15, 10, 8]
    for n in range(len(list_vertices)):
        tri_segm = np.zeros_like(image)
        vertices = list_vertices[n]

        # Obtains vertices (y,x) coordinates
        A, B, C = vertices
        y1, x1 = A
        y2, x2 = B
        y3, x3 = C

        # Calculate original triangle area
        area = abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))

        # Iterate over all coordinates in the image
        for i in range(shape[0]):
            if i == coord:
                for j in range(shape[1]):
                    for k in range(shape[2]):
                        # Calculate the area of each of the triangles formed by the point and each of the sides of the original triangle
                        sub_area1 = abs((x1 - k) * (y2 - j) - (x2 - k) * (y1 - j))
                        sub_area2 = abs((x2 - k) * (y3 - j) - (x3 - k) * (y2 - j))
                        sub_area3 = abs((x3 - k) * (y1 - j) - (x1 - k) * (y3 - j))

                        # Check whether the sum of the areas of the triangles formed by the point and the sides of the original triangle is equal to the area of the original triangle
                        if abs(sub_area1 + sub_area2 + sub_area3 - area) < 1e-6:
                            tri_segm[i, j, k] = 1

        # Expand segmentation to previous and next slice
        tri = tri_segm
        for l in range(-1, 0):
            seg_roll = np.roll(tri_segm, l, axis=0)
            tri = tri + seg_roll

        for l in range(1, 2):
            seg_roll = np.roll(tri_segm, l, axis=0)
            tri = tri + seg_roll

        # Analyse the resolution in triangular section #
        segme = np.copy(tri).astype('bool')
        image_copy = np.copy(image).astype('float')
        image_copy[~segme] = np.nan

        th_h = np.copy(th_holes).astype(float)  # holes' segmentation
        th_h[~segme] = 0

        # Contrast calculation
        c_background = np.nanmean(image_copy[image_copy < th_cyl])
        image_copy_c = np.copy(image_copy)
        thresh_water_dil_er = np.copy(th_water_filled).astype(bool)
        image_copy_c[~thresh_water_dil_er] = np.nan
        c_holes = np.nanmean(image_copy_c[image_copy_c >= th_cyl])
        contrast = c_holes - c_background

        contrast = abs(contrast)
        if contrast > 10:
            contrast = round(contrast)
        else:
            contrast = round(contrast, 1)
        contrast_total.append(contrast)
        section_names.append(str(n + 1))

        # Separate each hole of the section for INTENSITY PROFILE creation
        labeled_array, num_features = measure.label(th_h, background=0, return_num=True, connectivity=None)
        holes_total.append(num_features)
        rc_holes_total.append(num_features / theoretical_holes_total[n])

        dist_total = []
        centroid_total = []
        for i in range(1, num_features + 1):
            single_hole = np.zeros_like(th_h)
            single_hole[labeled_array == i] = 1

            z_coords, y_coords, x_coords = np.where(single_hole == 1)
            centroid_z, centroid_y, centroid_x = np.mean(z_coords), np.mean(y_coords), np.mean(x_coords)
            centroid = int(centroid_z), int(centroid_y), int(centroid_x)
            dist = np.sqrt((centroid_x - x1) ** 2 + (centroid_y - y1) ** 2)
            dist_total.append(dist)
            centroid_total.append(centroid)

        min_value_distance = min(dist_total)
        min_index = dist_total.index(min_value_distance)
        initial_centroid = centroid_total[min_index]  # we choose the hole closest to vertex A
        starting_coord = initial_centroid

        # Create intensity profile
        yx_vector = B - A
        yx_vector[0] = int(yx_vector[0])
        yx_vector[1] = int(yx_vector[1])
        xyz_vector = np.insert(yx_vector, 0, 0)

        xyz_final = xyz_vector + initial_centroid
        end_coord = tuple(xyz_final)

        # Calculate the original distance between the two points
        dist_original = math.sqrt((end_coord[1] - starting_coord[1]) ** 2 + (end_coord[2] - starting_coord[2]) ** 2)

        # Reduce the line by 30%
        dist_new = dist_original * 0.7

        # Calculate the new final point coordinates
        x_new = starting_coord[2] + (end_coord[2] - starting_coord[2]) * dist_new / dist_original
        y_new = starting_coord[1] + (end_coord[1] - starting_coord[1]) * dist_new / dist_original

        new_end_coord = tuple([end_coord[0], y_new, x_new])

        z_lin, y_lin, x_lin = np.linspace(starting_coord[0], new_end_coord[0], num=60), np.linspace(starting_coord[1], new_end_coord[1], num=60), np.linspace(starting_coord[2], new_end_coord[2], num=60)

        intensity_profile = [(image[int(round(z_lin[i])), int(round(y_lin[i])), int(round(x_lin[i]))]) for i in range(60)]
        max_intensity = max(intensity_profile)
        normalized_profile = [value / max_intensity for value in intensity_profile]

        fig, ax = plt.subplots(nrows=1, ncols=2)
        fig.suptitle(f'Resolution CT Section {n+1} (D={list_diameters[n]}mm)')
        ax[0].plot(normalized_profile)
        ax[0].set_title(f"Intensity Profile (From Image Center to Edges)")
        ax[0].set_ylabel("Normalized Intensity")

        ax[1].imshow(image[starting_coord[0], :, :], cmap='gray')
        ax[1].imshow(th_h[starting_coord[0], :, :], cmap='Blues', alpha=0.3)
        ax[1].plot([starting_coord[2], new_end_coord[2]], [starting_coord[1], new_end_coord[1]], color='red', linewidth=1)
        ax[1].set_title(f"Profile Line")

        fig.savefig(os.path.join(output, 'ResolutionCT_IntensityProfile_Section' + str(n+1) + '.png'))
        plt.close(fig)

    return contrast_total, holes_total, rc_holes_total, section_names


# Resolution analysis for MR T2 images
def resol_mr_t2(im, coord, vertices, output):
    image = sitk.GetArrayFromImage(im)
    shape = image.shape

    # Select all holes of the resolution insert
    otsu_filter = sitk.OtsuThresholdImageFilter()
    otsu_filter.SetInsideValue(0)
    otsu_filter.SetOutsideValue(1)
    otsu_result = otsu_filter.Execute(im[:, :, coord])
    otsu_threshold = otsu_filter.GetThreshold()

    thresh_water = np.where((image >= otsu_threshold), 1, 0)
    thresh_holes = np.where((image < otsu_threshold), 1, 0)

    labeled_array, num_features = ndimage.label(thresh_water[coord, :, :])
    sizes = ndimage.sum(thresh_water[coord, :, :], labeled_array, range(num_features + 1))
    largest_label = np.argmax(sizes[1:]) + 1
    water_region = np.zeros_like(thresh_water[coord, :, :])
    water_region[labeled_array == largest_label] = 1

    th_water_filled_slice = ndimage.binary_fill_holes(water_region)

    th_holes_slice = thresh_holes[coord, :, :] & th_water_filled_slice

    th_holes = np.zeros_like(image)
    th_water = np.zeros_like(image)
    th_holes[coord, :, :] = th_holes_slice
    th_water[coord, :, :] = water_region

    # Creates 6 triangular segmentations and analyses the resolution in each section #
    # Defines vertices coordinates of each triangular segmentation in the xy plane (z=coord)
    vertices6 = np.array(vertices[5])
    vertices5 = np.array(vertices[4])
    vertices4 = np.array(vertices[3])
    vertices3 = np.array(vertices[2])
    vertices2 = np.array(vertices[1])
    vertices1 = np.array(vertices[0])
    list_vertices = [vertices1, vertices2, vertices3, vertices4, vertices5, vertices6]
    list_diameters = [5, 7.5, 9, 11, 12, 15]

    contrast_total = []
    section_names = []
    holes_total = []
    rc_holes_total = []
    theoretical_holes_total = [55, 31, 20, 15, 10, 8]
    for n in range(len(list_vertices)):
        tri_segm = np.zeros_like(image)
        vertices = list_vertices[n]

        # Obtains vertices (y,x) coordinates
        A, B, C = vertices
        y1, x1 = A
        y2, x2 = B
        y3, x3 = C

        # Calculate original triangle area
        area = abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))

        # Iterate over all coordinates in the image
        for i in range(shape[0]):
            if i == coord:
                for j in range(shape[1]):
                    for k in range(shape[2]):
                        # Calculate the area of each of the triangles formed by the point and each of the sides of the original triangle
                        sub_area1 = abs((x1 - k) * (y2 - j) - (x2 - k) * (y1 - j))
                        sub_area2 = abs((x2 - k) * (y3 - j) - (x3 - k) * (y2 - j))
                        sub_area3 = abs((x3 - k) * (y1 - j) - (x1 - k) * (y3 - j))

                        # Check whether the sum of the areas of the triangles formed by the point and the sides of the original triangle is equal to the area of the original triangle
                        if abs(sub_area1 + sub_area2 + sub_area3 - area) < 1e-6:
                            tri_segm[i, j, k] = 1

        # Analyse the resolution in triangular section #
        segme = np.copy(tri_segm).astype(bool)
        image_copy = np.copy(image).astype(float)

        th_h = np.copy(th_holes).astype(bool)  # holes' segmentation
        th_h = th_h & segme
        th_w = np.copy(th_water).astype(bool)   # water segmentation
        th_w = th_w & segme

        # Contrast calculation
        c_background = np.nanmean(image_copy[th_w])
        c_holes = np.nanmean(image_copy[th_h])
        contrast = c_holes - c_background
        contrast = abs(contrast)
        if contrast > 10:
            contrast = round(contrast)
        else:
            contrast = round(contrast, 1)
        contrast_total.append(contrast)
        section_names.append(str(n + 1))

        # Separate each hole of the section for INTENSITY PROFILE creation
        labeled_array, num_features = measure.label(th_h, background=0, return_num=True, connectivity=None)
        holes_total.append(num_features)
        rc_holes_total.append(num_features / theoretical_holes_total[n])

        dist_total = []
        centroid_total = []
        for i in range(1, num_features + 1):
            single_hole = np.zeros_like(th_h)
            single_hole[labeled_array == i] = 1

            z_coords, y_coords, x_coords = np.where(single_hole == 1)
            centroid_z, centroid_y, centroid_x = np.mean(z_coords), np.mean(y_coords), np.mean(x_coords)
            centroid = int(centroid_z), int(centroid_y), int(centroid_x)
            dist = np.sqrt((centroid_x - x1) ** 2 + (centroid_y - y1) ** 2)
            dist_total.append(dist)
            centroid_total.append(centroid)

        min_value_distance = min(dist_total)
        min_index = dist_total.index(min_value_distance)
        initial_centroid = centroid_total[min_index]  # we choose the hole closest to vertex A
        starting_coord = initial_centroid

        # Create intensity profile
        yx_vector = B - A
        yx_vector[0] = int(yx_vector[0])
        yx_vector[1] = int(yx_vector[1])
        xyz_vector = np.insert(yx_vector, 0, 0)

        xyz_final = xyz_vector + initial_centroid
        end_coord = tuple(xyz_final)

        # Calculate the original distance between the two points
        dist_original = math.sqrt((end_coord[1] - starting_coord[1]) ** 2 + (end_coord[2] - starting_coord[2]) ** 2)

        # Reduce the line by 30%
        dist_new = dist_original * 0.7

        # Calculate the new final point coordinates
        x_new = starting_coord[2] + (end_coord[2] - starting_coord[2]) * dist_new / dist_original
        y_new = starting_coord[1] + (end_coord[1] - starting_coord[1]) * dist_new / dist_original

        new_end_coord = tuple([end_coord[0], y_new, x_new])

        z_lin, y_lin, x_lin = np.linspace(starting_coord[0], new_end_coord[0], num=60), np.linspace(starting_coord[1], new_end_coord[1], num=60), np.linspace(starting_coord[2], new_end_coord[2], num=60)

        intensity_profile = [(image[int(round(z_lin[i])), int(round(y_lin[i])), int(round(x_lin[i]))]) for i in range(60)]
        max_intensity = max(intensity_profile)
        normalized_profile = [value / max_intensity for value in intensity_profile]

        fig, ax = plt.subplots(nrows=1, ncols=2)
        fig.suptitle(f'Resolution MRT2 Section {n+1} (D={list_diameters[n]}mm)')
        ax[0].plot(normalized_profile)
        ax[0].set_title(f"Intensity Profile (From Image Center to Edges)")
        ax[0].set_ylabel("Normalized Intensity")

        ax[1].imshow(image[starting_coord[0], :, :], cmap='gray')
        ax[1].imshow(th_h[starting_coord[0], :, :], cmap='Blues', alpha=0.3)
        ax[1].plot([starting_coord[2], new_end_coord[2]], [starting_coord[1], new_end_coord[1]], color='red', linewidth=1)
        ax[1].set_title(f"Profile Line")

        fig.savefig(os.path.join(output, 'ResolutionMRT2_IntensityProfile_Section' + str(n + 1) + '.png'))
        plt.close(fig)

    return contrast_total, holes_total, rc_holes_total, section_names


# Resolution analysis for MR T1 images
def resol_mr_t1(im, coord, vertices, output, voxeldim):
    image = sitk.GetArrayFromImage(im)
    shape = image.shape

    # Obtain phantom's centre in xy plane
    otsu_filter = sitk.OtsuThresholdImageFilter()
    otsu_filter.SetInsideValue(0)
    otsu_filter.SetOutsideValue(1)
    otsu_result = otsu_filter.Execute(im[:, :, coord])
    otsu_threshold = otsu_filter.GetThreshold()
    thresh = np.where(image[coord, :, :] >= otsu_threshold, 1, 0)
    labeled_array, num_features = ndimage.label(thresh)
    sizes = ndimage.sum(thresh, labeled_array, range(num_features + 1))
    largest_label = np.argmax(sizes[1:]) + 1
    water_region = np.zeros_like(thresh)
    water_region[labeled_array == largest_label] = 1
    panal = ndimage.binary_fill_holes(water_region).astype(np.uint8)
    center_phantom = ndimage.measurements.center_of_mass(panal)     # phantom's centre (y,x)

    z, y, x = np.indices(shape)
    y = y - center_phantom[0]
    x = x - center_phantom[1]

    distance_to_center = np.sqrt((x*voxeldim[0]) ** 2 + (y*voxeldim[1]) ** 2)

    # Creates 6 triangular segmentations and analyses the resolution in each section #
    # Defines vertices coordinates of each triangular segmentation in the xy plane (z=coord)
    vertices6 = np.array(vertices[5])
    vertices5 = np.array(vertices[4])
    vertices4 = np.array(vertices[3])
    vertices3 = np.array(vertices[2])
    vertices2 = np.array(vertices[1])
    vertices1 = np.array(vertices[0])
    list_vertices = [vertices1, vertices2, vertices3, vertices4, vertices5, vertices6]
    list_diameters = [5, 7.5, 9, 11, 12, 15]

    # Generates each triangular section
    contrast_total = []
    section_names = []
    holes_total = []
    rc_holes_total = []
    theoretical_holes_total = [55, 31, 20, 15, 10, 8]
    for n in range(len(list_vertices)):
        tri_segm = np.zeros_like(image)
        vertices = list_vertices[n]

        # Obtains vertices (y,x) coordinates
        A, B, C = vertices
        y1, x1 = A
        y2, x2 = B
        y3, x3 = C

        # Calculate original triangle area
        area = abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))

        # Iterate over all coordinates in the image
        for i in range(shape[0]):
            if i == coord:
                for j in range(shape[1]):
                    for k in range(shape[2]):
                        # Calculate the area of each of the triangles formed by the point and each of the sides of the original triangle
                        sub_area1 = abs((x1 - k) * (y2 - j) - (x2 - k) * (y1 - j))
                        sub_area2 = abs((x2 - k) * (y3 - j) - (x3 - k) * (y2 - j))
                        sub_area3 = abs((x3 - k) * (y1 - j) - (x1 - k) * (y3 - j))

                        # Check whether the sum of the areas of the triangles formed by the point and the sides of the original triangle is equal to the area of the original triangle
                        if abs(sub_area1 + sub_area2 + sub_area3 - area) < 1e-6:
                            tri_segm[i, j, k] = 1

        # Analyse the resolution in triangular section #
        segme = np.copy(tri_segm).astype(bool)
        image_copy = np.copy(image).astype(float)
        image_copy[~segme] = np.nan  # image only in triangular section

        # Analyse resolution in each disk #
        disk_rad = 122.5    # max disk radius in mm
        holes_tri = []
        water_tri = []
        for i in range(5):
            disk_rad_min = disk_rad - 24.5
            disk_rad_max = disk_rad
            mask = (distance_to_center <= disk_rad_max) & (distance_to_center >= disk_rad_min)

            image_disk = np.copy(image_copy)
            image_disk[~mask] = np.nan
            max_value = np.nanmax(image_disk)  # max intensity value in triangular and disk section
            th40 = 0.7 * max_value  # 68% threshold of the max intensity value

            holes_disk = np.where(image_disk >= th40, 1, 0).astype(int)
            water_disk = np.where(image_disk < th40, 1, 0).astype(int)
            holes = ndimage.binary_fill_holes(holes_disk)
            holes_tri.append(holes)
            water_tri.append(water_disk)

            if i != 3:
                disk_rad -= 24.5
            elif i == 3:
                disk_rad -= 18

        sum_holes_tri = np.sum(holes_tri, axis=0)
        sum_holes_tri = np.where(sum_holes_tri > 1, 1, sum_holes_tri)   # holes' segmentation

        sum_water_tri = np.sum(water_tri, axis=0)
        sum_water_tri = np.where(sum_water_tri > 1, 1, sum_water_tri)
        sum_water_tri = sum_water_tri & panal   # water's segmentation

        th_h = np.copy(sum_holes_tri).astype(bool)
        th_w = np.copy(sum_water_tri).astype(bool)

        # Contrast calculation
        c_background = np.nanmean(image_copy[th_w])
        c_holes = np.nanmean(image_copy[th_h])
        contrast = c_holes - c_background
        contrast = abs(contrast)
        if contrast > 10:
            contrast = round(contrast)
        else:
            contrast = round(contrast, 1)
        contrast_total.append(contrast)
        section_names.append(str(n + 1))

        # Separate each hole of the section for INTENSITY PROFILE creation
        labeled_array, num_labels = measure.label(th_h, background=0, return_num=True, connectivity=None)
        holes_total.append(num_labels)
        rc_holes_total.append(num_features / theoretical_holes_total[n])

        dist_total = []
        centroid_total = []
        for i in range(1, num_labels + 1):
            single_hole = np.zeros_like(th_h)
            single_hole[labeled_array == i] = 1

            z_coords, y_coords, x_coords = np.where(single_hole == 1)
            centroid_z, centroid_y, centroid_x = np.mean(z_coords), np.mean(y_coords), np.mean(x_coords)
            centroid = int(centroid_z), int(centroid_y), int(centroid_x)
            dist = np.sqrt((centroid_x - x1) ** 2 + (centroid_y - y1) ** 2)
            dist_total.append(dist)
            centroid_total.append(centroid)

        min_value_distance = min(dist_total)
        min_index = dist_total.index(min_value_distance)
        initial_centroid = centroid_total[min_index]  # we choose the hole closest to vertex A
        starting_coord = initial_centroid

        # Create intensity profile
        yx_vector = B - A
        yx_vector[0] = int(yx_vector[0])
        yx_vector[1] = int(yx_vector[1])
        xyz_vector = np.insert(yx_vector, 0, 0)

        xyz_final = xyz_vector + initial_centroid
        end_coord = tuple(xyz_final)

        # Calculate the original distance between the two points
        dist_original = math.sqrt((end_coord[1] - starting_coord[1]) ** 2 + (end_coord[2] - starting_coord[2]) ** 2)

        # Reduce the line by 30%
        dist_new = dist_original * 0.7

        # Calculate the new final point coordinates
        x_new = starting_coord[2] + (end_coord[2] - starting_coord[2]) * dist_new / dist_original
        y_new = starting_coord[1] + (end_coord[1] - starting_coord[1]) * dist_new / dist_original

        new_end_coord = tuple([end_coord[0], y_new, x_new])

        z_lin, y_lin, x_lin = np.linspace(starting_coord[0], new_end_coord[0], num=60), np.linspace(starting_coord[1], new_end_coord[1], num=60), np.linspace(starting_coord[2], new_end_coord[2], num=60)

        intensity_profile = [(image[int(round(z_lin[i])), int(round(y_lin[i])), int(round(x_lin[i]))]) for i in range(60)]
        max_intensity = max(intensity_profile)
        normalized_profile = [value / max_intensity for value in intensity_profile]

        fig, ax = plt.subplots(nrows=1, ncols=2)
        fig.suptitle(f'Resolution MRT1 Section {n+1} (D={list_diameters[n]}mm)')
        ax[0].plot(normalized_profile)
        ax[0].set_title(f"Intensity Profile (From Image Center to Edges)")
        ax[0].set_ylabel("Normalized Intensity")

        ax[1].imshow(image[starting_coord[0], :, :], cmap='gray')
        ax[1].imshow(th_h[starting_coord[0], :, :], cmap='Blues', alpha=0.3)
        ax[1].plot([starting_coord[2], new_end_coord[2]], [starting_coord[1], new_end_coord[1]], color='red', linewidth=1)
        ax[1].set_title(f"Profile Line")

        fig.savefig(os.path.join(output, 'ResolutionMRT1_IntensityProfile_Section' + str(n + 1) + '.png'))
        plt.close(fig)

    return contrast_total, holes_total, rc_holes_total, section_names


# CO-REGISTRATION analysis #
# PET threshold to select the co-registration tube (only)
def pet_th(im, voxeldim):
    image = sitk.GetArrayFromImage(im)
    max_value = np.nanmax(image[35:, :, :])                     # max intensity value
    max_value70 = np.nanmean(image[image >= 0.7 * max_value])   # mean of the voxels with value >= 0.7max_value
    th40 = 0.4 * max_value70                                    # 40% of the max_value70

    # Create the segmentation
    thresh = np.where((image >= th40), 1, 0)
    thresh_pet_sitk = sitk.GetImageFromArray(thresh)
    thresh_pet_sitk.SetDirection(im.GetDirection())
    thresh_pet_sitk.SetSpacing(im.GetSpacing())
    thresh_pet_sitk.SetOrigin(im.GetOrigin())

    # Separate in individual volumes
    labeled_array, num_features = ndimage.label(thresh)
    centroids = ndimage.measurements.center_of_mass(thresh, labeled_array, range(1, num_features + 1))
    center_image = np.array(image.shape) // 2
    distances = [np.linalg.norm(np.array(center_image) - centroid) for centroid in centroids]
    closest_centroid_index = np.argmin(distances)
    chosen_segmentation = np.zeros_like(thresh)
    chosen_segmentation[labeled_array == closest_centroid_index + 1] = 1

    numvoxels = np.sum(chosen_segmentation)
    vol = numvoxels * (voxeldim / 1000)

    if vol > 150:
        voi_pet_sitk = sitk.GetImageFromArray(chosen_segmentation)
        voi_pet_sitk.SetDirection(im.GetDirection())
        voi_pet_sitk.SetSpacing(im.GetSpacing())
        voi_pet_sitk.SetOrigin(im.GetOrigin())

        ref_voi_pet = ndimage.binary_dilation(chosen_segmentation, iterations=4).astype(np.uint8)
        ref_voi_pet_sitk = sitk.GetImageFromArray(ref_voi_pet)
        ref_voi_pet_sitk.SetDirection(im.GetDirection())
        ref_voi_pet_sitk.SetSpacing(im.GetSpacing())
        ref_voi_pet_sitk.SetOrigin(im.GetOrigin())

    return voi_pet_sitk, ref_voi_pet_sitk


# CT and MR threshold in all image to select the co-registration tube (also selects other regions)
def ct_mr_th(im, system_type):
    image = sitk.GetArrayFromImage(im)

    # PET/CT system: CT threshold
    if "PET/CT" in system_type:
        threshold = np.where((image >= -40) & (image <= 30), 1, 0)
    # PET/MR system: MR threshold
    elif "PET/MR" in system_type:
        # MR T1 sequences' threshold
        if "T2" in system_type:
            vmax = np.nanmax(image)
            value_th = 0.055
            threshold = np.where((image >= value_th * vmax), 1, 0)
        # MR T2 sequences' threshold
        elif "T1" in system_type:
            vmax = np.nanmax(image)
            value_th = 0.055
            threshold = np.where((image >= value_th * vmax), 1, 0)

    seg = ndimage.binary_erosion(threshold)
    seg = ndimage.binary_dilation(seg)
    seg = ndimage.binary_closing(seg, structure=np.ones((3, 3, 3), dtype=int))
    seg = ndimage.binary_fill_holes(seg, structure=np.ones((3, 3, 3), dtype=int))
    seg = np.copy(seg).astype(int)

    seg_sitk = sitk.GetImageFromArray(seg)
    seg_sitk.SetDirection(im.GetDirection())
    seg_sitk.SetSpacing(im.GetSpacing())
    seg_sitk.SetOrigin(im.GetOrigin())

    return seg_sitk


# DISTORTION analysis #
# Get phantom's center in xy plane
def phantom_center(im, panal_slice, slice_spacing, modality):
    image = sitk.GetArrayFromImage(im)

    slice_cn = int(round(panal_slice - (15.6 / slice_spacing)))
    slice = image[slice_cn, :, :]

    if modality == 'CT':
        thresh = np.where(slice >= -644, 1, 0)
    elif modality == 'MR':
        otsu_filter = sitk.OtsuThresholdImageFilter()
        otsu_filter.SetInsideValue(0)
        otsu_filter.SetOutsideValue(1)
        otsu_result = otsu_filter.Execute(im[:, :, slice_cn])
        otsu_threshold = otsu_filter.GetThreshold()
        thresh = np.where(slice >= otsu_threshold, 1, 0)

    labeled_array, num_features = ndimage.label(thresh)
    sizes = ndimage.sum(thresh, labeled_array, range(num_features + 1))
    largest_label = np.argmax(sizes[1:]) + 1
    water_region = np.zeros_like(thresh)
    water_region[labeled_array == largest_label] = 1
    panal = ndimage.binary_fill_holes(water_region).astype(np.uint8)

    center_phantom = ndimage.measurements.center_of_mass(panal)

    return center_phantom


# Segmentation of cylindrical distortion insert: all modalities (CT, MR T1 and MR T2)
def th_cylinsert(im, in_sl, fin_sl, voxeldim, center_phantom, filled_cyl, modality):
    image = sitk.GetArrayFromImage(im)

    shape = image.shape
    z, y, x = np.indices(shape)

    if modality == 'CT':
        # Cylinder insert center calculation
        center_insert = [center_phantom[0] - (2 / voxeldim[1]), center_phantom[1] - (76 / voxeldim[0])]
        y = y - center_insert[0]
        x = x - center_insert[1]

        # Reduces threshold area
        max_radius = 24  # max radius in mm from cylinder's centre
        distance_to_center_mm = np.sqrt((x*voxeldim[0])**2 + (y*voxeldim[1])**2)
        mask = (distance_to_center_mm <= max_radius) & (z < fin_sl) & (in_sl < z)
        image_copy = np.copy(image).astype(float)
        image_copy[~mask] = np.nan
        if filled_cyl == 'no':
            threshold_str = np.where(image_copy >= -700, 1, 0)
        elif filled_cyl == 'yes':
            threshold_str = np.where(image_copy >= 14, 1, 0)

    elif "T2" in modality:
        # Cylinder insert center calculation
        center_insert = [center_phantom[0] - (2 / voxeldim[1]), center_phantom[1] - (71 / voxeldim[0])]
        y = y - center_insert[0]
        x = x - center_insert[1]

        # Reduces threshold area
        max_radius = 30     # max radius in mm from cylinder's centre
        distance_to_center_mm = np.sqrt((x*voxeldim[0])**2 + (y*voxeldim[1])**2)
        mask = (distance_to_center_mm <= max_radius) & (z < fin_sl) & (in_sl < z)
        image_copy = np.copy(image).astype(float)
        image_copy[~mask] = np.nan
        if filled_cyl == 'yes':
            otsu_filter = sitk.OtsuThresholdImageFilter()
            otsu_filter.SetInsideValue(0)
            otsu_filter.SetOutsideValue(1)
            otsu_result = otsu_filter.Execute(im[:, :, in_sl:fin_sl])
            otsu_threshold = otsu_filter.GetThreshold()
            threshold_mesh = np.where((image < otsu_threshold), 1, 0)
            threshold_water = np.where((image_copy >= otsu_threshold), 1, 0)
            closed_water = ndimage.binary_closing(threshold_water)
            threshold_str = threshold_mesh & closed_water

    elif "T1" in modality:
        # Cylinder insert center calculation
        center_insert = [center_phantom[0] - (2 / voxeldim[1]), center_phantom[1] - (71 / voxeldim[0])]
        y = y - center_insert[0]
        x = x - center_insert[1]

        # Reduces threshold area
        max_radius = 30     # max radius in mm from cylinder's centre
        distance_to_center_mm = np.sqrt((x*voxeldim[0])**2 + (y*voxeldim[1])**2)
        mask = (distance_to_center_mm <= max_radius) & (z < fin_sl) & (in_sl < z)
        if filled_cyl == 'yes':
            otsu_filter = sitk.OtsuMultipleThresholds(im, 3, 0, 128, False, True)
            otsu_array = sitk.GetArrayFromImage(otsu_filter)

            threshold_str = otsu_array == 1
            threshold_str = np.copy(threshold_str).astype(int)
            threshold_str = threshold_str & np.copy(mask).astype(int)


    segmentation_sitk = sitk.GetImageFromArray(threshold_str)
    segmentation_sitk.SetDirection(im.GetDirection())
    segmentation_sitk.SetSpacing(im.GetSpacing())
    segmentation_sitk.SetOrigin(im.GetOrigin())

    return segmentation_sitk


# CT segmentation of panal distortion insert
def ct_th_panalinsert(im, panal_slice, center_phantom, voxeldim):
    image = sitk.GetArrayFromImage(im)
    shape = image.shape

    threshold = np.where(image >= -590, 1, 0)

    z, y, x = np.indices(shape)
    y = y - center_phantom[0]
    x = x - center_phantom[1]
    distance_to_center_mm = np.sqrt((x*voxeldim[0])**2 + (y*voxeldim[1])**2)
    max_radius = 115  # max radius in mm
    mask = (distance_to_center_mm <= max_radius) & (z == panal_slice)

    new_segm = np.zeros_like(threshold, dtype=np.float)
    segmentation = np.copy(threshold).astype('bool')
    new_segm[mask] = 1
    new_segm[~segmentation] = 0

    segmentation_sitk = sitk.GetImageFromArray(new_segm)
    segmentation_sitk.SetDirection(im.GetDirection())
    segmentation_sitk.SetSpacing(im.GetSpacing())
    segmentation_sitk.SetOrigin(im.GetOrigin())

    return segmentation_sitk


# Choose xy points in a 3D image slice: returns [z,y,x] coordinates of the points
def choose_points(dicom_path, coord, image_name, direction):
    # Lectura de la imagen
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(dicom_path)
    reader.SetFileNames(dicom_names)
    im = reader.Execute()
    image = sitk.GetArrayFromImage(im)

    axial_slice = image[coord, :, :] + image[coord-1, :, :]

    fig, ax = plt.subplots()
    ax.imshow(axial_slice, cmap='gray', origin='lower')
    ax.set_title(f'Select Panal Insert Points in \n {image_name} image {direction} direction ')

    ax.invert_yaxis()

    coordinates = []

    def on_image_click(event):
        if event.inaxes == ax:
            x = int(event.xdata)
            y = int(event.ydata)
            coordinates.append([y, x])
            ax.plot(x, y, 'ro', markersize=3)
            plt.draw()

    # Conectar la función de clic a la figura
    fig.canvas.mpl_connect('button_press_event', on_image_click)
    fig.canvas.set_window_title('Select Points of Panal Insert for Distortion Analysis')

    # Mostrar la figura y esperar hasta que la ventana sea cerrada
    plt.show(block=True)

    return coordinates


# Calculates the difference between contiguous elements of a list
def calculate_differences(list):
    differences = []
    for i in range(len(list) - 1):
        difference = list[i+1] - list[i]
        differences.append(difference)
    return differences


# Function that returns panal walls coordinates
def panal_walls(pixel_list, coord_list):
    processed_list = []
    index_list = []
    start = None
    end = None

    for i, pixel_value in enumerate(pixel_list):
        if pixel_value == 1:
            if start is None:
                start = i
                end = i
            else:
                end = i
        else:
            if start is not None:
                index_list.append((start, end))
                start = None
                end = None

    # Calculate the averages for each sequence of 1s
    for start, end in index_list:
        coords_sum = [0, 0]  # Accumulator for coordinates sum
        count = 0  # Counter for coordinates in this sequence
        for i in range(start, end + 1):
            coords_sum[0] += coord_list[i][0]  # Adding y-coordinate
            coords_sum[1] += coord_list[i][1]  # Adding x-coordinate
            count += 1
        avg_coords = [coords_sum[0] / count, coords_sum[1] / count] if count > 0 else 0
        processed_list.append(avg_coords)

    return processed_list


# Calculates distances in mm between two consecutive coordinate pairs
def calculate_mm_distances(coord_list, y_spacing, x_spacing):
    distances_mm = []
    for i in range(1, len(coord_list)):
        diff_y = coord_list[i][0] - coord_list[i - 1][0]
        diff_x = coord_list[i][1] - coord_list[i - 1][1]

        distance = math.sqrt((diff_y * y_spacing) ** 2 + (diff_x * x_spacing) ** 2)
        distances_mm.append(distance)

    return distances_mm



# RADIOMICS analysis #
# Resample an image to fixed image dimensions
def resampling(fixed_image, moving_image):
    dimension = 3
    identity = sitk.Transform(dimension, sitk.sitkIdentity)
    moving_resampled = sitk.Resample(moving_image, fixed_image, identity, sitk.sitkNearestNeighbor, 0.0, moving_image.GetPixelID())

    return moving_resampled


# Register an image to a fixed/reference one
def register(fixed_image, moving_image):
    """
    Registers a moving image to a fixed image using SimpleITK.

    Parameters:
    - fixed_image: SimpleITK Image, the fixed image.
    - moving_image: SimpleITK Image, the moving image.

    Returns:
    - final_transform: SimpleITK Transform, the final transformation after registration.
    """

    initial_transform = sitk.CenteredTransformInitializer(fixed_image, moving_image, sitk.Euler3DTransform(), sitk.CenteredTransformInitializerFilter.MOMENTS)

    registration_method = sitk.ImageRegistrationMethod()

    # Similarity metric settings.
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=100)
    registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
    registration_method.SetMetricSamplingPercentage(0.01)

    registration_method.SetInterpolator(sitk.sitkLinear)

    # Optimizer settings.
    registration_method.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=300,
                                                      convergenceMinimumValue=1e-6, convergenceWindowSize=10)
    registration_method.SetOptimizerScalesFromPhysicalShift()

    # Setup for the multi-resolution framework.
    registration_method.SetShrinkFactorsPerLevel(shrinkFactors=[4, 2, 1])
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[2, 1, 0])
    registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    # Set the initial moving and optimized transforms.
    optimized_transform = sitk.TranslationTransform(fixed_image.GetDimension())
    registration_method.SetMovingInitialTransform(initial_transform)
    registration_method.SetInitialTransform(optimized_transform, inPlace=False)

    # Need to compose the transformations after registration.
    final_transform = sitk.CompositeTransform(
        [registration_method.Execute(fixed_image, moving_image), initial_transform])

    return final_transform


# Percent Inactive feature
def getPercentInactive(path_image, path_segm, thresh=0.6):
    raw, meta = nrrd.read(path_image)
    mask, metamask = nrrd.read(path_segm)
    mask = mask.astype(bool).astype('uint8')
    mask_bool = np.array(mask, dtype=bool)
    mask_bool = mask_bool.astype('uint8')

    ROI = raw.copy()
    ROI[mask_bool == 0] = 0
    ROI[mask_bool != 0] = raw[mask_bool != 0]

    maskInactive = np.zeros(np.shape(ROI))
    maskInactive[ROI > (thresh * np.max(ROI) ** 2)] = 1

    perim = bwperim(mask, 8)
    newMask = maskInactive + perim
    newMask[mask == 0] = 10
    newMask[newMask == 1] = 10
    newMask[newMask == 2] = 10
    newMask[newMask == 0] = 1
    newMask[newMask == 10] = 0

    all_labels = measure.label(newMask)
    numLabels = np.shape(np.unique(all_labels))[0]

    b = np.zeros([1, numLabels - 1])
    for l in np.arange(numLabels - 1):
        if np.shape(np.where(all_labels == l))[1] < 15:
            b[0, l] = 0
        else:
            b[0, l] = 1

    [row, col] = np.where(b > 0)

    sumInactive = 0
    for i in np.arange(np.size(col)):
        sumInactive = sumInactive + np.shape(np.where(all_labels == (i + 1)))[1]

    sumVolume = np.sum(mask)
    percentInactive = sumInactive / sumVolume * 100

    return percentInactive


# Eccentricity feature
def getEccentricity(segm_path):
    segm = sitk.ReadImage(segm_path)
    mask = sitk.GetArrayFromImage(segm)

    moments = np.array(regionprops(mask.astype(int))[0].inertia_tensor)

    eigenvalues = np.linalg.eigvals(moments)
    eigenvalues.sort()

    # Calcular la excentricidad
    eccentricity = np.sqrt(1 - (eigenvalues[0] ** 2) / (eigenvalues[2] * eigenvalues[1]))

    return eccentricity


def update_params_file(bin_width, params_file_path):
    with open(params_file_path, 'r') as file:
        params = yaml.safe_load(file)

    params['setting']['binWidth'] = bin_width

    with open(params_file_path, 'w') as file:
        yaml.dump(params, file)


# Obtains radiomics features
def get_radiomics(image_path, segm_path, modality, bin_width, params_path):
    if modality == 'PET':
        update_params_file(bin_width, params_path)
        extractor = featureextractor.RadiomicsFeatureExtractor(params_path)
    else:
        extractor = featureextractor.RadiomicsFeatureExtractor()

    result = extractor.execute(image_path, segm_path)

    df = pd.DataFrame.from_dict(dict(result), orient='index')
    df.transpose()

    mean = df.loc['original_firstorder_Mean'][0]
    median = df.loc['original_firstorder_Median'][0]
    variance = df.loc['original_firstorder_Variance'][0]
    skewness = df.loc['original_firstorder_Skewness'][0]
    kurtosis = df.loc['original_firstorder_Kurtosis'][0]
    entropy = df.loc['original_firstorder_Entropy'][0]
    energy = df.loc['original_firstorder_Energy'][0]
    percentInactive = getPercentInactive(image_path, segm_path)
    eccentricity = getEccentricity(segm_path)

    glcmEnergy = df.loc['original_glcm_JointEnergy'][0]
    glcmContrast = df.loc['original_glcm_Contrast'][0]
    glcmEntropy = df.loc['original_glcm_JointEnergy'][0]
    glcmHomogeinity = df.loc['original_glcm_Idm'][0]
    glcmCorrelation = df.loc['original_glcm_Correlation'][0]
    glcmVariance = df.loc['original_glcm_DifferenceVariance'][0]
    glcmDissimilarity = df.loc['original_glcm_DifferenceAverage'][0]
    glcmAutocorrelation = df.loc['original_glcm_Autocorrelation'][0]

    glszmSAE = df.loc['original_glszm_SmallAreaEmphasis'][0]
    glszmLAE = df.loc['original_glszm_LargeAreaEmphasis'][0]
    glszmGLN = df.loc['original_glszm_GrayLevelNonUniformity'][0]
    glszmSZN = df.loc['original_glszm_SizeZoneNonUniformity'][0]
    glszmZP = df.loc['original_glszm_ZonePercentage'][0]
    glszmLGLZE = df.loc['original_glszm_LowGrayLevelZoneEmphasis'][0]
    glszmHGLZE = df.loc['original_glszm_HighGrayLevelZoneEmphasis'][0]
    glszmSALGLE = df.loc['original_glszm_SmallAreaLowGrayLevelEmphasis'][0]
    glszmSAHGLE = df.loc['original_glszm_SmallAreaHighGrayLevelEmphasis'][0]
    glszmLALGLE = df.loc['original_glszm_LargeAreaLowGrayLevelEmphasis'][0]
    glszmLAHGLE = df.loc['original_glszm_LargeAreaHighGrayLevelEmphasis'][0]
    glszmGLV = df.loc['original_glszm_GrayLevelVariance'][0]
    glszmZV = df.loc['original_glszm_ZoneVariance'][0]

    glrlmSRE = df.loc['original_glrlm_ShortRunEmphasis'][0]
    glrlmLRE = df.loc['original_glrlm_LongRunEmphasis'][0]
    glrlmGLN = df.loc['original_glrlm_GrayLevelNonUniformity'][0]
    glrlmRLN = df.loc['original_glrlm_RunLengthNonUniformity'][0]
    glrlmRP = df.loc['original_glrlm_RunPercentage'][0]
    glrlmLGRE = df.loc['original_glrlm_LowGrayLevelRunEmphasis'][0]
    glrlmHGRE = df.loc['original_glrlm_HighGrayLevelRunEmphasis'][0]
    glrlmSRLGLE = df.loc['original_glrlm_ShortRunLowGrayLevelEmphasis'][0]
    glrlmSRHGLE = df.loc['original_glrlm_ShortRunHighGrayLevelEmphasis'][0]
    glrlmLRLGLE = df.loc['original_glrlm_LongRunLowGrayLevelEmphasis'][0]
    glrlmLRHGLE = df.loc['original_glrlm_LongRunHighGrayLevelEmphasis'][0]
    glrlmGLV = df.loc['original_glrlm_GrayLevelVariance'][0]
    glrlmRLV = df.loc['original_glrlm_RunVariance'][0]

    ngtdmCoarseness = df.loc['original_ngtdm_Coarseness'][0]
    ngtdmContrast = df.loc['original_ngtdm_Contrast'][0]
    ngtdmBusyness = df.loc['original_ngtdm_Busyness'][0]
    ngtdmComplexity = df.loc['original_ngtdm_Complexity'][0]
    ngtdmStrength = df.loc['original_ngtdm_Strength'][0]

    radiomics = {'Mean': mean, 'Median': median, 'Variance': variance, 'Skewness': skewness, 'Kurtosis': kurtosis,
                 'Entropy': entropy,
                 'Energy': energy, 'PercentInactive': percentInactive,
                 'Eccentricity': eccentricity, 'GLCMEnergy': glcmEnergy, 'GLCMContrast': glcmContrast,
                 'GLCMEntropy': glcmEntropy,
                 'GLCMHomogeinity': glcmHomogeinity, 'GLCMCorrelation': glcmCorrelation,
                 'GLCMVariance': glcmVariance, 'GLCMDissimilarity': glcmDissimilarity,
                 'GLCMAutocorrelation': glcmAutocorrelation,
                 'GLSZMSAE': glszmSAE, 'GLSZMLAE': glszmLAE, 'GLSZMGLN': glszmGLN, 'GLSZMSZN': glszmSZN,
                 'GLSZMZP': glszmZP, 'GLSZMLGLZE': glszmLGLZE, 'GLSZMHGLZE': glszmHGLZE,
                 'GLSZMSALGLE': glszmSALGLE,
                 'GLSZMSAHGLE': glszmSAHGLE, 'GLSZMLALGLE': glszmLALGLE, 'GLSZMLAHGLE': glszmLAHGLE,
                 'GLSZMGLV': glszmGLV, 'GLSZMZV': glszmZV, 'GLRLMSRE': glrlmSRE, 'GLRLMLRE': glrlmLRE,
                 'GLRLMGLN': glrlmGLN,
                 'GLRLMRLN': glrlmRLN, 'GLRLMRP': glrlmRP, 'GLRLMLGRE': glrlmLGRE, 'GLRLMHGRE': glrlmHGRE,
                 'GLRLMSRLGLE': glrlmSRLGLE, 'GLRLMSRHGLE': glrlmSRHGLE, 'GLRLMLRLGLE': glrlmLRLGLE,
                 'GLRLMLRHGLE': glrlmLRHGLE,
                 'GLRLMGLV': glrlmGLV, 'GLRLMRLV': glrlmRLV, 'NGTDMCoarseness': ngtdmCoarseness,
                 'NGTDMContrast': ngtdmContrast, 'NGTDMBusyness': ngtdmBusyness,
                 'NGTDMComplexity': ngtdmComplexity, 'NGTDMStrength': ngtdmStrength}

    return radiomics


# Main radiomic features analysis
def radiomics_analysis(path_image, path_image_ref, path_seg_ref, type_image, dir, output_path, params_path):
    # Convert DICOM image serie to NRRD format
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(path_image)
    reader.SetFileNames(dicom_names)
    serie = reader.Execute()
    nrrd_path = os.path.join(output_path, type_image + '_image.nrrd')
    sitk.WriteImage(serie, nrrd_path)

    # Radiomics analysis
    im = sitk.ReadImage(nrrd_path, sitk.sitkFloat32)  # analysis image
    image_array = sitk.GetArrayFromImage(im)
    unique_values = np.unique(image_array)
    data_range = np.max(unique_values) - np.min(unique_values)
    num_bins = 128
    bin_width = int(data_range / num_bins)

    im_ref = sitk.ReadImage(path_image_ref, sitk.sitkFloat32)  # reference image
    registration_transform = register(im, im_ref)

    segm_names = []
    all_resampled_segmentations = []
    all_resampled_segmentations_paths = []
    for seg_ref in path_seg_ref:
        segm_ref = sitk.ReadImage(seg_ref, sitk.sitkFloat32)  # reference segmentation
        name_segm_ref = os.path.splitext(os.path.basename(seg_ref))[0]
        segm_names.append(name_segm_ref)

        moving_seg_resampled = sitk.Resample(segm_ref, im, registration_transform, sitk.sitkLinear, 0.0, segm_ref.GetPixelID())
        # moving_seg_resampled = sitk.Resample(segm_ref, im, sitk.Transform(3, sitk.sitkIdentity), sitk.sitkNearestNeighbor, 0.0, segm_ref.GetPixelID())
        moving_seg_resampled_int = sitk.Cast(moving_seg_resampled, sitk.sitkUInt8)
        seg_resampled_path = os.path.join(output_path, name_segm_ref + str(type_image) + '_registered.nrrd')
        sitk.WriteImage(moving_seg_resampled_int, seg_resampled_path)

        all_resampled_segmentations.append(moving_seg_resampled_int)
        all_resampled_segmentations_paths.append(seg_resampled_path)

    # Visualize image and resampled segmentations
    visualize_slices_3d_with_segmentations(sitk.GetArrayFromImage(im), [sitk.GetArrayFromImage(seg) for seg in all_resampled_segmentations], type_image, dir)

    # Ask user if they want to continue with existing segmentations or define new ones
    user_input = input("Do you want to continue with existing segmentations (C) or define new ones (N)? ")

    if user_input.upper() == "C":
        all_segm_paths = all_resampled_segmentations_paths
    elif user_input.upper() == "N":
        print("Select the folder containing the new segmentations:")
        new_seg_folder = filedialog.askdirectory()
        if new_seg_folder:
            print("Selected folder:", new_seg_folder)
            new_seg_files = os.listdir(new_seg_folder)
            new_seg_paths = [os.path.join(new_seg_folder, file) for file in new_seg_files]
            all_segm_paths = new_seg_paths
        else:
            print("No folder selected.")
    else:
        print("Invalid input. Please enter C to continue with existing segmentations or N to define new ones.")

    # # Normalize PET image to first variability measure
    # original_im = sitk.ReadImage(nrrd_path, sitk.sitkFloat32)
    # original_im_array = sitk.GetArrayFromImage(original_im)
    # multiplier = 1.4211966947241974
    # normalized_im_array = original_im_array * multiplier

    # Create a new normalised image with the normalised matrix and the original position and orientation information
    # normalized_im = sitk.GetImageFromArray(normalized_im_array)
    # normalized_im.SetOrigin(original_im.GetOrigin())
    # normalized_im.SetSpacing(original_im.GetSpacing())
    # normalized_im.SetDirection(original_im.GetDirection())
    # normalized_nrrd_path = os.path.join(output_date_analysis_path, sub_folder + '_image_normalized.nrrd')
    # sitk.WriteImage(normalized_im, normalized_nrrd_path)

    # Get radiomic characteristics
    all_radiomics = []
    for i in range(len(all_segm_paths)):
        segm_radiomics = get_radiomics(nrrd_path, all_segm_paths[i], type_image, bin_width, params_path)
        # segm_radiomics = get_radiomics(normalized_nrrd_path, all_segm_paths[i], type_image, bin_width, params_path)
        all_radiomics.append(segm_radiomics)

    writer = pd.ExcelWriter(os.path.join(output_path, f'Radiomics{str(type_image)}.xlsx'))

    parametersRadiomicsDf = pd.DataFrame()
    for i in np.arange(len(path_seg_ref)):
        df = pd.DataFrame.from_records(all_radiomics[i], index=[segm_names[i]])
        parametersRadiomicsDf = pd.concat([parametersRadiomicsDf, df])
    parametersRadiomicsDf.sort_index(inplace=True)
    parametersRadiomicsDf.index.name = 'Segmentation'
    parametersRadiomicsDf.to_excel(writer, sheet_name='Radiomics')

    writer.save()


def visualize_slices_3d_with_segmentations(image_array, segmentation_arrays, image_name, date):
    slice_index = 0
    axial_slices = [image_array[i, :, :] for i in range(image_array.shape[0])]
    segmentations_on_slices = [seg for seg in segmentation_arrays]

    def update_axial_slice():
        nonlocal slice_index
        ax.clear()
        ax.imshow(axial_slices[slice_index], cmap='gray', origin='lower')
        ax.set_title(f'Axial Slice {slice_index} \n of {image_name} image dated {date}')
        ax.invert_yaxis()

        # Visualizar todas las segmentaciones con la misma intensidad
        for seg in segmentations_on_slices:
            ax.imshow(seg[slice_index] * 0.5, cmap='jet', origin='lower', alpha=0.3)

        fig.canvas.draw_idle()

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)

    update_axial_slice()

    slider_axial = Slider(ax=plt.axes([0.15, 0.1, 0.65, 0.03]), label='Axial Slice', valmin=0, valmax=image_array.shape[0]-1, valinit=slice_index, valstep=1)

    def on_slider_change(val):
        nonlocal slice_index
        slice_index = int(val)
        update_axial_slice()

    slider_axial.on_changed(on_slider_change)
    plt.show(block=True)


# Return an alphabetically sorted names list of NRRD or NIFTI images
def names(path):
    foodir = path
    names = list()
    for root, dirs, files in os.walk(foodir):
        for f in files:
            if os.path.splitext(f)[1].lower() == ".nrrd":
                g = f.replace(".nrrd", "")
                names.append(g)
            elif os.path.splitext(f)[1].lower() == ".nii":
                g = f.replace(".nii", "")
                names.append(g)
    names_sorted = sorted(names)
    return names_sorted


# Return an alphabetically sorted directories list of NRRD or NIFTI images
def directories(path):
    foodir = path
    barlist = list()

    for root, dirs, files in os.walk(foodir):
        for f in files:
            if os.path.splitext(f)[1].lower() == ".nrrd" or os.path.splitext(f)[1].lower() == ".nii":
                barlist.append(os.path.join(root, f))
    barlist_sorted = sorted(barlist)
    return barlist_sorted
