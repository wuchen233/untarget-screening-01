from numpy import genfromtxt
import numpy as np
import csv
import time
from functools import reduce
from statistics import median
import copy
from numpy import savetxt


################################################################################

ctrl = r'/content/PEAKScsvHRMS_G-RHA1_Ctrls.csv'
trtmnt = r'/content/PEAKScsvHRMS_G-RHA1_TRMT.csv'
output_file = 'novel_big_new_possibly_correct_mass_set.csv'
mass_error = .000005
internal_standard = 420.99373

################################################################################
start_time = time.time()


def average(lst):
    return reduce(lambda a, b: a + b, lst) / len(lst)


def cleaner(file_name):
    dirty = genfromtxt(file_name, delimiter=',', skip_header=8)
    l_listed = []
    dirty_shape = dirty.shape
    for x in range(int((dirty_shape[1] + 2) / 4)):
        l_listed.append(dirty[:, 4 * x].reshape(dirty_shape[0], 1))
        l_listed.append(dirty[:, 4 * x + 1].reshape(dirty_shape[0], 1))
    l = np.concatenate(l_listed, axis=1)
    dimensions = l.shape
    l = l.reshape((dimensions[0], int(dimensions[1] / 2), 2)).transpose((1, 0, 2))
    return l


def similar_peaks(peak, mass_error, matrix, size):
    fresh = list(np.where(abs(peak - matrix) / peak < mass_error))

    if fresh[0].size > size:
        coordinates = []
        for idx in range(len(fresh)):
            fresh[idx] = fresh[idx].tolist()
        for idx in range(len(fresh[0])):
            coordinates.append([fresh[0][idx], fresh[1][idx]])
        return coordinates


def get_trends(slot, lst):
    trends = []
    for time in range(1, len(lst[slot])):
        prev_pt = lst[slot][time - 1]
        time_pt = lst[slot][time]
        quotient = prev_pt / time_pt
        if quotient < .2:
            trends.append('Time ' + str(time - 1) + ' < ' + 'Time ' + str(time))
        elif quotient > 5:
            trends.append('Time ' + str(time - 1) + ' > ' + 'Time ' + str(time))
        elif quotient >= .2 and quotient <= 5:
            trends.append('Time ' + str(time - 1) + ' = ' + 'Time ' + str(time))
    if len(trends) > 1:
        trend_string = trends[0]
        for idx in range(1, len(trends)):
            trend_string += trends[idx][6:]
        return trend_string
    elif len(trends) == 1:
        return trends[0]
    else:
        return 'No Trend'


def get_relative_intensities(prism, internal_standard):
    test_prism = copy.deepcopy(prism)
    output_prism = copy.deepcopy(prism)
    tested = abs(test_prism - internal_standard) / internal_standard
    home = np.where(tested < mass_error)
    for layer_idx in range(prism.shape[0]):
        if layer_idx not in home[0]:
            return 'Internal standard not found in all samples'
    for idx_layer_idx in range(home[0].shape[0]):
        layer = home[0][idx_layer_idx]
        column = 1
        row = home[1][idx_layer_idx]
        output_prism[layer][:, 1] = output_prism[layer][:, 1] / prism[layer][row][column]
    return output_prism


def combine_sample_pairs(pair):
    scratch_prism = copy.deepcopy(pair)
    idx1 = 0
    mass_sample2 = scratch_prism[1][:, 0]
    inte_sample2 = scratch_prism[1][:, 1]
    mass_sample1 = scratch_prism[0][:, 0]
    inte_sample1 = scratch_prism[0][:, 1]
    new = []
    for idx in range(mass_sample1.shape[0]):
        element = mass_sample1[idx]
        comparison = np.where(abs(mass_sample2 - element) / element < mass_error)
        if comparison[0].tolist():
            new_row = np.reshape(np.array([(mass_sample2[comparison[0][0]] + element) / 2,
                                           (inte_sample2[comparison[0][0]] + inte_sample1[idx]) / 2]), (1, 2))
            new.append(new_row)
    newest = np.concatenate(new, axis=0)
    return newest


def doubler(relative):
    half_size = []
    for i in range(relative.shape[0]):
        if i % 2 == 0:
            half_size.append(combine_sample_pairs(relative[i:i + 2]))
    return half_size


def tripler(triple):
    mass_g1 = triple[0][:, 0]
    inte_g1 = triple[0][:, 1]
    mass_g2 = triple[1][:, 0]
    inte_g2 = triple[1][:, 1]
    mass_g3 = triple[2][:, 0]
    inte_g3 = triple[2][:, 1]

    new = []

    for idx in range(mass_g1.shape[0]):
        element = mass_g1[idx]
        inte_element = inte_g1[idx]
        comparison1 = np.where(abs(mass_g2 - element) / element < mass_error)
        comparison2 = np.where(abs(mass_g3 - element) / element < mass_error)
        if comparison1[0].tolist() and comparison2[0].tolist():
            new_row = np.reshape(np.array([(mass_g3[comparison2[0][0]] + mass_g2[comparison1[0][0]] + element) / 3, (
                        inte_g3[comparison2[0][0]] + inte_g2[comparison1[0][0]] + inte_element) / 3]), (1, 2))
            new.append(new_row)
        elif comparison1[0].tolist() and not comparison2[0].tolist():
            new_row = np.reshape(
                np.array([(mass_g2[comparison1[0][0]] + element) / 2, (inte_g2[comparison1[0][0]] + inte_element) / 2]),
                (1, 2))
            new.append(new_row)
        elif comparison2[0].tolist() and not comparison1[0].tolist():
            new_row = np.reshape(
                np.array([(mass_g3[comparison2[0][0]] + element) / 2, (inte_g3[comparison2[0][0]] + inte_element) / 2]),
                (1, 2))
            new.append(new_row)

    newer = [np.concatenate(new, axis=0)]
    mass_newer = newer[0][:, 0]

    for idx in range(mass_g2.shape[0]):
        element = mass_g2[idx]
        inte_element = inte_g2[idx]
        comparison1 = np.where(abs(mass_g3 - element) / element < mass_error)
        comparison2 = np.where(abs(mass_newer - element) / element < mass_error)
        if comparison1[0].tolist() and not comparison2[0].tolist():
            new_row = np.reshape(
                np.array([(mass_g3[comparison1[0][0]] + element) / 2, (inte_g3[comparison1[0][0]] + inte_element) / 2]),
                (1, 2))
            newer.append(new_row)

    newest = np.concatenate(newer, axis=0)
    return newest


def fixer(full_layer):
    masses = full_layer[:, 0].tolist()
    intensities = full_layer[:, 1].tolist()
    i = 0
    removals_masses, removals_intensities = [], []
    while i < len(masses) - 1:
        if abs(masses[i] - masses[i + 1]) / masses[i] < mass_error:
            masses[i] = average([masses[i], masses[i + 1]])
            intensities[i] = average([intensities[i], intensities[i + 1]])
            removals_masses.append(masses[i + 1])
            removals_intensities.append(intensities[i + 1])
        i += 1

    for i in range(len(removals_masses)):
        if removals_masses[i] in masses:
            masses.remove(removals_masses[i])
        else:
            if (removals_masses[i] + removals_masses[i + 1]) / 2 in masses:
                masses.remove((removals_masses[i] + removals_masses[i + 1]) / 2)
            else:
                print(removals_masses[i])

    for i in range(len(removals_intensities)):
        if removals_intensities[i] in intensities:
            intensities.remove(removals_intensities[i])
        else:
            if (removals_intensities[i] + removals_intensities[i + 1]) / 2 in intensities:
                intensities.remove((removals_intensities[i] + removals_intensities[i + 1]) / 2)
            else:
                print(removals_intensities[i])

    shortened_masses = np.reshape(np.asarray(masses), (len(masses), 1))
    shortened_intensities = np.reshape(np.asarray(intensities), (len(intensities), 1))
    return np.concatenate((shortened_masses, shortened_intensities), axis=1)


def get_all_sets(matrix):
    shape = matrix.shape
    bank = []
    full_mass_set = []
    full_intensity_set = []
    mass_matrix = matrix[:, :, 0]
    for layer_idx in range(shape[0]):
        for row_idx in range(mass_matrix.shape[1]):
            if [layer_idx, row_idx] not in bank:
                value = mass_matrix[layer_idx][row_idx]
                if str(value) != str(np.nan):
                    error_matrix = abs(mass_matrix - value) / value
                    positives = np.where(error_matrix < mass_error)
                    intensity_set = np.empty((shape[0]))
                    intensity_set[:] = np.NaN
                    intensity_set = intensity_set.tolist()
                    mass_set = []
                    for idx in range(positives[0].shape[0]):
                        bank.append([positives[0][idx], positives[1][idx]])
                        intensity_set[positives[0][idx]] = matrix[positives[0][idx]][positives[1][idx]][1]
                        mass_set.append(matrix[positives[0][idx]][positives[1][idx]][0])
                full_intensity_set.append(intensity_set)
                full_mass_set.append(average(mass_set))
    trtmnt_masses = []
    [trtmnt_masses.append(x) for x in full_mass_set if x not in trtmnt_masses]
    trtmnt_intensities = []
    [trtmnt_intensities.append(x) for x in full_intensity_set if x not in trtmnt_intensities]

    return trtmnt_masses, trtmnt_intensities


def fer(list_of_matrices):
    new = copy.deepcopy(list_of_matrices)
    for idx in range(len(list_of_matrices)):
        new[idx] = fixer(list_of_matrices[idx])
    return new


def condenser_reshaper(half_size):
    proper_condensed = []
    sizer = []
    for i in range(len(half_size)):
        if i % 3 == 0:
            combined_triple = tripler(half_size[i:i + 3])
            proper_condensed.append(combined_triple)
            sizer.append(combined_triple.shape[0])
    fixed = fer(fer(proper_condensed))
    comparable = []
    for matrix in fixed:
        max_rows = max(sizer)
        if matrix.shape[0] != max_rows:
            filler = np.full((max_rows - matrix.shape[0], 2), np.nan)
            resized = np.concatenate((matrix, filler))
            comparable.append(resized)
        else:
            comparable.append(matrix)
    comparable = np.concatenate(comparable, axis=1)
    header = []
    for i in range(int(comparable.shape[1] / 2)):
        header.append('Time ' + str(i))
        header.append('')
    array_header = np.reshape(np.asarray(header), (1, comparable.shape[1]))
    checkpt = np.concatenate((array_header, comparable))
    np.savetxt('checkpoint.csv', checkpt, delimiter = ',', fmt = '%s')
    dimensions = comparable.shape
    reshapen = comparable.reshape((dimensions[0], int(dimensions[1] / 2), 2)).transpose((1, 0, 2))
    return reshapen


def First(elem):
    return elem[0]


def final_writer(output_file, complete):
    with open(output_file, mode='w') as finished_file:
        final_data = csv.writer(finished_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        
        header = ['m/z', 'Trends']
        for sample_number in range(len(complete[0][-2])):
            header.append('Intensity (Time ' + str(sample_number) + ')')
        header.append('Control Presence')
        header.append('Control Frequency')
        final_data.writerow(header)
        for row in complete:
            if len(row) == 5:
                row[:-3] += row[2]
                row.pop(-3)
                final_data.writerow(row)
            else:
                row[:-2] += row[2]
                row.pop(-2)
                row.append('')
                final_data.writerow(row)


entire_prism = cleaner(trtmnt)

relative = get_relative_intensities(entire_prism, internal_standard)

half_size = doubler(relative)

condensed = condenser_reshaper(half_size)

trtmnt = get_all_sets(condensed)

control_entire_prism = get_relative_intensities(cleaner(ctrl), internal_standard)

control_mass_matrix = control_entire_prism[:, :, 0]

all_matches_ctrl = []

completely_absent_masses = []

completely_absent_intensities = []


for idx1 in range(len(trtmnt[0])):
    control_matches = similar_peaks(trtmnt[0][idx1], mass_error, control_mass_matrix, 0)
    if control_matches:
        all_matches_ctrl.append(control_matches)
    else:
        all_matches_ctrl.append([])
        completely_absent_masses.append(trtmnt[0][idx1])
        completely_absent_intensities.append(trtmnt[1][idx1])

control_indices = all_matches_ctrl


for idx1 in range(len(all_matches_ctrl)):
    if all_matches_ctrl[idx1]:
        for idx2 in range(len(all_matches_ctrl[idx1])):
            var = all_matches_ctrl[idx1][idx2]
            all_matches_ctrl[idx1][idx2] = control_entire_prism[var[0]][var[1]][1]

satisfied_intensities = []
for idx in range(len(all_matches_ctrl)):
    intensity_list = all_matches_ctrl[idx]
    if intensity_list:
        maximum = max(intensity_list)
        fake_intensities = []
        for intensity in trtmnt[1][idx]:
            if str(intensity) != str(np.nan):
                fake_intensities.append(intensity)

        if maximum < .0001 and min(fake_intensities) > 100 * maximum:
            satisfied_intensities.append([trtmnt[0][idx], get_trends(idx, trtmnt[1]), trtmnt[1][idx], 'Low Presence', len(intensity_list)])


complete = satisfied_intensities
[complete.append([completely_absent_masses[x], get_trends(x, completely_absent_intensities), completely_absent_intensities[x], 'No Presence']) for x in range(len(completely_absent_masses))]
complete.sort(key=First)

final_writer(output_file, complete)

print("0 to finish--- %s seconds ---" % (time.time() - start_time))
