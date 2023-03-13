from comparison_functions import cleaner, get_relative_intensities, doubler, condenser_reshaper, get_all_sets, similar_peaks, get_trends, First, final_writer
import time
import numpy as np

################################################################################

ctrl = r'C:\Users\milga\Documents\HRMS Project\Data\PEAKScsvHRMS_G-RHA1_Ctrls.csv'
trtmnt = r'C:\Users\milga\Documents\HRMS Project\Data\betaPEAKScsvHRMS_G-RHA1_TRMT.csv'
output_file = 'alphanovel_big_new_possibly_correct_mass_set.csv'
mass_error = .000005
internal_standard = 420.99373

################################################################################
start_time = time.time()


entire_prism = cleaner(trtmnt)

relative = get_relative_intensities(entire_prism, internal_standard, mass_error)

half_size = doubler(relative, mass_error)

condensed = condenser_reshaper(half_size, mass_error)

trtmnt = get_all_sets(condensed, mass_error)

control_entire_prism = get_relative_intensities(cleaner(ctrl), internal_standard, mass_error)

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
complete.sort( key = First)

final_writer(output_file, complete)

print("0 to finish--- %s seconds ---" % (time.time() - start_time))
