import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord#, match_coordinates_sky
from astropy import units as u
from astropy.table import Table, vstack
#import matplotlib.pyplot as plt
import sys
import multiprocessing as mp
import time
from collections import Counter

#def Repeat(x):
#    _size = len(x)
#    repeated = []
#    for i in range(_size):
#        k = i + 1
#        for j in range(k, _size):
#            if x[i] == x[j] and x[i] not in repeated:
#                repeated.append(x[i])
#    return repeated

def cross_match_2_tables(overlaps_final, overlaps_field):
        # Because we want to make sure the objects are the same, and we want to require that have compatible magnitudes, we save the difference between the matches found.
        #tolerable_difference = 2.0*(overlaps_final['GERR'] + overlaps_field['GERR'])
        tolerable_difference = [0.2]*len(overlaps_final['GERR'])
        mag_difference = abs(overlaps_final['G'] - overlaps_field['G'])
        gerr_final = overlaps_final['GERR']
        gerr_field = overlaps_field['GERR']
        counter = np.arange(len(gerr_final))

        #This is for sanity check purposes.
        stars_eliminated = 0
        remove_duplicate_f2 = 0
        remove_duplicate_f1 = 0
        remove_f2 = 0
        remove_f1 = 0
        two_objects = 0
        entries_to_remove_1 = []
        entries_to_remove_2 = []

        # Ok so now we start to go through the matches found to make some decisions.
        iter_time = []
        for i, diff, tol_diff, gerr1, gerr2 in zip(counter, mag_difference, tolerable_difference, gerr_final, gerr_field):
            start = time.time()
            #if gerr1 == 99.99:
            #    overlaps_final.remove_row(i-1-remove_f1-remove_duplicate_f1)
            #    stars_eliminated += 1
            #    remove_f1 += 1
            #    if gerr2 == 99.99:
            #        overlaps_field.remove_row(i-1-remove_f2-remove_duplicate_f2)
            #        stars_eliminated += 1
            #        remove_f2 += 1
            #    continue
            #if gerr2 == 99.99:
            #    overlaps_field.remove_row(i-1-remove_f2-remove_duplicate_f2)
            #    stars_eliminated += 1
            #    remove_f2 += 1
            #    if gerr1 == 99.99:
            #        overlaps_final.remove_row(i-1-remove_f1-remove_duplicate_f1)
            #        stars_eliminated += 1
            #        remove_f1 += 1
            #    continue

            # First of all, is important to check if the star has been accounted for in the previous loop. If it has, remove it from the previous row, unless it has already been removed.
            if overlaps_field['ID'][i] == overlaps_field['ID'][i-1]:
                if remove_field == False:
                    entries_to_remove_2.append(i)
                    #overlaps_field.remove_row(i-1-remove_f2-remove_duplicate_f2)
                    remove_duplicate_f2 += 1
                    stars_eliminated += 1
                else:
                    entries_to_remove_2.append(i)
                    #overlaps_field.remove_row(i-remove_f2-remove_duplicate_f2)
                    remove_duplicate_f2 += 1
                    two_objects += 1
                    stars_eliminated += 1
                    continue
            if overlaps_final['ID'][i] == overlaps_final['ID'][i-1]:
                if remove_final == False:
                    entries_to_remove_1.append(i)
                    #overlaps_final.remove_row(i-1-remove_f1-remove_duplicate_f1)
                    remove_duplicate_f1 += 1
                    stars_eliminated += 1
                else:
                    entries_to_remove_1.append(i)
                    #overlaps_final.remove_row(i-remove_f1-remove_duplicate_f1)
                    remove_duplicate_f1 += 1
                    two_objects += 1
                    stars_eliminated += 1
                    continue

            # We check if they are the same source by comparing magnitudes.
            if diff <= tol_diff:
                # If they are the same, we check which table has the best data by comparing the error in the magnitudes.
                if gerr1 < gerr2:
                    # If the initial table has the best data, we remove the other data row.
                    entries_to_remove_2.append(i)
                    #overlaps_field.remove_row(i-stars_eliminated)
                    remove_field = True
                    remove_final = False
                    remove_f2 += 1
                else:
                    # And viceversa.
                    entries_to_remove_1.append(i)
                    #overlaps_final.remove_row(i-stars_eliminated)
                    remove_final = True
                    remove_field = False
                    remove_f1 += 1
                stars_eliminated += 1
            else:
                # If the objects are not the same, we do not need to bother and we will include them both in the final catalogue.
                remove_final = False
                remove_field = False
                two_objects += 1
                continue
            stop = time.time()
            iter_time.append(stop-start)
        overlaps_field.remove_rows(entries_to_remove_2)
        overlaps_final.remove_rows(entries_to_remove_1)
        print('Hello there, I have taken '+str(np.mean(iter_time))+' seconds for each cross-match operation.')
        #print(overlaps_field.info)
        #print(overlaps_final.info)

        table = vstack([overlaps_final, overlaps_field])
        #print(table.info)
        #print(stars_eliminated, remove_f1, remove_f2, remove_duplicate_f1, remove_duplicate_f2, two_objects)
        return table



field_list = ['3', '4', '5', '6', '7', '9', '10', '11', '12', '14', '15', '16', '178'] #SMC fields
#field_list = ['3','4'] #test fields


data_path = '/scratch/pmassana/SMASH_DATA/allobj_files/'
# We read the first one to start of.
table = Table.read(data_path+'Field'+field_list[0]+'_combined_allobj.fits')
# Removing columns that will give problems afterwards if we leave them there and making the file lighter.
table.keep_columns(['ID','RA','DEC', 'RAERR', 'U', 'UERR', 'G', 'GERR', 'R', 'RERR', 'I', 'IERR', 'Z', 'ZERR', 'CHI', 'SHARP', 'PROB', 'EBV'])
# Removing rows with no photometry.
bad_phot_idcs = np.where(table['G']==99.99)
table.remove_rows(bad_phot_idcs)
#We want to make sure that all the rows have positional data...
bad_phot_idcs = np.where(table['RA']!=table['RA'])
table.remove_rows(bad_phot_idcs)

#Now we want to start to go through all the fields.
for field in field_list:
    # We want to skip the first one because we have already read it.
    if field==field_list[0]:
        continue

    # Now we read the next field in the list and do the same as in the previous one.
    table2 = Table.read(data_path+'Field'+field+'_combined_allobj.fits')
    table2.keep_columns(['ID', 'RA', 'DEC', 'RAERR', 'U', 'UERR', 'G', 'GERR', 'R', 'RERR', 'I', 'IERR', 'Z', 'ZERR', 'CHI', 'SHARP', 'PROB', 'EBV'])
    bad_phot_idcs = np.where(table2['G']==99.99)
    table2.remove_rows(bad_phot_idcs)
    #We want to make sure that all the rows have positional data...
    bad_phot_idcs = np.where(table2['RA']!=table2['RA'])
    table2.remove_rows(bad_phot_idcs)

    # Now is the time to read the coordinates for two.
    coordinates_field = SkyCoord(ra=table2['RA']*u.degree, dec=table2['DEC']*u.degree)
    coordinates_final = SkyCoord(ra=table['RA']*u.degree, dec=table['DEC']*u.degree)
    # We find the indeces and distances between all the objects that are within a distance less than the maximum error in both catalogues.
    #search_radius = 1.0 + 3.0*np.std(np.hstack([table['RAERR'], table2['RAERR']]))
    search_radius = 1.5
    idx_field, idx_final, d2d, d3d = coordinates_final.search_around_sky(coordinates_field, search_radius*u.arcsec)
    print('There are '+str(len(idx_field))+' pairs that are less than '+str(search_radius*u.arcsec)+' from each other.')
    counter_idx1 = Counter(idx_final)
    counter_idx2 = Counter(idx_field)
    #idx_field_repeated_list = Repeat(idx_field)
    #idx_final_repeated_list = Repeat(idx_final)
    #print('Finished creating rep. list for field '+field)

    duplicates1 = []
    duplicates2 = []
    single_pairs = []
    for idx1, idx2 in zip(idx_final, idx_field):
        if counter_idx1[idx1] > 1:
            duplicates1.append([idx1, idx2])
            continue
        if counter_idx2[idx2] > 1:
            duplicates2.append([idx1, idx2])
            continue
        else:
            single_pairs.append([idx1, idx2])
    duplicates1, duplicates2, single_pairs = np.array(duplicates1), np.array(duplicates2), np.array(single_pairs)
    if len(single_pairs) > 0:
        #print(duplicates1)
        duplicates1 = duplicates1[duplicates1[:,0].argsort()]
        duplicates2 = duplicates2[duplicates2[:,1].argsort()]
        #print(duplicates1.shape, duplicates2.shape, single_pairs.shape)
        #sys.exit()
        # Now we create four tables. Two that contain the overlapping stars from each respective incoming table. The other two contain the rest of the stars for which pairs have not been found.
        overlaps_final_single = table[single_pairs.T[0]]
        overlaps_field_single = table2[single_pairs.T[1]]
        overlaps_final_d1 = table[duplicates1.T[0]]
        overlaps_field_d1 = table2[duplicates1.T[1]]
        overlaps_final_d2 = table[duplicates2.T[0]]
        overlaps_field_d2 = table2[duplicates2.T[1]]
        table.remove_rows(idx_final)
        table2.remove_rows(idx_field)

        print('Starting cross-match for field '+field)
        print('Cross-matching '+str(len(single_pairs))+' single pairs.')
        overlaps_single = cross_match_2_tables(overlaps_final_single, overlaps_field_single)
        print('Cross-matching '+str(len(duplicates1))+' duplicates in the reference table.')
        overlaps_d1 = cross_match_2_tables(overlaps_final_d1, overlaps_field_d1)
        print('Cross-matching '+str(len(duplicates2))+' duplicates in the second table.')
        overlaps_d2 = cross_match_2_tables(overlaps_final_d2, overlaps_field_d2)
        table = vstack([table, table2, overlaps_single, overlaps_d1, overlaps_d2])
        print('Finished cross-match for field '+field)
    else:
        table = vstack([table, table2])
    print('Field '+field+' done')

print(table.info)
print('Starting to write into file...')
col1 = fits.Column(name='ID', format='14A', array=table['ID'])
col2 = fits.Column(name='RA', format='D', array=table['RA'])
col3 = fits.Column(name='DEC', format='D', array=table['DEC'])
col4 = fits.Column(name='RAERR', format='D', array=table['RAERR'])
col6 = fits.Column(name='U', format='D', array=table['U'])
col7 = fits.Column(name='UERR', format='D', array=table['UERR'])
col8 = fits.Column(name='G', format='D', array=table['G'])
col9 = fits.Column(name='GERR', format='D', array=table['GERR'])
col10 = fits.Column(name='R', format='D', array=table['R'])
col11 = fits.Column(name='RERR', format='D', array=table['RERR'])
col12 = fits.Column(name='I', format='D', array=table['I'])
col13 = fits.Column(name='IERR', format='D', array=table['IERR'])
col14 = fits.Column(name='Z', format='D', array=table['Z'])
col15 = fits.Column(name='ZERR', format='D', array=table['ZERR'])
col16 = fits.Column(name='CHI', format='D', array=table['CHI'])
col17 = fits.Column(name='SHARP', format='D', array=table['SHARP'])
col18 = fits.Column(name='PROB', format='D', array=table['PROB'])
col19 = fits.Column(name='EBV', format='D', array=table['EBV'])
hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19])
hdu.writeto(data_path+'combined_fields_allobj_largeG_02_error.fits', overwrite = True)
print('Finished writing into file.')
