#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 06:14:22 2023

@author: Alexandros Koukras
"""

import sunpy.map
import drms
from sunpy.net import Fido
from sunpy.net import Fido, attrs as a
import astropy.units as u
from datetime import datetime,timedelta
import pickle

from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy as np

#------------------------------------------------------------------------------
def download_aia_based_on_bmt_drms(bmt,wavelength,aiadir_name,aialist_name,
                              save_aialist=True):
    """ Similar to download_aia_based_on_bmt but it uses the drms package and
    not Fido to retrieve the AIA data.
    
    TODO:
        Adjust this function in order to be more generic. It should accept a date
        and download the closest AIA full disk image at that date (given a particular channel)
        
        It should also be able to differentiate between a single date or a list of
        dates and perfom the same task in both cases
    Notes
    -----
    There are 2 AIA products that are of interest, the 'aia.lev1' and the
    'aia.lev1_euv_12s'. It seems that with Fido I was downloading the first
    but with drms I can only download the second because this is the one 
    that has the wavelength as a prime key. I compared the two different files
    and seem almost identical (in the header)."""
    
    # Set the export parameters
    c = drms.Client(email='alexandros.koukras@gmail.com')
    aiaseries = 'aia.lev1_euv_12s'
    aiasegment = 'image'
    
    aiafile_list = []
    # Iterate on the list of bm times 
    for bmtime in bmt:
        print(bmtime)
        ds = '%s[%s][%s]{%s}'%(aiaseries,str(bmtime.value),str(wavelength),
                               aiasegment)
        # Search for the AIA file
        r = c.export(ds,method='url',protocol='fits')
        r.wait()
        # Download it at the specified path
        r.download(aiadir_name,verbose=True)
        # Save the location+name of the file
        aiafile_list.append(aiadir_name + '/' + r.data.iloc[0,1])
    
    # Save the list of AIA files in the same directory
    if save_aialist==True:
        save_aialist_name = aiadir_name+'/'+aialist_name+'.pkl'
        file = open(save_aialist_name,'wb')
        pickle.dump(aiafile_list,file)
        file.close()   
        print('The AIA list was writen in the file:')
        print(save_aialist_name)
        
#------------------------------------------------------------------------------
def y_shift_eis_map(eismap,filling_value=0.):
    """ Corrects for the 'grating tilt' by shifting an EISmap in the y-direction
    based on the observed wavelength. The y-shift is based on the He II 256 A line
    and is derived from the eispac.instr.ccd_offset routine. The user can specify
    the filling value but the default is 0.
    
    Parameters:
    -----------
        eismap: EISmap (sunpy.map.Map object)
            The intensity/velocity/width sunpy map for a certain EIS observation.
            This map is typically retrieved from the fit_res.get_map() routine.
            
    Returns:
    --------
        eismap_shifted: EISmap (sunpy.map.Map object)
            The original EIS map but shifted in the y-direction.
    
    Notes:
    ------
        Details about the ccd_offset routine can be found here:
        https://eispac.readthedocs.io/en/latest/api/eispac.instr.ccd_offset.html
    
    """
    
    from eispac.instr import ccd_offset
    
    y_offset = round(ccd_offset(eismap.wavelength.value)[0]) #in pixels
    
    if  y_offset != 0.:
        # Create an array of zeros (default) for the offset portion
        leny = eismap.data.shape[0]
        lenx = eismap.data.shape[1]
        if filling_value==0:
            arr_shift = np.zeros((abs(y_offset),lenx))
        else:
            arr_shift = np.full([abs(y_offset), lenx], filling_value)

        # Make the new 'shifted' eis map. 
        # The direction of the shift depends on the offset value 
        if y_offset>0.:
            tmp_arr = np.concatenate((eismap.data,arr_shift),axis=0) 
            eisdata_shifted = tmp_arr[y_offset:,:]
        else: 
            tmp_arr = np.concatenate((arr_shift,eismap.data),axis=0) 
            eisdata_shifted = tmp_arr[0:leny,:]

        eismap_shifted = sunpy.map.Map(eisdata_shifted,eismap.meta)
        
        return eismap_shifted
    
    else:
        print('Not offset was necessary for wavelength: %s'%eismap.wavelength.value)
        return eismap
    
#------------------------------------------------------------------------------
def replace_bad_pixels_eis_intmaps(intmap,option1=True,
                                  option1_upper_threshold=1000.,
                                  option1_lower_threshold=0.):
    """ Replaces 'bad' pixels in an EIS intensity map. 

    The bad pixels are replaced with zero values. The detection of the bad pixels
    is based on the different options.
    Option1: 
        Sets pixels above and below a threshold value to zero.
        
    TODO:
        - Add another option based on the uncertainty 
        - Add different options for filling the value of the 'bad' pixel 
        (e.g. interpolation from the neighboring pixels)
        
    Notes:
    ------
        There might be a problem with the uncertainty values when I perfom the 
        following steps. The 'fixed' maps I return have the uncertainty attribute 
        as NoneType.
    """
    
    print('For %s'%intmap.name)
    n_above_upper = intmap.data[intmap.data>option1_upper_threshold].size
    n_below_lower = intmap.data[intmap.data<option1_lower_threshold].size
    print('%d pixel above %.f'%(n_above_upper,option1_upper_threshold))
    print('%d pixel below %.f'%(n_below_lower,option1_lower_threshold))
    
    intmap.data[intmap.data>option1_upper_threshold] = 0. 
    intmap.data[intmap.data<option1_lower_threshold] = 0. 

    return intmap

#------------------------------------------------------------------------------
def get_gofnt_for_eismap(eismap,temp=None,dens=None):#,dens=1e+9):
    """ Retrieve the contribution function (G(n,T)) for the 'main' line in an
    EISmap.
    
    Parameters:
    -----------
        eismap : EISmap (sunpy.map.Map object)
            Most commonly an intensity map, that is retrieved from the 
            fit_res.get_map routine.
            
        temp : ndarray
            The temperature values to be used to calculate the contribution
            function. Default log(T):5.8-7.3 (step 0.05)
            
        dens : float, default:1e+9
            The density value to be used to calculate the contribution function.
            
    Returns:
    --------
        gofnt : array 
            The contribution function
    
    Notes:
    -----
        It's more practical to calculate the G(n,T) with the steps I follow,
        rather than using the gofnt function and the interactive GUI. This 
        is also mentioned in the ChiantiPy notes:
        https://chiantipy.readthedocs.io/en/latest/quick_start.html#g-n-t-function
        
        The roman library is used to convert from roman numeral to numbers. If 
        not available it can be install by: pip install roman.
        
    TODO:
    -----
        - Add the temperature and density as input parameters 
        
    """
    
    import ChiantiPy.core as ch
    import roman 
    
    if isinstance(temp,np.ndarray)!=True:
        temp = 10.**(5.8 + 0.05*np.arange(31.))  
    if isinstance(dens,np.float)!=True:
        dens = 1.e+9

    print(dens)    
    # Get the necessary information of the line in the map
    ion_name, ion_numeral, ion_wvl = eismap.meta['line_id'].split()
    
    # Convert them to a form compatible with ChiantiPy 
    ion_name_chianti = "_".join([ion_name.lower(),
                                 str(roman.fromRoman(ion_numeral))])

    # Set up the ion in ChiantiPy and get the G(n,T)
    ion_n = ch.ion(ion_name_chianti,temp,dens)
    ion_n.intensity()
    dist = np.abs(np.asarray(ion_n.Intensity['wvl']) - float(ion_wvl))
    idx = np.argmin(dist)
    gofnt = ion_n.Intensity['intensity'][:,idx]
    
    return gofnt


#------------------------------------------------------------------------------
def get_gofnt_for_pixel(line_id,temp=None,dens=None,use_sum_in_range=False,
                        use_max_in_range=False,wvl_sigma=0.01):
    """ Retrieve the contribution function (G(n,T)) for a pixel on an ion 
    intensity map. """
    
    import ChiantiPy.core as ch
    import roman 
    
    if isinstance(temp,np.ndarray)!=True:
        raise Exception('The temperature array must be provided')  
    if isinstance(dens,np.float)!=True:
        raise Exception('The density value must be provided')  

    ion_name, ion_numeral, ion_wvl = line_id.split()
    
    # Convert them to a form compatible with ChiantiPy 
    ion_name_chianti = "_".join([ion_name.lower(),
                                 str(roman.fromRoman(ion_numeral))])

    # Set up the ion in ChiantiPy and get the G(n,T)
    ion_n = ch.ion(ion_name_chianti,temp,dens)
    ion_n.intensity()

    
    if use_sum_in_range:
        # Return the sum of G(n,T)s in a specific wavelength range around the 
        # wavelength that is indicated in the file header.
        # This is often more appropriate as the wavelength of the ion in the header 
        # is not the same as the one in ChiantiPy.
        
        wvl_range = np.linspace(float(ion_wvl)-wvl_sigma,float(ion_wvl)+wvl_sigma,num=5,
                                endpoint=True)
        gofnt_range = []
        for wvl in wvl_range:
            dist = np.abs(np.asarray(ion_n.Intensity['wvl']) - float(wvl))
            idx = np.argmin(dist)
            gofnt_range.append(ion_n.Intensity['intensity'][:,idx])
            
        gofnt_sum = np.sum(gofnt_range,axis=0)
        return gofnt_sum

    elif use_max_in_range:
        # Return the largest G(n,T) inside the specified wavelength range
        wvl_range = np.linspace(float(ion_wvl)-wvl_sigma,float(ion_wvl)+wvl_sigma,num=5,
                                endpoint=True)
        gofnt_range = []
        max_list = [] #very hacky way, but it works for testing
        for wvl in wvl_range:
            dist = np.abs(np.asarray(ion_n.Intensity['wvl']) - float(wvl))
            idx = np.argmin(dist)
            gg = ion_n.Intensity['intensity'][:,idx]
            gofnt_range.append(gg)
            max_list.append(np.max(gg))
        
        return gofnt_range[np.argmax(max_list)]
        
    else:
        # Return the G(n,T) for the specific wavelength in the header
        dist = np.abs(np.asarray(ion_n.Intensity['wvl']) - float(ion_wvl))
        idx = np.argmin(dist)
        gofnt = ion_n.Intensity['intensity'][:,idx]
        
        return gofnt

def plot_fit_profile(data_cube,fit_res,iy,ix,axes):
    ''' Plot the fitting profile with all the components for a specific pixel.
    
    Given a fit_res object and coordinates it plots the real data and the 
    fitting profile with all the components.
    
    Parameters:
    -----------
        data_cube : 
        fit_res : 
            
        iy : int
            The index of the pixel in the y-direction.
        ix : int
            The index of the pixel in the x-direction.
    '''    
    
    data_x = data_cube.wavelength[iy, ix, :]
    data_y = data_cube.data[iy, ix, :]
    data_err = data_cube.uncertainty.array[iy, ix, :]
    # Remove bad data 
    data_y[data_y<0] = np.nan 
    data_err[data_err<0] = np.nan 
    
    fit_x, fit_y = fit_res.get_fit_profile(coords=[iy,ix], num_wavelengths=100)
    
    n_gauss = fit_res.fit['n_gauss']
    num_wvl = 100
    if n_gauss==1:
        main_comp_id = fit_res.fit['main_component']
        background_id = main_comp_id - 1
        c0_x, c0_y = fit_res.get_fit_profile(main_comp_id, coords=[iy,ix], 
                                             num_wavelengths=num_wvl)
        cb_x, cb_y = fit_res.get_fit_profile(background_id, coords=[iy,ix],
                                             num_wavelengths=num_wvl)
        
        # Plotting 
        axes.errorbar(data_x, data_y, yerr=data_err, ls='', marker='o', 
                      color='k')
        axes.plot(fit_x, fit_y, color='b', label='Fit profile')
        axes.plot(c0_x, c0_y, color='r', label=fit_res.fit['line_ids'][0])
        axes.plot(cb_x, cb_y, color='c', ls='--',label='Background')
        
    elif n_gauss>1:
        main_comp_id = fit_res.fit['main_component']
        if main_comp_id==0:
            # Same as before 
            background_id = main_comp_id - 1            
            c0_x, c0_y = fit_res.get_fit_profile(main_comp_id, coords=[iy,ix], 
                                                 num_wavelengths=num_wvl)
            cb_x, cb_y = fit_res.get_fit_profile(background_id, coords=[iy,ix],
                                             num_wavelengths=num_wvl)
        
            # Get the fits for the other components 
            cn_x, cn_y = [],[]
            for c_id in range(1,n_gauss):
                tmpx, tmpy = fit_res.get_fit_profile(c_id, coords=[iy,ix], 
                                                     num_wavelengths=num_wvl)
                cn_x.append(tmpx)
                cn_y.append(tmpy)
                
            # Plotting 
            axes.errorbar(data_x, data_y, yerr=data_err, ls='', marker='o', 
                          color='k')
            axes.plot(fit_x, fit_y, color='b', label='Fit profile')
            axes.plot(c0_x, c0_y, color='r', label=fit_res.fit['line_ids'][0])
            axes.plot(cb_x, cb_y, color='c', ls='--',label='Background')
            for i in range(n_gauss-1):
                axes.plot(cn_x[i], cn_y[i], label=fit_res.fit['line_ids'][i+1])
    
        else:
            print('The index of the main component in a multi-gaussian fit, was not 0')
            pass
            # To be continued 
    
    axes.set_xlabel('Wavelength [$\AA$]')
    axes.set_ylabel('%s'%data_cube.unit)
    axes.set_title('Fit profile for x:%.d and y:%.d'%(ix,iy))
    axes.legend(loc='upper left', frameon=False)
    
    return 

