#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Copy H5parm values from a reference station to all new stations, e.g., from ST001 to all core stations

Created on Tue Aug 28 2018

@author: Alexander Drabent (parts by Maaijke Mevius)
"""

import argparse
import logging
import os
import tables

from losoto.h5parm import h5parm
from losoto.lib_operations import reorderAxes

import numpy as np
import casacore.tables as ct


def makesolset(MS, data, solset_name, useh5coordinates):
    ''' Create a new solset.

    Args:
        MS (str): name of the input measurement set.
        data (h5parm): h5parm object to add the solset to.
        solset_name (str): name of the new solset.
    Returns:
        antennaNames (ndarray): array containing the names of antennas in the measurement set.
    '''

    solset = data.makeSolset(solset_name)

    antennaFile = MS + "::ANTENNA"
    logging.info('Collecting information from the ANTENNA table.')
    antennaTable = ct.table(antennaFile, ack=False)
    antennaNames = antennaTable.getcol('NAME')
    antennaPositions = antennaTable.getcol('POSITION')
    antennaTable.close()
    antennaTable = solset.obj._f_get_child('antenna')
    antennaTable.append(zip(*(antennaNames, antennaPositions)))

    fieldFile = MS + "::FIELD"
    logging.info('Collecting information from the FIELD table.')
    fieldTable = ct.table(fieldFile, ack=False)
    phaseDir = fieldTable.getcol('PHASE_DIR')
    pointing = phaseDir[0, 0, :]
    print ('pointing ', pointing)
    if pointing[0] < 0:
        pointing[0] = pointing[0] + 2.*np.pi
    fieldTable.close()

    sourceTable = solset.obj._f_get_child('source')
    # Add the field centre, that is also the direction for Gain and Common*
    print ('sourceTable ', sourceTable[:])
    
    if useh5coordinates:
        solset0 = data.getSolset('sol000')
        direction = solset0.obj.source[:]
        print ('useh5coordinates==True, so the dir should be taken from the given .h5 file at sol000.source ', direction)
        for i in range(len(direction)):
            if direction[i-1][1][0] < 0:
                direction[i-1][1][0] += 2.*np.pi
        #if direction[0] < 0:
        #    direction[0] = direction[0] + 2.*np.pi
        sourceTable.append(direction)
    else:
        sourceTable.append([('pointing', pointing)])
    print ('Updated sourcetable is ', sourceTable[:])

    return antennaNames


def main(h5parmfile, MSfiles, solset_in='sol000', solset_out='sol001', soltab_list='phase000,amplitude000', superstation='ST001', restrictToCS=True, matchPtg=False, useh5coordinates=True):
    ''' Copy the gains from the phased up core back to all core stations.

    Args:
        h5parmfile (str): input H5Parm
        MSfiles (list): list of (a) measurement set(s) from which the stations are read.
        solset_in (str): input solset to process.
        solset_out (str): output solset with the core stations added.
        soltab_list (list): list of strings, containing the soltabs to be processed.
        superstation (str): name of the phased up station.
        restrictToCS (bool): only do this operation for stations starting with CS. Default: True
        useh5coordinates (bool): copy the source coordinate from the input .h5 file at sol000

    Returns:
        0 or 1 (int): 0 if succesfull, 1 if an error occured.
    '''

    mslist = MSfiles.lstrip('[').rstrip(']').replace(' ', '').replace("'", "").split(',')

    soltab_list = soltab_list.split(',')
    #print (soltab_list)

    if len(mslist) == 0:
        logging.error("Did not find any existing directory in input MS list!")
        return 1
    else:
        MSfile = mslist[0]

    if not os.path.exists(h5parmfile):
        logging.error("H5parm file %s doesn't exist!" % h5parmfile)
        return 1

    if not os.path.exists(MSfile):
        logging.error("MS file %s doesn't exist!" % MSfile)
        return 1

    if solset_in == solset_out:
        logging.error("Output solset has to be different from input solset!")
        return 1

    # Open up the h5parm, get an example value
    H = tables.open_file(h5parmfile, 'r+')
    station_names = H.root.sol000.phase000.ant[:]
    print ('The input .h5 file has stations at sol000 shown as the following: ', station_names)
    H.close()
    
    data = h5parm(h5parmfile, readonly=False)

    # Create a new solset for the data
    if solset_out not in data.getSolsetNames():
        new_station_names = makesolset(MSfile, data, solset_out, useh5coordinates)
        

    print ('The output .h5 file has stations as the following: ', new_station_names)
    # loading solset
    solset = data.getSolset(solset_in)
    OutSolset = data.getSolset(solset_out)

    #station_names = solset.getAnt().keys()
        
    if superstation not in station_names:
        logging.error("Couldn't find station " + str(superstation))
        return 1

    # start copying
    for soltab_name in soltab_list:
        logging.info("Running copySTgains_toCS on: " + soltab_name)
        soltab = solset.getSoltab(soltab_name)
        STindex = soltab.getAxisValues('ant', ignoreSelection=True).tolist().index(superstation)

        if 'clock' in soltab_name or 'tec' in soltab_name:
            for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=['time', 'ant'], weight=True):
                vals = reorderAxes(vals, soltab.getAxesNames(), ['time', 'ant'])
                weights = reorderAxes(weights, soltab.getAxesNames(), ['time', 'ant'])
        elif 'amplitude' in soltab_name or 'phase' in soltab_name:
          if 'pol' in soltab.getAxesNames():  
            for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=['pol', 'ant', 'freq', 'time', 'dir'], weight=True):
                #vals = reorderAxes(vals, soltab.getAxesNames(), ['time', 'ant', 'freq', 'pol'])
                #weights = reorderAxes(weights, soltab.getAxesNames(), ['time', 'ant', 'freq', 'pol'])
                vals = reorderAxes(vals, soltab.getAxesNames(), ['time', 'freq', 'ant', 'dir', 'pol'])
                weights = reorderAxes(weights, soltab.getAxesNames(), ['time', 'freq', 'ant', 'dir', 'pol'])
                print ('HERE!', vals.shape)
          else: # in case we have have no polarization axis, so scalarphase-type      
            for vals, weights, coord, selection in soltab.getValuesIter(returnAxes=['ant', 'freq', 'time', 'dir'], weight=True):
                vals = reorderAxes(vals, soltab.getAxesNames(), ['time', 'freq', 'ant', 'dir'])
                weights = reorderAxes(weights, soltab.getAxesNames(), ['time', 'freq', 'ant', 'dir'])             
        else:
            logging.error('No phase or amplitude soltab has been found or specified.')
            return 1

        dimension = np.shape(vals)
        if 'clock' in soltab_name or 'tec' in soltab_name:
            new_vals = np.ndarray(shape=(dimension[0], len(new_station_names)))
            new_weights = np.ndarray(shape=(dimension[0], len(new_station_names)))
        if 'amplitude' in soltab_name or 'phase' in soltab_name:
            if 'pol' in soltab.getAxesNames():
                new_vals = np.ndarray(shape=(dimension[0], dimension[1], len(new_station_names), dimension[3], dimension[4]))
                new_weights = np.ndarray(shape=(dimension[0], dimension[1], len(new_station_names), dimension[3], dimension[4]))
            else:
                new_vals = np.ndarray(shape=(dimension[0], dimension[1], len(new_station_names), dimension[3]))
                new_weights = np.ndarray(shape=(dimension[0], dimension[1], len(new_station_names), dimension[3], dimension[4]))

        for i, new_station in enumerate(new_station_names):
            if new_station in station_names:
                ant_index = soltab.getAxisValues('ant', ignoreSelection=True).tolist().index(new_station)
                if 'clock' in soltab_name or 'tec' in soltab_name:
                    new_vals[:, i] = vals[:, ant_index]
                    new_weights[:, i] = weights[:, ant_index]
                if 'amplitude' in soltab_name or 'phase' in soltab_name:
                    if 'pol' in soltab.getAxesNames():
                        new_vals[:, :, i, :, :] = vals[:, :, ant_index, :, :]
                        new_weights[:, :, i, :, :] = weights[:, :, ant_index, :, :]
                    else:
                        new_vals[:, :, i, :] = vals[:, :, ant_index, :]
                        new_weights[:, :, i, :] = weights[:, :, ant_index, :]
            else:
                if restrictToCS and 'CS' not in new_station:
                    logging.info('RestrictToCS: Omitting station ' + new_station)
                    continue
                logging.info('Adding ' + str(soltab_name) + ' to ' + new_station)
                if 'clock' in soltab_name or 'tec' in soltab_name:
                    new_vals[:, i] = vals[:, STindex]
                    new_weights[:, i] = weights[:, STindex]
                if 'amplitude' in soltab_name or 'phase' in soltab_name:
                    if 'pol' in soltab.getAxesNames():
                        print (new_vals.shape)
                        new_vals[:, :, i, :, :] = vals[:, :, STindex, :, :]
                        new_weights[:, :, i, :, :] = weights[:, :, STindex, :, :]
                    else:
                        new_vals[:, :, i, :] = vals[:, :, STindex, :]
                        new_weights[:, :, i, :] = weights[:, :, STindex, :]

        if 'clock' in soltab_name:
            new_soltab = OutSolset.makeSoltab(soltype='clock', soltabName=soltab_name,
                                              axesNames=['time', 'ant'],
                                              axesVals=[soltab.time, new_station_names],
                                              vals=new_vals, weights=new_weights)
            pass
        elif 'tec' in soltab_name:
            new_soltab = OutSolset.makeSoltab(soltype='tec', soltabName=soltab_name,
                                              axesNames=['time', 'ant'],
                                              axesVals=[soltab.time, new_station_names],
                                              vals=new_vals, weights=new_weights)
        elif 'amplitude' in soltab_name:
            if 'pol' in soltab.getAxesNames():
                new_soltab = OutSolset.makeSoltab(soltype='amplitude', soltabName=soltab_name,
                                                  axesNames=['time', 'freq', 'ant', 'dir', 'pol'],
                                                  axesVals=[soltab.time, soltab.freq, new_station_names, soltab.dir, soltab.pol],
                                                  vals=new_vals, weights=new_weights)
            else:
                new_soltab = OutSolset.makeSoltab(soltype='amplitude', soltabName=soltab_name,
                                                  axesNames=['time', 'freq', 'ant', 'dir'],
                                                  axesVals=[soltab.time, soltab.freq, new_station_names, soltab.dir],
                                                  vals=new_vals, weights=new_weights)

        elif 'phase' in soltab_name:
            if 'pol' in soltab.getAxesNames():
                new_soltab = OutSolset.makeSoltab(soltype='phase', soltabName=soltab_name,
                                                  axesNames=['time', 'freq', 'ant', 'dir', 'pol'],
                                                  axesVals=[soltab.time, soltab.freq, new_station_names, soltab.dir, soltab.pol],
                                                  vals=new_vals, weights=new_weights)
            else:
                new_soltab = OutSolset.makeSoltab(soltype='phase', soltabName=soltab_name,
                                                  axesNames=['time', 'freq', 'ant', 'dir'],
                                                  axesVals=[soltab.time, soltab.freq, new_station_names, soltab.dir],
                                                  vals=new_vals, weights=new_weights)

        soltab = 0
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Copies calibration values from a reference station to a set of new stations, e.g., from the superterp (ST001) to all core stations.')

    parser.add_argument('h5parm', type=str,
                        help='H5parm to which this action should be performed .')
    parser.add_argument('MSfiles', type=str,
                        help='MS for which the new solset shall be created.')
    parser.add_argument('--solset_in', type=str, default='sol000',
                        help='Input solution set')
    parser.add_argument('--solset_out', type=str, default='sol001',
                        help='Output solution set (has to be different from input solution set)')
    parser.add_argument('--soltab_list', '--soltab_list', type=str, default='phase000,amplitude000',
                        help='Comma-separated list of soltabs to be copied')
    parser.add_argument('--superstation', type=str, default='ST001',
                        help='Reference station from which data should be copied')
    parser.add_argument('--restrictToCS',
                        help='Restrict the copy action to core stations only',
                        action='store_true', dest="restrictToCS")
    parser.add_argument('--match_ptg',
			help='Match the pointing direction of the MS, for transferring between directions',
			action='store_true', dest='match_ptg')
    parser.add_argument('--useh5coordinates', help='Copy the source coordinate from the input .h5 file at sol000', default=True)

    args = parser.parse_args()

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s", "%Y-%m-%d %H:%M:%S")
    format_file = logging.Formatter("%(asctime)s %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)

    MSfiles = args.MSfiles
    h5parmfile = args.h5parm

    main(h5parmfile, MSfiles, solset_in=args.solset_in, solset_out=args.solset_out, soltab_list=args.soltab_list, superstation=args.superstation, restrictToCS=args.restrictToCS, matchPtg=args.match_ptg, useh5coordinates=args.useh5coordinates)

