# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 09:39:59 2022

@author: Kim Degn Jensen
Spectro Inlets A/S
contact: support@spectroinlets.com
"""

# Call stuff
from pathlib import Path
import ixdat
import os
import glob

#Function aimed at easy data readout
def read_full_tsv_and_mpt_files(data_dir, tsv_function, tsv_reader_name, mpt_reader_name):

    
    # Define a nested function defning a parent directory for which you can and want to compiale multiple EC-MS experiements from. Note function is not selective in its input all .tsv and .mpt files iun the first subfolder layer is included 
    def read_tsv_and_mpt_files_as_experiment_per_directories(rootdir):
        paths=[]    #used to define all paths. has no usable function as of yet, but can be nice to also include MID data in a meaningful manner
        tsv_paths=[]    #used to define all paths for tsv files
        mpt_paths=[]    #used to define all paths for mpt files
        for path in Path(rootdir).iterdir():
            if path.is_dir():
                paths.append(path)
                os.chdir(path)
                for file in glob.glob("*.tsv"):
                    tsv_paths.append(path / file)
                for file in glob.glob("*.mpt"):
                    mpt_paths.append(path / file)    
        return(paths, tsv_paths, mpt_paths)
    
    # Combine tsv files w ixdat tuble of tsv paths using reader e.g. reader_name="Zilien" E.g. read_tsv_and_mpt_files_as_experiment_per_directories(Dir)[1] prodcues usable list of tsv files
    def combine_tsv_files_using_ixdat(tsv_paths, technique_name, reader_name):
        try:
            assert len(tsv_paths) > 0, "Empty list of .tsv files!"
        except AssertionError as msg:
            print(msg)
        combined = ixdat.Measurement.read(tsv_paths[0], technique=technique_name, reader=tsv_reader_name)
        for tsv_path in tsv_paths[1:]:
            combined=combined+ixdat.Measurement.read(tsv_path, technique=technique_name, reader=tsv_reader_name) #reader="zilien" or
        return(combined)
    
    
    # # Combine mpt files w ixdat tuble of mpt paths using reader e.g. reader_name="Zilien" E.g. read_tsv_and_mpt_files_as_experiment_per_directories(Dir)[2] prodcues usable list of mpt files
    def combine_mpt_files_using_ixdat(mpt_paths, reader_name):
        try:
            assert len(mpt_paths) > 1, "Empty list of .mpt files!"
        except AssertionError as msg:
            print(msg)
        combined = ixdat.Measurement.read(mpt_paths[0], reader=mpt_reader_name)
        for mpt_path in mpt_paths[1:]:
            combined=combined+ixdat.Measurement.read(mpt_path, reader=mpt_reader_name) #reader="zilien" or
        return(combined)
    MS_data=combine_tsv_files_using_ixdat(read_tsv_and_mpt_files_as_experiment_per_directories(data_dir)[1], tsv_function, tsv_reader_name)  
    EC_data=combine_mpt_files_using_ixdat(read_tsv_and_mpt_files_as_experiment_per_directories(data_dir)[2], mpt_reader_name)  
    ECMS_data=MS_data+EC_data
    
    return(ECMS_data, EC_data, MS_data)