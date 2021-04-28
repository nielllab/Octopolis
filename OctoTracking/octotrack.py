"""
octotrack.py
"""
import argparse, json, sys, os, subprocess, shutil, yaml
import cv2
import pandas as pd
import numpy as np
os.environ["DLClight"] = "True"
import deeplabcut

from utils.paths import find, list_subdirs
from utils.angles import get_angles, get_center, get_rotation

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, help='config options for the current session as a .yaml file')
    args = parser.parse_args()
    return args

def main(config_path):
    # open config file as dict
    with open(config_path, 'r') as infile:
        config = yaml.load(infile, Loader=yaml.FullLoader)
    # iterate through the videos
    session = list_subdirs(config['session_dir'])
    # iterate through each video
    for recording in session:
        # get video .mov and time .txt
        vidpath = find(recording+'.mov', config['session_dir'])
        timepath = find(recording+'.txt', config['session_dir'])
        # get props of the video filepath
        base_path, video_name = os.path.split(vidpath)
        # read in screen time data
        with open(timepath) as f:
            time_txt = f.read()
        screendict = json.loads(time_txt)
        # from the dictionary of frame numbers, create bool list for each frame where an entry is True whenever the simulus is moving to the right
        switch_frames = sorted(screendict['start_left'] + screendict['start_right'])
        first_onset_right = [True if np.min(switch_frames) in screendict['start_right'] else False][0]
        switches = np.digitize(range(num_frames), switch_frames) % 2 == 1 # True for all instances of the first type of stim presented
        if first_onset_right is True:
            # if the first stimulus was moving right, 'switches' will already label all right movements as True
            moving_right = switches
        else:
            # otherwise, flip the list so that right movement is still labeled as True
            moving_right = ~switches
        # run deeplabcut on the video
        deeplabcut.analyze_videos(config['pose_estimation']['network'], vidpath)
        # find and read in the DLC .h5 file that was saved out from deeplabcut.analyze_videos
        octodata = pd.read_hdf(path)
        # drop bad pts
        likcols = [i for i in octodata.columns.values if 'likelihood' in i]
        xcols = [i for i in octodata.columns.values if '_x' in i]
        ycols = [i for i in octodata.columns.values if '_y' in i]
        # iterate through each DLC point (each has three columns: x,y,likelihood)
        # set points below threshold to NaN
        for indnum in len(likcols):
            likcol = likcols[indnum]; xcol = xcols[indnum]; ycol = ycols[indnum]
            octodata[xcol][octodata[likcol]<config['tracking']['lik_thresh']] = np.nan
            octodata[ycol][octodata[likcol]<config['tracking']['lik_thresh']] = np.nan
        # add the direction of gratings to data as a new column
        octodata = pd.concat([octodata, pd.Series(moving_right).rename('moving_right')])
        # get angle of each eye
        octodata = get_angles(octodata, config['angles']['leye_names'])
        octodata = get_angles(octodata, config['angles']['reye_names'])
        # get centerpoint of each eye
        octodata = get_center(octodata, config['angles']['leye_names'])
        octodata = get_center(octodata, config['angles']['reye_names'])
        # get angle of head (angle between eyes)
        octodata = get_angles(octodata, config['angles']['head_names'])
        # get centerpoint of head (centerpoint between the center points of both eyes)
        octodata = get_center(octodata, config['angles']['head_names'])
        # get angle of body (angle between mantle and centerpoint between eyes)
        octodata = get_angle(octodata, config['angles']['body_names'])
        # get octo rotation around the center of the tube
        octodata = get_rotation(octodata)
        # save out as .h5 file in the recording directory
        save_path = os.path.join(base_path, 'tracking.h5')
        octodata.to_hdf(save_path, 'w')

if __name__ == '__main__':
    args = get_args()
    main(args.config)