"""
angles.py
"""
import pandas as pd
import numpy as np

def get_center(df, xy_names):
    """
    create new columns in dataframe for the center of existing points
    concretely, run this to get the center-point of each eye from the front and back pts
    of that eye, and the center of the head from those eye centerpoints
    INPUTS
        df: DeepLabCut .h5 file read in as a pandas dataframe (each row is a timepoint, each column is a labeled point)
        xy_names: dict of positions (i.e. 'x1') mapped to points labeled in the DLC network (i.e. 'front_left')
    OUPUTS
        pandas dataframe with new columns for the centerpoints in x and y
    """
    df[xynames['base_name']+'_cent_x'] = df.apply(lambda x: np.mean(x[xy_names['x1']], x[xy_names['x2']]))
    df[xynames['base_name']+'_cent_y'] = df.apply(lambda x: np.mean(x[xy_names['y1']], x[xy_names['y2']]))
    return df

def get_angles(df, xy_names)
    """
    get the angles of all frames in a video given a dictionary of pts to compare
    used for either eye
    INPUTS
        df: DeepLabCut .h5 file read in as a pandas dataframe (each row is a timepoint, each column is a labeled point)
        xy_names: dict of positions (i.e. 'x1') mapped to points labeled in the DLC network (i.e. 'front_left')
    OUTPUTS
        pandas dataframe with new column for the angle between points at each frame (in radians) 
    """
    df[xynames['base_name']+'_ang'] = df.apply(lambda x: x[xy_names['x1']] - x[xy_names['x2']])
    return df

def get_rotation(df):
    """
    get the angular velocity in rad/sec of animal using head angle at each frame and the known frame rate
    INPUTS
        df: dataframe with required column, 'head_ang', as angle of the animal head (angle between eyes) in radians
    OUTPUTS
        same as df, but with added column for head angular velocity
    """
    # angular velocity
    df['head_angular_velocity'] = df['head_ang'].diff().fillna(0).apply(lambda x: x / (1/config['tracking']['frame_rate']))
    return df