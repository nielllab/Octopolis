"""
stimlabel.py
"""
import tkinter as tk
from tkinter import filedialog
import cv2, os, json

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--count', type=int, help='number of videos to label')
    args = parser.parse_args()
    return args

def nothing(x):
    pass

def main(count):
    for step in range(count):
        # let user choose the video path
        root = tk.Tk()
        root.withdraw()
        video_path = filedialog.askopenfilename()
        # get save path for later use, and the name of the video
        base_path, video_name = os.path.split(video_path)
        # set up to play video
        cap = cv2.VideoCapture(video_path)
        frame_num = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
        cv2.namedWindow(video_name,cv2.WND_PROP_FULLSCREEN)
        cv2.setWindowProperty(video_name,cv2.WND_PROP_FULLSCREEN,cv2.WINDOW_FULLSCREEN)
        # tracking bars
        cv2.createTrackbar('frame',video_name,0,frame_num-1,nothing)
        font = cv2.FONT_HERSHEY_DUPLEX
        fontcol = (0,255,0); fontsize = 0.5; lastdisp = ''
        # set up lists of reversal times
        left_onset = []; right_onset = []
        # open window for labeling
        while(1):
            # get current positions of trackbar
            frame_pos = cv2.getTrackbarPos('frame',video_name)
            cap.set(cv2.CAP_PROP_POS_FRAMES, frame_pos)
            # make sure it read in correctly
            ret, frame = cap.read()
            # recording if any keys were pressed
            k = cv2.waitKey(1)
            # exit
            if k == 27: # key 27 is 'esc'
                break
            # reload GUI if it freezes
            elif k == 99: # key 99 is 'c'
                cv2.destroyAllWindows()
                cv2.namedWindow(video_name,cv2.WND_PROP_FULLSCREEN)
                cv2.setWindowProperty(video_name,cv2.WND_PROP_FULLSCREEN,cv2.WINDOW_FULLSCREEN)
                cv2.createTrackbar('frame',video_name,0,frame_num-1,nothing)
            # gratings started moving left
            elif k == 113: # key 113 is 'q'
                left_onset.append(frame_pos)
            # gratings started moving right
            elif k == 112: # key 112 is 'p'
                right_onset.append(frame_pos)
            # maybe add navigation through time at some point? they're short videos, so probably not needed?
            cv2.imshow(video_name,frame)
        cv2.destroyAllWindows()
        # put the reversal times into a dict
        timing_dict = {'right_onset':right_onset, 'left_onset': left_onset}
        # build the save path for the txt file
        save_path = os.path.join(base_path, os.path.splitext(video_name)+'.txt')
        with open(save_path, 'w') as f:
            f.write(json.dumps(timing_dict))

if __name__ == '__main__':
    main(args.count)