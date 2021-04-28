# OctoTracking

## Installation
1. Be sure that GPU drivers, CUDA, etc. are set up on the device.
2. Install anaconda., etc.
3. Clone the DeepLabCut repository locally, which can be found [on github](https://github.com/DeepLabCut/DeepLabCut).
4. Create the GPU version of the DeepLabCut anaconda environment, instructions [here](https://github.com/DeepLabCut/DeepLabCut/blob/master/conda-environments/README.md). Be sure to install it with `tensorflow-gpu`.
5. Run `pip install -r requirements.txt` in the directory of the FreelyMovingEphys repository to install any other dependencies.

## Data and Organization
### Raw Data
Before any analysis is run, the session directory should be structured so that there is a seperate directory for each .mov video, into which other outputs can be saved.

- session
    - recording1
        - recording1.mov
    - recording2
        - recording2.mov

### Output Files
A .txt file of stimulus times will be saved into each recording directory once `stimlabel.py` is run for that .mov file. Once `octotrack.py` is run on a session, each recording directory for that session will have a .h5 and .pickle written in it storing DeepLabCut pose estimates, as well as an .h5 file named 'tracking.h5', storing a single dataframe with raw points, eye and head angles, and animal angular velocity.

## Usage
Use the DeepLabCut environment (with added requirements), activated with `conda activate DLC-GPU2`, if 'DLC-GPU2' is the name of the built environment.

### Labeling Stimulus Directions
To launch the stimulus-labeling GUI, run `python stimlabel.py --count 5` if there are five videos to label in sequence.

To label the onset of leftwards movement on the stimulus, press 'q', and 'p' for rightwards movement. To reload the GUI (if it stops responding), press 'c'. To quit, press 'esc'.

### Octopus Tracking
To run octopus tracking, run `python octotrack.py --config /path/to/config.yaml`, where the argument for config is the custom .yaml for this specific session, **not** the DeepLabCut .yaml. The path to the data session is provided through this config, along with the path to the DeepLabCut config.yaml file for the network trained on the octopus tracking task. The names of the DLC points to expect, and how to use them, should be specified by the `angles` section of the .yaml. An example .yaml file can be found in this repository [here](/OctoTracking/config.yaml).

### Post-Tracking Analysis
TO DO:
* in jupyter, explore correlations between octopus activities with gratings direction
* func to write a .pdf of figures, once they're polished