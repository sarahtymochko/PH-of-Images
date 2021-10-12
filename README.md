# PH-of-Images

This repo contains code used in the paper:

T. Heiss, S. Tymochko, B. Story, A Garin, H. Bui, B. Bleile, and V. Robins. The impact of finite resolution on the persistent homology of images. In submission, 2021.


### Code

The following examples (some of which are used in the paper) are:
1. Concentric Rings - `Rings.ipynb`
2. Disk with Spots - 

The majority of the code is contained in the two .py files:
- generate_examples.py contains code to generate many example images 
    - `create_circles_image`: 2 dimensional disk (of any resolution)
    - `make_checkerboard`: Chessboard (of any resolution with any number of squares) 
    - `make_rings_example`: three concentric rings example from the paper (of any resolution)
    - `create_diskwithspot`: A disk with one spot in the middle (fixed resolution of 2048x2048)
    - `create_diskwithmanyspots`: A disk with many spots (fixed resolution of 2048x2048)
- utils.py contains functions related to:
    - Image processing: Compute DSEDT (discrete signed Euclidean distance transform), binarize an image based on a threshold, downsample using the averaging technique described in our paper
    - TDA: Use giotto-TDA to compute persistent homology and bottleneck distance. 
    
### Figures

The Figs folder contains figures used in the paper as well as the same plots with respect to pixel size (where the high resolution image has pixel size=1), rather than resolution.

### Citing

If you use code from this repo, please cite our paper for the methodology as:

bibtex here 



### Contact

Feel free to contact me with any questions or issues!

Sarah Tymochko
tymochko@egr.msu.edu

Contact: Sarah Tymochko (tymochko@egr.msu.edu)


# TITLE

 "Insert title name here" (insert link here).

