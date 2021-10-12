import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.morphology import distance_transform_edt
from sklearn.pipeline import Pipeline, make_pipeline, FeatureUnion, make_union
from gtda.homology import CubicalPersistence
from gtda.diagrams import PairwiseDistance
from skimage.measure import block_reduce
import skimage

from utils import make_binary, SEDT

def create_diskwithmanyspots():
    data = np.ones([2048,2048])
    bigdisk = create_circles_image(510*2,r=5.00001)
    data[2048//2-510:2048//2+510, 2048//2-510:2048//2+510] = bigdisk
    smalldisk = create_circles_image(5*2,r=5.1)

    for j in [i*85+5+75//2+2 for i in range(1,2048//85)]:
        for k in [i*85+5+75//2+2 for i in range(1,2048//85)]:
    #         if j > 554 and j < 1489 and k > 554 and k < 1489:
    #             if (j == 639 and (k == 724 or k == 1319)) or (j == 724 and (k == 639 or k == 1404)):
    #                 continue
    #             if (j== 1404 and (k == 724 or k == 1319)) or (j == 1319 and (k == 639 or k == 1404)):
    #                 continue
            if len(np.unique(data[j-5:j+5, k-5:k+5])) == 1:
                if np.unique(data[j-5:j+5, k-5:k+5])[0] == 0:
                    data[j-5:j+5, k-5:k+5] = data[j-5:j+5, k-5:k+5]+ -1*(smalldisk-1)

    data = -1*(data-1)
    
    return data

def create_diskwithspot():
    data = np.zeros([2048,2048])
    circle = -1*(make_binary(create_circles_image(1024, use_SEDT=False)) -1)
    data[2048//4: 3*2048//4, 2048//4:3*2048//4] = circle
    data[2048//2-10:2048//2+10, 2048//2-10:2048//2+10] = create_circles_image(20, r=5,use_SEDT=False)
    return data

def create_synthetic_image(num_pixels, use_SEDT = True):
    '''
    Creates image with specified number of pixels.
    Optionally thresholds image and uses SEDT to create greyscale image.
    '''
    
    def f(x, y):
        '''
        Creates a 2d image from sines and cosines
        '''
        return np.sin(x) ** 10 + np.cos(10 + y * x) * np.cos(x) 

    x1 = np.linspace(0, 10, num_pixels)
    y1 = np.linspace(0, 10, num_pixels)
    X1, Y1 = np.meshgrid(x1, y1)

    x2 = np.linspace(-5, 5, num_pixels)
    y2 = np.linspace(-5, 5, num_pixels)
    X2, Y2 = np.meshgrid(x2, y2)

    
    
    Z = f(X1, Y1) + f(X2, Y2)

    if use_SEDT:
        Z = SEDT(make_binary(Z))

    return Z

def create_circles_image(num_pixels, r=4,use_SEDT = True):
    '''
    Creates image with specified number of pixels.
    Optionally thresholds image and uses SEDT to create greyscale image.
    
    '''    
    
    x1 = np.linspace(0, 10, num_pixels)
    y1 = np.linspace(0, 10, num_pixels)
    X1, Y1 = np.meshgrid(x1, y1)
        
    Z = (X1-5)**2 + (Y1-5)**2 - r**2
    

    if use_SEDT:
        Z = SEDT(make_binary(Z))

    return make_binary(Z)

def create_circles_3D_image(num_pixels, r=4, use_SEDT = True):
    '''
    Creates image with specified number of pixels.
    Optionally thresholds image and uses SEDT to create greyscale image.
    '''    
    
    x1 = np.linspace(0, 10, num_pixels)
    y1 = np.linspace(0, 10, num_pixels)
    z1 = np.linspace(0, 10, num_pixels)
    X1, Y1, Z1 = np.meshgrid(x1, y1, z1)
        
    Z = (X1-5)**2 + (Y1-5)**2 + (Z1-5)**2 - r**2
    

    if use_SEDT:
        Z = SEDT(-1*(make_binary(Z)-1))

    return Z


def make_rings_example(num_pixels):

    Z1 = make_annulus(num_pixels, 0.1,0.15)*-1
    Z2 = make_annulus(num_pixels, 1,1.2)*-1
    Z3 = make_annulus(num_pixels, 3,3.5)*-1

    img = Z1+Z2+Z3
    HR_binary = np.copy(img)
    
    return HR_binary


def make_annulus(num_pixels, rmin, rmax, use_SEDT=False):
    x1 = np.linspace(0, 10, num_pixels)
    y1 = np.linspace(0, 10, num_pixels)
    X1, Y1 = np.meshgrid(x1, y1)
        
    Z1 = (X1-5)**2 + (Y1-5)**2 - rmin**2
    Z2 = (X1-5)**2 + (Y1-5)**2 - rmax**2
    Z = make_binary(Z2)-make_binary(Z1)
    
    if use_SEDT:
        Z = SEDT(Z)
        
    return Z



def make_checkerboard(onesqsize,nsqs=8):
    board = np.zeros([onesqsize*nsqs, onesqsize*nsqs])
    
    for i in range(nsqs):
        for j in range(nsqs):
            if i%2 == 0 and j%2 == 0:
                board[i*onesqsize:(i+1)*onesqsize, j*onesqsize:(j+1)*onesqsize] = np.ones([onesqsize,onesqsize])
            elif i%2 == 1 and j%2 == 1:
                board[i*onesqsize:(i+1)*onesqsize, j*onesqsize:(j+1)*onesqsize] = np.ones([onesqsize,onesqsize])
                
    return board