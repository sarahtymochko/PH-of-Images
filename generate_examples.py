import numpy as np
from utils import make_binary, SEDT

def create_diskwithmanyspots():
    '''
    Generate example from paper with one big disk and a grid of smaller disks inside.
    '''
    
    data = np.ones([2048,2048])
    bigdisk = create_circles_image(510*2,r=5.00001)
    data[2048//2-510:2048//2+510, 2048//2-510:2048//2+510] = bigdisk
    smalldisk = create_circles_image(5*2,r=5.1)

    for j in [i*85+5+75//2+2 for i in range(1,2048//85)]:
        for k in [i*85+5+75//2+2 for i in range(1,2048//85)]:
            if len(np.unique(data[j-5:j+5, k-5:k+5])) == 1:
                if np.unique(data[j-5:j+5, k-5:k+5])[0] == 0:
                    data[j-5:j+5, k-5:k+5] = data[j-5:j+5, k-5:k+5]+ -1*(smalldisk-1)

    data = -1*(data-1)
    
    return data

def create_diskwithspot():
    '''
    Generate example with one big disk and one small disk inside
    '''
    data = np.zeros([2048,2048])
    circle = -1*(make_binary(create_circles_image(1024, use_SEDT=False)) -1)
    data[2048//4: 3*2048//4, 2048//4:3*2048//4] = circle
    data[2048//2-10:2048//2+10, 2048//2-10:2048//2+10] = create_circles_image(20, r=5,use_SEDT=False)
    return data


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


def make_rings_example(num_pixels):
    '''
    Make concentric rings example from paper
    
    '''

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