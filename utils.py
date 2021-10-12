import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.morphology import distance_transform_edt, distance_transform_cdt
from sklearn.pipeline import Pipeline, make_pipeline, FeatureUnion, make_union
from gtda.homology import CubicalPersistence
from gtda.diagrams import PairwiseDistance
from skimage.measure import block_reduce
import skimage


def compute_pers(image, homology_dimensions=[0,1], dgm_format='ripser'):
    '''
    Parameters:
        image - np.array
        homology_dimensions - list of dimensions to compute (default is [0,1])
        dgm_format - string specifying what format you want persistence diagrams returned
            - 'ripser' - list of arrays like the output of ripser
            - 'giotto' - array where first two columns are birth and death and third column is the dimension
    '''
    
    cubicalpers = CubicalPersistence(homology_dimensions)
    diagrams = cubicalpers.fit_transform([image])
    
    if dgm_format == 'ripser':
        diagrams = reformat_dmgs(diagrams)[0]
    
    return diagrams

def calc_bottleneck(dgm_a, dgm_b, delta=0):
    '''
    Calculate the bottleneck distance between persistence diagrams using giotto
    
    Parameters:
        dgm_a - persistence diagram in giotto format (note only ONE PH dimension should be included)
        dgm_b - persistence diagram in giotto format (note only ONE PH dimension should be included)
        delta - approximation parameter for bottleneck computation
    '''
    
    bndst = PairwiseDistance(metric='bottleneck', metric_params={'delta':delta}) # if delta=0 then bottleneck distance is exact
    
    dim = int(dgm_a[0,2])
    
    if len(dgm_a)<len(dgm_b):
        dif = len(dgm_b) - len(dgm_a)
        dgm_a = np.concatenate([dgm_a, dim*np.ones([dif, 3])])
    elif len(dgm_b)<len(dgm_a):
        dif = len(dgm_a) - len(dgm_b)
        dgm_b = np.concatenate([dgm_b, dim*np.ones([dif, 3])])
    
    d = bndst.fit_transform([dgm_a, dgm_b])
    
    if np.mean(d) == 0:
        dis = 0
    else: 
        d[d==0] = np.nan
        dis = np.nanmean(d[d!=0])
    
    return dis

def make_binary(img, threshold=0):
    '''
    Makes a binary image by thresholding input image above specified threshold
    '''
    return 1*(img > threshold)

def SEDT(im):
    '''
    Pass in binary image where white (1 pixels) gets positive values and black (0 pixels) get negative values
    '''
    # white pixels
    dt0 = distance_transform_edt(im==0)

    # black pixels
    dt1 = distance_transform_edt(im==1)

    return dt1 - dt0


def reformat_dmgs(dgms):
    '''
    Transform giotto format for persistence diagrams into ripser style to use plotting
    
    '''
    new_dgms = []
    # List of multiple diagrams
    for dgm in dgms:
        new_ = []
        # Diagram dimension
        for j in np.unique(dgm[:,2]):
            new = dgm[dgm[:,2] == j]
            
            new_.append( new[:,:2] )
        new_dgms.append(new_)
    return new_dgms


def downsample_average(image, block_size, t, nothreshold=False):
    '''
    Downsample image using averaging method
    
    Parameters:
        image - np.array of dim 2 or 3
        block_size - size of kernel (kernel is always square of size block_size^dim)
        t - threshold to use to make it binary
    '''
    
    dim = len(np.shape(image))

    block_shape = tuple([block_size for _ in range(dim)])
        
    blocks = skimage.util.view_as_blocks(image, block_shape)  
    
    zeros_size = [np.shape(blocks)[i] for i in range(dim)]
    avgs = np.zeros( zeros_size )
    
    if dim == 2:
        for i in range(np.shape(blocks)[0]):
            for j in range(np.shape(blocks)[1]):
                b = np.copy(blocks[i,j])
                avgs[i,j] = np.mean(b) #,where = b!=fillval)
    else:
        for i in range(np.shape(blocks)[0]):
            for j in range(np.shape(blocks)[1]):
                for k in range(np.shape(blocks)[2]):
                    b = np.copy(blocks[i,j,k])
                    avgs[i,j,k] = np.mean(b) #,where = b!=fillval)
                    
    if nothreshold:
        return avgs
            
    th = avgs>t
                       
    return th*1








def largest_sphere_diameter(image, leash_length=1):
    max_kernel = int(np.amax(distance_transform_cdt(image))*2)
    print(max_kernel)
    
    for k in range(max_kernel,1,-1):
        m = find_leash_length(image, k)
        
        if m <= leash_length:
            return k

    print('Cound not find kernel size for given leash length')
    return 1


def find_leash_length(image, kernel, return_dt=False, plot=False):
    dim = len(np.shape(image))
    img_shape = [np.shape(image)[i] for i in range(dim)]

    is_hit = np.zeros(img_shape)
    kernel_block = [kernel for _ in range(dim)]

    if dim == 2:
        for i in range(np.shape(image)[0]):
            for j in range(np.shape(image)[1]):
                if (image[i:i+kernel, j:j+kernel] == 1).all():
                    is_hit[i:i+kernel, j:j+kernel] = np.ones(kernel_block)
    elif dim == 3:
        for i in range(np.shape(image)[0]):
            for j in range(np.shape(image)[1]):
                for j in range(np.shape(image)[2]):
                    if (image[i:i+kernel, j:j+kernel, k:k+kernel] == 1).all():
                        is_hit[i:i+kernel, j:j+kernel, k:k+kernel] = np.ones(kernel_block)
    else:
        print('Error. Data must be 2 or 3 dimensional.')
    
    # Where the kernel DOESNT hit
    dif = image-is_hit
    
    if plot:
        plt.imshow(dif)
    
    # Distance transform 
    dt = distance_transform_edt((1-is_hit)*-1)
    
    dt[dif == 0] = 0
    
    if return_dt:
        return np.amax(dt), dt
    else:
        return np.amax(dt)
    


