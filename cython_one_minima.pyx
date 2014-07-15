#cython: boundscheck=False
#cython: profile=False
#cython: cdivision = True
import numpy as np
cimport numpy as np
cimport cython


DTYPE1 = np.uint16
ctypedef np.uint16_t DTYPE_t
DTYPE2 = np.float
ctypedef np.float DTYPE_t2
DTYPE3 = np.uint32
ctypedef np.uint32_t DTYPE_t3

# fit simulated and experimental minima to calculate thickness
# problem: speed!! --> one main problem is the access of the large list sim_waves --> can this be written as an array?

cdef float _abs(float a): return a if a>=0 else -a


cdef list peakdetect(y_axis, x_axis = None, unsigned short lookahead_min=5, unsigned short lookahead_max=3, unsigned short delta = 0):
    
    # define output container
    #cdef list max_peaks=[]
    cdef list min_peaks = []
    cdef list dump = [] # used to pop the first hit which almost always is false
    cdef list y_axis_list = y_axis.tolist() # convert array to list, min() is faster for list

    # check input data --> this makes the algorithm 5 times slower
    #x_axis, y_axis = _datacheck_peakdetect(x_axis, y_axis) 
    
    # store data length for later use
    cdef unsigned int length = len(y_axis)

    #perform some checks
    #if lookahead < 1:
    #    raise ValueError, "Lookahead must be '1' or above in value"
    #if not (np.isscalar(delta) and delta >= 0):
    #    raise ValueError, "delta must be a positive number"

    #maxima and minima candidates are temporarily stored in 
    #mx and mn respectively
    cdef int mn = 100000
    cdef int mx = -100000
    cdef unsigned short x,y, index

    for index, (x,y) in enumerate(zip(x_axis[:-lookahead_min], y_axis[:-lookahead_min])):
        
        if y > mx:
            mx = y
            mxpos = x

        if y < mn:
            mn = y
            mnpos = x

        #### look for max ####
        
        if y < mx-delta and mx != 100000:
            #Maxima peak candidate found
            # lool ahead in signal to ensure that this is a peak and not jitter
            if max(y_axis_list[index:index+lookahead_max]) < mx:
                #max_peaks.append([mxpos, mx])
                dump.append(True)
                #set algorithm to only find minima now
                
                mx = 100000
                mn = 100000
                if index+lookahead_min >= length:
                    #end is within lookahead no more peaks can be found
                    break
                continue

        #### look for min ####    
        
        if y > mn+delta and mn != -100000:
            #Minima peak candidate found
            # look ahead in signal to ensure that this is a peak and not jitter
            if min(y_axis_list[index:index+lookahead_min]) > mn:
                min_peaks.append(mnpos)
                dump.append(False)
                #set algorithm to only find maximum now
                mn = -100000
                mx = -100000
                if index+lookahead_min >= length:
                    #end is within lookahead no more peaks can be found
                    break


    #Remove the false hit on the first value of the y_axis
    if len(dump)>0:
        if not dump[0]:
            min_peaks.pop(0)
        else:
            pass    
    
        #no peaks were found, should the function return empty lists?
    
    return min_peaks

    
# find reflection minima for every pixel

def c_Fit_Pixel(unsigned int start,unsigned int ende, np.ndarray[DTYPE_t, ndim=3] data,list waves, unsigned short lookahead_min,unsigned short lookahead_max, unsigned short delta, unsigned int one_minima, np.ndarray[DTYPE_t, ndim=2] min_reference):
    cdef unsigned int Image_width = len(data[0][0])
    cdef np.ndarray[DTYPE_t, ndim=2] thickness_ready = np.zeros((ende-start,Image_width),np.uint16 )
    cdef unsigned short spalte, zeile
    cdef np.ndarray[DTYPE_t, ndim=1] intensity
    cdef np.ndarray[double,ndim=1] minima_exp
    cdef unsigned int counter=start, result_minima


    print 'x ', len(data) 
    print 'y ', len(data[0])
    print 'z ', len(data[0][0])
    #print "using no thickness limits"
    for zeile in range(len(data[0])):
        print counter
        counter+=1
        for spalte in xrange(Image_width):
            intensity = data[:,zeile, spalte]
            minima_exp = np.array(peakdetect(intensity, waves, lookahead_min,lookahead_max, delta),dtype=np.float)
            if len(minima_exp)>1:
                if min_reference[zeile][spalte] != 0:
                    result_minima = minima_exp[abs(minima_exp - min_reference[zeile][spalte]).argmin()]
                else:
                    result_minima = 0
            else:
                result_minima = 0
            thickness_ready[zeile][spalte] = result_minima 


    return thickness_ready