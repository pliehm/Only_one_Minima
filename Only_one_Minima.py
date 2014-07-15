###############################################################################
### Script to calculate the minima position closest to a reference position ###
###############################################################################

#############
### INPUT ###
#############


# enter all the folders which contain subfolders with the image files, 
# the folder with the image files MUST be in another folder (its name has to be entered here)
# e.g. you have 5 folders with images: 1,2,3,4,5 --> all these 5 folders have to be e.g. in a folder
# named "data". "data" is what you would enter in the list below. You can enter more then one folder 
# in this list (e.g. for 3 different long-time measurementes)

data = ['Test_Irma']

# enter file name of the reference file

# chose wavelength range and step-width

wave_start = 600    # [nm]
wave_end = 700     # [nm]

# enter position to which the closest minimum should be found

one_minima = 650

# do you want to have a file header or not? Enter True if not, useful to import in ImageJ

no_header = True

# define parameters for minima detection  

lookahead_min = 5 # something like peak width for the minima
# delta = 7    # something like peak height


# enter True if you want to enable this smoothing
x_y_smooth = False
# enter sigma for the gaussian smoothing
x_y_sigma = 0.1

# parameters for printing
# color map is calculated like (mean_thickness - color_min, mean_thickness + color_max) 

color_min = 500
color_max = 500

# enter True if you want to enable this smoothing
x_y_smooth = False
# enter sigma for the gaussian smoothing
x_y_sigma = 0.1



############################
### END of INPUT SECTION ###
############################



#############################
#### start of the program ###
#############################

import cython_one_minima as Fit # import self-written cython code
import numpy as np
import time
import os 
import Image as im
from scipy import ndimage
import multiprocessing as mp
import matplotlib.pyplot as plt

version = 'BioCalc 2.0'

t_a_start = time.time() # start timer for runtime measurement

if __name__ == '__main__':

    for data_folder in data:
        folder_list = os.listdir(data_folder)
        folder_list.sort()
        folder_counter = 0
        for folder in folder_list:


            # enter number of cpu cores, this has to be an integer number!
            # number of physical cores is a good start, but you can try with a larger as well

            multi_p = True   # True for multiprocessing, False for single core (Windows)
            cores = 8

            
            # enter name of simulation_file

            #sim_file = 'Sim_0.5Cr_25Ag_50SiO2_Elastomer_RT601_25Au_500_760nm.txt'

            lookahead_max = lookahead_min-1 # for the maxima --> should not be larger than lookahead_min

            # make wavelength list

            wave_step = 1       # [nm]
            
            waves=[]

            waves=[wave_start + i*wave_step for i in xrange((wave_end-wave_start)/wave_step + 1)]

            ## read image data 
            dateien=os.listdir(data_folder+'/'+folder)
            dateien.sort()
            
            # get size, bit-depth of the images
            for i in range(len(dateien)):
                if dateien[i][-5:]=='.tiff' or dateien[i][-4:]=='.tif':
                    Img=im.open(data_folder + '/'+folder + '/' + dateien[i])
                    break
            Image_width = Img.size[0]
            Image_height = Img.size[1]
            Image_mode = Img.mode 
            if Image_mode == 'RGB' or Image_mode == 'P':
                Image_bit = '8'
            elif Image_mode == 'I;16' or Image_mode == 'I;12':
                Image_bit = '16'
            else:
                print 'Image mode is:', Img.mode
                Image_bit = '16B'
                #testttt = input('Unknown Image mode, ask Philipp for help, sorry for the inconvenience! Abort the programm with CTRL+C or CTRL+Z, maybe several times')     
            #generates an empty array --> image grey values 
            alle=np.zeros(((wave_end-wave_start)/wave_step + 1,Image_height,Image_width),np.uint16)
            alle_x_y_smooth = np.zeros(((wave_end-wave_start)/wave_step + 1,Image_height,Image_width),np.uint16)

            # create array for minima reference with zeros if first folder in calculation
            if folder_counter == 0:
                min_reference = np.zeros((Image_height,Image_width),np.uint16)
                min_reference = min_reference + one_minima
            folder_counter += 1

            # define function to convert the image-string to an array
            def image2array(Img):
                newArr= np.fromstring(Img.tostring(), np.uint16)
                newArr= np.reshape(newArr, (Image_height,Image_width))


            # read every image in folder and check if it is in the wavelength range --> 
            # write grey values into array
            
            counter=0
            print 'reading images from folder: ', folder
            if x_y_smooth == True:
                print 'smoothing of every wavelength image with sigma = ', x_y_sigma, ' before thickness calculation'
            else:
                print 'no smoothing of the wavelength images'    
            for i in xrange(len(dateien)):
                if dateien[i][-5:]=='.tiff' or dateien[i][-4:]=='.tif':
                    if int(dateien[i][:3]) >= wave_start and int(dateien[i][:3]) <= wave_end:
                        #print dateien[i]
                        #print counter
                        Img=im.open(data_folder + '/'+folder + '/' + dateien[i])
                        if Image_bit == '8':
                            Img = Img.convert('L')

                        alle[counter]=np.asarray(Img)
                        # smoothing x-y direction
                        if x_y_smooth == True:
                            Img_s = ndimage.gaussian_filter(Img, sigma=x_y_sigma)
                            alle_x_y_smooth[counter] = np.asarray(Img_s)
                        counter+= 1
    ##################################
    ##### Section for smoothing ######
    ##################################
            lambda_smooth = False
            lambda_sigma = 1
            if x_y_smooth == True:
                alle = alle_x_y_smooth.copy()

            if lambda_smooth == True:
                print 'smoothing over wavelength with sigma = ', lambda_sigma
                for zeile in range(Image_height):
                    #print zeile
                    for spalte in range(Image_width):
                        alle[:,zeile,spalte] = ndimage.gaussian_filter1d(alle[:,zeile,spalte],lambda_sigma)


    #####################################
    ##### END section for smoothing #####
    #####################################


            print 'perform the calculations'
            t1 = time.time()

####################################################
## Find new delta to deal with different dynamics ##
####################################################
               
            # use different delta scaling for 8,12,16 bit
            if Image_bit == '16':
                local_min = np.where(alle[:50]<=alle[:50].min())
                if len(local_min[0])>1:
                    local_min_temp = []
                    local_min_temp = [local_min[0][0],local_min[0][1],local_min[0][2]]
                    local_min = local_min_temp
                local_profile = alle[:,int(local_min[1]),int(local_min[2])][:50]
                new_delta = (local_profile.max() - local_profile.min())/2
                if new_delta > 7:
                    delta = new_delta

            if Image_bit == '8':
                delta = 7    


##############################################################
## Start calculations either with multiprocessing or single ##
##############################################################


            # define queue for the multiprocessing

            if multi_p == True:

                def put_into_queue(start,ende,que,alle, waves, lookahead_min, lookahead_max, delta, one_minima, min_reference):

                    que.put(Fit.c_Fit_Pixel(start,ende,alle, waves, lookahead_min, lookahead_max, delta,one_minima,min_reference)) # calls the C-Fit-function
                    #print 'Schlange ist satt'

                
                
                # devide the rows by the core-number --> to split it equally, assing the rest to the last process
                Zeile_Teil = Image_height/cores
                Zeile_Rest = Image_height%cores

                # start multiprocessing with queues

                Prozesse = []
                Queues = []

                for i in range(cores):
                    Queues.append(mp.Queue())

                for i in range(cores):
                    if i < cores-1:
                        Prozesse.append(mp.Process(target=put_into_queue,args=(i*Zeile_Teil,(i+1)*Zeile_Teil,Queues[i],alle[:,(i*Zeile_Teil):((i+1)*Zeile_Teil),:], waves, lookahead_min, lookahead_max, delta,one_minima, min_reference[(i*Zeile_Teil):((i+1)*Zeile_Teil),:])))
                    if i == cores-1:
                        Prozesse.append(mp.Process(target=put_into_queue,args=(i*Zeile_Teil,(i+1)*Zeile_Teil+Zeile_Rest,Queues[i],alle[:,(i*Zeile_Teil):((i+1)*Zeile_Teil),:], waves, lookahead_min, lookahead_max, delta,one_minima, min_reference[(i*Zeile_Teil):((i+1)*Zeile_Teil),:])))
                for i in range(cores):
                    Prozesse[i].start()
                    

                # initialise array for thicknesses
                dicke = np.ndarray((0,Image_width),dtype=np.uint16)

                for i in range(cores):
                    #print 'queuet', i
                    dicke = np.append(dicke,Queues[i].get(),axis=0)

                for i in range(cores):
                    #print 'joint', i
                    Prozesse[i].join()

            if multi_p == False:
                start = 0
                ende = Image_height

                dicke = Fit.c_Fit_Pixel(start,ende,alle, waves, lookahead_min, lookahead_max, delta, one_minima, min_reference)
            t2 = time.time()

            # assigne latest calculation to reference
            min_reference = dicke

            print t2-t1, 'seconds just for the calculation'

            # count not fitted values

            not_fitted = Image_height*Image_width - np.count_nonzero(dicke)
            not_fitted_percent = 100.0/(Image_height*Image_width)*not_fitted
            print 'not fitted values',not_fitted
            print 'in percent:', not_fitted_percent

            print 'write data to file'
            # use numpy function to save array to file, '0' and not '-' used for missing values
            HEADER = time.strftime('Version = ' + version + '\n' + "%d.%m.%Y at %H:%M:%S")+'\n' + 'folder with data = ' + folder + '\n' + 'wave_start = '+str(wave_start) + '\n' + 'wave_end = ' + str(wave_end) + '\n' + 'lookahead_min = ' + str(lookahead_min) + '\n'  + 'lookahead_max = ' + str(lookahead_max) + '\n' + 'delta = ' + str(delta) + ' delta was varied +-5'+ '\n' +  '\n' +  'not fitted values: ' + str(not_fitted) + ', percentage of whole image: ' + str(not_fitted_percent)  + '\n'
            if x_y_smooth == True:
                HEADER+= 'x_y_smoothing done with sigma = ' + str(x_y_sigma) + '\n'
            if lambda_smooth == True:
                HEADER+= 'lambda smoothing done with sigma = ' + str(lambda_sigma) + '\n'

            HEADER+= '\n'

            # If filename should be different for calculations with smoothing
            # if x_y_smooth == True and lambda_smooth == False:
            #     file_name = data_folder + '_' +folder + time.strftime("_%Y%m%d_%H%M%S")+'_x_y_sigma_' + str(x_y_sigma) + '_smoothed.txt'
                
            # if lambda_smooth == True and x_y_smooth == False:
            #     file_name = data_folder + '_' + folder + time.strftime("_%Y%m%d_%H%M%S")+'_lambda_sigma_' + str(lambda_sigma) + '_smoothed.txt'

            # if lambda_smooth == True and x_y_smooth == True:
            #     file_name = data_folder + '_' + folder +  time.strftime("_%Y%m%d_%H%M%S")+'_x_y_sigma_' + str(x_y_sigma) + '_lambda_sigma_' +str(lambda_sigma) + '_smoothed.txt'

            # if lambda_smooth == False and x_y_smooth == False: 
            #     file_name = data_folder + '_' + folder + time.strftime("_%Y%m%d_%H%M%S")+'.txt'

            if no_header == True:
                file_name = data_folder + '_' + folder + time.strftime("_%Y%m%d_%H%M%S")+'.txt'
                np.savetxt(data_folder + '/' + file_name,dicke,fmt='%d')

            if no_header == False:
                file_name = data_folder + '_' + folder + time.strftime("_%Y%m%d_%H%M%S")+'.txt'
                np.savetxt(data_folder + '/' + file_name,dicke,fmt='%d',header=HEADER )

            ### script to replace a certain string or string combination in a file

            #t_replace_1 = time.time()

            # do it twice because 0 0 would not be - - but - 0 since the empty character before the second 0 has already been read
            for i in range(2):
                p = open(data_folder + '/' + file_name,'r')
                if no_header == True and i == 0:
                    string = '\n'
                    string += p.read()

                p.close()
                if no_header == True and i == 1:
                   string = string[string.find('\n')+1:]
                string = string.replace(' 0 ', ' - ')
                string = string.replace('\n'+'0 ', '\n'+'- ')
                string = string.replace(' 0'+'\n', ' -'+'\n')

                p = open(data_folder + '/' + file_name,'w')
                p.write(string)

                p.close()
            #print (time.time()-t_replace_1), ' seconds for replaceing 0 with -'
            print (time.time()-t_a_start), ' seconds for the whole program'

            # plot datta
            
            plt.figure(folder)
            plt.imshow(dicke)
            plt.clim(dicke.mean()-color_min,dicke.mean()+color_max)
            plt.colorbar()
            #plt.show()



