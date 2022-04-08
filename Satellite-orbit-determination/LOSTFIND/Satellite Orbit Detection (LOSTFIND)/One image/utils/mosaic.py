import numpy as np
import time
from PIL import Image
from filtre import *

def get_raw_image(filename):
    print('passed 2')
    global_start = time.time()
    DATAPATH = './Processed/ImgCompressed/'
    filename2 = DATAPATH+filename
    raw_img = Image.open(filename)

    global_end = time.time()
    seconds = (global_end - global_start)
    raw_img = filtre(raw_img)
    global_end = time.time()
    seconds = (global_end - global_start)
    if seconds > 10:
        print('t2a > 10s: Error')
        return
    else: print('passed t2a')
    print('im_size = ', raw_img.size)
    compressed_img = raw_img.resize((1200,800))  ############ ERROR #############
    global_end = time.time()
    seconds = (global_end - global_start)
    if seconds > 10:
        print('t3a > 10s: Error')
        return
    compressed_img.save(filename2[:-4]+'compressed.jpg')
    global_end = time.time()
    seconds = (global_end - global_start)
    if seconds > 10:
        print('t4a > 10s: Error')
        return

    raw_img = np.asarray(raw_img)
    rgb_weights = [0.2989, 0.5870, 0.1140]
    raw_img = np.dot(raw_img[...,:3],rgb_weights)
    return raw_img


