from analyse_pic_from_folder import *
import utils.prologue as prologue
from satellite_trace_functions import *
from end_process import *
from utils.create_param_img import create_param_img

def main(args):
    os.chdir(args.folder)

    
    try:
        os.makedirs("./Processed/Lines")
    except OSError as error:
        pass
    try:
        os.makedirs("./Processed/Beginnings")
    except OSError as error:
        pass
    try:
        os.makedirs("./Processed/Middles")
    except OSError as error:
        pass
    try:
        os.makedirs("./Processed/Endings")
    except OSError as error:
        pass
    try:
        os.makedirs("./Processed/Sets")
    except OSError as error:
        pass
    try:
        os.makedirs("./Processed/ImgCompressed")
    except OSError as error:
        pass
    try:
        os.makedirs("./Processed/ImgLines")
    except OSError as error:
        pass
    try:
        os.makedirs("./Processed/ImgCenter")
    except OSError as error:
        pass
    try:
        os.makedirs("./Processed/IOD")
    except OSError as error:
        pass 
    '''
    #os.makedirs("./Processed/Lines")
    os.makedirs("./Processed/Beginnings")
    os.makedirs("./Processed/Middles")
    os.makedirs("./Processed/Endings")
    os.makedirs("./Processed/Sets")
    os.makedirs("./Processed/ImgCompressed")
    os.makedirs("./Processed/ImgLines")
    os.makedirs("./Processed/ImgCenter")
    os.makedirs("./Processed/IOD") '''
        
    print('Parameters are read')
    time, agvert, aghoriz, bouss, stnb, ststatut, texpos, focal, pixsize = read_param("parametres.txt")
    create_param_img(time, args.folder)

    print('Image Analysis')
    #analyse_pic_from_folder(args) # analyse all images in a Folder + save lines, middles, output imgages
    filenames = read_img_name("param_img.txt") #names of all images in the Folder
    print('End of the image analysis')

    print('Detection of the satellites')
    satellites = data_for_el_from_folder(filenames)

    print('Creation of the IOD files')
    IODs = creation_IOD_files(args, satellites)

    print('Creation of one IOD document')
    one_IOD(IODs)

    #print('Lancement de Elfind')
    #path = os.getcwd()
    #Run_elfind_on_sats(path)

if __name__ == '__main__':
    main(prologue.get_args())