from analyse_pic_from_folder import *
import utils.prologue as prologue
from satellite_trace_functions import *
from end_process import *
from utils.create_param_img import create_param_img

def main(args):
    os.chdir(args.folder)

    try:
        os.makedirs("./Processed/Lines")
        os.makedirs("./Processed/Middles")
        os.makedirs("./Processed/Sets")
        os.makedirs("./Processed/ImgCompressed")
        os.makedirs("./Processed/ImgLines")
        os.makedirs("./Processed/ImgCenter")
        os.makedirs("./Processed/IOD")
    except OSError as error:
        pass
        
    print('Lecture des parametres')
    times, agvert, aghoriz, bouss, stnb, ststatut, texpos, focal, pixsize = read_param("parametres.txt")
    create_param_img(times,args.folder)

    print('Image analysis')
    analyse_pic_from_folder(args) # analyse all images in a Folder + save lines, middles, output imgages
    filenames = read_img_name("param_img.txt") #names of all images in the Folder
    print('End of image analysis')

    print('Satellite detection in the different images')
    satellites = data_for_el_from_folder(filenames)

    print('IOd files are being created')
    IODs = creation_IOD_files(args, satellites)

    print('Creation of a complete IOD file')
    one_IOD(IODs)

if __name__ == '__main__':
    main(prologue.get_args())