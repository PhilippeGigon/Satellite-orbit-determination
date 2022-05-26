from photutils.background import Background2D, MedianBackground
from photutils.segmentation import make_source_mask
import numpy as np
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
#from astropy.visualization.mpl_normalize import ImageNormalize
#from astropy.visualization import SqrtStretch
#from astropy.stats import sigma_clipped_stats
#from astropy.stats import mad_std
#from astropy.stats import SigmaClip
import matplotlib.pyplot as plt
plt.style.use(astropy_mpl_style)

filename = '2022-04-27_234051_frame0037.fit'
img = fits.getdata(filename)

mask = make_source_mask(img, 3.0, npixels=100, dilate_size=11)
track = np.ma.array(img, mask = np.logical_not(mask), fill_value = 0)

plt.figure()
plt.imshow(1-mask, origin='lower', cmap='Greys', interpolation='none')
plt.axis('off')
plt.show()
plt.savefig('_noBackground.png')
print("here")

#fits.writeto(filename[:-4] + 'no_background.fits', data = mask, overwrite=True)
#fits.writeto(filename, data=mask, overwrite=True, checksum=True)
