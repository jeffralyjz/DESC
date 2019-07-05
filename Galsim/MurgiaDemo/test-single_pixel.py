import galsim

# Make blank image
im = galsim.ImageF(31,31, init_value=0)

# Set the coordinates so 0,0 refers to the central pixel.
#im.setCenter(0,0)
im.setCenter(0,0)

# Make 80K photons with position 0,0
photons = galsim.PhotonArray(800000)

    # """The PhotonArray class encapsulates the concept of a collection of photons incident on
    # a detector.
    # A PhotonArray object is not typically constructed directly by the user.  Rather, it is
    # typically constructed as the return value of the `GSObject.shoot` method.
    # At this point, the photons only have x,y,flux values.  Then there are a number of classes
    # that perform various modifications to the photons such as giving them wavelengths or
    # inclination angles or removing some due to fringing or vignetting.
    # TODO: fringing, vignetting, and angles are not implemented yet, but we expect them to
    # be implemented soon, so the above paragraph is a bit aspirational atm.
    # Attributes
    # ----------
    # A PhotonArray instance has the following attributes, each of which is a numpy array:
    # - x,y           the incidence positions at the top of the detector
    # - flux          the flux of the photons
    # - dxdz, dydz    the tangent of the inclination angles in each direction
    # - wavelength    the wavelength of the photons
    # Unlike most GalSim objects (but like Images), PhotonArrays are mutable.  It is permissible
    # to write values to the above attributes with code like
    #     >>> photon_array.x += numpy.random.random(1000) * 0.01
    #     >>> photon_array.flux *= 20.
    #     >>> photon_array.wavelength = sed.sampleWavelength(photonarray.size(), bandpass)
    #     etc.
    # All of these will update the existing numpy arrays being used by the photon_array instance.
    # Note about the flux attribute
    # -----------------------------
    # Normal photons have flux=1, but we allow for "fat" photons that combine the effect of several
    # photons at once for efficiency.  Also, some profiles need to use negative flux photons to
    # properly implement photon shooting (e.g. InterpolateImage, which uses negative flux photons to
    # get the interpolation correct).  Finally, when we "remove" photons, for better efficiency, we
    # actually just set the flux to 0 rather than recreate new numpy arrays.
    # Initialization
    # --------------
    # The initialization constructs a PhotonArray to hold N photons, but does not set the values of
    # anything yet.  The constructor allocates space for the x,y,flux arrays, since those are always
    # needed.  The other arrays are only allocated on demand if the user accesses these attributes.
    # @param N            The number of photons to store in this PhotonArray.  This value cannot be
    #                     changed.
    # @param x            Optionally, the initial x values. [default: None]
    # @param y            Optionally, the initial y values. [default: None]
    # @param flux         Optionally, the initial flux values. [default: None]
    # @param dxdz         Optionally, the initial dxdz values. [default: None]
    # @param dydz         Optionally, the initial dydz values. [default: None]
    # @param wavelength   Optionally, the initial wavelength values. [default: None]
    # """
photons.x = 0.5
photons.y = 0.5
photons.flux = 1.

# Accumulate these photons
sensor = galsim.SiliconSensor(strength=1, diffusion_factor=1, qdist=4, name='lsst_itl_32')


# """--------------------------------------------------

#     A model of a silicon-based CCD sensor that converts photons to electrons at a wavelength-
#     dependent depth (probabilistically) and drifts them down to the wells, properly taking
#     into account the repulsion of previously accumulated electrons (known as the brighter-fatter
#     effect).
#     There are currently three sensors shipped with GalSim, which you can specify as the `name`
#     parameter mentioned below.
#         lsst_itl_8      The ITL sensor being used for LSST, using 8 points along each side of the
#                         pixel boundaries.
#         lsst_itl_32     The ITL sensor being used for LSST, using 32 points along each side of the
#                         pixel boundaries.  (This is more accurate than the lsst_itl_8, but slower.)
#         lsst_etv_32     The ETV sensor being used for LSST, using 32 points along each side of the
#                         pixel boundaries.  (This file is still somewhat preliminary and may be
#                         updated in the future.)
#     The Silicon model is asymmetric in the behavior along rows and columns in the CCD.
#     The traditional meaning of (x,y) is (col,row), and the brighter-fatter effect is stronger
#     along the columns than across the rows, since charge flows more easily in the readout
#     direction.
#     There is also an option to include "tree rings" in the Silicon model, which add small
#     distortions to the sensor pixel positions due to non-uniform background doping in the silicon
#     sensor.  The tree rings are defined by a center and a radial amplitude function.  The radial
#     function needs to be a galsim.LookupTable instance.  Note that if you just want a simple cosine
#     radial function, you can use the helper class method `SiliconSensor.simple_treerings` to
#     build the LookupTable for you.
#     Note that there is an option to transpose the effect if your definition of the image is to
#     have the readout "columns" along the x direction.  E.g. to conform with the LSST Camera
#     Coordinate System definitions of x,y, which are transposed relative to the usual FITS meanings.
#     This only affects the direction of the brighter-fatter effect.  It does not change the meaning
#     of treering_center, which should still be defined in terms of the coordinate system of the
#     images being passed to `accumulate`.
#     @param name             The base name of the files which contains the sensor information,
#                             presumably calculated from the Poisson_CCD simulator, which may
#                             be specified either as an absolute path or as one of the above names
#                             that are in the `galsim.meta_data.share_dir/sensors` directory.
#                             name.cfg should be the file used to simulate the pixel distortions,
#                             and name.dat should containt the distorted pixel information.
#                             [default: 'lsst_itl_8']
#     @param strength         Set the strength of the brighter-fatter effect relative to the
#                             amount specified by the Poisson simulation results.  [default: 1]
#     @param rng              A BaseDeviate object to use for the random number generation
#                             for the stochastic aspects of the electron production and drift.
#                             [default: None, in which case one will be made for you]
#     @param diffusion_factor A factor by which to multiply the diffusion.  Use 0.0 to turn off the
#                             effect of diffusion entirely. [default: 1.0]
#     @param qdist            The maximum number of pixels away to calculate the distortion due to
#                             the charge accumulation. A large value will increase accuracy but
#                             take more time. If it is increased larger than 4, the size of the
#                             Poisson simulation must be increased to match. [default: 3]
#     @param nrecalc          The number of electrons to accumulate before recalculating the
#                             distortion of the pixel shapes. [default: 10000]
#     @param treering_func    A galsim.LookupTable giving the tree ring pattern f(r). [default: None]
#     @param treering_center  A PositionD object with the center of the tree ring pattern in pixel
#                             coordinates, which may be outside the pixel region. [default: None;
#                             required if treering_func is provided]
#     @param transpose        Transpose the meaning of (x,y) so the brighter-fatter effect is
#                             stronger along the x direction. [default: False]

# --------------------------------------------------"""


sensor.accumulate(photons, im)

# """--------------------------------------------------------
# Accumulate the photons incident at the surface of the sensor into the appropriate
#         pixels in the image.
#         Each photon has a position, which corresponds to the (x,y) position at the top of the
#         sensor.  In general, they may also have incidence directions and wavelengths, although
#         these are not used by the base class implementation.
#         The base class implementation simply accumulates the photons above each pixel into that
#         pixel.
#         @param photons      A PhotonArray instance describing the incident photons.
#         @param image        The image into which the photons should be accumulated.
#         @param orig_center  The position of the image center in the original image coordinates.
#                             [default: (0,0)]
#         @param resume       Resume accumulating on the same image as a previous call to accumulate.
#                             In the base class, this has no effect, but it can provide an efficiency
#                             gain for some derived classes. [default: False]
#         @returns the total flux that fell onto the image.
# --------------------------------------------------------"""

areas = sensor.calculate_pixel_areas(im)

# Presumably at this point you have some way to check what the sensor object
# has for the vertex locations.  It's not a public API feature, but if you
# think it would be useful, we could add it.

# Write it out.  All flux is in the center.
im.write('./output/imageitl32nod_800000_1-0_d1.fits')
areas.write('./output/areasitl32nod_800000_1-0_d1.fits')
