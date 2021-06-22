import numpy as np


""" This is the Wiggle class. It houses simulations of AO loops
with adaptive primaries. Go team.


EASTER EGGGGGGGGGS
------------------
- cn_squared is a vertical turbulence profile
-- we need it to calculate the structure function
-- we don't know how or why quite yet oops

- coherence cell size vs r0_scalar
-- this are the same
- NO LASER GUIDESTARS
-- NO TIPTILT ERRORS
-- NO ANISOPLANETISM 
- NO SHIRT NO SHOES NO LASER GUIDE STARS

"""


class Wiggle:
    """ 
    Parameters 
    ----------
    site : str
        Key for the site if a prebuilt set of conditions.
    telescope_diameter : float
        Diameter of the primary mirror in meters.
    science_wavelength : float
        Wavelength at which we do our ssscience in nanometers.
    coherence_cell_size : float
        Turbulence pocket scale length -- i.e. r0 NOT dependent 
        on zenith angle in meters.
    zenith_angle : float
        Angle at which we're observing -- assumed to be zero in degrees. 
    mean_wind_speed : float
        Mean speed of the wind in meters/s.
    controller_frequency : float
        Frequency at which the controller (i.e. AO system) can operate in Hz.
    cn_squared : is rude  
        Vertical turbulence profile. #FIXME
    science_field : float #FIXME unused
        Size of science field in arcseconds
    number_actuators : int
        Number of actuators on the adaptive mirror. 
    readnoise : float
        Readnoise of the wavefront sensor in electrons/pix.
    fitting_parameter : float # FIXME
        WHO IS SHE. (Between 0-1 -- set with ext reference from Hardy.)
    guide_star_mag : float
        Magnitude of star in band of WFS. 
    aperture : str #FIXME
        Sets propagator term something something equation something somethign
        WFS. Noll, JOSA 68, 1978
    dm : str
        Whether the DM is primary, secondary, or post focal plane.
    
    """


    def __init__(self, site=None, telescope_diameter=8, science_wavelength=1000, 
                 coherence_cell_size=0.1, zenith_angle=0, mean_wind_speed = 10,
                 controller_frequency=50, cn_squared=1, science_field=30, 
                 number_actuators=559, readnoise=1, fitting_parameter=0.18, 
                 guide_star_mag=4, aperture='square', dm='primary'):
        
        self.parameter_wavlength = 500 # DON'T CHANGE MEEEEEE
        # Site population
        sites = {'hamilton': {'r0': 10}}            
        if site is not None:
            self.site = site
            self.coherence_cell_size = sites[site]['coherence_cell_size'] # in meters
            self.mean_wind_speed = sites[site]['mean_wind_speed'] # m/s
            self.cn_squared = sites[site]['cn_squared'] # 
        

        # Overwriting/calculating site parameters -- FIXME
        self.coherence_cell_size = coherence_cell_size
        self.mean_wind_speed = mean_wind_speed # meters/2
        self.cn_squared = cn_squared

        # AO system parameters
        self.telescope_diameter = telescope_diameter # in meters
        self.science_field = science_field # in arcsec
        self.number_actuators = number_actuators # actuators
        self.sigma_readnoise = readnoise # in photons/pix
        self.controller_frequency = controller_frequency
        self.wavelength = science_wavelength # nanometers
        self.dm=dm
        self.aperture=aperture
        
        self.zenith_angle = zenith_angle # specify in radians 
        self.guide_star_mag = guide_star_mag 
        self._calculate_r0()
        self._calculate_structure_function()
        self._calculate_greenwood_frequency()
        
        self._calculate_diffraction_limit()
        self._calculate_actuators_across()
        self._calculate_spatial_frequency_cutoff()

        self.fitting_parameter = fitting_parameter
        
        self._calculate_fitting_error()
        self._calculate_measurement_error()
        self._calculate_anisoplatanism_error()
        self._calculate_bandwidth_error()
        self._calculate_high_order_wfe()
        self._calculate_strehl()
        
        
        # Detector parameters 


        #self.science_camera = 
        #self.wfs = 

        #self._calculate_influence_functions()
    

    def _calculate_structure_function(self):
        """ Calculates the structure function for a given site."""
        self.structure_function = self.cn_squared*(self.r0**(2/3))

    def _calculate_diffraction_limit(self):
        """ Calculates diffraction limit. """

        pass
    
    def _calculate_actuators_across(self):
        """ Calculates actuators across the diameter. I.e., 'unit actuators'. """
        
        # spherical mirrors only
        self.actuator_spacing =(self.telescope_diameter/2)**2*np.pi/self.number_actuators 
        self.unit_actuators = self.telescope_diameter / self.actuator_spacing

    def _calculate_spatial_frequency_cutoff(self):
        """ Calculates the highest spatial freq (smallest size) that's controllable. """

        self.spatial_cutoff = ( self.unit_actuators / self.telescope_diameter)/2 

    def _calculate_r0(self):
        """ Calcualtes r0 at the zenith angle. """

        self.r0 = self.coherence_cell_size * (np.cos(self.zenith_angle))**(3/5)

    def _calculate_greenwood_frequency(self):
        """ Calculates greendwood frequency for a site/target."""

        self.greenwood_frequency = 0.426 * (self.mean_wind_speed / self.r0)
   
    def _calculate_fitting_error(self):
        
        sigma_fitting_squared = self.fitting_parameter * (( (self.actuator_spacing)/(self.r0) )**(5/3)) *  (self.wavelength /(2*np.pi))
        self.sigma_fitting = np.sqrt(sigma_fitting_squared)

    def _calculate_measurement_error(self):
        """ Calculates measurement error for the AO system. """    
        
        # Calculate Hartmann Spot
        arc_seconds_across = 206265*5.89e-7
        term1, term2 = arc_seconds_across/self.actuator_spacing, arc_seconds_across/self.r0
        hartmann_spot = np.max([term1, term2])
        
        # Calculate SNR 
        n_pix=4
        sample_time = 1/(10*self.controller_frequency)
        n_photons = self.guide_star_mag*sample_time*((100*self.actuator_spacing)**2)
        snr = n_photons/np.sqrt(n_photons + n_pix*(self.sigma_readnoise)**2)
        
        # Calculate noise propagator 
        degrees_of_freedom = np.round( (np.pi/4) * (self.telescope_diameter/self.actuator_spacing) )
        propagation, factor = 0.0536, 0.0795  # FIXME
        if self.aperture == 'circular':
            proagation, factor = 0.0068, 0.0796
        noise_propagator = np.sqrt(2*(propagation + factor*np.log(degrees_of_freedom)))
        
        # Calculate close loop averaging
        controller_over_frame = 1/10
        close_loop_averaging = np.sqrt(2*controller_over_frame)*np.arctan(1/(2*controller_over_frame))

        sigma_measurement = noise_propagator * close_loop_averaging * self.actuator_spacing * (hartmann_spot / snr)
        self.sigma_measurement = sigma_measurement

    def _calculate_anisoplatanism_error(self):
        """ Calculates anisoplatanism error. Not important unless we're using a laser guide star. """

        self.sigma_anisoplatanism = 0
    
    def _calculate_bandwidth_error(self):
        """ Calculates bandwidth error. """

        self.sigma_bandwidth = (self.greenwood_frequency / self.controller_frequency)**(5/3)

    def _calculate_high_order_wfe(self):

        self.high_order_wfe = np.sqrt( self.sigma_fitting**2 + self.sigma_measurement**2 + 
                                       self.sigma_anisoplatanism**2 + self.sigma_bandwidth**2 )

    def _calculate_strehl(self):
        """ Calculates strehl given high order wfe. """

        self.strehl = np.exp( (0 - (2*np.pi / self.wavelength) * self.high_order_wfe)**2)


    def build_science_detector(self):

        self.science_detector = hcipy.NoiselessDectector()
        #.... 

    
    def _calculate_influence_functions(self):
        if self.dm == 'primary':
            pass

    
    def calculate_static_performance(self, verbose=False):
        
        self.calculate_measurement_error()


    def build_deformable_mirrr(self):

        self.deformable_mirror = hicpy.DeformableMirror(self.influence_functions)


    def run_simulation(self, verbose=False):

        # running through time???
        self.build_science_detector()
        self.build_wfs()
        self.build_deformable_mirror()
        
        # GO
        
    def plot_metrics(self):
        pass

        # idk plot some things


#shane_2_4 = Wiggle(site='hamilton', diameter=2.4)
#shane.run_simulation()

#actuators vs diameter:

if __name__ == "__main__":
    test = Wiggle(site=None, telescope_diameter=8, science_field=30,
        actuator_spacing=0.3, readnoise=.01, controller_frequency=30, 
        dm='primary', aperture='square', coherence_cell_size = 1,
        mean_wind_speed = 3, cn_squared = .5, r0_scalar = .02)
        
    print(test.strehl)

