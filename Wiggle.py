import numpy as np


""" This is the Wiggle class. It houses simulations of AO loops
with adaptive primaries. Go team."""


class Wiggle:
    """ 
    Parameters 
    ----------
    site : str

    diameter : 
    seeing : 

    """

    def __init__(self, site=None, telescope_diameter=8, science_field=30, 
                 actuator_spacing=0.3, readnoise=1, controller_frequency=30,
                 wavelength=100, zenith_angle=5, fitting_parameter=0.3, 
                 guide_star_mag=4, dm='primary', aperture='square', 
                 coherence_cell_size = 0.1, mean_wind_speed = 10, 
                 cn_squared = 0.1, r0_scalar = 0.1):      

        # Site population
        sites = {'hamilton': {'r0': 10}}            
        if site is not None:
            self.site = site
            self.coherence_cell_size = sites[site]['coherence_cell_size']
            self.mean_wind_speed = sites[site]['mean_wind_speed']
            self.cn_squared = sites[site]['cn_squared']
            self.r0_scalar = sites[site]['r0_scalar']
        

        # Overwriting/calculating site parameters -- FIXME
        self.coherence_cell_size = coherence_cell_size
        self.mean_wind_speed = mean_wind_speed # meters/2
        self.cn_squared = cn_squared
        self.r0_scalar = r0_scalar

        # AO system parameters
        self.telescope_diameter = telescope_diameter # in meters
        self.science_field = science_field # in arcsec
        self.actuator_spacing = actuator_spacing # in meters 
        self.sigma_readnoise = readnoise # in photons/pix
        self.controller_frequency = controller_frequency
        self.wavelength = wavelength 
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

        self.unit_actuators = self.telescope_diameter / self.actuator_spacing

    def _calculate_spatial_frequency_cutoff(self):
        """ Calculates the highest spatial freq (smallest size) that's controllable. """

        self.spatial_cutoff = ( self.unit_actuators / self.telescope_diameter)/2 

    def _calculate_r0(self):
        """ Calcualtes r0 at the zenith angle. """

        self.r0 = self.r0_scalar * (np.cos(self.zenith_angle))**(3/5)

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


test = Wiggle(site=None, telescope_diameter=8, science_field=30,
        actuator_spacing=0.3, readnoise=.01, controller_frequency=30, 
        dm='primary', aperture='square', coherence_cell_size = 1,
        mean_wind_speed = 3, cn_squared = .5, r0_scalar = .02)
        
print(test.strehl)

