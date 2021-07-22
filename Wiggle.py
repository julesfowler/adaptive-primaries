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
    isoplanatic_angle : float
        Max separation while moving through the same turbulence pocket in
        arcseconds. 
    science_object_separation : float
        Separation between science object and guidestar in arcseconds.
    mean_wind_speed : float
        Mean speed of the wind in meters/s.
    controller_frequency : float
        Frequency at which the controller (i.e. AO system) can operate in Hz.
    cn_squared : is rude  
        Vertical turbulence profile. #FIXME
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
                 coherence_cell_size=0.1, zenith_angle=0, isoplanatic_angle=3.87,
                 science_object_separation=0, mean_wind_speed = 10,
                 controller_frequency=50, cn_squared=1, 
                 actuator_spacing=0.3, readnoise=1, fitting_parameter=0.18, 
                 guide_star_mag=4, aperture='square', dm='primary'):
        
        self.parameter_wavelength = 500 # DON'T CHANGE MEEEEEE
        # Site population
        sites = {'hamilton': {'r0': 10}}            
        if site is not None:
            self.site = site if not hasattr(self, 'site') else self.site
            self.coherence_cell_size = sites[site]['coherence_cell_size'] # in meters
            self.mean_wind_speed = sites[site]['mean_wind_speed'] # m/s
            self.cn_squared = sites[site]['cn_squared'] # 
        

        # Overwriting/calculating site parameters -- FIXME
        self.coherence_cell_size = coherence_cell_size if not hasattr(self, 'coherence_cell_size') else self.coherence_cell_size
        self.mean_wind_speed = mean_wind_speed if not hasattr(self, 'mean_wind_speed') else self.mean_wind_speed # meters/2
        self.cn_squared = cn_squared if not hasattr(self, 'cn_squared') else self.cn_squared

        # AO system parameters
        self.telescope_diameter = telescope_diameter if not hasattr(self, 'telescope_diameter') else self.telescope_diameter# in meters
        self.actuator_spacing = actuator_spacing if not hasattr(self, 'actuator_spacing') else self.actuator_spacing # meters
        self.sigma_readnoise = readnoise if not hasattr(self, 'readnoise') else self.readnoise # in photons/pix
        self.controller_frequency = controller_frequency if not hasattr(self, 'controller_frequency') else self.controller_frequency 
        self.science_wavelength = science_wavelength if not hasattr(self, 'science_wavelength') else self.science_wavelength # nanometers
        self.dm = dm if not hasattr(self, 'dm') else self.dm 
        self.aperture = aperture if not hasattr(self, 'aperture') else self.aperture 
        
        self.zenith_angle = zenith_angle if not hasattr(self, 'zenith_angle') else self.zenith_angle 
        self.isoplanatic_angle = isoplanatic_angle if not hasattr(self, 'isoplanatic_angle') else self.isoplanatic_angle 
        self.science_object_separation = science_object_separation if not hasattr(self, 'science_object_separation') else self.science_object_separation
        self.guide_star_mag = guide_star_mag if not hasattr(self, 'guide_star_mag') else self.guide_star_mag 
        
        self.fitting_parameter = fitting_parameter if not hasattr(self, 'fitting_parameter') else self.fitting_parameter
        
        self.calculate_ao_error_terms()
        self.calculate_thickness_terms()

        
    def _calculate_structure_function(self):
        """ Calculates the structure function for a given site."""
        
        self.structure_function = self.cn_squared*(self.r0**(2/3))

    def _calculate_diffraction_limit(self):
        """ Calculates diffraction limit. """
        
        # return in mas
        self.diffraction_limit = 2.063e8 * self.science_wavelength/self.telescope_diameter
    
    def _calculate_actuators_across(self):
        """ Calculates actuators across the diameter. I.e., 'unit actuators'. """
        
        self.unit_actuators = self.telescope_diameter / self.actuator_spacing

    def _calculate_spatial_frequency_cutoff(self):
        """ Calculates the highest spatial freq (smallest size) that's controllable. """
        
        # unused 
        self.spatial_cutoff = (self.unit_actuators / self.telescope_diameter)/2 

    def _calculate_r0(self):
        """ Calcualtes r0 at the zenith angle. """

        self.r0 = self.coherence_cell_size * (np.cos(np.deg2rad(self.zenith_angle)))**(3/5)

    def _calculate_greenwood_frequency(self):
        """ Calculates greendwood frequency for a site/target."""
        
        # Where does this factor come from (pulled from spreadsheet) #FIXME
        self.greenwood_frequency = 0.426 * (self.mean_wind_speed / self.r0)
   
    def _calculate_fitting_error(self):
        
        sigma_fitting_squared = self.fitting_parameter*(((self.actuator_spacing)/(self.r0))**(5/3))
        self.sigma_fitting = np.sqrt(sigma_fitting_squared)*(self.parameter_wavelength /(2*np.pi))

    def _calculate_measurement_error(self):
        """ Calculates measurement error for the AO system. """    
        
        # Calculate Hartmann Spot
        # FIXME what are factor_1, factor_2 ???
        factor_1, factor_2 = 206265*5.89e-7, 206265*6.5e-7
        term1, term2 = factor_1/self.actuator_spacing, factor_2/self.r0
        hartmann_spot = np.max([term1, term2])
        
        # Calculate SNR 
        n_pix=4 # FIXME spreadsheet says not to change this idk why?
        sample_time = 1/(10*self.controller_frequency)
        brightness = (8.9e5)*10**((0-self.guide_star_mag)/2.5)
        n_photons = brightness*sample_time*((100*self.actuator_spacing)**2)
        snr = n_photons/np.sqrt(n_photons + n_pix*(self.sigma_readnoise)**2)

        # Calculate noise propagator 
        degrees_of_freedom = np.round((np.pi/4) * (self.telescope_diameter/self.actuator_spacing)**2)
        factor_1, factor_2 = 0.0536, 0.0795  # FIXME WHAT THE HECK IS THIS
        if self.aperture == 'circular':
            factor_1, factor_2 = 0.0068, 0.0796
        noise_propagator = np.sqrt(2*(factor_1 + factor_2*np.log(degrees_of_freedom)))

        # Calculate close loop averaging
        controller_over_frame = 1/10
        close_loop_averaging = np.sqrt(2*controller_over_frame)*np.arctan(1/(2*controller_over_frame))
        sigma_measurement = noise_propagator * close_loop_averaging * (self.actuator_spacing*1e9) * (hartmann_spot/snr*4.84814e-6)
        self.sigma_measurement = sigma_measurement # in nm

    def _calculate_anisoplatanism_error(self):
        """ Calculates anisoplatanism error. Not important unless we're using a laser guide star. """

        self.sigma_anisoplatanism = np.sqrt((self.science_object_separation/self.isoplanatic_angle)**(5/3))*(self.parameter_wavelength /(2*np.pi))
    
    def _calculate_bandwidth_error(self):
        """ Calculates bandwidth error. """

        self.sigma_bandwidth = np.sqrt((self.greenwood_frequency / self.controller_frequency)**(5/3))*(self.parameter_wavelength /(2*np.pi))

    def _calculate_high_order_wfe(self):

        self.high_order_wfe = np.sqrt( self.sigma_fitting**2 + self.sigma_measurement**2 + 
                                       self.sigma_anisoplatanism**2 + self.sigma_bandwidth**2 )

    def _calculate_strehl(self):
        """ Calculates strehl given high order wfe. """

        self.strehl = np.exp(-1*((2*np.pi/self.science_wavelength)*self.high_order_wfe)**2)

    def _calculate_matching_thickness(self):
        """ Calculates the thickness required to match Keck ASM specs."""

        asm_spacing = 0.039 # in meters
        asm_thickness = 3.55e-3 # in meters

        self.matching_thickness = asm_thickness*(self.actuator_spacing/asm_spacing)**2
        self.thickness = self.matching_thickness

    def _calculate_driven_mass(self):
        """ Calculates the mass each actuator needs to drive given a borofloat
        mirror."""

        borofloat_density = 2230 # in kg/m^3
        actuator_region = self.actuator_spacing**2 # for square regions
        self.driven_mass = borofloat_density*self.thickness*actuator_region

    def _calculate_resonant_frequency(self):
        """ Calculates the resonant frequency for a mirror. """
    
        spring_constant =  0.525e6 # in N/m #FIXME
        self.resonant_frequency = np.sqrt(spring_constant/self.driven_mass)
    
    def _calculate_influence_functions(self):
        if self.dm == 'primary':
            pass

    def calculate_quilting_error(self, resistance, thickness):
        """ Calculates gravity quilting error."""

        self.quilting_error = ((0.126*self.actuator_spacing**2)**2)/(resistance*thickness**2)
    
    def calculate_min_thickness(self, resistance, max_quilting_error):
        """ Calculates the thickness required to limit to a certain quilting error."""

        self.min_thickness = (0.126*self.actuator_spacing**2)/np.sqrt(resistance*max_quilting_error)
        self.thickness = self.min_thickness

    def calculate_ao_error_terms(self):
        """ Runs the calculations for AO error terms. """

        self._calculate_r0()
        self._calculate_structure_function() # WHO IS SHE??
        self._calculate_greenwood_frequency()
        
        self._calculate_diffraction_limit() # in mas
        self._calculate_actuators_across()
        self._calculate_spatial_frequency_cutoff()

        self._calculate_fitting_error()
        self._calculate_measurement_error()
        self._calculate_anisoplatanism_error()
        self._calculate_bandwidth_error()
        self._calculate_high_order_wfe()
        self._calculate_strehl()
    
    def calculate_thickness_terms(self):
        """ Calculates the terms related to facesheet thickness. """

        self._calculate_matching_thickness()
        self._calculate_driven_mass()
        self._calculate_resonant_frequency()
    
    def recalculate(self):
        """ Recalculates the AO error and thickness terms when we update
        params."""

        self.calculate_ao_error_terms()
        self._calculate_driven_mass()
        self._calculate_resonant_frequency()


if __name__ == "__main__":
    test = Wiggle()
        
    print(test.strehl)
    print(test.resonant_frequency/test.controller_frequency)
