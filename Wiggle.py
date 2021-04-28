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

    def __init__(self, telescope=None, diameter=None, seeing=None, dm='primary'):
        
        site = {'hamilton': {'r0': 10}}

        if site is not None:
            self.site = site
            self.r0 = ...
        
        self.r0 = r0 if r0 is not None
        
        self.science_camera = 
        self.wfs = 

        self._calculate_influence_functions()

    def some_math(self):
         
        # MATH MATH MATH

    def build_science_detector(self):

        self.science_detector = hcipy.NoiselessDectector()
        .... 

    
    def _calculate_influence_functions(self):
        if self.dm == 'primary':
            ...

    
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

        # idk plot some things


shane_2_4 = Wiggle(site='hamilton', diameter=2.4)
shane.run_simulation()

#actuators vs diameter:


    
        


