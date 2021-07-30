WIGGLE
======
Wavefront Improvements with GiGantic Light-deforming Experiments 

AKA ... adaptive primaries!

Required Packages
-----------------
- numpy


Setup
-----
```
>>> python setup.py install 
```

Example Use
-----------

```python
from wiggle.Wiggle import Wiggle

ap = Wiggle(site=None, telescope_diameter=1.4, science_wavelength=2500,
                 coherence_cell_size=0.15, zenith_angle=20, isoplanatic_angle=14,
                 science_object_separation=0, mean_wind_speed = 10,
                 controller_frequency=50,
                 actuator_spacing=0.14, readnoise=1, fitting_parameter=0.3,
                 guide_star_mag=8, aperture='square', dm='primary') 
print(ap.strehl)
```
