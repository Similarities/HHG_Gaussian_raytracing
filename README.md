# High Harmonic Gaussian Beam Raytracing
Two-lens treatment of a gaussian beam with high harmonic frequency (1/N in wavelength), via
Transfer-Matrix (q-vector, only rational part) treatment.

High harmonic generation e.g. via laser beams is developed for
high photon coherent sources with sub-fs pulse durations.
Commonly a laser with a fundamental wavelength is down focused into
/ onto a medium (gas or solid) where in the focus by the interaction
high harmonic radiation with mulitple frequencies of the fundamental freuqency is
emitted.

For the special case of solid dense - laser high harmonic creation, 
the interaction is based on the Relativistic Mirror Model, where a 
sufficient strong laser creates a plasma on the solid surface with which
the rest of the laser pulse interacts. Electrons in the density gradient
of the plasma surface are accelerated by the laser and form an 
relativistic oscillating dense surface upon which the laser is
reflected. By the doppler-shift the reflected light is once per
laser cycle upshifted to a higher photon energy (increase of the reflected
pulse bandwidth to a quasi continuum like distribution). This
creates a pulse train over the driving laser pulse, where each pulse
train member has a significant smaller pulse duration (< 0.5 laser cycle).
On a detector with a integration time that is much longer the driving
laser pulse duration (electronic switching can NOT resolve fs or as!!!)
each pulse of the pulse train is integrated, that by the uncertainty
relation leads to an interference on the detector (cf. SPIDER measurment
method) - where the interference can be described as modulation with
the laser cycle duration (once per cycle) that in the spectrum is found
as high harmonics of the fundamental frequency (or 1/N in wavelength).

This raytracing script tries to model high-harmonic gaussian beam creation close
and in the laser focal region. It focuses only on the spatial characteristics
of such a created beam. 
The created HHG are diverging in space after creation proportional
to the driving laser frequency *1/N (theoretical model). In case of
a plasma surface denting, the source of creation changes the HHG divergence 
like a curved mirror - in other words focus the beam. 
The mathematical treatment uses as a two lens - ray transfer matrix, 
where lens 1 is the initial driving laser focal lens and lens 2 is
the curved plasma surface, that can be either dented outwards or inwards. 
The model treats the case for a displacement of the plasma surface (target)
from the laser focuse, that leads to a bigger spotsize on the target as
well as a reduction of the intensity. The denting depth is depending on
both parameters scaling differently with either defocusing or additionaly
scaling the intensity by a pulse duration enlargement (cf. chirped pulses).
The model is based on Vincenti et. al, NatComm. 2014, but enlarged his model
by the defocusing range.

The script focuses on real experimental parameters, which are known for
the highfield laser community: M^2, Rayleigh-Length, hole boring velocity, 
w0 for gaussian beams, Denting depth of a plasma surface, scaling laws. 


Python 3.6 ... mathematical details will follow sooon

