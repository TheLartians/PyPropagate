# PyPropagate

A paraxial wave propagation framework for python 2.7 . 
Created by Lars Melchior for the [Institute for X-Ray Physics](http://www.roentgen.physik.uni-goettingen.de/) in Göttingen.
Funded by [SFB755](http://www.uni-goettingen.de/de/318955.html).

## Installation
    
    pip install pypropagate

## Dependencies

- Python 2.7
- gcc >= 5.0 or compatible
- boost python
- xraylib

Python dependencies are installed automatically via pip.

## Licence

PyPropagate is licenced under version 3 of the GNU General Public License.
Third party contributions are strongly encouraged.

## Citing PyPropagate

When using PyPropagate in your work it would be very appreciated to give us a citation. You can use the following BibTeX templates:

    @article{Melchior:17,
        author = {Lars Melchior and Tim Salditt},
        journal = {Opt. Express},
        keywords = {Pulses; Ultrafast phenomena; X-rays, soft x-rays, extreme ultraviolet (EUV); Propagation; Waves},
        number = {25},
        pages = {32090--32109},
        publisher = {OSA},
        title = {Finite difference methods for stationary and time-dependent X-ray propagation},
        volume = {25},
        month = {Dec},
        year = {2017},
        url = {http://www.opticsexpress.org/abstract.cfm?URI=oe-25-25-32090},
        doi = {10.1364/OE.25.032090},
        abstract = {We have generalized finite-difference (FD) simulations for time-dependent field propagation problems, in particular in view of ultra-short x-ray pulse propagation and dispersion. To this end, we first derive the stationary paraxial (parabolic) wave equation for the scalar field envelope in a more general manner than typically found in the literature. We then present an efficient FD implementation of propagators for different dimensionality for stationary field propagation, before we treat time-dependent problems by spectral decomposition, and suitable numerical sampling. We prove the validity of the numerical approach by comparison to analytical theory, using simple tractable propagation problems. Finally, we apply the framework to the problem of modal dispersion in X-ray waveguide. We show that X-ray waveguides can be considered as non-dispersive optical elements down to sub-femtosecond pulse width. Only when considering resonant absorption close to an X-ray absorption edge, we observe pronounced dispersion effects for experimentally achievable pulse widths. All code used for the work is made available as supplemental material.},
    }

