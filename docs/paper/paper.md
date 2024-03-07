---
title: 'PulPy: A Python Toolkit for MRI RF and Gradient Pulse Design'
tags:
  - Python
  - magnetic resonance imaging (MRI)
  - signal processing
  - RF pulse
  - gradient
  - radiofrequency pulse design
  - adiabatic pulse
  - parallel transmission (pTx)
  - multi-slice imaging
  - SLR pulse design
authors:
  - name: Jonathan B. Martin
    orcid: 0000-0002-9384-8056
    equal-contrib: false
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Heng Sun
    equal-contrib: false
    affiliation: 2
  - name: Madison Albert
    equal-contrib: false
    affiliation: 3
  - name: Kevin M. Johnson
    equal-contrib: false
    affiliation: 4
  - name: William A. Grissom
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 5
affiliations:
 - name: Vanderbilt University Institute of Imaging Science, Vanderbilt University Medical Center, Nashville, TN, USA
   index: 1
 - name: Department of Biomedical Engineering, Yale University, New Haven, CT, USA
   index: 2
 - name: Department of Biomedical Engineering, Vanderbilt University, Nashville, TN, USA
   index: 3
 - name: Department of Medical Physics and Radiology, University of Wisconsin School of Medicine and Public Health, Madison, WI, USA
   index: 4
 - name: Department of Biomedical Engineering, Case Western Reserve University, Cleveland, OH, USA
   index: 5
date: 26 Mar 2023
bibliography: paper.bib

---

# Summary
We present PulPy (Pulses in Python), an extensive set of open-source, Python-based tools for magnetic resonance imaging (MRI) radiofrequency (RF) and gradient pulse design. PulPy is a Python package containing implementations of a wide range of commonly used RF and gradient pulse design tools. Our implemented functions for RF pulse design include advanced Shinnar-LeRoux (SLR), multiband, adiabatic, optimal control, B$_1^+$-selective and small-tip parallel transmission (pTx) designers. Gradient waveform design functionality is included, providing the ability to design and optimize readout or excitation k-space trajectories [@Pauly1989].  Other useful tools such as vendor-specific waveform input/output, Bloch equation simulators, abstracted linear operators, and pulse reshaping functions are included.  This toolbox builds on the RF tools introduced previously in the SigPy.RF Python software package [@Martin2020a]. The current toolbox continues to leverage SigPy’s existing capabilities for GPU computation, iterative optimization, and powerful abstractions for linear operators and applications [@Ong2019]. The table below shows an outline of the implemented functions.

| Module               | Description                                                      |
|:--------------------:|:----------------------------------------------------------------:|
| .rf.adiabatic.py     | Adiabatic/frequency-swept RF pulses (e.g. [@Garwood2001])        |
| .rf.b1sel.py         | B1-selective pulses (e.g. [@Martin2022])                         |
| .rf.multiband.py     | Pulses for simultaneous multi-slice (e.g. [@Norris2011])         |
| .rf.optcont.py       | Large tip angle optimal control design (e.g. [@Connolly1986])    |
| .rf.ptx.py           | parallel transmit pulse designers (e.g. [@Grissom2006])          |
| .rf.shim.py          | parallel transmit RF shimming (e.g. [@Mao2006])                  |               
| .rf.slr.py           | Conventional SLR and variations (e.g. [@Pauly1991])              |
| .rf.util.py          | RF pulse design utilities                                        |
| .grad.waveform.py    | Gradient and trajectory designers (e.g. [@Kim2003])              |
| .grad.optim.py       | Gradient and trajectory optimization (e.g. [@Lustig2008])        |
| io.py                | Vendor-specific scanner input/output                             |
| linop.py             | Linear operators for pulse design (e.g. [@Grissom2006])          |
| sim.py               | 1-D/N-D/N-coil Bloch simulation (e.g. [@Mansfield1982])          |
| verse.py             | RF pulse/gradient reshaping tools                                |

Preliminary development of this toolbox was presented in reference [@Martin2020a]. The pulse design tools were initially implemented as a sub-package in the SigPy Python package for signal processing and image reconstruction [@Ong2019]. PulPy migrates those tools into a pulse design specific package, with SigPy as an external dependency. PulPy has been streamlined and expanded to include a larger collection of RF and gradient pulse design methods from the literature, as well as additional utility tools for I/O, pulse reshaping, and experimental B$_1^+$-selective pulse design algorithms. The toolbox has proved useful for prototyping novel pulse design algorithms, enabling the publication of Reference [@Martin2022] by the authors and several works from other groups [@Shin2021, @Wu2023]. Figure \ref{pulses_demo} shows an example of RF and gradient waveforms produced by PulPy.


# Statement of need
The field of magnetic resonance imaging is currently experiencing rapid growth in available open source imaging tools. Tools have been made freely available for MRI hardware development [@Amrein2022; @Anand2018], system simulation [@Villena2014; @Stocker2010], pulse sequence programming [@Layton2017], image reconstruction [@Ong2019; @Uecker2015], and post-processing and analysis [@Avants2014; @Duval2018; Soher2023]. The great increase in open-source tools has helped enable fully open-source imaging systems [@Arndt2017, @Artiges2024]. 
However, one critical aspect of the imaging pipeline which has seen limited open-source tool development is RF and gradient pulse design.  While RF pulse and gradient designers increasingly share code online
in independent repositories, there are few sets of common pulse design tools maintained in a rigorous and consistent manner with easy-to-read code and tutorials. This is despite the reality that in many cases, carefully designed or application-specific RF and gradient pulses are crucial to the success of MRI or NMR techniques. An open
source pulse design code library would facilitate the development and dissemination of
novel techniques and the comparison of approaches, similar to how BART [@Uecker2015] and SigPy
[@Ong2019] have made advanced parallel imaging and reconstruction methods widely accessible. To meet this need, we have developed a library of
MRI pulse design tools. We call this new package PulPy, short for Pulses in Python. 


![Example RF and gradient waveforms that PulPy can produce. Top left: 4-channel spokes RF pulse. Bottom left: associated 3-axis spokes gradient waveforms. Top right: PINS excitation RF pulse. Bottom right: associated slice-axis gradient \label{fig:pulses_demo}](pulpy_demo.png){ width=90% }

# Target Audience

The PulPy toolbox has been developed for use by MRI researchers who are interested in pulse sequence design, MRI physics, signal processing, and optimization. We believe that it will serve as an essential building block for more general image acquisition tools which require specialized RF pulses. The toolbox has already been incorporated into open-source sequence development software such as Pulseq [@Layton2017] and PyPulseq [@SravanRavi2019] to provide RF pulses critical to the performance of various pulse sequences. Finally, end-to-end optimization of MRI pulse sequences and reconstructions is 
being increasingly explored [@Radhakrishna2023; @Wang2022]; with the RF pulse and gradient waveform design functions
provided, the PulPy package could facilitate this research. 


Reproducibility and standardization are critical needs in MRI , and any method of reducing methodological variability is desirable. We believe that having centralized references for RF and gradient pulses will help promote consistency between studies by providing common code sources for the most widely used RF and gradient pulses. PulPy's predecessor toolbox, SigPy.RF, also served as a hands-on teaching aid for researchers and students (for example, see the educational ISMRM tutorial associated with [@Martin2020a]). This is a role that the PulPy toolbox will continue to fill. We have developed [several tutorials](https://github.com/jonbmartin/pulpy-tutorials), which are accessible to a wide audience with minimal prior MRI knowledge. 

# Availability and Use

The latest version of PulPy includes the latest stable release of the pulse
design tools and is available from 
[the main repository](https://github.com/jonbmartin/pulpy). It can be installed through pip- see the [documentation](https://pulpy.readthedocs.io/en/latest/) for more details. Jupyter notebook based [pulse design
tutorials](https://github.com/jonbmartin/pulpy-tutorials) for PulPy are also available, which demonstrate the toolbox being used for several classes of pulse design.


# Acknowledgements

We would like to particularly acknowledge that this toolbox is indebted to the MRI scientists who created the RF and gradient pulse design innovations showcased in this software. These are cited as much as possible in this paper and in the PulPy source code. We are particularly thankful for John Pauly's invaluable MATLAB SLR [pulse design toolbox](https://rsl.stanford.edu/research/software.html) [@Pauly1991], which helped inform the core of PulPy's SLR pulse design module. The EPI gradient waveform designer was based on Jeff Fessler's [MIRT](https://web.eecs.umich.edu/~fessler/code/) implementation [@Fessler]. Many other useful case-specific RF pulse design toolboxes not directly incorporated into this toolbox have been created, and we encourage PulPy users to investigate these toolboxes: 

- [Multiband-RF](https://github.com/mriphysics/Multiband-RF) (MATLAB-based) for advanced multiband RF pulse design, incorporating [@AboSeada2019]
- [Spectral-Spatial-RF-Pulse-Design](https://github.com/LarsonLab/Spectral-Spatial-RF-Pulse-Design) (MATLAB-based) for designing spectral-spatial RF pulses for MRS and MR imaging, incorporating [@Larson2008]
- [FastPtx](https://link.springer.com/article/10.1007/s10334-023-01134-7) (Python-based) for designing pTx RF and gradient pulses, from [@Bosch2023]


The authors gratefully acknowledge the assistance provided by Frank Ong, Jonathan Tamir, and Michael Lustig in developing the original SigPy.RF toolbox, which was foundational to this work. We thank the users of SigPy and SigPy.RF for contributing feedback, suggesting new features, and reporting bugs. We acknowledge the support this study received from NIH grants R01 EB 016695 and T32 EB 001628.

# References
