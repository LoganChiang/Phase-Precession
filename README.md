# Phase-Precession
Model code for phase precession in MATLAB based on Tsodyks et al 1996

Key parameters are contained within lines 3~25, LIF timestep loop begins line 80

Code for Figure 3 and Figure 4a are contained in different blocks of code within the loop, running one OR the other during a given runthrough of code to obtain a figure would save time. Full runthrough with both figures takes ~5 minutes.

Heat map code (for Figure 4a) is written in comments at line 184, use MATLAB Brush function to:

1. Highlight all points

2. Right-click and select 'Export Brushed'

3. Export points as data in the workspace

4. Replace 'XXX' in code with created workspace variable (must be single variable with n rows and 2 columns)

To produce sufficient data for the heatmap multiple runs may need to be saved and combined. Alternatively, the excitation/inhibition parameters or "neuronNumber" parameter (which randomly picks place cells to observe evolution of phase v time) can be tweaked to produce more conducive runs. This model differs from Tsodyks; baseline activity I0 is separated in this model into excitatory and inhibitory counterparts, I0E and I0I respectively on line 17.
