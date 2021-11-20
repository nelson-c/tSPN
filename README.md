# tSPN

A single compartment conductance-based thoracic sympathetic postganglionic neuron (tSPN) neuron model. 

## syapticGain.m
This is the main script used to generate the figures used in the document. Additional documentation is provided within the script. Each section (denoted with %%) generates a different figure. There are several subfunctions within the script which allow 
Calling the defaultSetup function returns a structure with fields that contain default values for a simulation. Different simulations can be run by first modifying the setup structure.  
Calling the tSPN_net function runs a complete simulation which returns all presynaptic event times, postsynaptic event times, synaptic conductance traces, and voltage traces. 
Calling the tSPN_sparseNet function runs a reduced simulation that is much faster, especially for low-firing rate simulations, but does not return as many variables. The sparse function only returns the number of presynaptic events and postsynaptic events and is only useful for calculating synaptic gain. It saves time by not running the simulation unless there are multiple synaptic events within a temporal window. The basic time-saving algorithm is to first calculate the threshold synaptic amplitude that would lead to an action potential. Then it splits the spike train into groups with at least 1 second between events. If there is a single event in a group (no other synaptic event occurred within a second before or within a second after the event) then simply compare the amplitude of the synaptic amplitude to the threshold. If it’s suprathreshold, then increase # spikes by 1. If it’s subthreshold, don’t increase the number of postsynaptic events. If there are 2 or more events within the group, run the simulation for the entire group and count spikes. If there are four events, each separated by 500ms, they would all be included in the same group because even though the total duration (2s) is longer than the temporal window (1s), there is less than 1s between each successive event. 

## tSPN.m
This is the function that runs the numerical simulation. Usage and documentation provided within the script.  In short, you provide: 
1)	a vector of conductance values
2)	a vector of injected current values
3)	optionally, you may also provide a starting value, y0 (you can initialize the simulation once, and then provide this value to prevent the need to initialize every time you run the simulation)
4)	optionally, a synaptic conductance trace (obviously necessary for doing any synaptic modeling). 

The function returns:
1)	A voltage trace
2)	A parameter matrix with many different dynamic variables. Documentation provided within code.

## tSPN_gsyn.m
Takes in a vector of event times and event weights and returns a synaptic conductance trace, gsyn, which can be used with the tSPN() function. 
