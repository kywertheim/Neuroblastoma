# Neuroblastoma_CellModel
Disclaimer: I developed the model and its source code under the supervision of Doctor Dawn Walker and in collaboration with Doctor Robert Chisholm. Although I wrote most of the code, the parser was coded by Robert.

About: The files in this repository constitute the first multicellular model of neuroblastoma, built for the PRIMAGE project. It is comprised of a continuous automaton, autonomous agents, and a centre-based mechanical model. Running the Python script named 'Simulator.py' will start a stochastic simulation of neuroblastoma progression. Note that the stochastic nature of the problem means multiple runs are needed to obtain robust averages.

Files:
1. The Python script named 'Simulator.py' should be implemented in Python 3.8.10.
2. The directory named 'PyNB' must be in the same directory as the Python script during a simulation.
3. The three files named 'Environment.py', 'Neuroblastoma.py', and 'Schwann.py' must be in the directory named 'PyNB' during a simulation.

Python Modules:
1. itertools.
2. numpy.
3. math.
4. statistics.
5. random.
6. copy.
7. matplotlib.
8. mpl_toolkits.mplot3d.
9. time.
10. argparse.

Configuration of the file named 'Environment.py' to parameterise and initialise the model to describe the problem under consideration:
1. Specify the initial volume of the tumour or tumoroid: V_tumour on line 56.
2. Specify the cellularity of the volume defined above: cellularity on line 70.
3. Specify the histology: histology_init on line 71. 0 is neuroblastoma, 1 is ganglioneuroblastoma, 2 is nodular ganglioneuroblastoma, 3 is intermixed ganglioneuroblastoma, 4 is ganglioneuroma, 5 is maturing ganglioneuroma, and 6 is mature ganglioneuroma.
4. Specify the grade of differentiation of the neuroblastic component of the volume: gradiff on line 84. 0 is undifferentiated, 1 is pooly differentiated, and 2 is differentiating.
5. Specify the initial level of oxygen, which is also a proxy for all types of nutrients: O2 on line 119.
6. Define the chemotherapy regimen: chemo_start on line 120, chemo_end on line 121, and chemo_effects on line 122.
7. Define the mutations present in the neuroblastic cells: MYCN_amp on line 133 (MYCN amplification), TERT_rarngm on line 134 (TERT rearrangement), ATRX_inact on line 135 (ATRX inactivation), ALT on line 136 (alternative lengthening of telomeres), and ALK on line 137 (0 means wild type, 1 means ALK amplification or activation, and 2 means other RAS mutations).

Configuration of the file named 'Simulator.py' to define the extent of the simulation:
1. Specify the number of one-hour steps making up the simulation: steps on line 1542.
