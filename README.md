# Ab_FcR_ODE
NOTE: Currently confirming permission to upload folder 4 from the Kirschner Lab

Arnold Lab, University of Michigan

Melissa Lemke, PhD Candidate

Matlab2019a - parallel computing toolbox

Windows 10

March 2020

Code used in "A personalized systems approach to evaluate 
antibody-Fc receptor activation post RV144 vaccination reveals 
a limited benefit of increasing IgG1 concentration", currently 
in submission 

The code is separated into 4 distinct folders containing all of
the necessary Matlab scripts to run each type of analysis.

	1. Baseline Simulation - Will run one simulation at all
		baseline parameters and plot a timecourse showing
		concentrations of all species and complexes
		note: Baseline IgG1-4 concentrations are the average
		of A244 vaccinee measurements (see Table S1). 
			Plots: Timecourse
			
	2. Personal Simulations - The main folder contents will
		run a baseline simulation for each individual,
		pulling their personal IgG1-4 concentrations from
		the .xlsx file of choice (A244 or BAL). 
			Plots: Validation of vaccinee 1-30 (Fig 2A),
			Bar graph of complex formation for each FcR
			for each person
			
		The subfolder "Personal Local Sensitivity Analysis" 
		will run a local sensitivity (alter each baseline 
		parameter by 0.004X-20X individually - 6 simulations
		per person) which can take ~30 min. 
			Plots: Local sensitivitiy heatmap (Fig 5A)
			Summary of low, med, high complex formation
			within local sensitivity (Fig 5B)
			
	3. Surface Simulation - Will run the 2500 simulations that
		create the surface in Fig 3A & 4A by altering initial
		IgG1 and IgG3 concentrations over 50 different values
		each (0.004X-20X baseline, uniformly spaced
		on a logarithmic scale). Once you move the data from
		the personal baseline simulations into this folder,
		the Make_Figure.m script can calculate the personal
		gradients as well.
			Plots: Move the personal baseline data here!
			then you can plot -> the surface with
			vaccinee dots overlaid (Fig 3A/4A), the IgG1
			& IgG3 gradients (Fig 3B)
			
	4. Global Sensitivity - Modified code provided by Paul Wolberg 
		in Denise Kirschner's Lab (University of Michigan). 
		Based off of:
			Marino, S., Hogue, I.B., Ray, C.J., and 
			Kirschner, D.E. (2008). A methodology for 
			performing global uncertainty and sensitivity 
			analysis in systems biology. J Theor Biol 254, 
			178-196
		https://www.ncbi.nlm.nih.gov/pubmed/18572196

		Our only additions/modifications appear in the 
		"lhs_ode_igg_subtype" folder which is model specific,
		a change of ODE solver in lines 135, 149 & 187 of 
		"lhs_ode_run_new.m" to ode113, and the addition of
		a subfolder "superbar" containing scripts we used
		to plot the PRCC from MathWorks File Exchange. 
		Copyright (c) 2016,  Scott C. Lowe 
		<scott.code.lowe@gmail.com>

		Detailed instructions/information is available in
		the Kirschner lab's own readme.txt file in 
		Global Sensitivity/matlab-ode-lhs. Below we give
		basic instructions on how to run the code as is and
		plot our results.

		Plots: PRCC values and significances (Fig 2B-c)

Instructions:		
Most folders (1-3) contain a 'Run_Me...m' file that you will run 
first to run the simulations and obtain the data. Some folders also 
contain a 'Make_figure' file used to plot the figures once you've 
obtained the data. The rest of the files within the folder are
functions used by these two scripts that you should not need to
directly interact with. The excpetion to this is the "Global 
Sensitivity" folder. 

Basic framework used in folders 1-3:

	Run_me.m -> you will run this file to run simulations
	
	Parameters...m -> will be used by Run_me to obtain the 
		correct intial parameters
		
	Simulate.m -> will be used by Run_me to run each individual
		simulation, sets initial conditions, time, ODE solver
		checks for steady state etc
		
	ODEs.m -> will be used by Simulate.m, contains the ODEs, 
		convservation equations etc.

Folder specific information:	

	1. Baseline Simulation: 
		Simulate.m here will plot the timecourse for you
		Takes ~1 sec
	2. Personal Simulations:
		Contains the personal IgG1-4 conc in .xlsx files
		Parameters...m file here will be pulling from this
		data for each vaccinee simulation
		Will create a data file you will need to copy to
		the "Surface Simulations"folder
		Takes ~20 sec

		Personal Local Sensitivity Analysis:
		USES PARALLEL COMPUTING TOOLBOX
		Also contains "Sensitivity_1D.m" which will be used
		by Run_me....m to run the simulation in parallel, and
		uses "Simulate.m" to do this
		Takes ~20-30 min
	3. Surface Simulations:
		USES PARALLEL COMPUTING TOOLBOX
		Contains "Sensitivity_2D.m" which will be used
		by Run_me....m to run the simulation in parallel, and
		uses "Simulate.m" to do this
		AGAIN, to make the figures, you must move the data
		generated in the main Personal Simulations folder
		into this folder
		Takes ~2 min

	4. Global sensitivity:
		See the readme.txt in 
		Global Sensitivity/matlab-ode-lhs for detailed info
		from the Kirschner lab (source of code)
		Simple instructions on how to run the code as we have
		provided in order to plot the PRCC resuts seen in
		Fig 2B-C:

		1. add the 'matlab-ode-lhs'  folder to your path
		2. Navigate to the 'lhs_ode_igg_subtype' folder
		3. Type in command window
		   >>lhs_ode_run_new(‘lhs_ode_igg_subtype_settings’)
		4. May take a few hours up to about a day
			lowering the sample size (line 165) in 
			"lhs_ode_subtype/lhs_ode_igg_settings.m"
			would decrease this time
		5. When prompted save the output file with an 
		informative name 

		To plot PRCC - code we have added to 
		"lhs_ode_igg_subtype"
			1. Make sure you have superbar.m and 
			supererror.m in your matlab path
			2. Open plot_prcc.m 
			(in “lhs_ode_igg_sutype” folder)
			3. Change: Load(‘your output file’)
			4. Run 

