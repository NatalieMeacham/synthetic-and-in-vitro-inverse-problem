# synthetic-and-in-vitro-inverse-problem

Forward problem: This section of the code defines the model and runs the forward problem for a given sensitivity distribution.

	-File to define different initial sensitivity distributions: DistFn2
 
	-File to solve forward problem and get aggregated tumor volume: ForwardFunctionN
 
	-File to generate forward problem results figures for given parameters: ForwardScript
 
	-File to visualize equilibrium analysis from supplementary material: SenseEqAnalysis

Synthetic data inverse problem: This section of the code generates synthetic tumor growth data with proportional error, then solves the inverse problem, comparing the known original curves with the recovered ones. 

	-File to compute error between data and model: Peerrorfn (used for both synthetic and in vitro datasets)
 
	-File to perform inverse problem: GLSInverseFnN
 
	-File to create synthetic data, run GLSInverseFnN on the data, and generate figures: GLSInverseScriptNPaper
 
	-File to generate figures comparing inverse problem results for synthetic data with different time meshes: timestepplots 
 
	-File to run different distributions, noise levels, and time meshes to generate figures for synthetic data inverse problem results: syntheticfiguremaker 
	
PhenoPop data inverse problem: This section of the code takes in in vitro tumor growth data, processes the data, and runs  the inverse problem on each dosage level of each dataset. 

	-Data from the PhenoPop paper is available at https://github.com/ocbe-uio/PhenoPop/blob/main/Phenopop_inference/data/BaF3-cell_line_data/README.txt
	
 	-File to visualize data: SeeData
	
 	-File to drop NaNs and adjust sizes of datasets: DropNaNsFn
	
 	-File to normalize data and average across dosage replicates: ScaleData
	
 	-Files to recover relevant growth and death rates: ErrorFnFindK, ErrorFnFindRho, ForwardFnFindK, FowardFnFindRho, InferMaxK, InferMaxRho
	
 	-File with core inverse problem function: GLSInverseFnNData
	
 	-File with inverse problem performed on a given dataset: GLSAllDataInclNaNsPaper
	
 	-File to generate figures with in vitro data inverse problem results: DataFigures 
