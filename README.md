## Gene Regulatory Network Generator in Julia (JuGRN)

### Introduction ###
JuGRN is a code generation system that transforms simple descriptions of the connectivity of gene regulatory networks into model code written in the [Julia](http://julialang.org) programming language. JuGRN has been used in the publications:

1. [Tasseff R, Jensen H, Congleton J, Dai W, Rogers K, Sagar A, Yen A and J. Varner (2017) An Effective Model of the Retinoic Acid Induced Differentiation Program, Sci Reports, 7:14327 doi:10.1038/s41598-017-14523-5](https://www.nature.com/articles/s41598-017-14523-5)
2. [Gould R, Bassen DM, Chakrabarti A, Varner JD and Butcher J (2016) Population Heterogeneity in the Epithelial to Mesenchymal Transition Is Controlled by NFAT and Phosphorylated Sp1. PLoS Comput Biol 12(12): e1005251. doi:10.1371/journal.pcbi.1005251](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005251)

to generate gene regulatory network code.

### Installation and Requirements
``JuGRN.jl`` is organized as a [Julia](http://julialang.org) package which 
can be installed in the ``package mode`` of Julia.

Start of the [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/index.html) and enter the ``package mode`` using the ``]`` key (to get back press the ``backspace`` or ``^C`` keys). Then, at the prompt enter:

    (v1.1) pkg> add https://github.com/varnerlab/JuGRN-Generator.git

This will install the ``JuGRN.jl`` package and the other required packages.
``JuGRN.jl`` requires Julia 1.3.x and above.

``JuGRN.jl`` is open source. 
You can download this repository as a zip file, or clone or pull it by using the command (from the command-line):

	$ git pull https://github.com/varnerlab/JuGRN-Generator.git

or

	$ git clone https://github.com/varnerlab/JuGRN-Generator.git

### How do I generate model code? ###
To generate a GRN model, issue the command ``make_julia_model.jl`` from the command line (outside of the REPL):

	$ julia make_julia_model.jl -m <input path> -o <output path> -s <host type>

The ``make_julia_model.jl`` command takes four command line arguments:

Argument | Required | Default | Description
--- | --- | --- | ---
-m | Yes	| none | Path to model input file (your \*.net file)
-o | No	| current directory | Path where files are written
-s | No	| bacterial | Host type (bacterial \| mammalian)

### Format for the GRN model input file ###
JuGRN-Generator transforms structured flat files into GRN model code. JuGRN-Generator takes flat files of the form:

~~~
// ----------------------------------------------------------------------- //
// JuGRN interactions -
//
// Record:
// actor {activate(*),induce(*) | inhibit(*),repress(*)} (target,...)
// ---------------------------------------------------------------------- //

// three gene memory network -
gene_1 induces (gene_2,gene_3)
gene_2 activates gene_3
gene_3 activates gene_2

~~~

The model specification file (by default given the filename `Network.net`) defines the biology of the model that gets generated.
JuGRN generates the files:

Filename | Description
--- | ---
``AdjDriver.jl`` | Driver function to solve the adjoint system of equations (sensitivity analysis)
``Balances.jl`` | Material balance equations for genes, mRNA and proteins in the GRN
``Control.jl`` | Encodes the control logic described in your GRN network file
``DataDictionary.jl`` | Encodes the model parameters e.g., initial conditions or promoter function parameters in a [Julia dictionary](https://docs.julialang.org/en/stable/stdlib/collections/#Base.Dict)
``Degradation.dat`` | Stoichiometric matrix for mRNA and protein degradation reactions
``Dilution.dat`` | Stoichiometric matrix for mRNA and protein dilution reactions
``Discrete.jl`` | Contains methods to evaluate the discrete mass balance equations (discretized using a zero-order-hold)
``Driver.jl`` | Example script that can be used to solve the continuous material balance equations
``Include.jl`` | Includes all the JuGRN files into the current workspace
``Kinetics.jl`` | Encodes the rate of transcription, translation and degradation for mRNA and protein species
``Network.dat`` | Stoichiometric array for the transcription and translation reactions
``SolveBalances.jl`` | Solves the material balance equations using ODE solvers from the [ODE package](https://github.com/JuliaDiffEq/ODE.jl)
``Utility.jl`` | Encodes utility functions required for sensitivity analysis (e.g., computation of the Jacobian)
