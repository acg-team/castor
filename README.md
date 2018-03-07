# miniJATI



### Usage
    Joint Alignment Tree Inference (miniJATI) 0.1.1 (refs/heads/develop afe0ca955db81c72c1ae705c41f9c7b71080e021, 07 Mar 2018, 18:07:17)
    Usage: miniJATI [options]


### Alphabet options

    alphabet={DNA|RNA|Protein)|Codon(letter={DNA|RNA},type={Standard|EchinodermMitochondrial|InvertebrateMitochondrial|VertebrateMitochondrial})}
                                            The alphabet to use when reading sequences. DNA and RNA alphabet can in addition take
                                            an argument:
                                            bangAsgap={bool}
                                            Tell is exclamation mark should be considered as a gap character. The default
                                            is to consider it as an unknown character such as 'N' or '?'.
    genetic_code={translation table}        Where ’translation table’ specifies the code to use, either as a text description,
                                            or as the NCBI number. The following table give the currently implemented codes
                                            with their corresponding names:
                                            Standard                    1
                                            VertebrateMitochondrial     2
                                            YeastMitochondrial          3
                                            MoldMitochondrial           4
                                            InvertebrateMitochondrial   5
                                            EchinodermMitochondrial     9
                                            AscidianMitochondrial       13
                                            The states of the alphabets are in alphabetical order.

### [Input] Reading sequences

    input.sequence.file={path}                                The sequence file to use. (These sequences can also be not aligned).
    input.sequence.format={format}                            The sequence file format.
    input.sequence.sites_to_use={all|nogap|complete}          Tells which sites to use
    input.sequence.remove_stop_codons={boolean}               Removes the sites where there is a stop codon (default: ’yes’)
    input.sequence.max_gap_allowed=100%                       It specifies the maximum amount of gap allowed per site.
    input.sequence.max_unresolved_allowed=100%                It specifies the maximum amount of unresolved states per site.
    input.site.selection={list of integers}                   Will only consider sites in the given list of positions, in extended format :
                                                              positions separated with ",", and "i-j" for all positions between i and j,
                                                              included.
    input.site.selection = {Sample(n={integer} [, replace={true}])}
                                                              Will consider {n} random sites, with optional replacement.

The following formats are currently supported:

    Fasta(extended={bool}, strictNames={bool})                The fasta format. The argument extended, default to 'no'
                                                              allows to enable the HUPO-PSI extension of the format.
                                                              The argument strict_names, default to 'no', specifies
                                                              that only the first word in the fasta header is used as
                                                              a sequence names, the rest of the header being considered
                                                              as comments.
    Mase(siteSelection={chars})                               The Mase format (as read by Seaview and Phylo_win for instance),
                                                              with an optional site selection name.
    Phylip(order={interleaved|sequential}, type={classic|extended}, split={spaces|tab})
                                                              The Phylip format, with several variations.
                                                              The argument order distinguishes between sequential and interleaved
                                                              format, while the option type distinguished between the plain old
                                                              Phylip format and the more recent extention allowing for sequence
                                                              names longer than 10 characters, as understood by PAML and PhyML.
                                                              Finally, the split argument specifies the type of character that
                                                              separates the sequence name from the sequence content.
                                                              The conventional option is to use one (classic) or more (extended)
                                                              spaces, but tabs can also be used instead.
    Clustal(extraSpaces={int})                                The Clustal format.
                                                              In its basic set up, sequence names do not have space characters,
                                                              and one space splits the sequence content from its name. The parser
                                                              can however be configured to allow for spaces in the sequence names,
                                                              providing a minimum number of space characters is used to split
                                                              the content from the name. Setting extraSpaces to 5 for instance, the
                                                              sequences are expected to be at least 6 spaces away for their names.
    Dcse()                                                    The DCSE alignment format. The secondary structure annotation will be ignored.
    Nexus()                                                   The Nexus alignment format. (Only very basic support is provided)
    GenBank()                                                 The GenBank not aligned sequences format.
                                                              Very basic support: only retrieves the sequence content for now,
                                                              all features are ignored.

### [Input] Reading trees

    input.tree.file={path}                      The phylogenetic tree file to use.
    input.tree.format={Newick|Nexus|NHX}        The format of the input tree file.

###  Alignment options

    alignment=<bool>                                     [requested]

###  Branch lengths initial values

    init.tree={user|random|distance}                  Set the method for the initial tree to use.
                                                      The user option allows you to use an existing file passed via input.tree.file
                                                      This file may have been built using another method like neighbor joining or
                                                      parsimony for instance. The random option picks a random tree, which is handy
                                                      to test convergence.  This may however slows down significantly the optimization  
                                                      process.
    init.distance.matrix.file={path}                  A distance matrix can be supplied instead of being computed from the alignment.
    init.distance.method={wpgma|upgma|nj|bionj}       When distance method is required, the user can specify which algorithm to use.
    init.brlen.method={method description}            Set how to initialize the branch lengths. Available methods include:

    Input(midpoint_root_branch={boolean})             Keep initial branch lengths as is. Additional argument specifies if the root
                                                      position should be moved to the midpoint position of the branch containing it.

    Equal(value={float>0})                            Set all branch lengths to the same value, provided as argumemt.

    Clock                                             Coerce to a clock tree.

    Grafen(height={{real>0}|input}, rho = {real>0})   Uses Grafen’s method to compute branch lengths.
                                                      In Grafen’s method, each node is given a weight equal to the number of underlying  
                                                      leaves. The length of each branch is then computed as the difference of the weights
                                                      of the connected nodes, and further divided by the number of leaves in the tree.
                                                      The height of all nodes are then raised to the power of ’rho’, a user specified value.
                                                      The tree is finally scaled to match a given total height, which can be the original
                                                      one (height=input), or fixed to a certain value (usually height=1). A value of
                                                      rho=0 provides a star tree, and the greater the value of rho, the more recent the
                                                      inner nodes.

    input.tree.check_root = {boolean}                 Tell if the input tree should be checked regarding to the presence of a root.

### Substitution model options

    model=<string>                                      A description of the substitution model to use, using the keyval syntax.

The following nucleotide models are currently available:

    JC69
    K80([kappa={real>0}])
    F84([kappa={real>0}, theta={real]0,1[}, theta1={real]0,1[},theta2={real]0,1[} ,"equilibrium frequencies"])
    HKY85([kappa={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,"equilibrium frequencies"])
    T92([kappa={real>0}, theta={real]0,1[} ,"equilibrium frequencies"])
    TN93([kappa1={real>0}, kappa2={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,"equilibrium frequencies"])
    GTR([a={real>0}, b={real>0}, c={real>0}, d={real>0}, e={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,"equilibrium frequencies"])
    L95([beta={real>0}, gamma={real>0}, delta={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,"equilibrium frequencies"])
    SSR([beta={real>0}, gamma={real>0}, delta={real>0}, theta={real]0,1[}])
    RN95([thetaR={real]0,1[}, thetaC={real]0,1[}, thetaG={real]0,1[}, kappaP={real[0,1[}, gammaP={real[0,1[}, sigmaP={real>1}, alphaP={real>1}])
    RN95s([thetaA={real]0,0.5[}, gamma={real]0,0.5[}, alphaP={real>1}])

The following protein models are currently available:

    JC69
    DSO78
    JTT92
    WAG01
    LG08
    LLG08_EX2([relrate1={real]0,1[}, relproba1={real]0,1[}])
    LLG08_EX3([relrate1={real]0,1[}, relrate2={real]0,1[}, relproba1={real]0,1[}, relproba2={real]0,1[}])
    LLG08_EHO([relrate1={real]0,1[}, relrate2={real]0,1[}, relproba1={real]0,1[}, relproba2={real]0,1[}])
    LLG08_UL2([relrate1={real]0,1[}, relproba1={real]0,1[}])
    LLG08_UL3([relrate1={real]0,1[}, relrate2={real]0,1[}, relproba1={real]0,1[}, relproba2={real]0,1[}])
    LGL08_CAT(nbCat={[10,20,30,40,50,60]}, [relrate1={real]0,1[}, relrate2={real]0,1[}, ..., relproba1={real]0,1[}, relproba2={real]0,1[}, ...] ))
    LGL08_CAT_C{[1,...,nbCat]}(nbCat={[10,20,30,40,50,60]})
    DSO78+F([theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ... ,"equilibrium frequencies"])
    JTT92+F([theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ..., "equilibrium frequencies"])
    WAG01+F([theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ..., "equilibrium frequencies"])
    LG08+F([theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ..., "equilibrium frequencies"])
    Empirical(name={chars}, file={path})
    Empirical+F(name={chars}, file={path}, [theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[}, ...,  "equilibrium frequencies"])

The following meta models are currently available:

    PIP(model={model description}, lambda={real>0}, mu={real>0} [, "equilibrium frequencies"])
    TS98(model={model description}, s1={real>0}, s2={real>0} [, "equilibrium frequencies"])
    G01(model={model description}, rdist={rate distribution description}, mu={real>0} [, "equilibrium frequencies"])
    RE08(model={model description}, lambda={real>0}, mu={real>0} [, "equilibrium frequencies"])

### Frequencies distribution sets

    The following frequencies distributions are available:
    Fixed()                                                                  All frequencies are fixed to their initial value and are not estimated.

    GC(theta={real]0,1[})                                                    For nucleotides only, set the G content equal to the C content.

    Full(theta1={real]0,1[}, theta2={real]0,1[}, ..., thetaN={real]0,1[})    Full parametrization. Contains N free parameters.


### Rate across site distribution

    rate_distribution={rate distribution description}                Specify the rate across sites distribution

The following distributions are currently available:

[Most used distributions]

    Constant                                                         Uses a constant rate across sites

    Gamma(n={int>=2}, alpha={float>0})                               A discretized gamma distribution of rates, with n classes, and a given shape,
                                                                     with mean 1 (scale=shape).

    Invariant(dist={rate distribution description}, p={real[0,1]})   A composite distribution allowing a special class of invariant site, with a
                                                                     probability p.

[Standard distributions]

    Beta(n={int>=2}, alpha={float>0}, beta={float>0})                A discretized beta distribution, with n classes, with standard parameters
                                                                     alpha and beta.

    Gamma(n={int>=2}, alpha={float>0}, beta={float>0})               A discretized gamma distribution, with n classes, a shape alpha and a rate
                                                                     beta, with parameters alpha and beta.

    Gaussian(n={int>=1}, mu={float}, sigma={float>0})                a discretized gaussian distribution, with n classes, a mean mu and a
                                                                     standard deviation sigma, with parameters mu and sigma.

    Exponential(n={int>=2}, lambda={float>0})                        a discretized exponential distribution, with n classes and parameter lambda.

    Simple(values={vector<double>}, probas={vector<double>} [, ranges={vector<parametername[min;max]>}])
                                                                     A discrete distribution with specific values (in values) and their respective
                                                                     non-negative probabibilities (in probas). The parameters are V1, V2, ..., Vn
                                                                     for all the values and the relative probabibility parameters are
                                                                     theta1, theta2, ..., thetan-1.
                                                                     Optional argument {ranges} sets the allowed ranges of values taken by the
                                                                     parameters; usage is like 'ranges=(V1[0.2;0.9],V2[1.1;999])'.

    TruncExponential(n={int>=2}, lambda={float>0}, tp={float>0})     A discretized truncated exponential distribution, with n classes, parameter
                                                                     lambda and a truncation point tp. The parameters are lambda and tp.

    Uniform(n={int>=1}, begin={float>0}, end={float>0})              A uniform distribution, with n classes in interval [begin,end].
                                                                     There are no parameters.

[Mixture Distributions]

    Invariant(dist={distribution description}, p={float>0})          A Mixture of a given discrete distributution and a 0 Dirac. p is the  
                                                                     probability of this 0 Dirac.

    Mixture(probas={vector<double>}, dist1={distribution description}, ..., distn={distribution description})
                                                                     A Mixture of discrete distributions with specific probabilities (in probas)
                                                                     and their respective desccriptions (in probas). The parameters are the relative
                                                                     probabibility parameters theta1, theta2, ..., thetan-1, and the parameters
                                                                     of the included distributions prefixed by Mixture.i_ where i is the order
                                                                     of the distribution.


### Likelihood computation options


### Numerical parameter optimisation options

This program allows to (re-)estimate numerical parameters, including
- Branch lengths
- Entries of the substitution matrices, included base frequencies values)
- Parameters of the rate distribution (currently shape parameter of the gamma law, proportion of invariant sites).


    optimization={method}

The following methods are currently available:

    None                                                         No optimization is performed, initial values are kept 'as is'.

    FullD(derivatives={Newton|Gradient})                         Full-derivatives method.
                                                                 Branch length derivatives are computed analytically, others numerically.
                                                                 The derivatives arguments specifies if first or second order derivatives
                                                                 should be used. In the first case, the optimization method used is the  
                                                                 so-called conjugate gradient method, otherwise the Newton-Raphson method
                                                                 will be used.

    D-Brent(derivatives={Newton|Gradient|BFGS}, nstep={int>0})   Branch lengths parameters are optimized using either the conjugate gradient,
                                                                 the Newton-Raphson method or BFGS, other parameters are estimated using the
                                                                 Brent method in one dimension. The algorithm then loops over all parameters
                                                                 until convergence. The nstep arguments allow to specify a number of
                                                                 progressive steps to perform during optimization. If nstep=3 and
                                                                 precision=E-6, a first optimization with precision=E-2, will be performed,
                                                                 then a round with precision set to E-4 and finally precision will be set to
                                                                 E-6. This approach generally increases convergence time.

    D-BFGS(derivatives={Newton|Gradient|BFGS}, nstep={int>0})    Branch lengths parameters are optimized using either the conjugate gradient,  
                                                                 the Newton-Raphson method or BFGS, other parameters are estimated using the
                                                                 BFGS method in one dimension. The algorithm then loops over all parameters
                                                                 until convergence. The nstep arguments allow to specify a number of
                                                                 progressive steps to perform during optimization. If nstep=3 and
                                                                 precision=E-6, a first optimization with precision=E-2, will be performed,
                                                                 then a round with precision set to E-4 and finally precision will be set to
                                                                 E-6. This approach generally increases convergence time.

    optimization.reparametrization=<bool>                        Tells if parameters should be transformed in order to remove constraints (for
                                                                 instance positivie-only parameters will be log transformed in order to obtain
                                                                 parameters defined from -inf to +inf). This may improve the optimization,
                                                                 particularly for parameter-rich models, but the likelihood calculations
                                                                 will take a bit more time.
    optimization.final={powell|simplex|bfgs}                     Optional final optimization step, useful if numerical derivatives are
                                                                 to be used. Leave the field empty in order to skip this step.
    optimization.profiler={{path}|std|none}                      A file where to dump optimization steps (a file path or std for standard
                                                                 output or none for no output).
    optimization.message_handler={{path}|std|none}               A file where to dump warning messages.
    optimization.max_number_f_eval=<int>0>                       The maximum number of likelihood evaluations to perform.
    optimization.tolerance=<float>0>                             The precision on the log-likelihood to reach.
    optimization.constrain_parameter={list<chars=interval>}      A list of parameters on which the authorized values are limited to a given
                                                                 interval. For example:
                                                                 optimization.constrain_parameter = YN98.omega = [-inf;1.9[, *theta* = [0.1;0.7[, BrLen*=[0.01;inf]
    optimization.ignore_parameter={list<chars>}                  A list of parameters to ignore during the estimation process. The parameter
                                                                 name should include there 'namespace', that is their model name, for
                                                                 instance K80.kappa,TN93.theta, GTR.a, Gamma.alpha, etc.'BrLen' will
                                                                 ignore all branch lengths and 'Model' will ignore all model parameters.
                                                                 The ’*’ wildcard can be used, as in *theta* for all the parameters whose
                                                                 name has theta in it.

### Topology optimisation options

    optimization.topology=<bool>                                                                Enable the tree topology estimation
    optimization.topology.algorithm={greedy|hillclimbing|swarm}                                 Algorithm to use for topology estimation
    optimization.topology.algorithm.operations={nni-search|spr-search|tbr-search|best-search}   Tree rearrangment operations to apply during tree-search
    optimization.topology.algorithm.maxcycles=<int>                                             Max number of tree-search cycles
    optimization.topology.likelihood={single,double}                                            Method for recomputing the likelihood after tree move
    optimization.topology.algorithm.hillclimbing.startnodes=<int>                               Number of starting random nodes to use during hc

### Output options

    output.tree.file={path}                         The phylogenetic tree file to write to.
    output.tree.format={Newick|Nexus|NHX}           The format of the output tree file.
    output.trees.file={path}                        The file that will contain multiple trees.
    output.trees.format={Newick|Nexus|NHX}          The format of the output tree file.
    output.infos={{path}|none}                      Alignment information log file (site specific rates, etc)
    output.estimates={{path}|none}                  Write parameter estimated values.
    output.estimates.alias={boolean}                Write the alias names of the aliased parameters instead of their values (default: true).

###  Verbosity

    export GLOG_v={0,1,2,3}
    export GLOG_minloglevel={INFO,FATAL,WARN}



### Examples

#### Analysis 1

- Dataset: sim_dna_5.fa (5 simulated nucleotide sequences with gaps) | 100nt (use all sites)
- Initial tree: distance based (bioNJ)
- Substitution model: PIP(lambda=0.1,mu=0.2)+GTR()+Gamma(cat=4)
- Optimisation= D-BFGS(derivatives=BFGS) + Topology(hillclimbing(n=4))


    ./miniJATI alphabet=DNA alignment=false input.sequence.file=../data/test_treesearch/sim_dna_5.fa input.sequence.sites_to_use=all model=PIP(model=GTR(),lambda=0.1,mu=0.2) rate_distribution=Gamma(n=4) init.tree=distance optimization=D-BFGS(derivatives=BFGS) optimization.message_handler=../out_opt.txt optimization.profiler=../out_opt_prof.txt optimization.reparametrization=false optimization.max_number_f_eval=1000 optimization.tolerance=0.001 optimization.final=bfgs optimization.topology=true optimization.topology.algorithm=hillclimbing optimization.topology.algorithm.hillclimbing.startnodes=4 optimization.topology.algorithm.operations=best-search optimization.topology.algorithm.maxcycles=100 optimization.topology.likelihood=bothways output.tree.file=../5leaves.nwk output.infos=../5leaves_info.log output.estimates=../5leaves_estimates.log

##### Results

- **Tree** ((jkqnk:10.7483,(zgedq:5.41497,gojyc:6.69881):0.996304):0.972286,(zoixi:3.22681,mcwps:5.07359):1e-06);<br>
- **Log likelihood :** -752.29059764981366242864<br>
- **Number of sites:** 100
- **Substitution model parameters: ** model=GTR+PIP(GTR.a=13198.746763264866,GTR.b=39108.651431432860,GTR.c=8549.078691481493,GTR.d=16777.297593105031,GTR.e=19546.980416169525,GTR.theta=0.537138911249,GTR.theta1=0.551576311231,GTR.theta2=0.544394277806,lambda=1.382997811983,mu=0.020728094422)
- **Rate distribution parameters: ** rate_distribution=Gamma(n=4,alpha=7722.244281054756, Gamma.beta=alpha)
- **Duration:** 27.000000s

### Wikipages

Documentation can be found at https://bitbucket.org/acg-team/minijati/
