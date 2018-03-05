# miniJATI



### Usage

    ./miniJATI alphabet=DNA alignment=false input_sequences=../data/test_treesearch/sim_dna_5.fa lkmove=bothways optimization=D-Brent(derivatives=Newton) optimization
    .message_handler=../out_opt.txt optimization.profiler=../out_opt_prof.txt optimization.reparametrization=false optimization.max_number_f_eval=1000 optimization.tolerance=0.001 optim_topology_algorithm=greedy optim_topology_numnodes=4 optim_topology_operations=best-search optim_topology_maxcycles=100 model_substitution=HKY85 model_indels=true model_setfreqsfromdata=false model_pip_lambda=0.2 model_pip_mu=0.1




## Verbose options

    export GLOG_v={0,1,2,3}
    export GLOG_minloglevel={INFO,FATAL,WARN}

### Arguments

    alphabet={DNA,RNA,Protein}                      (inherited from bpp documentation)      [requested]
    alignment={true,false}                                                                  [requested]  
    input_sequences=path/to/sequences.fa                                                    [requested]     
    input_tree=path/to/tree.nwk                                                             [not-requested] 
    lkmove="bothways"                                                                       [not-requested]
    optimization=D-Brent(derivatives=Newton)        (inherited from bpp documentation)      [requested]
    optimization.message_handler=../out_opt.txt     (inherited from bpp documentation)      [not-requested]
    optimization.profiler=../out_opt_prof.txt       (inherited from bpp documentation)      [not-requested]
    optimization.reparametrization=false            (inherited from bpp documentation)      [not-requested]
    optimization.max_number_f_eval=1000             (inherited from bpp documentation)      [not-requested]
    optimization.tolerance=0.001                    (inherited from bpp documentation)      [not-requested]
    optim_topology_algorithm="greedy"               {none, greedy, hillclimbing, swarm}     [requested]
    optim_topology_numnodes=4                       number of random nodes to use during the hillclimbing 
    optim_topology_operations="best-search"         {nni-search, spr-search, tbr-search, best-search}    
    optim_topology_maxcycles=100
    model_substitution="HKY85"                      (inherited from bpp documentation)  [requested]
    model_indels=true                               {true, false}
    model_setfreqsfromdata=false                    (inherited from bpp documentation)  [requested]
    model_pip_lambda=0.2                            (if models_indels=true)
    model_pip_mu=0.1                                (if models_indels=true)