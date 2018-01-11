# miniJATI



### Usage

    miniJATI --input_sequences=../data/MSA_all_combinations_of_gaps_5_leaves.fa --input_tree=../data/tree_5_leaves_r_bl.nwk --lkmove_bothways=true --optim_topology="full-search" 
    --model_substitution="GTR" --model_indels=true

## Verbose options

    export GLOG_v={0,1,2,3}
    export GLOG_minloglevel={INFO,FATAL,WARN}

### Flags

    -input_sequences (File path containing the sequences [FASTA]) type: string
      default: "path/to/fasta.fa"
    -input_tree (File path containing the tree file [NWK]) type: string
      default: "path/to/newick.nwk"
    -lkmove_bothways (Compute likelihood of the model for each topology
      rearrangment (apply and revert)) type: bool default: false
    -model_indels (Extend the substitution model to include Insertion &
      Deletion events (PIP)) type: bool default: false
    -model_substitution (Substitution model to apply on the input data)
      type: string default: "GTR"
    -optim_topology (Topology optimisation under a predefined scheme)
      type: string default: "full-search,smart-search,nni-search,spr-search"
