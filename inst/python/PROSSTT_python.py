def PROSSTT_sim_Python(newick_string,
                       gene_num,
                       seed,
                       modules = 15,
                       alpha = 0.2,
                       beta = 3):
 
  import numpy as np
  from numpy import random
  import newick
  from prosstt import tree
  from prosstt import simulation as sim
  from prosstt import sim_utils as sut
  # set random seed
  random.seed(seed)
  newick_string = newick_string
  newick_tree = newick.loads(newick_string)
  lineage = tree.Tree.from_newick(newick_string, modules = int(modules), genes = gene_num, density = None)
  uMs, Ws, Hs = sim.simulate_lineage(lineage, a = 0.05)
  gene_scale = sut.simulate_base_gene_exp(lineage, uMs)
  lineage.add_genes(uMs, gene_scale)
  random.seed(seed)
  X1, labs1, brns1, scalings1 = sim.sample_whole_tree(lineage, 1, alpha = alpha, beta = beta)
  return X1
