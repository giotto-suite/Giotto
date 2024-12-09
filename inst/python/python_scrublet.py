import scrublet as scr

def python_scrublet(counts_matrix, expected_doublet_rate, min_counts, min_cells, min_gene_variability_pctl, n_prin_comps, seed_number=1234):
  
  min_counts = int(min_counts)
  min_cells = int(min_cells)
  min_gene_variability_pctl = int(min_gene_variability_pctl)
  n_prin_comps = int(n_prin_comps)
  random_state = int(seed_number)
  
  scrub = scr.Scrublet(counts_matrix=counts_matrix, expected_doublet_rate=expected_doublet_rate, random_state=random_state)
  doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=min_counts,min_cells=min_cells, min_gene_variability_pctl=min_gene_variability_pctl, n_prin_comps=n_prin_comps)
  
  return_list = []
  return_list.append(doublet_scores)
  return_list.append(predicted_doublets)
  
  return(return_list)
