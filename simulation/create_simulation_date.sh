# Create simulation_files

Rscript simulate_sites.R \
  --chr_len 10000000 \
  --site_lambda 20000 \
  --block_frac 0.05 \
  --delta_grid 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0 \
  --block_len_lo 300 --block_len_hi 2000 \
  --groups WT,MT \
  --rep_per_group 4 \
  --rho 0.3 \
  --logcov_mu 3.0 \
  --logcov_sd 1.2 \
  --dmr_site_prop 0.5 \
  --mu_site_sd 0.05 \
  --miss_rate 0.30 \
  --min_cov 5 \
  --seed 202509
