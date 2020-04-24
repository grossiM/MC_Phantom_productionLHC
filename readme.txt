Readme about all command to generate phantom event:
file,command used:
1) Firt step - generate GRIDS:
 python submit_phantom_final.py ww_semileptonic_mu_unpol_nob.cfg --verbose
NOTE: ok works

2) second step: CHECK grid generated and calculate xsection:
   python gridPack_verification.py ww_semileptonic_mu_trans.cfg

3) Third step: once generated  gridpack, generate events:
  python submit_phantom_final.py ww_semileptonic_mu_unpol_nob.cfg --produce 20000 100 > submission_unpol_list.log

