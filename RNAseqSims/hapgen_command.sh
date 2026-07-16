# After using the following command to get the recombination map for GBR 
# wget ftp://ftp.1000genomes.ebi.ac.uk:21/vol1/ftp/technical/working/ 20130507_omni_recombination_rates/GBR_omni_recombination_20130507.tar

# Get first SNP position
FIRST_POS=$(awk 'NR==2 {print $2}' gbr_subset.legend)

./hapgen2 \
  -m GBR/GBR-22-final.txt \
  -l gbr_subset.legend \
  -h gbr_subset.hap \
  -o gbr_simulated \
  -n 60 0 \
  -dl $FIRST_POS 0 0 0 \
  -no_haps_output

