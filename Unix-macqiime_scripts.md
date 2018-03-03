# The following Unix, perl, and macqiime scripts were used to analyze the Ion Torrent data.
## Quality control & chimera removal
This was done using macqiime running in the Terminal.

`$split_libraries.py -m mapping.txt -f seqs.fasta -q quality.qual -r -l 200 -L350 -M 1 -w 50 -g -b variable_length -o split_libraries -a 0 -H 6`

Chimeras were detected using Usearch61. 
`$ ~/usearch6.1.544_i86osx32  ~/ITXseqs.fasta -uchime_ref -db ~/fungiLSU_train_041714.fasta  -strand plus -uchimeout results_ushime150721.txt -chimeras chimeras.fasta -nonchimeras nonchimeras.fasta`

Chimeras were removed using macqiime.

`$ filter_fasta.py -f ~/ITXseqs.fasta -a ~/chimeras.fasta -n -o chimera_filtered_ITXseq.fasta`

## OTU picking & taxonomy

`$ pick_otus.py -i ~/chimera_filtered_ITXseq.fasta -m uclust -s 0.95 -o picked_otus_95`

Made a representative set of sequences.

`$pick_rep_set.py -i ~/chimera_filtered_ITXseq_otus.txt -m most_abundant -o rep_set.fasta`

### Taxonomy with uclust

`$ assign_taxonomy.py -i ~/rep_set.fasta -r ~/LSU_Silva_accessions.fasta -t ~/LSU_Silva_rdp.txt -m uclust -o uclust_SILVA_assigned_taxonomy/`

### Taxonomy with BLAST + MEGAN

`$ blastn -db nt -query ~/rep_set.fasta -out results/`

### Taxonomy with RDP-NBC

`$java -Xmx1g -jar ~/rdp_classifier_2.7/dist/classifier.jar ~/mytrained/rRNAClassifier.properties -o classify.txt ~/rep_set.fasta`

The output file rep_set_tax_assignments.txt generated with assign_taxonomy.py was edited to reflect the consensus of the three methods.

## OTU table generation & filtering

`$ make_otu_table.py -i ~/chimera_filtered_ITXseq_otus.txt -t ~/rep_set_tax_assignments.txt -o all_otus_tax.biom`

Potential contamination was removed using the positive control sample. First, an OTU table with only the OTUs appearing in the positive control was created.

`$ filter_samples_from_otu_table.py -i ~/all_otus_tax.biom -m mapping.txt -s ‘Treatment:Control’ -o otu_table_control_only.biom`

The control OTU table was filtered to remove an OTU with an abundance <1.

`$ filter_otus_from_otu_table.py -i otu_table_control_only.biom -n 1 -o filtered_otu_table_control.biom`

The control OTU table was converted to a list of OTUs.

`$ biom convert -b -i filtered_otu_table_control.biom -o otus_to_remove.txt`

The list was then used to remove potential contaminants from the full OTU set.

`$ filter_otus_from_otu_table.py -i ~/ all_otus_tax.biom -e otus_to_remove.txt -o filtered_all_otus.biom` 

The dataset was reduced to just the temporary forest pond October and December samples.

`$ filter_samples_from_otu_table.py -i ~/filtered_all_otus.biom -m mapping.txt -s ‘DOB:oct,dec’ -o filtered_otus_oct_dec.biom`

The OTU was reduced to just the fungal sequences.

`$ filter_taxa_from_otu_table.py -i ~/filtered_otus_oct_dec.biom -p D_2_Fungi -o filtered_fungal_otus_oct_dec.biom`

A third OTU table was created containing only chytrid OTUs.

`$ filter_taxa_from_otu_table.py -i ~/filtered_otus_oct_dec.biom -p D_3__Chytridiomycota -o filtered_chytrid_otus_oct_dec.biom`

Both were filtered to remove OTUs with an abundance <1. 

``

## Generation of phylogenetic tree
The three OTU tables were converted into fasta files using the following scripts. Only the full OTU table is #shown.

`$biom convert -b -i ~/ filtered_otus_oct_dec.biom -o oct_dec_otus.txt`

`$ cat oct_dec_otus.txt | cut -f 1| grep “^denovo” | tr “\n” “ ” > oct_dec_otus2.txt`

`$ perl (-ne ‘if(/^>(\S+)/){$c=grep{/^$1/}qw( [space delimited list of OTUs copy & pasted from oct_dec_otus.txt]>)} print if $c’ ~/rep_set.fasta > ./oct_dec_otus.fasta`

Sequences were then aligned, and a phylogeny was inferred. (Showing only full OTU table.)

`$ align_seqs.py -i ~/oct_dec_otus.fasta -m pynast -t ~/ref_align.fasta -o Alignment/`

`$ filter_alignment.py -i ~/oct_dec_otus_aligned.fasta -o filtered_algn_oct_dec_otus/`

`$make_phylogeny.py -i ~/[alignment filename] -l tree_log.txt -o Tree.tre`

## Taxonomy plots

## Alpha diversity indices

Alpha diversity indices were calculated with the OTU data.

`$ alpha_diversity.py -i ~/filtered_oct_dec_otus.biom -m Shannon, simpson_reciprocal, PD_whole_tree, chao1, observed_species -t Tree.tre -o alpha_indices.txt`

Rarefaction curves were calculated as well. A parameter file listing the alpha diversity indices used above was created first.

`$alpha_rarefaction.py -i ~/filtered_oct_dec_otus.biom -m mapping.txt -t Tree.tre -p alphaparams.txt -o arare/`

## Beta diversity indices & analyses

`beta_diversity_through_plots.py -i oct_dec_otus_wo_ITX3.biom -m Mapping_ID_wm4_corrected.txt -p beta_params -t Tree.tre -o all_otu_beta_diversity_plots`

`beta_diversity_through_plots.py -i chytrid_oct_dec_otus.biom -m Mapping_ID_wm4_corrected.txt -p beta_params -t Tree.tre -o chytrid_otu_beta_diversity_plots`

`beta_diversity.py -i oct_dec_otus_wo_ITX3.biom -m unweighted_unifrac,weighted_unifrac,bray_curtis,bray_curtis_faith -t Tree.tre -o all_beta`

`beta_diversity.py -i chytrid_oct_dec_otus.biom -m unweighted_unifrac,weighted_unifrac,bray_curtis,bray_curtis_faith -t Tree.tre -o chytrid_beta`

`compare_categories.py --method adonis -i ./all_beta/bray_curtis_oct_dec_otus_wo_ITX3.txt -m Mapping_ID_wm4_corrected.txt -c DOB -o all_adonis_bray_curtis_DOB`

`compare_categories.py --method anosim -i ./all_beta/bray_curtis_oct_dec_otus_wo_ITX3.txt -m Mapping_ID_wm4_corrected.txt -c DOB -o all_anosim_bray_curtis_DOB`

