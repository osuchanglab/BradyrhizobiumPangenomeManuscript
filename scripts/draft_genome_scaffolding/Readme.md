## Scripts for the scaffolding of draft genomes against completely assembled genomes, as well as extracting regions corresponding to symICEs.

#### Align contigs from draft genomes to the most similar complete/hybrid assembled genome using contiguator:

map_draft_to_complete.sh

#### Select the best reference strain for each draft assembly based on visual inspection of contiguator output and/or the most closely related symbiosis genes. Save to a table file.

#### Make a list of contigs in the order and direction they map to the complete genome. Uses the above mentioned table:

getcontigorder.sh

#### Make a scaffolded version of the draft assembly based on the contiguator alignment:

make_combo_gbk.sh

#### Identify the predicted att sites for each symICE type using BLASTN and save the filtered results to a table ("SI_att_position_table") with the columns:

strain	ICEtype	startposition	endposition

The second column is used for labeling output files.

#### Extract each symICE element from the scaffolded genomes:

extract_islands.sh


