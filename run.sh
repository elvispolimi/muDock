# add  ncu --set full -f -o my_rep to profile 

./build/application/muDock --protein data/1fkb/1fkb_protein.pdb --ligand ../helpersMuDock/the-molecular-docking-dataset/out_sweep/out_10000/dataset.adtmol2 --use CUDA:GPU:0
