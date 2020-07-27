#!/usr/bin/env python

source activate imc
python hcs_data_processing.py \
    -i test_data/HCS-input-form.xlsx \
    -r test_data/PlateMaps_Hamilton_June\ 19\,\ 2020.xls \
    -m spectramax test_data/SBB_20200624_assay13A_final_ctg_reading.txt \
    -m phenix-day1 test_data/SBB_Assay_13A_20200619_calcein_viol_SBB_20200720_day1_texture_identified_orgs.txt CM,CO,CR \
    -m phenix-day5 test_data/SBB_assay13A_20200624_day5\ _\ SBB_20200709_class_A_B_.txt EQ,ES,IU,IW \
    -f -v \
    -o testfile.csv \
    -p testfile_plot.pdf

python gr50.py \
  -i testfile.csv \
  -g0 'phenix-day1.*Sum.*' \
  -g1 'phenix-day5.*Sum.*' \
  -p testgr50
