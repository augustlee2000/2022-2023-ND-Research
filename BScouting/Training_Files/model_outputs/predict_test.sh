#!/bin/bash

weaver \
	--predict \
    --data-test "/hadoop/store/user/alee43/bscouting_slimmed/TTBar_120240209_1435/b_tagging/output_1.root" \
	--data-config /afs/crc.nd.edu/user/a/alee43/BScouting/CMSSW_13_1_0_pre1/src/Run3ScoutingJetTagging/Training/AK4/flavour/data/flavour.yaml \
    --network-config /afs/crc.nd.edu/user/a/alee43/BScouting/CMSSW_13_1_0_pre1/src/Run3ScoutingJetTagging/Training/AK4/flavour/networks/flavour.py \
	--model-prefix /afs/crc.nd.edu/user/a/alee43/BScouting/CMSSW_13_1_0_pre1/src/Run3ScoutingJetTagging/Training/AK4/flavour/output/ak4_flavour_20240314-123152_flavour_ranger_lr0.005_batch512/net_best_epoch_state.pt \
	--gpus '' --batch-size 512 \
	--predict-output /afs/crc.nd.edu/user/a/alee43/BScouting/CMSSW_13_1_0_pre1/src/Run3ScoutingJetTagging/Training/AK4/flavour/output/ak4_flavour_20240314-123152_flavour_ranger_lr0.005_batch512/pred.root