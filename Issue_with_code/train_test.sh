weaver \
    --data-train test.root \
    --data-config /afs/cern.ch/user/a/aulee/BScouting/CMSSW_13_1_0_pre1/src/Run3ScoutingJetTagging/Training/AK4/flavour/data/flavour.yaml \
    --network-config /afs/cern.ch/user/a/aulee/BScouting/CMSSW_13_1_0_pre1/src/Run3ScoutingJetTagging/Training/AK4/flavour/networks/flavour.py \
    --model-prefix 'output/ak4_flavour_{auto}/net' \
    --gpus '' --demo --num-epochs 1 \
    --log 'output/ak4_flavour_{auto}.log'
