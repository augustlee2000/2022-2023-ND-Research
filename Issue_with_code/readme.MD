# Issues with Btagging especially with weaver
I am following almost exactly this [GitHub](https://github.com/alintulu/Run3ScoutingJetTagging/blob/0ee5c8785642df8464692de27055ccbf105a9043/README.md) the goal is to be able to take AODSIM and successfully btag our jets produced by the hlt


## Creating nTuples

I followed almost exactly the GitHub the only differece is that I added muon information from the HLT to the flat tree that is produced. I have already check and having the muon information does not cause (The error shows up with or without the muon data). I have also checked if the dataset I was using was causing the error, it does not affect the results, I check on a TTbar sample and a QCD (the same QCD example provided in the orginal github). I have provided the ttbar, QCD, AK4.py and AK$JetNTupleProducer.cc which are the only files changed to make the current .root files

## Using weaver 

This is where the problems starts. I now take the two different .root files and run them locally using this command script (train_test.sh) after some issues with gpu, this method is suppose to not need a gpu. I should note that a active conda enviroment is need to run these commands, I followed this [GitHub](https://github.com/hqucms/weaver-core/#set-up-your-environment) to set up my environment.
