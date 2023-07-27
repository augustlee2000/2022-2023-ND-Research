# Issues with Btagging especially with weaver
I am following almost exactly this [GitHub](https://github.com/alintulu/Run3ScoutingJetTagging/blob/0ee5c8785642df8464692de27055ccbf105a9043/README.md) the goal is to be able to take AODSIM and successfully btag our jets produced by the hlt


## Creating nTuples

I followed almost exactly the GitHub the only differece is that I added muon information from the HLT to the flat tree that is produced. I have already check and having the muon information does not cause (The error shows up with or without the muon data). I have also checked if the dataset I was using was causing the error, it does not affect the results, I check on a TTbar sample and a QCD (the same QCD example provided in the orginal github). I have provided the ttbar 
