# Issues with Btagging especially with weaver
I am following almost exactly this [GitHub](https://github.com/alintulu/Run3ScoutingJetTagging/blob/0ee5c8785642df8464692de27055ccbf105a9043/README.md) the goal is to be able to take AODSIM and successfully btag our jets produced by the hlt


## Creating nTuples

I followed almost exactly the orginal GitHub the only differece is that I added muon information from the HLT to the flat tree. I have already check that having the muon information does not cause the error (The error shows up with or without the muon data). I have also checked if the dataset I am using was causing the error, it does not affect the results, I check on a TTbar sample and a QCD (the same QCD example provided in the orginal github). I have provided the ttbar.root, QCD.root, AK4.py and AK4JetNTupleProducer.cc which are the only files that were changed to make the current .root files

## Using weaver 

This is where the problems starts. I now take the two different .root files and run them locally using this command script (train_test.sh) (we had issues running gpu using HTCondor), this method is not supose to need a gpu. I should note that an active conda enviroment is need to run these commands, I followed this [GitHub](https://github.com/hqucms/weaver-core/#set-up-your-environment) to set up my environment. When I run the test_TTBar.root file uing the .sh file I get this [error log](https://github.com/augustlee2000/2022-2023-ND-Research/blob/main/Issue_with_code/test_TTBAR.root_error_log.txt). When I run QCD we get this [error log](https://github.com/augustlee2000/2022-2023-ND-Research/blob/main/Issue_with_code/test_QCD.root_error_log.txt).

### Where I think what the error is

In both error logs we have two different errors "IndexError: list index out of range" and "RuntimeError: Inconsistent label definition: some of the entries are assigned to multiple classes!". I have been focusing my efforts on the "RuntimeError: Inconsistent label definition: some of the entries are assigned to multiple classes!" After following the Traceback I belive that the code has an internal check to make sure an event only has one label associated with, and from my understanding the label are defined in the .yaml file used by weaver where is [here]() 
