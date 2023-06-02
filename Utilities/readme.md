This folder containts all utilities scrips that I or anyone might need to use ROOT

Instructions on how to run each file and their uses:

# **CheckFile.py**

This program is made to take a mini or aod root file and provided all the tree that your could possible need
## Steps

1) Launch up your favorite instance of CMSSW
2) type in the command: edmDumpEventContent example.root
3) find the Header and Lable from the edmDump
4) Now go into the file and replace line 58 59 with the header and label you want to look at
5) If you are only intrested in the possible branches everything after line 64 can be commented out
6) If you are intrested in a very basic histogram line 63 can be commented out and you can replace line 83 from hist_array[0].Fill(ig.eta()) to hist_array[0].Fill(ig.YOURPARAMETER())
7) If you want both do not comment anything
8) Save your updated python script
9) type in the command: python3 CheckFile.py example.root
10) Enjoy

