
evgenConfig.description = "Pythia8 dijet samples with a fake taus filter"
evgenConfig.keywords = ["tau", "performance", "jets", "tau", "QCD"]
evgenConfig.process  = ""
evgenConfig.generators =  ["Pythia8"]
evgenConfig.contact  = [ "quentin.buat@cern.ch" ]

include("MC15JobOptions/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")

genSeq.Pythia8.Commands += [
    "HardQCD:all = on",
    # "PromptPhoton:qg2qgamma = on",
    # "PromptPhoton:qqbar2ggamma = on",
    # "WeakSingleBoson:all = on",
    "Top:gg2ttbar = on",
    "Top:qqbar2ttbar = on",
    "PhaseSpace:pTHatMin = 15",         # PtHat and mHat thresholds interpolated from DC14 numbers
    "PhaseSpace:mHatMin = 30",
]          # to be used for performance studies only



if not hasattr( filtSeq, "FakeTauFilter" ):
    from GeneratorFilters.GeneratorFiltersConf import FakeTauFilter
    filtSeq += FakeTauFilter()
    pass

filtSeq.FakeTauFilter.UseDr = False
filtSeq.FakeTauFilter.MinTracksCore = 0
