#add_subdirectory(NewSubdirectory)
add_subdirectory(ndkAna)
#add_subdirectory(ndkReco)


# Respect the order specified in lar_build, which is:
#
#     "SoftRelTools",
#     "SRT_LAR",
#     "SimpleTypesAndConstants",
#     "larcore_Geometry",
#     "larcore_SummaryData",
#     "lardata_RawData",
#     "lardata_RecoBase",
#     "lardata_OpticalDetectorData",
#     "lardata_AnalysisBase",
#     "lardata_Utilities",
#     "larsim_Simulation",
#     "larsim_MCCheater",
#     "lardata_AnalysisAlg",
#     "lardata_RecoObjects",
#     "larreco_RecoAlg",
#     "larsim_LArG4",
#     "larsim_PhotonPropagation",
#     "larsim_DetSim",
#     "larsim_EventGenerator",
#     "larevt_Filters",
#     "CalData",
#     "larreco_ClusterFinder",
#     "larreco_Genfit",
#     "larreco_HitFinder",
#     "larreco_TrackFinder",
#     "larreco_VertexFinder",
#     "larreco_EventFinder",
#     "larsim_TriggerAlgo",
#     "larana_OpticalDetector",
#     "larana_Calorimetry", 
#     "larreco_ShowerFinder",
#     "lareventdisplay_EventDisplay",
#     "larana_ParticleIdentification",
#     "larexamples_AnalysisExample",
#     "LArPandoraAlgorithms",
#     "larpandora_LArPandoraInterface"
