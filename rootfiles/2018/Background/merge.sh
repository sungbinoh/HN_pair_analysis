rm HN_pair_all_SkimTree_LRSMHighPt_TT_powheg.root
rm HN_pair_all_SkimTree_LRSMHighPt_VV.root
rm HN_pair_all_SkimTree_LRSMHighPt_DYJets_MG_HT.root
rm HN_pair_all_SkimTree_LRSMHighPt_WJets_MG_HT.root
rm HN_pair_all_SkimTree_LRSMHighPt_fake.root

cp HN_pair_all_SkimTree_LRSMHighPt_TTLL_powheg.root HN_pair_all_SkimTree_LRSMHighPt_TT_powheg.root
hadd HN_pair_all_SkimTree_LRSMHighPt_fake.root HN_pair_all_SkimTree_LRSMHighPt_TTLJ_powheg.root HN_pair_all_SkimTree_LRSMHighPt_WJets_MG_HT-*.root HN_pair_all_SkimTree_LRSMHighPt_QCD.root
hadd HN_pair_all_SkimTree_LRSMHighPt_DYJets_MG_HT.root HN_pair_all_SkimTree_LRSMHighPt_DYJets_MG_HT-*.root
hadd HN_pair_all_SkimTree_LRSMHighPt_VV.root HN_pair_all_SkimTree_LRSMHighPt_ZZ_pythia.root HN_pair_all_SkimTree_LRSMHighPt_WZ_pythia.root HN_pair_all_SkimTree_LRSMHighPt_WW_pythia.root HN_pair_all_SkimTree_LRSMHighPt_WWW.root HN_pair_all_SkimTree_LRSMHighPt_WWZ.root HN_pair_all_SkimTree_LRSMHighPt_WZZ.root HN_pair_all_SkimTree_LRSMHighPt_ZZZ.root HN_pair_all_SkimTree_LRSMHighPt_ttZ.root HN_pair_all_SkimTree_LRSMHighPt_ttW.root HN_pair_all_SkimTree_LRSMHighPt_SingleTop_tW_top_NoFullyHad.root HN_pair_all_SkimTree_LRSMHighPt_SingleTop_tW_antitop_NoFullyHad.root HN_pair_all_SkimTree_LRSMHighPt_SingleTop_tch_top_Incl.root HN_pair_all_SkimTree_LRSMHighPt_SingleTop_tch_antitop_Incl.root HN_pair_all_SkimTree_LRSMHighPt_SingleTop_sch_Lep.root


#hadd HN_pair_all_SkimTree_LRSMHighPt_TT_powheg.root HN_pair_all_SkimTree_LRSMHighPt_TTLL_powheg.root HN_pair_all_SkimTree_LRSMHighPt_TTLJ_powheg.root #HN_pair_all_SkimTree_LRSMHighPt_TTJJ_powheg.root
#hadd HN_pair_all_SkimTree_LRSMHighPt_VV.root HN_pair_all_SkimTree_LRSMHighPt_ZZ_pythia.root HN_pair_all_SkimTree_LRSMHighPt_WZ_pythia.root HN_pair_all_SkimTree_LRSMHighPt_WW_pythia.root HN_pair_all_SkimTree_LRSMHighPt_WWW.root HN_pair_all_SkimTree_LRSMHighPt_WWZ.root HN_pair_all_SkimTree_LRSMHighPt_WZZ.root HN_pair_all_SkimTree_LRSMHighPt_ZZZ.root HN_pair_all_SkimTree_LRSMHighPt_ttZ.root HN_pair_all_SkimTree_LRSMHighPt_ttW.root HN_pair_all_SkimTree_LRSMHighPt_SingleTop_tW_top_NoFullyHad.root HN_pair_all_SkimTree_LRSMHighPt_SingleTop_tW_antitop_NoFullyHad.root HN_pair_all_SkimTree_LRSMHighPt_SingleTop_tch_top_Incl.root HN_pair_all_SkimTree_LRSMHighPt_SingleTop_tch_antitop_Incl.root HN_pair_all_SkimTree_LRSMHighPt_SingleTop_sch_Lep.root #HN_pair_all_SkimTree_LRSMHighPt_TTG.root
#hadd HN_pair_all_SkimTree_LRSMHighPt_DYJets_MG_HT.root HN_pair_all_SkimTree_LRSMHighPt_DYJets_MG_HT-*.root
#hadd HN_pair_all_SkimTree_LRSMHighPt_WJets_MG_HT.root HN_pair_all_SkimTree_LRSMHighPt_WJets_MG_HT-*.root
