// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest, Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/ClusterMassTracesByPrecursor.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusMap.h>


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_ClusterMassTracesByPrecursor ClusterMassTracesByPrecursor

@brief Identifies precursor mass traces and tries to correlate them with fragment ion mass traces in SWATH maps.

This algorithm will try to correlate the masstraces to find co-eluting traces and cluster them.

This program looks at mass traces in a precursor MS1 map and tries to
correlate them with features found in the corresponding MS2 map based on
their elution profile. It uses

 - the mass traces from the MS1 in consensusXML format [note this is an unintended use of the consesusXML format to also store intensities]
 - the mass traces from the MS2 (SWATH map)

 It does a separate correlation analysis on the MS1 and the MS2 map,
 both produces a set of pseudo spectra.
 In a second (optional) step, the MS2 pseudo spectra are correlated with
 the MS1 traces and the most likely precursor is assigned to the pseudo
 spectrum.
  
It is based on the following papers:
ETISEQ -- an algorithm for automated elution time ion sequencing of concurrently fragmented peptides for mass spectrometry-based proteomics
  BMC Bioinformatics 2009, 10:244 doi:10.1186/1471-2105-10-244 ; http://www.biomedcentral.com/1471-2105/10/244
  they use FFT to correlate and then use lag of at least 1 scan and pearson correlation of 0.7 to assign precursors to product ions
  If one fragment matches to multiple precursors, it is assigned to all of them. If it doesn't match any, it is assigned to all

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ClusterMassTracesByPrecursor.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ClusterMassTracesByPrecursor.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace std;
using namespace OpenMS;

class TOPPCorrelateMasstraces
  : public TOPPBase, 
    public ProgressLogger
{

 public:

  TOPPCorrelateMasstraces()
    : TOPPBase("ClusterMassTracesByPrecursor", "Correlate precursor masstraces with fragment ion masstraces in SWATH maps based on their elution profile.")
  {
  }

 protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in_ms1","<file>","","MS1 mass traces");
    setValidFormats_("in_ms1",ListUtils::create<String>("consensusXML"));

    registerInputFile_("in_swath","<file>","","MS2 / SWATH mass traces");
    setValidFormats_("in_swath",ListUtils::create<String>("consensusXML"));

    registerOutputFile_("out","<file>","","output file");
    setValidFormats_("out",ListUtils::create<String>("mzML"));

    // registerFlag_("ms1_centric","MS1 centric - find MS1 features first and then add MS2s (MSE like)");
    registerFlag_("assign_unassigned_to_all","Assign unassigned MS2 fragments to all precursors (only for ms1_centrif)");

    registerDoubleOption_("min_pearson_correlation", "<double>", 0.7, "Minimal pearson correlation score to match elution profiles to each other.", false); // try 0.3, 0.5 and 0.7
    registerIntOption_("max_lag", "<number>", 1, "Maximal lag (e.g. by how many spectra the peak may be shifted at most). This parameter will depend on your chromatographic setup but a number between 1 and 3 is usually sensible.", false);
    registerIntOption_("min_nr_ions", "<number>", 3, "Minimal number of ions to report a spectrum.", false);
    registerDoubleOption_("max_rt_apex_difference", "<double>", 5.0, "Maximal difference of the apex in retention time (in seconds). This is a hard parameter, all profiles further away will not be considered at all.", false);

    registerDoubleOption_("swath_lower", "<double>", 0.0, "Swath lower isolation window", false);
    registerDoubleOption_("swath_upper", "<double>", 0.0, "Swath upper isolation window", false);
  }

 public:

  ExitCodes main_(int , const char**) override
  {
    setLogType(log_type_);

    String ms1 = getStringOption_("in_ms1");
    String in_swath = getStringOption_("in_swath");
    String out = getStringOption_("out");

    double swath_lower = getDoubleOption_("swath_lower");
    double swath_upper = getDoubleOption_("swath_upper");

    // Load input:
    // - MS1 feature map containing the MS1 mass traces
    // - MS2 feature map containing the MS2 (SWATH) mass traces
    ConsensusMap MS1_feature_map;
    ConsensusMap MS2_feature_map;
    FileHandler().loadConsensusFeatures(ms1, MS1_feature_map, {FileTypes::CONSENSUSXML}, log_type_);
    FileHandler().loadConsensusFeatures(in_swath, MS2_feature_map, {FileTypes::CONSENSUSXML}, log_type_);
    OPENMS_LOG_INFO << "Loaded consensus maps: " << MS1_feature_map.size() << " MS1 traces, "
                    << MS2_feature_map.size() << " MS2 traces" << endl;

    // Set up clustering parameters
    Param cluster_param;
    cluster_param.setValue("min_pearson_correlation", getDoubleOption_("min_pearson_correlation"));
    cluster_param.setValue("max_lag", getIntOption_("max_lag"));
    cluster_param.setValue("max_rt_apex_difference", getDoubleOption_("max_rt_apex_difference"));
    cluster_param.setValue("min_nr_ions", getIntOption_("min_nr_ions"));
    cluster_param.setValue("assign_unassigned_to_all", getFlag_("assign_unassigned_to_all") ? "true" : "false");

    // Run clustering
    MSExperiment pseudo_spectra;
    ClusterMassTracesByPrecursor clusterFragments;
    clusterFragments.setParameters(cluster_param);
    clusterFragments.setLogType(log_type_);
    clusterFragments.run(MS1_feature_map, MS2_feature_map, swath_lower, swath_upper, pseudo_spectra);

    // Store output
    FileHandler().storeExperiment(out, pseudo_spectra, {FileTypes::MZML}, log_type_);

    return EXECUTION_OK;
  }
};

int main( int argc, const char** argv )
{
  TOPPCorrelateMasstraces tool;
  return tool.main(argc,argv);
}

///@endcond
