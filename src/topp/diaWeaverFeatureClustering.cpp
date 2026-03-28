// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/DiaWeaverFeatureClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_diaWeaverFeatureClustering diaWeaverFeatureClustering

@brief Aligns diaWeaver FeatureMaps across experiments using pose clustering.

<CENTER>
  <table>
    <tr>
      <th ALIGN = "center"> potential predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=2> &rarr; diaWeaverFeatureClustering &rarr;</td>
      <th ALIGN = "center"> potential successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_diaWeaver </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeledQT </td>
    </tr>
  </table>
</CENTER>

Takes diaWeaver pseudo MS2 spectra represented as featureXML files (one per
experiment) and aligns their retention time scales using pose clustering. The
algorithm selects a reference map (largest by default), then computes an affine
RT transformation for every other map to align it to the reference.

Each feature's (RT, m/z, charge) corresponds directly to a diaWeaver pseudo
spectrum precursor. After alignment, features with matching m/z and charge
across experiments will have comparable RTs, enabling cross-experiment grouping
with FeatureLinkerUnlabeledQT.

<B>The command line parameters of this tool are:</B> @n
@verbinclude TOPP_diaWeaverFeatureClustering.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_diaWeaverFeatureClustering.html
*/

/// @cond TOPPCLASSES

class TOPPdiaWeaverFeatureClustering :
  public TOPPMapAlignerBase
{

public:
  TOPPdiaWeaverFeatureClustering() :
    TOPPMapAlignerBase("diaWeaverFeatureClustering",
      "Aligns diaWeaver FeatureMaps across experiments using pose clustering.")
  {}

protected:
  void registerOptionsAndFlags_() override
  {
    TOPPMapAlignerBase::registerOptionsAndFlagsMapAligners_("featureXML", REF_RESTRICTED);
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& section) const override
  {
    if (section == "algorithm")
    {
      DiaWeaverFeatureClustering algo;
      return algo.getParameters();
    }
    return Param();
  }

  ExitCodes main_(int, const char**) override
  {
    ExitCodes ret = TOPPMapAlignerBase::checkParameters_();
    if (ret != EXECUTION_OK) return ret;

    DiaWeaverFeatureClustering algorithm;
    Param algo_params = getParam_().copy("algorithm:", true);
    algorithm.setParameters(algo_params);
    algorithm.setLogType(log_type_);

    StringList in_files  = getStringList_("in");
    StringList out_files = getStringList_("out");
    StringList out_trafos = getStringList_("trafo_out");

    if (in_files.size() == 1)
    {
      OPENMS_LOG_WARN << "Only one input file provided to diaWeaverFeatureClustering." << std::endl;
    }

    Size reference_index = getIntOption_("reference:index");
    String reference_file = getStringOption_("reference:file");

    String ref_file;
    if (!reference_file.empty())
    {
      ref_file = reference_file;
      reference_index = in_files.size(); // invalid — points past end
    }
    else if (reference_index > 0)
    {
      ref_file = in_files[--reference_index]; // 1-based in params, 0-based here
    }
    else // reference_index == 0: auto-select largest map
    {
      OPENMS_LOG_INFO << "Picking reference by feature count ..." << std::flush;
      Size max_count(0);
      FeatureXMLFile f;
      for (Size i = 0; i < in_files.size(); ++i)
      {
        Size s = f.loadSize(in_files[i]);
        if (s > max_count)
        {
          max_count = s;
          reference_index = i;
        }
      }
      OPENMS_LOG_INFO << " done" << std::endl;
      ref_file = in_files[reference_index];
    }

    // Load and set reference map
    {
      FileHandler fh;
      fh.getFeatOptions().setLoadConvexHull(false);
      fh.getFeatOptions().setLoadSubordinates(false);
      FeatureMap map_ref;
      fh.loadFeatures(ref_file, map_ref, {FileTypes::FEATUREXML}, log_type_);
      algorithm.setReference(map_ref);
    }

    FileHandler f_fxml;
    if (out_files.empty())
    {
      f_fxml.getFeatOptions().setLoadConvexHull(false);
      f_fxml.getFeatOptions().setLoadSubordinates(false);
    }

    ProgressLogger plog;
    plog.setLogType(log_type_);

    vector<TransformationDescription> transformations(in_files.size());

    plog.startProgress(0, in_files.size(), "Aligning diaWeaver feature maps");
    Size progress(0);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int i = 0; i < static_cast<int>(in_files.size()); ++i)
    {
      TransformationDescription trafo;
      FeatureMap map;

      FileHandler f_fxml_tmp;
      f_fxml_tmp.getFeatOptions() = f_fxml.getFeatOptions();
      f_fxml_tmp.loadFeatures(in_files[i], map);

      if (i == static_cast<int>(reference_index))
      {
        trafo.fitModel("identity");
      }
      else
      {
        try
        {
          algorithm.align(map, trafo);
        }
        catch (Exception::IllegalArgument& e)
        {
          OPENMS_LOG_ERROR << "Aligning " << in_files[i] << " failed. "
                           << "No transformation applied (RT unchanged)." << endl;
          writeLogError_("Illegal argument (" + String(e.getName()) + "): " + String(e.what()) + ".");
          trafo.fitModel("identity");
        }
      }

      if (!out_files.empty())
      {
        MapAlignmentTransformer::transformRetentionTimes(map, trafo);
        addDataProcessing_(map, getProcessingInfo_(DataProcessing::ALIGNMENT));
        f_fxml_tmp.storeFeatures(out_files[i], map, {FileTypes::FEATUREXML}, log_type_);
      }

      transformations[i] = trafo;

      if (!out_trafos.empty())
      {
        FileHandler().storeTransformations(out_trafos[i], trafo, {FileTypes::TRANSFORMATIONXML});
      }

#ifdef _OPENMP
#pragma omp critical (diaWeaverFC_Progress)
#endif
      {
        plog.setProgress(++progress);
      }
    }

    plog.endProgress();

    // Transform optional spectra files
    StringList in_spectra  = getStringList_("in_spectra_files");
    StringList out_spectra = getStringList_("out_spectra_files");
    transformSpectraFiles_(in_spectra, out_spectra, transformations, false);

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPdiaWeaverFeatureClustering tool;
  return tool.main(argc, argv);
}

/// @endcond
