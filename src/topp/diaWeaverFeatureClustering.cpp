// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/DiaWeaverFeatureClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

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

@brief Aligns diaWeaver pseudo-spectrum mzML files across experiments using pose clustering.

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

Takes diaWeaver pseudo-spectrum mzML files (one per experiment) and aligns their
retention time scales using pose clustering.

Each MS2 spectrum is converted to a feature using its precursor RT, m/z, charge,
and intensity. If the input contains ion mobility data, the precursor drift time
is stored as feature metadata ("ion_mobility").

A reference map is selected automatically (largest by MS2 count) or can be
specified explicitly. Every other map is then aligned to the reference independently.

Output is one featureXML per input file (aligned precursor features) and
optionally one trafoXML per file.

<B>The command line parameters of this tool are:</B> @n
@verbinclude TOPP_diaWeaverFeatureClustering.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_diaWeaverFeatureClustering.html
*/

/// @cond TOPPCLASSES

class TOPPdiaWeaverFeatureClustering :
  public TOPPBase
{

public:
  TOPPdiaWeaverFeatureClustering() :
    TOPPBase("diaWeaverFeatureClustering",
      "Aligns diaWeaver pseudo-spectrum mzML files across experiments using pose clustering.",
      false)
  {}

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(),
      "Input diaWeaver pseudo-spectrum mzML files (one per experiment).", true);
    setValidFormats_("in", {"mzML"});

    registerOutputFileList_("out", "<files>", StringList(),
      "Output featureXML files with aligned precursor features (one per input). "
      "Must match the number of input files if provided.", false);
    setValidFormats_("out", {"featureXML"});

    registerOutputFileList_("trafo_out", "<files>", StringList(),
      "Output trafoXML files describing the RT transformations (one per input). "
      "Must match the number of input files if provided.", false);
    setValidFormats_("trafo_out", {"trafoXML"});

    registerIntOption_("reference:index", "<index>", 0,
      "1-based index of the reference file from the 'in' list. "
      "0 (default) selects the file with the most MS2 spectra automatically.", false);
    setMinInt_("reference:index", 0);

    registerInputFile_("reference:file", "<file>", "",
      "External reference mzML file. Overrides reference:index if provided.", false);
    setValidFormats_("reference:file", {"mzML"});

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
    StringList in_files   = getStringList_("in");
    StringList out_files  = getStringList_("out");
    StringList out_trafos = getStringList_("trafo_out");

    if (in_files.empty())
    {
      writeLogError_("No input files provided.");
      return ILLEGAL_PARAMETERS;
    }
    if (!out_files.empty() && out_files.size() != in_files.size())
    {
      writeLogError_("Number of 'out' files must match number of 'in' files.");
      return ILLEGAL_PARAMETERS;
    }
    if (!out_trafos.empty() && out_trafos.size() != in_files.size())
    {
      writeLogError_("Number of 'trafo_out' files must match number of 'in' files.");
      return ILLEGAL_PARAMETERS;
    }
    if (in_files.size() == 1)
    {
      OPENMS_LOG_WARN << "Only one input file provided to diaWeaverFeatureClustering." << std::endl;
    }

    DiaWeaverFeatureClustering algorithm;
    Param algo_params = getParam_().copy("algorithm:", true);
    algorithm.setParameters(algo_params);
    algorithm.setLogType(log_type_);

    // -----------------------------------------------------------------------
    // Reference selection
    // -----------------------------------------------------------------------
    Size reference_index = getIntOption_("reference:index");
    String reference_file = getStringOption_("reference:file");

    if (!reference_file.empty())
    {
      // explicit external reference
      FeatureMap map_ref;
      DiaWeaverFeatureClustering::loadMzMLAsFeatureMap(reference_file, map_ref);
      algorithm.setReference(map_ref);
      reference_index = in_files.size(); // sentinel: no in-list reference
    }
    else if (reference_index > 0)
    {
      --reference_index; // convert 1-based parameter to 0-based index
      FeatureMap map_ref;
      DiaWeaverFeatureClustering::loadMzMLAsFeatureMap(in_files[reference_index], map_ref);
      algorithm.setReference(map_ref);
    }
    else // reference_index == 0: auto-select largest map by MS2 count
    {
      OPENMS_LOG_INFO << "Picking reference by MS2 spectrum count ..." << std::flush;
      Size max_count = 0;
      for (Size i = 0; i < in_files.size(); ++i)
      {
        FeatureMap m;
        DiaWeaverFeatureClustering::loadMzMLAsFeatureMap(in_files[i], m);
        if (m.size() > max_count)
        {
          max_count = m.size();
          reference_index = i;
        }
      }
      OPENMS_LOG_INFO << " done (selected index " << reference_index << ")" << std::endl;

      FeatureMap map_ref;
      DiaWeaverFeatureClustering::loadMzMLAsFeatureMap(in_files[reference_index], map_ref);
      algorithm.setReference(map_ref);
    }

    // -----------------------------------------------------------------------
    // Align each input map to the reference
    // -----------------------------------------------------------------------
    ProgressLogger plog;
    plog.setLogType(log_type_);

    vector<TransformationDescription> transformations(in_files.size());

    plog.startProgress(0, in_files.size(), "Aligning diaWeaver feature maps");
    Size progress = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int i = 0; i < static_cast<int>(in_files.size()); ++i)
    {
      FeatureMap map;
      DiaWeaverFeatureClustering::loadMzMLAsFeatureMap(in_files[i], map);

      TransformationDescription trafo;
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
        FileHandler().storeFeatures(out_files[i], map, {FileTypes::FEATUREXML}, log_type_);
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
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPdiaWeaverFeatureClustering tool;
  return tool.main(argc, argv);
}

/// @endcond
