// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/APPLICATIONS/diaWeaver.h>

using namespace OpenMS;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_diaWeaver diaWeaver

@brief Splits a DIA mzML file into per-window MS1 and MS2 mzML files.

This tool extracts DIA windows from ion mobility DIA data and writes separate mzML files
for each precursor isolation window. It processes both MS1 and MS2 spectra, filtering
peaks by m/z and ion mobility ranges defined by the DIA acquisition windows.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_diaWeaver.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_diaWeaver.html
*/

class TOPPDiaWeaver :
  public TOPPBase
{
public:
  TOPPDiaWeaver() :
    TOPPBase(
      "diaWeaver",
      "Split a DIA mzML into per-window MS1 and MS2 mzML files",
      false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_(
      "in",
      "<file>",
      "",
      "Input DIA mzML file",
      true);
    setValidFormats_("in", {"mzML"});

    registerOutputPrefix_(
      "out",
      "<prefix>",
      "",
      "Output directory prefix (default: <input>_diaWindows)");

    registerFlag_(
      "save_unfragmented_precursors",
      "If set, save peaks within precursor isolation window in MS2 and apply precursor detection algorithm");
  }

  ExitCodes main_(int, const char**) override
  {
    const String in = getStringOption_("in");
    const bool save_precursors = getFlag_("save_unfragmented_precursors");

    String out = getStringOption_("out");
    if (out.empty())
    {
      out = File::path(in) + "/" +
            File::basename(in) + "_diaWindows";
    }
    File::makeDir(out);

    // ------------------------------
    // Load input mzML using on-disk access
    // ------------------------------
    OnDiscMSExperiment raw;
    if (!raw.openFile(in))
    {
      OPENMS_LOG_ERROR << "Failed to open file as indexed mzML." << std::endl;
      return INPUT_FILE_NOT_FOUND;
    }

    // ------------------------------
    // Determine DIA windows and IM info
    // ------------------------------
    DiaWeaver::WindowMap windows;
    DiaWeaver::determineWindows(raw, windows);

    // Determine IM info once upfront
    DiaWeaver::IMInfo im_info = DiaWeaver::determineIMInfo(raw, windows);

    OPENMS_LOG_INFO << "Processing " << windows.size() << " DIA windows..." << std::endl;

    // ------------------------------
    // Process each window incrementally: extract, write, release memory
    // ------------------------------
    MzMLFile mzml;
    MSExperiment ms2_exp;
    MSExperiment ms1_exp;
    MSExperiment precursor_exp;

    Size window_idx = 0;
    for (const auto& it : windows)
    {
      const DiaWeaver::DIAWindow& w = it.first;
      const std::vector<Size>& indices = it.second;
      ++window_idx;

      OPENMS_LOG_INFO << "Processing window " << window_idx << "/" << windows.size()
                      << " (m/z: " << w.lower_mz << "-" << w.upper_mz << ")" << std::endl;

      // Build filename base
      String fname_base =
        String(w.lower_mz) + "-" + String(w.upper_mz) + "_" +
        String(w.lower_im) + "-" + String(w.upper_im) + ".mzML";
      fname_base.substitute(".", "-");

      // Extract and write MS2
      DiaWeaver::extractSingleMS2Window(raw, w, indices, im_info, ms2_exp,
                                         save_precursors ? &precursor_exp : nullptr);
      if (!ms2_exp.empty())
      {
        mzml.store(out + "/ms2_" + fname_base, ms2_exp);
      }

      // Write precursors if requested
      if (save_precursors && !precursor_exp.empty())
      {
        mzml.store(out + "/precursor_" + fname_base, precursor_exp);
      }

      // Extract and write MS1
      DiaWeaver::extractSingleMS1Window(raw, w, im_info, ms1_exp);
      if (!ms1_exp.empty())
      {
        mzml.store(out + "/ms1_" + fname_base, ms1_exp);
      }
    }

    OPENMS_LOG_INFO << "Finished processing all windows." << std::endl;

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPDiaWeaver tool;
  return tool.main(argc, argv);
}
