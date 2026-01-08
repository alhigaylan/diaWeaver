// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
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
    // Load input mzML
    // ------------------------------
    MSExperiment raw;
    MzMLFile().load(in, raw);

    // ------------------------------
    // Run DiaWeaver algorithm
    // ------------------------------
    DiaWeaver::WindowMap windows;
    DiaWeaver::determineWindows(raw, windows);

    DiaWeaver::WindowedExperiments ms2_windows;
    DiaWeaver::WindowedExperiments ms1_windows;
    DiaWeaver::WindowedExperiments precursor_windows;

    DiaWeaver::extractMS2Windows(raw, windows, ms2_windows,
                                  save_precursors ? &precursor_windows : nullptr);
    DiaWeaver::extractMS1Windows(raw, windows, ms1_windows);

    // ------------------------------
    // Write outputs
    // ------------------------------
    MzMLFile mzml;

    for (const auto& it : ms2_windows)
    {
      const auto& w = it.first;
      const MSExperiment& exp = it.second;

      String fname = "ms2_" +
        String(w.lower_mz) + "-" + String(w.upper_mz) + "_" +
        String(w.lower_im) + "-" + String(w.upper_im) + ".mzML";

      fname.substitute(".", "-");

      mzml.store(out + "/" + fname, exp);
    }

    for (const auto& it : ms1_windows)
    {
      const auto& w = it.first;
      const MSExperiment& exp = it.second;

      String fname = "ms1_" +
        String(w.lower_mz) + "-" + String(w.upper_mz) + "_" +
        String(w.lower_im) + "-" + String(w.upper_im) + ".mzML";

      fname.substitute(".", "-");

      mzml.store(out + "/" + fname, exp);
    }

    // Write unfragmented precursor files if requested
    if (save_precursors)
    {
      for (const auto& it : precursor_windows)
      {
        const auto& w = it.first;
        const MSExperiment& exp = it.second;

        String fname = "precursor_" +
          String(w.lower_mz) + "-" + String(w.upper_mz) + "_" +
          String(w.lower_im) + "-" + String(w.upper_im) + ".mzML";

        fname.substitute(".", "-");

        mzml.store(out + "/" + fname, exp);
      }
    }

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPDiaWeaver tool;
  return tool.main(argc, argv);
}
