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

#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

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

    registerIntOption_(
      "threads",
      "<n>",
      1,
      "Number of threads to use for parallel window processing (default: 1)",
      false);
    setMinInt_("threads", 1);
  }

  ExitCodes main_(int, const char**) override
  {
    const String in = getStringOption_("in");
    const bool save_precursors = getFlag_("save_unfragmented_precursors");
#ifdef _OPENMP
    const int num_threads = getIntOption_("threads");
#endif

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

    // Convert map to vector for OpenMP indexed access
    std::vector<std::pair<DiaWeaver::DIAWindow, std::vector<Size>>> window_vec(
      windows.begin(), windows.end());
    const Size total_windows = window_vec.size();

    OPENMS_LOG_INFO << "Processing " << total_windows << " DIA windows";
#ifdef _OPENMP
    OPENMS_LOG_INFO << " using " << num_threads << " thread(s)";
#endif
    OPENMS_LOG_INFO << "..." << std::endl;

    // ------------------------------
    // Process windows in parallel
    // Each thread opens its own file handle for thread-safe reading
    // ------------------------------
    Size processed = 0;

#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (SignedSize idx = 0; idx < static_cast<SignedSize>(total_windows); ++idx)
    {
      const DiaWeaver::DIAWindow& w = window_vec[idx].first;
      const std::vector<Size>& indices = window_vec[idx].second;

      // Each thread needs its own file reader for thread-safe access
      OnDiscMSExperiment thread_raw;
      thread_raw.openFile(in);

      // Thread-local buffers
      MzMLFile mzml;
      MSExperiment ms2_exp;
      MSExperiment ms1_exp;
      MSExperiment precursor_exp;

      // Build human-readable filename: mz<lower>-<upper>_im<lower>-<upper>.mzML
      String fname_base = "mz" + String(static_cast<int>(std::round(w.lower_mz))) +
                          "-" + String(static_cast<int>(std::round(w.upper_mz)));
      if (w.hasIonMobility())
      {
        String im_lower = String::number(w.lower_im, 2);
        String im_upper = String::number(w.upper_im, 2);
        im_lower.substitute(".", "p");
        im_upper.substitute(".", "p");
        fname_base += "_im" + im_lower + "-" + im_upper;
      }
      fname_base += ".mzML";

      // Extract and write MS2
      DiaWeaver::extractSingleMS2Window(thread_raw, w, indices, im_info, ms2_exp,
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
      DiaWeaver::extractSingleMS1Window(thread_raw, w, im_info, ms1_exp);
      if (!ms1_exp.empty())
      {
        mzml.store(out + "/ms1_" + fname_base, ms1_exp);
      }

      // Progress logging (thread-safe)
#ifdef _OPENMP
#pragma omp critical (progress_log)
#endif
      {
        ++processed;
        OPENMS_LOG_INFO << "Processed window " << processed << "/" << total_windows
                        << " (m/z: " << w.lower_mz << "-" << w.upper_mz << ")" << std::endl;
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
