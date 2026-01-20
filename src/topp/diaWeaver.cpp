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
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>

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

    registerSubsection_("PeakPickerIM", "Parameters for ion mobility peak picking (used when input has IM data)");

    registerSubsection_("PeakPickerHiRes", "Parameters for high-resolution peak picking (used when input has no IM data)");

    registerIntOption_(
      "threads",
      "<n>",
      1,
      "Total number of threads to use for processing",
      false);
    setMinInt_("threads", 1);

    registerIntOption_(
      "threads_outer_loop",
      "<n>",
      -1,
      "Number of threads for the outer loop (over DIA windows). Remaining threads are used for "
      "inner loop (peak picking within each window). Set to -1 to use all threads in outer loop only (no nested parallelism). "
      "Example: with 24 total threads and 4 outer threads, each window gets 6 threads for peak picking.",
      false);
  }

  Param getSubsectionDefaults_(const String& name) const override
  {
    if (name == "PeakPickerIM")
    {
      PeakPickerIM ppim;
      return ppim.getDefaults();
    }
    if (name == "PeakPickerHiRes")
    {
      PeakPickerHiRes pphr;
      return pphr.getDefaults();
    }
    return Param();
  }

  ExitCodes main_(int, const char**) override
  {
    const String in = getStringOption_("in");
    const bool save_precursors = getFlag_("save_unfragmented_precursors");
    const Param ppim_params = getParam_().copy("PeakPickerIM:", true);
    const Param pphr_params = getParam_().copy("PeakPickerHiRes:", true);
#ifdef _OPENMP
    const int num_threads = getIntOption_("threads");
    const int threads_outer_loop = getIntOption_("threads_outer_loop");
#endif

    String out = getStringOption_("out");
    if (out.empty())
    {
      out = File::path(in) + "/" +
            File::basename(in) + "_diaWindows";
    }
    File::makeDir(out);

    // ------------------------------
    // Step 1: Use OnDiscMSExperiment for memory-efficient metadata access
    // This determines DIA windows and IM info without loading all peak data
    // ------------------------------
    OPENMS_LOG_INFO << "Opening file for metadata access..." << std::endl;
    OnDiscMSExperiment on_disc;
    if (!on_disc.openFile(in))
    {
      OPENMS_LOG_ERROR << "Failed to open file as indexed mzML." << std::endl;
      return INPUT_FILE_NOT_FOUND;
    }

    // Determine DIA windows from MS2 metadata (efficient - no peak data loaded)
    DiaWeaver::WindowMap windows;
    DiaWeaver::determineWindows(on_disc, windows);

    // Determine IM info (loads only representative spectra)
    DiaWeaver::IMInfo im_info = DiaWeaver::determineIMInfo(on_disc, windows);

    // Convert map to vector for OpenMP indexed access
    std::vector<std::pair<DiaWeaver::DIAWindow, std::vector<Size>>> window_vec(
      windows.begin(), windows.end());
    const Size total_windows = window_vec.size();

    // Note: OnDiscMSExperiment provides memory-efficient on-demand spectrum loading
    // Each thread will get its own copy via firstprivate (creates separate file handles)

    OPENMS_LOG_INFO << "Processing " << total_windows << " DIA windows";
#ifdef _OPENMP
    OPENMS_LOG_INFO << " using " << num_threads << " thread(s)";
#endif
    OPENMS_LOG_INFO << " with on-disc random access..." << std::endl;

    if (im_info.available)
    {
      OPENMS_LOG_INFO << "Ion mobility data detected. Using PeakPickerIM (mobilogram method)." << std::endl;
    }
    else
    {
      OPENMS_LOG_INFO << "No ion mobility data detected. Using PeakPickerHiRes." << std::endl;
    }

    // ------------------------------
    // Step 3: Process windows in parallel with nested parallelism
    // Outer loop: over DIA windows
    // Inner loop: over spectra within each window (for peak picking)
    // ------------------------------
    Size processed = 0;

#ifdef _OPENMP
    // Store total number of threads available
    const int total_threads = num_threads;

    // Calculate outer and inner thread counts
    int outer_threads = total_threads;
    int inner_threads = 1;

    if (threads_outer_loop > 0)
    {
      // User specified nested parallelism
      outer_threads = std::min(threads_outer_loop, total_threads);
      inner_threads = std::max(1, total_threads / outer_threads);
      omp_set_nested(1);
      omp_set_dynamic(0);
      OPENMS_LOG_INFO << "Using nested parallelism: " << outer_threads << " outer threads x "
                      << inner_threads << " inner threads for peak picking." << std::endl;
    }
    else
    {
      OPENMS_LOG_INFO << "Using " << outer_threads << " threads for window processing (no nested parallelism)." << std::endl;
    }

    omp_set_num_threads(outer_threads);
#pragma omp parallel for schedule(dynamic, 1) firstprivate(on_disc)
#endif
    for (SignedSize idx = 0; idx < static_cast<SignedSize>(total_windows); ++idx)
    {
      const DiaWeaver::DIAWindow& w = window_vec[idx].first;
      const std::vector<Size>& indices = window_vec[idx].second;

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

      // Extract MS2 (on-demand from disk - each thread has its own file handle)
      DiaWeaver::extractSingleMS2Window(on_disc, w, indices, im_info, ms2_exp,
                                         save_precursors ? &precursor_exp : nullptr);

      // Apply peak picking to MS2 spectra (inner parallel loop)
      // Use parallel region to create thread-private pickers (more efficient than per-iteration)
#ifdef _OPENMP
#pragma omp parallel num_threads(inner_threads)
      {
        PeakPickerIM picker_im;
        PeakPickerHiRes picker_hr;
        if (im_info.available)
        {
          picker_im.setParameters(ppim_params);
        }
        else
        {
          picker_hr.setParameters(pphr_params);
        }

#pragma omp for schedule(dynamic, 1)
        for (SignedSize s = 0; s < static_cast<SignedSize>(ms2_exp.size()); ++s)
        {
          if (im_info.available)
          {
            picker_im.pickIMTraces(ms2_exp[s]);
          }
          else
          {
            MSSpectrum picked;
            picker_hr.pick(ms2_exp[s], picked);
            ms2_exp[s] = std::move(picked);
          }
        }
      }
#else
      for (SignedSize s = 0; s < static_cast<SignedSize>(ms2_exp.size()); ++s)
      {
        PeakPickerIM picker_im;
        PeakPickerHiRes picker_hr;
        if (im_info.available)
        {
          picker_im.setParameters(ppim_params);
          picker_im.pickIMTraces(ms2_exp[s]);
        }
        else
        {
          picker_hr.setParameters(pphr_params);
          MSSpectrum picked;
          picker_hr.pick(ms2_exp[s], picked);
          ms2_exp[s] = std::move(picked);
        }
      }
#endif

      if (!ms2_exp.empty())
      {
        mzml.store(out + "/ms2_" + fname_base, ms2_exp);
      }

      // Apply peak picking to precursors (inner parallel loop)
      if (save_precursors && !precursor_exp.empty())
      {
#ifdef _OPENMP
#pragma omp parallel num_threads(inner_threads)
        {
          PeakPickerIM picker_im;
          PeakPickerHiRes picker_hr;
          if (im_info.available)
          {
            picker_im.setParameters(ppim_params);
          }
          else
          {
            picker_hr.setParameters(pphr_params);
          }

#pragma omp for schedule(dynamic, 1)
          for (SignedSize s = 0; s < static_cast<SignedSize>(precursor_exp.size()); ++s)
          {
            if (im_info.available)
            {
              picker_im.pickIMTraces(precursor_exp[s]);
            }
            else
            {
              MSSpectrum picked;
              picker_hr.pick(precursor_exp[s], picked);
              precursor_exp[s] = std::move(picked);
            }
          }
        }
#else
        for (SignedSize s = 0; s < static_cast<SignedSize>(precursor_exp.size()); ++s)
        {
          PeakPickerIM picker_im;
          PeakPickerHiRes picker_hr;
          if (im_info.available)
          {
            picker_im.setParameters(ppim_params);
            picker_im.pickIMTraces(precursor_exp[s]);
          }
          else
          {
            picker_hr.setParameters(pphr_params);
            MSSpectrum picked;
            picker_hr.pick(precursor_exp[s], picked);
            precursor_exp[s] = std::move(picked);
          }
        }
#endif
        mzml.store(out + "/precursor_" + fname_base, precursor_exp);
      }

      // Extract MS1 (on-demand from disk - each thread has its own file handle)
      DiaWeaver::extractSingleMS1Window(on_disc, w, im_info, ms1_exp);

      // Apply peak picking to MS1 spectra (inner parallel loop)
#ifdef _OPENMP
#pragma omp parallel num_threads(inner_threads)
      {
        PeakPickerIM picker_im;
        PeakPickerHiRes picker_hr;
        if (im_info.available)
        {
          picker_im.setParameters(ppim_params);
        }
        else
        {
          picker_hr.setParameters(pphr_params);
        }

#pragma omp for schedule(dynamic, 1)
        for (SignedSize s = 0; s < static_cast<SignedSize>(ms1_exp.size()); ++s)
        {
          if (im_info.available)
          {
            picker_im.pickIMTraces(ms1_exp[s]);
          }
          else
          {
            MSSpectrum picked;
            picker_hr.pick(ms1_exp[s], picked);
            ms1_exp[s] = std::move(picked);
          }
        }
      }
#else
      for (SignedSize s = 0; s < static_cast<SignedSize>(ms1_exp.size()); ++s)
      {
        PeakPickerIM picker_im;
        PeakPickerHiRes picker_hr;
        if (im_info.available)
        {
          picker_im.setParameters(ppim_params);
          picker_im.pickIMTraces(ms1_exp[s]);
        }
        else
        {
          picker_hr.setParameters(pphr_params);
          MSSpectrum picked;
          picker_hr.pick(ms1_exp[s], picked);
          ms1_exp[s] = std::move(picked);
        }
      }
#endif

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

#ifdef _OPENMP
    // Restore total thread count if nested parallelism was used
    if (threads_outer_loop > 0)
    {
      omp_set_num_threads(total_threads);
    }
#endif

    OPENMS_LOG_INFO << "Finished processing all windows." << std::endl;

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPDiaWeaver tool;
  return tool.main(argc, argv);
}
