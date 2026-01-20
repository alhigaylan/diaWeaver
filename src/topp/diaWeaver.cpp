// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/CachedMzML.h>
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
      "Number of threads to use for parallel window processing (default: 1)",
      false);
    setMinInt_("threads", 1);
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

    // ------------------------------
    // Step 2: Create binary cache for fast parallel I/O
    // Load the full file once and cache it for parallel access
    // ------------------------------
    OPENMS_LOG_INFO << "Loading and caching data for fast parallel access..." << std::endl;
    MzMLFile mzml_loader;
    MSExperiment full_exp;
    mzml_loader.load(in, full_exp);

    const String cache_file = out + "/.diaWeaver_cache.mzML";
    CachedmzML::store(cache_file, full_exp);

    // Clear full_exp to free memory - we'll read from cache now
    full_exp.clear(true);

    OPENMS_LOG_INFO << "Processing " << total_windows << " DIA windows";
#ifdef _OPENMP
    OPENMS_LOG_INFO << " using " << num_threads << " thread(s)";
#endif
    OPENMS_LOG_INFO << " with fast binary cache..." << std::endl;

    if (im_info.available)
    {
      OPENMS_LOG_INFO << "Ion mobility data detected. Using PeakPickerIM (mobilogram method)." << std::endl;
    }
    else
    {
      OPENMS_LOG_INFO << "No ion mobility data detected. Using PeakPickerHiRes." << std::endl;
    }

    // ------------------------------
    // Step 3: Process windows in parallel
    // Each thread opens its own CachedmzML for thread-safe fast reading
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

      // Each thread opens its own CachedmzML for thread-safe fast binary access
      CachedmzML thread_cache(cache_file);

      // Thread-local peak picker instances (avoids lock contention)
      PeakPickerIM peak_picker_im;
      PeakPickerHiRes peak_picker_hr;
      if (im_info.available)
      {
        peak_picker_im.setParameters(ppim_params);
      }
      else
      {
        peak_picker_hr.setParameters(pphr_params);
      }

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

      // Extract MS2 (using fast binary cache)
      DiaWeaver::extractSingleMS2Window(thread_cache, w, indices, im_info, ms2_exp,
                                         save_precursors ? &precursor_exp : nullptr);

      // Apply peak picking to MS2 spectra
      for (auto& spec : ms2_exp)
      {
        if (im_info.available)
        {
          peak_picker_im.pickIMTraces(spec);
        }
        else
        {
          MSSpectrum picked;
          peak_picker_hr.pick(spec, picked);
          spec = std::move(picked);
        }
      }

      if (!ms2_exp.empty())
      {
        mzml.store(out + "/ms2_" + fname_base, ms2_exp);
      }

      // Apply peak picking to precursors and write if requested
      if (save_precursors && !precursor_exp.empty())
      {
        for (auto& spec : precursor_exp)
        {
          if (im_info.available)
          {
            peak_picker_im.pickIMTraces(spec);
          }
          else
          {
            MSSpectrum picked;
            peak_picker_hr.pick(spec, picked);
            spec = std::move(picked);
          }
        }
        mzml.store(out + "/precursor_" + fname_base, precursor_exp);
      }

      // Extract MS1 (using fast binary cache)
      DiaWeaver::extractSingleMS1Window(thread_cache, w, im_info, ms1_exp);

      // Apply peak picking to MS1 spectra
      for (auto& spec : ms1_exp)
      {
        if (im_info.available)
        {
          peak_picker_im.pickIMTraces(spec);
        }
        else
        {
          MSSpectrum picked;
          peak_picker_hr.pick(spec, picked);
          spec = std::move(picked);
        }
      }

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

    // Clean up cache file
    File::remove(cache_file);
    File::remove(cache_file + ".cached");

    OPENMS_LOG_INFO << "Finished processing all windows." << std::endl;

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPDiaWeaver tool;
  return tool.main(argc, argv);
}
