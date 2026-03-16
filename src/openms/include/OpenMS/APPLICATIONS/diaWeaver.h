// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/CachedMzML.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/OpenMSConfig.h>
#include <cmath>
#include <map>
#include <vector>


namespace OpenMS
{
/**
  @brief An application to generate pseudo MS2 fragmentation spectra from DIA data
*/
  class OPENMS_DLLAPI DiaWeaver
  {
  public:
    /// Sentinel value indicating ion mobility is not set
    static constexpr double IM_NOT_SET = -1.0;

    /// Default tolerance for comparing window boundaries
    static constexpr double WINDOW_TOLERANCE = 1e-6;

    /**
      @brief Data structure representing a DIA isolation window

      Stores m/z boundaries (lower, upper, center) and optional ion mobility
      boundaries. Ion mobility values are set to IM_NOT_SET (-1) when not present.
    */
    struct DIAWindow
    {
      double lower_mz{0.0};
      double upper_mz{0.0};
      double center_mz{0.0};
      double lower_im{IM_NOT_SET};
      double upper_im{IM_NOT_SET};

      /// Comparison operator for use in std::map
      bool operator<(const DIAWindow& rhs) const;

      /// Tolerance-based equality check (recommended for window matching)
      bool isEqual(const DIAWindow& other, double tolerance = WINDOW_TOLERANCE) const
      {
        return (std::fabs(lower_mz - other.lower_mz) < tolerance) &&
               (std::fabs(upper_mz - other.upper_mz) < tolerance) &&
               (std::fabs(center_mz - other.center_mz) < tolerance) &&
               (std::fabs(lower_im - other.lower_im) < tolerance) &&
               (std::fabs(upper_im - other.upper_im) < tolerance);
      }

      /// Check if ion mobility boundaries are set
      bool hasIonMobility() const
      {
        return lower_im > IM_NOT_SET && upper_im > IM_NOT_SET;
      }
    };

    using WindowMap = std::map<DIAWindow, std::vector<Size>>;
    using WindowedExperiments = std::map<DIAWindow, MSExperiment>;

    /// Build a TSV filename for a given DIA window, e.g. "mz500-525_im068-079.tsv".
    /// An optional tag (e.g. "ms1", "precursor") is appended before the extension
    /// to distinguish files when clustering is run more than once per window.
    static String windowTSVFilename(const DIAWindow& w, const String& tag = "")
    {
      String fname = "mz" + String(static_cast<int>(std::round(w.lower_mz)))
                   + "-" + String(static_cast<int>(std::round(w.upper_mz)));
      if (w.hasIonMobility())
      {
        // Encode IM as zero-padded 3-digit integer (e.g. 0.68 -> "068")
        auto fmt = [](double v) -> String {
          int vi = static_cast<int>(std::round(v * 100));
          String s = String(vi);
          while (s.size() < 3) s = "0" + s;
          return s;
        };
        fname += "_im" + fmt(w.lower_im) + "-" + fmt(w.upper_im);
      }
      if (!tag.empty()) fname += "_" + tag;
      return fname + ".tsv";
    }

    /**
      @brief Information about ion mobility data availability in the dataset

      Used to determine once whether ion mobility filtering should be applied
      and which float data array indices contain the IM data for MS1 and MS2 spectra.
    */
    struct IMInfo
    {
      bool available{false};     ///< True if IM filtering should be applied
      Size ms1_im_index{0};      ///< Index of IM data in MS1 FloatDataArrays (valid only if available=true)
      Size ms2_im_index{0};      ///< Index of IM data in MS2 FloatDataArrays (valid only if available=true)
      DriftTimeUnit unit{DriftTimeUnit::NONE};  ///< Unit of ion mobility data (for proper output naming)

      /// Get the appropriate array name for the detected IM unit
      const std::string& getIMArrayName() const;
    };

    /**
      @brief Determine if ion mobility filtering should be used for this dataset

      Checks three conditions:
      1. At least one DIA window has ion mobility bounds
      2. MS1 spectra contain ion mobility data arrays
      3. MS2 spectra contain ion mobility data arrays

      All conditions must be true for ion mobility filtering to be enabled.
      The IM array index is assumed to be consistent across all spectra.

      @param raw Input MSExperiment containing DIA data
      @param window_map Map of DIA windows to check for IM bounds
      @return IMInfo struct with availability flag and array index
    */
    static IMInfo determineIMInfo(
      const MSExperiment& raw,
      const WindowMap& window_map);

    /**
      @brief Determine DIA windows from MS2 spectra

      Analyzes all MS2 spectra in the experiment and groups them by their
      precursor isolation window. Uses tolerance-based matching to handle
      floating-point precision issues.

      @param raw Input MSExperiment containing DIA data
      @param window_map Output map of DIAWindow to spectrum indices
    */
    static void determineWindows(
      const MSExperiment& raw,
      WindowMap& window_map);

    /**
      @brief Extract MS2 spectra for each window

      Creates separate MSExperiment objects for each DIA window containing
      the corresponding MS2 spectra. Peaks within the precursor isolation
      window are always excluded from the output (fragment ions only).

      @param raw Input MSExperiment containing DIA data
      @param window_map Map of windows to spectrum indices (from determineWindows)
      @param out_ms2 Output map of DIAWindow to MSExperiment (fragment ions only)
      @param out_precursors Optional output map for unfragmented precursor peaks
             (peaks within isolation window). If nullptr, precursor peaks are discarded.
    */
    static void extractMS2Windows(
      const MSExperiment& raw,
      const WindowMap& window_map,
      WindowedExperiments& out_ms2,
      WindowedExperiments* out_precursors = nullptr);

    /**
      @brief Extract MS1 spectra for each window

      Filters MS1 spectra by m/z and ion mobility ranges defined by each
      DIA window. Only peaks within the window boundaries are retained.

      @param raw Input MSExperiment containing DIA data
      @param window_map Map of windows (used for window definitions)
      @param out_ms1 Output map of DIAWindow to filtered MS1 MSExperiment
    */
    static void extractMS1Windows(
      const MSExperiment& raw,
      const WindowMap& window_map,
      WindowedExperiments& out_ms1);

    // ========================================================================
    // OnDiscMSExperiment overloads for memory-efficient processing
    // ========================================================================

    /**
      @brief Determine DIA windows from MS2 spectra (on-disk version)

      Uses only metadata from OnDiscMSExperiment without loading peak data.
      This is memory-efficient for large files as only precursor information
      is needed to determine windows.

      @param raw Input OnDiscMSExperiment containing DIA data
      @param window_map Output map of DIAWindow to spectrum indices
    */
    static void determineWindows(
      OnDiscMSExperiment& raw,
      WindowMap& window_map);

    /**
      @brief Determine if ion mobility filtering should be used (on-disk version)

      @param raw Input OnDiscMSExperiment containing DIA data
      @param window_map Map of DIA windows to check for IM bounds
      @return IMInfo struct with availability flag and array index
    */
    static IMInfo determineIMInfo(
      OnDiscMSExperiment& raw,
      const WindowMap& window_map);

    /**
      @brief Extract MS2 spectra for each window (on-disk version)

      Loads spectra on-demand from disk, processing one spectrum at a time
      to minimize memory usage.

      @param raw Input OnDiscMSExperiment containing DIA data
      @param window_map Map of windows to spectrum indices (from determineWindows)
      @param out_ms2 Output map of DIAWindow to MSExperiment (fragment ions only)
      @param out_precursors Optional output map for unfragmented precursor peaks
    */
    static void extractMS2Windows(
      OnDiscMSExperiment& raw,
      const WindowMap& window_map,
      WindowedExperiments& out_ms2,
      WindowedExperiments* out_precursors = nullptr);

    /**
      @brief Extract MS1 spectra for each window (on-disk version)

      Loads MS1 spectra on-demand from disk, processing one at a time.

      @param raw Input OnDiscMSExperiment containing DIA data
      @param window_map Map of windows (used for window definitions)
      @param out_ms1 Output map of DIAWindow to filtered MS1 MSExperiment
    */
    static void extractMS1Windows(
      OnDiscMSExperiment& raw,
      const WindowMap& window_map,
      WindowedExperiments& out_ms1);

    // ========================================================================
    // Single-window extraction for incremental/streaming processing
    // ========================================================================

    /**
      @brief Extract MS2 spectra for a single window (on-disk version)

      Processes only the specified window, allowing incremental processing
      where each window can be written to disk immediately after extraction.

      @param raw Input OnDiscMSExperiment containing DIA data
      @param window The DIA window to extract
      @param indices Spectrum indices belonging to this window
      @param im_info Pre-computed ion mobility info (from determineIMInfo)
      @param out_ms2 Output MSExperiment for fragment ions
      @param out_precursor Optional output MSExperiment for unfragmented precursors
    */
    static void extractSingleMS2Window(
      OnDiscMSExperiment& raw,
      const DIAWindow& window,
      const std::vector<Size>& indices,
      const IMInfo& im_info,
      MSExperiment& out_ms2,
      MSExperiment* out_precursor = nullptr);

    /**
      @brief Extract MS1 spectra for a single window (on-disk version)

      Processes only the specified window, allowing incremental processing
      where each window can be written to disk immediately after extraction.

      @param raw Input OnDiscMSExperiment containing DIA data
      @param window The DIA window to extract
      @param im_info Pre-computed ion mobility info (from determineIMInfo)
      @param out_ms1 Output MSExperiment for filtered MS1 spectra
    */
    static void extractSingleMS1Window(
      OnDiscMSExperiment& raw,
      const DIAWindow& window,
      const IMInfo& im_info,
      MSExperiment& out_ms1);

    // ========================================================================
    // CachedmzML overloads for fast binary I/O with parallel processing
    // Each thread should create its own CachedmzML instance for thread safety
    // ========================================================================

    /**
      @brief Determine if ion mobility filtering should be used (cached version)

      @param cache Input CachedmzML containing DIA data
      @param window_map Map of DIA windows to check for IM bounds
      @return IMInfo struct with availability flag and array index
    */
    static IMInfo determineIMInfo(
      CachedmzML& cache,
      const WindowMap& window_map);

    /**
      @brief Extract MS2 spectra for a single window (cached version)

      Fast binary I/O version. Each thread should have its own CachedmzML instance.

      @param cache Input CachedmzML containing DIA data
      @param window The DIA window to extract
      @param indices Spectrum indices belonging to this window
      @param im_info Pre-computed ion mobility info (from determineIMInfo)
      @param out_ms2 Output MSExperiment for fragment ions
      @param out_precursor Optional output MSExperiment for unfragmented precursors
    */
    static void extractSingleMS2Window(
      CachedmzML& cache,
      const DIAWindow& window,
      const std::vector<Size>& indices,
      const IMInfo& im_info,
      MSExperiment& out_ms2,
      MSExperiment* out_precursor = nullptr);

    /**
      @brief Extract MS1 spectra for a single window (cached version)

      Fast binary I/O version. Each thread should have its own CachedmzML instance.

      @param cache Input CachedmzML containing DIA data
      @param window The DIA window to extract
      @param im_info Pre-computed ion mobility info (from determineIMInfo)
      @param out_ms1 Output MSExperiment for filtered MS1 spectra
    */
    static void extractSingleMS1Window(
      CachedmzML& cache,
      const DIAWindow& window,
      const IMInfo& im_info,
      MSExperiment& out_ms1);

  };
} // namespace OpenMS
