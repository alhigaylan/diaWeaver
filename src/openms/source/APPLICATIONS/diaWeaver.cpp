// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/diaWeaver.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace OpenMS;

// ----------------------------------------------------------------------
// IMInfo helper - get appropriate array name for the detected unit
// ----------------------------------------------------------------------
const std::string& DiaWeaver::IMInfo::getIMArrayName() const
{
  if (unit == DriftTimeUnit::VSSC)
  {
    return Constants::UserParam::INVERSE_REDUCED_ION_MOBILITY;
  }
  return Constants::UserParam::ION_MOBILITY;
}

// ----------------------------------------------------------------------
// DIAWindow comparator
// ----------------------------------------------------------------------
bool DiaWeaver::DIAWindow::operator<(const DIAWindow& rhs) const
{
  return std::tie(lower_mz, upper_mz, center_mz, lower_im, upper_im) <
         std::tie(rhs.lower_mz, rhs.upper_mz, rhs.center_mz,
                  rhs.lower_im, rhs.upper_im);
}

// ----------------------------------------------------------------------
// Determine DIA windows from MS2 spectra
// ----------------------------------------------------------------------
void DiaWeaver::determineWindows(
  const MSExperiment& raw,
  WindowMap& window_map)
{
  window_map.clear();

  // Store known windows for tolerance-based matching
  std::vector<DIAWindow> known_windows;

  for (Size i = 0; i < raw.size(); ++i)
  {
    const MSSpectrum& spec = raw[i];
    if (spec.getMSLevel() != 2) continue;

    // Validate precursor information
    if (spec.getPrecursors().empty())
    {
      OPENMS_LOG_WARN << "MS2 spectrum at index " << i
                      << " has no precursor information, skipping." << std::endl;
      continue;
    }

    const Precursor& p = spec.getPrecursors()[0];

    // Validate isolation window offsets
    if (p.getIsolationWindowLowerOffset() == 0 && p.getIsolationWindowUpperOffset() == 0)
    {
      OPENMS_LOG_WARN << "MS2 spectrum at index " << i
                      << " has zero isolation window offsets." << std::endl;
    }

    // Build window from precursor info
    DIAWindow candidate;
    candidate.center_mz = p.getMZ();
    candidate.lower_mz = p.getMZ() - p.getIsolationWindowLowerOffset();
    candidate.upper_mz = p.getMZ() + p.getIsolationWindowUpperOffset();

    // Check for ion mobility metadata (use sentinel if not present)
    if (spec.metaValueExists("ion mobility lower limit"))
    {
      candidate.lower_im = spec.getMetaValue("ion mobility lower limit");
    }
    if (spec.metaValueExists("ion mobility upper limit"))
    {
      candidate.upper_im = spec.getMetaValue("ion mobility upper limit");
    }

    // Find existing window using tolerance-based matching
    bool found = false;
    for (auto& known : known_windows)
    {
      if (candidate.isEqual(known))
      {
        window_map[known].push_back(i);
        found = true;
        break;
      }
    }

    // If no matching window found, add as new window
    if (!found)
    {
      known_windows.push_back(candidate);
      window_map[candidate].push_back(i);
    }
  }

  OPENMS_LOG_INFO << "Determined " << window_map.size()
                  << " unique DIA windows from " << raw.size() << " spectra." << std::endl;
}

// ----------------------------------------------------------------------
// Determine ion mobility info for the dataset
// ----------------------------------------------------------------------
DiaWeaver::IMInfo DiaWeaver::determineIMInfo(
  const MSExperiment& raw,
  const WindowMap& window_map)
{
  IMInfo info;

  // Check 1: Do any windows have ion mobility bounds?
  bool windows_have_im = false;
  for (const auto& it : window_map)
  {
    if (it.first.hasIonMobility())
    {
      windows_have_im = true;
      break;
    }
  }

  if (!windows_have_im)
  {
    OPENMS_LOG_INFO << "DIA windows do not have ion mobility bounds. "
                    << "Filtering by m/z only." << std::endl;
    return info; // available = false
  }

  // Check 2: Do spectra contain ion mobility data arrays?
  // Check first MS1 spectrum
  bool ms1_has_im = false;
  Size ms1_im_index = 0;
  DriftTimeUnit im_unit = DriftTimeUnit::NONE;
  for (const MSSpectrum& spec : raw)
  {
    if (spec.getMSLevel() != 1) continue;
    if (spec.containsIMData())
    {
      const auto [idx, unit] = spec.getIMData();
      ms1_im_index = idx;
      im_unit = unit;
      ms1_has_im = true;
    }
    break;
  }

  // Check first MS2 spectrum
  bool ms2_has_im = false;
  Size ms2_im_index = 0;
  for (const MSSpectrum& spec : raw)
  {
    if (spec.getMSLevel() != 2) continue;
    if (spec.containsIMData())
    {
      const auto [idx, unit] = spec.getIMData();
      ms2_im_index = idx;
      // Use MS2 unit if MS1 didn't have IM
      if (im_unit == DriftTimeUnit::NONE)
      {
        im_unit = unit;
      }
      ms2_has_im = true;
    }
    break;
  }

  // All conditions must be met for IM filtering
  if (!ms1_has_im || !ms2_has_im)
  {
    OPENMS_LOG_WARN << "DIA windows have ion mobility bounds but spectra lack ion mobility data. "
                    << "MS1 has IM: " << (ms1_has_im ? "yes" : "no")
                    << ", MS2 has IM: " << (ms2_has_im ? "yes" : "no")
                    << ". Filtering by m/z only." << std::endl;
    return info;
  }

  // Log if IM array indices differ between MS1 and MS2
  if (ms1_im_index != ms2_im_index)
  {
    OPENMS_LOG_WARN << "Ion mobility array index differs between MS1 (" << ms1_im_index
                    << ") and MS2 (" << ms2_im_index << ")." << std::endl;
  }

  info.available = true;
  info.ms1_im_index = ms1_im_index;
  info.ms2_im_index = ms2_im_index;
  info.unit = im_unit;

  OPENMS_LOG_INFO << "Ion mobility filtering enabled. MS1 IM index: "
                  << info.ms1_im_index << ", MS2 IM index: "
                  << info.ms2_im_index << std::endl;

  return info;
}

// ----------------------------------------------------------------------
// Extract MS2 windows
// ----------------------------------------------------------------------
void DiaWeaver::extractMS2Windows(
  const MSExperiment& raw,
  const WindowMap& window_map,
  DiaWeaver::WindowedExperiments& out_ms2,
  DiaWeaver::WindowedExperiments* out_precursors)
{
  out_ms2.clear();

  const bool save_precursors = (out_precursors != nullptr);
  if (save_precursors)
  {
    out_precursors->clear();
  }

  // Determine ion mobility info for proper handling of float data arrays
  const IMInfo im_info = determineIMInfo(raw, window_map);

  for (const auto& it : window_map)
  {
    const DIAWindow& window = it.first;
    const std::vector<Size>& scan_indices = it.second;

    MSExperiment fragment_exp;
    MSExperiment precursor_exp;

    for (Size idx : scan_indices)
    {
      const MSSpectrum& orig_spec = raw[idx];

      // Fragment spectrum (peaks outside precursor window)
      MSSpectrum frag_spec;
      frag_spec.setRT(orig_spec.getRT());
      frag_spec.setMSLevel(1);  // required by downstream tools

      // Precursor spectrum (peaks inside precursor window)
      MSSpectrum prec_spec;
      if (save_precursors)
      {
        prec_spec.setRT(orig_spec.getRT());
        prec_spec.setMSLevel(1);
      }

      // Prepare IM arrays if available
      MSSpectrum::FloatDataArray frag_im_fda;
      MSSpectrum::FloatDataArray prec_im_fda;
      const MSSpectrum::FloatDataArray* im_array = nullptr;

      if (im_info.available)
      {
        frag_im_fda.setName(im_info.getIMArrayName());
        im_array = &orig_spec.getFloatDataArrays()[im_info.ms2_im_index];

        if (save_precursors)
        {
          prec_im_fda.setName(im_info.getIMArrayName());
        }
      }

      // Separate peaks into fragments and precursors
      for (Size i = 0; i < orig_spec.size(); ++i)
      {
        const double mz = orig_spec[i].getMZ();
        const bool in_precursor_window = (mz >= window.lower_mz && mz <= window.upper_mz);

        if (in_precursor_window)
        {
          // Peak is within precursor isolation window (unfragmented precursor)
          if (save_precursors)
          {
            prec_spec.push_back(orig_spec[i]);
            if (im_array)
            {
              prec_im_fda.push_back((*im_array)[i]);
            }
          }
        }
        else
        {
          // Peak is outside precursor window (fragment ion)
          frag_spec.push_back(orig_spec[i]);
          if (im_array)
          {
            frag_im_fda.push_back((*im_array)[i]);
          }
        }
      }
      
      // Add fragment spectrum
      if (!frag_spec.empty())
      {
        if (im_array)
        {
          frag_spec.getFloatDataArrays().push_back(std::move(frag_im_fda));
        }
        frag_spec.sortByPosition();
        fragment_exp.addSpectrum(frag_spec);
      }

      // Add precursor spectrum
      if (save_precursors && !prec_spec.empty())
      {
        if (im_array)
        {
          prec_spec.getFloatDataArrays().push_back(std::move(prec_im_fda));
        }
        prec_spec.sortByPosition();
        precursor_exp.addSpectrum(prec_spec);
      }
    }

    fragment_exp.sortSpectra();
    out_ms2.emplace(window, std::move(fragment_exp));

    if (save_precursors && !precursor_exp.empty())
    {
      precursor_exp.sortSpectra();
      out_precursors->emplace(window, std::move(precursor_exp));
    }
  }
}

// ----------------------------------------------------------------------
// Extract MS1 windows
// ----------------------------------------------------------------------
void DiaWeaver::extractMS1Windows(
  const MSExperiment& raw,
  const WindowMap& window_map,
  DiaWeaver::WindowedExperiments& out_ms1)
{
  out_ms1.clear();

  // Determine ion mobility availability once for the entire dataset
  const IMInfo im_info = determineIMInfo(raw, window_map);

  for (const auto& it : window_map)
  {
    const DIAWindow& window = it.first;
    MSExperiment exp;

    for (const MSSpectrum& spec : raw)
    {
      if (spec.getMSLevel() != 1) continue;

      MSSpectrum new_spec;
      new_spec.setRT(spec.getRT());

      MSSpectrum::FloatDataArray im_fda;
      if (im_info.available)
      {
        im_fda.setName(im_info.getIMArrayName());
      }

      const MSSpectrum::FloatDataArray* im_array = nullptr;
      if (im_info.available)
      {
        im_array = &spec.getFloatDataArrays()[im_info.ms1_im_index];
      }

      for (Size i = 0; i < spec.size(); ++i)
      {
        const double mz = spec[i].getMZ();

        // Always filter by m/z
        if (mz < window.lower_mz || mz > window.upper_mz)
        {
          continue;
        }

        // Optionally filter by ion mobility
        if (im_array)
        {
          const double im = (*im_array)[i];
          if (im < window.lower_im || im > window.upper_im)
          {
            continue;
          }
          im_fda.push_back(im);
        }

        new_spec.push_back(spec[i]);
      }

      if (new_spec.empty()) continue;

      if (im_array)
      {
        new_spec.getFloatDataArrays().push_back(std::move(im_fda));
      }
      new_spec.sortByPosition();

      exp.addSpectrum(new_spec);
    }

    if (!exp.empty())
    {
      exp.sortSpectra();
      out_ms1.emplace(window, std::move(exp));
    }
  }
}

// ========================================================================
// OnDiscMSExperiment implementations for memory-efficient processing

// Determine DIA windows from MS2 spectra (on-disk version)
// Uses metadata only - no peak data loading required
// ========================================================================
void DiaWeaver::determineWindows(
  OnDiscMSExperiment& raw,
  WindowMap& window_map)
{
  window_map.clear();

  // Get metadata - this is loaded once and cached
  auto meta = raw.getMetaData();
  if (!meta)
  {
    OPENMS_LOG_ERROR << "Failed to load metadata from OnDiscMSExperiment." << std::endl;
    return;
  }

  // Store known windows for tolerance-based matching
  std::vector<DIAWindow> known_windows;

  for (Size i = 0; i < meta->size(); ++i)
  {
    const MSSpectrum& spec = (*meta)[i];
    if (spec.getMSLevel() != 2) continue;

    // Validate precursor information
    if (spec.getPrecursors().empty())
    {
      OPENMS_LOG_WARN << "MS2 spectrum at index " << i
                      << " has no precursor information, skipping." << std::endl;
      continue;
    }

    const Precursor& p = spec.getPrecursors()[0];

    // Validate isolation window offsets
    /// TODO: is this acceptable? is default MSSpectrum() set to 0 offset if not DIA data?
    if (p.getIsolationWindowLowerOffset() == 0 && p.getIsolationWindowUpperOffset() == 0)
    {
      OPENMS_LOG_WARN << "MS2 spectrum at index " << i
                      << " has zero isolation window offsets." << std::endl;
    }

    // Build window from precursor info
    DIAWindow candidate;
    candidate.center_mz = p.getMZ();
    candidate.lower_mz = p.getMZ() - p.getIsolationWindowLowerOffset();
    candidate.upper_mz = p.getMZ() + p.getIsolationWindowUpperOffset();

    // Check for ion mobility metadata (use sentinel if not present)
    if (spec.metaValueExists("ion mobility lower limit"))
    {
      candidate.lower_im = spec.getMetaValue("ion mobility lower limit");
    }
    if (spec.metaValueExists("ion mobility upper limit"))
    {
      candidate.upper_im = spec.getMetaValue("ion mobility upper limit");
    }

    // Find existing window using tolerance-based matching
    bool found = false;
    for (auto& known : known_windows)
    {
      if (candidate.isEqual(known))
      {
        window_map[known].push_back(i);
        found = true;
        break;
      }
    }

    // If no matching window found, add as new window
    if (!found)
    {
      known_windows.push_back(candidate);
      window_map[candidate].push_back(i);
    }
  }

  OPENMS_LOG_INFO << "Determined " << window_map.size()
                  << " unique DIA windows from " << raw.getNrSpectra() << " spectra." << std::endl;
}

// ----------------------------------------------------------------------
// Determine ion mobility info for the dataset
// Loads only the first MS1 and MS2 spectra to check for IM data
// ----------------------------------------------------------------------
DiaWeaver::IMInfo DiaWeaver::determineIMInfo(
  OnDiscMSExperiment& raw,
  const WindowMap& window_map)
{
  IMInfo info;

  // Check 1: Do any windows have ion mobility bounds?
  bool windows_have_im = false;
  for (const auto& it : window_map)
  {
    if (it.first.hasIonMobility())
    {
      windows_have_im = true;
      break;
    }
  }

  if (!windows_have_im)
  {
    OPENMS_LOG_INFO << "DIA windows do not have ion mobility bounds. "
                    << "Filtering by m/z only." << std::endl;
    return info; // available = false
  }

  // Get metadata for MS level info
  auto meta = raw.getMetaData();
  if (!meta)
  {
    OPENMS_LOG_WARN << "Failed to load metadata. Disabling IM filtering." << std::endl;
    return info;
  }

  // Check 2: Do spectra contain ion mobility data arrays?
  // Check first MS1 spectrum by loading it from disk
  bool ms1_has_im = false;
  Size ms1_im_index = 0;
  DriftTimeUnit im_unit = DriftTimeUnit::NONE;
  for (Size i = 0; i < meta->size(); ++i)
  {
    if ((*meta)[i].getMSLevel() != 1) continue;
    // Load the actual spectrum to check for IM data
    MSSpectrum spec = raw.getSpectrum(i);
    if (spec.containsIMData())
    {
      const auto [idx, unit] = spec.getIMData();
      ms1_im_index = idx;
      im_unit = unit;
      ms1_has_im = true;
    }
    break;
  }

  // Check first MS2 spectrum
  bool ms2_has_im = false;
  Size ms2_im_index = 0;
  for (Size i = 0; i < meta->size(); ++i)
  {
    if ((*meta)[i].getMSLevel() != 2) continue;
    // Load the actual spectrum to check for IM data
    MSSpectrum spec = raw.getSpectrum(i);
    if (spec.containsIMData())
    {
      const auto [idx, unit] = spec.getIMData();
      ms2_im_index = idx;
      // Use MS2 unit if MS1 didn't have IM
      if (im_unit == DriftTimeUnit::NONE)
      {
        im_unit = unit;
      }
      ms2_has_im = true;
    }
    break;
  }

  // All conditions must be met for IM filtering
  if (!ms1_has_im || !ms2_has_im)
  {
    OPENMS_LOG_WARN << "DIA windows have ion mobility bounds but spectra lack ion mobility data. "
                    << "MS1 has IM: " << (ms1_has_im ? "yes" : "no")
                    << ", MS2 has IM: " << (ms2_has_im ? "yes" : "no")
                    << ". Filtering by m/z only." << std::endl;
    return info;
  }

  // Log if IM array indices differ between MS1 and MS2
  if (ms1_im_index != ms2_im_index)
  {
    OPENMS_LOG_WARN << "Ion mobility array index differs between MS1 (" << ms1_im_index
                    << ") and MS2 (" << ms2_im_index << ")." << std::endl;
  }

  info.available = true;
  info.ms1_im_index = ms1_im_index;
  info.ms2_im_index = ms2_im_index;
  info.unit = im_unit;

  OPENMS_LOG_INFO << "Ion mobility filtering enabled. MS1 IM index: "
                  << info.ms1_im_index << ", MS2 IM index: "
                  << info.ms2_im_index << std::endl;

  return info;
}

// ----------------------------------------------------------------------
// Extract MS2 windows
// Loads spectra on-demand from disk
// ----------------------------------------------------------------------
void DiaWeaver::extractMS2Windows(
  OnDiscMSExperiment& raw,
  const WindowMap& window_map,
  DiaWeaver::WindowedExperiments& out_ms2,
  DiaWeaver::WindowedExperiments* out_precursors)
{
  out_ms2.clear();

  const bool save_precursors = (out_precursors != nullptr);
  if (save_precursors)
  {
    out_precursors->clear();
  }

  // Determine ion mobility info for proper handling of float data arrays
  const IMInfo im_info = determineIMInfo(raw, window_map);

  for (const auto& it : window_map)
  {
    const DIAWindow& window = it.first;
    const std::vector<Size>& scan_indices = it.second;

    MSExperiment fragment_exp;
    MSExperiment precursor_exp;

    for (Size idx : scan_indices)
    {
      // Load spectrum from disk on-demand
      MSSpectrum orig_spec = raw.getSpectrum(idx);

      // Fragment spectrum (peaks outside precursor window)
      MSSpectrum frag_spec;
      frag_spec.setRT(orig_spec.getRT());
      frag_spec.setMSLevel(1);  // required by downstream tools

      // Precursor spectrum (peaks inside precursor window)
      MSSpectrum prec_spec;
      if (save_precursors)
      {
        prec_spec.setRT(orig_spec.getRT());
        prec_spec.setMSLevel(1);
      }

      // Prepare IM arrays if available
      MSSpectrum::FloatDataArray frag_im_fda;
      MSSpectrum::FloatDataArray prec_im_fda;
      const MSSpectrum::FloatDataArray* im_array = nullptr;

      if (im_info.available && orig_spec.getFloatDataArrays().size() > im_info.ms2_im_index)
      {
        frag_im_fda.setName(im_info.getIMArrayName());
        im_array = &orig_spec.getFloatDataArrays()[im_info.ms2_im_index];

        if (save_precursors)
        {
          prec_im_fda.setName(im_info.getIMArrayName());
        }
      }

      // Separate peaks into fragments and precursors
      for (Size i = 0; i < orig_spec.size(); ++i)
      {
        const double mz = orig_spec[i].getMZ();
        const bool in_precursor_window = (mz >= window.lower_mz && mz <= window.upper_mz);

        if (in_precursor_window)
        {
          // Peak is within precursor isolation window (unfragmented precursors)
          if (save_precursors)
          {
            prec_spec.push_back(orig_spec[i]);
            if (im_array)
            {
              prec_im_fda.push_back((*im_array)[i]);
            }
          }
        }
        else
        {
          // Peak is outside precursor window (fragment ions)
          frag_spec.push_back(orig_spec[i]);
          if (im_array)
          {
            frag_im_fda.push_back((*im_array)[i]);
          }
        }
      }

      // Add fragment spectrum
      if (!frag_spec.empty())
      {
        if (im_array)
        {
          frag_spec.getFloatDataArrays().push_back(std::move(frag_im_fda));
        }
        frag_spec.sortByPosition();
        fragment_exp.addSpectrum(frag_spec);
      }

      // Add precursor spectrum
      if (save_precursors && !prec_spec.empty())
      {
        if (im_array)
        {
          prec_spec.getFloatDataArrays().push_back(std::move(prec_im_fda));
        }
        prec_spec.sortByPosition();
        precursor_exp.addSpectrum(prec_spec);
      }
    }

    fragment_exp.sortSpectra();
    out_ms2.emplace(window, std::move(fragment_exp));

    if (save_precursors && !precursor_exp.empty())
    {
      precursor_exp.sortSpectra();
      out_precursors->emplace(window, std::move(precursor_exp));
    }
  }
}

// ----------------------------------------------------------------------
// Extract MS1 windows
// Loads MS1 spectra on-demand from disk
// ----------------------------------------------------------------------
void DiaWeaver::extractMS1Windows(
  OnDiscMSExperiment& raw,
  const WindowMap& window_map,
  DiaWeaver::WindowedExperiments& out_ms1)
{
  out_ms1.clear();

  // Determine ion mobility availability once for the entire dataset
  const IMInfo im_info = determineIMInfo(raw, window_map);

  // Get metadata for MS level checking
  auto meta = raw.getMetaData();
  if (!meta)
  {
    OPENMS_LOG_ERROR << "Failed to load metadata from OnDiscMSExperiment." << std::endl;
    return;
  }

  for (const auto& it : window_map)
  {
    const DIAWindow& window = it.first;
    MSExperiment exp;

    for (Size i = 0; i < meta->size(); ++i)
    {
      if ((*meta)[i].getMSLevel() != 1) continue;

      // Load MS1 spectrum from disk
      MSSpectrum spec = raw.getSpectrum(i);

      MSSpectrum new_spec;
      new_spec.setRT(spec.getRT());

      MSSpectrum::FloatDataArray im_fda;
      if (im_info.available)
      {
        im_fda.setName(im_info.getIMArrayName());
      }

      const MSSpectrum::FloatDataArray* im_array = nullptr;
      if (im_info.available && spec.getFloatDataArrays().size() > im_info.ms1_im_index)
      {
        im_array = &spec.getFloatDataArrays()[im_info.ms1_im_index];
      }

      for (Size j = 0; j < spec.size(); ++j)
      {
        const double mz = spec[j].getMZ();

        // Always filter by m/z
        if (mz < window.lower_mz || mz > window.upper_mz)
        {
          continue;
        }

        // Optionally filter by ion mobility
        if (im_array)
        {
          const double im = (*im_array)[j];
          if (im < window.lower_im || im > window.upper_im)
          {
            continue;
          }
          im_fda.push_back(im);
        }

        new_spec.push_back(spec[j]);
      }

      if (new_spec.empty()) continue;

      if (im_array)
      {
        new_spec.getFloatDataArrays().push_back(std::move(im_fda));
      }
      new_spec.sortByPosition();

      exp.addSpectrum(new_spec);
    }

    if (!exp.empty())
    {
      exp.sortSpectra();
      out_ms1.emplace(window, std::move(exp));
    }
  }
}

// ========================================================================
// Single-window extraction for incremental/streaming processing
// ========================================================================

// ----------------------------------------------------------------------
// Extract MS2 spectra for a single window
// ----------------------------------------------------------------------
void DiaWeaver::extractSingleMS2Window(
  OnDiscMSExperiment& raw,
  const DIAWindow& window,
  const std::vector<Size>& indices,
  const IMInfo& im_info,
  MSExperiment& out_ms2,
  MSExperiment* out_precursor)
{
  out_ms2.clear(true);
  const bool save_precursors = (out_precursor != nullptr);
  if (save_precursors)
  {
    out_precursor->clear(true);
  }

  for (Size idx : indices)
  {
    // Load spectrum from disk on-demand
    MSSpectrum orig_spec = raw.getSpectrum(idx);

    // Fragment spectrum (peaks outside precursor window)
    MSSpectrum frag_spec;
    frag_spec.setRT(orig_spec.getRT());
    frag_spec.setMSLevel(1);  // required by downstream tools

    // Precursor spectrum (peaks inside precursor window)
    MSSpectrum prec_spec;
    if (save_precursors)
    {
      prec_spec.setRT(orig_spec.getRT());
      prec_spec.setMSLevel(1);
    }

    // Prepare IM arrays if available
    MSSpectrum::FloatDataArray frag_im_fda;
    MSSpectrum::FloatDataArray prec_im_fda;
    const MSSpectrum::FloatDataArray* im_array = nullptr;

    if (im_info.available && orig_spec.getFloatDataArrays().size() > im_info.ms2_im_index)
    {
      frag_im_fda.setName(im_info.getIMArrayName());
      im_array = &orig_spec.getFloatDataArrays()[im_info.ms2_im_index];

      if (save_precursors)
      {
        prec_im_fda.setName(im_info.getIMArrayName());
      }
    }

    // Separate peaks into fragments and precursors
    for (Size i = 0; i < orig_spec.size(); ++i)
    {
      const double mz = orig_spec[i].getMZ();
      const bool in_precursor_window = (mz >= window.lower_mz && mz <= window.upper_mz);

      if (in_precursor_window)
      {
        if (save_precursors)
        {
          prec_spec.push_back(orig_spec[i]);
          if (im_array)
          {
            prec_im_fda.push_back((*im_array)[i]);
          }
        }
      }
      else
      {
        frag_spec.push_back(orig_spec[i]);
        if (im_array)
        {
          frag_im_fda.push_back((*im_array)[i]);
        }
      }
    }

    // Add fragment spectrum
    if (!frag_spec.empty())
    {
      if (im_array)
      {
        frag_spec.getFloatDataArrays().push_back(std::move(frag_im_fda));
      }
      frag_spec.sortByPosition();
      out_ms2.addSpectrum(frag_spec);
    }

    // Add precursor spectrum
    if (save_precursors && !prec_spec.empty())
    {
      if (im_array)
      {
        prec_spec.getFloatDataArrays().push_back(std::move(prec_im_fda));
      }
      prec_spec.sortByPosition();
      out_precursor->addSpectrum(prec_spec);
    }
  }

  out_ms2.sortSpectra();
  if (save_precursors)
  {
    out_precursor->sortSpectra();
  }
}

// ----------------------------------------------------------------------
// Extract MS1 spectra for a single window
// ----------------------------------------------------------------------
void DiaWeaver::extractSingleMS1Window(
  OnDiscMSExperiment& raw,
  const DIAWindow& window,
  const IMInfo& im_info,
  MSExperiment& out_ms1)
{
  out_ms1.clear(true);

  // Get metadata for MS level checking
  auto meta = raw.getMetaData();
  if (!meta)
  {
    OPENMS_LOG_ERROR << "Failed to load metadata from OnDiscMSExperiment." << std::endl;
    return;
  }

  for (Size i = 0; i < meta->size(); ++i)
  {
    if ((*meta)[i].getMSLevel() != 1) continue;

    // Load MS1 spectrum from disk
    MSSpectrum spec = raw.getSpectrum(i);

    MSSpectrum new_spec;
    new_spec.setRT(spec.getRT());

    MSSpectrum::FloatDataArray im_fda;
    if (im_info.available)
    {
      im_fda.setName(im_info.getIMArrayName());
    }

    const MSSpectrum::FloatDataArray* im_array = nullptr;
    if (im_info.available && spec.getFloatDataArrays().size() > im_info.ms1_im_index)
    {
      im_array = &spec.getFloatDataArrays()[im_info.ms1_im_index];
    }

    for (Size k = 0; k < spec.size(); ++k)
    {
      const double mz = spec[k].getMZ();

      // Always filter by m/z
      if (mz < window.lower_mz || mz > window.upper_mz)
      {
        continue;
      }

      // Optionally filter by ion mobility
      if (im_array)
      {
        const double im = (*im_array)[k];
        if (im < window.lower_im || im > window.upper_im)
        {
          continue;
        }
        im_fda.push_back(im);
      }

      new_spec.push_back(spec[k]);
    }

    if (new_spec.empty()) continue;

    if (im_array)
    {
      new_spec.getFloatDataArrays().push_back(std::move(im_fda));
    }
    new_spec.sortByPosition();

    out_ms1.addSpectrum(new_spec);
  }

  out_ms1.sortSpectra();
}
