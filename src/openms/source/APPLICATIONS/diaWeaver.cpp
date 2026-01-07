// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/diaWeaver.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace OpenMS;

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
// Extract MS2 windows
// ----------------------------------------------------------------------
void DiaWeaver::extractMS2Windows(
  const MSExperiment& raw,
  const WindowMap& window_map,
  DiaWeaver::WindowedExperiments& out_ms2)
{
  out_ms2.clear();

  for (const auto& it : window_map)
  {
    const DIAWindow& window = it.first;
    const std::vector<Size>& scan_indices = it.second;

    MSExperiment exp;

    for (Size idx : scan_indices)
    {
      MSSpectrum spec = raw[idx];  // copy
      spec.setMSLevel(1);          // required by downstream tools
      exp.addSpectrum(spec);
    }

    exp.sortSpectra();
    out_ms2.emplace(window, std::move(exp));
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

  for (const auto& it : window_map)
  {
    const DIAWindow& window = it.first;
    MSExperiment exp;

    // Check if this window has ion mobility filtering
    const bool filter_by_im = window.hasIonMobility();

    for (const MSSpectrum& spec : raw)
    {
      if (spec.getMSLevel() != 1) continue;

      // If filtering by IM, we need the float data array
      const MSSpectrum::FloatDataArray* im_array = nullptr;
      if (filter_by_im)
      {
        if (spec.getFloatDataArrays().empty())
        {
          continue; // Skip spectra without IM data when IM filtering is needed
        }
        im_array = &spec.getFloatDataArrays()[0];
      }

      MSSpectrum new_spec;
      new_spec.setRT(spec.getRT());

      MSSpectrum::FloatDataArray im_fda;
      if (filter_by_im)
      {
        im_fda.setName("Ion Mobility");
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
        if (filter_by_im)
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

      // Only add IM array if we were filtering by IM
      if (filter_by_im)
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
