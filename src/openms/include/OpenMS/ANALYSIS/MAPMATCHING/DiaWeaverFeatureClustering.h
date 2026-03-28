// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConversionHelper.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  /**
    @brief Feature map alignment for diaWeaver pseudo spectra, based on pose clustering.

    This is a copy of MapAlignmentAlgorithmPoseClustering, dedicated to aligning
    diaWeaver FeatureMaps (built from pseudo MS2 spectra precursor information).
    Keeping it separate allows diaWeaver-specific modifications without affecting
    the upstream algorithm.

    Pose clustering analyzes pair distances to find the most probable affine
    transformation of retention times between a scene map and a reference map.

    The algorithm selects the top N most intense features per map (controlled by
    @p max_num_peaks_considered). Set to -1 to use all features.

    @ingroup MapAlignment
  */
  class OPENMS_DLLAPI DiaWeaverFeatureClustering :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Default constructor
    DiaWeaverFeatureClustering();

    /// Destructor
    ~DiaWeaverFeatureClustering() override;

    void align(const FeatureMap& map, TransformationDescription& trafo);
    void align(const PeakMap& map, TransformationDescription& trafo);
    void align(const ConsensusMap& map, TransformationDescription& trafo);

    /**
      @brief Load a diaWeaver pseudo-spectrum mzML file as a FeatureMap.

      Each MS2 spectrum in the file becomes one Feature:
        - RT        ← spectrum retention time
        - m/z       ← precursor m/z
        - charge    ← precursor charge state
        - intensity ← precursor intensity assigned by diaWeaver,
                      used by the superimposer to weight alignment votes

      Spectra with no precursor or zero precursor m/z are skipped.

      @param mzml_file  Path to a diaWeaver pseudo-spectrum mzML file.
      @param feature_map  Output FeatureMap, one feature per MS2 spectrum.
    */
    static void loadMzMLAsFeatureMap(const String& mzml_file, FeatureMap& feature_map);

    /// Sets the reference map for alignment
    template <typename MapType>
    void setReference(const MapType& map)
    {
      MapType map2 = map; // avoid const issue in MapConversion::convert
      MapConversion::convert(0, map2, reference_, max_num_peaks_considered_);
    }

protected:

    void updateMembers_() override;

    PoseClusteringAffineSuperimposer superimposer_;

    StablePairFinder pairfinder_;

    ConsensusMap reference_;

    Int max_num_peaks_considered_;

private:

    /// Copy constructor intentionally not implemented -> private
    DiaWeaverFeatureClustering(const DiaWeaverFeatureClustering&);
    /// Assignment operator intentionally not implemented -> private
    DiaWeaverFeatureClustering& operator=(const DiaWeaverFeatureClustering&);
  };

} // namespace OpenMS
