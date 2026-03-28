// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Mohammed Alhigaylan $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/DiaWeaverFeatureClustering.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Feature.h>

using namespace std;

namespace OpenMS
{

  DiaWeaverFeatureClustering::DiaWeaverFeatureClustering() :
    DefaultParamHandler("DiaWeaverFeatureClustering"),
    ProgressLogger(),
    max_num_peaks_considered_(0)
  {
    defaults_.insert("superimposer:", PoseClusteringAffineSuperimposer().getParameters());
    defaults_.insert("pairfinder:", StablePairFinder().getParameters());
    defaults_.setValue("max_num_peaks_considered", 1000,
      "The maximal number of peaks/features to be considered per map. To use all, set to '-1'.");
    defaults_.setMinInt("max_num_peaks_considered", -1);

    defaultsToParam_();
  }

  DiaWeaverFeatureClustering::~DiaWeaverFeatureClustering() = default;

  void DiaWeaverFeatureClustering::updateMembers_()
  {
    superimposer_.setParameters(param_.copy("superimposer:", true));
    superimposer_.setLogType(getLogType());

    pairfinder_.setParameters(param_.copy("pairfinder:", true));
    pairfinder_.setLogType(getLogType());

    max_num_peaks_considered_ = param_.getValue("max_num_peaks_considered");
  }

  void DiaWeaverFeatureClustering::align(const FeatureMap& map, TransformationDescription& trafo)
  {
    ConsensusMap map_scene;
    MapConversion::convert(1, map, map_scene, max_num_peaks_considered_);
    align(map_scene, trafo);
  }

  void DiaWeaverFeatureClustering::align(const PeakMap& map, TransformationDescription& trafo)
  {
    ConsensusMap map_scene;
    PeakMap map2(map);
    MapConversion::convert(1, map2, map_scene, max_num_peaks_considered_);
    align(map_scene, trafo);
  }

  void DiaWeaverFeatureClustering::align(const ConsensusMap& map, TransformationDescription& trafo)
  {
    const ConsensusMap& map_model = reference_;
    ConsensusMap map_scene = map;

    // Run superimposer to find the global affine transformation
    TransformationDescription si_trafo;
    superimposer_.run(map_model, map_scene, si_trafo);

    // Apply global transformation to consensus features and their handles
    for (Size j = 0; j < map_scene.size(); ++j)
    {
      double rt = map_scene[j].getRT();
      rt = si_trafo.apply(rt);
      map_scene[j].setRT(rt);
      map_scene[j].begin()->asMutable().setRT(rt);
    }

    // Run pairfinder to find feature pairs between reference and (globally transformed) scene
    ConsensusMap result;
    std::vector<ConsensusMap> input(2);
    input[0] = map_model;
    input[1] = map_scene;
    pairfinder_.run(input, result);

    // Build local transformation from matched pairs, undoing the global transform
    si_trafo.invert();
    TransformationDescription::DataPoints data;
    for (ConsensusFeature& cfeature : result)
    {
      if (cfeature.size() == 2)
      {
        ConsensusFeature::iterator feat_it = cfeature.begin();
        double y = feat_it->getRT();
        double x = si_trafo.apply((++feat_it)->getRT());
        if (feat_it->getMapIndex() != 0)
        {
          data.push_back(make_pair(x, y));
        }
        else
        {
          data.push_back(make_pair(y, x));
        }
      }
    }

    trafo = TransformationDescription(data);
    trafo.fitModel("linear");
  }

  void DiaWeaverFeatureClustering::loadMzMLAsFeatureMap(const String& mzml_file, FeatureMap& feature_map)
  {
    feature_map.clear(true);

    MSExperiment exp;
    MzMLFile mzml;
    mzml.load(mzml_file, exp);

    for (const MSSpectrum& spec : exp)
    {
      if (spec.getMSLevel() != 2) continue;
      if (spec.getPrecursors().empty()) continue;

      const Precursor& prec = spec.getPrecursors()[0];
      Feature f;
      f.setRT(spec.getRT());
      f.setMZ(prec.getMZ());
      f.setCharge(prec.getCharge());
      f.setIntensity(prec.getIntensity());

      feature_map.push_back(std::move(f));
    }

    feature_map.updateRanges();
  }

} // namespace OpenMS
