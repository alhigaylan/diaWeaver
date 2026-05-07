// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/RedeemBatchPredictor.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/CsvFile.h>

#include <algorithm>
#include <array>
#include <fstream>
#include <initializer_list>
#include <map>
#include <optional>
#include <vector>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_RedeemTransitionPredictor RedeemTransitionPredictor

@brief Predicts RT, MS2 fragment intensities, and optional ion mobility from a peptide TSV using ReDeem and writes an OpenSWATH-style transition TSV.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; RedeemTransitionPredictor &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> Digestor, external peptide TSV generation </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> TargetedFileConverter, OpenSwathAssayGenerator, OpenSwathDecoyGenerator </td>
        </tr>
    </table>
</CENTER>

The input is a peptide-centric TSV with one row per precursor. The tool requires
`precursor_charge` plus either `modified_peptide_sequence` or `peptide_sequence`.
Optional metadata columns such as `protein_accession`, `gene_name`, `nce`, and
`instrument` are carried through when present.

Optionally, a second transition-level TSV can be supplied via `fine_tune_tsv` to
fine-tune the loaded ReDeem models before prediction. Fine-tuning is performed
in-place on the loaded RT, CCS, and MS2 models unless individual model types are
disabled. Fine-tuned checkpoints can optionally be saved to explicit output paths.

The output is an OpenSWATH-style transition TSV. The RT column is written as
`RetentionTime` rather than `NormalizedRetentionTime`, since ReDeem currently
predicts local chromatographic RT rather than iRT. Predicted ion mobility is
written to `PrecursorIonMobility` as Bruker-style `1/k0` converted from the
predicted CCS.

MS2 fragment rows follow the same mapping used in easypqp for the non-neutral-loss
channels: `b+1`, `b+2`, `y+1`, `y+2`. Neutral-loss channels remain available in
the underlying ReDeem API but are not written into the transition list in v1.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_RedeemTransitionPredictor.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_RedeemTransitionPredictor.html
*/

/// @cond TOPPCLASSES

namespace
{
  struct PeptideRecord
  {
    Size input_line = 0;
    String stripped_sequence;
    String canonical_modified_peptide;
    Int precursor_charge = 0;
    Int nce = 30;
    String instrument = "Lumos";
    vector<String> protein_ids;
    String gene_name;
    String peptide_group_label;
    String transition_group_id;
  };

  String normalizeHeader_(const String& value)
  {
    String normalized = value;
    normalized.trim();
    return normalized.toLower();
  }

  optional<Size> findHeaderIndex_(
    const map<String, Size>& header_map,
    initializer_list<const char*> candidates)
  {
    for (const auto* candidate : candidates)
    {
      const auto it = header_map.find(String(candidate));
      if (it != header_map.end())
      {
        return it->second;
      }
    }
    return nullopt;
  }

  String getCell_(const StringList& row, const optional<Size>& index)
  {
    if (!index.has_value() || *index >= row.size())
    {
      return {};
    }
    String value = row[*index];
    value.trim();
    return value;
  }

  Int parseRequiredIntCell_(const String& value, const String& column_name, Size input_line)
  {
    if (value.empty())
    {
      throw Exception::IllegalArgument(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        "Missing required integer column '" + column_name + "' on input line " + String(input_line) + ".");
    }

    try
    {
      return value.toInt();
    }
    catch (const Exception::ConversionError&)
    {
      throw Exception::IllegalArgument(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        "Could not parse integer column '" + column_name + "' with value '" + value + "' on input line " + String(input_line) + ".");
    }
  }

  Int parseOptionalIntCell_(
    const String& value,
    const String& column_name,
    Size input_line,
    Int fallback)
  {
    if (value.empty())
    {
      return fallback;
    }

    try
    {
      return value.toInt();
    }
    catch (const Exception::ConversionError&)
    {
      throw Exception::IllegalArgument(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        "Could not parse integer column '" + column_name + "' with value '" + value + "' on input line " + String(input_line) + ".");
    }
  }

  vector<String> splitListField_(const String& raw)
  {
    vector<String> values;
    if (raw.empty())
    {
      return values;
    }

    StringList split_values;
    if (raw.hasSubstring(";"))
    {
      raw.split(';', split_values);
    }
    else
    {
      split_values.push_back(raw);
    }

    for (const auto& value : split_values)
    {
      String trimmed = value;
      trimmed.trim();
      if (!trimmed.empty())
      {
        values.push_back(trimmed);
      }
    }

    return values;
  }

  String makeTransitionGroupId_(
    const String& canonical_modified_peptide,
    Int precursor_charge,
    map<String, Size>& seen_group_ids)
  {
    const String base_id = canonical_modified_peptide + "/" + String(precursor_charge);
    Size& count = seen_group_ids[base_id];
    ++count;
    if (count == 1)
    {
      return base_id;
    }
    return base_id + "#" + String(count);
  }

  String makeAnnotation_(const RedeemFragmentPrediction& fragment)
  {
    String annotation = fragment.ion_type + String(fragment.ordinal);
    if (fragment.charge > 1)
    {
      annotation += "^" + String(fragment.charge);
    }
    return annotation;
  }

  vector<RedeemPrediction> predictBatchWithDiagnostics_(
    const RedeemBatchPredictor& predictor,
    const vector<PeptideRecord>& records,
    Size begin,
    Size end)
  {
    vector<RedeemPrecursorInput> batch_inputs;
    batch_inputs.reserve(end - begin);
    for (Size i = begin; i < end; ++i)
    {
      batch_inputs.push_back(
        RedeemPrecursorInput{
          records[i].canonical_modified_peptide,
          records[i].precursor_charge,
          records[i].nce,
          records[i].instrument
        });
    }

    try
    {
      return predictor.predict(batch_inputs);
    }
    catch (const Exception::BaseException& batch_error)
    {
      for (Size i = 0; i < batch_inputs.size(); ++i)
      {
        try
        {
          predictor.predict({batch_inputs[i]});
        }
        catch (const Exception::BaseException& single_error)
        {
          throw Exception::IllegalArgument(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "Redeem prediction failed on input line " + String(records[begin + i].input_line)
              + " for precursor '" + records[begin + i].canonical_modified_peptide
              + "' (charge " + String(records[begin + i].precursor_charge) + "): "
              + String(single_error.what()));
        }
      }

      throw Exception::IllegalArgument(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        "Redeem batch prediction failed for input lines " + String(records[begin].input_line)
          + "-" + String(records[end - 1].input_line) + ": " + String(batch_error.what()));
    }
    catch (const exception& batch_error)
    {
      for (Size i = 0; i < batch_inputs.size(); ++i)
      {
        try
        {
          predictor.predict({batch_inputs[i]});
        }
        catch (const exception& single_error)
        {
          throw Exception::IllegalArgument(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "Redeem prediction failed on input line " + String(records[begin + i].input_line)
              + " for precursor '" + records[begin + i].canonical_modified_peptide
              + "' (charge " + String(records[begin + i].precursor_charge) + "): "
              + String(single_error.what()));
        }
      }

      throw Exception::IllegalArgument(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        "Redeem batch prediction failed for input lines " + String(records[begin].input_line)
          + "-" + String(records[end - 1].input_line) + ": " + String(batch_error.what()));
    }
  }

  void writeRow_(ofstream& os, const vector<String>& row)
  {
    for (Size i = 0; i < row.size(); ++i)
    {
      if (i != 0)
      {
        os << '\t';
      }
      os << row[i];
    }
    os << '\n';
  }
}

class TOPPRedeemTransitionPredictor :
  public TOPPBase,
  public ProgressLogger
{
public:
  TOPPRedeemTransitionPredictor() :
    TOPPBase("RedeemTransitionPredictor", "Predicts RT, MS2 fragment intensities, and optional ion mobility from a peptide TSV using ReDeem.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_(
      "in",
      "<file>",
      "",
      "Input peptide TSV with one row per precursor. Requires 'precursor_charge' and either 'modified_peptide_sequence' or 'peptide_sequence'.");
    setValidFormats_("in", {"tsv"});

    registerOutputFile_("out", "<file>", "", "Output OpenSWATH-style transition TSV.");
    setValidFormats_("out", {"tsv"});

    registerIntOption_("batch_size", "<int>", 512, "Number of precursors per ReDeem prediction batch.", false);
    setMinInt_("batch_size", 1);

    registerIntOption_("nce", "<int>", 30, "Default normalized collision energy when the input TSV does not provide an NCE column.", false);
    setMinInt_("nce", 0);
    registerStringOption_("instrument", "<str>", "Lumos", "Default instrument label when the input TSV does not provide an instrument column.", false);
    registerStringOption_("device", "<str>", "cpu", "ReDeem device preference passed to the Rust bridge. Supported values are 'cpu', 'auto', 'cuda', or 'cuda:N'.", false);
    registerFlag_("disable_ccs", "Disable CCS prediction and omit PrecursorIonMobility output.", false);
    registerDoubleOption_("library_intensity_scale", "<float>", 10000.0, "Scale factor applied to predicted library intensities before writing output TSV rows.", false);
    setMinFloat_("library_intensity_scale", 0.0);
    registerDoubleOption_("min_library_intensity", "<float>", 0.0, "Skip fragment rows whose predicted library intensity is below this threshold.", false);
    setMinFloat_("min_library_intensity", 0.0);

    registerInputFile_("rt_model", "<file>", "", "Optional explicit RT model path.", false, true);
    registerInputFile_("ccs_model", "<file>", "", "Optional explicit CCS model path.", false, true);
    registerInputFile_("ms2_model", "<file>", "", "Optional explicit MS2 model path.", false, true);

    registerTOPPSubsection_("fine_tune", "Advanced options for optional ReDeem fine-tuning before prediction.");
    registerInputFile_(
      "fine_tune:tsv",
      "<file>",
      "",
      "Optional transition-level TSV used to fine-tune the loaded ReDeem models before prediction. Expected format is one row per fragment transition, e.g. an OpenSWATH-style assay table. Required columns are a peptide column such as 'sequence' or 'modified_peptide_sequence' (using UniMod:X annotation) plus 'precursor_charge'. RT fine-tuning additionally requires 'retention_time'. CCS fine-tuning requires either 'ccs' or both 'ion_mobility' and 'precursor_mz'. MS2 fine-tuning requires fragment columns such as 'fragment_type', 'fragment_series_number', 'product_charge', and 'intensity'.",
      false,
      true);
    setValidFormats_("fine_tune:tsv", {"tsv", "csv"});
    registerInputFile_(
      "fine_tune:validation_tsv",
      "<file>",
      "",
      "Optional validation TSV used during ReDeem fine-tuning. Must use the same transition-level schema and required columns as fine_tune:tsv. When unset, the tool uses an internal validation split from fine_tune:tsv.",
      false,
      true);
    setValidFormats_("fine_tune:validation_tsv", {"tsv", "csv"});
    registerDoubleOption_("fine_tune:validation_fraction", "<float>", 0.2, "Validation fraction used when fine_tune:validation_tsv is not provided.", false, true);
    setMinFloat_("fine_tune:validation_fraction", 0.0);
    setMaxFloat_("fine_tune:validation_fraction", 0.95);
    registerIntOption_("fine_tune:batch_size", "<int>", 64, "Batch size used during ReDeem fine-tuning.", false, true);
    setMinInt_("fine_tune:batch_size", 1);
    registerIntOption_("fine_tune:epochs", "<int>", 10, "Maximum number of fine-tuning epochs.", false, true);
    setMinInt_("fine_tune:epochs", 1);
    registerIntOption_("fine_tune:early_stopping_patience", "<int>", 5, "Early stopping patience in epochs during fine-tuning.", false, true);
    setMinInt_("fine_tune:early_stopping_patience", 1);
    registerDoubleOption_("fine_tune:learning_rate", "<float>", 1e-4, "Learning rate used during ReDeem fine-tuning.", false, true);
    setMinFloat_("fine_tune:learning_rate", 0.0);
    registerDoubleOption_("fine_tune:warmup_fraction", "<float>", 0.0, "Optional warmup fraction used during ReDeem fine-tuning.", false, true);
    setMinFloat_("fine_tune:warmup_fraction", 0.0);
    setMaxFloat_("fine_tune:warmup_fraction", 1.0);
    registerFlag_("fine_tune:disable_rt", "Skip RT fine-tuning even when fine_tune:tsv is provided.", true);
    registerFlag_("fine_tune:disable_ccs", "Skip CCS fine-tuning even when fine_tune:tsv is provided.", true);
    registerFlag_("fine_tune:disable_ms2", "Skip MS2 fine-tuning even when fine_tune:tsv is provided.", true);
    registerStringOption_("fine_tune:rt_model_out", "<file>", "", "Optional output path for the fine-tuned RT model (.safetensors).", false, true);
    registerStringOption_("fine_tune:ccs_model_out", "<file>", "", "Optional output path for the fine-tuned CCS model (.safetensors).", false, true);
    registerStringOption_("fine_tune:ms2_model_out", "<file>", "", "Optional output path for the fine-tuned MS2 model (.safetensors).", false, true);
  }

  ExitCodes main_(int, const char**) override
  {
    const String in = getStringOption_("in");
    const String out = getStringOption_("out");
    const Int batch_size = getIntOption_("batch_size");
    const Int default_nce = getIntOption_("nce");
    const String default_instrument = getStringOption_("instrument");
    const String device = getStringOption_("device");
    const bool enable_ccs = !getFlag_("disable_ccs");
    const double library_intensity_scale = getDoubleOption_("library_intensity_scale");
    const double min_library_intensity = getDoubleOption_("min_library_intensity");
    const String fine_tune_tsv = getStringOption_("fine_tune:tsv");
    const String fine_tune_validation_tsv = getStringOption_("fine_tune:validation_tsv");
    const String fine_tuned_rt_model = getStringOption_("fine_tune:rt_model_out");
    const String fine_tuned_ccs_model = getStringOption_("fine_tune:ccs_model_out");
    const String fine_tuned_ms2_model = getStringOption_("fine_tune:ms2_model_out");

    if (fine_tune_tsv.empty() && (!fine_tune_validation_tsv.empty() || !fine_tuned_rt_model.empty() || !fine_tuned_ccs_model.empty() || !fine_tuned_ms2_model.empty()))
    {
      throw Exception::IllegalArgument(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        "fine_tune:validation_tsv and fine_tune:*_model_out options require fine_tune:tsv.");
    }

    CsvFile peptide_table(in, '\t');
    if (peptide_table.rowCount() < 2)
    {
      writeLogError_("Input TSV does not contain any peptide rows.\n");
      return PARSE_ERROR;
    }

    StringList header;
    peptide_table.getRow(0, header);
    map<String, Size> header_map;
    for (Size i = 0; i < header.size(); ++i)
    {
      header_map[normalizeHeader_(header[i])] = i;
    }

    const optional<Size> modified_peptide_col = findHeaderIndex_(
      header_map,
      {"modified_peptide_sequence", "modifiedpeptidesequence", "modifiedsequence", "fullpeptidename", "fullunimodpeptidename"});
    const optional<Size> peptide_col = findHeaderIndex_(
      header_map,
      {"peptide_sequence", "peptidesequence", "sequence", "strippedsequence"});
    const optional<Size> precursor_charge_col = findHeaderIndex_(
      header_map,
      {"precursor_charge", "precursorcharge", "charge"});
    const optional<Size> protein_col = findHeaderIndex_(
      header_map,
      {"protein_accession", "proteinaccession", "proteinid", "proteinname"});
    const optional<Size> gene_col = findHeaderIndex_(
      header_map,
      {"gene_name", "genename"});
    const optional<Size> nce_col = findHeaderIndex_(
      header_map,
      {"nce", "collisionenergy", "ce"});
    const optional<Size> instrument_col = findHeaderIndex_(
      header_map,
      {"instrument"});

    if (!precursor_charge_col.has_value())
    {
      writeLogError_("Input TSV is missing a precursor charge column. Expected one of: precursor_charge, PrecursorCharge, charge.\n");
      return PARSE_ERROR;
    }
    if (!modified_peptide_col.has_value() && !peptide_col.has_value())
    {
      writeLogError_("Input TSV is missing a peptide sequence column. Expected one of: modified_peptide_sequence, ModifiedPeptideSequence, peptide_sequence, PeptideSequence.\n");
      return PARSE_ERROR;
    }

    vector<PeptideRecord> records;
    records.reserve(peptide_table.rowCount() - 1);
    map<String, Size> seen_group_ids;
    map<String, Size> duplicate_precursors;

    for (Size row_index = 1; row_index < peptide_table.rowCount(); ++row_index)
    {
      StringList row;
      peptide_table.getRow(row_index, row);
      const Size input_line = row_index + 1;

      const String modified_value = getCell_(row, modified_peptide_col);
      const String stripped_value = getCell_(row, peptide_col);
      const String peptide_value = !modified_value.empty() ? modified_value : stripped_value;

      if (peptide_value.empty())
      {
        throw Exception::IllegalArgument(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Missing peptide sequence on input line " + String(input_line) + ".");
      }

      const Int precursor_charge = parseRequiredIntCell_(
        getCell_(row, precursor_charge_col),
        "precursor_charge",
        input_line);
      if (precursor_charge <= 0)
      {
        throw Exception::IllegalArgument(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Precursor charge must be positive on input line " + String(input_line) + ".");
      }

      const Int nce = parseOptionalIntCell_(
        getCell_(row, nce_col),
        "nce",
        input_line,
        default_nce);
      const String instrument = [&]()
      {
        const String value = getCell_(row, instrument_col);
        return value.empty() ? default_instrument : value;
      }();

      AASequence sequence;
      try
      {
        sequence = AASequence::fromString(peptide_value);
      }
      catch (const Exception::BaseException& e)
      {
        throw Exception::IllegalArgument(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Could not parse peptide sequence '" + peptide_value + "' on input line " + String(input_line) + ": " + String(e.what()));
      }

      if (sequence.size() < 2)
      {
        throw Exception::IllegalArgument(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Peptides shorter than length 2 are not supported for MS2 prediction on input line " + String(input_line) + ".");
      }

      String canonical_modified_peptide;
      try
      {
        canonical_modified_peptide = RedeemBatchPredictor::canonicalizeModifiedPeptide(peptide_value);
      }
      catch (const Exception::BaseException& e)
      {
        throw Exception::IllegalArgument(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Could not canonicalize peptide '" + peptide_value + "' on input line " + String(input_line) + ": " + String(e.what()));
      }

      const String duplicate_key = canonical_modified_peptide + "/" + String(precursor_charge);
      ++duplicate_precursors[duplicate_key];

      PeptideRecord record;
      record.input_line = input_line;
      record.stripped_sequence = sequence.toUnmodifiedString();
      record.canonical_modified_peptide = canonical_modified_peptide;
      record.precursor_charge = precursor_charge;
      record.nce = nce;
      record.instrument = instrument;
      record.protein_ids = splitListField_(getCell_(row, protein_col));
      record.gene_name = getCell_(row, gene_col);
      record.peptide_group_label = record.stripped_sequence;
      record.transition_group_id = makeTransitionGroupId_(
        record.canonical_modified_peptide,
        record.precursor_charge,
        seen_group_ids);
      records.push_back(std::move(record));
    }

    Size duplicate_precursor_count = 0;
    for (const auto& [_, count] : duplicate_precursors)
    {
      if (count > 1)
      {
        ++duplicate_precursor_count;
      }
    }
    if (duplicate_precursor_count > 0)
    {
      writeLogWarn_(
        "Input TSV contains " + String(duplicate_precursor_count)
        + " duplicated precursor keys (modified peptide + charge). The tool keeps all rows and disambiguates TransitionGroupId values.\n");
    }

    RedeemBatchPredictor::Config config;
    if (!getStringOption_("rt_model").empty())
    {
      config.rt_model_path = getStringOption_("rt_model");
    }
    if (!getStringOption_("ccs_model").empty())
    {
      config.ccs_model_path = getStringOption_("ccs_model");
    }
    if (!getStringOption_("ms2_model").empty())
    {
      config.ms2_model_path = getStringOption_("ms2_model");
    }
    config.enable_ccs = enable_ccs;
    config.device_preference = device;

    RedeemBatchPredictor predictor(config);

    if (!fine_tune_tsv.empty())
    {
      RedeemFineTuneConfig fine_tune_config;
      fine_tune_config.training_tsv_path = fine_tune_tsv;
      if (!fine_tune_validation_tsv.empty())
      {
        fine_tune_config.validation_tsv_path = fine_tune_validation_tsv;
      }
      fine_tune_config.validation_fraction = getDoubleOption_("fine_tune:validation_fraction");
      fine_tune_config.batch_size = static_cast<Size>(getIntOption_("fine_tune:batch_size"));
      fine_tune_config.epochs = static_cast<Size>(getIntOption_("fine_tune:epochs"));
      fine_tune_config.early_stopping_patience = static_cast<Size>(getIntOption_("fine_tune:early_stopping_patience"));
      fine_tune_config.learning_rate = getDoubleOption_("fine_tune:learning_rate");
      fine_tune_config.warmup_fraction = getDoubleOption_("fine_tune:warmup_fraction");
      fine_tune_config.default_nce = default_nce;
      fine_tune_config.default_instrument = default_instrument;
      fine_tune_config.enable_rt = !getFlag_("fine_tune:disable_rt");
      fine_tune_config.enable_ccs = enable_ccs && !getFlag_("fine_tune:disable_ccs");
      fine_tune_config.enable_ms2 = !getFlag_("fine_tune:disable_ms2");

      if (!fine_tuned_rt_model.empty())
      {
        fine_tune_config.rt_model_output_path = fine_tuned_rt_model;
      }
      if (!fine_tuned_ccs_model.empty())
      {
        fine_tune_config.ccs_model_output_path = fine_tuned_ccs_model;
      }
      if (!fine_tuned_ms2_model.empty())
      {
        fine_tune_config.ms2_model_output_path = fine_tuned_ms2_model;
      }

      if (!fine_tune_config.enable_rt && !fine_tune_config.enable_ccs && !fine_tune_config.enable_ms2)
      {
        throw Exception::IllegalArgument(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "fine_tune_tsv was provided, but all fine-tuning targets were disabled.");
      }

      writeLogInfo_(
        "Fine-tuning ReDeem models from " + fine_tune_tsv
        + " before transition prediction.\n");
      predictor.fineTuneFromTransitionTsv(fine_tune_config);
      writeLogInfo_("Finished ReDeem fine-tuning.\n");
    }

    ofstream os(out.c_str());
    if (!os)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out);
    }
    os.precision(writtenDigits(double()));

    const array<String, 29> header_names{
      "PrecursorMz",
      "ProductMz",
      "PrecursorCharge",
      "ProductCharge",
      "LibraryIntensity",
      "RetentionTime",
      "PeptideSequence",
      "ModifiedPeptideSequence",
      "PeptideGroupLabel",
      "LabelType",
      "CompoundName",
      "SumFormula",
      "SMILES",
      "Adducts",
      "ProteinId",
      "UniprotId",
      "GeneName",
      "FragmentType",
      "FragmentSeriesNumber",
      "Annotation",
      "CollisionEnergy",
      "PrecursorIonMobility",
      "TransitionGroupId",
      "TransitionId",
      "Decoy",
      "DetectingTransition",
      "IdentifyingTransition",
      "QuantifyingTransition",
      "Peptidoforms"
    };
    writeRow_(os, vector<String>(header_names.begin(), header_names.end()));

    startProgress(0, records.size(), "Predicting ReDeem transition rows");
    Size written_rows = 0;
    for (Size batch_begin = 0; batch_begin < records.size(); batch_begin += static_cast<Size>(batch_size))
    {
      const Size batch_end = min(batch_begin + static_cast<Size>(batch_size), records.size());
      const auto predictions = predictBatchWithDiagnostics_(predictor, records, batch_begin, batch_end);

      for (Size offset = 0; offset < predictions.size(); ++offset)
      {
        const auto& record = records[batch_begin + offset];
        const auto& prediction = predictions[offset];
        const String protein_ids = ListUtils::concatenate(record.protein_ids, ";");
        const String mobility = prediction.mobility_1k0.has_value() ? String(*prediction.mobility_1k0) : String();

        for (const auto& fragment : prediction.annotated_fragments)
        {
          const double scaled_intensity = static_cast<double>(fragment.intensity) * library_intensity_scale;
          if (scaled_intensity < min_library_intensity)
          {
            continue;
          }

          const String annotation = makeAnnotation_(fragment);
          const String transition_id = record.transition_group_id + "_" + annotation;

          writeRow_(
            os,
            {
              String(prediction.precursor_mz),
              String(fragment.mz),
              String(record.precursor_charge),
              String(fragment.charge),
              String(scaled_intensity),
              String(prediction.rt),
              record.stripped_sequence,
              record.canonical_modified_peptide,
              record.peptide_group_label,
              "",
              "",
              "",
              "",
              "",
              protein_ids,
              protein_ids,
              record.gene_name,
              fragment.ion_type,
              String(fragment.ordinal),
              annotation,
              String(record.nce),
              mobility,
              record.transition_group_id,
              transition_id,
              "0",
              "1",
              "0",
              "1",
              ""
            });
          ++written_rows;
        }
      }

      setProgress(batch_end);
    }
    endProgress();

    if (written_rows == 0)
    {
      writeLogWarn_("No transition rows were written. Check the input peptides and the intensity threshold.\n");
    }
    else
    {
      writeLogInfo_("Wrote " + String(written_rows) + " transition rows for " + String(records.size()) + " precursors.\n");
    }

    return EXECUTION_OK;
  }
};

/// @endcond

int main(int argc, const char** argv)
{
  TOPPRedeemTransitionPredictor tool;
  return tool.main(argc, argv);
}
