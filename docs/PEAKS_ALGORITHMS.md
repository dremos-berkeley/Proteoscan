# PEAKS Studio DDA identification: published algorithms, exports, and reimplementation notes

Prepared 2026-09-25 for an open-source reimplementation that imitates the PEAKS Studio (Bioinformatics Solutions Inc., BSI) DDA workflow: de novo, then PEAKS DB, PEAKS PTM and SPIDER, then protein inference and export.

**Evidence labels used throughout**

- **[D]** Documented. Stated in a peer-reviewed paper, an official BSI manual or web page, or an official BSI brochure. The source is cited next to the claim.
- **[O]** Observed. Seen in real PEAKS export files that third parties reproduce in open-source test data or parsers. This shows what PEAKS emits, not why it emits it.
- **[I]** Inferred. My reading or recommendation. It is not stated by BSI. Treat it as a hypothesis to test against real PEAKS output.

No PEAKS binary was downloaded or inspected. All manuals were read as PDFs from bioinfor.com. There is no public PEAKS Studio 13/13.5 manual: the PEAKS 12 manual (June 2024) is the newest one online, and URL guesses for 12.5, 13 and 13.5 manuals all returned 404.

---

## Source key

| Key | Source |
|---|---|
| S1 | Ma B, Zhang K, Hendrie C, Liang C, Li M, Doherty-Kirby A, Lajoie G. *PEAKS: powerful software for peptide de novo sequencing by tandem mass spectrometry.* Rapid Commun Mass Spectrom 2003;17:2337-2342. doi:10.1002/rcm.1196. https://pubmed.ncbi.nlm.nih.gov/14558135/ |
| S2 | Zhang J, Xin L, Shan B, Chen W, Xie M, Yuen D, Zhang W, Zhang Z, Lajoie GA, Ma B. *PEAKS DB: de novo sequencing assisted database search for sensitive and accurate peptide identification.* Mol Cell Proteomics 2012;11:M111.010587. https://pmc.ncbi.nlm.nih.gov/articles/PMC3322562/ |
| S3 | Ma B, Johnson R. *De novo sequencing and homology searching.* Mol Cell Proteomics 2012;11:O111.014902. https://pmc.ncbi.nlm.nih.gov/articles/PMC3277775/ |
| S4 | Han X, He L, Xin L, Shan B, Ma B. *PeaksPTM: mass spectrometry-based identification of peptides with unspecified modifications.* J Proteome Res 2011;10:2930-2936. doi:10.1021/pr200153k |
| S5 | Han Y, Ma B, Zhang K. *SPIDER: software for protein identification from sequence tags with de novo sequencing error.* J Bioinform Comput Biol 2005;3:697-716. doi:10.1142/s0219720005001247 |
| S6 | PEAKS Studio 12 User Manual (June 24, 2024). https://www.bioinfor.com/wp-content/uploads/2024/07/PEAKSStudio12_UserManual.pdf |
| S7 | PEAKS Studio 11 User Manual (Jan 4, 2023). https://www.bioinfor.com/wp-content/uploads/2023/01/PEAKSStudio11_UserManual.pdf |
| S8 | PEAKS Studio 8.5 User Manual. https://www.bioinfor.com/wp-content/uploads/2017/08/PEAKS_Studio_85_Manual.pdf |
| S9 | PEAKS Online 11 User Manual. https://www.bioinfor.com/wp-content/uploads/2023/09/PEAKSOnline11UserManual.pdf |
| S10 | BSI FDR / decoy-fusion tutorial. https://www.bioinfor.com/fdr-tutorial/ |
| S11 | Xin L, Qiao R, Chen X, et al., Li M. *A streamlined platform for analyzing tera-scale DDA and DIA mass spectrometry data enables highly sensitive immunopeptidomics.* Nat Commun 2022;13:3108. https://pmc.ncbi.nlm.nih.gov/articles/PMC9174175/ |
| S12 | PEAKS version history. https://www.bioinfor.com/versions/ |
| S13 | PEAKS Studio 13.5 release note. https://www.bioinfor.com/peaks-studio-13-5-release/ |
| S14 | PEAKS 11 product page. https://www.bioinfor.com/peaks-11/ |
| S15 | PEAKS Studio 13 brochure. https://www.bioinfor.com/wp-content/uploads/2025/06/PEAKS-Studio-13-Brochure.pdf |
| S16 | PEAKS CMD Linux User Guide (Dec 2015). https://www.bioinfor.com/wp-content/uploads/2016/12/PEAKSCMD-Linux-UserGuide.pdf |
| S17 | *Benchmarking Software for DDA-PASEF Immunopeptidomics.* Mol Cell Proteomics 2025, doi:10.1016/j.mcpro.2025.101492 (PMC13085059; bioRxiv 10.1101/2025.05.28.656277) |
| S18 | *Maximizing Immunopeptidomics-Based Bacterial Epitope Discovery by Multiple Search Engines and Rescoring.* J Proteome Res 2025, doi:10.1021/acs.jproteome.4c00864. https://pmc.ncbi.nlm.nih.gov/articles/PMC11976845/ |
| S19 | *Protein Identification Improvement in Complex Samples Using Higher Frequency MS Acquisition and PEAKS Software.* Appl Sci 2025;15(2):666. https://www.mdpi.com/2076-3417/15/2/666 (the publisher page returned HTTP 403; figures taken from search-engine abstracts) |
| S20 | ContraNovo, arXiv:2312.11584 (Table 1: nine-species benchmark including PEAKS). https://arxiv.org/abs/2312.11584 |
| S21 | Tran NH, Zhang X, Xin L, Shan B, Li M. *De novo peptide sequencing by deep learning.* PNAS 2017;114:8247. https://pmc.ncbi.nlm.nih.gov/articles/PMC5547637/ |
| S22 | Mao Z, Zhang R, Xin L, Li M. *Mitigating the missing fragmentation problem in de novo peptide sequencing with a two stage graph-based deep learning model (GraphNovo).* 2023, doi:10.21203/rs.3.rs-2593528/v1 |
| S23 | rusteomics/mzcore PEAKS parser and real-export test fixtures: `mzident/src/formats/peaks.rs`, `mzident/src/formats/peaks_tests.rs`, `fuzz/in_peaks/*.csv`. https://github.com/rusteomics/mzcore |
| S24 | patrick-willems/immunopeptidomics (`core.py`, `README.md`): PEAKS 12 `db.psms.csv` to Percolator conversion. https://github.com/patrick-willems/immunopeptidomics |
| S25 | zack-kirsch/timsplot `app.py` (PEAKS `db.psms.csv` import). https://github.com/zack-kirsch/timsplot |
| S26 | Parker R et al. *The Choice of Search Engine Affects Sequencing Depth and HLA Class I Allele-Specific Peptide Repertoires.* Mol Cell Proteomics 2021, doi:10.1016/j.mcpro.2021.100124 (PMC8724928) |
| S27 | PEAKS Studio product page (licensing, system requirements). https://www.bioinfor.com/peaks-studio/ |

---

## 1. PEAKS de novo sequencing, ALC and local confidence

### 1.1 Core algorithm (PEAKS "classic" de novo)

- **[D] Model (S1, S3).** PEAKS computes "the best peptide sequences whose fragment ions can best interpret the peaks in the MS/MS spectrum". It reports a confidence score for the whole sequence and "an additional novel positional scoring scheme for portions of the sequences" (S1 abstract).
- **[D] Dynamic programming over mass values (S3).** Unlike spectrum-graph methods, "PEAKS' algorithm performs dynamic programming on the mass values regardless of the presence of an observed fragment ion". The recurrence is BestScore(m) = f(m) + max_x BestScore(m − m(x)) over residue masses m(x). This lets PEAKS score through missing fragment ions.
- **[D] Two-round scoring (S3).** "A simple score function is used in the first round, which computes 10,000 sequences with the top matching scores. These candidates are further evaluated by a more sophisticated scoring function."
- **[D] Prefix and suffix odds functions (S2, describing the same scorer reused as the PEAKS DB "quick scorer").** A spectrum is converted into f_N(m), the odds that the correct peptide has an N-terminal prefix of residue mass m, and f_C(m), the same for a C-terminal suffix. The odds are estimated from fragment ions:
  - CID: a, b, c, y, z, b-H2O, y-H2O, y-NH3. This is the verbatim S2 list, even though c and z are unusual for CID.
  - ETD: a, b, c, c-H, y, z, z+H.
- **[D] Instrument-specific internal parameters (S6 §2.3, S8).** The "MS/MS Product Scan" type (e.g. IT, FT, LIT/FT) "will help PEAKS decide which internal parameters (for weighing fragments and amount of noise) to use during PEAKS auto de novo sequencing and a PEAKS DB search". The weights themselves are not published.
- **[D] Candidates reported.** The top 5 de novo candidates per spectrum are kept ("De novo all candidates"). The displayed candidate is the one with the highest ALC (S6 §4.4.3; S9 §3.2).
- **[D] Spectral pairs and triplets.** CID/HCD/ETD spectra of the same precursor are sequenced jointly into one row (S8 ch. 8 §6).

### 1.2 Local confidence and ALC

- **[D] Local confidence (LC)** is a per-residue value of 0-99% (S8) or 0-100 in exports. It indicates "how confident the algorithm considers a particular amino acid as correctly sequenced" (S8).
  - S3: "PEAKS computes a confidence score by comparing the top de novo sequencing peptide with the suboptimal peptides."
  - S9 (DeepNovo section, same concept): "Local confidence scores are calculated by permutation of the possible peptide sequences that can match the precursor mass and comparing the sequence alignments. The local confidence represents the likelihood that the amino acid reported at the position is correct compared to the other possible matches."
- **[I] Practical reimplementation of LC.** Take the top-N candidates (N of order 10^3-10^4 from round 1), convert their scores to weights (e.g. softmax of the round-2 score), and set LC(i) to the weight-fraction of candidates that place the same residue at the same mass interval as the top candidate. Mass-aligned positional agreement, not index agreement, is what "comparing the sequence alignments" implies.
- **[D] ALC** = (sum of residue LC) / peptide length (S6 §4.4.1).
  - The default ALC display threshold is 50% (S6 §4.3).
  - The DB-search report step also defaults "De Novo ALC (%)" for de novo-only results to 50% (S6 §6.1).
- **[D] Colour bins:** >90% red, 80-90% purple, 60-80% blue, <60% black (S6).
- **[D] Tag, mass tags and Tag Length.**
  - Residues whose LC is below a user threshold are collapsed into mass tags. The remaining residues form sequence tags.
  - **Tag Length** is "the length of the longest sequence tag in the peptide" at the current threshold (S8 §3.2.1; S9 §3.4). The `tag(>=X%)` export column shows the tagged sequence at threshold X.
  - **[O]** In PEAKS 13 exports the threshold in that header is 0.0% (`tag(>=0.0%)`), so Tag Length equals Length in those rows (S23 test data).
- **[D] Newer engines.**
  - DeepNovo (deep learning, S21) appears as "Deep Novo Score (%)" in PEAKS 12 de novo exports [O, S23].
  - PEAKS 12 adds GraphNovo, a two-stage graph neural network that first finds an optimal path, then predicts the sequence (S22; S6 §1.3). It requires a CUDA GPU with compute capability ≥ 8.
  - PEAKS 12 introduced de novo FDR, "the first of its kind", using decoy spectra (S6 ch. 12; S15). PEAKS 13.5 adds "DDA de novo FDR" with decoy generation (S12).
- **[D] DeepNovo "Score (%)".** "A peptide confidence score which adds penalty to ALC (%) based on the difference between predicted RT and observed RT" (S9). This is the precedent for RT-aware rescoring in PEAKS.

---

## 2. PEAKS DB (S2 unless noted)

### 2.1 Pipeline

1. **De novo sequencing** of every spectrum. Only the top de novo peptide per spectrum is used.
2. **Tag construction.** "Each stretch of adjacent amino acid residues with <30% confidence is replaced by a 'mass segment' that is equal to the total mass of the residues." [D]
3. **Protein shortlisting.**
   - Each tag is aligned against database peptides so that "the mass of each aligned block … is equal for both sequences". The **CAA score** is the number of common amino acids in that alignment.
   - During shortlisting a modified de novo residue may match an unmodified database residue. During final scoring it must match the same residue with the same modification.
   - Proteins are ranked by their best peptide CAA, with ties broken by the 2nd and then 3rd best.
   - **The top 7,000 proteins** form the shortlist. This number "can be adjusted in the configuration file of PEAKS DB". [D]
4. **Peptide shortlisting.**
   - All in-silico peptides from the shortlist are enumerated, including all combinations of the user's variable PTMs.
   - Each is scored against spectra with matching precursor mass using the quick scorer: the ion-match score is the sum of f_N(m) over prefix masses plus f_C(m′) over suffix masses.
   - A priority queue keeps **the top 512 candidates per spectrum**. [D]
5. **Peptide scoring (rerank).**
   - The ion-match score is normalised as **s′ = (s − μ)/σ**, where μ is the mean of the top 10 candidates and σ is the standard deviation of the top 150 candidates. [D]
   - Nine further features [D]:
     1. CAA score against the de novo tag.
     2. Protein feature: each protein scores the sum of its 3 highest peptide CAA scores; a peptide's feature is the maximum over the proteins that contain it.
     3. Peptide length.
     4. Average sequence length per missed cleavage.
     5. Average sequence length per PTM.
     6. Precursor mass error.
     7. Charge state.
     8. Maximum length of a consecutively matched fragment-ion series.
     9. Number of termini violating the enzyme rule.
   - "These nine features, together with the normalized ion match score, are combined with a weighted sum." The weights are trained per instrument type "to maximize the area on the left of the 1% FDR curve" and "do not change from experiment to experiment." [D]
   - **[D] (S11, PEAKS Online, a later generation):** "18 hand-selected features … A pre-trained linear discriminant function is applied on the 18 features". Also: "a small set of high quality spectra are searched against the filtered protein set. The match results are used to calibrate the mass error and train the scoring model", and "the overlap of the de novo sequencing result and peptide … added as a feature in the LDA model … significantly boosts the performance". The tag search "starts with the seeds constructed from de novo sequences then tries to extend the seed to the whole peptide", and is fast enough to filter NCBI-scale databases.
6. **Conversion to −10lgP.**
   - "The weighted sum score is converted to a p value … the probability that a false identification in the current search achieves the same or better matching score." The significance score is **−10·log10(P)**. [D]
   - So 20 corresponds to p = 0.01 and 15 to p ≈ 0.03 (S2; S8).
   - BSI stresses that p is a false-positive rate, not FDR (S8).
   - **[I]** The paper does not say how P is estimated. A faithful imitation fits the null distribution of the best-decoy or second-best score per search, e.g. an empirical CDF over decoy-fusion hits, a Gumbel or exponential tail fit, or the lower-ranked candidates of the 512-candidate queue. Then define −10lgP = −10·log10(1 − F_null(score)). PEAKS DB -10lgP values commonly reach 100-130 for excellent PSMs (observed maxima of about 133 in S23 fixtures) [O].
7. **Result validation by decoy fusion** (§2.2).
8. **Protein inference** (§4).

### 2.2 Decoy fusion

- **[D] Decoy construction.**
  - "The decoy sequences were produced by randomly shuffling the amino acids in each protein" (S2 Methods).
  - The target and the decoy are **concatenated into one entry per protein**, so the database has the same number of entries, each twice as long.
  - A separator **'J'** is inserted between them, and the digester cleaves on both sides of J. This keeps the target C-terminal peptide when the last target residue is not a cleavage site.
  - After the search, hits are labelled decoy or target by which half of the fused sequence they fall in.
  - **[I]** S2 does not say whether shuffling preserves K/R positions or protein termini. Plain random shuffling of the whole protein appears to be what was done.
- **[D] FDR formula.** FDR(T) = #decoy hits / #target hits at score ≥ T. It is not 2·D/(T+D) (S2; S10).
- **[D] Why decoy fusion is needed (S2; S10).**
  - Shortlisting selects "more target proteins than the decoy proteins", so later false hits land preferentially in targets.
  - The protein feature (feature 2) gives "bonuses" that favour targets.
  - Fusing equalises both effects: "the same amount of bonus will be equally applied to the target and decoy parts of the same fused sequence".
  - S2 Fig. 5: conventional TDA on PEAKS DB "generated fewer decoy hits than the false target hits … which may lead to FDR underestimation"; decoy fusion did not. For Mascot, SEQUEST and Percolator the two methods did not differ.
- **[D] Decoy identifiers in exports.** With the hidden `export-decoy` option (§8), decoy PSMs carry accessions prefixed `#DECOY#` and a boolean `decoy` column [O, S24].
- **[D] FDR levels (S6 §6.1, §6.3.6).**
  - **PSM FDR:** decoy assignments to spectra / target assignments to spectra. "PEAKS Studio keeps one peptide per spectrum at most" when counting PSMs.
  - **Peptide FDR:** "Peptides are scored based on their top PSM". FDR = decoy unique sequences / target unique sequences. The report step lets you choose PSM FDR, Peptide FDR or a raw PSM −10lgP threshold.
  - **Protein-group FDR:** "decoy protein groups where the top hit in the protein group is a decoy … relative to … target" groups.
  - Decoys are excluded from all displayed counts.
  - **[I]** No q-value monotonisation is described. Imitate it by scanning thresholds from high to low and reporting the lowest score threshold whose running D/T ≤ α.
- **[D] Independent check (S17).** An *Arabidopsis* entrapment test in DDA-PASEF immunopeptidomics found the target-decoy FDR of PEAKS 11 (with deep learning boost) and the other engines "unbiased and not overestimating", and not optimistic.

### 2.3 Published numbers (S2)

- **Settings:** 15 ppm precursor, 0.8 Da fragment (LTQ-Orbitrap CID/ETD), up to 3 missed cleavages, semi-specific allowed. Fixed carbamidomethyl-C; variable deamidation (NQ), Met-oxidation and pyro-Glu (Q).
- **P. aeruginosa CID at 1% FDR (PSMs):**

  | Engine | PSMs |
  |---|---|
  | PEAKS DB | 10,668 |
  | Mascot+Percolator | 9,969 |
  | SEQUEST | 8,236 |
  | Mascot | 7,515 |

  PEAKS DB found about 42% more than Mascot.
- **Yeast ETD (iPRG 2011) at 1% FDR:**

  | Engine | PSMs |
  |---|---|
  | PEAKS DB | 3,652 |
  | Mascot+Percolator | 2,702 |
  | Mascot | 2,398 |
  | SEQUEST | 2,233 |

---

## 3. PEAKS PTM and SPIDER (brief)

### PEAKS PTM (S4; S6 §6.1, §6.10; S8)

- **[D]** Unspecified-modification search on spectra that have good de novo tags but no confident PEAKS DB hit.
- By default it searches all "Common" and "Uncommon" built-in PTMs, which covers the natural modifications in Unimod: "313" per S15, "over 650 PTMs and mutations" per the PEAKS PTM web page.
- It restricts the search space to avoid exponential growth, and it can take a de novo ALC threshold.
- **[D]** PEAKS DB PTMs are automatically included, and the maximum number of variable PTMs per peptide defaults to the DB-search value (<4 recommended).
- **[D]** Results are merged with PEAKS DB under a single −10lgP scale.
- **[D]** Since PEAKS 12:
  - **Confident PTM:** at least X consecutive fragment ions at the site (default X = 3, or 3+1 at a terminus) with relative intensity ≥ Y (e.g. 2%) (S6 §2.5).
  - **Potential signature ions:** ≥ 30% of modified PSMs, and at least 30 percentage points more prevalent in modified than unmodified PSMs, are used for rescoring.
- **[D]** AScore = −10·log10(p) of site localisation versus alternative positions. It is exported as e.g. `C6:Carbamidomethylation:1000.00` (S6; S23).

### SPIDER (S5; S8)

- **[D]** Homology search of de novo tags that "deals with the possible overlaps between the de novo sequencing errors and the homology mutations" (S3). A segment of de novo residues may be replaced by a mass-equivalent segment while the search also allows point mutations.
- **[D]** Runs on spectra with ALC > 15% that remain unidentified, where PEAKS 8.5 used −10lgP < 30 as the "not confidently identified" cut (S8).
- **[D]** "For each spectrum, the better sequence constructed by SPIDER or found by PEAKS DB will be used."
- **[D]** Confident mutations use the same consecutive-ion rule as confident PTMs (S6 §2.5).

### PEAKS 13.5 additions (S13)

- **[D]** "Open PTM Search": a wide precursor tolerance on unmatched high-quality spectra.
- **[D]** "Multi-Round Search": re-searches unassigned spectra that have good de novo tags.

---

## 4. Protein inference, grouping, scores and the "Top" rule

- **[D] Protein −10lgP (S2; S6 §6.3).**
  - Remove redundant peptides. The same peptide seen in several spectra keeps only its best score. Peptides that differ only in PTM location count as the same; different sequences or different PTMs count as different.
  - Sort the scores s1 ≥ s2 ≥ … ≥ sk.
  - **Protein score = Σ_k s_k / k.**
  - A group's score equals the score of its dominating protein.
  - The default report filter is protein −10lgP ≥ 15 (S6). S2 uses T = 15 in grouping.
- **[D] Grouping by domination (S2).**
  - Given T = 15 (p ≈ 0.03), protein X dominates Y if every peptide of Y with score ≥ T is also in X. A dominated Y joins X's group.
  - Several proteins may dominate one another (equivalent peptide sets).
  - The UI can show one dominating protein, all dominating proteins, or all proteins.
- **[D] UI and export semantics (S6 §6.4).**
  - "Proteins are grouped based on parsimony."
  - **Top** proteins are "supported by the most unique peptides in the group". Proteins supported by a subset ("subproteins") are hidden in Top view.
  - A protein "supported by a peptide not supporting the top protein … will be added to a new group".
  - `Cluster` or `Protein Group` is the group identifier; `Top` is a flag.
- **[D] Unique peptide.** "A peptide that can be mapped to only one protein group." The minimum number of unique peptides is a report filter; UI defaults are commonly 1, but the manual does not state the value [I].
  - `#Peptides` and `#Unique` count "two peptides with the same starting and ending positions … as one, regardless of their PTM forms" (MCP guidelines).
  - `Unique` in `protein-peptides.csv` is `Y` or `N` [O, S23].
- **[D] Protein area.** "Total of all peptide features from unique supporting peptides" per sample (S6 §6.4.1).
- **[D] Protein FDR.** Groups whose top hit is a decoy, divided by groups whose top hit is a target.
- **[D] Coverage** counts supporting peptides only, not de novo tags. Low-confidence peptides are shown light blue.
- **[I] Reimplementation.**
  1. Build a peptide-to-protein bipartite graph from peptides that pass the peptide FDR, using peptides with score ≥ 15 for domination.
  2. Compute connected dominance groups.
  3. Pick the Top protein by most unique peptides, with ties broken by protein score.
  4. Compute uniqueness at the group level.
  5. Compute the 1/k-weighted score.

---

## 5. Feature detection, precursor correction and chimeras

- **[D] Data Refine (S6 §3.5; S9 §2.1).** The options are:
  - "Mass Correction": "By examining the isotope shapes in the corresponding MS scans, this function can accurately adjust the precursor m/z to be equal to the monoisotopic ion m/z."
  - "Associate feature with Chimera".
  - "Preprocessed Data - Skip Data Refinement".
  - "Filter Features".
- **[D] Older refinement (S8 ch. 7).**
  - Merge scans within an RT window and m/z tolerance; CID and HCD are not merged unless forced.
  - Correct precursor "Mass only" or "Mass and Charge states", with min and max charge.
  - Low-quality spectrum removal (quality threshold 0.65 recommended).
  - Centroiding, deisotoping and charge deconvolution of MS/MS, applied to all data.
  - Feature detection for LFQ used "an Expectation-Maximization (EM) based algorithm … for feature detection, deconvolution and refinement".
- **[D] Current generation (S11).** "A feature-based principle": "precursor features are first detected from LC-MS map; MS/MS spectra are then grouped with their corresponding precursors based on their co-elution profiles". Deep learning is used for "MS1 isotope feature detection" (DeepIso and the PointIso-type models by Zohora et al.).
- **[D] Chimera handling.**
  - With chimera association on, "PEAKS will assume that a tandem scan may contain two or more peptides… by searching for peptides that match the monoisotopic mass of any peptide feature in the acquisition window" (S6).
  - So one MS2 scan is searched once per candidate feature in its isolation window, and each feature supplies its own precursor mass and charge.
  - Summary statistics separate `#MS/MS` (excludes chimeric scans), `#Chimera`, `Identified #PSMs`, `Identified #Scans` and `Identified #Features` (S6 §6.3.3).
  - Older DB PSM exports carry a `from Chimera` Yes/No column [O, S23].
  - **[I]** A PSM therefore belongs to a (scan, feature) pair. #PSMs can exceed #Scans, and "one peptide per spectrum" applies per (scan, precursor) instance.
- **[D] Fake features.** "Features with Feature Id = 0 are known as fake features and are constructed to enhance identification" (S6 §6.5.3). This covers MS2 scans with no detected MS1 feature, which use the instrument-reported precursor. PEAKS DIA likewise "creates a fake peptide feature" (S11).
- **[D] Feature Id** is "an arbitrary number provided by PEAKS during data refinement for each identifiable feature". Older exports used `F<fraction>:<id>` (S6 §4.5.1; S23).
- **[D] Area** is "the area under the curve of the peptide feature found at the same m/z and retention time as the MS/MS scan" (S6).
  - Peptide `Area` per sample is "the combined area for all features in that sample associated with this peptide".
  - PEAKS Online can switch the abundance measure to feature intensity instead of area (S9).
- **[D] Calibration.** PEAKS "corrects m/z value for each PSM to correct systematic biases and instrument measurement drift" (S6 §6.3.7). The ppm column is 10^6 × (precursor mass − peptide mass) / peptide mass.
  - **[O]** Newer `db.psms.csv` files carry both `m/z` and `raw m/z`, and `RT` and `raw RT` (S25). The reported m/z is the calibrated or feature value; the raw value is the instrument's.
- **[I] Isotope-error correction.** No PEAKS manual offers an isotope-error search parameter such as MSFragger's −1/0/+1/+2. PEAKS relies on feature-based monoisotopic assignment ("Mass Correction") instead. To imitate it, assign the precursor from the MS1 feature's monoisotopic peak by isotope-envelope fitting, and search a single mass. Optionally allow 0/+1 only when no feature is found.
- **[D] Precursor charge.** Charge comes from the MS1 feature, or from the survey-scan isotope spacing (instrument "MS Precursor Scan" resolution setting, S6 §2.3). It is a scoring feature (#7) in PEAKS DB.

---

## 6. Deep Learning Boost, CAA, Tag Length and Positional Confidence

- **[D] Introduced in PEAKS Studio / Online 11 (Jan 2023)** as "Deep learning-boosted ID for PEAKS DB" (S12).
  - The manual text is only: "Check this option to use PEAKS enhanced deeplearning algorithm to better identify PSMs in the DB Search process" (S6 §6.1; S7).
  - Marketing claims: "can increase the identification over 10% for DDA data analysis" (S14; S9 §1).
  - The PEAKS 13 brochure: "deep learning boost for peptide confidence re-scoring … by using algorithms to predict retention time, fragment ion intensity, mass-to-charge ratios, and ion mobility" (S15).
- **[D] Features exposed per PSM or peptide (S6 §5.6.3, shown for the Peptide-search result but present in DB exports too [O]):**
  - **ΔRT / `Delta RT`:** "The difference between predicted RT and RT of peptide." [O] Values are in minutes, signed, e.g. −12.8 to +11.8 (S23).
  - **`MS2 Correlation`:** "The correlation of the real spectrum versus the predicted spectrum for the best PSM." [O] Range 0-1, typically 0.5-0.95. The metric is not stated. S11 used Pearson, cosine and spectral angle as three Percolator features in PEAKS Online rescoring, so **[I]** Pearson or cosine on matched fragment intensities is most likely.
  - **`Delta 1/k0`** for timsTOF: "difference between average predicted 1/k0 and 1/k0 of peptide".
  - **`CAA (%)`:** "The percent of number of confident amino acids in the sequence." It depends on the report-step parameter **Confident Amino Acid Threshold (%)**, which "sets the threshold for each amino acid local confidence for it to be considered confident when calculating the CAA column". The manual does not state the default value.
  - **`Tag Length`:** longest run of residues above the tag threshold. **[I]** In DB exports this is the de novo-style tag computed against the matched peptide's positional confidences.
  - **`Positional Confidence`:** [O] a space-separated per-residue integer vector in PEAKS 13 exports, e.g. `0 0 2 8 11 9 0 0 0 0 0 0` (S23, PEAKS 13 DIA export).
    - Its definition is not in any public manual. **[I]** It is probably a per-residue fragment-evidence confidence for the database peptide, analogous to de novo LC.
    - CAA does not equal the fraction of non-zero positions in all rows, so CAA likely uses a separate threshold or ion-evidence rule.
- **[D] Earlier precedents.**
  - PEAKS X+ de novo exports had `Denovo Score` and `Predict RT` [O, S23].
  - The DeepNovo Score penalises ALC by |ΔRT| (S9).
  - PEAKS Online rescoring added predicted-spectrum similarity features and used Percolator (S11). In three HLA DDA datasets, standard PEAKS Online beat MaxQuant 1.7-4.1× and its rescoring beat MaxQuant+Prosit 1.0-1.4×. It was comparable to MS²Rescore on 2 of 3 sets; MS²Rescore found 17% more on A*02:04.
- **[D] How the features feed rescoring.** BSI says only that PEAKS "uses internally calculated correlations between predicted and empirical properties (including delta RT and MS2 correlation) to re-score PSM −10log10(P) scores" (paraphrased in S19). S18 notes "the built-in 'Deep learning boost' already leverages spectrum and RT scoring features, explaining the high identification numbers prior MS²Rescore". The boosted score is still reported as `-10LgP`, so FDR is computed on the rescored value.
  - **[D] Reproducibility:** PEAKS Studio 11 and PEAKS Online 11 gave identical results, "Only minor variations … when activating the deep learning boost … attributable to the inherent randomness of neural network models" (S17). **[I]** This suggests per-run training, e.g. RT-model fine-tuning or a semi-supervised rescorer, rather than a fixed LDA.
- **[I] Recommended imitation.**
  1. Base features: PEAKS-DB-style LDA features (§2.1), de novo CAA and Tag Length, ppm, charge and length.
  2. Add ΔRT from an RT predictor (e.g. DeepLC or AlphaPeptDeep, calibrated per run on high-confidence PSMs, as |observed − aligned predicted|).
  3. Add MS2 correlation (Pearson on matched b/y intensities against Prosit, MS²PIP or AlphaPeptDeep predictions).
  4. Add Δ1/k0 for timsTOF data.
  5. Rescore with a semi-supervised linear SVM or LDA (Percolator- or mokapot-like) trained on decoy-fusion decoys.
  6. Map the output back to a −10lgP-like scale.

  S24 builds exactly this for PEAKS 12 exports: a mokapot PIN built from `-10LgP`, `Tag Length`, `Delta RT`, `MS2 Correlation`, `ppm` and `Delta 1/k0`.
- **[D] Benchmarks of the boost itself.**
  - More than 10% (S14).
  - Protein IDs +22.9% (CHO) and +23.2% (E. coli), i.e. about 20-26% and 18-28% ranges across conditions (S19).
  - PEAKS 11 Online (with boost) versus PEAKS X Pro gave +6-11% peptides for MHC-I and MHC-II (S17).

---

## 7. Export formats and column meanings

File names differ by version. PEAKS 11/12 auto-export and the Export node write `<analysis>.<kind>.csv`, and the DB-search node prefix is typically `db.`. The manual names the kinds `proteins.csv`, `peptides.csv`, `protein-peptides.csv`, `psms.csv`, `denovoOnly` and `features`. The names `db.ptmProfile.csv` and `db.modifiedPeptide.csv` correspond to the PTM Profiling CSV and the Modified Peptide view added in updated PEAKS 12 (S6 §6.6; S9 §3.2). **[I]** The exact file names for those two are not in the manual.

Columns come from S6, S9 and S8 definitions [D] plus real headers [O] from S23, S24 and S25.

### 7.1 `db.psms.csv` (one row per PSM, top PSM per spectrum/precursor; I/L variants split into separate rows [D S6 §6.9.1])

**[O] Headers seen:**
- PEAKS X / X Pro DB: `Peptide,-10lgP,Mass,Length,ppm,m/z,Z,RT,Area,Fraction,Id,Scan,from Chimera,Source File,Accession,PTM,AScore,Found By`
- PEAKS 12/13 with DL boost: `Peptide, -10LgP, Mass, Tag Length, CAA (%), Length, Delta RT, MS2 Correlation, ppm, m/z, raw m/z, z, RT, raw RT, Area, Scan, Feature Id, Source File, Accession, PTM, AScore, Positional Confidence, tag(>=0.0%)`, plus `Found By`, and for timsTOF `Precursor Id`, `1/k0 Range`, `Delta 1/k0`. The column set is from S25, but the order is partly inferred.
- With `export-decoy`, a `decoy` column is added (S24).

| Column | Meaning |
|---|---|
| Peptide | Sequence. Modifications are shown as `(+57.02)` mass deltas; SPIDER substitutions as `(sub X)` [O] |
| -10lgP / -10LgP | PSM significance, −10·log10 p (after DL boost when enabled). Capitalisation changed in 12+ [O] |
| Mass | Monoisotopic peptide mass |
| Length | Residues |
| ppm | 10^6·(precursor mass − peptide mass)/peptide mass (calibrated) |
| m/z, z | Precursor m/z (calibrated or feature) and charge; `raw m/z` and `raw RT` are the instrument values [O] |
| RT | Retention time (minutes) of the MS2 (or feature) |
| Area | Area of the associated feature; 0 or blank when there is no feature (fake feature) [O] |
| Fraction, Id | Fraction index; internal PSM or spectrum id (X era) [O] |
| Scan | Scan identifier: `F<frac>:<scan>` in X era; for mzML input in 12+, **spectrum index + 1**, not the native scan number [O, S24] |
| Feature Id | MS1 feature id; 0 means fake feature [D] |
| from Chimera | Yes/No, whether the PSM came from chimera association [O] |
| Source File | Raw or mzML filename |
| Accession | Protein accession(s), `;`/`:` separated; decoys use `#DECOY#` [O] |
| PTM | Modification names |
| AScore | `<AA><pos>:<PTM name>:<AScore>` entries (1000.00 is used for unambiguous sites) [O] |
| Found By | PEAKS DB / PEAKS PTM / SPIDER (12 peptidome: DB Search, DeepNovo, Homolog, Validated DeepNovo) [D] |
| Tag Length, CAA (%), Positional Confidence, tag(>=x%) | §6 |
| Delta RT, MS2 Correlation, Delta 1/k0 | §6 |

### 7.2 `db.peptides.csv` (one row per peptidoform, best PSM; same sequence with different PTMs is a separate row [D])

- **[O] X Pro header:** `Peptide,-10lgP,Mass,Length,ppm,m/z,RT,Area <sample…>,Fraction,Scan,Source File,#Feature,#Feature <sample…>,Accession,PTM,AScore,Found By`
- **[D] 12 adds:** CAA (%), ΔRT, MS2 Correlation, Gene, Database (Target / Non-Canonical / Contaminant), and #Spec in Online.
- **Meanings:**
  - Peptide −10lgP is the max over its PSMs.
  - `Area <sample>` is the sum of feature areas in that sample.
  - `#Feature` is the number of features assigned (or used in feature vectors for multi-sample runs).
  - `Scan` is the best-PSM scan.
  - `#Spec` is the number of MS2 spectra assigned.

### 7.3 `db.proteins.csv`

- **[D] Columns:** Accession, Cluster/Protein Group, Top, −10lgP, Coverage (%), #Peptides, #Unique, PTM, Avg. Mass, Description.
- **Per-sample optional columns** (always exported): Coverage <sample>, #Spec <sample>, Area <sample>.
- Filtered by the Summary protein filters: −10lgP, FDR, unique peptides.

### 7.4 `db.protein-peptides.csv` (one row per protein × supporting peptide)

- **[O] Header:** `Protein Group,Protein ID,Protein Accession,Peptide,Unique,-10lgP,Mass,Length,ppm,m/z,z,RT,Area <samples>,Fraction,Scan,Source File,#Feature,#Feature <samples>,Start,End,PTM,AScore,Found By`
- `Peptide` includes flanking residues (`K.PEPTIDE.R`).
- `Unique` is Y/N.
- `Start` and `End` are 1-based protein positions.
- It "usually contains more entries than the Peptide table" because shared peptides repeat [D].

### 7.5 `db.features.csv` (per sample; ignores view filters [D])

- **[O] PEAKS 11 de novo features header:** `Source File,Feature Id,m/z,z,RT,RT Begin,RT End,Area,Quality,Denovo Peptide,ALC (%),ppm`
- **[D] Meanings:**
  - `RT Begin` and `RT End` are the feature elution bounds.
  - `Quality` "may have quality score of 0 … no scoring" in 12.
  - The de novo peptide(s) are listed; multiple MS2 per feature are all listed.
  - A "Denovo Score (%)" combines ALC and RT difference.
- **[I]** The DB version replaces or augments these with the DB peptide and −10lgP of the best PSM on the feature.

### 7.6 `db.denovoOnly.csv` (per sample; spectra with ALC ≥ threshold and no passing DB/PTM/SPIDER hit [D])

- Same columns as the de novo table.
- **[O] PEAKS 11 header:** `Source File,Scan,Peptide,Tag length,ALC (%),Length,m/z,z,RT,Area,Mass,ppm,PTM,local confidence (%),mode,tag(>=0.0%),Feature Id`
- **[O] PEAKS 12** adds `Deep Novo Score (%)`.
- `local confidence (%)` is space-separated per-residue LC; `mode` is CID/HCD/ETHCD etc.
- The report parameters "Denovo Only Tag Sharing" and "Fully Matched" control coverage display, not this file [D].

### 7.7 `db.ptmProfile.csv` (PTM Profiling CSV [D S6 §6.4.3.3; S9])

- Grouped by protein, then PTM type, then site.
- **Columns:** Protein (accession), Protein Position (e.g. `S123`), PTM, −10lgP (best modified peptide), AScore (best), Ion Intensity (best), `<sample> Modified` area, `<sample> Unmodified` area, **Modified Peptide Sequence Window** (±10 aa, padded), Sample Best AScore.
- With "Fully digested": Peptide and Peptide Position.
- Only confident PTMs are included, using the Summary-page AScore and ion-intensity filters. Fixed PTMs are excluded.

### 7.8 `db.modifiedPeptide.csv` (Modified Peptide view, updated PEAKS 12 [D S6 §6.6])

- **[I]** Columns mirror the Peptide table, pairing each modified peptide with its unmodified counterpart. PTMs are flagged confident or unconfident by the AScore/ion-intensity filters (green or red in the UI).

---

## 8. CLI, batch mode, API and project format

- **[D] PEAKS Studio 11/12: GUI only (Windows 10+ 64-bit).**
  - "Batch submission" in the Report step "runs a separate analysis with the same workflow parameters for each sample".
  - An "Export" checkbox auto-exports all CSVs (S6 §1.3).
  - Workflows are saved and reused. "Modify Analysis" reruns with changed parameters (S6 §1.6).
- **[D] PEAKS Studio 13.5 (June 12, 2026)** adds "Command Line Interface (CLI) support enabling project submission" and "Integration of CLI into PEAKS Studio allows for automated analyses pipelines" (S12; S13). No syntax is published.
- **[D] PEAKS Online** has a Web client and a "Client CLI" for pipelines (bioinfor.com/peaks-online). It is a distributed Akka/JVM master-worker system (S11).
  - Archived Online projects can be imported into Studio 11/12 but cannot be re-analysed "due to different project structures" (S6 §1.7).
- **[D] PEAKS CMD (2015, PEAKS 7/8 era), Windows and Linux (S16):**
  - `peakscmd -c template.xml (-i file | -s samplelist.txt | -d dir_of_sample_lists)`.
  - The XML template is built with a GUI `TemplateEditor` and covers loading, de novo or DB search (+PTM/SPIDER), filters and exports.
  - Defaults: `$HOME/PEAKS_CMD_TEMPLATE/CMDTemplate.xml`, projects in `$HOME/PEAKS_CMD_PROJECT`, exports in `$HOME/PeaksCMD_export/<project>`.
  - Accepts PKL, MGF, mzXML and mzML.
  - Also ships `ProjectConverter`, `MassCalculator` and `PerformanceConfiguration`.
- **[D] "PEAKS Desktop" is not a product.** BSI sells PEAKS Studio under a "Desktop License" (24 threads) or a "Workstation License" (48 threads) (S27). The product family is PEAKS Studio, Online, Client, Viewer, AB and GlycanFinder.
- **[D] Project on disk.**
  - "A project is saved as a directory that contains multiple files, including the compressed spectral data and the analysis results. To share a PEAKS project, you can copy the whole project directory" (S8). The default is `[USER_HOME]/PeaksProjects`.
  - Older versions used an embedded Derby DB (FAQ: delete `defaultproject`, `defaultTaskDB`, `derbyServer` folders).
  - 11+ uses a "Default Project Directory" setting. Configured databases and saved workflows live in the install directory's `db` and `es` folders ("copy over the db and es folders" when upgrading, S6 §1.4.2).
  - A JSON-style `peaks.conf` exists with hidden keys `"show-debug-panel"`, `"export-decoy"` and `"detach-service"` (S18; S24).
  - **[I]** The internal file layout and serialisation are undocumented and proprietary. No public API or SDK exists; interoperability is via CSV, pepXML 1.22 and mzIdentML 1.2.0 exports (S6 §6.9.3).

---

## 9. Defaults and scoring details

- **[D] Peptide length:** default **6-45** (S6 §6.1; S9).
- **[D] Parameters.**
  - Variable PTMs per peptide: <4 recommended.
  - Missed cleavages: user-set.
  - Digest mode: Specific / Semispecific / Unspecific / No digestion.
  - Enzyme "None" is recommended for immunopeptides.
- **[D] Defaults and thresholds.**
  - De novo ALC 50%.
  - SPIDER ALC > 15%.
  - Protein −10lgP ≥ 15.
  - Without FDR, a starting peptide threshold of −10lgP 20 (p = 0.01) is suggested (S8).
- **[D] Tolerances.**
  - Precursor in ppm or Da. Fragment is in Da in 11/12 de novo/DB dialogs; "Units can be selected" for DeepNovo.
  - Reported defaults in studies: 10 ppm / 0.02 Da for Q Exactive and 15 ppm / 0.03 Da for timsTOF (S18); 15 ppm / 0.8 Da for LTQ (S2).
  - "A fragment ion is considered matched if the calculated m/z is within the specified error tolerance" (S6).
  - For annotation, a match also needs relative intensity ≥ 2%, with 5% and 10% options (S6 §6.5.5).
  - **[I]** At 0.5 Da the scoring bins f_N and f_C by tolerance: an ion contributes if any peak lies within ±tol. There is no evidence of intensity-weighted mass-error scoring for low-res data. Setting instrument type to IT/LIT selects different internal fragment weights and noise models, and ion-trap setups may use **average precursor mass** ("For ion-trap instruments, it is usually beneficial to allow PEAKS DB search to use an average mass", S8).
- **[D] Charge.** Precursor charge comes from features or refinement. Older refinement could correct charge within min/max bounds. Charge is a scoring feature. Fragment charge states are not documented.
  - **[I]** Imitate with fragment charges up to min(z−1, 2) for CID/HCD.
- **[D] Ion types for DB scoring:** see §1.1 (CID a, b, c, y, z plus neutral losses; ETD a, b, c, c-H, y, z, z+H). EThcD is supported (S8).

---

## 10. Published benchmarks involving PEAKS

| Study | Comparison | Result |
|---|---|---|
| S2 (2012) | PEAKS DB vs Mascot, Mascot+Percolator, SEQUEST (1% FDR PSMs) | CID: 10,668 vs 7,515 / 9,969 / 8,236. ETD: 3,652 vs 2,398 / 2,702 / 2,233 (+52% vs Mascot) |
| S26 (MCP 2021) | PEAKS vs MaxQuant, Comet, MS-GF+ (HLA-I) | PEAKS reports "~30-70%" more peptides at <1% FDR, with more L/I/P/V (hydrophobic) bias |
| S11 (Nat Commun 2022) | PEAKS Online vs MaxQuant(+Prosit), MS²Rescore (HLA DDA) | 1.7-4.1× MaxQuant; rescoring 1.0-1.4× MaxQuant+Prosit; MS²Rescore +17% on one of 3 sets. DIA: 5-30% more precursors than DIA-NN/Spectronaut |
| S17 (MCP 2025, timsTOF DDA-PASEF, JY cells) | PEAKS X Pro, PEAKS 11 Online (+DL boost), FragPipe 18 (+MSBooster), MaxQuant 2.4.4, MHCquant 2.6 | PEAKS 11 highest: +6-11% over PEAKS X Pro and FragPipe. MHCquant −6-9%; MaxQuant 50-60% fewer. PEAKS lost IDs as the database grew (UniProt+transcriptome+de novo), unlike FragPipe and MHCquant. Entrapment showed unbiased FDR for all tools. PEAKS reported the most I/L-variant "SAVs" |
| S18 (JPR 2025, Listeria HLA) | PEAKS 12 (+DL boost), MSFragger 4.1, Comet 2023.01, Sage 0.14.7, each ± MS²Rescore | MS²Rescore added 30-89% peptides to Sage, Comet and MSFragger. PEAKS was already high "prior MS²Rescore" because of the DL boost. Combined 4 engines: +27% average. 57-69% of strong binders found by all four |
| S19 (Appl Sci 2025) | PEAKS DL boost on/off | Protein IDs +22.9% CHO, +23.2% E. coli |
| S14 (BSI) | DL boost on/off | ">10%" more IDs |
| S20 (ContraNovo, 9-species V1) | De novo peptide precision | PEAKS 0.322 average peptide precision (AA precision 0.663). DeepNovo 0.376, PointNovo 0.439, Casanovo 0.448, CasanovoV2 0.539, ContraNovo 0.621 |
| S21 (PNAS 2017) | DeepNovo vs PEAKS, Novor, PepNovo | DeepNovo 38.1-64.0% higher peptide-level and 7.7-22.9% higher AA-level accuracy than the best competitor (PEAKS or Novor were usually second) |

No peer-reviewed head-to-head of PEAKS DB against Sage on standard tryptic DDA was found. S18 is the closest, and it compares Sage in HLA data after rescoring. ProteoBench (bioRxiv 10.64898/2025.12.09.692895) lists PEAKS submissions: in the Astral module, PEAKS quantified about 125,026 precursors versus DIA-NN's 125,114 at k = 3, per search snippet. I could not verify the full text (HTTP 429).

---

## 11. Checklist for imitating the PEAKS DDA workflow

1. **Refine.** Detect MS1 isotope features. Assign each MS2 scan to every feature in its isolation window (chimera mode). Correct the precursor to the feature's monoisotopic m/z and charge. Create fake features (Feature Id 0) for unassigned MS2. Record Area as the feature XIC area.
2. **De novo.** Run DP over mass with f_N/f_C odds and the top 10,000, then rescore and take the top 5. Compute per-residue LC by candidate agreement, and ALC as mean LC. Build tags by collapsing LC < 30% residues into mass segments.
3. **Tag-filter the database.** Compute CAA by mass-block alignment. Shortlist the top 7,000 proteins (configurable).
4. **Decoy fusion.** For each protein, build `target + 'J' + shuffle(target)`. Digest with J as an always-cleave site. Label a hit decoy if it lies in the second half.
5. **Candidates and scoring.** Take the top 512 per spectrum by quick score. Rerank with an LDA over the normalised ion score ((s−μ10)/σ150) plus the 9 features (or about 18 in the newer engine), including CAA and the protein feature. Convert to p against a null and report −10lgP.
6. **Optional DL boost.** Add ΔRT, MS2 correlation, Δ1/k0, CAA and Tag Length. Semi-supervised rescoring on fused decoys; report the rescored −10lgP.
7. **FDR.** Use D/T at the PSM level (one per spectrum-precursor) and at the peptide level (best PSM per peptidoform, with sequences deduplicated ignoring I/L as PEAKS counts #Peptides).
8. **Proteins.** Group by domination at T = 15. Choose Top by most unique peptides. Score = Σ s_k/k over non-redundant peptides. Protein-group FDR is based on the top-hit decoy status. Defaults: −10lgP ≥ 15 and ≥ 1 unique [I].
9. **PTM and SPIDER passes** on unexplained spectra with ALC above threshold.
10. **Exports.** Write the column sets in §7, using PEAKS conventions: `(+57.02)` mods, `#DECOY#` prefix, `-10LgP` capitalisation (12+), `Scan` as spectrum index + 1 for mzML, and `K.PEPTIDE.R` flanks in protein-peptides.
