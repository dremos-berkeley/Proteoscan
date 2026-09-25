# How PEAKS Studio preprocesses DDA spectra: evidence from its MGF export

Run: DrewR_Adam_91526_DDA_OT_4 (S4, Jurkat), with a check against DrewR_Adam_91526_DDA_OT_1 (S1, panel).
Everything below is inferred from the data only (PEAKS MGF, features.csv, db.psms.csv, denovoOnly.csv,
ThermoRawFileParser 2.0 mzML); no binaries were inspected. Numbers are S4 unless marked S1.

Terms:
- "mzML precursor": the ThermoRawFileParser selected-ion m/z. It equals the Thermo trailer "Monoisotopic M/Z" in 99.99% of scans.
- "Isolation target": the m/z in the filter string, i.e. the window centre.
- "Master MS1": the MS1 named by precursor/@spectrumRef. On this Lumos it is always the MS1 before the most recent one.
- "Feature-derived" vs "featureless": an entry is feature-derived if its RTINSECONDS differs from the scan RT by more than 0.01 s. This proxy agrees 99.7% with the known Feature Id on 109,227 entries that have a PSM or de novo row.

## Summary
1. Each MGF entry is one (scan, precursor) pair, of two kinds.
   - Feature-derived (39.1%; S1 50.8%). PEPMASS, CHARGE and RT are the feature's monoisotopic m/z, charge (1-9) and apex RT, exactly as in features.csv. A feature is linked when the scan RT is within [RT Begin - 0.02, RT End + 0.02] min and some isotope of the feature is inside the window. The monoisotopic peak is the one inside in 90-92% of cases.
   - Featureless (60.9%; S1 49.2%). CHARGE is the instrument charge (99.99%) and RT is the scan RT (within 10 ms). PEPMASS comes from PEAKS's own monoisotope correction in the master MS1. It is usually an exact centroid there (91%; S1 77%), often shifted to a lower isotope. The entry is emitted unless a feature entry in the same scan has the same charge within about 12 ppm (zero exceptions at 10 ppm; smallest gap seen 12.0 ppm in both samples).
2. Entries per scan: 72.6% have 1, 21.2% have 2, 4.9% have 3, 1.3% have 4 or more (max 7). S1: 66.7 / 23.8 / 6.8 / 2.7% (max 11). All 194,508 MS2 scans are present.
3. Peak lists are not the mzML centroids. About 0.71x the peak count (S1 0.75x), no entry has an identical count. About 77% of MGF peaks are verbatim centroids (m/z rounded to 0.001, same intensity). The rest come from merging neighbours <= 0.5 Th apart (summed intensity); folding singly-charged isotopes (+1.003 +/- 0.2 Th) into the lighter peak; ~1% charge-reduced peaks; removing most of the precursor peak. Processing depends on the precursor charge. No top-N filter, no normalisation, no truncation at the precursor mass. A simple re-implementation matched 63% of peaks, versus 73% for raw centroids: the exact algorithm was not recovered.
4. PSM table columns: m/z equals PEPMASS (100%), RT is RTINSECONDS/60, raw RT is the scan start time. raw m/z is the mzML precursor, but filled only for the one representative PSM per peptide (3,472 non-zero rows = 3,472 peptides).

## Q1. MGF entries per raw MS2 scan
| entries/scan | S4 scans | S4 frac | S1 scans | S1 frac |
|---|---|---|---|---|
| 1 | 141,245 | 0.7262 | 62,144 | 0.6673 |
| 2 | 41,187 | 0.2117 | 22,159 | 0.2379 |
| 3 | 9,607 | 0.0494 | 6,313 | 0.0678 |
| 4 | 1,997 | 0.0103 | 1,785 | 0.0192 |
| 5 | 399 | 0.0021 | 483 | 0.0052 |
| 6 | 66 | 0.0003 | 143 | 0.0015 |
| 7 | 7 | <0.0001 | 52 | 0.0006 |
| 8-11 | 0 | 0 | 46 | 0.0005 |
| total | 194,508 scans / 262,868 entries | | 93,125 / 136,562 | |

Scan composition in S4. Rows: number of featureless entries; columns: number of feature-derived entries.

| featureless \ feature | 0 | 1 | 2 | 3 | 4 | 5 | 6+ | All |
|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 24,489 | 7,729 | 1,788 | 358 | 63 | 17 | 34,444 |
| 1 | 116,756 | 33,435 | 7,814 | 1,636 | 336 | 51 | 5 | 160,033 |
| 2 | 23 | 5 | 3 | 0 | 0 | 0 | 0 | 31 |

- A scan has 0 or 1 featureless entry plus 0-7 feature entries (S1: up to 11).
- The 31 scans with 2 featureless entries are almost certainly feature entries whose apex equals the scan RT.

Pairwise relations between entries of the same scan:

| relation | S4 pairs | S1 pairs |
|---|---|---|
| unrelated species, different z | 49,572 | 38,478 |
| unrelated species, same z | 20,398 | 8,884 |
| isotopes of each other (same z, delta = k*1.00335/z) | 9,533 | 9,805 |
| same m/z (<0.01), different z | 5,844 | 3,301 |
| same m/z, same z (duplicate features) | 1,770 | 830 |

73% of multi-entry scans (S1 74%) contain more than one charge.

Each entry vs the mzML precursor (S4):

| relation | multi/feature | multi/featureless | single/feature | single/featureless | All |
|---|---|---|---|---|---|
| identical (<2 ppm), same z | 7,986 | 27,799 | 15,919 | 74,662 | 126,366 |
| same z, 2-20 ppm off | 1,208 | 1,064 | 1,616 | 1,692 | 5,580 |
| same z, isotope shift, abs(k) 1-3 | 9,373 | 13,195 | 6,486 | 38,322 | 67,376 |
| same m/z, charge reassigned | 5,143 | 5 | 0 | 0 | 5,148 |
| other species, m/z in window | 44,172 | 1,019 | 397 | 1,746 | 47,334 |
| other, m/z outside window (an isotope is inside) | 10,402 | 257 | 71 | 334 | 11,064 |

## Q2. Single-entry scans
PEPMASS class vs the mzML precursor:

| class | single/feature (24,489) | single/featureless (116,756) |
|---|---|---|
| identical <2 ppm | 15,919 (65.0%) | 74,662 (63.9%) |
| 2-20 ppm, same z | 1,616 | 1,692 |
| isotope shift k = -3 / -2 / -1 | 5 / 586 / 4,124 | 180 / 5,794 / 26,569 |
| isotope shift k = +1 / +2 / +3 | 1,433 / 255 / 83 | 4,830 / 798 / 151 |
| charge reassignment | 0 | 0 |
| other | 468 | 2,080 |

- Positive k means the instrument's trailer mono is lower than PEAKS's choice. The trailer mono differs from the isolation target in 21% of S4 scans (S1 43%).
- Relative to the isolation target, PEAKS never goes up: k > 0 in 8 of 160,095 featureless entries.

Featureless PEPMASS vs the instrument's two m/z values (single-entry scans, 2 ppm):

| sample | = mono & = target | = mono only | = target only | neither |
|---|---|---|---|---|
| S4 | 66,204 | 8,458 | 6,558 | 35,536 |
| S1 | 16,608 | 5,728 | 3,141 | 17,631 |

Where featureless PEPMASS comes from (S4 / S1):
- A. Exact centroid of the master MS1: 91.2% / 76.8%.
- B. Arithmetic: the centroid nearest the target plus k*1.0034/z, used when the monoisotopic peak is not observed. 8.2% / 19.4%, mostly above 2 kDa. Implied spacing median 1.00341.
- C. Instrument trailer mono: 0.6% / 3.8%.

Isotope shift relative to the target peak (k_filt) by neutral mass, S4, row fractions:

| M (Da) | k = -4 | -3 | -2 | -1 | 0 | n |
|---|---|---|---|---|---|---|
| <1000 | 0 | 0 | .084 | .149 | .768 | 42,101 |
| 1000-1500 | 0 | 0 | .034 | .320 | .647 | 53,576 |
| 1500-2000 | 0 | 0 | .035 | .374 | .591 | 24,949 |
| 2000-2500 | 0 | 0 | .052 | .460 | .488 | 14,467 |
| 2500-3000 | 0 | 0 | .199 | .319 | .482 | 12,187 |
| 3000-4000 | 0 | .004 | .328 | .246 | .422 | 7,977 |
| 4000-6000 | .018 | .170 | .274 | .053 | .484 | 4,461 |

PEAKS shifts more often than the instrument, and the extra shifts look partly wrong:
- Where the instrument said "the isolated peak is monoisotopic", PEAKS still shifted 30,977 scans by -1 and 7,540 by -2.
- PSM rate per featureless entry falls with the shift: 1.57% (k = 0), 1.08% (-1), 0.26% (-2), 0% (-3).
- The shift is not a simple intensity-ratio cut. Median I(-1)/I(0) is 0.49 when PEAKS kept k = 0 and 0.68 when it shifted to -1, with heavy overlap.

CHARGE:
- Single-entry scans: identical to the mzML charge in 141,245 of 141,245 (S1 likewise). Featureless entries: 99.99% identical.
- All charge changes come from features. E.g. 40,715 z = 1 entries sit on instrument z = 2 (33,762), 3 (6,043) and 4-7 (910).
- The mzML has no charge 0/unknown scans (instrument charge filter 2-7: 2: 145,441; 3: 41,717; 4: 2,573; 5: 2,607; 6: 1,319; 7: 851), so PEAKS's handling of unknown charge cannot be observed here.

## Q3. MGF peak lists vs mzML centroids
Sample: 2,676 entries from 2,015 scans (S1: 1,507 entries from 967 scans).

| metric, median [5-95%] | S4 | S1 |
|---|---|---|
| mzML peaks | 719 [315-1200] | 676 [314-1723] |
| MGF peaks | 509 [246-720] | 516 [274-814] |
| count ratio | 0.71 [0.55-0.86] | 0.75 [0.46-0.93] |
| identical count | 0.07% | 0% |
| MGF peaks within 0.01 Th of an mzML peak | 0.82 | 0.87 |
| mzML peaks with an MGF peak within 0.01 Th | 0.58 | 0.65 |
| intensity of matched peaks, MGF/mzML | 1.000 (85% exact; rest up to ~2x from folded isotopes) | 1.000 (83% exact) |
| TIC ratio, all entries | 0.95 [0.83-1.00], never >1 | 0.92 [0.72-1.00] |
| base-peak ratio | 0.82 | 0.79 |

Origin of MGF peaks: verbatim centroid (m/z rounded to 0.001, same intensity) 76.8% (S1 82.0%); charge-reduced copy z*m - (z-1)*1.00728: 1.4% (S1 0.2%); other (merged or re-centred): 21.9% (S1 17.8%).

Fate of mzML centroids (S4): 52% kept verbatim, 19% folded into the lighter peak as its +1.00335 isotope, 29% absent (merged away or removed as precursor).

Established:
- Merging. Neighbouring centroids 0.50 Th apart or closer are never both kept (0 of 187,684 pairs); above 0.5 Th, 9-33% are. The merged peak has the summed intensity; for simple pairs its m/z is the intensity-weighted mean, e.g. 175.152 (18.7) + 175.521 (20.4) -> 175.344 (39.1). Dense clusters are re-centred in a way a weighted mean does not explain (looks like profile re-picking). S4's ion trap ran at turbo rate ("t d"), so close centroids are common; S1 ran at rapid rate ("r d") where they are rare.
- z = 1 deisotoping. When the heavier peak sits 0.80-1.20 Th above a kept peak, it is folded into the lighter one (intensities added) in 27-30% of cases; outside that range 0%. Folding depends on heavy/light intensity ratio: 64% when ratio < 0.25, 7% at 1-1.5, 0% above 4.
- Precursor removal. Median 87% of the lost intensity (S1 93%) lies within target +/- (half-width + 0.5 Th). Median 72% of that region's intensity is removed (S1 89%); at most 3-4 peaks remain there. The mzML base peak sits in that region in 63% of entries (S1 72%).
- Charge dependence. Entries of the same scan with the same charge have identical peak lists in 88% of cases (S1 85%); with different charges 9% (S1 34%).
- Not applied: top-N filtering (up to 955 peaks kept, down to 0.03% of base peak); intensity scaling; removal of peaks above the singly-charged precursor mass (83k such peaks kept on z = 1 entries).
- Formatting: peak m/z rounded to 0.001, stored as float32, printed with 5 decimals (129.10699 is float32(129.107)). Intensity is the mzML float32 printed with 5 decimals. PEPMASS printed %.5f with trailing zeros stripped, no intensity field. Peaks sorted by m/z.
- Simulation (simulate_peaks.py): chain-merging at <= 0.5 Th plus z = 1 deisotoping (+/-0.2 Th) reproduced only 63% of MGF peaks exactly, versus 73% for raw centroids.

## Q4. Features
PSMs:
- 7,522 of 7,524 PSMs match an MGF entry by (Scan = RAWSCANS, m/z = PEPMASS); the 2 misses are duplicate rows. S1: 834 of 834.
- All 5,413 PSMs with a Feature Id (S1: 530) are found in features.csv; their m/z, z and RT equal the feature's in 100% of cases.
- PSM RT equals RTINSECONDS/60 in 99.9%; raw RT equals the scan RT.
- 2,111 PSMs (28%) have Feature Id 0.

| PSM origin (S4) | PSMs | MGF entries | PSM rate |
|---|---|---|---|
| feature-derived | 5,417 | 102,773 | 5.2% |
| featureless, k = 0 | 1,581 | 100,187 | 1.6% |
| featureless, shifted | 524 | 59,908 | 0.9% |

Entries with a detected feature: 102,773 (39.1%) vs 160,095 without (60.9%). S1: 50.8% / 49.2%. In S4, 60% of scans have no feature entry at all.

features.csv is only a subset of the feature map: it holds only features with a PSM or de novo result >= 50 ALC (23,540 rows, Feature Ids up to 80,874). The MGF references 50,296 distinct features. 46,712 feature-derived entries (S1: 23,011) could be tied to a features.csv row by (m/z, z, apex).

For those tied entries:
- RT coverage: scan RT inside [RT Begin, RT End] in 88.5% (S1 93.0%), within 0.02 min of it in the rest; furthest outside 0.0201 min.
- Window position: monoisotopic m/z inside the window in 89.6% (S1 92.1%); otherwise only isotope 1, 2 or 3 is inside: 8.7 / 1.5 / 0.2%.
- RT reported: RTINSECONDS equals features.csv RT x 60 in 98.3% (S1 99.0%), equals RT Begin in only 1%.
- Feature width: median 0.36 min, 95th percentile 6.6 min, capped at 10.0 min (so RTINSECONDS can be up to +/-301 s from the scan RT).

Which candidate features get linked? (features.csv subset; precision is a lower bound.)

| rule | S4 precision | S4 recall | S1 precision | S1 recall |
|---|---|---|---|---|
| A: feature = mzML precursor (10 ppm, same z) | 0.936 | 0.387 | 0.958 | 0.518 |
| B: A, or feature is an isotope (abs(k) <= 3) of the precursor | 0.730 | 0.558 | 0.798 | 0.654 |
| C: any isotope of the feature in the window | 0.231 | 1.000 | 0.292 | 1.000 |
| D: monoisotopic m/z in the window | 0.515 | 0.899 | 0.545 | 0.923 |
| G: B, or (D and window-current share >= 5% and width <= 0.55 min) | 0.678 | 0.740 | 0.746 | 0.765 |

Co-isolated features that are other species get linked only 12% of the time. Their share of the window ion current barely matters (~17% linked for any share > 0; 3% when share is 0). Width matters: 31% for features under 0.22 min wide, 5% above 1.3 min. Position matters: ~39% when the monoisotopic m/z is in the window, under 1% when it is more than 1.5 Th below. The remaining selection criterion was not identified.

When does a featureless entry sit alongside feature entries? (S4, 77,729 scans with >= 1 feature entry)

| best feature relation to the instrument precursor | no featureless | + featureless |
|---|---|---|
| same m/z and z | 26,844 | 2,273 (in 83% of these PEAKS shifted its own precursor to another isotope) |
| same m/z, different z | 24 | 4,930 |
| isotope, same z | 6,880 | 3,940 |
| none | 696 | 32,142 |

No featureless entry has a same-z feature entry within 10 ppm in the same scan; minimum gap 12.0 ppm (S4 and S1).

## Q5. Isolation window
- Constant and symmetric: +/-0.65 Th (1.3 Th) in every S4 scan and +/-0.8 Th (1.6 Th) in every S1 scan. Centred on the filter-string target, not on the trailer mono.
- Featureless PEPMASS is always the target peak or a lower isotope of it (-2.03 to +0.52 Th).
- Feature PEPMASS ranges from -2.58 to +0.65 Th of the target; the monoisotopic m/z can fall below the window but some isotope is always inside.

First isotope j of each entry's envelope inside target +/- h (S4):

| group | j = 0 | 1 | 2 | 3 | 4 | none |
|---|---|---|---|---|---|---|
| multi/feature | 65,638 | 10,781 | 1,497 | 271 | 94 | 3 |
| multi/featureless | 40,131 | 3,182 | 24 | 1 | 1 | 0 |
| single/feature | 23,241 | 1,234 | 14 | 0 | 0 | 0 |
| single/featureless | 107,892 | 8,428 | 436 | 0 | 0 | 0 |

## Q6. RTINSECONDS vs scan start time
- Featureless: RTINSECONDS - scan RT ranges -0.010 to +0.010 s (median -0.0035 s): the scan RT with ms-level offset/rounding.
- Feature-derived: difference ranges -301 to +301 s (median -0.57 s, IQR -8.9 to +4.4 s): the feature apex RT.
- The MGF is sorted by RTINSECONDS, not by scan number. SCANS is the 1-based entry index and equals the TITLE suffix. SPECGROUPID is <one uuid>_<SCANS>.

## Q7. Other observations
- No PEPMASS intensity field; CHARGE is never missing.
- Charge 1: 40,715 entries (15.5%), all feature-derived, never the only entry of a scan. 4,485 have the same m/z as the instrument's z >= 2 precursor. PSM rate 0.48% (S1 0.34%) vs 3.7% for z = 2.
- Charges 8-9: 117 entries (S1 445), all feature-derived, 0 PSMs.
- S1 charges 4-7: 36,562 entries with 1 PSM.
- Master MS1 matters: the featureless PEPMASS is an exact centroid of the spectrumRef MS1, not the latest MS1 (2% match in the latest vs 91% in the referenced one).
- Duplicate features: same m/z and z linked twice to one scan in 1,770 pairs (S4).
- Trailer mono field: pyopenms drops the trailer mono userParam; extract_trailer.sh recovers it.

## Rules to imitate (S4 values; S1 in brackets)
1. Per MS2 scan, read: scan RT; isolation target t (filter string or mzML isolation window target); half-width h = 0.65 [0.8]; instrument charge z_i and trailer mono; the master MS1 from precursor/@spectrumRef (the MS1 before the latest).
2. Feature entries. Candidates are features with RT_begin - 0.02 min <= RT <= RT_end + 0.02 min and some isotope m + j*1.00335/z (j = 0-4) inside [t - h, t + h]. To approximate PEAKS, keep candidates whose monoisotopic m/z is inside the window (recall 0.90, precision >= 0.52). Always keep one that equals the instrument precursor, or an isotope of it, at the same z. For more precision, drop co-isolated features wider than ~0.55 min or with < 5% of the window current (precision 0.68, recall 0.74). Emit one entry per linked feature: PEPMASS = feature m/z, CHARGE = feature z, RT = apex. Do not deduplicate features against each other.
3. Featureless precursor. Take the master-MS1 centroid nearest t (always within 10 ppm). With z = z_i, step down by 1.00335/z while the isotope pattern supports it. PEAKS's outcomes: k = 0 in 63%, -1 in 28%, -2 in 8%, -3/-4 in 0.7%. Shifts of -3 appear only above ~3 kDa and -4 above ~4 kDa. It never shifts upward. A more conservative averagine test is probably better, since shifted entries identify less often. PEPMASS is the observed centroid of the chosen isotope (91% [77%]); if absent, centroid(t) - k*1.00335/z (8% [19%]); otherwise the trailer mono. CHARGE = z_i, never changed. RT = scan RT.
4. Emitting the featureless entry. Emit it unless a feature entry in the scan has the same z and m/z within 12 ppm. Emit it even when a feature has the same m/z but a different z. A scan with no feature gets exactly one entry. Expect ~1.35 entries per scan, 39% feature-derived.
5. Output order. Sort all entries by RT, number them SCANS = 1..N, keep RAWSCANS as the native scan number.
6. MS2 peaks (approximate only). Round m/z to 0.001; no normalisation, no top-N. Remove the precursor region (window +/- 0.5 Th, leaving ~3 peaks or fewer). Merge centroids <= 0.5 Th apart (summed intensity, weighted m/z). Fold +1.00335 +/- 0.20 Th isotopes into the lighter peak when the heavier one is weaker. Allow charge-aware deconvolution up to the precursor charge (1-2% charge-reduced peaks). Expect ~0.7x the peak count and >= 95% of the TIC. Because the exact algorithm was not recovered, the engine's own deisotoping on raw centroids plus precursor-region removal is a reasonable substitute.
7. Charge 1 and 8-9 precursors yield almost nothing (0.3-0.5% and 0% PSM rates; 194 of 7,524 PSMs are z = 1). They can be omitted with little loss.

## Scripts
All in this directory; run with `conda activate openms` and `python <script> s4|s1`.

| script | purpose |
|---|---|
| parse_mgf.py | streamed MGF -> header table plus peak lists for raw scan % 97 == 0 |
| parse_mzml.py | pyopenms MS2 metadata, sampled MS2 peaks, all MS1 peaks |
| extract_trailer.sh | recovers the trailer "Monoisotopic M/Z" |
| common.py | loaders, joins, markdown helper |
| tables.py | Q1/Q2/Q4/Q6/Q7 tables -> tables_<s>.md |
| multi_entry.py | relations within a scan, window table |
| featureless_ms1.py, featureless_summary.py | featureless PEPMASS source |
| featureless_rule.py | featureless vs feature entries, 12 ppm check |
| feature_match.py, feature_rule.py, feature_intensity.py, rule_eval.py | feature-linking analyses |
| charge1.py | charge-1 and charge-8/9 entries |
| peaks_compare.py, peaks_explain.py, merge_gap.py, precursor_removal.py, intensity_loss.py, entry_peak_diff.py, simulate_peaks.py | Q3 peak-list analyses |

Cached .pkl intermediates take about 1.2 GB. The a2/ subdirectory was created by another process during this analysis and was left untouched.
