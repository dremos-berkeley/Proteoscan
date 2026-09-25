# Brief for the Windows-side Claude session (PEAKS desktop)

Paste this as the first message of a Code session running locally on the Windows machine
where PEAKS is installed (Computer use enabled in Settings > Desktop app).

---

PEAKS Studio is installed on this Windows machine. A Linux-side session is building an
open-source replica of its DDA identification pipeline and needs reference material from
PEAKS. Do not decompile or reverse-engineer the PEAKS binaries; only use documented
behavior, its GUI, its CLI if it has one, and its exports. Work through these tasks in
order, write everything into C:\peaks-parity-exports\ (create it), and keep a
NOTES.md there describing exactly what you did and what each file is.

## 1. Inventory
- Find the PEAKS installation directory, product name and exact version (Help > About).
- List the executables, any command-line tools, batch/automation options, config files,
  and where projects are stored on disk. Copy any user manual / release notes PDFs you
  find into the export folder.
- Screenshot the main search parameter dialogs (database search, de novo, DL boost,
  FDR / filter settings, export options) into the export folder.

## 2. Existing project "DrewR_Adam_91526_DDA_OT" (4 runs)
If the PEAKS project for the four DrewR_Adam_91526_DDA_OT_{1..4}.raw runs exists on this
machine, open it; otherwise create it with EXACTLY these parameters (they are the
facility's original settings and must be copied verbatim, even the fixed carbamidomethyl):
  Enzyme Trypsin, semi-specific, max 2 missed cleavages; precursor 15 ppm, fragment 0.5 Da,
  monoisotopic; fixed Carbamidomethylation (C); variable Oxidation (M), Deamidation (NQ);
  max 3 variable PTM per peptide; peptide length 6-45; database = the same "Homo_sapiens"
  database with 20597 entries; Deep Learning Boost ON; FDR estimation ON; PEAKS PTM OFF;
  SPIDER OFF.
Then produce these exports, each into its own subfolder:
  a. dlboost_on_default:   default filters (1% peptide FDR). All CSVs, mzIdentML, MGF,
                            features, de novo only, PTM profile, summary HTML.
  b. dlboost_on_unfiltered: same result, but set the peptide -10lgP filter to 0 (or FDR
                            100%) and the protein -10lgP filter to 0, unique peptides 0,
                            de novo ALC 0, so that ALL PSMs including low-scoring ones are
                            exported. If decoy hits can be shown/exported, include them.
  c. dlboost_off_default:  re-run the identical search with Deep Learning Boost OFF,
                            default filters.
  d. dlboost_off_unfiltered: as (b) for the DL-boost-off run.
  e. Export the exact FASTA database used (Configuration > Database, or copy the file it
     points to) with its md5.
  f. Export the de novo results (all spectra, ALC >= 0) as CSV for each run.
  g. Export the "all PSMs / peptide candidates" if PEAKS offers a lower-rank candidate
     export; note if it does not.

## 3. Public HeLa benchmark PXD004940 (Lumos)
Download these 9 raw files from PRIDE (ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2017/03/PXD004940/):
  151118_Q_Hela_03_0{1,2,3}_50cm_HCD_IT_EasySpray_90min_500ng.raw
  151118_Q_Hela_05_0{1,2,3}_50cm_HCD_IT_EasySpray_90min_1000ng.raw
  151019_Q_Hela_03_0{1,2,3}_1ug_50cm_HCD_FT_EasySpray_90min.raw
Create a PEAKS project with one sample per file. Search with the same parameters as in
section 2 (CAM fixed is chemically correct for these; keep 0.5 Da fragment tolerance for
the HCD_IT files; for the three HCD_FT files make a second project with fragment
tolerance 0.02 Da). Export (a) and (b) as above, plus the FASTA.

## 4. Hand-off
Zip C:\peaks-parity-exports\ and copy it to the Linux server at
/drive_a/drew/vendor-exports/ (scp/rsync/WinSCP; Drew has the credentials), or tell Drew
it is ready.
