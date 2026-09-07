# Steps that need a browser (cannot be done from the command line)

1. **Zenodo integration** — **DONE 2026-09-08.** The repository was made
   public, the Zenodo GitHub switch for `dmai287/csa-audit` was flipped, and
   release `v0.1.0` was cut, which minted the DOIs below. The concept DOI is
   recorded in `CITATION.cff` and in the manuscript's data availability
   statement.

   - Concept DOI (always the latest version): `10.5281/zenodo.22648224`
   - Version DOI for `v0.1.0`: `10.5281/zenodo.22648225`

   Every later release is archived automatically; cut a fresh tagged release
   before submission so the archived code matches the reported analysis.

2. **fastMRI data-use agreement.** Request access at
   https://fastmri.med.nyu.edu/ (brain multi-coil). The download links arrive
   by e-mail and expire; store the data outside this repository (`data/` is
   git-ignored). Status 2026-09-08: links received (valid to 2026-12-07); a
   resumable downloader runs on the external T9 drive at
   `/Volumes/T9/fastMRI_brain/` (val, test_full and train batch 0 first; the
   full set is 1.51 TB and exceeds the 1 TB drive, see its README.txt).
   fastMRI+ annotations are public at
   https://github.com/microsoft/fastmri-plus (CSV files).

3. **OSF pre-registration.** Create the registration at https://osf.io from
   the analysis plan (paper Section 4.8) and `configs/thresholds.yaml` once the
   thresholds are chosen; save the PDF and its SHA-256 under `prereg/`.

4. **Ethics** for the optional reader spot-check, through the RMIT HREC portal.
