# Steps that need a browser (cannot be done from the command line)

1. **Zenodo integration** (needed for the concept DOI in the paper's data
   availability statement). Log in at https://zenodo.org with the GitHub
   account `dmai287`, open https://zenodo.org/account/settings/github/, flip
   the switch for `dmai287/csa-audit`. The repository must be public at the
   time of the first GitHub release; Zenodo mints the DOI from that release.
   Then paste the concept DOI into `CITATION.cff` (`doi:`) and the manuscript.

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
