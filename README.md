# CDPKit-Utils

This repository provides a collection of Python scripts for free use by researchers in the field of computer-aided drug design.
All scripts make heavy use of CDPKit functionality and thus require a proper CDPKit installation to run (see https://cdpkit.org/master/installation.html
for information about different CDPKit installation options).

## 1. Generation of Decoy Data Sets for Docking and Pharmacophore Screening

The script *gen_decoys.py* as closely as possible implements the workflow that was used for the selection of target-specific decoys
in the DUDE-Z data set [[Stein et al.]](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.0c00598).

Decoys for a given set of input molecules (e.g. a set of known actives specified by the *-i* option) are retrieved from a large 
database of diverse drug-like molecules (e.g. a ZINC subset specified by the *-d* option). In principle, any compound database can
be used as a decoy source. However, it should be large enough (in the range of several millions of molecules) to be able find enough
structurally diverse decoys for each provided input molecule (default: 50 per input molecule).

Decoy candidates are primarily selected by input molecule property matching. One of these properties is the net formal charge and the
default tolerance for matching this property is zero. Thus, for finding decoy candidates to be successful, possible protonation states
of the input and decoy database molecules should be enumerated (e.g. using a tool like *LigPrep* available from Schrodinger Inc.)
before performing the decoy search. For the decoy database this needs to be done only once.

In order to speed up the search for decoys all relevant properties of the molecules in the decoy source database need to be precalculated
and stored in a CSV-file for later use.

For precalculating these properties the *gen_decoys.py* script provided a special mode of operation that is enabled by
specifying the *--calc-props* flag on the command line:

```console
$ python3 gen_decoys.py --calc-props -d <DECOY-DB> -p <OUT-CSV-FILE>
```

This will precalculate all required properties of the molecules in the decoy database specified by `<DECOY-DB>` and writes
the results in CSV format into the output file `<OUT-CSV-FILE>`. The precalculation procedure needs to be carried out
only once for the used decoy database and the generated output file can be reused by later decoy generation runs.
