# CDPKit-Utils

This repository provides a collection of Python scripts for free use by researchers in the field of computer-aided drug design.
All scripts make heavy use of CDPKit functionality and thus require a proper CDPKit installation to run without errors
(see https://cdpkit.org/master/installation.html for information about different CDPKit installation options).

## 1. Generation of Decoy Data Sets for Docking and Pharmacophore Screening

The script *gen_decoys.py* (as closely as possible) implements the workflow that is being used for the selection of target-specific decoys
in the DUDE-Z data set [[Stein et al.]](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.0c00598).

Decoys for a given set of input molecules (e.g. a set of known actives specified by the *-i* option) are retrieved from a large 
database of diverse drug-like molecules (e.g. a suitable ZINC subset; specified by the *-d* option). In principle, any compound database can
be used as decoy source. However, the database should be large enough (in the order of several millions of molecules) so that it can be made
sure that a sufficient amount of structurally diverse decoys is found for each input molecule (default: 50 output decoys per molecule).

Decoy candidates are primarily selected by an input molecule property matching procedure. One of these properties is the net formal charge and the
tolerance for matching this property should be quite small (default: zero). Thus, to be able to find a sufficient number of decoy candidates,
possible protonation states of the input and decoy database molecules either need to be standardized or have to be enumerated (the preferred option;
can be done with a tool like *Epik* available from Schrodinger Inc.) before performing any decoy search. For the used decoy database this needs to be
done only once.

In order to speed up the search for decoys, all relevant properties of the molecules in the decoy database need to be precalculated
and stored in a CSV-file for later use.

For the precalculation of these properties the *gen_decoys.py* script offers a special mode of operation that is enabled by
the *--calc-props* command line flag:

```console
$ python3 gen_decoys.py --calc-props -d <DECOY-DB> -p <PROP-CSV-FILE>
```

When the *gen_decoys.py* script is executed in this way it will calculate all required properties of the molecules in the database file `<DECOY-DB>`
and write the results in CSV format to the file `<PROP-CSV-FILE>`. The precalculation procedure needs to be carried out
only once for the used decoy database and the generated property CSV file can then be used by all later decoy searches (note:
if the `<DECOY-DB>` file gets altered by removing or adding molecules the property calculation procedure needs to be carried out
again!).

### Generating decoys for set of molecules

```console
$ python3 gen_decoys.py -i <INPUT-MOLS> -o <GEN-DECOYS> -d <DECOY-DB> -p <PROP-CSV-FILE>
```

The input molecules (e.g. known ligands for the target of interest) are specified via the *-i* option and the found decoys (by default 50
for each input molecule) will be written to the file `<GEN-DECOYS>` specified by the *-d* option. `<DECOY-DB>` denotes the decoy source database
and `<PROP-CSV-FILE>` the associated property data file that was generated previously (see above).

To some extent the decoy search procedure can be fine-tuned by modifying the tolerance range for each matched property type. Tolerances
for each property can be changed by a corresponding command line option. The list of available options can be displayed when executing *gen_decoys.py*
with the *-h* flag.

### Implemented extension of the original protocol

The used decoy database may contain molecules that are active towards the biological target of the ligands for which the decoy search is carried out.
In order to further reduce the risk that active molecules end up in the output decoy set (rest assured, the new DUDE-Z protocol takes several measures
that greatly diminish the likeliness that this will happen), the *gen_decoys.py* script supports the specification of an 'excluded' molecule file
(*-x* option). Any decoy candidate that is contained in the set of excluded molecules will not be considered for further processing and thus not enter the
final decoy set. A reasonable way of using this option is to provide a file containing all known actives for the target of interest (e.g. retrieved
from the [ChEMBL](https://www.ebi.ac.uk/chembl/) database). Doing so will then ensure that at least none of the known active molecules is among the generated decoys.
