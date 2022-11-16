# Standardising Names In Datasets

## Installation

Run:
`pip install git+https://github.com/alrichardbollans/automatchnames.git#egg=automatchnames`

## Usage

```python
import pandas as pd
from wcvp_name_matching import get_accepted_info_from_names_in_column

data_csv = 'path_to_data.csv'

your_data_df = pd.read_csv(data_csv)  # Data to use
name_col = 'taxa'  # Name of column in data with names to check
# Names of families in your data 
# This is optional, but the program is much faster if specified
families_in_occurrences = ['Apocynaceae', 'Rubiaceae']
# Manual resolutions are optional and included by specifying a csv file, in the same format as
# the `manual_match_template.csv` file.
manual_resolution_csv = 'manual_match_template.csv'

# Match level specifies how conservative to be. One of ['full', 'weak', 'knms']
# weak: only include direct matches to wcvp
# knms: Include direct matches to wcvp and matches from KNMS
# full: include both of the above, and autoresolution step
match_level = 'full'
data_with_accepted_information = get_accepted_info_from_names_in_column(your_data_df, name_col,
                                                                        families_of_interest=families_in_occurrences,
                                                                        manual_resolution_csv=manual_resolution_csv,
                                                                        match_level=match_level)
```

## Steps

In the first step, to avoid the program spending time trying to find names we know to be problematic we do
some manual matching. Manual resolutions are optional and included by specifying a csv file, in the same
format as the `manual_match_template.csv` file. Tag= 'manual'

Once manual matches have been found, we try to match names directly to taxa in WCVP. This finds taxa in WCVP
which match our submitted names exactly. In cases where multiple taxa are returned for a given submission,
taxa labelled as 'Accepted' are prioritised. Tag= 'direct_wcvp'

Submitted names which aren't found in these first steps are then matched to names using KNMS, which contains
multiple steps. Firstly, in simple cases where KNMS returns a single match for a submitted name we use the
match IPNI ID to find accepted information from WCVP. Tag = 'knms_single'

Frequently however, submissions will be matched to multiple names in KNMS. In these cases we attempt to find
the 'best' match. To do this, first we find accepted info for each of the matches using the match IPNI ID and
WCVP. In cases where the accepted name for a given match is the same as the submitted name, we use this
match (Tag = 'knms_multiple_1').
Next, in cases where a given submitted name matches (to many) names which all have the same accepted name, we
use this accepted name (Tag = 'knms_multiple_2').

Next, for submissions which have been matched in KNMS but haven't been resolved so far we look for matches
where the accepted name is contained in the submitted name. This is useful for catching instances where author
names have been provided and so the submission may be unresolved in the previous step. In some cases, for a
single submitted name this may return multiple matches, in which case we take the most specific match (i.e. "
Subspecies" > "Variety" > "Species"> "Genus"). Tag = 'knms_multiple_3'

Once we have tried to resolve submitted names through KNMS in the above, we may still have some names left
over. In these cases we first try to do some automated resolution. In this step we search through WCVP for
taxa where the taxon name is contained in the submitted name. This is similar to the previous step but is much
slower as many more names must be checked (specifying families of interest really helps here). For each
submitted name, we then have a list (possibly empty) of taxa where the taxon name is contained in the
submitted name. We want to prioritise accepted taxa over synonyms etc.. so a given submitted name is resolved
to the best taxonomic status i.e. "Accepted" > "Synonym" > ...

As before, generic names may be contained in more specific names, so we must account for this somehow. This
is achieved for a given submitted name by resolving it to the most precise matched taxa i.e. "Subspecies" > "
Variety" > "Species"> "Genus". Moreover, say a species has been submitted where the species part of the name
has been misspelled e.g. "Neonauclea observifolia"; in these cases the genus will be the only match to the
name and this resolution would be incorrect. Furthermore, genera names can be shared across family names (
e.g. **Condylocarpus**). Therefore, when families have not been specified we don't match submissions to genera
where the submitted name contains a space and the genera are known to be contained in multiple families. Note
that this is conservative and will cause some good matches to not be matched, in particular genera given with
authors. Tag= 'autoresolution'

Finally, the resolutions are recompiled and an updated dataframe is returned. Submitted names which haven't
been matched at any point are output to a csv file for you to check. Note that unmatched submissions are
included in the output dataframe without any accepted information.

A rough diagram is given below.

![pipe](pipe.svg)

### Notes on outputs

* Output dataframe is the same as the input, with additional columns providing resolved accepted name
  information. Where names are unresolved, values in these columns are empty.
* `matched_by` column specifies how the name has been resolved. One of:
    * 'direct_wcvp': name resolved directly matching to WCVP
    * 'knms_single': where KNMS provides a single matching name
    * 'knms_multiple_1': where KNMS provides multiple matches for the submitted name, but the submitted name
      is exactly the same as the accepted name for one of the matches
    * 'knms_multiple_2': where KNMS provides multiple matches for the submitted name, but matches are all
      synonyms of the same accepted name
    * 'knms_multiple_3': where KNMS provides multiple matches for the submitted name, picks matches where
      the accepted name is contained in the submitted name. Where there are multiple such cases, the most
      specific resolutions are picked.
    * 'autoresolution': Resolutions found in autoresolution step

## Name Formatting

The program does some automatic formatting (not affecting output)

* Remove whitespace
* Capitalisations

Some formatting issues are hard to deal with and are subject of improvements. If cases are unmatched:

* Spelling errors
* Infraspecific ranks (var., subsp.) should in general end with '.' e.g. (**Psychotria guadalupensis subsp.
  grosourdieana** not **Psychotria guadalupensis subsp grosourdieana**)
* Authors and publication info end with '.'

## Notes on KNMS

* KNMS may not return anything if you submit too many names and/or requests. We mitigate this by only checking
  names in
  KNMS which can't be found in WCVP. Also, results from KNMS are stored for reuse in
  a `name matching temp outputs`
  folder.
* KNMS does not appear to account for spelling errors e.g. 'Neonauclea observifolia' returns no info (it
  should be '
  Neonauclea obversifolia').
* KNMS does not always find matches for correctly spelled accepted names. Some examples are given
  in `knms_unmatched_accepted_names.csv`.
* KNMS does not handle captilisation particularly well. For example, 'PALICOUREA GRACILENTA' is unmatched
  and 'ROTHMANIA
  ENGLERIANA (K. SCHUM.) KEAV' and 'ROTHMANIA ENGLERIANA (K. Schum) Keav' match the genus 'Rothmannia Kongl.
  Vetensk.
  Acad. Handl. 37: 63 (1776) Thunb. 1776'. Moreover, uncapitalised authors cause no matches e.g. 'Acokanthera
  deflersii
  schweinf. ex lewin' returns no match.

## Notes on WCVP

* Using most up to date version of WCVP
* Some records in WCVP are not given accepted information e.g. 'Psychotria guadalupensis subsp.
  grosourdieana', '
  Asperula nitida' or '
  Urtica angustifolia'
* Some records have ranks differing from their accepted taxa e.g. the Variety 'Diodia teres var. hirsutior' is
  a synonym
  of the Species 'Hexasepalum teres'
* Some times POWO and WCVP don't agree (mostly due to short lag in POWO updates?) as of writing **Gunnessia**
  , **Oistonema** and **Gentingia** are accepted in POWO but are synonyms in WCVP

## Notes on Kew Reconciliation Service

* KRS relies a little on manually matching unknown samples/multiple matches.
* If we try to include KRS by only including records with single matches, we may still get some errors.
* When doing automatic matching as in `open_reconciling.py` each epithet needs extracting and adding to the
  query
  otherwise e.g. 'Vaccinium vitis-idaea L.' is matched to its genus '
  Vaccinium L.'
* I've created an implementation which includes KRS but so far it is very slow (possibly because extracting
  epithets for
  lots of samples is slow).

## Possible Improvements

* Levels of strictness (specify how conservative you want to be about matching)
* Test homotypic synonyms
* Resolution of hybrid genera
* Resolution of capitalised names
* Compare to reconciliation service/include as initial step

### Known Issues

Hard cases which fail are given in `examples_to_fix.csv` in `unit_tests`.

## Sources

WCVP (2022). World Checklist of Vascular Plants. Facilitated by the Royal Botanic Gardens, Kew.
Published
on the Internet
http://wcvp.science.kew.org/
Retrieved XX/XX/XX.

KNMS (2022). Kew Names Matching Service.
http://namematch.science.kew.org/

Kew Reconciliation Service

gnparser Mozzherin, D.Y., Myltsev, A.A. & Patterson, D.J. “gnparser”: a powerful parser for scientific names
based on
Parsing Expression Grammar. BMC Bioinformatics 18, 279 (2017).https://doi.org/10.1186/s12859-017-1663-3
