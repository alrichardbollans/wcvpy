# Standardising Names In Datasets

## Steps

Run `get_accepted_info_from_names_in_column` with Pandas dataframe containing names to standardise in `name_col`. The
program will be much quicker if you specify which families to include
e.g. `families_of_interest=['Rubiaceae','Apocynaceae']`.

In the first step, to avoid the program spending time trying to find names we know to be problematic we do some manual
matching. Manual resolutions are included by editing the given `manual_match.csv` file.

Once manual matches have been found, we try to match names directly to taxa in WCVP. This finds taxa in WCVP which match
our submitted names exactly. In cases where multiple taxa are returned for a given submission, taxa labelled as '
Accepted' are prioritised.

Submitted names which aren't found in these first steps are then matched to names using KNMS, which contains multiple
steps. Firstly, in simple cases where KNMS returns a single match for a submitted name we use the match IPNI ID to find
accepted information from WCVP.

Frequently however, submissions will be matched to multiple names in KNMS. In these cases we attempt to find the 'best'
match. To do this, first we find accepted info for each of the matches using the match IPNI ID and WCVP. In cases where
the accepted name for a given match is the same as the submitted name, we use this match. Then in cases where a given
submitted name matches (to many) names which all have the same accepted name, we use this accepted name.

Next, for submissions which have been matched in KNMS but haven't been resolved so far we look for matches where the
accepted name is contained in the submitted name. This is useful for catching instances where author names have been
provided and so the submission may be unresolved in the last step. Note: this is only applied to submissions where the
returned KNMS matches are all of the same rank. In doing this we hope to avoid submissions being incorrectly matched to
more generic taxa (though I think this is unlikely anyway).

Once we have tried to resolve submitted names through KNMS in the above, we may still have some names left over. In
these cases we first try to do some automated resolution. In this step we search through WCVP for accepted taxa where
the taxon name is contained in the submitted name. This is similar to the previous step but is much slower due as many
more names must be checked (specifying families of interest really helps here). For each submitted name, we then have a
list (possibly empty) of accepted taxa where the taxon name is contained in the submitted name. As before, generic names
may be contained in more specific names and so we must account for this somehow. This is achieved for a given submitted
name by resolving it to the most precise matched taxa i.e. "Subspecies" > "Variety" > "Species"> "Genus". Moreover, say
a species has been submitted where the species part of the name has been misspelled e.g. "Neonauclea observifolia"; in
these cases the genus will be the only match to the same and this resolution would be incorrect. We therefore don't
match submissions to genera where the submitted name contains a space. Note that this is conservative and will cause
some good matches to not be matched, in particular hybrid genera or genera given with authors.

Finally, the resolutions are recompiled and an updated dataframe is returned. Submitted names which haven't been matched
at any point are output to a csv file for you to check. Note that unmatched submissions are included in the output
dataframe without out any accepted information.

## Notes on KNMS

* KNMS may not return anything if you submit too many names and/or requests. We mitigate this by only checking names in
  KNMS which can't be found in WCVP. Also, results from KNMS are stored for reuse in a `name matching temp outputs`
  folder.
* KNMS does not appear to account for spelling errors e.g. 'Neonauclea observifolia' returns no info (it should be '
  Neonauclea obversifolia').
* KNMS does not always find matches for correctly spelled accepted names. Some examples are given
  in `knms_unmatched_accepted_names.csv`.
* KNMS does not handle captilisation particularly well. For example, 'PALICOUREA GRACILENTA' is unmatched and 'ROTHMANIA
  ENGLERIANA (K. SCHUM.) KEAV' and 'ROTHMANIA ENGLERIANA (K. Schum) Keav' match the genus 'Rothmannia Kongl. Vetensk.
  Acad. Handl. 37: 63 (1776) Thunb. 1776'. Moreover, uncapitalised authors cause no matches e.g. 'Acokanthera deflersii
  schweinf. ex lewin' returns no match.

## Notes on WCVP

* Some records in WCVP are not given accepted information e.g. 'Psychotria guadalupensis subsp. grosourdieana'
* Some records have ranks differing from their accepted taxa e.g. the Variety 'Diodia teres var. hirsutior' is a synonym
  of the Species 'Hexasepalum teres'

## Notes on Kew Reconciliation Service

* KRS relies a little on manually matching unknown samples/multiple matches.
* If we try to include KRS by only including records with single matches, we may still get some errors.
* When doing automatic matching as in `open_reconciling.py` each epithet needs extracting and adding to the query
  otherwise e.g. 'Vaccinium vitis-idaea L.' is matched to its genus '
  Vaccinium L.'
* I've created an implementation which includes KRS but so far it is very slow (possibly because extracting epithets for
  lots of samples is slow).

## Possible Improvements

* Levels of strictness (specify how conserative you want to be about matching)
* Test homotypic synonyms
* Resolution of hybrid genera
* Resolution of capitalised names
* Provide more info in output e.g. genus/parent
* Compare to reconciliation service/include as initial step

## Sources

WCVP (2022). World Checklist of Vascular Plants, version 2.0. Facilitated by the Royal Botanic Gardens, Kew. Published
on the Internet
http://wcvp.science.kew.org/
Retrieved 20/01/2022.

KNMS (2022). Kew Names Matching Service.
http://namematch.science.kew.org/

Kew Reconciliation Service

gnparser Mozzherin, D.Y., Myltsev, A.A. & Patterson, D.J. “gnparser”: a powerful parser for scientific names based on
Parsing Expression Grammar. BMC Bioinformatics 18, 279 (2017).https://doi.org/10.1186/s12859-017-1663-3
