Methods in this package are based on: https://gist.github.com/nickynicolson/11fe9e57a198d31fa010fb3feaa65d94

# Collected Issues
## GUI vs. API

Using the GUI version 3.7.2 on linux with and reconciling
to  https://data1.kew.org/reconciliation/reconcile/IpniName service
get slightly different results than with this API set up:

* 'Palicourea gracilenta (Müll. Arg.) Delprete & J. H. Kirkbr.' doesn't match with API but matches with GUI

## Capitalisation

* 'ROTHMANIA ENGLERIANA (K. SCHUM.) KEAV' only matches to the genus in but 'Rothmania engleriana (k. schum.)
  keav' matches correctly.

## Unmatched examples
Some examples of names where matches aren't found are given in [unmatched_in.csv](unittests%2Ftest_inputs%2Funmatched_in.csv)