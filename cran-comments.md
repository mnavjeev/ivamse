## Test environments

* local macOS, R 4.1.3

## R CMD check results

0 errors | 0 warnings | 2 notes

* "New submission" -- this is a first submission.
* "unable to verify current time" -- the check machine could not reach a time
  server; unrelated to the package.

## Notes

Examples and tests run in a few seconds. The test suite includes an optional
comparison against the replication code for the paper the package implements,
which is skipped unless the environment variable IVAMSE_PAPER_REPO points at a
checkout of that repository.
