## R CMD check results

On Arch Linux (x86_64) with R 4.6.2 Patched, the source tarball returned
`Status: OK` from `R CMD check --no-tests --no-manual`. The full suite
ran separately with 33,820 passing expectations and zero failures,
warnings or skips. The PDF reference manual also built locally.
The network-enabled URL check passed for the package sources.

`R CMD check --as-cran --no-tests --no-manual` reported 0 errors,
0 warnings and 2 notes with remote incoming checks enabled. The local
run with remote incoming checks disabled reported the same notes.

- The CRAN incoming check identifies a new submission and suggests
  appending `/issues` to `BugReports`. The configured GitLab work-item
  URL returned HTTP 200. The suggested path redirected to sign-in.
- The compilation flags named in the note come from this Arch Linux R
  installation's default `CFLAGS`.

## Test environments

* Local: Arch Linux, R 4.6.2
* GitHub Actions: macOS-latest (R-release)
* GitHub Actions: windows-latest (R-release)
* GitHub Actions: ubuntu-latest (R-release, R-devel, R-oldrel-1)

## Downstream dependencies

None (new package).
