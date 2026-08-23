## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

  The same NOTE reports the `BugReports` URL as possibly invalid and
  suggests

      https://gitlab.com/dickoa/svyplan/-/work_items/issues

  The check heuristic appends `/issues` to the `BugReports` value. This
  package is hosted on GitLab, where the value given already is the issue
  tracker, so the suggestion is that tracker plus a further path segment the
  project does not serve. The URL is correct as given.

## Test environments

* Local: Arch Linux, R 4.6.1
* GitHub Actions: macOS-latest (R-release)
* GitHub Actions: windows-latest (R-release)
* GitHub Actions: ubuntu-latest (R-release, R-devel, R-oldrel-1)

## Downstream dependencies

None (new package).
