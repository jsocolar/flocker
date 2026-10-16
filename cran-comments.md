## R CMD check results

0 errors | 0 warnings | 0 notes

## Windows tests

On CRAN Windows only, tests that require runtime compilation of Stan models are
skipped because the current rstan/StanHeaders Windows toolchain fails during C++
compilation before those tests can run. Tests that do not require fitted model
fixtures continue to run on CRAN Windows. The complete test suite runs on macOS
and Linux and in continuous integration.

## URL checks

`urlchecker::url_check()` reports `403 Forbidden` responses for the following
three references:

- <https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13897>
- <https://doi.org/10.1002/ecy.3867>
- <https://onlinelibrary.wiley.com/doi/10.1002/ece3.9328>

These are valid links to published articles. The destinations are hosted by
Wiley Online Library, which rejects the automated requests made by the URL
checker. The `403` responses therefore reflect the publisher's automated-access
controls rather than broken or unavailable URLs.
