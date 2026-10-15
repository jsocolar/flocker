## R CMD check results

0 errors | 0 warnings | 0 notes

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
