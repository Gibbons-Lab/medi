# Why don't we simply provide the DB for download?

Believe us, we would love nothing more than to be able to do this. Unfortunately, the
MEDI database is very large (~500GB) and compressing it does not help much (~350GB).
Right now there is no public data repository that allows data sets of this size.

For the sake of completeness, here is what we tried up to now.

## Zenodo

Zenodo is a great data repository. However, maximum data uploads are limited to
50GB with a possible extension to an absolute maximum of 200GB.

We considered splitting up the database into several ZENODO depositions but this is
specifically forbidden under the [ZENODO fair-use regulations](https://support.zenodo.org/help/en-gb/1-upload-deposit/80-what-are-the-size-limitations-of-zenodo).

## Self-sharing using AWS S3 or other cloud providers

AWS S3 has no size limitations and actual storage costs are manageable for large data sets
(~20-50 USD per month). However, S3 charges for egress, meaning for transferring data to S3
to the internet. This means we need to pay for every download of the data set and this is
expensive (~50-100 USD for every *single* download). So if the database is downloaded 1000 times per month the projected costs would be 100,000 USD per month, which is prohibitively expensive for us.

Pricing for competing cloud providers is similar.

## Sharing using Google Drive

Google Drive can also host large files and does not charge for download. However, there
are strict API limits for publically shared files. In our tests sharing the database per
an open link made it virtually impossible to download as there would be repeated errors.

> [!NOTE]
> We did find a workaround with personalized sharing (sharing with specific Google accounts)
> only, in particular when using the rclone tool. So if you need the version of the DB from
> the publication please contact us at lab[at]cdiener.com with your Gmail address and we
> can share it with you.