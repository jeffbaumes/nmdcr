# nmdcr

* GitHub: https://github.com/jeffbaumes/nmdcr
* Documentation: https://jeffbaumes.github.io/nmdcr

This module provides functions to interact with the NMDC API,
including retrieving metadata into data frames and linking across schema objects.

## Getting started

To install:

```r
install.packages("devtools")
library(devtools)
install_github("jeffbaumes/nmdcr")
```

To use:

```r
library(nmdcr)
client <- NmdcClient$new()
studies <- client$find(nmdcr::Study)
```

## User guide

This module provides functions to interact with the NMDC API,
including retrieving metadata into data frames and linking across schema objects.

The `NmdcClient$find` and `NmdcClient$lookup` functions return data as a data frame.
See the `NMDC schema documentation <https://microbiomedata.github.io/nmdc-schema/>`_
for more information on NMDC schema classes.

The `NmdcClient$merge_related` function
merges related objects into one data frame.
This uses the more low-level `NmdcClient$related_ids` function to
find related objects of a given type from a single ID.

When you load `nmdcr`, it fetches the latest NMDC schema
and provides R representations of schema classes as `nmdcr::<class>`.
Names of classes and slots can be discovered by autocompleting `nmdc::` in your IDE, or by
visiting the `NMDC schema documentation <https://microbiomedata.github.io/nmdc-schema/>`_.

To interact with the NMDC API, create an `NmdcClient` instance.

```r
client <- NmdcClient$new()
```

The following retrieves all studies with Wrighton as principal investigator as a data frame.
Queries are specified using
`MongoDB query syntax <https://www.mongodb.com/docs/manual/tutorial/query-documents/>`_.

```r
query <- list(`principal_investigator.has_raw_value` = list(`$regex` = "Wrighton"))
studies <- client$find(nmdcr::Study, query = query)
```

The following retrieves a specific biosample by its identifier and associates it
with study metadata. The first time `NmdcClient$merge_related` is called,
it will fetch all linkages between objects in the NMDC schema and cache them in
the client object. This takes a few minutes.
You can force an update of this cache at a later time with
by calling `NmdcClient$fetch_links` with `force` set to `TRUE`.

```r
ids <- "nmdc:bsm-13-7qxjvr77"
biosample <- client$lookup(nmdcr::Biosample, ids, fields = c("id", "name"))
add_study <- client$merge_related(biosample, "id", nmdcr::Study, fields=["id", "title"])
```

The following retrieves 10 samples then merges them with all associated workflow executions.

```r
samp <- client$find(nmdcr::Biosample, fields = c("id", "type"), limit = 10)
exec <- client$merge_related(samp, "id", nmdcr::WorkflowExecution, fields=c("id", "type"))
```
