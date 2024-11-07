# Load necessary libraries
library(httr)
library(jsonlite)
library(dplyr)
library(yaml)

schema <- yaml.load_file("https://raw.githubusercontent.com/microbiomedata/nmdc-schema/refs/heads/main/nmdc_schema/nmdc_materialized_patterns.yaml")
for (symbol in names(schema$classes)) {
  assign(symbol, schema$classes[[symbol]])
}

class_descendants = function(cls) {
  descendants <- c(cls$name)
  for (symbol in names(schema$classes)) {
    if (!is.null(schema$classes[[symbol]]$is_a) && schema$classes[[symbol]]$is_a == cls$name) {
      descendants <- c(descendants, class_descendants(schema$classes[[symbol]]))
    }
  }
  return(descendants)
}

class_ancestors = function(cls) {
  ancestors <- c(cls$name)
  while (!is.null(cls$is_a)) {
    ancestors <- c(ancestors, cls$is_a)
    cls <- schema$classes[[cls$is_a]]
  }
  return(ancestors)
}

database_collections_for_type = function(cls) {
  descendants <- class_descendants(cls)
  collections = c()
  for (collection in schema$classes$Database$slots) {
    slot <- schema$slots[[collection]]
    collection_descendants <- class_descendants(schema$classes[[slot$range]])
    common_elements <- intersect(descendants, collection_descendants)
    if (length(common_elements) > 0) {
      collections <- c(collections, collection)
    }
  }
  return(collections)
}

NmdcClient <- setRefClass(
  "NmdcClient",
  fields = list(
    base_url = "character",
    max_page_size = "numeric",
    sleep_seconds = "numeric",
    timeout = "numeric",
    verbose = "logical",
    links = "ANY",
    codes = "ANY"
  ),
  methods = list(
    initialize = function(base_url = "https://api.microbiomedata.org", max_page_size = 1000, sleep_seconds = 0.5, timeout = 30, verbose = FALSE, links = NULL) {
      .self$base_url <- base_url
      .self$max_page_size <- max_page_size
      .self$sleep_seconds <- sleep_seconds
      .self$timeout <- timeout
      .self$verbose <- verbose
      .self$links <- links
      .self$codes <- NULL
    },

    fetch_links = function(force = FALSE) {
      if (!force && !is.null(.self$links)) {
        return(.self$links)
      }

      # # Check if the RDS file exists and load it
      # if (file.exists("links.rds")) {
      #   .self$links <- readRDS("links.rds")
      #   return(.self$links)
      # }

      row_chunk_size <- 10000
      links_env <- new.env()
      links_env$links <- data.frame(source = character(row_chunk_size), target = character(row_chunk_size))
      links_env$current_row <- 1

      add_link <- function(env, source, target) {
        if (env$current_row > nrow(env$links)) {
          env$links <- bind_rows(env$links, data.frame(source = character(row_chunk_size), target = character(row_chunk_size)))
        }
        env$links[env$current_row, ] <- list(source, target)
        env$current_row <- env$current_row + 1
      }

      studies <- .self$find("study_set", fields=c("id", "part_of"))
      apply(studies, 1, function(x) {
        if (!is.null(x$part_of)) {
          for (input_id in x$part_of) {
            add_link(links_env, input_id, x$id)
          }
        }
      })

      material_processings <- .self$find("material_processing_set", fields=c("id", "has_input", "has_output"))
      apply(material_processings, 1, function(x) {
        for (input_id in x$has_input) {
          add_link(links_env, input_id, x$id)
        }
        for (output_id in x$has_output) {
          add_link(links_env, x$id, output_id)
        }
      })

      data_generations <- .self$find("data_generation_set", fields=c("id", "has_input", "has_output", "associated_studies"))
      apply(data_generations, 1, function(x) {
        for (input_id in x$has_input) {
          add_link(links_env, input_id, x$id)
        }
        if (!is.null(x$has_output)) {
          for (output_id in x$has_output) {
            add_link(links_env, x$id, output_id)
          }
        }
        for (study_id in x$associated_studies) {
          add_link(links_env, study_id, x$id)
        }
      })

      workflow_executions <- .self$find("workflow_execution_set", fields=c("id", "was_informed_by", "has_input", "has_output"))
      apply(workflow_executions, 1, function(x) {
        for (informed_id in x$was_informed_by) {
          add_link(links_env, informed_id, x$id)
        }
        for (input_id in x$has_input) {
          add_link(links_env, input_id, x$id)
        }
        for (output_id in x$has_output) {
          add_link(links_env, x$id, output_id)
        }
      })

      biosamples <- .self$find("biosample_set", fields=c("id", "associated_studies"))
      apply(biosamples, 1, function(x) {
        for (study_id in x$associated_studies) {
          add_link(links_env, study_id, x$id)
        }
      })

      # # Cache the links
      # saveRDS(links_env$links, file = "links.rds")

      .self$links <- links_env$links[1:(links_env$current_row - 1), ]
    },

    fetch_codes = function(force = FALSE) {
      if (!force && !is.null(.self$codes)) {
        return(.self$codes)
      }

      url <- paste0(.self$base_url, "/nmdcschema/typecodes")
      if (.self$verbose) {
        print(url)
      }

      response <- GET(url, timeout(.self$timeout))
      data <- content(response, as = "text", encoding = "UTF-8")
      parsed <- fromJSON(data, flatten = TRUE)

      .self$codes <- parsed
    },

    find = function(cls, query = list(a = NULL)[-1], fields = NULL, limit = NULL) {
      if (is.null(fields)) {
        fields <- list()
      }

      collection <- NULL
      if (is.character(cls)) {
        collection <- cls
      } else {
        collections <- database_collections_for_type(cls)

        if(length(collections) == 0) {
          print(paste0("Note: No collection found for ", cls$name))
          return (data.frame())
        }
        if(length(collections) > 1) {
          print(paste0("Note: Multiple collections found for ", cls$name, ": ", collections, ". Using the first one."))
        }
        collection <- collections[1]

        # Add the class to the filter (for heterogeneous collections)
        query[["type"]] <- list("$in" = I(paste0("nmdc:", class_descendants(cls))))
      }

      url <- paste0(.self$base_url, "/nmdcschema/", collection, "?filter=", toJSON(query, auto_unbox = TRUE), "&projection=", paste(fields, collapse = ","), "&max_page_size=", .self$max_page_size)
      if (.self$verbose) {
        print(url)
      }

      response <- GET(url, timeout(.self$timeout))
      data <- content(response, as = "text", encoding = "UTF-8")
      parsed <- fromJSON(data, flatten = TRUE)
      results <- parsed$resources

      while (!is.null(parsed$next_page_token) && (is.null(limit) || nrow(results) < limit)) {
        Sys.sleep(.self$sleep_seconds)
        page_url <- paste0(url, "&page_token=", parsed$next_page_token)
        if (.self$verbose) {
          print(page_url)
        }

        response <- GET(page_url, timeout(.self$timeout))
        data <- content(response, as = "text", encoding = "UTF-8")
        parsed <- fromJSON(data, flatten = TRUE)
        results <- bind_rows(results, parsed$resources)
      }

      if (is.null(limit)) {
        return(results)
      } else {
        return(head(results, limit))
      }
    },

    lookup = function(cls, ids, fields = NULL) {
      chunked_ids <- split(ids, ceiling(seq_along(ids) / 100))
      results <- NULL

      for (chunk in chunked_ids) {
        data <- .self$find(cls, query = list(id = list(`$in` = I(chunk))), fields = fields)
        if (is.null(results)) {
          results <- data
        } else {
          results <- bind_rows(results, data)
        }
      }

      return(results)
    },

    follow_links = function(root_id, direction, steps = NULL) {
      if (is.null(.self$links)) {
        .self$fetch_links()
      }

      result <- c()
      frontier <- c(root_id)
      while ((is.null(steps) || step < steps) && length(frontier) > 0) {
        result <- c(result, frontier)
        if (direction == 1) {
          frontier <- .self$links %>% filter(target %in% frontier) %>% pull(source)
        } else {
          frontier <- .self$links %>% filter(source %in% frontier) %>% pull(target)
        }
      }

      return(unique(result))
    },

    related_ids = function(ids, cls) {
      descendants <- class_descendants(cls)
      if (is.null(.self$codes)) {
        .self$fetch_codes()
      }
      descendant_codes <- client$codes %>% filter(schema_class %in% paste0("nmdc:", descendants)) %>% pull(name)
      source <- c()
      target <- c()
      for (source_id in ids) {
        linked <- unique(c(.self$follow_links(source_id, 1), .self$follow_links(source_id, 2)))
        for (d in linked) {
          if (any(sapply(descendant_codes, function(code) startsWith(d, paste0("nmdc:", code))))) {
            source <- c(source, source_id)
            target <- c(target, d)
          }
        }
      }
      return(data.frame(source = source, target = target))
    },

    merge_related = function(data, id_column, cls, fields = NULL) {
      # Ensure the ID column is included in the fields
      if (!is.null(fields) && !("id" %in% fields)) {
        fields <- c(fields, "id")
      }

      if (is.null(fields)) {
        fields <- c()
      }

      related <- .self$related_ids(data[[id_column]], cls)
      related_id_column <- paste0(id_column, "_", cls$name)
      related_df <- setNames(data.frame(a = related$source, b = related$target), c(id_column, related_id_column))

      with_ids_df <- merge(data, related_df, by = id_column, all.y = TRUE)

      related_deduplicated <- unique(related$target)
      related_data_df <- .self$lookup(cls, related_deduplicated, fields)
      names(related_data_df)[names(related_data_df) == "id"] <- related_id_column

      result <- merge(with_ids_df, related_data_df, by = related_id_column, all.x = TRUE, suffixes = c("", paste0("_", cls$name)))

      return(result)
    }

  )
)
