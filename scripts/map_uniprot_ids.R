library(httr)

isJobReady <- function(jobId) {
    pollingInterval = 5
    nTries = 20
    for (i in 1:nTries) {
        url <- paste("https://rest.uniprot.org/idmapping/status/", jobId, sep = "")
        r <- GET(url = url, accept_json())
        status <- content(r, as = "parsed")
        if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
            return(TRUE)
        }
        if (!is.null(status[["messages"]])) {
            print(status[["messages"]])
            return (FALSE)
        }
        Sys.sleep(pollingInterval)
    }
    return(FALSE)
}

getResultsURL <- function(redirectURL) {
    if (grepl("/idmapping/results/", redirectURL, fixed = TRUE)) {
        url <- gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
    } else {
        url <- gsub("/results/", "/results/stream/", redirectURL)
    }
}

myids <- readRDS('/home/fstruebi/projects/SNCA_MAPT_network/export/uniprot_ids_list.rds') %>% 
    lapply(., function(x) {
        unlist(strsplit(x, split = '\\|')[1])[1]
    }) %>% unlist(.) %>% str_remove(., pattern = '\\-[0-9]')

files = list(
    ids = paste(myids, collapse = ','),
    from = "UniProtKB_AC-ID",
    to = "Gene_Name"
)

r <- POST(url = "https://rest.uniprot.org/idmapping/run", body = files, encode = "multipart", accept_json())
submission <- httr::content(r, as = "parsed")

if (isJobReady(submission[["jobId"]])) {
    url <- paste("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]], sep = "")
    r <- GET(url = url, accept_json())
    details <- httr::content(r, as = "parsed")
    url <- getResultsURL(details[["redirectURL"]])
    # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available
    url <- paste(url, "?format=tsv", sep = "")
    r <- GET(url = url, accept_json())
    resultsTable = read.table(text = httr::content(r), sep = "\t", header=TRUE)
    # print(resultsTable)
}

