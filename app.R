#### app.R
# Shiny app for DADA2 (ASV + SILVA) and Kraken2 classification
# Local single-core app. Supports single-end and paired-end uploads (user-provided pairing).
# Includes shinyFiles browse buttons, full SILVA/Kraken validation, and pair preview/confirm.

library(shiny)
library(shinyFiles)
library(DT)
library(fs)   # for path helpers

# Increase file upload size (500 MB per file)
options(shiny.maxRequestSize = 500 * 1024^2)  # 500 MB

# ---- UI ----
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
    .indicator {
      font-size: 18px;
      font-weight: bold;
      margin-left: 8px;
    }
  "))
  ),
  titlePanel("DADA2 + SILVA and Kraken2 — Local Shiny App (Validated)"),
  sidebarLayout(
    sidebarPanel(
      radioButtons(
        "mode", "Choose mode:",
        choices = c("DADA2 + SILVA", "Kraken2 (fast)"),
        selected = "DADA2 + SILVA"
      ),
      radioButtons(
        "layout", "Read layout:",
        choices = c("Single-end", "Paired-end"),
        selected = "Paired-end"
      ),
      
      # Dynamic file upload UI (single / paired)
      uiOutput("fastq_ui"),
      
      hr(),
      
      # SILVA FASTA path + Browse (file)
      fluidRow(
        column(10, textInput("silva_path", "SILVA FASTA (for DADA2)", value = "")),
        column(2, htmlOutput("silva_indicator"))
      ),
      shinyFilesButton("silva_file", "Browse SILVA FASTA", title = "Select SILVA FASTA file", multiple = FALSE),
      actionButton("validate_silva", "Validate SILVA now"),
      verbatimTextOutput("silva_validation_log"),
      br(),
      
      # Kraken2 DB path + Browse (directory)
      fluidRow(
        column(10, textInput("kraken_db", "Kraken2 DB directory", value = "")),
        column(2, htmlOutput("kraken_indicator"))
      ),
      shinyDirButton("kraken_dir", "Browse Kraken2 DB folder", title = "Select Kraken2 DB folder"),
      actionButton("validate_kraken", "Validate Kraken2 DB now"),
      verbatimTextOutput("kraken_validation_log"),
      br(),
      
      hr(),
      h4("DADA2 trimming parameters"),
      numericInput("truncLenF", "truncLenF (0 = no truncation)", value = 240, min = 0),
      numericInput("truncLenR", "truncLenR (0 = no truncation)", value = 200, min = 0),
      numericInput("maxEE", "maxEE", value = 2, min = 0),
      numericInput("truncQ", "truncQ", value = 2, min = 0),
      checkboxInput("remove_chimera", "Remove chimeras (DADA2)", value = TRUE),
      
      hr(),
      checkboxInput("run_kraken_report", "Generate Kraken2 report & summary (Kraken2 mode)", value = TRUE),
      
      hr(),
      # The Run button will be prevented from running when paired-end pairs are not confirmed
      actionButton("run", "Run Pipeline", class = "btn-primary"),
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Pair Preview",
                 h4("Paired-end file preview (confirm before running)"),
                 DT::dataTableOutput("pair_preview"),
                 actionButton("confirm_pairs", "Confirm pairs (allow Run)"),
                 verbatimTextOutput("pair_confirm_log")
        ),
        tabPanel("Logs", verbatimTextOutput("log")),
        tabPanel("DADA2 Outputs", DT::dataTableOutput("dada_tax")),
        tabPanel("Kraken2 Outputs", DT::dataTableOutput("kraken_tbl")),
        tabPanel("Download", uiOutput("download_ui"))
      ),
      width = 9
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  rv <- reactiveValues(
    log = "",
    dada_tax = NULL,
    kraken = NULL,
    outdir = NULL,
    zipfile = NULL,
    pair_df = NULL,
    pairs_confirmed = FALSE,
    silva_validation = NULL,
    kraken_validation = NULL
  )
  
  append_log <- function(text) {
    rv$log <- paste(rv$log, paste0(Sys.time(), ": ", text), sep = "\n")
  }
  
  # Setup shinyFiles roots (volumes)
  roots <- c(Home = path_home(), shinyFiles::getVolumes()())
  
  # Configure shinyFiles
  shinyFileChoose(input, "silva_file", roots = roots, filetypes = c("fa", "fasta", "gz", "txt"))
  shinyDirChoose(input, "kraken_dir", roots = roots)
  
  # When user picks a SILVA file via the Browse button, update the text input
  observeEvent(input$silva_file, {
    req(input$silva_file)
    fileinfo <- parseFilePaths(roots, input$silva_file)
    if (nrow(fileinfo) > 0) {
      upd <- as.character(fileinfo$datapath[1])
      updateTextInput(session, "silva_path", value = upd)
      append_log(paste("SILVA file selected:", upd))
    }
  })
  
  # When user picks a Kraken2 directory via the Browse button, update the text input
  observeEvent(input$kraken_dir, {
    req(input$kraken_dir)
    dirpath <- parseDirPath(roots, input$kraken_dir)
    if (!is.null(dirpath) && length(dirpath) > 0) {
      upd <- as.character(dirpath)
      updateTextInput(session, "kraken_db", value = upd)
      append_log(paste("Kraken2 DB folder selected:", upd))
    }
  })
  
  # Dynamic upload UI
  output$fastq_ui <- renderUI({
    if (is.null(input$layout)) return(NULL)
    if (input$layout == "Paired-end") {
      tagList(
        fileInput(
          "fastqF",
          "Forward Reads (R1) — upload one file per sample (multiple allowed)",
          multiple = TRUE,
          accept = c(".fastq", ".fq", ".fastq.gz", ".fq.gz", "application/gzip")
        ),
        fileInput(
          "fastqR",
          "Reverse Reads (R2) — upload corresponding reverse files in the same order",
          multiple = TRUE,
          accept = c(".fastq", ".fq", ".fastq.gz", ".fq.gz", "application/gzip")
        )
      )
    } else {
      fileInput(
        "fastq_single",
        "Single-end Reads — upload one file (or multiple single-end samples)",
        multiple = TRUE,
        accept = c(".fastq", ".fq", ".fastq.gz", ".fq.gz", "application/gzip")
      )
    }
  })
  
  # -----------------------
  # SILVA Validation (A3)
  # -----------------------
  observeEvent(input$validate_silva, {
    rv$silva_validation <- NULL
    silva_path <- input$silva_path
    log_msgs <- character()
    
    if (is.null(silva_path) || silva_path == "" || !file.exists(silva_path)) {
      log_msgs <- c(log_msgs, "SILVA path is empty or file does not exist.")
      rv$silva_validation <- list(ok = FALSE, messages = log_msgs)
      output$silva_validation_log <- renderText(paste(log_msgs, collapse = "\n"))
      return()
    }
    
    # Basic checks
    fi <- file.info(silva_path)
    log_msgs <- c(log_msgs, paste("SILVA file found. Size (MB):", round(fi$size / 1024^2, 2)))
    
    # size check (warn if small)
    if (fi$size < 10 * 1024^2) {
      log_msgs <- c(log_msgs, "Warning: SILVA filesize < 10 MB — file may be incomplete or wrong.")
    }
    
    # check for FASTA headers '>' in first 200 lines
    con <- gzfile(silva_path, open = "rt")
    first_lines <- tryCatch(readLines(con, n = 200), error = function(e) { close(con); character(0) }, finally = { close(con) })
    if (length(first_lines) == 0) {
      log_msgs <- c(log_msgs, "Could not read SILVA file — not readable or empty.")
      rv$silva_validation <- list(ok = FALSE, messages = log_msgs)
      output$silva_validation_log <- renderText(paste(log_msgs, collapse = "\n"))
      return()
    }
    has_header <- any(grepl("^>", first_lines))
    if (!has_header) {
      log_msgs <- c(log_msgs, "Warning: No '>' FASTA headers detected in first 200 lines. The file may not be a FASTA.")
    } else {
      log_msgs <- c(log_msgs, "FASTA header '>' detected.")
    }
    
    # Final - test assignTaxonomy on a dummy short sequence to ensure compatibility (A3)
    # This requires dada2 package; run only if available
    if (!requireNamespace("dada2", quietly = TRUE)) {
      log_msgs <- c(log_msgs, "Package 'dada2' not installed. Cannot perform assignTaxonomy test. Run setup.R to install it.")
      rv$silva_validation <- list(ok = FALSE, messages = log_msgs)
      output$silva_validation_log <- renderText(paste(log_msgs, collapse = "\n"))
      return()
    }
    
    append_log("Running small assignTaxonomy test against provided SILVA file (this may take a few seconds)...")
    test_seq <- c("ACGTACGTACGTACGTACGT")  # tiny dummy sequence
    tax_test <- tryCatch({
      suppressWarnings({
        # assignTaxonomy expects character vector of sequences
        dada2::assignTaxonomy(test_seq, silva_path, multithread = FALSE)
      })
    }, error = function(e) e)
    
    if (inherits(tax_test, "error")) {
      log_msgs <- c(log_msgs, paste("assignTaxonomy test failed:", tax_test$message))
      rv$silva_validation <- list(ok = FALSE, messages = log_msgs)
    } else {
      log_msgs <- c(log_msgs, "assignTaxonomy test succeeded. SILVA file appears compatible with dada2::assignTaxonomy().")
      rv$silva_validation <- list(ok = TRUE, messages = log_msgs)
    }
    
    output$silva_validation_log <- renderText(paste(log_msgs, collapse = "\n"))
  })
  
  output$silva_indicator <- renderUI({
    if (is.null(rv$silva_validation)) {
      HTML("<span class='indicator' style='color: grey;'>&#9675;</span>")  # grey circle
    } else if (rv$silva_validation$ok) {
      HTML("<span class='indicator' style='color: green;'>&#10004;</span>") # green check
    } else {
      HTML("<span class='indicator' style='color: red;'>&#10006;</span>")   # red cross
    }
  })
  
  # -----------------------
  # Kraken2 Validation (B3)
  # -----------------------
  observeEvent(input$validate_kraken, {
    rv$kraken_validation <- NULL
    kraken_path <- input$kraken_db
    log_msgs <- character()
    
    if (is.null(kraken_path) || kraken_path == "" || !dir.exists(kraken_path)) {
      log_msgs <- c(log_msgs, "Kraken2 DB path is empty or does not exist.")
      rv$kraken_validation <- list(ok = FALSE, messages = log_msgs)
      output$kraken_validation_log <- renderText(paste(log_msgs, collapse = "\n"))
      return()
    }
    
    # Check for .k2d files (a standard Kraken2 DB file extension)
    k2d_files <- list.files(kraken_path, pattern = "\\.k2d$", recursive = TRUE, full.names = TRUE)
    if (length(k2d_files) == 0) {
      log_msgs <- c(log_msgs, "No .k2d files found in the selected directory (expected for Kraken2 DB).")
    } else {
      log_msgs <- c(log_msgs, paste(length(k2d_files), ".k2d files found."))
    }
    
    # Folder size heuristic
    total_size <- sum(file.info(list.files(kraken_path, full.names = TRUE, recursive = TRUE))$size, na.rm = TRUE)
    log_msgs <- c(log_msgs, paste("Kraken DB total size (MB):", round(total_size / 1024^2, 2)))
    if (total_size < 50 * 1024^2) {
      log_msgs <- c(log_msgs, "Warning: Kraken2 DB size < 50 MB. This is unusually small for a Kraken2 DB and may be incomplete.")
    }
    
    # Check kraken2 binary availability and try a DB-specific dry-run if possible
    kraken_check <- tryCatch({
      # Try kraken2 --db <path> --version (some builds accept it); capture status
      cmd_db_version <- sprintf("kraken2 --db '%s' --version 2>&1", kraken_path)
      dbver <- system(cmd_db_version, intern = TRUE, ignore.stdout = FALSE, ignore.stderr = FALSE)
      list(success = TRUE, output = dbver)
    }, error = function(e) {
      # fallback: try kraken2 --version
      tryCatch({
        v <- system("kraken2 --version 2>&1", intern = TRUE)
        list(success = TRUE, output = v)
      }, error = function(e2) {
        list(success = FALSE, output = paste0("kraken2 not found or failed: ", e2$message))
      })
    })
    
    if (!isTRUE(kraken_check$success)) {
      log_msgs <- c(log_msgs, paste("Kraken2 binary test failed:", kraken_check$output))
      rv$kraken_validation <- list(ok = FALSE, messages = log_msgs)
      output$kraken_validation_log <- renderText(paste(log_msgs, collapse = "\n"))
      return()
    } else {
      log_msgs <- c(log_msgs, "Kraken2 binary seems available. Output of version test:")
      log_msgs <- c(log_msgs, paste(kraken_check$output, collapse = "\n"))
    }
    
    # Very small runtime test: try running kraken2 with --help to ensure it runs (no DB required)
    kraken_help <- tryCatch({
      h <- system("kraken2 --help 2>&1", intern = TRUE)
      list(ok = TRUE, out = h)
    }, error = function(e) {
      list(ok = FALSE, out = paste0("kraken2 help failed: ", e$message))
    })
    
    if (!kraken_help$ok) {
      log_msgs <- c(log_msgs, paste("kraken2 execution test failed:", kraken_help$out))
      rv$kraken_validation <- list(ok = FALSE, messages = log_msgs)
    } else {
      log_msgs <- c(log_msgs, "kraken2 executed --help successfully.")
      rv$kraken_validation <- list(ok = TRUE, messages = log_msgs)
    }
    
    output$kraken_validation_log <- renderText(paste(log_msgs, collapse = "\n"))
  })
  
  output$kraken_indicator <- renderUI({
    if (is.null(rv$kraken_validation)) {
      HTML("<span class='indicator' style='color: grey;'>&#9675;</span>")  # grey circle
    } else if (rv$kraken_validation$ok) {
      HTML("<span class='indicator' style='color: green;'>&#10004;</span>") # green check
    } else {
      HTML("<span class='indicator' style='color: red;'>&#10006;</span>")   # red cross
    }
  })
  
  # -----------------------
  # Paired-end preview & confirmation (C3)
  # -----------------------
  observe({
    # Build pair_df whenever files are uploaded (Paired-end mode)
    if (is.null(input$layout) || input$layout != "Paired-end") {
      rv$pair_df <- NULL
      rv$pairs_confirmed <- FALSE
      output$pair_preview <- DT::renderDataTable(NULL)
      return()
    }
    
    
    filesF <- input$fastqF
    filesR <- input$fastqR
    if (is.null(filesF) || is.null(filesR)) {
      rv$pair_df <- NULL
      rv$pairs_confirmed <- FALSE
      output$pair_preview <- DT::renderDataTable(NULL)
      return()
    }
    
    # If numbers mismatch, still create a table with NA for missing entries
    n <- max(nrow(filesF), nrow(filesR))
    f_names <- if (nrow(filesF) >= 1) filesF$name else rep(NA, n)
    r_names <- if (nrow(filesR) >= 1) filesR$name else rep(NA, n)
    
    # Truncate/prefix logic: attempt to infer sample prefix for each pair
    # We will compute a simple common-prefix before _R1/_R2 or before first underscore
    get_prefix <- function(x) {
      if (is.na(x) || x == "") return(NA_character_)
      # remove common R1/R2 patterns
      y <- gsub("(_R?1|_R?2|_1|_2)(\\.|_).*", "", x)
      # also remove file extension
      y <- sub("(\\.f(ast)?q(\\.gz)?)$", "", y, ignore.case = TRUE)
      y
    }
    prefixesF <- sapply(f_names, get_prefix)
    prefixesR <- sapply(r_names, get_prefix)
    
    match_flag <- mapply(function(a, b) {
      if (is.na(a) || is.na(b)) return(FALSE)
      tolower(a) == tolower(b)
    }, prefixesF, prefixesR, SIMPLIFY = TRUE)
    
    pair_df <- data.frame(
      index = seq_len(n),
      forward = if (length(f_names) < n) c(f_names, rep(NA, n - length(f_names))) else f_names,
      reverse = if (length(r_names) < n) c(r_names, rep(NA, n - length(r_names))) else r_names,
      prefixF = if (length(prefixesF) < n) c(prefixesF, rep(NA, n - length(prefixesF))) else prefixesF,
      prefixR = if (length(prefixesR) < n) c(prefixesR, rep(NA, n - length(prefixesR))) else prefixesR,
      matched = match_flag,
      stringsAsFactors = FALSE
    )
    
    rv$pair_df <- pair_df
    rv$pairs_confirmed <- FALSE
    output$pair_preview <- DT::renderDataTable({
      df <- pair_df
      df$matched <- ifelse(df$matched, "OK", "MISMATCH")
      DT::datatable(df[, c("index", "forward", "reverse", "matched")], options = list(pageLength = 10, scrollX = TRUE))
    })
    output$pair_confirm_log <- renderText("Please confirm pairs if they look correct (click Confirm pairs). The Run button will be disabled until confirmation.")
  })
  
  observeEvent(input$confirm_pairs, {
    if (is.null(rv$pair_df)) {
      rv$pairs_confirmed <- FALSE
      output$pair_confirm_log <- renderText("No pairs to confirm.")
    } else {
      # if any mismatches exist, warn user but still allow confirm after warning
      if (any(!rv$pair_df$matched)) {
        output$pair_confirm_log <- renderText(paste0("Warning: Some pairs do not match by filename prefix. Confirmed anyway; ensure ordering is correct. (Mismatches: ", sum(!rv$pair_df$matched), ")"))
      } else {
        output$pair_confirm_log <- renderText("Pairs confirmed OK.")
      }
      rv$pairs_confirmed <- TRUE
    }
  })
  
  # -----------------------
  # Main Run observer
  # -----------------------
  observeEvent(input$run, {
    # Prevent run when paired mode pairs not confirmed
    if (input$layout == "Paired-end" && !isTRUE(rv$pairs_confirmed)) {
      append_log("Paired-end pairs are not confirmed. Please confirm pairs in the 'Pair Preview' tab before running.")
      return()
    }
    
    # Reset reactive outputs/logs for run
    rv$log <- ""
    rv$dada_tax <- NULL
    rv$kraken <- NULL
    rv$outdir <- tempfile("dada2_app_")
    dir.create(rv$outdir, recursive = TRUE)
    append_log(paste("Created output directory:", rv$outdir))
    
    # Validate uploads depending on layout and save uploaded files
    if (input$layout == "Single-end") {
      files <- input$fastq_single
      if (is.null(files) || nrow(files) == 0) {
        append_log("No single-end FASTQ uploaded. Aborting.")
        return()
      }
      # Save files
      for (i in seq_len(nrow(files))) {
        dest <- file.path(rv$outdir, files$name[i])
        file.copy(files$datapath[i], dest)
      }
      append_log(paste0("Saved ", nrow(files), " single-end file(s) to output directory."))
      # Path list for downstream
      fnFs <- file.path(rv$outdir, files$name)
      fnRs <- NULL
      
    } else { # Paired-end
      filesF <- input$fastqF
      filesR <- input$fastqR
      if (is.null(filesF) || is.null(filesR) || nrow(filesF) == 0 || nrow(filesR) == 0) {
        append_log("Forward and/or reverse files are missing. Please upload both R1 and R2 files.")
        return()
      }
      if (nrow(filesF) != nrow(filesR)) {
        append_log("Number of forward and reverse files does not match. Please upload matching pairs in the same order.")
        return()
      }
      # Save files (keep order)
      for (i in seq_len(nrow(filesF))) {
        file.copy(filesF$datapath[i], file.path(rv$outdir, filesF$name[i]))
        file.copy(filesR$datapath[i], file.path(rv$outdir, filesR$name[i]))
      }
      append_log(paste0("Saved ", nrow(filesF), " forward and ", nrow(filesR), " reverse files to output directory."))
      fnFs <- file.path(rv$outdir, filesF$name)
      fnRs <- file.path(rv$outdir, filesR$name)
    }
    
    # Start selected mode
    if (input$mode == "DADA2 + SILVA") {
      append_log("Starting DADA2 pipeline (single-core).")
      
      tryCatch({
        # load packages
        if (!requireNamespace("dada2", quietly = TRUE)) stop("Package 'dada2' not installed. Run setup.R to install dependencies.")
        if (!requireNamespace("ShortRead", quietly = TRUE)) stop("Package 'ShortRead' not installed.")
        
        library(dada2)
        library(ShortRead)
        
        # Log files to process
        append_log(paste("Files to process (F):", paste(basename(fnFs), collapse = ", ")))
        if (!is.null(fnRs)) append_log(paste("Files to process (R):", paste(basename(fnRs), collapse = ", ")))
        
        # Make filtered directory and set filtered file paths
        filt_dir <- file.path(rv$outdir, "filtered")
        dir.create(filt_dir, showWarnings = FALSE)
        filtFs <- file.path(filt_dir, basename(fnFs))
        if (!is.null(fnRs)) filtRs <- file.path(filt_dir, basename(fnRs))
        
        append_log("Running filterAndTrim...")
        if (is.null(fnRs)) {
          out <- filterAndTrim(fnFs, filtFs, truncQ = input$truncQ, maxEE = input$maxEE, multithread = FALSE)
        } else {
          out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                               truncLen = c(input$truncLenF, input$truncLenR),
                               maxEE = c(input$maxEE, input$maxEE),
                               truncQ = input$truncQ,
                               multithread = FALSE)
        }
        append_log(paste("filterAndTrim completed. Summary:\n", paste(capture.output(print(out)), collapse = "\n")))
        
        append_log("Learning error rates...")
        errF <- learnErrors(filtFs, multithread = FALSE)
        if (!is.null(fnRs)) errR <- learnErrors(filtRs, multithread = FALSE)
        append_log("Error models learned.")
        
        append_log("Running dada() inference...")
        if (is.null(fnRs)) {
          dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
          seqtab <- makeSequenceTable(dadaFs)
        } else {
          dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
          dadaRs <- dada(filtRs, err = errR, multithread = FALSE)
          mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
          seqtab <- makeSequenceTable(mergers)
        }
        append_log("Sequence table constructed.")
        
        if (input$remove_chimera) {
          append_log("Removing chimeras (consensus)...")
          seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE)
        } else {
          seqtab.nochim <- seqtab
        }
        
        # Assign taxonomy
        silva_path <- input$silva_path
        if (is.null(silva_path) || silva_path == "" || !file.exists(silva_path)) {
          append_log("SILVA FASTA not provided or file not found. Please provide a valid path in the UI.")
          stop("SILVA not provided")
        }
        append_log(paste("Assigning taxonomy using SILVA file:", silva_path))
        taxa <- assignTaxonomy(seqtab.nochim, silva_path, multithread = FALSE)
        append_log("Taxonomy assignment complete.")
        
        # Save outputs
        asv_seqs <- colnames(seqtab.nochim)
        asv_headers <- paste0(">ASV", seq_along(asv_seqs))
        asv_fasta <- file.path(rv$outdir, "asvs.fasta")
        cat(paste(rbind(asv_headers, asv_seqs), collapse = "\n"), file = asv_fasta)
        
        asv_table <- t(seqtab.nochim)
        colnames(asv_table) <- paste0("Sample_", seq_len(ncol(asv_table)))
        taxa_df <- as.data.frame(taxa, stringsAsFactors = FALSE)
        taxa_df$ASV <- rownames(taxa_df)
        
        write.csv(asv_table, file = file.path(rv$outdir, "asv_table.csv"), row.names = TRUE)
        write.csv(taxa_df, file = file.path(rv$outdir, "taxonomy.csv"), row.names = FALSE)
        
        rv$dada_tax <- taxa_df
        append_log("DADA2 pipeline finished. Outputs saved.")
        
      }, error = function(e) {
        append_log(paste("Error in DADA2 pipeline:", e$message))
      })
      
    } else if (input$mode == "Kraken2 (fast)") {
      append_log("Starting Kraken2 classification (single-core).")
      
      tryCatch({
        # Find fastq files in outdir (we saved files earlier)
        fq <- list.files(rv$outdir, pattern = "(fastq|fq)(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
        if (length(fq) == 0) stop("No FASTQ files found in upload directory to run Kraken2.")
        
        kraken_db <- input$kraken_db
        if (is.null(kraken_db) || kraken_db == "" || !dir.exists(kraken_db)) stop("Kraken2 DB path not provided or invalid.")
        
        out_report <- file.path(rv$outdir, "kraken_report.txt")
        out_labels <- file.path(rv$outdir, "kraken_labels.txt")
        
        for (f in fq) {
          append_log(paste("Running kraken2 on", basename(f)))
          cmd <- sprintf("kraken2 --db '%s' --report '%s' --output '%s' '%s'", kraken_db, out_report, out_labels, f)
          ret <- system(cmd)
          append_log(paste("kraken2 exit status:", ret))
        }
        
        if (file.exists(out_report)) {
          rpt <- tryCatch(read.table(out_report, sep = "\t", quote = "", fill = TRUE, stringsAsFactors = FALSE),
                          error = function(e) NULL)
          rv$kraken <- rpt
          if (!is.null(rpt)) write.csv(rpt, file.path(rv$outdir, "kraken_report.csv"), row.names = FALSE)
        }
        append_log("Kraken2 run finished. Outputs saved.")
        
      }, error = function(e) {
        append_log(paste("Error in Kraken2 pipeline:", e$message))
      })
    } else {
      append_log("Unknown mode selected.")
    }
    
    # Create downloadable zip of outputs
    tryCatch({
      files_to_zip <- list.files(rv$outdir, full.names = TRUE, recursive = TRUE)
      rv$zipfile <- file.path(rv$outdir, "results.zip")
      if (length(files_to_zip) > 0) {
        if (requireNamespace("zip", quietly = TRUE)) {
          zip::zip(rv$zipfile, files_to_zip)
        } else {
          oldwd <- getwd()
          setwd(rv$outdir)
          on.exit(setwd(oldwd), add = TRUE)
          utils::zip(basename(rv$zipfile), list.files(".", recursive = TRUE))
        }
        append_log(paste("Created results zip at", rv$zipfile))
      } else {
        append_log("No output files found to zip.")
      }
    }, error = function(e) {
      append_log(paste("Error creating zip:", e$message))
    })
    
    # Update UI outputs
    output$download_ui <- renderUI({
      if (!is.null(rv$zipfile) && file.exists(rv$zipfile)) {
        tagList(downloadButton("downloadZip", "Download results.zip"))
      } else {
        tags$span("No results available for download yet.")
      }
    })
    
    output$dada_tax <- DT::renderDataTable({
      if (is.null(rv$dada_tax)) return(NULL)
      DT::datatable(rv$dada_tax, options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$kraken_tbl <- DT::renderDataTable({
      if (is.null(rv$kraken)) return(NULL)
      DT::datatable(rv$kraken, options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$log <- renderText({ rv$log })
    
    output$downloadZip <- downloadHandler(
      filename = function() {
        paste0("dada2_kraken_results_", Sys.Date(), ".zip")
      },
      content = function(file) {
        file.copy(rv$zipfile, file)
      },
      contentType = "application/zip"
    )
  }) # end run observer
  
} # end server

# Run the app
shinyApp(ui, server)