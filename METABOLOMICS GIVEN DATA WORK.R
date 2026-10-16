############################################################
# MSc Bioinformatics - structToolbox Practical
# Adapted for: ralfdata.csv
# Same workflow as provided Rmd
############################################################

# -------------------------------
# 1) Load Libraries
# -------------------------------

library(structToolbox)
library(ggplot2)

# Install Bioconductor manager (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install pmp (required by structToolbox)
BiocManager::install("pmp")

# Reinstall structToolbox to be safe
BiocManager::install("structToolbox")

# -------------------------------
# 2) Import Data
# -------------------------------

# Import CSV file
X <- read.csv("ralfdata.csv",
              header = TRUE,
              check.names = FALSE,
              stringsAsFactors = FALSE)

# View structure
head(X)
dim(X)


# -------------------------------
# 3) Prepare Sample Meta & Peak Table
# -------------------------------

# Sample meta = first two columns
SM <- X[, 1:2]

# Peak table = remaining columns
PT <- X[, 3:ncol(X)]


# -------------------------------
# 4) Convert Label to Factor
# -------------------------------

SM$Label = factor(SM$Label)

# show counts
table(SM$Label)

SM <- X[, 1:2]
SM$Label <- factor(SM$Label)

# -------------------------------
# 3) Prepare Sample Meta & Peak Table (FIXED for your ralfdata.csv)
# -------------------------------

# Sample IDs are in the first column (named "mz" in your file, but it's sample names)
sample_id <- as.character(X[[1]])

# Create Label from sample_id (structToolbox practical expects a "Label" factor)
Label <- ifelse(grepl("QC", sample_id, ignore.case = TRUE), "QC",
                ifelse(grepl("_C", sample_id, ignore.case = TRUE), "C",
                       ifelse(grepl("_S", sample_id, ignore.case = TRUE), "S", "Sample")))

SM <- data.frame(
  id = sample_id,
  Label = factor(Label),
  stringsAsFactors = FALSE
)

# Peak table = all remaining columns
PT <- X[, -1, drop = FALSE]

# Convert to numeric (important)
PT[] <- lapply(PT, function(v) as.numeric(as.character(v)))

# (Optional but common) Treat zeros as missing values
PT[PT == 0] <- NA

# Variable metadata
VM <- data.frame(
  id = colnames(PT),
  stringsAsFactors = FALSE
)

# Quick checks
table(SM$Label)
dim(PT)



# -------------------------------
# 5) Prepare Variable Metadata 
# -------------------------------

VM <- data.frame(
  id = colnames(PT),
  stringsAsFactors = FALSE
)


# -------------------------------
# 6) Create DatasetExperiment Object
# -------------------------------

DE <-
  DatasetExperiment(
    name = "Workshop data",
    description = "Data for the MS-based metabolomics practical",
    data = PT,
    sample_meta = SM,
    variable_meta = VM
  )

# Show summary
show(DE)


# -------------------------------
# 7) Build PCA Workflow
# -------------------------------

WF <-
  mv_feature_filter(
    threshold = 50,      # Same as original (0.5 = 50%)
    qc_label = "QC",
    method = "across",
    factor_name = "Label"
  ) +
  rsd_filter(
    rsd_threshold = 30,
    qc_label = "QC",
    factor_name = "Label"
  ) +
  knn_impute(
    neighbours = 5,
    by = "features"
  ) +
  pqn_norm(
    qc_label = "QC",
    factor_name = "Label"
  ) +
  log_transform(
    base = 10
  ) +
  mean_centre() +
  PCA(
    number_components = 8
  )


# -------------------------------
# 8) Apply Workflow
# -------------------------------

WF <- model_apply(WF, DE)


# -------------------------------
# 9) PCA Charts
# -------------------------------

# Scree plot
C1 <- pca_scree_plot()

# Scores plot
C2 <- pca_scores_plot(
  factor_name = "Label",
  components = c(1, 2)
)


# Plot PCA (Step 7 = PCA)
chart_plot(C1, WF[7]) +
  theme(panel.background = element_rect(fill = "white"))

chart_plot(C2, WF[7]) +
  theme(panel.background = element_rect(fill = "white"))


# -------------------------------
# 10) Workflow Excluding QCs
# -------------------------------

WF2 <-
  mv_feature_filter(
    threshold = 50,
    qc_label = "QC",
    method = "across",
    factor_name = "Label"
  ) +
  rsd_filter(
    rsd_threshold = 30,
    qc_label = "QC",
    factor_name = "Label"
  ) +
  knn_impute(
    neighbours = 5,
    by = "features"
  ) +
  pqn_norm(
    qc_label = "QC",
    factor_name = "Label"
  ) +
  log_transform(
    base = 10
  ) +
  
  # Exclude QCs
  filter_smeta(
    mode = "exclude",
    levels = "QC",
    factor_name = "Label"
  ) +
  
  mean_centre() +
  PCA(
    number_components = 8
  )


# -------------------------------
# 11) Apply Second Workflow
# -------------------------------

WF2 <- model_apply(WF2, DE)


# -------------------------------
# 12) PCA After QC Removal
# -------------------------------

# PCA is now Step 8

chart_plot(C1, WF2[8]) +
  theme(panel.background = element_rect(fill = "white"))

chart_plot(C2, WF2[8]) +
  theme(panel.background = element_rect(fill = "white"))


# -------------------------------
# 13) View PCA Loadings
# -------------------------------

C3 <- pca_loadings_plot(
  components = c(1, 2)
)

chart_plot(C3, WF2[8]) +
  theme(panel.background = element_rect(fill = "white"))


############################################################
# END SCRIPT
############################################################

############################################################
# Demo-script adapted to YOUR data (ralfdata.csv)
# Same workflow as your demo:
# UNI = mv_feature_filter + rsd_filter + knn_impute + filter_smeta
# Threshold 0.5 => 50 (%)
############################################################

# -------------------------------
# 0) Packages
# -------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("pmp", quietly = TRUE)) BiocManager::install("pmp", ask = FALSE)
if (!requireNamespace("structToolbox", quietly = TRUE)) BiocManager::install("structToolbox", ask = FALSE)

library(structToolbox)

# -------------------------------
# 1) Read YOUR CSV
# -------------------------------
input_csv <- "ralfdata.csv"

X <- read.csv(input_csv, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# -------------------------------
# 2) Build SM / PT exactly for your file
#    (Your file: first column = sample names; no Label column provided)
# -------------------------------
sample_id <- as.character(X[[1]])

# Create Label the same way structToolbox practical expects (QC / C / S / Sample)
Label <- ifelse(grepl("QC", sample_id, ignore.case = TRUE), "QC",
                ifelse(grepl("_C", sample_id, ignore.case = TRUE), "C",
                       ifelse(grepl("_S", sample_id, ignore.case = TRUE), "S", "Sample")))

SM <- data.frame(
  id = sample_id,
  Label = factor(Label),
  stringsAsFactors = FALSE
)

# Peak table = everything except first column
PT <- X[, -1, drop = FALSE]

# Ensure numeric
PT[] <- lapply(PT, function(v) as.numeric(as.character(v)))

# Treat 0 as missing (common in MS exports)
PT[PT == 0] <- NA

# Row names = sample IDs
rownames(PT) <- SM$id

# Variable metadata (minimum needed)
VM <- data.frame(
  id = colnames(PT),
  stringsAsFactors = FALSE
)

# Quick checks
print(table(SM$Label))
print(dim(PT))

# -------------------------------
# 3) Create DatasetExperiment (DE)
# -------------------------------
DE <- DatasetExperiment(
  name = "ralfdata",
  description = "ralfdata.csv adapted to demo UNI workflow",
  data = PT,
  sample_meta = SM,
  variable_meta = VM
)

# -------------------------------
# 4) DEMO WORKFLOW (same as yours)
# -------------------------------
UNI =
  mv_feature_filter(
    threshold = 50,       # 0.5 => 50%
    qc_label = "QC",
    method = "across",
    factor_name = "Label"
  ) +
  rsd_filter(
    rsd_threshold = 30,
    qc_label = "QC",
    factor_name = "Label"
  ) +
  knn_impute(
    neighbours = 5,
    by = "features"
  ) +
  filter_smeta(
    mode = "exclude",
    levels = "QC",
    factor_name = "Label"
  )

# apply
UNI <- model_apply(UNI, DE)

DE <- DatasetExperiment(
  name = "ralfdata",
  description = "ralfdata.csv adapted to demo UNI workflow",
  data = PT,
  sample_meta = SM,
  variable_meta = VM
)

UNI <-
  mv_feature_filter(
    threshold = 50,
    qc_label = "QC",
    method = "across",
    factor_name = "Label"
  ) +
  rsd_filter(
    rsd_threshold = 30,
    qc_label = "QC",
    factor_name = "Label"
  ) +
  knn_impute(
    neighbours = 5,
    by = "features"
  ) +
  filter_smeta(
    mode = "exclude",
    levels = "QC",
    factor_name = "Label"
  )


# -------------------------------
# 5) Export the processed matrix (after QC removal)
# -------------------------------
last_step <- UNI[length(UNI)]
processed_X <- last_step$data

write.csv(processed_X, "ralfdata_UNI_processed.csv", row.names = TRUE)

message("✅ Done. Output saved as: ralfdata_UNI_processed.csv")

D_final <- predicted(UNI)        # DatasetExperiment after the whole workflow
processed_X <- D_final$data      # processed matrix
write.csv(processed_X, "ralfdata_UNI_processed.csv", row.names = TRUE)

# --------------------------------
# Get filtered DatasetExperiment
# --------------------------------
DEfilt <- predicted(UNI)

# Check available groups (important)
table(DEfilt$sample_meta$Label)


# --------------------------------
# Create Fold Change Model
# --------------------------------
FC <-
  fold_change(
    factor_name = "Label",
    paired = FALSE,
    control_group = "S",   # Change if needed (see table above)
    method = "mean"
  )


# --------------------------------
# Apply Fold Change Model
# --------------------------------
FC <- model_apply(FC, DEfilt)


# --------------------------------
# Extract Fold Change Results
# --------------------------------
fc_step <- FC[length(FC)]

# View structure (optional, for debugging)
str(fc_step)

# Most structToolbox versions store FC here:
fc_results <- fc_step$fold_change

# Preview
head(fc_results)

# Save to file
write.csv(fc_results,
          "ralfdata_fold_change_results.csv",
          row.names = FALSE)


message("✅ Fold change analysis complete. Results saved.")

dim(PT)
length(SM$id)
class(fc_step)
str(fc_step)

# Apply fold change model
FC <- model_apply(FC, DEfilt)

# Create fc_step (last step of the workflow output)
fc_step <- FC[length(FC)]

# Check what fc_step contains
str(fc_step)

# ---------- PRE-REQ: UNI must already be applied ----------
# UNI <- model_apply(UNI, DE)

# 1) Get filtered dataset
DEfilt <- predicted(UNI)

# Sanity check: make sure Label exists and has groups
stopifnot("Label" %in% colnames(DEfilt$sample_meta))
print(table(DEfilt$sample_meta$Label))
print(levels(DEfilt$sample_meta$Label))

# 2) Create fold change model (MAKE SURE control group exists!)
control_grp <- "S"
if (!(control_grp %in% levels(DEfilt$sample_meta$Label))) {
  stop(paste0("control_group '", control_grp, "' not found. Available: ",
              paste(levels(DEfilt$sample_meta$Label), collapse = ", ")))
}

FC_model <- fold_change(
  factor_name = "Label",
  paired = FALSE,
  control_group = control_grp,
  method = "mean"
)

# 3) Apply fold change model
FC <- model_apply(FC_model, DEfilt)

# 4) NOW create fc_step (this is the line you’re missing / not reaching)
stopifnot(exists("FC"))
fc_step <- FC[[length(FC)]]

# 5) Inspect and extract results (slot name depends on version)
str(fc_step)

# Try common slots safely
fc_results <- NULL
if (!is.null(fc_step$fold_change)) fc_results <- fc_step$fold_change
if (is.null(fc_results) && !is.null(fc_step$results)) fc_results <- fc_step$results
if (is.null(fc_results) && !is.null(fc_step$output)) fc_results <- fc_step$output

if (is.null(fc_results)) {
  stop("No fold-change table found inside fc_step. Look at str(fc_step) output to see the correct slot name.")
}

# 6) Save
head(fc_results)
write.csv(fc_results, "ralfdata_fold_change_results.csv", row.names = FALSE)
message("✅ Saved: ralfdata_fold_change_results.csv")


UNI <-
  mv_feature_filter(
    threshold = 50,
    qc_label = "QC",
    method = "across",
    factor_name = "Label"
  ) +
  rsd_filter(
    rsd_threshold = 30,
    qc_label = "QC",
    factor_name = "Label"
  ) +
  knn_impute(
    neighbours = 5,
    by = "features"
  ) +
  filter_smeta(
    mode = "exclude",
    levels = "QC",
    factor_name = "Label"
  )

# apply
UNI <- model_apply(UNI, DE)

DEfilt = predicted(UNI)

# create a fold change object
FC = 
  fold_change(
    factor_name = 'Label',
    paired = FALSE,
    control_group = 'S',
    method = 'mean')

# apply
FC = model_apply(FC,DEfilt)

fc = FC$fold_change

# log transform
fc = log(fc,base = 2)

# create an id column
fc$id=rownames(fc)

# create an index column
fc$idx = as.numeric(1:nrow(fc))

# create a important column
fc$important = 'Not important'
fc$important[fc$`C/S` < -log2(2)] = 'Important (Down)'
fc$important[fc$`C/S` > log2(2)] = 'Important (Up)'

# create ggplot
g = ggplot(
  data = fc,
  aes(x = idx,
      y = `C/S`,
      colour = important,
      fill = important)) +
  geom_point(
    alpha = 0.5,
    shape = 21,
    size = 3) +
  scale_color_manual(
    values=c(
      'Important (Down)' = 'blue',
      'Important (Up)' = 'red',
      'Not important'='grey')) +
  scale_fill_manual(
    values = c(
      'Significant (Down)' = 'blue',
      'Significant (Up)' = 'red',
      'Not significant' = 'grey')) +
  scale_x_continuous(
    limits = c(1,nrow(fc)),
    breaks = seq(0,nrow(fc),by = 100)) +
  theme(
    panel.background = element_rect(fill='white'),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black")) +
  ylab('Log2(FC)')+
  xlab('Feature index') 


imp = subset(fc,important != 'Not important')
head(imp)

fc = FC$fold_change

# log transform
fc = log(fc,base = 2)

# create an id column
fc$id=rownames(fc)

# create an index column
fc$idx = as.numeric(1:nrow(fc))

# create a important column
fc$important = 'Not important'
fc$important[fc$`C/S` < -log2(2)] = 'Important (Down)'
fc$important[fc$`C/S` > log2(2)] = 'Important (Up)'

# create ggplot
g = ggplot(
  data = fc,
  aes(x = idx,
      y = `C/S`,
      colour = important,
      fill = important)) +
  geom_point(
    alpha = 0.5,
    shape = 21,
    size = 3) +
  scale_color_manual(
    values=c(
      'Important (Down)' = 'blue',
      'Important (Up)' = 'red',
      'Not important'='grey')) +
  scale_fill_manual(
    values = c(
      'Significant (Down)' = 'blue',
      'Significant (Up)' = 'red',
      'Not significant' = 'grey')) +
  scale_x_continuous(
    limits = c(1,nrow(fc)),
    breaks = seq(0,nrow(fc),by = 100)) +
  theme(
    panel.background = element_rect(fill='white'),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black")) +
  ylab('Log2(FC)')+
  xlab('Feature index') 
g


library(ggplot2)

# -----------------------------
# 1) Get fold-change results safely
# -----------------------------
# FC is the output of: FC <- model_apply(FC_model, DEfilt)

fc_step <- FC[[length(FC)]]

# Try common slot names (depends on structToolbox version)
fc_tbl <- NULL
if (!is.null(fc_step$fold_change)) fc_tbl <- fc_step$fold_change
if (is.null(fc_tbl) && !is.null(fc_step$results)) fc_tbl <- fc_step$results
if (is.null(fc_tbl) && !is.null(fc_step$output)) fc_tbl <- fc_step$output

if (is.null(fc_tbl)) {
  stop("Couldn't find a fold-change table in FC. Run: str(FC[[length(FC)]]) and check slot names.")
}

# Make sure it's a data.frame
fc <- as.data.frame(fc_tbl, check.names = FALSE)

# Add feature id (from rownames if present; otherwise make one)
if (is.null(rownames(fc)) || all(rownames(fc) == "")) {
  fc$id <- paste0("feature_", seq_len(nrow(fc)))
} else {
  fc$id <- rownames(fc)
}

# Index column
fc$idx <- seq_len(nrow(fc))


# -----------------------------
# 2) Identify the fold-change column in YOUR results
# -----------------------------
# Prefer a column like "C/S" (your demo), otherwise guess common names
candidate_cols <- c("C/S", "S/C", "C_vs_S", "S_vs_C", "FC", "fold_change", "log2FC", "logFC")
fc_col <- candidate_cols[candidate_cols %in% colnames(fc)][1]

if (is.na(fc_col) || is.null(fc_col)) {
  stop(
    paste0(
      "Couldn't find a fold-change column. Available columns are:\n",
      paste(colnames(fc), collapse = ", ")
    )
  )
}

# -----------------------------
# 3) Convert to numeric + log2 transform if needed
# -----------------------------
fc[[fc_col]] <- as.numeric(as.character(fc[[fc_col]]))

# If the column is NOT already log2 scale, convert to log2.
# Heuristic: if values are mostly positive and not centered around 0, it's likely a ratio FC.
# If it already looks like log2FC (both positive and negative), keep as-is.
looks_like_log <- any(fc[[fc_col]] < 0, na.rm = TRUE)

if (!looks_like_log) {
  # Avoid log2(0) or negative
  fc[[fc_col]][fc[[fc_col]] <= 0] <- NA
  fc$log2FC <- log2(fc[[fc_col]])
} else {
  fc$log2FC <- fc[[fc_col]]
}

# -----------------------------
# 4) Mark important features (same thresholds as your demo)
#    log2(2) = 1; -log2(2) = -1
# -----------------------------
fc$important <- "Not important"
fc$important[fc$log2FC < -1] <- "Important (Down)"
fc$important[fc$log2FC >  1] <- "Important (Up)"

# -----------------------------
# 5) Plot (same style as your demo, fixed labels)
# -----------------------------
g <- ggplot(
  data = fc,
  aes(x = idx, y = log2FC, colour = important, fill = important)
) +
  geom_point(alpha = 0.5, shape = 21, size = 3) +
  scale_color_manual(values = c(
    "Important (Down)" = "blue",
    "Important (Up)"   = "red",
    "Not important"    = "grey"
  )) +
  scale_fill_manual(values = c(
    "Important (Down)" = "blue",
    "Important (Up)"   = "red",
    "Not important"    = "grey"
  )) +
  scale_x_continuous(
    limits = c(1, nrow(fc)),
    breaks = seq(0, nrow(fc), by = 100)
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ylab("Log2(FC)") +
  xlab("Feature index")

g

DEfilt = predicted(UNI)

# create a fold change object
TS = 
  wilcox_test(
    factor_names = 'Label',
    paired = FALSE,
    mtc = 'none'
  )

# apply
TS = model_apply(TS,DEfilt)

pv = TS$p_value

# apply FDR
pv$fdr = p.adjust(pv$p_value,method = 'fdr')

# log transform
pv$log_p_value = -log(pv$p_value,base = 10)

# create an id column
pv$id=rownames(pv)

# create an index column
pv$idx = as.numeric(1:nrow(pv))

# create a important column
pv$significant = 'Not significant'
pv$significant[pv$fdr < 0.05] = 'Significant'

# create ggplot
g = 
  ggplot(
    data = pv,
    aes(x = idx,
        y = log_p_value,
        colour = significant,
        fill = significant)) +
  geom_point(
    alpha = 0.5,
    shape = 21,
    size =3) +
  scale_color_manual(
    values = c(
      'Significant' = 'green',
      'Not significant' = 'lightgrey'))+
  scale_color_manual(
    values = c(
      'Significant' = 'purple',
      'Not significant' = 'yellow'))+
  scale_x_continuous(
    limits = c(1,nrow(pv)),
    breaks = seq(0,nrow(pv),by = 100)) +
  theme(
    panel.background = element_rect(fill='black'),
    panel.grid.major = element_line(color = "blue"),
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "red")) +
  ylab('-log10(p-value)')+
  xlab('Feature index') 

library(ggplot2)

# ----------------------------
# 1) Extract p-values (S4)
# ----------------------------
pv <- as.data.frame(TS@p_value, check.names = FALSE)

# If the p-value column isn't literally named "p_value", detect it
pcol <- c("p_value", "p.value", "pval", "p", "P.Value")[c("p_value","p.value","pval","p","P.Value") %in% colnames(pv)][1]
if (is.na(pcol)) {
  stop(paste0("No p-value column found. Available columns: ", paste(colnames(pv), collapse = ", ")))
}

# Make numeric
pv[[pcol]] <- as.numeric(as.character(pv[[pcol]]))

# ----------------------------
# 2) FDR + -log10(p)
# ----------------------------
pv$fdr <- p.adjust(pv[[pcol]], method = "fdr")

# avoid -log10(0) = Inf
pv[[pcol]][pv[[pcol]] <= 0] <- .Machine$double.xmin
pv$log_p_value <- -log10(pv[[pcol]])

# ----------------------------
# 3) id + index + significance
# ----------------------------
pv$id  <- rownames(pv)
pv$idx <- seq_len(nrow(pv))

pv$significant <- "Not significant"
pv$significant[pv$fdr < 0.05] <- "Significant"

# ----------------------------
# 4) Plot
# ----------------------------
g <-
  ggplot(
    data = pv,
    aes(
      x = idx,
      y = log_p_value,
      colour = significant,
      fill = significant
    )
  ) +
  geom_point(alpha = 0.5, shape = 21, size = 3) +
  scale_color_manual(values = c(
    "Significant" = "purple",
    "Not significant" = "lightgrey"
  )) +
  scale_fill_manual(values = c(
    "Significant" = "purple",
    "Not significant" = "lightgrey"
  )) +
  scale_x_continuous(
    limits = c(1, nrow(pv)),
    breaks = seq(0, nrow(pv), by = 100)
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ylab("-log10(p-value)") +
  xlab("Feature index")

g


library(ggplot2)

# ----------------------------
# 1) Extract p-values (S4 object)
# ----------------------------
pv <- as.data.frame(TS@p_value, check.names = FALSE)

# Check column names (for debugging)
print(colnames(pv))

library(ggplot2)

# --- helper: convert struct::entity (or similar S4) -> data.frame ---
entity_to_df <- function(x) {
  # If it's already a data.frame
  if (is.data.frame(x)) return(x)
  
  # If it's a matrix
  if (is.matrix(x)) return(as.data.frame(x, check.names = FALSE))
  
  # If it's an S4 object (like struct::entity)
  if (isS4(x)) {
    sn <- slotNames(x)
    
    # Most common: @data
    if ("data" %in% sn) {
      return(as.data.frame(slot(x, "data"), check.names = FALSE))
    }
    
    # Sometimes: @values
    if ("values" %in% sn) {
      return(as.data.frame(slot(x, "values"), check.names = FALSE))
    }
    
    # Sometimes: nested object that contains data
    for (nm in sn) {
      obj <- try(slot(x, nm), silent = TRUE)
      if (!inherits(obj, "try-error")) {
        if (is.data.frame(obj)) return(obj)
        if (is.matrix(obj)) return(as.data.frame(obj, check.names = FALSE))
        if (isS4(obj) && "data" %in% slotNames(obj)) {
          return(as.data.frame(slot(obj, "data"), check.names = FALSE))
        }
      }
    }
  }
  
  stop("Couldn't extract a data table from this object. Use: str(TS@p_value) and slotNames(TS@p_value).")
}

# ----------------------------
# 1) Extract p-values table
# ----------------------------
pv <- entity_to_df(TS@p_value)

# quick look
print(colnames(pv))


# ----------------------------
# 2) Find p-value column safely
# ----------------------------
possible_cols <- c("p_value", "p.value", "pval", "p", "P.Value")

pcol <- possible_cols[possible_cols %in% colnames(pv)][1]

# Stop if not found
if (is.na(pcol) || is.null(pcol)) {
  stop(
    paste0(
      "No p-value column found. Available columns: ",
      paste(colnames(pv), collapse = ", ")
    )
  )
}

# Make numeric
pv[[pcol]] <- as.numeric(as.character(pv[[pcol]]))

# ----------------------------
# 3) FDR + -log10(p)
# ----------------------------
pv$fdr <- p.adjust(pv[[pcol]], method = "fdr")

# Avoid log10(0)
pv[[pcol]][pv[[pcol]] <= 0] <- .Machine$double.xmin

pv$log_p_value <- -log10(pv[[pcol]])

# ----------------------------
# 4) id + index + significance
# ----------------------------
pv$id  <- rownames(pv)
pv$idx <- seq_len(nrow(pv))

pv$significant <- "Not significant"
pv$significant[pv$fdr < 0.05] <- "Significant"

# ----------------------------
# 5) Plot
# ----------------------------
g <-
  ggplot(
    pv,
    aes(
      x = idx,
      y = log_p_value,
      colour = significant,
      fill = significant
    )
  ) +
  geom_point(alpha = 0.5, shape = 21, size = 3) +
  scale_color_manual(values = c(
    "Significant" = "purple",
    "Not significant" = "lightgrey"
  )) +
  scale_fill_manual(values = c(
    "Significant" = "purple",
    "Not significant" = "lightgrey"
  )) +
  scale_x_continuous(
    limits = c(1, nrow(pv)),
    breaks = seq(0, nrow(pv), by = 100)
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_blank(),
    axis.line = element
    
    
    g = 
      ggplot(
        data = pv,
        aes(x = idx,
            y = log_p_value,
            colour = significant,
            fill = significant)) +
      geom_point(
        alpha = 0.5,
        shape = 21,
        size =3) +
      scale_color_manual(
        values = c(
          'Significant' = 'purple',
          'Not significant' = 'lightgrey'))+
      scale_color_manual(
        values = c(
          'Significant' = 'purple',
          'Not significant' = 'lightgrey'))+
      scale_x_continuous(
        limits = c(1,nrow(pv)),
        breaks = seq(0,nrow(pv),by = 100)) +
      theme(
        panel.background = element_rect(fill='white'),
        panel.grid.major = element_line(color = "lightgrey"),
        panel.grid.minor = element_blank(),  
        axis.line = element_line(color = "black")) +
      ylab('-log10(p-value)')+
      xlab('Feature index') 
    
    
    g <-
      ggplot(
        data = pv,
        aes(
          x = idx,
          y = log_p_value,
          colour = significant,
          fill = significant
        )
      ) +
      geom_point(
        alpha = 0.5,
        shape = 21,
        size = 3
      ) +
      scale_color_manual(
        values = c(
          "Significant" = "purple",
          "Not significant" = "lightgrey"
        )
      ) +
      scale_fill_manual(
        values = c(
          "Significant" = "purple",
          "Not significant" = "lightgrey"
        )
      ) +
      scale_x_continuous(
        limits = c(1, nrow(pv)),
        breaks = seq(0, nrow(pv), by = 100)
      ) +
      theme(
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "lightgrey"),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
      ) +
      ylab("-log10(p-value)") +
      xlab("Feature index")
    
    # Show plot
    g
    
    
    
    g <-
      ggplot(
        data = pv,
        aes(
          x = idx,
          y = log_p_value,
          colour = significant,
          fill = significant
        )
      ) +
      geom_point(
        alpha = 0.5,
        shape = 21,
        size = 3
      ) +
      scale_color_manual(
        values = c(
          "Significant" = "purple",
          "Not significant" = "lightgrey"
        )
      ) +
      scale_fill_manual(
        values = c(
          "Significant" = "purple",
          "Not significant" = "lightgrey"
        )
      ) +
      scale_x_continuous(
        limits = c(1, nrow(pv)),
        breaks = seq(0, nrow(pv), by = 100)
      ) +
      theme(
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "lightgrey"),
        panel.grid.minor = element_blank(),
        axis.lin
        
    
        
        
        pv = TS$p_value
        
        # apply FDR
        pv$fdr = p.adjust(pv$p_value,method = 'fdr')
        
        # log transform
        pv$log_p_value = -log(pv$p_value,base = 10)
        
        # create an id column
        pv$id=rownames(pv)
        
        # create an index column
        pv$idx = as.numeric(1:nrow(pv))
        
        # create a important column
        pv$significant = 'Not significant'
        pv$significant[pv$fdr < 0.05] = 'Significant'
        
        # create ggplot
        g = 
          ggplot(
            data = pv,
            aes(x = idx,
                y = log_p_value,
                colour = significant,
                fill = significant)) +
          geom_point(
            alpha = 0.5,
            shape = 21,
            size =3) +
          scale_color_manual(
            values = c(
              'Significant' = 'purple',
              'Not significant' = 'lightgrey'))+
          scale_color_manual(
            values = c(
              'Significant' = 'purple',
              'Not significant' = 'lightgrey'))+
          scale_x_continuous(
            limits = c(1,nrow(pv)),
            breaks = seq(0,nrow(pv),by = 100)) +
          theme(
            panel.background = element_rect(fill='white'),
            panel.grid.major = element_line(color = "lightgrey"),
            panel.grid.minor = element_blank(),  
            axis.line = element_line(color = "black")) +
          ylab('-log10(p-value)')+
          xlab('Feature index') 



