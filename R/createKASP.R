#'Create KASP Function
#'This allows users to create KASP seqeunces using a reference sequence and variant information
#'Currently applicable for single SNP and INDELS
#'Requires phsyical positions of flanking sequence start and end, the physical position of the variant, the variant type, the variant itself and the reference sequence
#'All input data to be stored as R dataframe object
#'Computes simple stats and checks on KASP sequences inlcuding flank lengths and GC content. 

createKASP <- function(data,
                       f.start.pos = "f.start.pos",
                       f.end.pos = "f.end.pos",
                       var.type = "var.type",
                       var.pos = "var.pos",
                       ref.seq = "ref.seq",
                       var = "var",
                       seq.name = NULL) {
  
  # ---- 1. Package check ----
  if (!requireNamespace("stringr", quietly = TRUE)) {
    install.packages("stringr")
  }
  library(stringr)
  
  # ---- 2. Column checks ----
  required_cols <- c(f.start.pos, f.end.pos, var.type, var.pos, ref.seq, var)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # ---- 3. Variant type check ----
  allowed_types <- c("snp", "insertion", "deletion")
  bad_types <- setdiff(tolower(unique(data[[var.type]])), allowed_types)
  if (length(bad_types) > 0) {
    stop("Invalid var.type(s): ", paste(bad_types, collapse = ", "),
         "\nAllowed types: SNP, insertion, deletion")
  }
  
  # ---- 4. Assign seq.name if not given ----
  if (is.null(seq.name) || !(seq.name %in% names(data))) {
    data$seq.name <- paste0("kasp.seq.", seq_len(nrow(data)))
    seq.name <- "seq.name"
  }
  
  # ---- 5. Processing ----
  output_list <- list()
  message("Starting KASP sequence generation for ", nrow(data), " variants...")
  
  for (i in seq_len(nrow(data))) {
    row <- data[i, ]
    message("Processing ", row[[seq.name]], " (", row[[var.type]], ")")
    
    # Calculate relative position of variant within the reference sequence
    rel_var_pos <- row[[var.pos]] - row[[f.start.pos]] + 1
    
    # Extract flanking sequences based on variant type
    if (tolower(row[[var.type]]) == "snp") {
      # For SNPs: split at the variant position, replace the base
      left_flank <- substr(row[[ref.seq]], 1, rel_var_pos - 1)
      right_flank <- substr(row[[ref.seq]], rel_var_pos + 1, nchar(row[[ref.seq]]))
      ref_base <- substr(row[[ref.seq]], rel_var_pos, rel_var_pos)
      
      kasp_seq <- paste0(left_flank, "[", ref_base, "/", row[[var]], "]", right_flank)
      alt_seq  <- paste0(left_flank, row[[var]], right_flank)
      
    } else if (tolower(row[[var.type]]) == "insertion") {
      # For insertions: insert the new base AFTER the reference position
      # The variant position indicates where to insert (after this position)
      left_flank <- substr(row[[ref.seq]], 1, rel_var_pos)
      right_flank <- substr(row[[ref.seq]], rel_var_pos + 1, nchar(row[[ref.seq]]))
      
      # KASP notation: [inserted_base/-] means insert vs no insert
      kasp_seq <- paste0(left_flank, "[", row[[var]], "/-]", right_flank)
      alt_seq  <- paste0(left_flank, row[[var]], right_flank)
      
    } else if (tolower(row[[var.type]]) == "deletion") {
      # For deletions: remove the base at the variant position
      left_flank <- substr(row[[ref.seq]], 1, rel_var_pos - 1)
      right_flank <- substr(row[[ref.seq]], rel_var_pos + 1, nchar(row[[ref.seq]]))
      ref_base <- substr(row[[ref.seq]], rel_var_pos, rel_var_pos)
      
      # KASP notation: [deleted_base/-] means delete vs keep
      kasp_seq <- paste0(left_flank, "[", ref_base, "/-]", right_flank)
      alt_seq  <- paste0(left_flank, right_flank)  # Base is deleted
    }
    
    # Calculate GC content for both reference and alternative sequences
    # Count G and C bases (case insensitive) and calculate percentage
    gc_ref <- round(100 * str_count(row[[ref.seq]], "[GCgc]") / nchar(row[[ref.seq]]), 2)
    gc_alt <- round(100 * str_count(alt_seq, "[GCgc]") / nchar(alt_seq), 2)
    
    # Calculate flank lengths for quality assessment
    left_len <- nchar(left_flank)
    right_len <- nchar(right_flank)
    
    # Store results for this variant
    output_list[[i]] <- data.frame(
      seq.name = row[[seq.name]],
      kasp.seq = kasp_seq,
      ref.seq.full = row[[ref.seq]],
      alt.seq.full = alt_seq,
      total.length.ref = nchar(row[[ref.seq]]),
      total.length.alt = nchar(alt_seq),
      left.flank.length = left_len,
      right.flank.length = right_len,
      GC.content.ref = gc_ref,
      GC.content.alt = gc_alt,
      stringsAsFactors = FALSE
    )
    
    # Quality warnings for short flanks (recommended minimum for primer design)
    if (left_len < 50) warning("Left flank <50bp for ", row[[seq.name]], " (", left_len, "bp)")
    if (right_len < 50) warning("Right flank <50bp for ", row[[seq.name]], " (", right_len, "bp)")
    
    # Additional quality checks
    total_len <- left_len + right_len + 1  # +1 for variant position
    if (total_len > 300) warning("Very long sequence (", total_len, "bp) for ", row[[seq.name]], " - may affect PCR efficiency")
    if (gc_ref < 30 || gc_ref > 70) warning("GC content outside optimal range (30-70%) for ", row[[seq.name]], ": ", gc_ref, "%")
  }
  
  # Combine original data with generated KASP sequences
  output_df <- cbind(data, do.call(rbind, output_list))
  
  # ---- 6. Comprehensive QC checks ----
  message("\nRunning comprehensive QC validation...")
  
  # Debug: Check if KASP sequences were actually generated
  message("Debug: First few KASP sequences generated:")
  for (i in 1:min(3, nrow(output_df))) {
    message("  ", output_df[i, ][[seq.name]], ": ", substr(output_df$kasp.seq[i], 1, 100), "...")
  }
  
  qc_issues <- c()
  
  # Check that sequence length changes match variant type expectations
  insertion_issues <- output_df[tolower(output_df[[var.type]]) == "insertion" &
                                  output_df$total.length.alt != (output_df$total.length.ref + 1), ][[seq.name]]
  deletion_issues <- output_df[tolower(output_df[[var.type]]) == "deletion" &
                                 output_df$total.length.alt != (output_df$total.length.ref - 1), ][[seq.name]]
  snp_issues <- output_df[tolower(output_df[[var.type]]) == "snp" &
                            output_df$total.length.alt != output_df$total.length.ref, ][[seq.name]]
  
  # Check for KASP notation format issues in the generated sequences
  kasp_format_issues <- output_df[!grepl("\\[[A-Z/-]+\\]", output_df$kasp.seq), ][[seq.name]]
  
  # Integrated validation checks from helper function
  for (i in 1:nrow(output_df)) {
    row <- output_df[i, ]
    var_type <- tolower(row[[var.type]])
    seq_name <- row[[seq.name]]
    
    # Extract variant notation from KASP sequence
    # Pattern should match [A/T], [G/-], [T/-], etc.
    variant_match <- regmatches(row$kasp.seq, regexpr("\\[[A-Z/-]+\\]", row$kasp.seq))
    
    if (length(variant_match) == 0) {
      qc_issues <- c(qc_issues, paste("No variant notation found in", seq_name))
      next
    }
    
    # Check notation format matches variant type
    if (var_type == "snp" && grepl("/-", variant_match)) {
      qc_issues <- c(qc_issues, paste("SNP", seq_name, "has indel notation:", variant_match))
    }
    if (var_type %in% c("insertion", "deletion") && !grepl("/-", variant_match)) {
      qc_issues <- c(qc_issues, paste("Indel", seq_name, "missing proper indel notation:", variant_match))
    }
    
    # Check that variant base matches what's expected
    if (var_type == "snp") {
      # For SNPs, check if the variant allele is in the KASP notation
      alleles <- gsub("\\[|\\]", "", variant_match)
      allele_list <- unlist(strsplit(alleles, "/"))
      if (length(allele_list) == 2 && !row[[var]] %in% allele_list) {
        qc_issues <- c(qc_issues, paste("Variant allele", row[[var]], "not found in KASP notation for", seq_name))
      }
    }
    
    # Check for unusual flank imbalances (>70% difference)
    flank_ratio <- min(row$left.flank.length, row$right.flank.length) / max(row$left.flank.length, row$right.flank.length)
    if (flank_ratio < 0.3) {
      qc_issues <- c(qc_issues, paste("Highly unbalanced flanks for", seq_name, 
                                      paste0("(", row$left.flank.length, "bp vs ", row$right.flank.length, "bp)")))
    }
  }
  
  # Report all QC failures
  if (length(insertion_issues) > 0) {
    warning("Insertion length QC fail (alt should be ref+1): ", paste(insertion_issues, collapse = ", "))
  }
  if (length(deletion_issues) > 0) {
    warning("Deletion length QC fail (alt should be ref-1): ", paste(deletion_issues, collapse = ", "))
  }
  if (length(snp_issues) > 0) {
    warning("SNP length QC fail (alt should equal ref): ", paste(snp_issues, collapse = ", "))
  }
  if (length(kasp_format_issues) > 0) {
    warning("KASP notation format issues: ", paste(kasp_format_issues, collapse = ", "))
  }
  if (length(qc_issues) > 0) {
    warning("Additional validation issues found:")
    for (issue in qc_issues) {
      message("  - ", issue)
    }
  } else {
    message("✓ All validation checks passed!")
  }
  
  # ---- 7. Summary statistics and reporting ----
  type_counts <- table(tolower(output_df[[var.type]]))
  
  # Flank length summary
  flank_summary <- data.frame(
    flank.side = c("Left", "Right"),
    min.length = c(min(output_df$left.flank.length), min(output_df$right.flank.length)),
    max.length = c(max(output_df$left.flank.length), max(output_df$right.flank.length)),
    mean.length = c(round(mean(output_df$left.flank.length), 1), 
                    round(mean(output_df$right.flank.length), 1))
  )
  
  # GC content summary
  mean_gc_ref <- round(mean(output_df$GC.content.ref), 2)
  mean_gc_alt <- round(mean(output_df$GC.content.alt), 2)
  gc_range_ref <- paste0(round(min(output_df$GC.content.ref), 1), "-", 
                         round(max(output_df$GC.content.ref), 1), "%")
  gc_range_alt <- paste0(round(min(output_df$GC.content.alt), 1), "-", 
                         round(max(output_df$GC.content.alt), 1), "%")
  
  # Print comprehensive summary
  message("\n===== KASP Processing Summary =====")
  message("Total variants processed: ", nrow(output_df))
  message("\nVariant type breakdown:")
  for (t in names(type_counts)) {
    message("  ", toupper(t), ": ", type_counts[[t]])
  }
  message("\nFlanking sequence statistics (bp):")
  print(flank_summary, row.names = FALSE)
  message("\nGC content summary:")
  message("  Reference sequences: ", mean_gc_ref, "% (range: ", gc_range_ref, ")")
  message("  Alternative sequences: ", mean_gc_alt, "% (range: ", gc_range_alt, ")")
  
  # Quality assessment summary
  short_flanks <- sum(output_df$left.flank.length < 50 | output_df$right.flank.length < 50)
  long_seqs <- sum((output_df$left.flank.length + output_df$right.flank.length) > 300)
  gc_issues <- sum(output_df$GC.content.ref < 30 | output_df$GC.content.ref > 70)
  
  message("\nQuality assessment:")
  message("  Sequences with short flanks (<50bp): ", short_flanks)
  message("  Very long sequences (>300bp): ", long_seqs)
  message("  GC content outside 30-70% range: ", gc_issues)
  message("===================================\n")
  
  return(output_df)
}
