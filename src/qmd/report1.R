# ============================================================================
# COMPUTATIONAL BIOLOGY FOUNDATIONS (CSB195) - REPORT 1
# Standard Genetic Code Optimization Analysis
# Evaluates whether the SGC is optimized for tolerance to point mutations
# ============================================================================

set.seed(42)

# Load the official aaSim() function from course repository
baseURL <- "https://raw.githubusercontent.com/hyginn/CSB195/main/"
fn <- "./dat/aaSim.4.1.Rds"
download.file(paste0(baseURL, fn), fn)
aaSim <- readRDS("./dat/aaSim.4.1.Rds")

# ============================================================================
# 1. DEFINE THE STANDARD GENETIC CODE
# ============================================================================

SGC <- c(
  "UUU" = "F", "UUC" = "F", "UUA" = "L", "UUG" = "L",
  "UCU" = "S", "UCC" = "S", "UCA" = "S", "UCG" = "S",
  "UAU" = "Y", "UAC" = "Y", "UAA" = "*", "UAG" = "*",
  "UGU" = "C", "UGC" = "C", "UGA" = "*", "UGG" = "W",
  "CUU" = "L", "CUC" = "L", "CUA" = "L", "CUG" = "L",
  "CCU" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
  "CAU" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
  "CGU" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
  "AUU" = "I", "AUC" = "I", "AUA" = "I", "AUG" = "M",
  "ACU" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
  "AAU" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
  "AGU" = "S", "AGC" = "S", "AGA" = "R", "AGG" = "R",
  "GUU" = "V", "GUC" = "V", "GUA" = "V", "GUG" = "V",
  "GCU" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
  "GAU" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
  "GGU" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"
)

# ============================================================================
# 2. HELPER FUNCTIONS
# ============================================================================

# Generate all 9 single-nucleotide neighbours of a codon
generate_neighbours <- function(codon) {
  nucleotides <- c("U", "C", "A", "G")
  neighbours <- character(9)
  idx <- 1

  for (pos in 1:3) {
    for (nuc in nucleotides) {
      if (nuc != substr(codon, pos, pos)) {
        neighbours[idx] <- paste0(
          substr(codon, 1, pos - 1),
          nuc,
          substr(codon, pos + 1, 3)
        )
        idx <- idx + 1
      }
    }
  }
  return(neighbours)
}

# Compute mutation-tolerance score for a genetic code
compute_score <- function(code) {
  total_score <- 0

  for (codon in names(code)) {
    aa1 <- code[codon]
    neighbours <- generate_neighbours(codon)

    for (neighbour in neighbours) {
      aa2 <- code[neighbour]
      similarity <- aaSim(aa1, aa2)
      total_score <- total_score + similarity
    }
  }
  return(total_score)
}

# Generate a random valid genetic code
# Constraints: maps all 64 codons to 20 amino acids + stop; includes ≥1 stop codon
generate_random_code <- function() {
  # All 20 standard amino acids plus stop codon
  amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*")

  # All 64 codons
  codons <- names(SGC)

  # Randomly assign each codon to an amino acid
  assignments <- sample(amino_acids, size = 64, replace = TRUE)

  # Create the code vector
  code <- assignments
  names(code) <- codons

  # Ensure at least one stop codon (marked as "*")
  # If no stop codon exists, replace a random codon with "*"
  if (!("*" %in% code)) {
    idx <- sample(1:64, 1)
    code[idx] <- "*"
  }

  return(code)
}

# ============================================================================
# 3. VERIFY SGC BENCHMARK
# ============================================================================

sgc_score <- compute_score(SGC)
cat("Standard Genetic Code Score: ", sgc_score, "\n")
cat("Expected Score: 9856.116\n")
cat("Difference: ", abs(sgc_score - 9856.116), "\n\n")

# ============================================================================
# 4. GENERATE 1000 RANDOM VALID CODES AND COMPUTE SCORES
# ============================================================================

n_random <- 1000
random_scores <- numeric(n_random)

cat("Generating", n_random, "random codes...\n")
for (i in 1:n_random) {
  random_code <- generate_random_code()
  random_scores[i] <- compute_score(random_code)

  # Progress indicator every 100 iterations
  if (i %% 100 == 0) {
    cat("  Completed:", i, "/", n_random, "\n")
  }
}

# ============================================================================
# 5. COMPUTE STATISTICS
# ============================================================================

mean_random <- mean(random_scores)
sd_random <- sd(random_scores)
min_random <- min(random_scores)
max_random <- max(random_scores)

# Calculate percentile and z-score of SGC within the random distribution
percentile_sgc <- sum(random_scores <= sgc_score) / n_random * 100
z_score_sgc <- (sgc_score - mean_random) / sd_random

cat("\n=== SUMMARY STATISTICS ===\n")
cat("Standard Genetic Code Score:", sgc_score, "\n")
cat("Random Codes Mean:           ", mean_random, "\n")
cat("Random Codes SD:             ", sd_random, "\n")
cat("Random Codes Min:            ", min_random, "\n")
cat("Random Codes Max:            ", max_random, "\n")
cat("SGC Percentile:              ", percentile_sgc, "%\n")
cat("SGC Z-score:                 ", z_score_sgc, "\n")

# ============================================================================
# 6. VISUALIZATION: HISTOGRAM WITH SGC MARKED
# ============================================================================

png("mutation_tolerance_distribution.png", width = 800, height = 600)

hist(random_scores,
     breaks = 50,
     main = "Mutation-Tolerance Scores: Random vs Standard Genetic Code",
     xlab = "Quality Score",
     ylab = "Frequency",
     col = "lightblue",
     border = "gray30",
     cex.main = 1.3,
     cex.lab = 1.1)

# Add a vertical line marking the SGC score
abline(v = sgc_score, col = "red", lwd = 3, lty = 2)

# Add a legend
legend("topleft",
       legend = c(paste("SGC Score:", round(sgc_score, 2)),
                  paste("Mean Random:", round(mean_random, 2)),
                  paste("Percentile:", round(percentile_sgc, 1), "%")),
       col = c("red", "black", "black"),
       lty = c(2, 1, 0),
       lwd = c(3, 1, 1),
       cex = 1.0)

dev.off()

cat("\n✓ Histogram saved as 'mutation_tolerance_distribution.png'\n")

# ============================================================================
# 7. ADDITIONAL ANALYSIS TABLE
# ============================================================================

# Compare SGC with quantiles of random distribution
quantiles <- quantile(random_scores, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

cat("\n=== QUANTILE ANALYSIS ===\n")
cat("Random Code Quantiles:\n")
print(quantiles)
cat("\nSGC Score Position:\n")
cat("  SGC > 95th percentile?", sgc_score > quantiles["95%"], "\n")
cat("  SGC > 75th percentile?", sgc_score > quantiles["75%"], "\n")
cat("  SGC > median?", sgc_score > quantiles["50%"], "\n")

# ============================================================================
# END OF ANALYSIS
# ============================================================================
