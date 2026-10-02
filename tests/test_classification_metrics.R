# Unit tests for the sensitivity/PPV/F1 metrics in R/peak_classification_metrics.R.

source("R/peak_classification_metrics.R")

## Known TP/FP/FN counts -> hand-computed sensitivity/PPV/F1.
classified <- data.frame(class = c(rep("TP", 6), rep("FP", 3), rep("FN", 2)))

expected_sensitivity <- 6 / (6 + 2) # 0.75
expected_ppv <- 6 / (6 + 3)         # 0.6667
expected_f1 <- 2 * 6 / (2 * 6 + 3 + 2) # 0.8276

stopifnot(abs(sensitivity(classified) - expected_sensitivity) < 1e-9)
stopifnot(abs(PPV(classified) - expected_ppv) < 1e-9)
stopifnot(abs(F1(classified) - expected_f1) < 1e-9)

## Perfect classification -> sensitivity = PPV = F1 = 1.
perfect <- data.frame(class = rep("TP", 5))
stopifnot(sensitivity(perfect) == 1, PPV(perfect) == 1, F1(perfect) == 1)

cat("test_classification_metrics.R: all tests passed\n")
