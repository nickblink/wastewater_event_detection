# Minimal base-R test runner (no testthat -- see R/00_install_packages.R for why).
# Run from the repo root: Rscript tests/run_tests.R
# Exits with status 1 if any test file fails, for use in CI.

test_files <- sort(Sys.glob("tests/test_*.R"))
failures <- character(0)

for (f in test_files) {
  cat(sprintf("---- Running %s ----\n", f))
  result <- tryCatch({
    source(f, local = new.env())
    "PASS"
  }, error = function(e) {
    cat(sprintf("FAIL: %s\n", conditionMessage(e)))
    "FAIL"
  })
  if (result == "FAIL") failures <- c(failures, f)
}

cat("\n==== Summary ====\n")
cat(sprintf("%d/%d test files passed\n", length(test_files) - length(failures), length(test_files)))

if (length(failures) > 0) {
  cat("Failed:", paste(failures, collapse = ", "), "\n")
  quit(status = 1)
}
