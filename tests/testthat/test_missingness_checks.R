# Unit tests for missingness_checks functions
# These are general-purpose missing data detection functions
# For gene alias handling, see na_check() in test-na_check.R

# Test detected_missingness() ------------------------------------------

test_that("detected_missingness() handles NULL input", {
  expect_equal(detected_missingness(NULL), TRUE)
})

test_that("detected_missingness() handles empty vector", {
  expect_equal(detected_missingness(character(0)), TRUE)
  expect_equal(detected_missingness(numeric(0)), TRUE)
})

test_that("detected_missingness() detects NA values", {
  x <- c("value", NA, "another")
  result <- detected_missingness(x)
  expect_equal(result, c(FALSE, TRUE, FALSE))
})

test_that("detected_missingness() detects empty strings", {
  x <- c("value", "", "another")
  result <- detected_missingness(x)
  expect_equal(result, c(FALSE, TRUE, FALSE))
})

test_that("detected_missingness() detects whitespace-only strings", {
  x <- c("value", "  ", "   ", "\t", "\n", "another")
  result <- detected_missingness(x)
  expect_equal(result, c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE))
})

test_that("detected_missingness() handles mixed cases", {
  x <- c("human", "mouse", NA, "", "  ", "rat")
  result <- detected_missingness(x)
  expect_equal(result, c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE))
})

test_that("detected_missingness() treats string 'NA' as valid data", {
  # This is the key difference from na_check()
  expect_false(detected_missingness("NA"))
  
  x <- c("a", NA, "", "b", "NA", "NANA", "na")
  result <- detected_missingness(x)
  expect_equal(result, c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
  
  # String "NA" should be kept
  filtered <- x[!detected_missingness(x)]
  expect_equal(filtered, c("a", "b", "NA", "NANA", "na"))
})

test_that("detected_missingness() handles numeric vectors", {
  x <- c(1, 2, NA, 3)
  result <- detected_missingness(x)
  expect_equal(result, c(FALSE, FALSE, TRUE, FALSE))
})

test_that("detected_missingness() handles factor vectors", {
  x <- factor(c("human", "mouse", NA, "rat"))
  result <- detected_missingness(x)
  expect_equal(result, c(FALSE, FALSE, TRUE, FALSE))
})


# Test any_missingness() -----------------------------------------------

test_that("any_missingness() returns TRUE when NULL", {
  expect_true(any_missingness(NULL))
})

test_that("any_missingness() returns TRUE when empty vector", {
  expect_true(any_missingness(character(0)))
})

test_that("any_missingness() returns TRUE when any NA present", {
  expect_true(any_missingness(c("value", NA, "another")))
})

test_that("any_missingness() returns TRUE when any empty string present", {
  expect_true(any_missingness(c("value", "", "another")))
})

test_that("any_missingness() returns TRUE when any whitespace present", {
  expect_true(any_missingness(c("value", "  ", "another")))
})

test_that("any_missingness() returns FALSE when all values valid", {
  expect_false(any_missingness(c("human", "mouse", "rat")))
})

test_that("any_missingness() treats string 'NA' as valid", {
  expect_false(any_missingness(c("NA", "NANA", "na")))
})

test_that("any_missingness() returns TRUE when all values missing", {
  expect_true(any_missingness(c(NA, "", "  ")))
})


# Test all_truthy() ----------------------------------------------------

test_that("all_truthy() returns FALSE for NULL", {
  expect_false(all_truthy(NULL))
})

test_that("all_truthy() returns FALSE for empty vector", {
  expect_false(all_truthy(character(0)))
})

test_that("all_truthy() returns FALSE when any NA present", {
  expect_false(all_truthy(c("human", NA, "rat")))
})

test_that("all_truthy() returns FALSE when any empty string present", {
  expect_false(all_truthy(c("human", "", "rat")))
})

test_that("all_truthy() returns FALSE when any whitespace present", {
  expect_false(all_truthy(c("human", "  ", "rat")))
})

test_that("all_truthy() returns TRUE when all values valid", {
  expect_true(all_truthy(c("human", "mouse", "rat")))
})

test_that("all_truthy() treats string 'NA' as truthy", {
  # This is the key difference from na_check()
  expect_true(all_truthy(c("NA", "NANA", "na")))
})

test_that("all_truthy() returns FALSE when all values missing", {
  expect_false(all_truthy(c(NA, "", "  ")))
})

test_that("all_truthy() handles single valid value", {
  expect_true(all_truthy("human"))
})

test_that("all_truthy() handles single NA value", {
  expect_false(all_truthy(NA))
})

test_that("all_truthy() handles numeric vectors", {
  expect_true(all_truthy(c(1, 2, 3)))
  expect_false(all_truthy(c(1, NA, 3)))
  
})

test_that("all_truthy() handles factor vectors", {
  expect_true(all_truthy(factor(c("human", "mouse", "rat"))))
  expect_false(all_truthy(factor(c("human", NA, "rat"))))
})


# Test any_truthy() ----------------------------------------------------

test_that("any_truthy() returns FALSE for NULL", {
  expect_false(any_truthy(NULL))
})

test_that("any_truthy() returns FALSE for empty vector", {
  expect_false(any_truthy(character(0)))
})

test_that("any_truthy() returns TRUE when any valid value present", {
  expect_true(any_truthy(c(NA, "", "  ", "rat")))
})

test_that("any_truthy() returns FALSE when all values missing", {
  expect_false(any_truthy(c(NA, "", "  ")))
})

test_that("any_truthy() returns TRUE when all values valid", {
  expect_true(any_truthy(c("human", "mouse", "rat")))
})

test_that("any_truthy() treats string 'NA' as truthy", {
  expect_true(any_truthy(c(NA, "", "NA")))
})

test_that("any_truthy() handles single valid value", {
  expect_true(any_truthy("human"))
})

test_that("any_truthy() handles single NA value", {
  expect_false(any_truthy(NA))
})

test_that("any_truthy() handles numeric vectors", {
  expect_true(any_truthy(c(1, 2, 3)))
  expect_true(any_truthy(c(NA, NA, 3)))
  expect_false(any_truthy(c(NA, NA, NA)))
})


# Test edge cases ------------------------------------------------------

test_that("functions handle vectors with special characters", {
  x <- c("test@example.com", "user-name", "value_123")
  expect_false(any_missingness(x))
  expect_true(all_truthy(x))
})

test_that("functions handle very long strings", {
  long_string <- paste(rep("a", 1000), collapse = "")
  x <- c(long_string, "short")
  expect_false(any_missingness(x))
  expect_true(all_truthy(x))
})

test_that("functions handle mixed whitespace types", {
  x <- c(" ", "\t", "\n", "\r", "  \t\n  ")
  expect_true(all(detected_missingness(x)))
  expect_true(any_missingness(x))
  expect_false(all_truthy(x))
  expect_false(any_truthy(x))
})

test_that("functions handle strings with leading/trailing whitespace", {
  # Note: strings with ANY non-whitespace character should be truthy
  x <- c("  value  ", " another ", "test")
  expect_false(any_missingness(x))
  expect_true(all_truthy(x))
})


# Test use cases -------------------------------------------------------

test_that("functions work for user input validation", {
  # Project names where "NA" might be valid
  projects <- c("Project Alpha", "Project NA", NA, "  ", "")
  clean <- projects[!detected_missingness(projects)]
  expect_equal(clean, c("Project Alpha", "Project NA"))
})

test_that("functions work for species validation", {
  
  report_species <- c("human", "mouse", "rat")
  expect_true(all_truthy(report_species))
  
  bad_species <- c("human", NA, "rat")
  expect_false(all_truthy(bad_species))
  
})

test_that("functions work for filtering operations", {
  study_vcs <- c("vc1", "vc2", NA, "", "vc3")
  valid_vcs <- study_vcs[!detected_missingness(study_vcs)]
  expect_equal(valid_vcs, c("vc1", "vc2", "vc3"))
})


# Test return types ----------------------------------------------------

test_that("detected_missingness() returns logical vector of correct length", {
  x <- c("a", "b", "c")
  result <- detected_missingness(x)
  expect_type(result, "logical")
  expect_length(result, 3)
})

test_that("any_missingness() returns single logical value", {
  result <- any_missingness(c("a", NA, "c"))
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("all_truthy() returns single logical value", {
  result <- all_truthy(c("a", "b", "c"))
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("any_truthy() returns single logical value", {
  result <- any_truthy(c("a", "b", "c"))
  expect_type(result, "logical")
  expect_length(result, 1)
})


# Comparison with shiny::isTruthy() -------------------------------------

test_that("all_truthy() is vector-aware unlike shiny::isTruthy()", {
  # shiny::isTruthy only checks first element
  # all_truthy checks ALL elements
  
  # All valid - should be TRUE
  expect_true(all_truthy(c("a", "b", "c")))
  
  # First valid, second NA - all_truthy should be FALSE
  # (different from shiny::isTruthy which would be TRUE)
  expect_false(all_truthy(c("a", NA)))
  
  # This is the key difference - all_truthy is vector-aware
  # Solves the problem from: if(!all(shiny::isTruthy(species)))
})


# Test NaN handling ----------------------------------------------------

test_that("detected_missingness() detects NaN in numeric vectors", {
  x <- c(1, NA, 3, NaN)
  result <- detected_missingness(x)
  expect_equal(result, c(FALSE, TRUE, FALSE, TRUE))
})

test_that("detected_missingness() handles NaN coerced to character", {
  # When NaN is in character vector, it becomes "NaN" string (valid data)
  x <- c("1", NA, "3", NaN)
  result <- detected_missingness(x)
  expect_equal(result, c(FALSE, TRUE, FALSE, FALSE))
})

test_that("all_truthy() returns FALSE when NaN present in numeric", {
  expect_false(all_truthy(c(1, 2, NaN)))
  expect_false(all_truthy(c(NaN, NA, 3)))
})

test_that("any_missingness() returns TRUE when NaN present in numeric", {
  expect_true(any_missingness(c(1, 2, NaN)))
  expect_true(any_missingness(c(NaN)))
})
