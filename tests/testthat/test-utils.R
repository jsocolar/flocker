test_that("log_inv_logit handles scalar input", {
  expect_equal(log_inv_logit(0), log(0.5))
  expect_equal(log_inv_logit(10), log(1 / (1 + exp(-10))))
})

test_that("log_inv_logit handles vector input", {
  input <- c(-3, 0, 3)
  expected_output <- log(1 / (1 + exp(-input)))
  expect_equal(log_inv_logit(input), expected_output)
})

test_that("log_inv_logit handles extreme values", {
  expect_equal(log_inv_logit(100), log(1 / (1 + exp(-100))))
  expect_equal(log_inv_logit(-100), log(1 / (1 + exp(100))))
  expect_false(identical(log_inv_logit(100), log(1 / (1 + exp(-100)))))
  expect_true(is.finite(log_inv_logit(-10000)))
})

test_that("log_inv_logit throws error for non-numeric input", {
  expect_error(log_inv_logit("abc"), "x must be numeric")
  expect_error(log_inv_logit(TRUE), "x must be numeric")
})


test_that("log1m_inv_logit returns correct output", {
  r <- rnorm(100)
  expect_identical(log_inv_logit(r), log1m_inv_logit(-r))
})

test_that("array utils work properly", {
  testmat <- matrix(c(NA, 1:3), nrow = 2)
  testdf <- as.data.frame(testmat)
  testarray <- array(1:8, dim = c(2,2,2))
  
  # expand_matrix
  expect_error(expand_matrix(1:4))
  expect_identical(expand_matrix(testmat), c(NA, 1:3))
  expect_identical(expand_matrix(testdf), c(NA, 1:3))
  
  # expand_array_3D
  expect_error(expand_array_3D(testarray[,,1]))
  expect_equal(expand_array_3D(testarray), matrix(c(1,2,5,6,3,4,7,8), ncol = 2))
  
  # nslice
  expect_error(nslice(testmat))
  expect_equal(nslice(testarray), 2)
  
  # stack_matrix
  expect_error(stack_matrix(testarray, 2))
  expect_equal(stack_matrix(testmat, 2), matrix(c(NA, 1, NA, 1, 2, 3, 2, 3), ncol = 2))
  expect_equivalent(as.matrix(stack_matrix(testdf, 2)), matrix(c(NA, 1, NA, 1, 2, 3, 2, 3), ncol = 2))
})


test_that("bookkeeping works properly", {
  # flocker_col_names
  expect_true(all(grepl("^ff_", flocker_col_names(3, 4))))
  expect_equal(7, length(flocker_col_names(3, 4)) - length(flocker_col_names()))
  
  # flocker_reserved
  expect_true(all(grepl(flocker_reserved()[1], flocker_col_names())))
  expect_true(all(grepl(flocker_reserved()[2], paste0(".", c(".", "foo", 1:2)))))
  
  # flocker_model_types
  expect_true(all(grepl("^single|^augmented|^multi", flocker_model_types())))
  
  # flocker_data_input_types
  expect_true(
    all(
      grepl(
        paste0("^", 
               paste(flocker_data_input_types(), collapse = "|^")
               ),
        flocker_model_types()
        )
      )
    )
  for(i in seq_along(flocker_data_input_types())){
    expect_false(
      all(
        grepl(
          paste0("^", 
                 paste(flocker_data_input_types()[-i], collapse = "|^")
          ),
          flocker_model_types()
        )
      )
    )
  }
  
  # flocker_data_output_types
  expect_true(all(grepl("^single|^augmented|^multi", flocker_data_output_types())))
})

test_that("fdtl function returns expected dataframe", {
  # Call the fdtl function
  result <- fdtl()
  
  # Check if the result is a dataframe
  expect_is(result, "data.frame")
  
  # Check if the result has the correct column names
  expect_named(result, c("model_type", "data_output_type", "data_input_type"))
  
  # Check if the result has the correct number of rows (assuming 10 model types)
  expect_equal(nrow(result), 10)
  
  # Check if the result has the correct number of columns
  expect_equal(ncol(result), 3)
  
  # Check if the model_type column contains the expected values
  expect_identical(result$model_type, flocker_model_types())
  
  # Check if the data_output_type and data_input_type columns contain the expected values
  expected_data_input_types <- c(
    "single", "single", "augmented", "multi", "multi", "multi", "multi",
    "single", "multi", "multi"
  )
  
  expected_data_output_types <- c(
    "single", "single_C", "augmented", "multi", "multi", "multi", "multi",
    "single", "multi", "multi"
  )
  
  
  expect_identical(result$data_output_type, expected_data_output_types)
  expect_identical(result$data_input_type, expected_data_input_types)
})

test_that("is_flocker_fit works", {
  expect_true(is_flocker_fit(example_flocker_model_single))
  expect_false(is_flocker_fit("foo"))
  expect_false(is_flocker_fit(NULL))
  expect_false(is_flocker_fit(list(f = example_flocker_model_single)))
})

test_that("type_flocker_fit function returns expected string", {
  expect_identical(type_flocker_fit(example_flocker_model_single), "single")
  
  # Create a dummy flocker_fit object
  dummy_flocker_fit <- structure(
    list(),
    class = "flocker_fit",
    data_type = "multi",
    multiseason = "colex",
    multi_init = "equilibrium"
  )
  
  # Call the type_flocker_fit function
  result <- type_flocker_fit(dummy_flocker_fit)
  
  # Check if the result is a character string
  expect_is(result, "character")
  
  # Check if the result has the correct value
  expected_value <- "multi_colex_eq"
  expect_identical(result, expected_value)
  
  # Check if the function throws an error for non-flocker_fit objects
  non_flocker_fit <- list()
  expect_error(type_flocker_fit(non_flocker_fit), "x must be a flocker_fit object")
  
  # Check if the function throws an error for objects with missing or altered attributes
  corrupted_flocker_fit <- structure(
    list(),
    class = "flocker_fit",
    data_type = "MT1"
  )
  expect_error(type_flocker_fit(corrupted_flocker_fit), "the attributes of the flocker_fit object have been altered or corrupted")
})


test_that("emission_likelihood function returns expected output", {
  # Test cases for state 0
  obs1 <- matrix(c(0, 0, 0, 0, NA), nrow = 1)
  det1 <- matrix(c(0.3, 0.4, 0.5, 0.6, 0.7), nrow = 1)
  expected_output1 <- 1
  
  obs2 <- rbind(obs1, c(0, 1, 0, 1, 1))
  det2 <- rbind(det1, c(0.3, 0.4, 0.5, 0.6, 0.7))
  expected_output2 <- c(1,0)
  
  # Test cases for state 1
  obs3 <- matrix(c(1, 1, 0, 1, NA), nrow = 1)
  det3 <- matrix(c(0.3, 0.4, 0.5, 0.6, 0.7), nrow = 1)
  expected_output3 <- 0.3 * 0.4 * (1 - 0.5) * 0.6
  
  obs4 <- matrix(c(0, 0, 1, 1, 0), nrow = 1)
  det4 <- matrix(c(0.3, 0.4, 0.5, 0.6, 0.7), nrow = 1)
  expected_output4 <- (1 - 0.3) * (1 - 0.4) * 0.5 * 0.6 * (1 - 0.7)
  
  # Test the function with the test cases
  result1 <- emission_likelihood(0, obs1, det1)
  expect_identical(result1, expected_output1)
  
  result2 <- emission_likelihood(0, obs2, det2)
  expect_identical(result2, expected_output2)
  
  result3 <- emission_likelihood(1, obs3, det3)
  expect_identical(result3, expected_output3)
  
  result4 <- emission_likelihood(1, obs4, det4)
  expect_equal(result4, expected_output4)
  
  # Test the function with invalid inputs
  obs_invalid1 <- matrix(c(0, 1, -1, 1, 0), nrow = 1)
  obs_invalid2 <- matrix(c(0, 1, 2, 1, 0), nrow = 1)
  det_invalid1 <- matrix(c(0.5, 0.5, -0.5, 0.5, 0.5), nrow = 1)
  det_invalid2 <- matrix(c(0.5, 0.5, 1.5, 0.5, 0.5), nrow = 1)
  det_invalid3 <- matrix(c(NA, .5, .5, .5, .5), nrow = 1)
  
  expect_error(emission_likelihood(0, obs_invalid1, det1), "all\\(obs")
  expect_error(emission_likelihood(0, obs_invalid2, det1), "all\\(obs")
  expect_error(emission_likelihood(0, obs1, det_invalid1), "all\\(det")
  expect_error(emission_likelihood(0, obs1, det_invalid2), "all\\(det")
  expect_error(emission_likelihood(0, obs1, det_invalid3))
  expect_error(emission_likelihood(1, obs1, det_invalid3))
})

test_that("formula_error works", {
  result <- formula_error("x")
  expect_is(result, "character")
  expect_identical(result, "Formula error: x formula has incorrect syntax.")
})

test_that("is_formula and is_flocker_formula work", {
  expect_true(is_formula(stats::formula(a ~ b)))
  expect_true(is_formula(stats::formula(~ b)))
  expect_true(is_formula(stats::formula(~1)))
  expect_true(is_formula(stats::formula(~ a ~ b)))
  
  expect_false(is_formula(1))
  expect_false(is_formula(list(stats::formula(~1), stats::formula(a ~ b))))
  
  expect_false(is_flocker_formula(stats::formula(a ~ b)))
  expect_false(is_flocker_formula(stats::formula(~ a ~ b)))
  expect_true(is_flocker_formula(stats::formula(~ b)))
  expect_true(is_flocker_formula(stats::formula(~ b + (1 || c))))
})

test_that("is_flocker_data works", {
  expect_false(is_flocker_data(example_flocker_data))
  
  expect_true(
    make_flocker_data(example_flocker_data$obs, example_flocker_data$unit_covs) |>
      is_flocker_data()
  )
})  



test_that("is_named_list returns TRUE for a named list with unique names", {
  named_list <- list(a = 1, b = 2, c = 3)
  expect_true(is_named_list(named_list))
})

test_that("is_named_list returns FALSE for an unnamed list", {
  unnamed_list <- list(1, 2, 3)
  expect_false(is_named_list(unnamed_list))
})

test_that("is_named_list returns FALSE for a partially named list", {
  partially_named_list <- list(a = 1, 2, c = 3)
  expect_false(is_named_list(partially_named_list))
})

test_that("is_named_list returns FALSE for a named list with duplicate names", {
  duplicate_named_list <- list(a = 1, b = 2, a = 3)
  expect_false(is_named_list(duplicate_named_list))
})

test_that("is_named_list returns FALSE for an empty list", {
  empty_list <- list()
  expect_false(is_named_list(empty_list))
})

test_that("is_named_list returns FALSE for a non-list object", {
  non_list_object <- 42
  expect_false(is_named_list(non_list_object))
})

test_that("is_one_logical returns TRUE for a single logical value (TRUE)", {
  single_logical_true <- TRUE
  expect_true(is_one_logical(single_logical_true))
})

test_that("is_one_logical returns TRUE for a single logical value (FALSE)", {
  single_logical_false <- FALSE
  expect_true(is_one_logical(single_logical_false))
})

test_that("is_one_logical returns FALSE for a numeric value", {
  numeric_value <- 42
  expect_false(is_one_logical(numeric_value))
})

test_that("is_one_logical returns FALSE for a character value", {
  character_value <- "test"
  expect_false(is_one_logical(character_value))
})

test_that("is_one_logical returns FALSE for a vector of logical values", {
  logical_vector <- c(TRUE, FALSE)
  expect_false(is_one_logical(logical_vector))
})

test_that("is_one_logical returns FALSE for a NULL value", {
  null_value <- NULL
  expect_false(is_one_logical(null_value))
})

test_that("is_one_pos_int returns TRUE for a single positive integer greater than m", {
  single_positive_int <- 5
  m <- 3
  expect_true(is_one_pos_int(single_positive_int, m))
})

test_that("is_one_pos_int returns FALSE for a single positive integer equal to m", {
  single_positive_int <- 5
  m <- 5
  expect_false(is_one_pos_int(single_positive_int, m))
})

test_that("is_one_pos_int returns FALSE for a single positive integer less than m", {
  single_positive_int <- 3
  m <- 5
  expect_false(is_one_pos_int(single_positive_int, m))
})

test_that("is_one_pos_int returns FALSE for a single negative integer", {
  negative_int <- -5
  m <- 0
  expect_false(is_one_pos_int(negative_int, m))
})

test_that("is_one_pos_int returns FALSE for a single non-integer numeric value", {
  non_integer_numeric <- 2.5
  m <- 0
  expect_false(is_one_pos_int(non_integer_numeric, m))
})

test_that("is_one_pos_int returns FALSE for a character value", {
  character_value <- "test"
  m <- 0
  expect_false(is_one_pos_int(character_value, m))
})

test_that("is_one_pos_int returns FALSE for a vector of integers", {
  integer_vector <- c(2, 5, 8)
  m <- 0
  expect_false(is_one_pos_int(integer_vector, m))
})

test_that("is_one_pos_int returns FALSE for a NULL value", {
  null_value <- NULL
  m <- 0
  expect_false(is_one_pos_int(null_value, m))
})

test_that("shared_elements returns correct shared elements for two non-empty vectors", {
  vec1 <- c(1, 2, 3, 4, 5, 6)
  vec2 <- c(4, 5, 6, 7, 8, 9)
  expected_output <- c(4, 5, 6)
  expect_equal(shared_elements(vec1, vec2), expected_output)
})

test_that("shared_elements returns an empty vector when there are no shared elements", {
  vec1 <- c(1, 2, 3)
  vec2 <- c(4, 5, 6)
  expected_output <- integer(0)
  expect_equal(shared_elements(vec1, vec2), expected_output)
})

test_that("shared_elements handles vectors with duplicate elements correctly", {
  vec1 <- c(1, 1, 2, 2, 3, 3)
  vec2 <- c(2, 2, 3, 3, 4, 4)
  expected_output <- c(2, 3)
  expect_equal(shared_elements(vec1, vec2), expected_output)
})

test_that("shared_elements returns an empty vector when both input vectors are empty", {
  vec1 <- integer(0)
  vec2 <- integer(0)
  expected_output <- integer(0)
  expect_equal(shared_elements(vec1, vec2), expected_output)
})

test_that("shared_elements returns an empty vector when one input vector is empty", {
  vec1 <- c(1, 2, 3)
  vec2 <- integer(0)
  expected_output <- integer(0)
  expect_equal(shared_elements(vec1, vec2), expected_output)
})

test_that("shared_elements works with character vectors", {
  vec1 <- c("apple", "banana", "cherry")
  vec2 <- c("banana", "cherry", "date")
  expected_output <- c("banana", "cherry")
  expect_equal(shared_elements(vec1, vec2), expected_output)
})

test_that("shared_elements works with mixed-type vectors", {
  vec1 <- c(1, "apple", 2, "banana")
  vec2 <- c("banana", 2, "date", 3)
  expected_output <- c(2, "banana")
  expect_equal(shared_elements(vec1, vec2), expected_output)
})

test_that("max_position_not_na works correctly", {
  # Test with no NAs and no -99 values
  expect_equal(max_position_not_na(c(1, 2, 3, 4, 5)), 5)
  
  # Test with NAs only
  expect_equal(max_position_not_na(c(NA, NA, NA)), 0)
  
  # Test with -99 values treated as NAs
  expect_equal(max_position_not_na(c(1, 2, -99, 4, 5), treat_m99_NA = TRUE), 5)
  expect_equal(max_position_not_na(c(1, 2, -99, 4, -99), treat_m99_NA = TRUE), 4)
  
  # Test with NAs and -99 values at the end of the vector
  expect_equal(max_position_not_na(c(1, 2, 3, NA, -99), treat_m99_NA = TRUE), 3)
  expect_equal(max_position_not_na(c(1, 2, 3, NA, -99), treat_m99_NA = FALSE), 5)
  
  # Test with NAs and -99 values at the beginning of the vector
  expect_equal(max_position_not_na(c(NA, -99, 1, 2, 3), treat_m99_NA = TRUE), 5)
  expect_equal(max_position_not_na(c(NA, -99, 1, 2, 3), treat_m99_NA = FALSE), 5)
  
  # Test with non-numeric vector
  expect_equal(max_position_not_na(c("A", "B", "C", "D", "E")), 5)
  expect_equal(max_position_not_na(c("A", "B", "C", "D", "NA")), 5)
  expect_equal(max_position_not_na(c("A", "B", "C", "D", "-99"), treat_m99_NA = TRUE), 4)
})


test_that("remove_rownames works correctly", {
  # Test with a matrix
  m1 <- matrix(1:9, nrow = 3, ncol = 3, dimnames = list(c("r1", "r2", "r3"), c("c1", "c2", "c3")))
  m1_expected <- matrix(1:9, nrow = 3, ncol = 3, dimnames = list(c(), c("c1", "c2", "c3")))
  expect_equal(remove_rownames(m1), m1_expected)
  
  # Test with a data.frame
  df1 <- data.frame(a = 1:3, b = 4:6, row.names = c("r1", "r2", "r3"))
  df1_expected <- data.frame(a = 1:3, b = 4:6)
  expect_equal(remove_rownames(df1), df1_expected)
  
  # Test with a tibble
  if (requireNamespace("tibble", quietly = TRUE)) {
    tb1 <- tibble::tibble(a = 1:3, b = 4:6)
    suppressWarnings(rownames(tb1) <- c("r1", "r2", "r3"))
    tb1_expected <- tibble::tibble(a = 1:3, b = 4:6)
    expect_equal(remove_rownames(tb1), tb1_expected)
  }
})
