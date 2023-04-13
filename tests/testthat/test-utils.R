test_that("math utils work properly", {
  expect_error(log_inv_logit("foo"))
  expect_equal(log_inv_logit(0), log(.5))
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
  
  
  # TODO: test get_positions_single (requires creating an example flocker fit in data)
  
})
