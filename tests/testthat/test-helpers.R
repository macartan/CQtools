context("helpers")

testthat::test_that(
	desc = "zero_range",
	code = {
		expect_true(!zero_range(1:2))
		expect_true(zero_range(rep(2,2)))
		expect_true(zero_range(3))
	}
)

testthat::test_that(
	desc = "combine_lists",
	code = {
		a <- combine_lists(list(B = 2, A = 1), list(A = 3, c = 4))

		expect_true(length(a) == 3)
	}
)
