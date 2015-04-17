#include "cutest.h"

void test_checks() {
        TEST_CHECK(5 == 5);
        TEST_CHECK(2 != 5);
}

void test_verbose_checks() {
        TEST_CHECK_(3 != 5, "%d != %d is very true!", 3, 5);
}

TEST_LIST = {
        {"test_checks", test_checks},
        {"test_verbose_checks", test_verbose_checks},
        { 0 }
};