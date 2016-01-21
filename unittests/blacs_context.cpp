#include <iostream>
#include "catch.hpp"

#include "Types.h"
#include "BlacsContext.h"

using namespace optimet;

TEST_CASE("Creates a blacs context 1x1") {
  BlacsContext context{1, 1};
}
