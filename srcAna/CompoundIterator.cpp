// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//
// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.

#include "CompoundIterator.h"

#include <cmath>

CompoundIterator::CompoundIterator() {
  first = 0;
  second = 0;
  forwardMap();
}

CompoundIterator::CompoundIterator(int compound_) {
  compound = compound_;
  backwardMap();
}

CompoundIterator::CompoundIterator(int first_, int second_) {
  first = first_;
  second = second_;
  forwardMap();
}

void CompoundIterator::forwardMap() {
  compound = first * (first + 1) - second - 1;
}

void CompoundIterator::backwardMap() {
  first = (int)std::sqrt(compound + 1.0);
  second = -(compound + 1) + first * (first + 1);
}

CompoundIterator::~CompoundIterator() {

}

void CompoundIterator::init(int compound_) {
  compound = compound_;
  backwardMap();
}

void CompoundIterator::init(int first_, int second_) {
  first = first_;
  second = second_;
  forwardMap();
}

int CompoundIterator::max(int first_) { return first_ * first_ + 2 * first_; }

void CompoundIterator::operator++() {
  compound++;
  backwardMap();
}

void CompoundIterator::operator--() {
  if (compound > 0) {
    compound--;
    backwardMap();
  }
}

CompoundIterator CompoundIterator::operator+(int increase_) {
  return CompoundIterator(compound + increase_);
}

CompoundIterator CompoundIterator::operator-(int decrease_) {
  if (compound >= decrease_) {
    return CompoundIterator(compound - decrease_);
  }

  return CompoundIterator(0);
}

void CompoundIterator::operator=(int compound_) {
  compound = compound_;
  backwardMap();
}

bool CompoundIterator::operator<(int compound_) { return compound < compound_; }

bool CompoundIterator::operator>(int compound_) { return compound > compound_; }

bool CompoundIterator::operator<=(int compound_) {
  return compound <= compound_;
}

bool CompoundIterator::operator>=(int compound_) {
  return compound >= compound_;
}

void CompoundIterator::operator++(int) {
  compound++;
  backwardMap();
}

void CompoundIterator::operator--(int) {
  if (compound >= 0) {
    compound--;
    backwardMap();
  }
}

bool CompoundIterator::operator==(int compound_) {
  return compound == compound_;
}

CompoundIterator::operator optimet::t_uint() const { return compound; }
