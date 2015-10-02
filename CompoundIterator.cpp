#include "CompoundIterator.h"

using std::sqrt;

CompoundIterator::CompoundIterator()
{
	first = 0;
	second = 0;
	forwardMap();
}

CompoundIterator::CompoundIterator(int compound_)
{
	compound = compound_;
	backwardMap();
}

CompoundIterator::CompoundIterator(int first_, int second_)
{
	first = first_;
	second = second_;
	forwardMap();
}

void CompoundIterator::forwardMap()
{
	compound = first*(first+1)-second-1;
}

void CompoundIterator::backwardMap()
{
	first = (int) sqrt(compound + 1.0);
	second = -(compound+1)+first*(first+1);
}

CompoundIterator::~CompoundIterator()
{
	//
}

void CompoundIterator::init(int compound_)
{
	compound = compound_;
	backwardMap();
}

void CompoundIterator::init(int first_, int second_)
{
	first = first_;
	second = second_;
	forwardMap();
}

int CompoundIterator::max(int first_)
{
	return first_*first_ + 2*first_;
}

void CompoundIterator::operator ++()
{
	compound++;
	backwardMap();
}

void CompoundIterator::operator --()
{
	if(compound > 0)
	{
		compound--;
		backwardMap();
	}
}

CompoundIterator CompoundIterator::operator +(int increase_)
{
	return CompoundIterator(compound + increase_);
}

CompoundIterator CompoundIterator::operator -(int decrease_)
{
	if(compound >= decrease_)
	{
		return CompoundIterator(compound - decrease_);
	}

	return CompoundIterator(0);
}

void CompoundIterator::operator =(int compound_)
{
	compound = compound_;
	backwardMap();
}

bool CompoundIterator::operator <(int compound_)
{
	return compound < compound_;
}

bool CompoundIterator::operator >(int compound_)
{
	return compound > compound_;
}

bool CompoundIterator::operator <=(int compound_)
{
	return compound <= compound_;
}

bool CompoundIterator::operator >=(int compound_)
{
	return compound >= compound_;
}

void CompoundIterator::operator ++(int compound_)
{
	compound ++;
	backwardMap();
}

void CompoundIterator::operator --(int compound_)
{
	if(compound >= 0)
	{
		compound --;
		backwardMap();
	}
}

bool CompoundIterator::operator ==(int compound_)
{
	return compound == compound_;
}

CompoundIterator::operator int()
{
	return compound;
}

CompoundIterator::operator long()
{
	return compound;
}



