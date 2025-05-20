#include "ActivePoint.h"

ActivePoint::ActivePoint() {
}

ActivePoint::~ActivePoint() {
}

void ActivePoint::set(int n, short e, short p) {
	node = n;
	idx = p;
	edge = e;
}

void ActivePoint::reset() {
	node = 0;
	idx = 0;
	edge = 0;
	remainder = 0;
}
