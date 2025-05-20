#ifndef ACTIVEPOINT_H_
#define ACTIVEPOINT_H_

class ActivePoint {
public:
	ActivePoint();
	virtual ~ActivePoint();

	int node = -1;
	short edge = -1;
	short idx = -1;
	short remainder = -1;

	void set(int n , short e, short p);
	void reset();
};

#endif /* ACTIVEPOINT_H_ */
