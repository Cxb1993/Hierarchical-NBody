#ifndef Point_h
#define Point_h

class Point {
	
public:

	/// Coordinate: x component is real part, y component is imaginary part
	Complex coord;
	
	/// Global index of the point, for book-keeping
	int index;

	/// Potential evaluated at this point
	double potential;

	/// Constructor
	Point(const Complex& coord, const int index) 
	: coord(coord), index(index)
	{
		potential = 0;
	}

};

#endif