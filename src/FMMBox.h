#ifndef Box_h
#define Box_h

#include "GeneralUtilities.h"
#include "Point.h"

class Box {

public:

	/// Level of the box
	int level;
	/// Index of the box
	int index;
	
	/// Size of the box
	double size;
	/// Coordinate of the box center
	Complex center;
	
	/// Truncation number
	int degree;	
	
	/// External Multipole expansion coefficients
	ComplexVec externalMultipoleCoeffs;
	/// Local Multipole expansion coefficients
	ComplexVec localMultipoleCoeffs;
	/// Local Multipole expansion coefficients (temporary)
	ComplexVec localMultipoleCoeffsTilde;

	/// Collection of sources inside this box
	std::vector<Point*> sources;
	/// Collection of targets inside this box
	std::vector<Point*> targets;


	/// Constructor
	Box(const int level, const int index, const int degree) 
	: level(level), index(index), degree(degree) 
	{
		size = pow(2.0, -level);
		center = (uninterleave(index, level) + Complex(0.5, 0.5)) * size;
		externalMultipoleCoeffs.resize(degree, Complex(0,0));
		localMultipoleCoeffs.resize(degree, Complex(0,0));
		localMultipoleCoeffsTilde.resize(degree, Complex(0,0));
	}

	// Destructor
	~Box() 
	{ 

	}

	/// Overload comparison operator
	inline bool operator==(const Box& other) const 
	{
		return (this->level == other.level) && (this->index == other.index);
	}

	inline bool operator!=(const Box& other) const 
	{
		return !(*this == other);
	}

	/// Add a source point to this box
	inline void AddSource(Point* source) 
	{
		sources.push_back(source);
	}

	/// Add a target point to this box
	inline void AddTarget(Point* target)
	{
		targets.push_back(target);
	}

};

#endif