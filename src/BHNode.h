#ifndef BHNode_h
#define BHNode_h

#include "GeneralUtilities.h"
#include "Point.h"

/// Barnes-Hut Treecode, Adaptive Quadtree, 2D Coulomb Potential
class BHNode {

public:
	/// Center of this node
	Complex center;
	/// Size of this node
	Complex size;
	/// Collection of sources in this node
	std::vector<Point*> sources;
	/// Pointers to children
	BHNode *children[4];
	/// Flag for existing children
	bool hasChild[4];
	/// Flag for any existing children
	bool hasChildren;
	/// Depth in the tree
	int depth;
	/// Maximum depth of tree
	int maxDepth;
	/// Total charge in this node
	double charge;
	/// Center of total charge
	Complex centerOfCharge;
	/// FLOP counter
	static long flops;
	/// Constructor
	BHNode(const Complex& center, const Complex& size, const int depth, int maxDepth);
	/// Destructor
	~BHNode();
	/// Get index of quadrant for a coordinate
	int GetQuadrant(const Complex& coord);
	/// Get index of quadrant for a point
	inline int GetQuadrant(const Point* p) { return GetQuadrant(p->coord); }
	/// Create a child at a given quadrant
	void CreateChild(const int quadrant);
	/// Add a source to this node (recursive)
	int AddSource(Point* source);
	/// Compute the charge distribution of this node (recursive)
	void ComputeChargeDistribution();
	/// Compute the approximate potential due to the charge distribution of this node (recursive)
	double ComputePotential(const Point* target, const double theta);
	/// Compute the potential directoy due to a collection of sources
	double ComputePotentialDirect(const std::vector<Point*>& sources, const Point* target);
};

#endif