#ifndef MLFMM_h
#define MLFMM_h

#include "GeneralUtilities.h"
#include "Point.h"
#include "FMMPotential.h"
#include "FMMBox.h"

class MLFMM {

public:

	/// Number of levels in the tree
	int levels;
	/// Index of the deepest level in the tree
	int maxLevel; 
	
	/// Collection of sources
	std::vector<Point*> sources;
	/// Collection of targets
	std::vector<Point*> targets;
	
	/// Hierarchical tree structure
	std::vector<std::vector<Box*>> structure; 
	
	/// Multipole potential
	Potential* potential; 
	
	/// FLOP counter
	long flops;

	/// Constructor 
	MLFMM(const int levels, Potential& potential);

	/// Destructor 
	~MLFMM();

	/// Initialize FMM tree structure
	void InitializeStructure();

	/// Add a source to the FMM tree
	void AddSource(Point* source);

	/// Add a target to the FMM tree
	void AddTarget(Point* target);

	/// Solve by direct evaluation of the potential
	void DirectSolve();

	/// Solve using the Fast Multipole Method
	void Solve();

	/// Multipole expansion
	void MultipoleExpansion();

	/// Multipole-to-multipole translation
	void MultipoleToMultipoleTranslation();

	/// Multipole-to-local translation
	void MultipoleToLocalTranslation();

	/// Local-to-local translation
	void LocalToLocalTranslation();

	/// Local expansion
	void LocalExpansion();

	/// Get the index of box from a coordinate and level
	int GetBoxIndex(const Complex& coord, const int level);

	/// Get the parent of a box
	Box* GetParent(Box* box);

	/// Get the children of a box as a std::vector of boxes
	std::vector<Box*> GetChildren(Box* box);
	
	/// Get the neighbors of a box as a std::vector of boxes
	std::vector<Box*> GetNeighbors(Box* box);

	/// Check if a candidate is a neighbor of box
	bool IsNeighbor(Box* box, Box* candidate);

	/// Get the interaction list of a box as a std::vector of boxes
	std::vector<Box*> GetInteractionList(Box* box);

};

#endif