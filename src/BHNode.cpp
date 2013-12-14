#include <cstdio>
#include <cassert>
#include "BHNode.h"

BHNode::BHNode(const Complex& center, const Complex& size, const int depth, int maxDepth)
: center(center), size(size), hasChildren(false), depth(depth), maxDepth(maxDepth)
{
	if (depth == 0) 
		BHNode::flops = 0;
	for (int quadrant = 0; quadrant < 4; quadrant++) 
		hasChild[quadrant] = false; 
}

BHNode::~BHNode() 
{
	for (int quadrant = 0; quadrant < 4; quadrant++)
		if (hasChild[quadrant])
			delete children[quadrant];
}

int BHNode::GetQuadrant(const Complex& coord) 
{
  if (real(coord) >= real(center) && imag(coord) >= imag(center)) return 0;
  if (real(coord) <  real(center) && imag(coord) >= imag(center)) return 1;
  if (real(coord) <  real(center) && imag(coord) <  imag(center)) return 2;
  if (real(coord) >= real(center) && imag(coord) <  imag(center)) return 3;
  return -1;
}

void BHNode::CreateChild(const int quadrant) 
{
	Complex size(real(this->size) / 2, imag(this->size) / 2);
	int depth = this->depth + 1;
	Complex center = this->center;
	switch (quadrant) {
		case 0: center += size; break;
		case 1: center -= conj(size); break;
		case 2: center -= size; break;
		case 3: center += conj(size); break;
	}
	children[quadrant] = new BHNode(center, size, depth, maxDepth);
	hasChild[quadrant] = true;
	hasChildren = true;
}

int BHNode::AddSource(Point* source) 
{
	// we are at maximum depth, insert here
	if (depth >= maxDepth) {
		sources.push_back(source);
		return depth;
	}
	// push an existing particle down to child level, 
	// then insert new particle at child level
	if (sources.size() == 1) {
		Point* existing = sources.back();
		sources.pop_back();
		assert(sources.size() == 0);
		int quadrant = GetQuadrant(existing);
		if (!hasChild[quadrant])
			CreateChild(quadrant);
		children[quadrant]->AddSource(existing);
		quadrant = GetQuadrant(source);
		if (!hasChild[quadrant])
			CreateChild(quadrant);
		return children[quadrant]->AddSource(source);
	}
	// otherwise insert normally
	if (sources.size() == 0) {
		if (hasChildren) {
			int quadrant = GetQuadrant(source);
			if (!hasChild[quadrant]) 
				CreateChild(quadrant);
			return children[quadrant]->AddSource(source);
		} else {
			sources.push_back(source);
			assert(sources.size() == 1);
			return depth;
		}
	}
	return -1;
}

void BHNode::ComputeChargeDistribution() 
{
	charge = 0;
	centerOfCharge = Complex(0, 0);
	if (sources.size() > 0) {
		for (auto &source : sources) {
			charge += 1.0;
			centerOfCharge += Complex(source->coord);
			BHNode::flops++;
		}
	} else {
		for (int quadrant = 0; quadrant < 4; quadrant++) {
			if (hasChild[quadrant]) {
				children[quadrant]->ComputeChargeDistribution();
				charge += children[quadrant]->charge;
				centerOfCharge += (children[quadrant]->centerOfCharge * children[quadrant]->charge);
				BHNode::flops++;
			}
		}
	}
	centerOfCharge /= charge;
}

double BHNode::ComputePotential(const Point* target, const double theta)
{
	if (sources.size() > 0)
		return ComputePotentialDirect(sources, target);

	double distance = norm(target->coord - centerOfCharge);
	double ratio = sqrt(distance / norm(size));

	if (ratio > theta) {
		BHNode::flops++;
		return charge * real(log(target->coord - centerOfCharge));
	} else {
		double potential = 0;
		for (int quadrant = 0; quadrant < 4; quadrant++)
			if (hasChild[quadrant])
				potential += children[quadrant]->ComputePotential(target, theta);
		return potential;
	}
}

double BHNode::ComputePotentialDirect(const std::vector<Point*>& sources, const Point* target) {
	double potential = 0;
	for (auto &source : sources) {
		if (source->coord != target->coord) {
			potential += real(log(target->coord - source->coord));
			BHNode::flops++;
		}
	}
	return potential;
}