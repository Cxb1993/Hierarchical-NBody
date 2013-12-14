#include "MLFMM.h"

MLFMM::MLFMM(const int levels, Potential& potential) 
: levels(levels), flops(0) 
{
	maxLevel = levels - 1;
	this->potential = &potential;
	InitializeStructure();
}

MLFMM::~MLFMM() 
{
	for (int level = maxLevel; level >= 0; level--)
		for (int index = 0; index < structure[level].size(); index++)
			delete structure[level][index]; 
}

void MLFMM::Solve() 
{
	MultipoleExpansion();
	MultipoleToMultipoleTranslation();
	MultipoleToLocalTranslation();
	LocalToLocalTranslation();
	LocalExpansion();
}

void MLFMM::DirectSolve() 
{
	for (auto &target : targets) {
		target->potential = 0.0;
		for (auto &source : sources){
			if (source->coord != target->coord) {
				target->potential += potential->DirectEvaluate(target->coord, source->coord);
			}
		}
	}
}

void MLFMM::InitializeStructure() 
{
	structure.resize(levels, std::vector<Box*>());
	for (int level = 0; level < levels; level++) {
		structure[level].resize((int)pow(4, level));
		for (int index = 0; index < structure[level].size(); index++) {
			structure[level][index] = new Box(level, index, potential->degree);
		}
	}
}

void MLFMM::AddSource(Point* source) {
	sources.push_back(source);
	int index = GetBoxIndex(source->coord, maxLevel);
	structure[maxLevel][index]->AddSource(source);
}

void MLFMM::AddTarget(Point* target) {
	targets.push_back(target);
	int index = GetBoxIndex(target->coord, maxLevel);
	structure[maxLevel][index]->AddTarget(target);
}

int MLFMM::GetBoxIndex(const Complex& coord, const int level) 
{
	return interleave(
		(int)floor(real(coord) * pow(2, level)), 
		(int)floor(imag(coord) * pow(2, level)), level);
}

Box* MLFMM::GetParent(Box* box) 
{
	return structure[box->level - 1][box->index >> 2];
}

std::vector<Box*> MLFMM::GetChildren(Box* box) 
{
	std::vector<Box*> children;
	for (int i = 0; i < 4; i++) {
		children.push_back(structure[box->level + 1][(box->index << 2) + i]);
	}
	return children;
}

std::vector<Box*> MLFMM::GetNeighbors(Box* box) 
{
	std::vector<Box*> neighbors;
	Complex location = uninterleave(box->index, box->level);
	int x = (int)real(location);
	int y = (int)imag(location);
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			// keep in bounds
			if ((i != 0 || j != 0)  
				&& (x + i >= 0)
				&& (y + j >= 0)
				&& (x + i < pow(2, box->level)) 
				&& (y + j < pow(2, box->level))) {
				neighbors.push_back(structure[box->level][interleave(x + i, y + j, box->level)]);
			}
		}
	}
	return neighbors;
}

bool MLFMM::IsNeighbor(Box* box, Box* candidate) 
{
	for (auto &neighbor : GetNeighbors(box)) 
		if (neighbor == candidate) 
			return true;
	return false;
}

std::vector<Box*> MLFMM::GetInteractionList(Box* box) 
{
	std::vector<Box*> interactionList;
	std::vector<Box*> neighbors = GetNeighbors(box);
	for (auto &parentsNeighbor : GetNeighbors(GetParent(box)))
		for (auto &candidate : GetChildren(parentsNeighbor))
			if (box != candidate && !IsNeighbor(box, candidate))
				interactionList.push_back(candidate);
	return interactionList;
}

void MLFMM::MultipoleExpansion() 
{
	for (auto &box : structure[maxLevel]){
		for (auto &source : box->sources){
			box->externalMultipoleCoeffs += potential->GetMultipoleCoeffs(source->coord, box->center);
			flops += potential->degree; 
		}
	}
}

void MLFMM::MultipoleToMultipoleTranslation() 
{
	for (int level = maxLevel; level >= 2; level--) {
		for (auto &box : structure[level]) {
			Box* parent = GetParent(box);
			parent->externalMultipoleCoeffs += potential->MultipoleToMultipole(box->center, parent->center, box->externalMultipoleCoeffs);
			flops += potential->degree * potential->degree; 
		}
	}
}

void MLFMM::MultipoleToLocalTranslation() 
{
	for (int level = 2; level <= maxLevel; level++) {
		for (auto &box : structure[level]) {
			for (auto &neighbor : GetInteractionList(box)) {
				box->localMultipoleCoeffsTilde += potential->MultipoleToLocal(neighbor->center, box->center, neighbor->externalMultipoleCoeffs);
				flops += potential->degree * potential->degree;
			}
		}
	}
}

void MLFMM::LocalToLocalTranslation() 
{
	for (auto &box : structure[2]) {
		box->localMultipoleCoeffs += box->localMultipoleCoeffsTilde;
		flops += potential->degree;
	}
	for (int level = 2; level < maxLevel; level++) {
		for (auto &box : structure[level]) {
			for (auto &child : GetChildren(box)) {
				child->localMultipoleCoeffs += child->localMultipoleCoeffsTilde;
				child->localMultipoleCoeffs += potential->LocalToLocal(box->center, child->center, box->localMultipoleCoeffs);
				flops += potential->degree * potential->degree + potential->degree;
			}
		}
	}
}

void MLFMM::MLFMM::LocalExpansion() 
{
	for (auto &box : structure[maxLevel]) {
		for (auto &target : box->targets) {
			double outsidePotential = 0.0;
			ComplexVec localVector = potential->GetLocalCoeffs(target->coord, box->center);
			for (int k = 0; k < potential->degree; k++) {
				outsidePotential += real(box->localMultipoleCoeffs[k] * localVector[k]);
			}
			flops += potential->degree;
			double insidePotential = 0.0;
			for (auto &source : box->sources) {
				if (source->coord != target->coord){
					insidePotential += potential->DirectEvaluate(target->coord, source->coord);
					flops += 1;
				}
			}
			for (auto &neighbor : GetNeighbors(box)) {
				for (auto &source : neighbor->sources) {
					if (source->coord != target->coord) {
						insidePotential += potential->DirectEvaluate(target->coord, source->coord);
						flops += 1;
					}
				}
			}
			target->potential = outsidePotential + insidePotential;
			flops += 1;
		}
	}
}
