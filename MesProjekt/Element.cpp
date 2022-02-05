#include "Element.h"
//Element::Element(): ID[0](0), ID[1](0), ID[2](0), ID[3](0)
//{
//
//}
//Element::Element(int id1, int id2, int id3, int id4)
//{
//	this->ID[0] = id1;
//	this->ID[1] = id2;
//	this->ID[2] = id3;
//	this->ID[3] = id4;
//}
void Element::print()
{
	std::cout << "(" << ID[0] << " " << ID[1] << " " << ID[2] << " " << ID[3] << ")\n";
}