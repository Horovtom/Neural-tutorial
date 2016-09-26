// Neural-net-tutorial.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include <iostream>

class Net {
public: 
	Net(const std::vector<unsigned> &topology);
	void feedForward(const std::vector<double> &inputVals);
	void backProp(const std::vector<double> &targetVals);
	//getResults(resultVals) is a function that only reads the result, does not modify the Net object at all. Thats why its const
	void getResults(std::vector<double> &resultVals) const;

private:

};

int main()
{
	// e.g., {3, 2, 1} (3 layers, 1st - 3 nodes, 2nd - 2 nodes, 3rd - 1 node
	std::vector<unsigned> topology;
	//Construct
	Net myNet(topology);
	
	std::vector<double> inputVals;
	//Train
	myNet.feedForward(inputVals);

	std::vector<double> targetVals;
	myNet.backProp(targetVals);
	//Save results

	std::vector<double> resultVals;
	myNet.getResults(resultVals);
}

