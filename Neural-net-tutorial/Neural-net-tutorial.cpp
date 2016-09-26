// Neural-net-tutorial.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

class Net {
public: 
	Net(topology);
	void feedForward(inputVals);
	void backProp(targetVals);
	//getResults(resultVals) is a function that only reads the result, does not modify the Net object at all. Thats why its const
	void getResults(resultVals) const;

private:


};

int main()
{
	//Construct
	Net myNet(topology);
	
	//Train
	myNet.feedForward(inputVals);
	myNet.backProp(targetVals);
	//Save results
	myNet.getResults(resultVals);
}

