// Neural-net-tutorial.cpp : Defines the entry point for the console application.
//TODO: DONT FORGET BIAS

#include "stdafx.h"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <cstdlib>


struct Connection {
	double weight, deltaWeight;
};

class Neuron;
typedef std::vector<Neuron> Layer;


//****************** class Neuron ******************
class Neuron {
public:
	Neuron(unsigned numOutputs);


private:
	static double randomWeight(void) { return rand() / double(RAND_MAX); } //This returns random value between 0 and 1
	double m_outputVal;
	double m_bias;
	std::vector<Connection> m_outputWeights;

};

Neuron::Neuron(unsigned numOutputs) {
	for (unsigned c = 0; c < numOutputs; ++c) {
		m_outputWeights.push_back(Connection());
		m_outputWeights.back().weight = randomWeight();

	}
};
//****************** class Net ******************
class Net {
public:
	Net(const std::vector<unsigned> &topology);
	void feedForward(const std::vector<double> &inputVals) {};
	void backProp(const std::vector<double> &targetVals) {};
	//getResults(resultVals) is a function that only reads the result, does not modify the Net object at all. Thats why its const
	void getResults(std::vector<double> &resultVals) const {};

private:
	std::vector<Layer> m_layers; //m_layers[layerNum][neuronNum]
};

//Define the Net constructor:
Net::Net(const std::vector<unsigned> &topology) {
	unsigned numLayers = topology.size();
	for (unsigned layerNum = 0; layerNum < numLayers; ++layerNum) {
		//append
		m_layers.push_back(Layer());
		unsigned numOutputs = layerNum == topology.size() - 1 ? 0 : topology[layerNum + 1];

		//We made layer: add neurons to layer
		for (unsigned neuronNum = 0; neuronNum < topology[layerNum]; ++neuronNum) {

			//get last element in m_layers (Layer) and append a new neuron to it.
			m_layers.back().push_back(Neuron(numOutputs));
			std::cout << "Made a Neuron!" << std::endl;
		}

	}

}

int main()
{
	// e.g., {3, 2, 1} (3 layers, 1st - 3 nodes, 2nd - 2 nodes, 3rd - 1 node
	std::vector<unsigned> topology;
	topology.push_back(3);
	topology.push_back(2);
	topology.push_back(1);

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

