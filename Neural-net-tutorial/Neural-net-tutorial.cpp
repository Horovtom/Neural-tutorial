// Neural-net-tutorial.cpp : Defines the entry point for the console application.
//TODO: DONT FORGET BIAS

#include "stdafx.h"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <cassert>
#include <cmath>

struct Connection {
    double weight, deltaWeight;
};

class Neuron;

typedef std::vector<Neuron> Layer;


//****************** class Neuron ******************
class Neuron {
public:
    Neuron(unsigned numOutputs, unsigned myIndex);

    void setOutputVal(const double val) { m_outputVal = val; }

    double getOutputVal(void) const { return m_outputVal; }

    void feedForward(const Layer &prevLayer);

private:
    static double randomWeight(void) { return rand() / double(RAND_MAX); } //This returns random value between 0 and 1
    double m_outputVal;

    static double activationFunction(double x);

    static double activationFunctionDerivative(double x);

    double m_bias;
    std::vector<Connection> m_outputWeights;
    unsigned m_myIndex;
};

//This is gonna be a hyperbolic tangent function
double Neuron::activationFunction(double x) {
    // tanh - output range [-1,0...1,0]
    return tanh(x);
}

double Neuron::activationFunctionDerivative(double x) {
    // tanh derivative
    double tan = tanh(x);
    return 1.0 - tan * tan;
}

void Neuron::feedForward(const Layer &prevLayer) {
    double sum = 0.0;

    //Sum outputs from previous layer and add bias to it
    for (unsigned n = 0; n < prevLayer.size(); ++n) {
        sum += prevLayer[n].getOutputVal() * prevLayer[n].m_outputWeights[m_myIndex].weight;
    }

    sum += m_bias;

    //Applying a static activationfunction
    m_outputVal = Neuron::activationFunction(sum);
}

Neuron::Neuron(unsigned numOutputs, unsigned myIndex) {
    m_myIndex = myIndex;
    for (unsigned c = 0; c < numOutputs; ++c) {
        m_outputWeights.push_back(Connection());
        m_outputWeights.back().weight = randomWeight();

    }
}

//****************** class Net ******************
class Net {
public:
    Net(const std::vector<unsigned> &topology);

    void feedForward(const std::vector<double> &inputVals);

    void backProp(const std::vector<double> &targetVals);

    //getResults(resultVals) is a function that only reads the result, does not modify the Net object at all. Thats why its const
    void getResults(std::vector<double> &resultVals) const {};

private:
    std::vector<Layer> m_layers; //m_layers[layerNum][neuronNum]
    double m_error;
    double m_recentAverageError;
    double m_recentAverageSmoothingFactor;
};

void Net::backProp(const std::vector<double> &targetVals) {
    //Calculate overall net error (RMS - root mean square)
    Layer &outputLayer = m_layers.back();
    m_error = 0.0;

    for (unsigned i = 0; i < outputLayer.size(); ++i) {
        double delta = targetVals[i] - outputLayer[i].getOutputVal();
        m_error += delta * delta;
    }
    m_error /= outputLayer.size();
    m_error = sqrt(m_error);

    //Implement a recent average measurement:

    m_recentAverageError = (m_recentAverageError * m_recentAverageSmoothingFactor + m_error) / (m_recentAverageSmoothingFactor + 1.0);

    //Calculate output layer gradients

    //Calculat gradients on hidden layers

    //For all layers from outputs to first hidden layer update connection weights
}

void Net::feedForward(const std::vector<double> &inputVals) {
    assert(inputVals.size() == m_layers[0].size());

    //feed the first layer (input neurons)
    for (unsigned i = 0; i < inputVals.size(); ++i) {
        m_layers[0][i].setOutputVal(inputVals[i]);
    }

    //Forward propagate
    for (unsigned layerNum = 1; layerNum < m_layers.size(); ++layerNum) {
        Layer &prevLayer = m_layers[layerNum - 1];
        for (unsigned n = 0; n < m_layers[layerNum].size(); ++n) {
            m_layers[layerNum][n].feedForward(prevLayer);
        }
    }
}

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
            m_layers.back().push_back(Neuron(numOutputs, neuronNum));
            std::cout << "Made a Neuron!" << std::endl;
        }

    }

}

int main() {
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

